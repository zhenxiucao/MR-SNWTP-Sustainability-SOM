# -*- coding: utf-8 -*-

import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests


INPUT_FILE = Path(r"D:\00实验数据-重要\00CAS2511\00CAS\09南水北调\03实验\13som聚类\06SOM_EWD\02python\00somDriver1kdayu0(ewd.xlsx")
OUTPUT_DIR = INPUT_FILE.parent
OUTPUT_NAME = "multinomial_SOM1_spatial_cluster_robust.xlsx"

MAX_PER_GROUP = 3000
RANDOM_STATE = 2026
DISTANCE_IN_METERS = True
YEARS = (2015, 2020)
REFERENCE_SOM = "1"
MIN_CLUSTER_WARNING = 30

REQUIRED_COLUMNS = [
    "pixel_id", "year", "SOM_type", "EWD", "distance", "urban_type", "spa_blo_ID"
]
FORMULA = "EWD_z + distance_z + EWD_z:distance_z + C(urban_type) + C(year)"


def standardize_columns(frame):
    df = frame.copy()
    df.columns = df.columns.astype(str).str.strip()
    rename = {
        "som_type": "SOM_type", "SOM": "SOM_type", "type": "SOM_type",
        "Type": "SOM_type", "som": "SOM_type", "ewd": "EWD",
        "Ewd": "EWD", "Distance": "distance", "dist": "distance",
        "DIST": "distance", "Urban_type": "urban_type",
        "urban": "urban_type", "Urban": "urban_type",
        "urbanclass": "urban_type", "urban_class": "urban_type",
        "Year": "year", "spatial_block_id": "spa_blo_ID",
        "SPA_BLO_ID": "spa_blo_ID", "Pixel_ID": "pixel_id",
        "PIXEL_ID": "pixel_id"
    }
    df = df.rename(columns=rename)
    if df.columns.duplicated().any():
        raise ValueError("Duplicate column names after normalization; check the Excel headers.")
    return df


def load_data(input_file):
    if not input_file.is_file():
        raise FileNotFoundError(f"Input file does not exist: {input_file}")
    frames = []
    suffix = input_file.suffix.lower()
    if suffix in (".xlsx", ".xls"):
        xls = pd.ExcelFile(input_file)
        print("Excel sheets:", xls.sheet_names)
        for sheet in xls.sheet_names:
            frame = standardize_columns(pd.read_excel(xls, sheet_name=sheet))
            if "year" not in frame.columns:
                try:
                    frame["year"] = int(sheet)
                except ValueError as exc:
                    raise ValueError(f"Sheet '{sheet}' has no year column or numeric year name.") from exc
            frames.append(frame)
    elif suffix == ".csv":
        frames.append(standardize_columns(pd.read_csv(input_file)))
    else:
        raise ValueError("Input must be .xlsx, .xls, or .csv")
    data = pd.concat(frames, ignore_index=True)
    missing = sorted(set(REQUIRED_COLUMNS) - set(data.columns))
    if missing:
        raise ValueError(f"Missing required fields: {missing}. Required field is spa_blo_ID.")
    return data[REQUIRED_COLUMNS].copy()


def prepare_data(raw):
    df = raw.copy()
    for field in ["pixel_id", "SOM_type", "urban_type", "spa_blo_ID"]:
        df[field] = df[field].astype("string").str.strip()
    df["SOM_type"] = (
        df["SOM_type"].str.upper()
        .str.replace(r"^(SOM|TYPE)\s*", "", regex=True)
        .str.replace(r"\.0+$", "", regex=True)
    )
    df["year"] = pd.to_numeric(df["year"], errors="coerce")
    df["EWD"] = pd.to_numeric(df["EWD"], errors="coerce")
    df["distance"] = pd.to_numeric(df["distance"], errors="coerce")
    df = df.loc[df["year"].isin(YEARS) & df["SOM_type"].isin([str(x) for x in range(1, 7)])].copy()
    df = df.loc[df["urban_type"] != "水体"].copy()
    for field in ["pixel_id", "urban_type", "spa_blo_ID"]:
        df.loc[df[field].isin(["", "nan", "NaN", "None", "<NA>"]), field] = pd.NA
    df = df.dropna(subset=REQUIRED_COLUMNS).copy()
    if df.empty:
        raise ValueError("No valid observations after filtering the required variables.")
    df["year"] = df["year"].astype(int)
    df["distance_km"] = df["distance"] / 1000.0 if DISTANCE_IN_METERS else df["distance"]
    if df.duplicated(["pixel_id", "year"]).any():
        examples = df.loc[df.duplicated(["pixel_id", "year"], keep=False), ["pixel_id", "year"]].head()
        raise ValueError(f"Duplicate (pixel_id, year) combinations found:\n{examples}")
    consistency = df.groupby("pixel_id")["spa_blo_ID"].nunique(dropna=False)
    if (consistency > 1).any():
        examples = consistency.loc[consistency > 1].head().index.tolist()
        raise ValueError(f"The same pixel_id has different spa_blo_ID across years: {examples}")

    counts_full = (
        df.groupby(["year", "SOM_type"], observed=True).size().rename("n_pixels").reset_index()
    )
    counts_full["pct_within_year"] = (
        100.0 * counts_full["n_pixels"] /
        counts_full.groupby("year")["n_pixels"].transform("sum")
    )

    ewd_low, ewd_high = df["EWD"].quantile([0.01, 0.99])
    distance_low, distance_high = df["distance_km"].quantile([0.01, 0.99])
    df = df.loc[
        df["EWD"].between(ewd_low, ewd_high) &
        df["distance_km"].between(distance_low, distance_high)
    ].copy()
    print("Observations after filtering and 1%-99% trimming:", len(df))
    if df.empty:
        raise ValueError("No observations left after 1%-99% trimming.")
    return df, counts_full, {
        "EWD_trim_lower": float(ewd_low), "EWD_trim_upper": float(ewd_high),
        "distance_km_trim_lower": float(distance_low),
        "distance_km_trim_upper": float(distance_high)
    }


def stratified_sample(df):
    parts = []
    for (_, _), group in df.groupby(["year", "SOM_type"], observed=True, sort=True):
        parts.append(group.sample(n=min(len(group), MAX_PER_GROUP), random_state=RANDOM_STATE))
    sample = pd.concat(parts, ignore_index=True)
    if sample["SOM_type"].nunique() != 6:
        raise ValueError("The model subsample must include all six SOM categories.")
    if REFERENCE_SOM not in sample["SOM_type"].unique():
        raise ValueError("SOM1 is not present in the sample.")
    for src, target in (("EWD", "EWD_z"), ("distance_km", "distance_z")):
        sd = sample[src].std()
        if not np.isfinite(sd) or sd == 0:
            raise ValueError(f"Cannot standardize {src}; zero or non-finite standard deviation.")
        sample[target] = (sample[src] - sample[src].mean()) / sd
    sample["year"] = sample["year"].astype(str)
    sample["urban_type"] = sample["urban_type"].astype(str)
    sample["spa_blo_ID"] = sample["spa_blo_ID"].astype(str)
    categories = [REFERENCE_SOM] + [str(x) for x in range(1, 7) if str(x) != REFERENCE_SOM]
    y = pd.Categorical(sample["SOM_type"], categories=categories).codes
    if np.any(y < 0):
        raise ValueError("Some SOM categories were not coded successfully.")
    sample_counts = (
        sample.groupby(["year", "SOM_type"], observed=True).size().rename("n_sample").reset_index()
    )
    return sample, y, categories, sample_counts


def fit_cluster_robust(sample, y):
    X = patsy.dmatrix(FORMULA, data=sample, return_type="dataframe", NA_action="raise")
    groups = sample["spa_blo_ID"].to_numpy()
    n_clusters = pd.Series(groups).nunique()
    if n_clusters < 10:
        raise ValueError(f"Only {n_clusters} spatial clusters: too few for credible cluster-robust inference.")
    if n_clusters < MIN_CLUSTER_WARNING:
        warnings.warn(
            f"Only {n_clusters} spatial clusters. Standard cluster-robust inference may be unreliable; "
            "consider cluster-level resampling or alternative spatial-block sensitivity checks.",
            RuntimeWarning,
        )
    model = sm.MNLogit(y, X)
    fit = model.fit(
        method="lbfgs", maxiter=400, disp=False,
        cov_type="cluster",
        cov_kwds={"groups": groups, "use_correction": True, "df_correction": True}
    )
    if not bool(fit.mle_retvals.get("converged", False)):
        raise RuntimeError("MNLogit optimizer did not converge. Do not report these estimates.")
    if not np.isfinite(fit.params.to_numpy()).all():
        raise RuntimeError("Non-finite model estimates; check the data and convergence.")
    return fit, X, n_clusters


def coefficient_table(fit, X, categories):
    params = fit.params
    ses = fit.bse
    pvalues = fit.pvalues
    zcrit = norm.ppf(0.975)
    rows = []
    for j, column in enumerate(params.columns):
        target_category = categories[j + 1]
        for variable in X.columns:
            coef = float(params.loc[variable, column])
            se = float(ses.loc[variable, column])
            pvalue = float(pvalues.loc[variable, column])
            lo, hi = coef - zcrit * se, coef + zcrit * se
            with np.errstate(over="ignore"):
                rrr, rrr_lo, rrr_hi = np.exp([coef, lo, hi])
            rows.append({
                "SOM_type": f"SOM{target_category}", "reference": "SOM1",
                "variable": variable, "coef": coef, "cluster_SE": se,
                "coef_CI_lower_95": lo, "coef_CI_upper_95": hi,
                "RRR": float(rrr), "RRR_CI_lower_95": float(rrr_lo),
                "RRR_CI_upper_95": float(rrr_hi),
                "p_cluster": pvalue
            })
    table = pd.DataFrame(rows)
    table["q_BH_all_nonintercept_coefficients"] = np.nan
    mask = table["variable"] != "Intercept"
    table.loc[mask, "q_BH_all_nonintercept_coefficients"] = multipletests(
        table.loc[mask, "p_cluster"], method="fdr_bh"
    )[1]
    return table


def robust_wald_tests(fit):
    idx = list(fit.cov_params().index)
    tests = {
        "EWD_joint_including_interaction": lambda v: v == "EWD_z" or v == "EWD_z:distance_z",
        "distance_joint_including_interaction": lambda v: v == "distance_z" or v == "EWD_z:distance_z",
        "EWD_distance_interaction": lambda v: v == "EWD_z:distance_z",
        "urban_type": lambda v: v.startswith("C(urban_type)"),
        "year": lambda v: v.startswith("C(year)"),
    }
    rows = []
    for test_name, selector in tests.items():
        positions = [i for i, (_, v) in enumerate(idx) if selector(v)]
        if not positions:
            continue
        R = np.eye(len(idx))[positions]
        wt = fit.wald_test(R, scalar=True)
        rows.append({
            "test": test_name, "Wald_chi2": float(wt.statistic),
            "df": int(len(positions)), "p_cluster_Wald": float(wt.pvalue),
            "covariance": "spatial-block cluster-robust",
        })
    if not rows:
        raise RuntimeError("No matching terms found for robust Wald tests.")
    table = pd.DataFrame(rows)
    table["q_BH_five_tests"] = multipletests(table["p_cluster_Wald"], method="fdr_bh")[1]
    return table


def build_model_info(fit, sample, categories, n_clusters, cutoffs, input_file):
    info = {
        "input_file": str(input_file),
        "n_observations": int(fit.nobs),
        "n_unique_pixels_in_sample": int(sample["pixel_id"].nunique()),
        "n_spatial_clusters": int(n_clusters),
        "n_SOM_types": len(categories),
        "reference_SOM_type": "SOM1",
        "SOM_category_order": ", ".join("SOM" + c for c in categories),
        "log_likelihood": float(fit.llf),
        "pseudo_R2_McFadden": float(fit.prsquared),
        "AIC": float(fit.aic),
        "BIC": float(fit.bic),
        "covariance_type": fit.cov_type,
        "covariance_groups": "spa_blo_ID (same block across years)",
        "cluster_small_sample_correction": True,
        "confidence_interval_type": "95% normal-approximation cluster-robust CI",
        "multiple_comparisons": "Benjamini-Hochberg FDR in coefficient and joint Wald outputs",
        "sampling_method": "stratified by year and SOM_type",
        "max_per_year_SOM_type": MAX_PER_GROUP,
        "random_seed": RANDOM_STATE,
        "outlier_trim": "EWD and distance 1st-99th percentile before sampling",
        "distance_original_units": "metres" if DISTANCE_IN_METERS else "kilometres",
        "years": ", ".join(map(str, YEARS)),
        "model_formula": FORMULA,
        "converged": bool(fit.mle_retvals.get("converged", False)),
        "model_fit_scope": "stratified subsample; not full-grid population",
    }
    info.update(cutoffs)
    return pd.DataFrame({"item": list(info.keys()), "value": list(info.values())})


def main():
    input_file = Path(sys.argv[1]) if len(sys.argv) > 1 else INPUT_FILE
    output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else (input_file.parent if len(sys.argv) > 1 else OUTPUT_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)
    out_file = output_dir / OUTPUT_NAME

    raw = load_data(input_file)
    print("Raw rows:", len(raw))
    trimmed, counts_full, cutoffs = prepare_data(raw)
    sample, y, categories, sample_counts = stratified_sample(trimmed)
    print("Sample observations:", len(sample))
    print("Reference category:", "SOM1 (code 0)")
    fit, X, n_clusters = fit_cluster_robust(sample, y)
    print("Spatial clusters:", n_clusters)
    print("McFadden pseudo-R2:", round(float(fit.prsquared), 6))

    coefs = coefficient_table(fit, X, categories)
    tests = robust_wald_tests(fit)
    class_mapping = pd.DataFrame({
        "code": range(len(categories)),
        "SOM_type": ["SOM" + c for c in categories],
        "reference": [c == REFERENCE_SOM for c in categories]
    })
    block_summary = (
        sample.groupby("spa_blo_ID").size().rename("sampled_n").reset_index()
    )
    model_info = build_model_info(fit, sample, categories, n_clusters, cutoffs, input_file)
    with pd.ExcelWriter(out_file, engine="openpyxl") as writer:
        model_info.to_excel(writer, sheet_name="model_info", index=False)
        class_mapping.to_excel(writer, sheet_name="class_mapping", index=False)
        coefs.to_excel(writer, sheet_name="coefficients_cluster", index=False)
        tests.to_excel(writer, sheet_name="Wald_cluster", index=False)
        counts_full.to_excel(writer, sheet_name="SOM_input_counts", index=False)
        sample_counts.to_excel(writer, sheet_name="sample_counts", index=False)
        block_summary.to_excel(writer, sheet_name="cluster_sizes", index=False)
    print("Cluster-robust Wald tests:")
    print(tests.to_string(index=False))
    print("Saved:", out_file)


if __name__ == "__main__":
    main()
