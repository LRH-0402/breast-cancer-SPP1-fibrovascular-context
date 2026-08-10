#!/usr/bin/env python3
"""Reviewer-requested spatial sensitivity and leave-one-P4-gene audits."""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
import statsmodels.api as sm
import yaml


ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "results" / "tables"
SEED = 20260722
N_PERMUTATIONS = 250


def zscore(values: np.ndarray) -> np.ndarray:
    sd = np.std(values, ddof=1)
    return np.zeros_like(values, dtype=float) if not np.isfinite(sd) or sd == 0 else (
        values - np.mean(values)
    ) / sd


def load_units() -> list[tuple[str, str, pd.DataFrame, str]]:
    units = []
    wu = pd.read_csv(ROOT / "metadata" / "Wu_Visium_frozen_program_spot_scores.tsv.gz", sep="\t")
    for unit, data in wu.groupby("patientid", observed=True):
        units.append(("Wu discovery", str(unit), data.reset_index(drop=True), "pathology"))
    for path in sorted((ROOT / "data" / "derived" / "figshare_21591429_scores").glob("P*_scores.tsv")):
        data = pd.read_csv(path, sep="\t")
        units.append(("Figshare validation", path.stem.split("_")[0], data, "malignant_quintile"))
    paths = sorted(
        (ROOT / "data" / "derived" / "tnbc94_annotated_spot_scores").glob("TNBC*_annotated_scores.tsv"),
        key=lambda p: int(p.stem.split("_")[0][4:]),
    )
    for path in paths:
        data = pd.read_csv(path, sep="\t")
        data = data.loc[data["nFeature_RNA"].ge(200) & ~data["Classification"].eq("Artefacts")].reset_index(drop=True)
        if len(data) >= 50:
            units.append(("TNBC94 extension", str(data["patient"].iloc[0]), data, "pathology"))
    return units


def strata_indices(data: pd.DataFrame, mode: str, stratified: bool) -> list[np.ndarray]:
    if not stratified:
        return [np.arange(len(data))]
    if mode == "pathology":
        labels = data["Classification"].astype(str).copy()
        counts = labels.value_counts()
        labels = labels.where(labels.map(counts).ge(20), "RARE")
    else:
        ranks = data["malignant_epithelial"].rank(method="first")
        labels = pd.qcut(ranks, q=5, labels=False, duplicates="drop").astype(str)
    return [np.flatnonzero(labels.to_numpy() == label) for label in labels.unique()]


def spatial_excess(
    data: pd.DataFrame,
    exposure: np.ndarray,
    mode: str,
    *,
    k: int = 6,
    high_fraction: float = 0.25,
    covariates: tuple[str, ...] = ("macrophage_identity", "hypoxia", "log1p_umi"),
    stratified: bool = True,
    seed: int = SEED,
) -> float:
    design_columns = {}
    for covariate in covariates:
        if covariate == "log1p_umi":
            design_columns[covariate] = np.log1p(data["nCount_RNA"].to_numpy())
        elif covariate == "fibroblast_ecm":
            design_columns[covariate] = data["fibroblast_ecm"].to_numpy()
        else:
            design_columns[covariate] = data[covariate].to_numpy()
    design = sm.add_constant(pd.DataFrame(design_columns), has_constant="add")
    residual = sm.OLS(exposure, design).fit().resid
    target = (zscore(data["fibroblast_ecm"].to_numpy()) + zscore(data["endothelial_angiogenic"].to_numpy())) / 2
    coordinates = data[["array_row", "array_col"]].to_numpy(float) if "array_row" in data else data[["x", "y"]].to_numpy(float)
    neighbors = NearestNeighbors(n_neighbors=min(k + 1, len(data))).fit(coordinates).kneighbors(
        coordinates, return_distance=False
    )[:, 1:]
    neighbor_target = target[neighbors].mean(axis=1)
    cutoff = np.quantile(residual, 1 - high_fraction)
    high = residual >= cutoff
    observed = neighbor_target[high].mean() - neighbor_target[~high].mean()
    strata = strata_indices(data, mode, stratified)
    rng = np.random.default_rng(seed)
    null = np.empty(N_PERMUTATIONS)
    for permutation in range(N_PERMUTATIONS):
        permuted = residual.copy()
        for indices in strata:
            permuted[indices] = rng.permutation(permuted[indices])
        permuted_high = permuted >= np.quantile(permuted, 1 - high_fraction)
        null[permutation] = neighbor_target[permuted_high].mean() - neighbor_target[~permuted_high].mean()
    return float(observed - np.median(null))


def summarize(rows: list[dict], group_columns: list[str]) -> pd.DataFrame:
    data = pd.DataFrame(rows)
    return data.groupby(group_columns, observed=True, as_index=False).agg(
        units=("unit", "nunique"),
        positive_units=("excess_over_null", lambda x: int((x > 0).sum())),
        positive_fraction=("excess_over_null", lambda x: float((x > 0).mean())),
        median_excess=("excess_over_null", "median"),
        mean_excess=("excess_over_null", "mean"),
    )


def main() -> None:
    with (ROOT / "config" / "locked_macrophage_programs_v1.yml").open() as handle:
        genes = yaml.safe_load(handle)["primary_program"]["genes"]
    units = load_units()
    missing = [(cohort, unit, gene) for cohort, unit, data, _ in units for gene in genes if f"P4gene_{gene}" not in data]
    if missing:
        raise RuntimeError(f"Individual P4 gene scores are missing; rerun spatial scoring first: {missing[:5]}")

    settings = []
    for k in (4, 6, 8, 12):
        settings.append((f"neighbors_k{k}", k, 0.25, ("macrophage_identity", "hypoxia", "log1p_umi"), True))
    for fraction in (0.20, 0.25, 0.33):
        settings.append((f"high_fraction_{fraction:.2f}", 6, fraction, ("macrophage_identity", "hypoxia", "log1p_umi"), True))
    settings.extend([
        ("residual_macrophage_depth", 6, 0.25, ("macrophage_identity", "log1p_umi"), True),
        ("residual_depth_only", 6, 0.25, ("log1p_umi",), True),
        ("residual_plus_fibroblast", 6, 0.25, ("macrophage_identity", "hypoxia", "log1p_umi", "fibroblast_ecm"), True),
        ("unstratified_permutation", 6, 0.25, ("macrophage_identity", "hypoxia", "log1p_umi"), False),
    ])
    sensitivity_rows = []
    for unit_index, (cohort, unit, data, mode) in enumerate(units):
        exposure = data["P4_primary"].to_numpy(float)
        for setting, k, fraction, covariates, stratified in settings:
            sensitivity_rows.append({
                "cohort": cohort,
                "unit": unit,
                "setting": setting,
                "neighbors": k,
                "high_fraction": fraction,
                "covariates": ";".join(covariates),
                "permutation_stratified": stratified,
                "excess_over_null": spatial_excess(
                    data, exposure, mode, k=k, high_fraction=fraction, covariates=covariates,
                    stratified=stratified, seed=SEED + unit_index,
                ),
            })
    pd.DataFrame(sensitivity_rows).to_csv(TABLES / "spatial_parameter_sensitivity_by_unit.tsv", sep="\t", index=False)
    summarize(sensitivity_rows, ["cohort", "setting", "neighbors", "high_fraction", "covariates", "permutation_stratified"]).to_csv(
        TABLES / "spatial_parameter_sensitivity_summary.tsv", sep="\t", index=False
    )

    loo_rows = []
    for unit_index, (cohort, unit, data, mode) in enumerate(units):
        for omitted in ["none", *genes]:
            retained = genes if omitted == "none" else [gene for gene in genes if gene != omitted]
            exposure = data[[f"P4gene_{gene}" for gene in retained]].mean(axis=1).to_numpy(float)
            loo_rows.append({
                "cohort": cohort,
                "unit": unit,
                "omitted_gene": omitted,
                "genes_retained": len(retained),
                "excess_over_null": spatial_excess(data, exposure, mode, seed=SEED + unit_index),
            })
    pd.DataFrame(loo_rows).to_csv(TABLES / "spatial_P4_leave_one_gene_out_by_unit.tsv", sep="\t", index=False)
    summarize(loo_rows, ["cohort", "omitted_gene", "genes_retained"]).to_csv(
        TABLES / "spatial_P4_leave_one_gene_out_summary.tsv", sep="\t", index=False
    )
    print(f"Audited {len(units)} spatial units across three cohorts")


if __name__ == "__main__":
    main()
