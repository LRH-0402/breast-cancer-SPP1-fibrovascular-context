#!/usr/bin/env python3
"""Post hoc spatial T-cell-state boundary analysis across three cohorts."""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest
from sklearn.neighbors import NearestNeighbors
import statsmodels.api as sm


ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "results" / "tables"
WU = ROOT / "metadata" / "Wu_Visium_frozen_program_spot_scores.tsv.gz"
FIGSHARE = ROOT / "data" / "derived" / "figshare_21591429_scores"
TNBC94 = ROOT / "data" / "derived" / "tnbc94_annotated_spot_scores"
TARGETS = [
    "cd8_cytotoxic",
    "tcell_exhaustion",
    "tcell_ifng_response",
    "tcell_progenitor",
]
N_PERMUTATIONS = 2_000
SEED = 20260723


def bh(values):
    values = np.asarray(values, float)
    order = np.argsort(values)
    adjusted = values[order] * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.minimum(adjusted, 1)
    return output


def residualize(scores):
    design = sm.add_constant(pd.DataFrame({
        "macrophage_identity": scores["macrophage_identity"],
        "hypoxia": scores["hypoxia"],
        "log1p_umi": np.log1p(scores["nCount_RNA"]),
    }), has_constant="add")
    return sm.OLS(scores["P4_primary"], design).fit().resid


def analyze_unit(scores, coordinate_columns, strata_labels, rng):
    residual = np.asarray(residualize(scores))
    high = residual >= np.quantile(residual, 0.75)
    coordinates = scores[list(coordinate_columns)].to_numpy(float)
    neighbors = NearestNeighbors(n_neighbors=7).fit(coordinates).kneighbors(
        coordinates, return_distance=False
    )[:, 1:]
    neighbor_targets = np.column_stack([
        scores[target].to_numpy(float)[neighbors].mean(axis=1) for target in TARGETS
    ])
    observed = neighbor_targets[high].mean(axis=0) - neighbor_targets[~high].mean(axis=0)

    labels = pd.Series(strata_labels).astype(str)
    counts = labels.value_counts()
    labels = labels.where(labels.map(counts).ge(20), "RARE")
    strata = [np.flatnonzero(labels.to_numpy() == value) for value in labels.unique()]
    null = np.empty((N_PERMUTATIONS, len(TARGETS)))
    for permutation in range(N_PERMUTATIONS):
        permuted = residual.copy()
        for indices in strata:
            permuted[indices] = rng.permutation(permuted[indices])
        permuted_high = permuted >= np.quantile(permuted, 0.75)
        null[permutation] = (
            neighbor_targets[permuted_high].mean(axis=0)
            - neighbor_targets[~permuted_high].mean(axis=0)
        )
    return observed, null, int(high.sum())


def cohort_units():
    wu = pd.read_csv(WU, sep="\t")
    for patient, scores in wu.groupby("patientid", sort=True):
        yield "Wu discovery", str(patient), scores.reset_index(drop=True), (
            "array_row", "array_col"
        ), scores["Classification"].reset_index(drop=True)

    for path in sorted(FIGSHARE.glob("P*_scores.tsv"), key=lambda p: int(p.stem.split("_")[0][1:])):
        scores = pd.read_csv(path, sep="\t")
        ranks = scores["malignant_epithelial"].rank(method="first")
        strata = pd.qcut(ranks, q=5, labels=False)
        yield "Figshare validation", path.stem.split("_")[0], scores, ("x", "y"), strata

    paths = sorted(TNBC94.glob("TNBC*_annotated_scores.tsv"),
                   key=lambda p: int(p.stem.split("_")[0][4:]))
    for path in paths:
        scores = pd.read_csv(path, sep="\t")
        scores = scores.loc[
            scores["nFeature_RNA"].ge(200) & ~scores["Classification"].eq("Artefacts")
        ].reset_index(drop=True)
        if len(scores) < 50:
            continue
        yield "TNBC94 extension", str(scores["patient"].iloc[0]), scores, (
            "x", "y"
        ), scores["Classification"]


def main():
    rows = []
    nulls = {}
    unit_counts = {}
    for index, (cohort, unit, scores, coordinates, strata) in enumerate(cohort_units()):
        missing = [target for target in TARGETS if target not in scores]
        if missing:
            raise KeyError(f"{cohort} {unit} missing scores: {missing}")
        observed, null, high_count = analyze_unit(
            scores, coordinates, strata, np.random.default_rng(SEED + index)
        )
        nulls.setdefault(cohort, []).append(null)
        unit_counts[cohort] = unit_counts.get(cohort, 0) + 1
        for target_index, target in enumerate(TARGETS):
            center = np.median(null[:, target_index])
            p_value = (
                1
                + np.sum(
                    np.abs(null[:, target_index] - center)
                    >= abs(observed[target_index] - center)
                )
            ) / (N_PERMUTATIONS + 1)
            rows.append({
                "cohort": cohort,
                "unit": unit,
                "target": target,
                "spots": len(scores),
                "high_P4_spots": high_count,
                "observed_neighbor_effect": observed[target_index],
                "null_median": center,
                "excess_over_null": observed[target_index] - center,
                "permutation_p": p_value,
            })

    by_unit = pd.DataFrame(rows)
    by_unit["fdr_within_unit"] = by_unit.groupby(
        ["cohort", "unit"], observed=True
    )["permutation_p"].transform(bh)
    by_unit.to_csv(TABLES / "spatial_tcell_state_boundary_by_unit.tsv", sep="\t", index=False)

    summary_rows = []
    for cohort in unit_counts:
        cohort_null = np.stack(nulls[cohort], axis=1)
        for target_index, target in enumerate(TARGETS):
            selected = by_unit.loc[
                by_unit["cohort"].eq(cohort) & by_unit["target"].eq(target)
            ]
            values = selected["observed_neighbor_effect"].to_numpy()
            excess = selected["excess_over_null"].to_numpy()
            synchronized = cohort_null[:, :, target_index].mean(axis=1)
            center = np.median(synchronized)
            observed_mean = values.mean()
            p_value = (
                1
                + np.sum(np.abs(synchronized - center) >= abs(observed_mean - center))
            ) / (N_PERMUTATIONS + 1)
            positives = int((excess > 0).sum())
            summary_rows.append({
                "cohort": cohort,
                "target": target,
                "spatial_units": len(selected),
                "positive_units": positives,
                "positive_fraction": positives / len(selected),
                "mean_excess_over_null": excess.mean(),
                "cohort_permutation_p": p_value,
                "directional_sign_test_p": binomtest(
                    positives, len(selected), 0.5
                ).pvalue,
            })
    summary = pd.DataFrame(summary_rows)
    summary["fdr_within_cohort"] = summary.groupby("cohort", observed=True)[
        "cohort_permutation_p"
    ].transform(bh)
    summary.to_csv(TABLES / "spatial_tcell_state_boundary_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
