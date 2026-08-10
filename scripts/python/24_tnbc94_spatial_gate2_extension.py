#!/usr/bin/env python3
"""Extend the frozen Gate 2 neighborhood test to 94 annotated TNBC arrays."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest
from sklearn.neighbors import NearestNeighbors
import statsmodels.api as sm


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "data" / "derived" / "tnbc94_annotated_spot_scores"
TABLES = ROOT / "results" / "tables"
METADATA = ROOT / "metadata"
TARGETS = [
    "fibrovascular", "fibroblast_ecm", "endothelial_angiogenic",
    "malignant_epithelial", "cd8_cytotoxic", "treg",
]
N_PERMUTATIONS = 2_000
SEED = 20260721


def zscore(values: np.ndarray) -> np.ndarray:
    sd = np.std(values, ddof=1)
    return np.zeros_like(values, dtype=float) if not np.isfinite(sd) or sd == 0 else (
        values - np.mean(values)
    ) / sd


def bh(values: pd.Series) -> np.ndarray:
    p = values.to_numpy(float)
    order = np.argsort(p)
    adjusted = p[order] * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.minimum(adjusted, 1.0)
    return output


def permuted_high_matrix(high: np.ndarray, strata: list[np.ndarray], rng) -> np.ndarray:
    # Float64 avoids platform-specific overflow warnings observed in mixed
    # float32/float64 BLAS matrix multiplication during the initial dry run.
    output = np.zeros((N_PERMUTATIONS, len(high)), dtype=np.float64)
    for indices in strata:
        count_high = int(high[indices].sum())
        if count_high == 0:
            continue
        if count_high == len(indices):
            output[:, indices] = 1
            continue
        random_values = rng.random((N_PERMUTATIONS, len(indices)))
        selected = np.argpartition(random_values, count_high - 1, axis=1)[:, :count_high]
        rows = np.arange(N_PERMUTATIONS)[:, None]
        output[rows, indices[selected]] = 1
    return output


def main() -> None:
    effect_rows = []
    null_by_patient = []
    spot_outputs = []
    files = sorted(INPUT.glob("TNBC*_annotated_scores.tsv"), key=lambda p: int(p.stem.split("_")[0][4:]))
    for patient_index, path in enumerate(files):
        scores = pd.read_csv(path, sep="\t")
        scores = scores.loc[
            scores["nFeature_RNA"].ge(200) & ~scores["Classification"].eq("Artefacts")
        ].reset_index(drop=True)
        if len(scores) < 50:
            continue
        design = sm.add_constant(pd.DataFrame({
            "macrophage_identity": scores["macrophage_identity"],
            "hypoxia": scores["hypoxia"],
            "log1p_umi": np.log1p(scores["nCount_RNA"]),
        }), has_constant="add")
        model = sm.OLS(scores["P4_primary"], design).fit()
        scores["P4_residual"] = model.resid
        scores["fibrovascular"] = (
            zscore(scores["fibroblast_ecm"].to_numpy())
            + zscore(scores["endothelial_angiogenic"].to_numpy())
        ) / 2
        coordinates = scores[["x", "y"]].to_numpy(float)
        neighbors = NearestNeighbors(n_neighbors=7).fit(coordinates).kneighbors(
            coordinates, return_distance=False
        )[:, 1:]
        neighbor_targets = np.column_stack([
            scores[target].to_numpy()[neighbors].mean(axis=1) for target in TARGETS
        ])
        residual = scores["P4_residual"].to_numpy()
        high = residual >= np.quantile(residual, 0.75)
        observed = neighbor_targets[high].mean(axis=0) - neighbor_targets[~high].mean(axis=0)

        labels = scores["Classification"].astype(str).copy()
        label_counts = labels.value_counts()
        labels = labels.where(labels.map(label_counts).ge(20), "RARE")
        strata = [np.flatnonzero(labels.to_numpy() == label) for label in labels.unique()]
        rng = np.random.default_rng(SEED + patient_index)
        high_permuted = permuted_high_matrix(high, strata, rng)
        high_count = high_permuted.sum(axis=1)[:, None]
        # Use an explicit contraction rather than the platform BLAS matmul,
        # which emits spurious floating-point warnings for these thin matrices.
        high_sum = np.einsum("pn,nt->pt", high_permuted, neighbor_targets, optimize=False)
        total_sum = neighbor_targets.sum(axis=0)[None, :]
        null = high_sum / high_count - (total_sum - high_sum) / (len(high) - high_count)
        null_by_patient.append(null)

        patient = scores["patient"].astype(str).iloc[0]
        for target_index, target in enumerate(TARGETS):
            null_center = np.median(null[:, target_index])
            if target == "fibrovascular":
                p_value = (1 + np.sum(null[:, target_index] >= observed[target_index])) / (
                    N_PERMUTATIONS + 1
                )
                alternative = "greater"
            else:
                p_value = (
                    1 + np.sum(
                        np.abs(null[:, target_index] - null_center)
                        >= abs(observed[target_index] - null_center)
                    )
                ) / (N_PERMUTATIONS + 1)
                alternative = "two-sided"
            effect_rows.append({
                "patient": patient, "array": scores["array"].iloc[0], "target": target,
                "spots": len(scores), "high_P4_spots": int(high.sum()),
                "observed_neighbor_effect": observed[target_index], "null_median": null_center,
                "excess_over_null": observed[target_index] - null_center,
                "permutation_p": p_value, "alternative": alternative,
                "P4_residual_r_squared": model.rsquared,
            })
        spot_outputs.append(scores.assign(high_P4_residual=high))

    effects = pd.DataFrame(effect_rows)
    effects["fdr_within_target"] = effects.groupby("target", observed=True)["permutation_p"].transform(bh)
    effects.to_csv(TABLES / "TNBC94_spatial_neighbor_effects_by_patient.tsv", sep="\t", index=False)

    null_stack = np.stack(null_by_patient, axis=1)
    cohort_rows = []
    for target_index, target in enumerate(TARGETS):
        patient_values = effects.loc[effects["target"].eq(target), "observed_neighbor_effect"].to_numpy()
        excess_values = effects.loc[effects["target"].eq(target), "excess_over_null"].to_numpy()
        cohort_null = null_stack[:, :, target_index].mean(axis=1)
        observed_mean = patient_values.mean()
        null_center = np.median(cohort_null)
        if target == "fibrovascular":
            p_value = (1 + np.sum(cohort_null >= observed_mean)) / (N_PERMUTATIONS + 1)
            alternative = "greater"
        else:
            p_value = (
                1 + np.sum(np.abs(cohort_null - null_center) >= abs(observed_mean - null_center))
            ) / (N_PERMUTATIONS + 1)
            alternative = "two-sided"
        positives = int((excess_values > 0).sum())
        cohort_rows.append({
            "target": target, "patients": len(patient_values),
            "mean_patient_effect": observed_mean, "cohort_null_median": null_center,
            "mean_excess_over_null": excess_values.mean(), "patients_positive_excess": positives,
            "sign_test_p": binomtest(positives, len(patient_values), 0.5).pvalue,
            "cohort_permutation_p": p_value, "alternative": alternative,
        })
    cohort = pd.DataFrame(cohort_rows)
    secondary = ~cohort["target"].eq("fibrovascular")
    cohort.loc[secondary, "secondary_fdr"] = bh(cohort.loc[secondary, "cohort_permutation_p"])
    cohort.to_csv(TABLES / "TNBC94_spatial_neighbor_effects_cohort.tsv", sep="\t", index=False)
    pd.concat(spot_outputs, ignore_index=True).to_csv(
        METADATA / "TNBC94_annotated_frozen_program_spot_scores.tsv.gz",
        sep="\t", index=False, compression="gzip",
    )
    print(cohort.to_string(index=False))


if __name__ == "__main__":
    main()
