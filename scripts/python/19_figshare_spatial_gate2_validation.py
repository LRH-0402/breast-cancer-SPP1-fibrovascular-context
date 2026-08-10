#!/usr/bin/env python3
"""Validate the frozen Gate 2 spatial test in eight independent Visium sections."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
import statsmodels.api as sm


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "data" / "derived" / "figshare_21591429_scores"
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
SECTIONS = [f"P{index}" for index in range(1, 9)]
TARGETS = [
    "fibrovascular", "fibroblast_ecm", "endothelial_angiogenic",
    "malignant_epithelial", "cd8_cytotoxic", "treg",
]
N_PERMUTATIONS = 2_000
SEED = 20260721


def zscore(values: np.ndarray) -> np.ndarray:
    sd = np.std(values, ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros_like(values, dtype=float)
    return (values - np.mean(values)) / sd


def bh(values: pd.Series) -> np.ndarray:
    p = values.to_numpy(float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.minimum(adjusted, 1.0)
    return output


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    effect_rows: list[dict] = []
    null_by_section: dict[str, np.ndarray] = {}
    plot_data = None

    for section_index, section in enumerate(SECTIONS):
        scores = pd.read_csv(INPUT / f"{section}_scores.tsv", sep="\t")
        residual_x = sm.add_constant(
            pd.DataFrame({
                "macrophage_identity": scores["macrophage_identity"],
                "hypoxia": scores["hypoxia"],
                "log1p_umi": np.log1p(scores["nCount_RNA"]),
            }),
            has_constant="add",
        )
        model = sm.OLS(scores["P4_primary"], residual_x).fit()
        scores["P4_residual"] = model.resid
        scores["fibrovascular"] = (
            zscore(scores["fibroblast_ecm"].to_numpy())
            + zscore(scores["endothelial_angiogenic"].to_numpy())
        ) / 2

        coordinates = scores[["x", "y"]].to_numpy(float)
        neighbors = (
            NearestNeighbors(n_neighbors=7)
            .fit(coordinates)
            .kneighbors(coordinates, return_distance=False)[:, 1:]
        )
        neighbor_targets = np.column_stack([
            scores[target].to_numpy()[neighbors].mean(axis=1) for target in TARGETS
        ])
        residual = scores["P4_residual"].to_numpy()
        high = residual >= np.quantile(residual, 0.75)
        observed = neighbor_targets[high].mean(axis=0) - neighbor_targets[~high].mean(axis=0)

        # Locked fallback for sections without pathology labels: preserve the
        # malignant-expression composition by shuffling only within quintiles.
        malignant_rank = scores["malignant_epithelial"].rank(method="first")
        strata_labels = pd.qcut(malignant_rank, q=5, labels=False)
        strata = [np.flatnonzero(strata_labels.to_numpy() == label) for label in range(5)]
        rng = np.random.default_rng(SEED + section_index)
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
        null_by_section[section] = null

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
                "section": section,
                "target": target,
                "spots": len(scores),
                "high_P4_spots": int(high.sum()),
                "observed_neighbor_effect": observed[target_index],
                "null_median": null_center,
                "permutation_p": p_value,
                "alternative": alternative,
                "P4_residual_r_squared": model.rsquared,
            })

        if section == "P8":
            plot_data = scores.assign(high_P4_residual=high)

    effects = pd.DataFrame(effect_rows)
    effects["fdr_within_target"] = effects.groupby("target", observed=True)[
        "permutation_p"
    ].transform(lambda values: bh(values))
    effects.to_csv(
        TABLES / "Figshare_spatial_neighbor_effects_by_section.tsv", sep="\t", index=False
    )

    cohort_rows = []
    for target_index, target in enumerate(TARGETS):
        observed_values = effects.loc[
            effects["target"].eq(target), "observed_neighbor_effect"
        ].to_numpy()
        cohort_observed = observed_values.mean()
        cohort_null = np.column_stack([
            null_by_section[section][:, target_index] for section in SECTIONS
        ]).mean(axis=1)
        null_center = np.median(cohort_null)
        if target == "fibrovascular":
            p_value = (1 + np.sum(cohort_null >= cohort_observed)) / (N_PERMUTATIONS + 1)
            alternative = "greater"
        else:
            p_value = (
                1 + np.sum(
                    np.abs(cohort_null - null_center) >= abs(cohort_observed - null_center)
                )
            ) / (N_PERMUTATIONS + 1)
            alternative = "two-sided"
        cohort_rows.append({
            "target": target,
            "sections": len(SECTIONS),
            "mean_section_effect": cohort_observed,
            "sections_positive": int((observed_values > 0).sum()),
            "cohort_null_median": null_center,
            "cohort_permutation_p": p_value,
            "alternative": alternative,
        })
    cohort = pd.DataFrame(cohort_rows)
    secondary = ~cohort["target"].eq("fibrovascular")
    cohort.loc[secondary, "secondary_fdr"] = bh(
        cohort.loc[secondary, "cohort_permutation_p"]
    )
    cohort.to_csv(
        TABLES / "Figshare_spatial_neighbor_effects_cohort.tsv", sep="\t", index=False
    )

    if plot_data is not None:
        figure, axes = plt.subplots(1, 3, figsize=(12, 3.7), constrained_layout=True)
        for axis, column, title in zip(
            axes,
            ["P4_residual", "fibrovascular", "high_P4_residual"],
            ["P4 residual", "Fibrovascular score", "Top-quartile P4"],
        ):
            plotted = axis.scatter(
                plot_data["y"], plot_data["x"], c=plot_data[column],
                s=5, cmap="viridis", linewidths=0,
            )
            axis.invert_yaxis()
            axis.set_aspect("equal")
            axis.set_title(title)
            axis.axis("off")
            figure.colorbar(plotted, ax=axis, fraction=0.035, pad=0.01)
        figure.savefig(FIGURES / "Figshare_P8_P4_fibrovascular_validation.pdf", dpi=300)
        plt.close(figure)

    print("Cohort validation")
    print(cohort.to_string(index=False))
    print("\nPrimary target by section")
    print(effects.loc[effects["target"].eq("fibrovascular")].to_string(index=False))


if __name__ == "__main__":
    main()
