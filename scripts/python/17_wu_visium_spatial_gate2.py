#!/usr/bin/env python3
"""Run the frozen Gate 2 discovery test in the Wu Visium cohort."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from sklearn.neighbors import NearestNeighbors
import statsmodels.api as sm
import yaml


ROOT = Path(__file__).resolve().parents[2]
MATRIX_ROOT = ROOT / "data" / "interim" / "spatial_audit" / "wu_visium" / "filtered_count_matrices"
SPATIAL_ROOT = ROOT / "data" / "interim" / "spatial_audit" / "wu_visium" / "spatial"
METADATA_ROOT = ROOT / "data" / "interim" / "spatial_audit" / "wu_metadata" / "metadata"
LOCK_FILE = ROOT / "config" / "locked_macrophage_programs_v1.yml"
SPEC_FILE = ROOT / "config" / "spatial_analysis_spec_v1.yml"
TCELL_SPEC_FILE = ROOT / "config" / "tcell_state_sensitivity_v1.yml"
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
METADATA = ROOT / "metadata"
LOGS = ROOT / "logs"
PATIENTS = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]
N_PERMUTATIONS = 2_000
SEED = 20260721


def read_plain_lines(path: Path) -> list[str]:
    with path.open() as handle:
        return [line.rstrip("\n") for line in handle]


def score_module(matrix: sparse.csr_matrix, gene_lookup: dict[str, int], genes: list[str]):
    present = [gene for gene in genes if gene in gene_lookup]
    if not present:
        raise ValueError(f"No genes present for {genes}")
    columns = [gene_lookup[gene] for gene in present]
    return np.asarray(matrix[:, columns].mean(axis=1)).ravel(), present


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
    with LOCK_FILE.open() as handle:
        lock = yaml.safe_load(handle)
    with SPEC_FILE.open() as handle:
        specification = yaml.safe_load(handle)
    with TCELL_SPEC_FILE.open() as handle:
        tcell_specification = yaml.safe_load(handle)

    modules = {
        "P4_primary": lock["primary_program"]["genes"],
        "P3_secondary": lock["secondary_program"]["genes"],
        "macrophage_identity": ["C1QA", "C1QB", "C1QC", "CD68", "MSR1"],
        "hypoxia": [
            "HIF1A", "EGLN1", "CA9", "VEGFA", "LDHA", "PDK1", "SLC2A1",
            "BNIP3", "NDRG1", "ENO1", "ALDOA", "PGK1", "HK2", "PFKP",
        ],
        **specification["spatial_signatures"],
        **tcell_specification["tcell_state_signatures"],
    }
    target_names = [
        "fibrovascular", "fibroblast_ecm", "endothelial_angiogenic",
        "malignant_epithelial", "cd8_cytotoxic", "treg",
    ]

    patient_effect_rows = []
    spot_tables = []
    null_by_patient: dict[str, np.ndarray] = {}
    coverage_rows = []

    for patient_index, patient in enumerate(PATIENTS):
        matrix_dir = MATRIX_ROOT / f"{patient}_filtered_count_matrix"
        with (matrix_dir / "matrix.mtx.gz").open("rb") as handle:
            counts_gene_by_spot = mmread(handle).tocsr()
        genes = read_plain_lines(matrix_dir / "features.tsv.gz")
        barcodes = read_plain_lines(matrix_dir / "barcodes.tsv.gz")
        counts = counts_gene_by_spot.T.tocsr().astype(np.float64)
        gene_lookup = {gene: index for index, gene in enumerate(genes)}

        metadata = pd.read_csv(METADATA_ROOT / f"{patient}_metadata.csv", index_col=0)
        metadata = metadata.reindex(barcodes)
        positions = pd.read_csv(
            SPATIAL_ROOT / f"{patient}_spatial" / "tissue_positions_list.csv",
            header=None,
            names=["barcode", "in_tissue", "array_row", "array_col", "pixel_row", "pixel_col"],
        ).set_index("barcode").reindex(barcodes)
        if metadata.isna().all(axis=1).any() or positions[["array_row", "array_col"]].isna().any().any():
            raise ValueError(f"Metadata or coordinates failed to align for {patient}")

        keep = metadata["nFeature_RNA"].ge(specification["spot_filter"]["minimum_detected_genes"])
        for label in specification["spot_filter"]["exclude_pathology_labels_containing"]:
            keep &= ~metadata["Classification"].astype(str).str.contains(label, case=False, regex=False)
        keep_indices = np.flatnonzero(keep.to_numpy())
        counts = counts[keep_indices]
        metadata = metadata.iloc[keep_indices].copy()
        positions = positions.iloc[keep_indices].copy()

        library = np.asarray(counts.sum(axis=1)).ravel()
        normalized = counts.multiply((10_000 / library)[:, None]).tocsr()
        normalized.data = np.log1p(normalized.data)

        scores = metadata[["patientid", "subtype", "Classification", "nCount_RNA", "nFeature_RNA"]].copy()
        scores.index.name = "barcode"
        for module, requested_genes in modules.items():
            values, present = score_module(normalized, gene_lookup, requested_genes)
            scores[module] = values
            coverage_rows.append(
                {
                    "patient": patient,
                    "module": module,
                    "requested_genes": len(requested_genes),
                    "present_genes": len(present),
                    "coverage": len(present) / len(requested_genes),
                    "genes_present": ";".join(present),
                }
            )

        for gene in lock["primary_program"]["genes"]:
            values, _ = score_module(normalized, gene_lookup, [gene])
            scores[f"P4gene_{gene}"] = values

        residual_x = sm.add_constant(
            pd.DataFrame(
                {
                    "macrophage_identity": scores["macrophage_identity"],
                    "hypoxia": scores["hypoxia"],
                    "log1p_umi": np.log1p(scores["nCount_RNA"]),
                },
                index=scores.index,
            ),
            has_constant="add",
        )
        residual_model = sm.OLS(scores["P4_primary"], residual_x).fit()
        scores["P4_residual"] = residual_model.resid
        scores["fibrovascular"] = (
            zscore(scores["fibroblast_ecm"].to_numpy())
            + zscore(scores["endothelial_angiogenic"].to_numpy())
        ) / 2

        coordinates = positions[["array_row", "array_col"]].to_numpy(float)
        # Query the fitted coordinates explicitly so that the first returned
        # index is the spot itself; removing it leaves the six locked neighbours.
        neighbors = (
            NearestNeighbors(n_neighbors=7)
            .fit(coordinates)
            .kneighbors(coordinates, return_distance=False)[:, 1:]
        )
        neighbor_targets = np.column_stack(
            [scores[target].to_numpy()[neighbors].mean(axis=1) for target in target_names]
        )

        residual = scores["P4_residual"].to_numpy()
        high = residual >= np.quantile(residual, 0.75)
        observed = neighbor_targets[high].mean(axis=0) - neighbor_targets[~high].mean(axis=0)

        pathology = scores["Classification"].astype(str).copy()
        counts_by_class = pathology.value_counts()
        pathology = pathology.where(pathology.map(counts_by_class).ge(20), "RARE")
        strata = [np.flatnonzero(pathology.to_numpy() == label) for label in pathology.unique()]
        rng = np.random.default_rng(SEED + patient_index)
        null = np.empty((N_PERMUTATIONS, len(target_names)))
        for permutation in range(N_PERMUTATIONS):
            permuted = residual.copy()
            for indices in strata:
                permuted[indices] = rng.permutation(permuted[indices])
            permuted_high = permuted >= np.quantile(permuted, 0.75)
            null[permutation] = (
                neighbor_targets[permuted_high].mean(axis=0)
                - neighbor_targets[~permuted_high].mean(axis=0)
            )
        null_by_patient[patient] = null

        for target_index, target in enumerate(target_names):
            null_center = np.median(null[:, target_index])
            if target == "fibrovascular":
                p_value = (1 + np.sum(null[:, target_index] >= observed[target_index])) / (
                    N_PERMUTATIONS + 1
                )
                alternative = "greater"
            else:
                p_value = (
                    1
                    + np.sum(
                        np.abs(null[:, target_index] - null_center)
                        >= abs(observed[target_index] - null_center)
                    )
                ) / (N_PERMUTATIONS + 1)
                alternative = "two-sided"
            patient_effect_rows.append(
                {
                    "patient": patient,
                    "subtype": scores["subtype"].iloc[0],
                    "target": target,
                    "spots": scores.shape[0],
                    "high_P4_spots": int(high.sum()),
                    "observed_neighbor_effect": observed[target_index],
                    "null_median": null_center,
                    "permutation_p": p_value,
                    "alternative": alternative,
                    "P4_residual_r_squared": residual_model.rsquared,
                }
            )

        scores = scores.join(positions[["array_row", "array_col", "pixel_row", "pixel_col"]])
        scores["high_P4_residual"] = high
        spot_tables.append(scores.reset_index())

    patient_effects = pd.DataFrame(patient_effect_rows)
    patient_effects["fdr_within_target"] = patient_effects.groupby("target", observed=True)[
        "permutation_p"
    ].transform(lambda values: bh(values))
    patient_effects.to_csv(
        TABLES / "Wu_Visium_spatial_neighbor_effects_by_patient.tsv", sep="\t", index=False
    )

    cohort_rows = []
    for target_index, target in enumerate(target_names):
        observed_values = patient_effects.loc[
            patient_effects["target"].eq(target), "observed_neighbor_effect"
        ].to_numpy()
        cohort_observed = observed_values.mean()
        cohort_null = np.column_stack(
            [null_by_patient[patient][:, target_index] for patient in PATIENTS]
        ).mean(axis=1)
        null_center = np.median(cohort_null)
        if target == "fibrovascular":
            p_value = (1 + np.sum(cohort_null >= cohort_observed)) / (N_PERMUTATIONS + 1)
            alternative = "greater"
        else:
            p_value = (
                1 + np.sum(np.abs(cohort_null - null_center) >= abs(cohort_observed - null_center))
            ) / (N_PERMUTATIONS + 1)
            alternative = "two-sided"
        cohort_rows.append(
            {
                "target": target,
                "patients": len(PATIENTS),
                "mean_patient_effect": cohort_observed,
                "patients_positive": int((observed_values > 0).sum()),
                "cohort_null_median": null_center,
                "cohort_permutation_p": p_value,
                "alternative": alternative,
            }
        )
    cohort_results = pd.DataFrame(cohort_rows)
    secondary_mask = ~cohort_results["target"].eq("fibrovascular")
    cohort_results.loc[secondary_mask, "secondary_fdr"] = bh(
        cohort_results.loc[secondary_mask, "cohort_permutation_p"]
    )
    cohort_results.to_csv(
        TABLES / "Wu_Visium_spatial_neighbor_effects_cohort.tsv", sep="\t", index=False
    )

    spots = pd.concat(spot_tables, ignore_index=True)
    spots.to_csv(
        METADATA / "Wu_Visium_frozen_program_spot_scores.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    pd.DataFrame(coverage_rows).to_csv(
        TABLES / "Wu_Visium_spatial_signature_coverage.tsv", sep="\t", index=False
    )

    representative = spots.loc[spots["patientid"].eq("1160920F")]
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.2))
    for axis, column, title in [
        (axes[0], "P4_residual", "Residualized P4"),
        (axes[1], "fibrovascular", "Fibrovascular target"),
    ]:
        scatter = axis.scatter(
            representative["array_col"], -representative["array_row"],
            c=representative[column], cmap="coolwarm", s=7, linewidths=0,
        )
        axis.set_title(title)
        axis.set_axis_off()
        fig.colorbar(scatter, ax=axis, fraction=0.045, pad=0.02)
    fig.suptitle("Wu Visium patient 1160920F (prespecified representative)")
    fig.tight_layout()
    fig.savefig(FIGURES / "Wu_Visium_P4_fibrovascular_representative.pdf")
    plt.close(fig)

    with (LOGS / "17_wu_visium_spatial_gate2.log").open("w") as handle:
        handle.write(f"patients\t{len(PATIENTS)}\n")
        handle.write(f"spots_analyzed\t{spots.shape[0]}\n")
        handle.write(f"permutations\t{N_PERMUTATIONS}\n")
        handle.write("statistical_unit\tpatient\n")
        handle.write("pathology_stratified_permutations\ttrue\n")

    print("Cohort results")
    print(cohort_results.to_string(index=False))
    print("\nPrimary target by patient")
    print(patient_effects.loc[patient_effects["target"].eq("fibrovascular")].to_string(index=False))


if __name__ == "__main__":
    main()
