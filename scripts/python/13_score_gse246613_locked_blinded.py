#!/usr/bin/env python3
"""Create outcome-blinded frozen-program scores for GSE246613."""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import yaml


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "data" / "raw" / "GSE246613_PembroRT_immune_R100_final.h5ad"
LOCK_FILE = ROOT / "config" / "locked_macrophage_programs_v1.yml"
DISCOVERY = ROOT / "results" / "tables" / "GSE176078_macrophage_nmf_program_genes_pilot.tsv"
TABLES = ROOT / "results" / "tables"
METADATA = ROOT / "metadata"
LOGS = ROOT / "logs"
CHUNK_SIZE = 5_000

MACROPHAGE_LIKE_CLUSTERS = {
    "myeloid_00", "myeloid_01", "myeloid_02", "myeloid_04",
    "myeloid_05", "myeloid_06", "myeloid_10", "myeloid_11",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    with LOCK_FILE.open() as handle:
        lock = yaml.safe_load(handle)
    discovery = pd.read_csv(DISCOVERY, sep="\t")
    modules = {
        "P4_primary": lock["primary_program"]["genes"],
        "P3_secondary": lock["secondary_program"]["genes"],
        "P1_inflammatory": discovery.loc[
            discovery["program"].eq("P1") & discovery["rank"].le(30), "gene"
        ].astype(str).tolist(),
        "P6_interferon": discovery.loc[
            discovery["program"].eq("P6") & discovery["rank"].le(30), "gene"
        ].astype(str).tolist(),
        "hypoxia": [
            "HIF1A", "EGLN1", "CA9", "VEGFA", "LDHA", "PDK1", "SLC2A1",
            "BNIP3", "NDRG1", "ENO1", "ALDOA", "PGK1", "HK2", "PFKP",
        ],
        "macrophage_identity": ["C1QA", "C1QB", "C1QC", "CD68", "MSR1"],
        "phagolysosome": ["FCER1G", "TYROBP", "LST1", "AIF1", "CTSB", "LAMP1"],
    }

    adata = sc.read_h5ad(INPUT, backed="r")
    gene_to_col = {str(gene): index for index, gene in enumerate(adata.var_names)}
    present = {name: [gene for gene in genes if gene in gene_to_col] for name, genes in modules.items()}
    union = sorted({gene for genes in present.values() for gene in genes})
    union_cols = np.array([gene_to_col[gene] for gene in union], dtype=int)
    union_lookup = {gene: index for index, gene in enumerate(union)}

    obs = adata.obs
    selected_mask = (
        obs["celltype"].astype(str).eq("myeloid")
        & obs["subcluster"].astype(str).isin(MACROPHAGE_LIKE_CLUSTERS)
    ).to_numpy()
    score_chunks = []
    for start in range(0, adata.n_obs, CHUNK_SIZE):
        stop = min(start + CHUNK_SIZE, adata.n_obs)
        local = selected_mask[start:stop]
        if not local.any():
            continue
        block = adata.X[start:stop, union_cols][local].tocsr().astype(np.float64)
        cell_obs = obs.iloc[start:stop].loc[local]
        library = cell_obs["n_counts"].astype(float).to_numpy()
        block = block.multiply((10_000 / library)[:, None]).tocsr()
        block.data = np.log1p(block.data)

        scores = cell_obs[
            ["cohort", "patient_treatment", "treatment", "subcluster", "n_counts", "n_genes_by_counts"]
        ].copy()
        scores.index.name = "cell_id"
        for name, genes in present.items():
            columns = [union_lookup[gene] for gene in genes]
            scores[name] = np.asarray(block[:, columns].mean(axis=1)).ravel()
        p4_genes = present["P4_primary"]
        for omitted in p4_genes:
            columns = [union_lookup[gene] for gene in p4_genes if gene != omitted]
            scores[f"P4_loo_{omitted}"] = np.asarray(block[:, columns].mean(axis=1)).ravel()
        score_chunks.append(scores.reset_index())

    cell_scores = pd.concat(score_chunks, ignore_index=True)
    cell_path = METADATA / "GSE246613_locked_macrophage_program_cell_scores_blinded.tsv.gz"
    cell_scores.to_csv(cell_path, sep="\t", index=False, compression="gzip")

    loo_columns = [f"P4_loo_{gene}" for gene in present["P4_primary"]]
    score_columns = list(modules) + loo_columns
    patient_timepoint = (
        cell_scores.groupby(["cohort", "patient_treatment", "treatment"], observed=True)[score_columns]
        .median()
        .reset_index()
    )
    selected_counts = (
        cell_scores.groupby(["cohort", "patient_treatment", "treatment"], observed=True)
        .agg(
            macrophage_like_cells=("cell_id", "size"),
            median_n_counts=("n_counts", "median"),
            median_n_genes=("n_genes_by_counts", "median"),
        )
        .reset_index()
    )
    patient_timepoint = patient_timepoint.merge(
        selected_counts, on=["cohort", "patient_treatment", "treatment"], how="left"
    )

    all_counts = (
        obs.groupby(["cohort", "patient_treatment", "treatment"], observed=True)
        .size()
        .rename("all_immune_cells")
        .reset_index()
    )
    myeloid_counts = (
        obs.loc[obs["celltype"].astype(str).eq("myeloid")]
        .groupby(["cohort", "patient_treatment", "treatment"], observed=True)
        .size()
        .rename("all_myeloid_cells")
        .reset_index()
    )
    patient_timepoint = patient_timepoint.merge(
        all_counts, on=["cohort", "patient_treatment", "treatment"], how="left"
    ).merge(
        myeloid_counts, on=["cohort", "patient_treatment", "treatment"], how="left"
    )
    patient_timepoint["macrophage_like_fraction_of_immune"] = (
        patient_timepoint["macrophage_like_cells"] / patient_timepoint["all_immune_cells"]
    )
    patient_timepoint["macrophage_like_fraction_of_myeloid"] = (
        patient_timepoint["macrophage_like_cells"] / patient_timepoint["all_myeloid_cells"]
    )

    patient_path = TABLES / "GSE246613_locked_program_by_patient_timepoint_blinded.tsv"
    patient_timepoint.to_csv(patient_path, sep="\t", index=False)
    coverage = pd.DataFrame(
        {
            "module": list(modules),
            "requested_genes": [len(modules[name]) for name in modules],
            "present_genes": [len(present[name]) for name in modules],
            "coverage": [len(present[name]) / len(modules[name]) for name in modules],
            "genes_present": [";".join(present[name]) for name in modules],
        }
    )
    coverage.to_csv(
        TABLES / "GSE246613_locked_program_module_coverage.tsv", sep="\t", index=False
    )

    with (LOGS / "13_score_gse246613_locked_blinded.log").open("w") as handle:
        handle.write("outcome_columns_read_or_exported\tfalse\n")
        handle.write(f"lock_sha256\t{lock['source_table_sha256']}\n")
        handle.write(f"selected_cells\t{cell_scores.shape[0]}\n")
        handle.write(f"patient_timepoints\t{patient_timepoint.shape[0]}\n")
        handle.write(f"patient_table_sha256\t{sha256(patient_path)}\n")
        handle.write(f"cell_table_sha256\t{sha256(cell_path)}\n")

    print(coverage.to_string(index=False))
    print(f"selected_cells={cell_scores.shape[0]}; patient_timepoints={patient_timepoint.shape[0]}")
    print(f"patient_table_sha256={sha256(patient_path)}")


if __name__ == "__main__":
    main()
