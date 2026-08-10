#!/usr/bin/env python3
"""Project frozen discovery programs into GSE246613 myeloid cells.

Clinical response labels are intentionally not used in this script. Scores are
computed from log1p(CP10K) expression and summarized at the patient-timepoint
level so that the patient, rather than the cell, remains the inferential unit.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "data" / "raw" / "GSE246613_PembroRT_immune_R100_final.h5ad"
PROGRAMS = ROOT / "results" / "tables" / "GSE176078_macrophage_nmf_program_genes_pilot.tsv"
TABLES = ROOT / "results" / "tables"
METADATA = ROOT / "metadata"
LOGS = ROOT / "logs"
TOP_N = 30
CHUNK_SIZE = 5_000


def main() -> None:
    genes = pd.read_csv(PROGRAMS, sep="\t")
    genes = genes.loc[genes["rank"].le(TOP_N)].copy()
    program_order = sorted(genes["program"].unique())

    adata = sc.read_h5ad(INPUT, backed="r")
    var_names = pd.Index(adata.var_names.astype(str))
    gene_to_col = {gene: i for i, gene in enumerate(var_names)}

    program_genes = {
        program: group.sort_values("rank")["gene"].astype(str).tolist()
        for program, group in genes.groupby("program", observed=True)
    }
    present = {
        program: [gene for gene in values if gene in gene_to_col]
        for program, values in program_genes.items()
    }
    union_genes = sorted({gene for values in present.values() for gene in values})
    union_cols = np.array([gene_to_col[gene] for gene in union_genes], dtype=int)
    union_lookup = {gene: i for i, gene in enumerate(union_genes)}

    coverage_rows = []
    for program in program_order:
        coverage_rows.append(
            {
                "program": program,
                "requested_genes": len(program_genes[program]),
                "present_genes": len(present[program]),
                "coverage": len(present[program]) / len(program_genes[program]),
                "genes_present": ";".join(present[program]),
            }
        )
    pd.DataFrame(coverage_rows).to_csv(
        TABLES / "GSE246613_program_projection_gene_coverage.tsv", sep="\t", index=False
    )

    obs = adata.obs
    myeloid_mask = obs["celltype"].astype(str).eq("myeloid").to_numpy()
    output_chunks: list[pd.DataFrame] = []

    for start in range(0, adata.n_obs, CHUNK_SIZE):
        stop = min(start + CHUNK_SIZE, adata.n_obs)
        local_mask = myeloid_mask[start:stop]
        if not local_mask.any():
            continue

        block = adata.X[start:stop, union_cols]
        block = block[local_mask].tocsr().astype(np.float64)
        cell_obs = obs.iloc[start:stop].loc[local_mask]
        library_size = cell_obs["n_counts"].astype(float).to_numpy()
        if np.any(library_size <= 0):
            raise ValueError("Non-positive n_counts encountered in myeloid cells")

        block = block.multiply((10_000.0 / library_size)[:, None]).tocsr()
        block.data = np.log1p(block.data)

        scored = cell_obs[
            ["cohort", "patient_treatment", "treatment", "subcluster", "batch"]
        ].copy()
        scored.index.name = "cell_id"
        for program in program_order:
            local_cols = [union_lookup[gene] for gene in present[program]]
            if not local_cols:
                raise ValueError(f"No genes available for {program}")
            scored[f"{program}_projection"] = np.asarray(
                block[:, local_cols].mean(axis=1)
            ).ravel()
        output_chunks.append(scored.reset_index())

    scores = pd.concat(output_chunks, ignore_index=True)
    scores.to_csv(
        METADATA / "GSE246613_myeloid_program_projection.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    score_cols = [f"{program}_projection" for program in program_order]
    patient_timepoint = (
        scores.groupby(["cohort", "patient_treatment", "treatment"], observed=True)[score_cols]
        .median()
        .reset_index()
    )
    cell_counts = (
        scores.groupby(["cohort", "patient_treatment", "treatment"], observed=True)
        .size()
        .rename("myeloid_cells")
        .reset_index()
    )
    patient_timepoint = patient_timepoint.merge(
        cell_counts, on=["cohort", "patient_treatment", "treatment"], how="left"
    )
    patient_timepoint.to_csv(
        TABLES / "GSE246613_program_projection_by_patient_timepoint_blinded.tsv",
        sep="\t",
        index=False,
    )

    subclusters = (
        scores.groupby(["subcluster"], observed=True)[score_cols].median().reset_index()
    )
    subcluster_counts = scores.groupby("subcluster", observed=True).size().rename("cells").reset_index()
    subclusters = subclusters.merge(subcluster_counts, on="subcluster", how="left")
    subclusters.to_csv(
        TABLES / "GSE246613_program_projection_by_myeloid_subcluster.tsv",
        sep="\t",
        index=False,
    )

    with (LOGS / "05_project_programs_gse246613.log").open("w") as handle:
        handle.write(f"myeloid_cells\t{scores.shape[0]}\n")
        handle.write(f"patient_timepoints\t{patient_timepoint.shape[0]}\n")
        handle.write(f"top_genes_per_program\t{TOP_N}\n")
        handle.write("response_labels_used\tfalse\n")
        for row in coverage_rows:
            handle.write(
                f"{row['program']}_gene_coverage\t{row['present_genes']}/{row['requested_genes']}\n"
            )

    print(patient_timepoint.groupby("treatment", observed=True).size().to_string())
    print(pd.DataFrame(coverage_rows)[["program", "present_genes", "coverage"]].to_string(index=False))


if __name__ == "__main__":
    main()
