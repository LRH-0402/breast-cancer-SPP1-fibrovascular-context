#!/usr/bin/env python3
"""Audit GSE114725 and profile immune clusters using a marker-only read."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "data" / "raw" / "GSE114725_rna_raw.csv.gz"
PROGRAMS = ROOT / "results" / "tables" / "GSE176078_macrophage_nmf_program_genes_pilot.tsv"
TABLES = ROOT / "results" / "tables"
METADATA = ROOT / "metadata"

META = ["patient", "tissue", "replicate", "cluster", "cellid"]
MARKERS = {
    # Lineage gate excludes discovery-program genes to avoid circularly
    # selecting cells for the independent program analysis.
    "macrophage_lineage": ["C1QA", "C1QB", "C1QC", "CD68", "MSR1"],
    "macrophage_core": ["C1QA", "C1QB", "C1QC", "APOE", "CD68", "CTSD", "LGMN", "MSR1"],
    "lipid_lysosomal": ["APOC1", "GPNMB", "LIPA", "LPL", "TREM2", "FABP5", "PLIN2", "CD36"],
    "spp1_matrix": ["SPP1", "FN1", "MMP9", "MARCO", "VEGFA", "HIF1A"],
    "monocyte": ["S100A8", "S100A9", "FCN1", "VCAN", "CTSS", "LYZ"],
    "dendritic": ["FCER1A", "CD1C", "CLEC10A", "CLEC9A", "BATF3", "XCR1"],
    "t_cell": ["CD3D", "CD3E", "TRAC", "IL7R", "CD8A"],
    "b_cell": ["CD79A", "MS4A1", "CD37", "CD74", "CD22"],
}


def main() -> None:
    discovery = pd.read_csv(PROGRAMS, sep="\t")
    discovery = discovery.loc[discovery["rank"].le(30)]
    program_sets = {
        f"discovery_{program}": group.sort_values("rank")["gene"].astype(str).tolist()
        for program, group in discovery.groupby("program", observed=True)
    }
    modules = {**MARKERS, **program_sets}
    requested_genes = sorted({gene for genes in modules.values() for gene in genes})

    header = pd.read_csv(INPUT, nrows=0)
    available = set(header.columns)
    present = {name: [gene for gene in genes if gene in available] for name, genes in modules.items()}
    usecols = META + sorted({gene for genes in present.values() for gene in genes})
    frame = pd.read_csv(INPUT, usecols=usecols, low_memory=False)

    score_frame = frame[META].copy()
    for name, genes in present.items():
        values = np.log1p(frame[genes].to_numpy(dtype=np.float32))
        score_frame[f"{name}_mean_log1p_count"] = values.mean(axis=1)
        score_frame[f"{name}_detected_fraction"] = (values > 0).mean(axis=1)

    counts = (
        score_frame.groupby(["patient", "tissue", "cluster"], observed=True)
        .size()
        .rename("cells")
        .reset_index()
    )
    counts.to_csv(TABLES / "GSE114725_cells_by_patient_tissue_cluster.tsv", sep="\t", index=False)

    score_cols = [column for column in score_frame if column not in META]
    cluster_profiles = (
        score_frame.groupby("cluster", observed=True)[score_cols].mean().reset_index()
    )
    cluster_profiles = cluster_profiles.merge(
        score_frame.groupby("cluster", observed=True).size().rename("cells").reset_index(),
        on="cluster",
        how="left",
    )
    cluster_profiles.to_csv(TABLES / "GSE114725_cluster_marker_profiles.tsv", sep="\t", index=False)

    patient_tissue = (
        score_frame.groupby(["patient", "tissue"], observed=True).size().rename("cells").reset_index()
    )
    patient_tissue.to_csv(METADATA / "GSE114725_patient_tissue_counts.tsv", sep="\t", index=False)

    coverage = pd.DataFrame(
        {
            "module": list(modules),
            "requested_genes": [len(modules[name]) for name in modules],
            "present_genes": [len(present[name]) for name in modules],
            "genes_present": [";".join(present[name]) for name in modules],
        }
    )
    coverage.to_csv(TABLES / "GSE114725_marker_program_coverage.tsv", sep="\t", index=False)

    selected_cell_metadata = score_frame[META + score_cols].copy()
    selected_cell_metadata.to_csv(
        METADATA / "GSE114725_marker_program_scores_raw.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    summary = pd.DataFrame(
        {
            "metric": ["cells", "patients", "clusters", "tumor_patients", "normal_patients"],
            "value": [
                frame.shape[0],
                frame["patient"].nunique(),
                frame["cluster"].nunique(),
                frame.loc[frame["tissue"].eq("TUMOR"), "patient"].nunique(),
                frame.loc[frame["tissue"].eq("NORMAL"), "patient"].nunique(),
            ],
        }
    )
    summary.to_csv(TABLES / "GSE114725_audit_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))
    print(cluster_profiles.sort_values("macrophage_core_mean_log1p_count", ascending=False).to_string(index=False))


if __name__ == "__main__":
    main()
