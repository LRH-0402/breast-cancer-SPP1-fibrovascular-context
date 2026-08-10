#!/usr/bin/env python3
"""Profile lineage marker modules across GSE246613 author-defined myeloid clusters."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "data" / "raw" / "GSE246613_PembroRT_immune_R100_final.h5ad"
TABLES = ROOT / "results" / "tables"
CHUNK_SIZE = 5_000

MARKERS = {
    "macrophage_core": ["C1QA", "C1QB", "C1QC", "APOE", "CD68", "CTSD", "LGMN", "MSR1"],
    "lipid_lysosomal": ["APOC1", "GPNMB", "LIPA", "LPL", "TREM2", "FABP5", "PLIN2", "CD36"],
    "spp1_matrix": ["SPP1", "FN1", "MMP9", "MARCO", "VEGFA", "HIF1A"],
    "classical_monocyte": ["S100A8", "S100A9", "S100A10", "FCN1", "VCAN", "CTSS", "LYZ"],
    "dendritic": ["FCER1A", "CD1C", "CLEC10A", "CLEC9A", "BATF3", "XCR1"],
    "neutrophil": ["FCGR3B", "CSF3R", "CXCR2", "MNDA", "FPR1", "SELL"],
    "antigen_presentation": ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "CIITA"],
}


def main() -> None:
    adata = sc.read_h5ad(INPUT, backed="r")
    gene_to_col = {str(gene): i for i, gene in enumerate(adata.var_names)}
    present = {name: [gene for gene in genes if gene in gene_to_col] for name, genes in MARKERS.items()}
    union = sorted({gene for genes in present.values() for gene in genes})
    union_cols = np.array([gene_to_col[gene] for gene in union], dtype=int)
    lookup = {gene: i for i, gene in enumerate(union)}
    myeloid = adata.obs["celltype"].astype(str).eq("myeloid").to_numpy()

    chunks = []
    for start in range(0, adata.n_obs, CHUNK_SIZE):
        stop = min(start + CHUNK_SIZE, adata.n_obs)
        local = myeloid[start:stop]
        if not local.any():
            continue
        block = adata.X[start:stop, union_cols][local].tocsr().astype(np.float64)
        cell_obs = adata.obs.iloc[start:stop].loc[local]
        lib = cell_obs["n_counts"].astype(float).to_numpy()
        block = block.multiply((10_000 / lib)[:, None]).tocsr()
        block.data = np.log1p(block.data)

        scores = pd.DataFrame({"subcluster": cell_obs["subcluster"].astype(str).to_numpy()})
        for name, genes in present.items():
            indices = [lookup[gene] for gene in genes]
            scores[name] = np.asarray(block[:, indices].mean(axis=1)).ravel()
        chunks.append(scores)

    cell_scores = pd.concat(chunks, ignore_index=True)
    summary = cell_scores.groupby("subcluster", observed=True).agg(["median", "mean"])
    summary.columns = [f"{module}_{stat}" for module, stat in summary.columns]
    summary = summary.reset_index()
    counts = cell_scores.groupby("subcluster", observed=True).size().rename("cells").reset_index()
    summary = summary.merge(counts, on="subcluster", how="left")
    summary.to_csv(TABLES / "GSE246613_myeloid_cluster_marker_modules.tsv", sep="\t", index=False)

    coverage = pd.DataFrame(
        {
            "module": list(MARKERS),
            "requested_genes": [len(MARKERS[name]) for name in MARKERS],
            "present_genes": [len(present[name]) for name in MARKERS],
            "genes_present": [";".join(present[name]) for name in MARKERS],
        }
    )
    coverage.to_csv(TABLES / "GSE246613_myeloid_marker_module_coverage.tsv", sep="\t", index=False)
    print(summary.sort_values("macrophage_core_median", ascending=False).to_string(index=False))


if __name__ == "__main__":
    main()
