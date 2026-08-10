#!/usr/bin/env python3
"""Outcome-blinded independent NMF of baseline GSE246613 myeloid cells.

Use ``--scope broad`` as a negative-control decomposition or the default
``--scope macrophage_like`` after marker-based lineage restriction.
"""

from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.optimize import linear_sum_assignment
from sklearn.decomposition import NMF


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "data" / "raw" / "GSE246613_PembroRT_immune_R100_final.h5ad"
DISCOVERY = ROOT / "results" / "tables" / "GSE176078_macrophage_nmf_program_genes_pilot.tsv"
TABLES = ROOT / "results" / "tables"
INTERIM = ROOT / "data" / "interim"
LOGS = ROOT / "logs"
SEED = 20260721
MAX_CELLS_PER_PATIENT = 100
N_HVG = 500
RANKS = (5, 6, 7)
N_STARTS = 5
TOP_N = 50
MACROPHAGE_LIKE_CLUSTERS = {
    "myeloid_00", "myeloid_01", "myeloid_02", "myeloid_04",
    "myeloid_05", "myeloid_06", "myeloid_10", "myeloid_11",
}


def component_similarity(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    a_norm = a / np.maximum(np.linalg.norm(a, axis=1, keepdims=True), 1e-12)
    b_norm = b / np.maximum(np.linalg.norm(b, axis=1, keepdims=True), 1e-12)
    similarity = np.nan_to_num(a_norm @ b_norm.T, nan=0.0, posinf=1.0, neginf=0.0)
    rows, cols = linear_sum_assignment(-similarity)
    return float(similarity[rows, cols].mean())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scope", choices=["broad", "macrophage_like"], default="macrophage_like")
    args = parser.parse_args()
    prefix = f"GSE246613_baseline_{args.scope}_independent_nmf"

    rng = np.random.default_rng(SEED)
    backed = sc.read_h5ad(INPUT, backed="r")
    obs = backed.obs
    eligible = obs["celltype"].astype(str).eq("myeloid") & obs["treatment"].astype(str).eq("Base")
    if args.scope == "macrophage_like":
        eligible &= obs["subcluster"].astype(str).isin(MACROPHAGE_LIKE_CLUSTERS)
    eligible_positions = np.flatnonzero(eligible.to_numpy())

    chosen: list[int] = []
    eligible_obs = obs.iloc[eligible_positions]
    for _, group in eligible_obs.groupby("cohort", observed=True):
        positions = obs.index.get_indexer(group.index)
        size = min(MAX_CELLS_PER_PATIENT, len(positions))
        chosen.extend(rng.choice(positions, size=size, replace=False).tolist())
    chosen = sorted(chosen)

    adata = backed[chosen, :].to_memory()
    adata.obs = adata.obs.drop(columns=["pCR", "RCB"], errors="ignore")
    sc.pp.normalize_total(adata, target_sum=10_000)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=N_HVG)
    adata = adata[:, adata.var["highly_variable"]].copy()
    x = adata.X.tocsr().astype(np.float64)

    fits: dict[int, list[tuple[NMF, np.ndarray]]] = {}
    metric_rows = []
    for rank in RANKS:
        rank_fits: list[tuple[NMF, np.ndarray]] = []
        for start in range(N_STARTS):
            model = NMF(
                n_components=rank,
                init="random",
                solver="cd",
                max_iter=1_000,
                random_state=SEED + 100 * rank + start,
                tol=1e-4,
            )
            model.fit_transform(x)
            rank_fits.append((model, model.components_.copy()))
        similarities = [
            component_similarity(first[1], second[1])
            for first, second in combinations(rank_fits, 2)
        ]
        errors = [fit[0].reconstruction_err_ for fit in rank_fits]
        metric_rows.append(
            {
                "rank": rank,
                "starts": N_STARTS,
                "mean_component_stability": np.mean(similarities),
                "min_component_stability": np.min(similarities),
                "mean_reconstruction_error": np.mean(errors),
                "min_reconstruction_error": np.min(errors),
            }
        )
        fits[rank] = rank_fits

    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(TABLES / f"{prefix}_metrics.tsv", sep="\t", index=False)
    selected_rank = int(
        metrics.sort_values(
            ["mean_component_stability", "mean_reconstruction_error"],
            ascending=[False, True],
        ).iloc[0]["rank"]
    )
    selected_model, components = min(
        fits[selected_rank], key=lambda value: value[0].reconstruction_err_
    )

    program_rows = []
    genes = adata.var_names.to_numpy()
    for component_index, loadings in enumerate(components, start=1):
        order = np.argsort(loadings)[::-1][:TOP_N]
        for gene_rank, position in enumerate(order, start=1):
            program_rows.append(
                {
                    "validation_program": f"B{component_index}",
                    "rank": gene_rank,
                    "gene": genes[position],
                    "loading": loadings[position],
                }
            )
    validation = pd.DataFrame(program_rows)
    validation.to_csv(
        TABLES / f"{prefix}_program_genes.tsv",
        sep="\t",
        index=False,
    )

    discovery = pd.read_csv(DISCOVERY, sep="\t")
    discovery = discovery.loc[discovery["rank"].le(TOP_N)]
    comparisons = []
    for discovery_program, dgroup in discovery.groupby("program", observed=True):
        dgenes = set(dgroup["gene"].astype(str))
        for validation_program, vgroup in validation.groupby("validation_program", observed=True):
            vgenes = set(vgroup["gene"].astype(str))
            overlap = sorted(dgenes & vgenes)
            comparisons.append(
                {
                    "discovery_program": discovery_program,
                    "validation_program": validation_program,
                    "overlap_n": len(overlap),
                    "jaccard": len(overlap) / len(dgenes | vgenes),
                    "overlap_coefficient": len(overlap) / min(len(dgenes), len(vgenes)),
                    "overlapping_genes": ";".join(overlap),
                }
            )
    comparison = pd.DataFrame(comparisons)
    comparison.to_csv(
        TABLES / f"{prefix}_vs_GSE176078_program_overlap.tsv", sep="\t", index=False
    )
    best = (
        comparison.sort_values(
            ["discovery_program", "overlap_coefficient"], ascending=[True, False]
        )
        .groupby("discovery_program", observed=True)
        .head(1)
    )
    best.to_csv(
        TABLES / f"{prefix}_vs_GSE176078_best_matches.tsv", sep="\t", index=False
    )

    adata.write_h5ad(INTERIM / f"{prefix}_balanced_hvg.h5ad", compression="gzip")
    with (LOGS / f"07_{prefix}.log").open("w") as handle:
        handle.write(f"response_labels_used\tfalse\n")
        handle.write(f"scope\t{args.scope}\n")
        if args.scope == "macrophage_like":
            handle.write(f"included_subclusters\t{';'.join(sorted(MACROPHAGE_LIKE_CLUSTERS))}\n")
        handle.write(f"patients\t{adata.obs['cohort'].nunique()}\n")
        handle.write(f"cells\t{adata.n_obs}\n")
        handle.write(f"highly_variable_genes\t{adata.n_vars}\n")
        handle.write(f"selected_rank_by_stability\t{selected_rank}\n")

    print(metrics.to_string(index=False))
    print(f"selected_rank={selected_rank}")
    print(best.to_string(index=False))


if __name__ == "__main__":
    main()
