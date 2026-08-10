#!/usr/bin/env python3
"""Independent, lineage-gated macrophage NMF in GSE114725 tumor samples."""

from __future__ import annotations

from itertools import combinations
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.optimize import linear_sum_assignment
from scipy.stats import hypergeom
from sklearn.decomposition import NMF
from sklearn.metrics.pairwise import cosine_similarity


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw" / "GSE114725_rna_raw.csv.gz"
SCORES = ROOT / "metadata" / "GSE114725_marker_program_scores_raw.tsv.gz"
DISCOVERY = ROOT / "results" / "tables" / "GSE176078_macrophage_nmf_program_genes_pilot.tsv"
TABLES = ROOT / "results" / "tables"
INTERIM = ROOT / "data" / "interim"
LOGS = ROOT / "logs"

SEED = 20260721
MAX_CELLS_PER_PATIENT = 150
LINEAGE_DETECTION_THRESHOLD = 0.4
N_HVG = 500
RANKS = (5, 6, 7)
N_STARTS = 5
TOP_N = 50
META = ["patient", "tissue", "replicate", "cluster", "cellid"]


def matched_component_similarity(a: np.ndarray, b: np.ndarray) -> float:
    similarity = cosine_similarity(a, b)
    rows, cols = linear_sum_assignment(-similarity)
    return float(similarity[rows, cols].mean())


def benjamini_hochberg(p_values: pd.Series) -> np.ndarray:
    values = p_values.to_numpy(dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.minimum(adjusted, 1.0)
    return output


def main() -> None:
    marker_scores = pd.read_csv(SCORES, sep="\t")
    eligible = (
        marker_scores["tissue"].astype(str).eq("TUMOR")
        & marker_scores["macrophage_lineage_detected_fraction"].ge(LINEAGE_DETECTION_THRESHOLD)
        & marker_scores["t_cell_detected_fraction"].lt(0.3)
        & marker_scores["b_cell_detected_fraction"].lt(0.5)
    )
    candidates = marker_scores.loc[eligible].copy()
    sampled = (
        candidates.groupby("patient", observed=True, group_keys=False)
        .apply(
            lambda group: group.sample(
                n=min(group.shape[0], MAX_CELLS_PER_PATIENT), random_state=SEED
            ),
            include_groups=False,
        )
        .sort_index()
    )
    selected_rows = set(sampled.index.astype(int))

    # The marker-score table preserves raw CSV row order. Skipping unselected
    # rows avoids materializing the 2.8-GB dense CSV while retaining all genes.
    raw = pd.read_csv(
        RAW,
        skiprows=lambda row_number: row_number > 0 and (row_number - 1) not in selected_rows,
        low_memory=False,
    )
    gene_columns = [column for column in raw.columns if column not in META]
    counts = sparse.csr_matrix(raw[gene_columns].to_numpy(dtype=np.float32))
    obs = raw[META].copy()
    obs.index = [f"GSE114725_{i}" for i in range(obs.shape[0])]
    var = pd.DataFrame(index=pd.Index(gene_columns, name="gene"))
    adata = ad.AnnData(X=counts, obs=obs, var=var)
    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=10_000)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=N_HVG)
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    x = adata_hvg.X.tocsr().astype(np.float64)

    fits: dict[int, list[NMF]] = {}
    metrics = []
    for rank in RANKS:
        rank_fits = []
        for start in range(N_STARTS):
            model = NMF(
                n_components=rank,
                init="random",
                solver="cd",
                max_iter=1_000,
                tol=1e-4,
                random_state=SEED + rank * 100 + start,
            )
            model.fit(x)
            rank_fits.append(model)
        similarities = [
            matched_component_similarity(first.components_, second.components_)
            for first, second in combinations(rank_fits, 2)
        ]
        errors = [model.reconstruction_err_ for model in rank_fits]
        metrics.append(
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

    metrics = pd.DataFrame(metrics)
    metrics.to_csv(TABLES / "GSE114725_macrophage_independent_nmf_metrics.tsv", sep="\t", index=False)
    selected_rank = int(
        metrics.sort_values(
            ["mean_component_stability", "mean_reconstruction_error"], ascending=[False, True]
        ).iloc[0]["rank"]
    )
    selected_model = min(fits[selected_rank], key=lambda model: model.reconstruction_err_)
    exposures = selected_model.transform(x)
    components = selected_model.components_

    validation_rows = []
    hvg_names = adata_hvg.var_names.to_numpy()
    for index, component in enumerate(components, start=1):
        order = np.argsort(component)[::-1][:TOP_N]
        for gene_rank, position in enumerate(order, start=1):
            validation_rows.append(
                {
                    "validation_program": f"A{index}",
                    "rank": gene_rank,
                    "gene": hvg_names[position],
                    "loading": component[position],
                }
            )
    validation = pd.DataFrame(validation_rows)
    validation.to_csv(
        TABLES / "GSE114725_macrophage_independent_nmf_program_genes.tsv", sep="\t", index=False
    )

    discovery = pd.read_csv(DISCOVERY, sep="\t")
    discovery = discovery.loc[discovery["rank"].le(TOP_N)]
    universe = set(hvg_names.astype(str))
    component_percentiles = {}
    for index, component in enumerate(components, start=1):
        descending = np.argsort(component)[::-1]
        percentile = 1.0 - np.arange(len(descending)) / max(len(descending) - 1, 1)
        component_percentiles[f"A{index}"] = {
            str(hvg_names[position]): float(percentile[rank_position])
            for rank_position, position in enumerate(descending)
        }
    comparisons = []
    for discovery_program, dgroup in discovery.groupby("program", observed=True):
        dgenes = set(dgroup["gene"].astype(str))
        comparable_discovery = dgenes & universe
        for validation_program, vgroup in validation.groupby("validation_program", observed=True):
            vgenes = set(vgroup["gene"].astype(str))
            overlap = sorted(comparable_discovery & vgenes)
            hypergeom_p = hypergeom.sf(
                len(overlap) - 1,
                len(universe),
                len(comparable_discovery),
                len(vgenes),
            )
            comparisons.append(
                {
                    "discovery_program": discovery_program,
                    "validation_program": validation_program,
                    "overlap_n": len(overlap),
                    "jaccard": len(overlap) / len(dgenes | vgenes),
                    "overlap_coefficient": len(overlap) / min(len(dgenes), len(vgenes)),
                    "universe_genes": len(universe),
                    "discovery_genes_in_universe": len(comparable_discovery),
                    "hypergeom_p": hypergeom_p,
                    "mean_loading_percentile": np.mean(
                        [component_percentiles[validation_program][gene] for gene in comparable_discovery]
                    ),
                    "overlapping_genes": ";".join(overlap),
                }
            )
    comparison = pd.DataFrame(comparisons)
    comparison["hypergeom_fdr"] = benjamini_hochberg(comparison["hypergeom_p"])
    comparison.to_csv(
        TABLES / "GSE114725_vs_GSE176078_nmf_program_overlap.tsv", sep="\t", index=False
    )
    best = (
        comparison.sort_values(
            ["discovery_program", "hypergeom_fdr", "mean_loading_percentile"],
            ascending=[True, True, False],
        )
        .groupby("discovery_program", observed=True)
        .head(1)
    )
    best.to_csv(
        TABLES / "GSE114725_vs_GSE176078_nmf_best_matches.tsv", sep="\t", index=False
    )

    exposure_frame = adata.obs[["patient", "tissue"]].copy()
    for index in range(selected_rank):
        exposure_frame[f"A{index + 1}_exposure"] = exposures[:, index]
    patient_exposures = exposure_frame.groupby("patient", observed=True).median(numeric_only=True).reset_index()
    patient_exposures.to_csv(
        TABLES / "GSE114725_macrophage_independent_nmf_patient_exposures.tsv", sep="\t", index=False
    )

    adata.write_h5ad(INTERIM / "GSE114725_balanced_tumor_macrophages.h5ad", compression="gzip")
    with (LOGS / "10_independent_nmf_gse114725.log").open("w") as handle:
        handle.write("response_labels_used\tfalse\n")
        handle.write("cell_selection_uses_discovery_genes\tfalse\n")
        handle.write(f"lineage_detection_threshold\t{LINEAGE_DETECTION_THRESHOLD}\n")
        handle.write(f"patients\t{adata.n_obs and adata.obs['patient'].nunique()}\n")
        handle.write(f"cells\t{adata.n_obs}\n")
        handle.write(f"selected_rank\t{selected_rank}\n")

    print(metrics.to_string(index=False))
    print(f"selected_rank={selected_rank}; cells={adata.n_obs}; patients={adata.obs['patient'].nunique()}")
    print(best.to_string(index=False))


if __name__ == "__main__":
    main()
