#!/usr/bin/env python3
"""Finalize background-corrected cross-cohort program matching and core genes."""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import hypergeom


ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "results" / "tables"
DISCOVERY_FILE = TABLES / "GSE176078_macrophage_nmf_program_genes_pilot.tsv"


def bh(p_values: pd.Series) -> np.ndarray:
    values = p_values.to_numpy(float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.minimum(adjusted, 1.0)
    return output


def match_programs(
    universe: set[str], validation: pd.DataFrame, discovery: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for dp, dgroup in discovery.groupby("program", observed=True):
        dgenes = set(dgroup["gene"].astype(str)) & universe
        for vp, vgroup in validation.groupby("validation_program", observed=True):
            vgenes = set(vgroup["gene"].astype(str))
            overlap = sorted(dgenes & vgenes)
            rows.append(
                {
                    "discovery_program": dp,
                    "validation_program": vp,
                    "overlap_n": len(overlap),
                    "universe_genes": len(universe),
                    "discovery_genes_in_universe": len(dgenes),
                    "hypergeom_p": hypergeom.sf(
                        len(overlap) - 1, len(universe), len(dgenes), len(vgenes)
                    ),
                    "overlapping_genes": ";".join(overlap),
                }
            )
    comparison = pd.DataFrame(rows)
    comparison["hypergeom_fdr"] = bh(comparison["hypergeom_p"])
    best = (
        comparison.sort_values(
            ["discovery_program", "hypergeom_fdr", "overlap_n"],
            ascending=[True, True, False],
        )
        .groupby("discovery_program", observed=True)
        .head(1)
    )
    return comparison, best


def main() -> None:
    discovery = pd.read_csv(DISCOVERY_FILE, sep="\t").query("rank <= 50")

    h5ad = ad.read_h5ad(
        ROOT / "data" / "interim" / "GSE246613_baseline_macrophage_like_independent_nmf_balanced_hvg.h5ad"
    )
    gse246613_validation = pd.read_csv(
        TABLES / "GSE246613_baseline_macrophage_like_independent_nmf_program_genes.tsv", sep="\t"
    )
    gse246613_comparison, gse246613_best = match_programs(
        set(h5ad.var_names.astype(str)), gse246613_validation, discovery
    )
    gse246613_comparison.to_csv(
        TABLES / "GSE246613_formal_program_matching.tsv", sep="\t", index=False
    )
    gse246613_best.to_csv(
        TABLES / "GSE246613_formal_program_best_matches.tsv", sep="\t", index=False
    )

    best_tables = {
        "GSE161529": pd.read_csv(TABLES / "GSE161529_formal_program_best_matches.tsv", sep="\t"),
        "GSE246613": gse246613_best,
        "GSE114725": pd.read_csv(TABLES / "GSE114725_vs_GSE176078_nmf_best_matches.tsv", sep="\t"),
    }
    summary = pd.concat(
        [table.assign(cohort=cohort) for cohort, table in best_tables.items()],
        ignore_index=True,
    )
    retained = [
        "cohort", "discovery_program", "validation_program", "overlap_n",
        "universe_genes", "discovery_genes_in_universe", "hypergeom_p",
        "hypergeom_fdr", "overlapping_genes",
    ]
    summary[retained].to_csv(
        TABLES / "gate1_cross_cohort_program_matching_summary.tsv", sep="\t", index=False
    )

    validation_program_files = {
        "GSE161529": TABLES / "GSE161529_macrophage_independent_nmf_program_genes.tsv",
        "GSE246613": TABLES / "GSE246613_baseline_macrophage_like_independent_nmf_program_genes.tsv",
        "GSE114725": TABLES / "GSE114725_macrophage_independent_nmf_program_genes.tsv",
    }
    core_rows = []
    for program in ["P3", "P4"]:
        discovery_genes = set(discovery.loc[discovery["program"].eq(program), "gene"].astype(str))
        gene_sets = [discovery_genes]
        matched = summary.loc[summary["discovery_program"].eq(program)]
        for row in matched.itertuples(index=False):
            validation = pd.read_csv(validation_program_files[row.cohort], sep="\t")
            genes = set(
                validation.loc[
                    validation["validation_program"].eq(row.validation_program), "gene"
                ].astype(str)
            )
            gene_sets.append(genes)
        support = Counter(gene for genes in gene_sets for gene in genes)
        for gene in sorted(discovery_genes):
            if support[gene] >= 3:
                core_rows.append(
                    {
                        "program": program,
                        "gene": gene,
                        "datasets_supporting_top50_membership": support[gene],
                        "locking_rule": "discovery_top50_and_present_in_at_least_2_of_3_independent_top50_factors",
                    }
                )
    core = pd.DataFrame(core_rows)
    core.to_csv(TABLES / "locked_macrophage_program_core_genes_v1.tsv", sep="\t", index=False)

    print(summary.loc[summary["discovery_program"].isin(["P3", "P4"]), retained].to_string(index=False))
    print("\nLocked core genes")
    print(core.groupby("program", observed=True)["gene"].apply(lambda x: ";".join(x)).to_string())


if __name__ == "__main__":
    main()
