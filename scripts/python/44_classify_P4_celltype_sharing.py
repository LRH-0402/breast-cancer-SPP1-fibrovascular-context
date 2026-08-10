#!/usr/bin/env python3
"""Classify locked P4 genes by empirically observed cell-compartment sharing."""

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "results" / "tables" / "GSE176078_P4_gene_expression_by_celltype.tsv"
OUTPUT = ROOT / "results" / "tables" / "GSE176078_P4_gene_celltype_classification.tsv"


def classify(row):
    dominant = row["dominant_celltype"]
    ratio = row["macrophage_to_strongest_other_ratio"]
    if dominant == "Macrophage" and ratio >= 1.5:
        return "macrophage-dominant"
    if dominant in {"Macrophage", "Other myeloid"}:
        return "myeloid-enriched/shared"
    if dominant == "Lymphocyte":
        return "lymphocyte-weighted"
    return "multicellular/stromal"


def main():
    data = pd.read_csv(INPUT, sep="\t")
    expression = data.pivot(
        index="gene", columns="analysis_celltype", values="mean_expression"
    )
    detection = data.pivot(
        index="gene", columns="analysis_celltype", values="fraction_detected"
    )
    rows = []
    for gene, values in expression.iterrows():
        dominant = values.idxmax()
        strongest_other = values.drop("Macrophage").max()
        rows.append({
            "gene": gene,
            "dominant_celltype": dominant,
            "macrophage_mean_expression": values["Macrophage"],
            "strongest_other_mean_expression": strongest_other,
            "macrophage_to_strongest_other_ratio": (
                values["Macrophage"] / strongest_other if strongest_other > 0 else float("inf")
            ),
            "macrophage_fraction_detected": detection.loc[gene, "Macrophage"],
        })
    output = pd.DataFrame(rows)
    output["empirical_expression_class"] = output.apply(classify, axis=1)
    output = output.sort_values(["empirical_expression_class", "gene"])
    output.to_csv(OUTPUT, sep="\t", index=False)
    print(output.to_string(index=False))


if __name__ == "__main__":
    main()
