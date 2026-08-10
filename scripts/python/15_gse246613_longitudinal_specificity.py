#!/usr/bin/env python3
"""Exploratory specificity audit of longitudinal GSE246613 module changes."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import wilcoxon


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "results" / "tables" / "GSE246613_locked_program_clinical_join.tsv"
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"

MODULES = [
    "P4_primary", "P3_secondary", "P1_inflammatory", "P6_interferon",
    "hypoxia", "macrophage_identity", "phagolysosome",
]
CONTRASTS = [("Base", "PD1"), ("PD1", "RTPD1"), ("Base", "RTPD1")]


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
    data = pd.read_csv(INPUT, sep="\t")
    rows = []
    for module in MODULES:
        wide = data.pivot(index="cohort", columns="treatment", values=module)
        reference_sd = data.loc[data["treatment"].eq("Base"), module].std(ddof=1)
        for start, end in CONTRASTS:
            paired = wide[[start, end]].dropna()
            difference = paired[end] - paired[start]
            test = wilcoxon(paired[end], paired[start], alternative="two-sided", method="auto")
            rows.append(
                {
                    "module": module,
                    "contrast": f"{end}_minus_{start}",
                    "paired_patients": paired.shape[0],
                    "patients_increased": int((difference > 0).sum()),
                    "median_change": difference.median(),
                    "median_change_in_baseline_sd": difference.median() / reference_sd,
                    "wilcoxon_statistic": test.statistic,
                    "p_value": test.pvalue,
                    "analysis_status": "exploratory specificity audit",
                }
            )
    results = pd.DataFrame(rows)
    results["fdr_across_modules_and_contrasts"] = bh(results["p_value"])
    results.to_csv(
        TABLES / "GSE246613_longitudinal_module_specificity.tsv", sep="\t", index=False
    )

    heat = results.pivot(
        index="module", columns="contrast", values="median_change_in_baseline_sd"
    ).reindex(index=MODULES, columns=[f"{end}_minus_{start}" for start, end in CONTRASTS])
    annotations = results.pivot(
        index="module", columns="contrast", values="fdr_across_modules_and_contrasts"
    ).reindex(index=heat.index, columns=heat.columns)
    labels = annotations.map(lambda value: "***" if value < 0.001 else "**" if value < 0.01 else "*" if value < 0.05 else "")

    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    sns.heatmap(
        heat,
        annot=labels,
        fmt="",
        cmap="vlag",
        center=0,
        linewidths=0.5,
        cbar_kws={"label": "Median change / baseline SD"},
        ax=ax,
    )
    ax.set_xlabel(None)
    ax.set_ylabel(None)
    ax.set_title("Longitudinal module-change specificity audit")
    fig.tight_layout()
    fig.savefig(FIGURES / "GSE246613_longitudinal_module_specificity.pdf")
    plt.close(fig)
    print(results.sort_values("fdr_across_modules_and_contrasts").to_string(index=False))


if __name__ == "__main__":
    main()
