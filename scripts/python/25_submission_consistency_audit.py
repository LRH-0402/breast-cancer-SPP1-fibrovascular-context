#!/usr/bin/env python3
"""Fail fast when headline manuscript claims drift from frozen result tables."""

from pathlib import Path
import re
import zipfile

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "results" / "tables"
MANUSCRIPT = ROOT / "manuscript" / "manuscript_npj_breast_cancer.md"


def close(actual: float, expected: float, tolerance: float = 1e-10) -> None:
    assert abs(float(actual) - expected) <= tolerance, (actual, expected)


def row(path: str, column: str, value: str) -> pd.Series:
    table = pd.read_csv(TABLES / path, sep="\t")
    selected = table.loc[table[column].astype(str) == value]
    assert len(selected) == 1, (path, column, value, len(selected))
    return selected.iloc[0]


def main() -> None:
    gate1 = pd.read_csv(
        TABLES / "gate1_cross_cohort_program_matching_summary.tsv", sep="\t"
    )
    p4 = gate1.loc[gate1["discovery_program"] == "P4"].set_index("cohort")
    expected_fdr = {
        "GSE161529": 1.5912831306686003e-10,
        "GSE246613": 6.140483952936966e-11,
        "GSE114725": 2.383347127407293e-05,
    }
    for cohort, expected in expected_fdr.items():
        close(p4.loc[cohort, "hypergeom_fdr"], expected)

    spatial_specs = [
        ("Wu_Visium_spatial_neighbor_effects_cohort.tsv", "patients", 6, "patients_positive", 5),
        ("Figshare_spatial_neighbor_effects_cohort.tsv", "sections", 8, "sections_positive", 7),
        ("TNBC94_spatial_neighbor_effects_cohort.tsv", "patients", 94, "patients_positive_excess", 82),
    ]
    for filename, n_column, n, positive_column, positive in spatial_specs:
        spatial = row(filename, "target", "fibrovascular")
        assert int(spatial[n_column]) == n
        assert int(spatial[positive_column]) == positive
        close(spatial["cohort_permutation_p"], 0.0004997501249375312)

    tnbc94 = row("TNBC94_spatial_neighbor_effects_cohort.tsv", "target", "fibrovascular")
    close(tnbc94["sign_test_p"], 5.6223687374767304e-14)

    response = row("GSE246613_P4_primary_pCR_analysis.tsv", "predictor", "z_P4_primary")
    assert int(response["n"]) == 34
    close(response["odds_ratio_per_sd"], 1.639585142836826)
    close(response["likelihood_ratio_p"], 0.2073729660565187)

    survival = pd.read_csv(TABLES / "TNBC94_ecosystem_survival_models.tsv", sep="\t")
    survival = survival.loc[
        (survival["endpoint"] == "DRFS") & (survival["model"] == "adjusted")
    ].iloc[0]
    assert int(survival["n"]) == 92 and int(survival["events"]) == 22
    close(survival["hazard_ratio"], 0.756466172355962)
    close(survival["wald_p"], 0.160587363600399)

    text = MANUSCRIPT.read_text()
    abstract = text.split("## Abstract", 1)[1].split("## Introduction", 1)[0]
    abstract_words = len(re.findall(r"\b[\w–+-]+\b", abstract.split("**Keywords:**", 1)[0]))
    assert abstract_words <= 250, abstract_words
    assert "[@" not in text

    required = [
        ROOT / "manuscript" / "figure_legends.md",
        ROOT / "manuscript" / "key_resources_table.tsv",
        ROOT / "supplementary" / "table_index.md",
        ROOT / "manuscript" / "Nature_Portfolio_Reporting_Summary_draft.md",
        ROOT / "output" / "doc" / "npj_Breast_Cancer_main_manuscript.docx",
        ROOT / "output" / "pdf" / "npj_Breast_Cancer_main_manuscript.pdf",
        ROOT / "output" / "doc" / "npj_Breast_Cancer_review_manuscript_with_figures.docx",
        ROOT / "output" / "pdf" / "npj_Breast_Cancer_review_manuscript_with_figures.pdf",
        ROOT / "output" / "tables" / "Supplementary_Tables_S1-S10.xlsx",
        ROOT / "output" / "pdf" / "Supplementary_Information.pdf",
        ROOT / "results" / "figures" / "supplementary" / "FigureS2_factor_matching_and_confounders.pdf",
        ROOT / "results" / "figures" / "supplementary" / "FigureS4_Wu_spatial_complete_effects.pdf",
        ROOT / "results" / "figures" / "supplementary" / "FigureS5_Figshare_spatial_complete_effects.pdf",
        ROOT / "results" / "figures" / "supplementary" / "FigureS6_TNBC94_spatial_QC_and_complete_effects.pdf",
    ]
    required += [ROOT / "results" / "figures" / "main" / f"{stem}.pdf" for stem in [
        "Figure1_cross_cohort_and_cell_specificity",
        "Figure2_three_cohort_spatial_replication",
        "Figure3_spatial_context_and_robustness",
        "Figure4_treatment_boundary",
        "Figure5_survival_and_evidence_model",
    ]]
    required += [
        ROOT / "results" / "figures" / "supplementary" / name
        for name in [
            "FigureS1_discovery_QC_and_NMF_rank.pdf",
            "FigureS3_tumor_normal_boundary_robustness.pdf",
            "FigureS7_response_and_longitudinal_controls.pdf",
            "FigureS8_survival_controls_and_robustness.pdf",
            "FigureS9_celltype_specificity.pdf",
            "FigureS10_spatial_sensitivity.pdf",
        ]
    ]
    missing = [str(path.relative_to(ROOT)) for path in required if not path.exists()]
    assert not missing, f"Missing submission files: {missing}"

    illustrated = ROOT / "output" / "doc" / "npj_Breast_Cancer_review_manuscript_with_figures.docx"
    with zipfile.ZipFile(illustrated) as archive:
        embedded_images = [name for name in archive.namelist() if name.startswith("word/media/")]
    assert len(embedded_images) == 5, embedded_images
    for callout in ["(Fig. 1a–c", "(Fig. 2a–e", "(Fig. 3c,d", "(Fig. 4 ", "(Fig. 5 "]:
        assert callout in text, f"Missing main-text figure callout: {callout}"

    print("PASS: frozen headline values match source tables")
    print(f"PASS: Abstract contains {abstract_words} tokenized words (limit 250)")
    print("PASS: numbered submission source contains no unresolved citation keys")
    print("PASS: main figures, core supplementary figures, and submission files exist")
    print("PASS: illustrated manuscript contains 5 embedded images and main-text callouts")


if __name__ == "__main__":
    main()
