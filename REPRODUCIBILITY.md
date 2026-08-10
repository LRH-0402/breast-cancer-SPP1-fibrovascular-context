# Breast cancer fibrovascular macrophage ecosystem project

This repository contains the public-data analysis, frozen specifications, machine-readable results, figures, supplementary material, and manuscript for a study of a recurrent SPP1–matrix–lipid macrophage program in breast cancer.

The current submission target is *npj Breast Cancer*. The central conclusion is deliberately bounded: the P4 program recurs across independently factorized single-cell cohorts and occupies reproducible fibrovascular, CD8/cytotoxic-infiltrated spatial neighborhoods, but it is not a universal marker of tumor induction, immune exclusion, treatment resistance, or adverse prognosis.

## Current manuscript

Title: “A recurrent SPP1 matrix macrophage program marks fibrovascular immune infiltrated breast cancer neighborhoods”

The journal-formatted source is `manuscript/manuscript_npj_breast_cancer.md`. Submission-ready Word and PDF files are in `output/doc/` and `output/pdf/`. The original iScience-format files are retained as historical working outputs and are not the current target version.

## Repository structure

```text
config/              Frozen genes, covariates, null models, and clinical specifications
data/raw/            Immutable public inputs; excluded from a portable code release
data/interim/        Unpacked and intermediate objects
data/derived/        Analysis-ready score tables too large for the manuscript package
environment/         Python and R dependency specifications
metadata/            Dataset audits and frozen cell/spot score metadata
scripts/R/           R analysis stages
scripts/python/      Python analysis, figure, document, and audit stages
results/tables/      Machine-readable statistical outputs
results/figures/     Main and supplementary artwork
manuscript/          Article source, references, legends, cover letter, and compliance report
supplementary/       Supplementary legends, table index, and limited validation plan
output/              Submission files, tables, and reproducibility manifests
logs/                Session information and auditable stage logs
```

## Evidence gates and final decisions

1. Cross-cohort macrophage-program recurrence: passed. P4 independently matched factors in GSE161529, GSE246613, and GSE114725 after discovery in GSE176078.
2. Composition-preserving spatial replication: passed. The frozen fibrovascular test reproduced in the Wu, eight-section Figshare, and 94-patient TNBC spatial series.
3. Causal upstream mechanism: not established. No velocity, ligand–receptor prediction, or cross-sectional association is presented as causal evidence.
4. Universal adverse clinical interpretation: not supported. Frozen response and survival analyses were retained as negative boundary tests rather than used to optimize a prognostic signature.

The preregistration and all dated amendments are in `analysis_preregistration.md`.

## Software environments

Analysis environment:

```bash
python3 -m venv .venv
.venv/bin/pip install -r environment/requirements.txt
```

Manuscript and PDF environment:

```bash
.venv/bin/pip install -r environment/requirements-manuscript.txt
```

R 4.5.3 and the direct R dependencies listed in `environment/R_requirements.txt` were used. Complete `sessionInfo()` outputs are retained in `logs/*sessionInfo.log`.

## Input data

All inputs are public. Dataset accession, scientific role, access URL, status, and audited sample counts are recorded in `data_manifest.tsv`. Raw files must be placed at the paths expected by the scripts; the authoritative local paths and SHA-256 values are written to `output/reproducibility/raw_data_sha256.tsv` by the manifest builder.

The full raw-data footprint is approximately 19 GB. Large public objects are not duplicated in the portable submission package.

## Reproduction levels

### Verify frozen claims and submission files

This level does not rerun the large single-cell or spatial computations.

```bash
.venv/bin/python scripts/python/25_submission_consistency_audit.py
python scripts/python/32_prepare_npj_source.py
python scripts/python/34_audit_npj_submission.py
```

### Rebuild figures, supplementary files, and the npj submission

Run after the machine-readable tables in `results/tables/` are available.

```bash
.venv/bin/python scripts/python/21_build_main_figures.py
.venv/bin/python scripts/python/28_build_core_supplementary_figures.py
.venv/bin/python scripts/python/30_build_remaining_supplementary_figures.py
python scripts/python/29_build_supplementary_tables_xlsx.py
python scripts/python/31_build_supplementary_information_pdf.py
python scripts/python/32_prepare_npj_source.py
python scripts/python/33_build_npj_submission.py
python scripts/python/34_audit_npj_submission.py
python scripts/python/35_build_reproducibility_manifests.py
python scripts/python/36_build_npj_submission_package.py
python scripts/python/38_build_public_repository_release.py
```

The document scripts require LibreOffice or another DOCX-to-PDF converter for the final PDF export. DOCX files themselves are generated directly by Python.

### Recompute the complete analysis from public inputs

Run from the repository root. Large stages can require substantial memory and several hours.

```bash
Rscript scripts/R/01_prepare_gse176078.R
Rscript scripts/R/02_nmf_macrophage_pilot.R
Rscript scripts/R/03_project_programs_gse161529.R
.venv/bin/python scripts/python/04_audit_gse246613.py
.venv/bin/python scripts/python/05_project_programs_gse246613.py
Rscript scripts/R/06_independent_nmf_gse161529.R
.venv/bin/python scripts/python/07_independent_nmf_gse246613_baseline.py
.venv/bin/python scripts/python/07_independent_nmf_gse246613_baseline.py --scope broad
.venv/bin/python scripts/python/08_profile_gse246613_myeloid_markers.py
.venv/bin/python scripts/python/09_audit_gse114725.py
.venv/bin/python scripts/python/10_independent_nmf_gse114725.py
Rscript scripts/R/11_formal_match_gse161529.R
.venv/bin/python scripts/python/11_finalize_gate1_matching.py
.venv/bin/python scripts/python/12_gse114725_paired_locked_programs.py
.venv/bin/python scripts/python/13_score_gse246613_locked_blinded.py
.venv/bin/python scripts/python/14_analyze_gse246613_clinical.py
.venv/bin/python scripts/python/15_gse246613_longitudinal_specificity.py
Rscript scripts/R/16_audit_spatial_cohorts.R
.venv/bin/python scripts/python/17_wu_visium_spatial_gate2.py
Rscript scripts/R/18_score_figshare_spatial_sections.R
.venv/bin/python scripts/python/19_figshare_spatial_gate2_validation.py
Rscript scripts/R/20_tnbc94_ecosystem_survival.R
Rscript scripts/R/22_gate1_patient_confounder_audit.R
Rscript scripts/R/23_score_tnbc94_annotated_spatial_arrays.R
.venv/bin/python scripts/python/24_tnbc94_spatial_gate2_extension.py
```

Then run the figure and submission commands from the preceding section.

## Reproducibility safeguards

- Outcome labels were withheld during factor discovery and gene locking.
- Frozen YAML files record genes, covariates, thresholds, permutation strata, statistical units, and prohibited inferences.
- Patients or sections, rather than cells or spots, are the independent inferential units.
- Spatial nulls preserve pathology class or malignant-expression composition.
- Negative findings and leave-one-gene-out checks are retained.
- `scripts/python/25_submission_consistency_audit.py` verifies headline values against source tables.
- `scripts/python/34_audit_npj_submission.py` verifies journal structure, word and figure limits, artwork properties, and output integrity.
- `scripts/python/35_build_reproducibility_manifests.py` records raw-input and release-file SHA-256 values.
- `scripts/python/38_build_public_repository_release.py` creates a deterministic, path-sanitized GitHub/Zenodo archive with an independent standard-library verifier.

## Remaining author-supplied information

The scientific and formatting files are complete, but submission still requires real author names and affiliations, corresponding-author details, CRediT contributions and guarantor, funding information, confirmation of competing interests and author approval, and a permanent public code-repository URL.

Enter these once in a private copy of `manuscript/author_metadata_template.yml`, following `manuscript/author_metadata_instructions.md`, then run `scripts/python/39_finalize_submission_metadata.py`. The validator refuses incomplete declarations or unresolved placeholders and generates separate final Markdown, DOCX, and PDF files without overwriting the audited templates.
