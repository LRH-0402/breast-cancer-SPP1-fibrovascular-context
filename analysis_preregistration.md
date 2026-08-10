# Analysis preregistration v0.1

Date frozen: 2026-07-21, before joining GSE246613 scores to pCR/RCB outcomes  
Primary target: iScience or a journal of comparable scope  
Primary disease setting: triple-negative breast cancer  
Secondary setting: pan-breast cancer generalization

## Study question

Does a patient-reproducible lipid-processing tumor-associated macrophage program form a spatially organized multicellular niche with malignant, fibroblast, endothelial, and T-cell states, and does this niche associate with breast cancer progression or treatment response beyond established clinical and immune features?

## Primary hypothesis

Across independent breast cancer single-cell cohorts, a continuous macrophage program enriched for lipid uptake, handling, storage, or catabolism is positively associated with an epithelial lipid-associated program at the patient level. In independent spatial cohorts, the macrophage program is non-randomly enriched near extracellular-matrix-producing fibroblasts and aggressive malignant regions and is associated with reduced local cytotoxic T-cell access.

## Primary outcomes

The primary single-cell outcome is the cross-cohort reproducibility of a macrophage gene program estimated without using clinical outcome labels. The primary spatial outcome is patient-level enrichment of the predefined macrophage program within a predefined multicellular neighborhood relative to a patient-stratified spatial null model. The primary clinical outcome is pathological complete response or residual cancer burden in a TNBC neoadjuvant cohort, conditional on data availability and sample size.

## Discovery and validation separation

GSE176078/SCP1039 is the initial discovery reference. At least one independent single-cell cohort will be used for program replication without re-estimating the program direction. Spatial discovery and spatial validation will be assigned before model fitting after completion of the metadata audit. Bulk score genes, directions, weights, and thresholds will be frozen before external validation.

## Statistical unit

Patients are the independent units for single-cell and bulk inference. Spatial spots or cells are nested within sections, and sections are nested within patients. Donors or independently treated culture batches are the independent units for in vitro experiments. Technical replicates, image fields, cells, and regions of interest do not increase biological sample size.

## Core models

Single-cell differential analyses will use patient-level pseudobulk or mixed models with patient as a random effect. Cross-dataset effects will be combined by random-effects meta-analysis when at least three cohorts are available. Spatial enrichment will be tested with patient-stratified label permutations or point-process null models that preserve tissue compartment and cell abundance. Clinical outcomes will be modeled with logistic regression for binary response, ordinal regression when valid for residual cancer burden, and Cox regression for time-to-event outcomes. All models will report effect sizes and 95% confidence intervals.

## Prespecified confounders

The minimum prespecified adjustment set includes molecular subtype, stage or tumor burden where available, treatment arm, tumor purity, broad immune infiltration, and dataset/platform. Hypoxia, necrosis, and general phagocytic activity will be evaluated as biological alternative explanations for the putative lipid-processing program.

## Multiple testing

Gene- and pathway-level discovery analyses will control the Benjamini–Hochberg false-discovery rate. A single primary test will be declared for each gate before execution. Secondary and exploratory results will be labeled accordingly.

## Gate 1 decision rule

Gate 1 passes when the program shows directionally concordant patient-level enrichment in at least three cohorts or in two cohorts plus one orthogonal validation dataset; the meta-analytic FDR is below 0.05; leave-one-dataset-out analysis retains at least 60% of the locked core genes; and the result is not abolished after adjustment for sequencing depth, hypoxia, and broad macrophage abundance.

## Gate 2 decision rule

Gate 2 passes when a predefined multicellular neighborhood is enriched in two independent spatial cohorts, remains significant after controlling for tissue compartment and cell abundance, and shows a consistent effect direction in patient-level estimates.

## Gate 3 decision rule

Gate 3 passes when one tumor-intrinsic candidate program shows cross-cohort association, spatial precedence or distance dependence, an experimentally plausible mediator, and a feasible perturbation/rescue design. Trajectory inference or ligand–receptor prediction alone cannot satisfy this gate.

## Gate 4 decision rule

Gate 4 passes when the frozen ecosystem score retains association in an external cohort and improves at least one prespecified performance measure beyond subtype, purity, and broad immune infiltration. Optimal cut points will not be selected in validation data.

## Prohibited inference

Cross-sectional transcriptomic association cannot establish lineage, metabolite flux, cell–cell communication, or causality. RNA velocity will only be used when raw spliced and unspliced counts are available and will be interpreted as state dynamics, not lineage proof. CellChat, NicheNet, LIANA, SCENIC, and related tools are candidate-prioritization methods rather than experimental validation.

## Amendment policy

Any change after the date is frozen will be recorded in an amendment table with date, rationale, affected hypothesis, and whether the change was made before or after viewing the relevant outcome.

| Date | Amendment | Rationale | Outcome viewed? |
|---|---|---|---|
| 2026-07-21 | Separate discovery P3 and P4 axes; designate the 13-gene P4 core as primary and the eight-gene P3 core as secondary | Independent outcome-blinded NMF showed a reproducible lysosomal-associated axis and a partially distinct SPP1–matrix–lipid axis | No GSE246613 score–outcome association viewed; only response-class counts were known |
| 2026-07-21 | Interpret GSE114725 tumor–normal comparison as a boundary analysis rather than proof of universal tumor induction | Only two of four paired patients had higher tumor P4 scores | Not applicable; tissue identity is the tested exposure and the negative result is retained |
| 2026-07-21 | For spatial validation sections lacking pathology annotations, stratify permutations by within-section malignant-epithelial score quintile, pooling strata with fewer than 20 spots | The validation RDS objects contain coordinates and counts but no harmonized pathology labels; this preserves a frozen proxy for tumor–stroma composition | No Figshare P1–P8 program score or spatial association viewed |
| 2026-07-21 | Refine the supported spatial phenotype to a fibrovascular niche and retain the failed immune-exclusion prediction | The frozen primary fibrovascular test reproduced, whereas malignant-neighbor enrichment was null or negative and CD8/cytotoxic association was positive in all three spatial cohorts | Yes; this is a result interpretation, not an analysis change |
| 2026-07-21 | Define the TNBC94 patient-level ecosystem score as the equal-weight mean of standardized P4 residual and fibrovascular expression; prespecify continuous adjusted DRFS Cox analysis | The 94-patient spatial pseudobulk and clinical data permit a non-spatial patient-level external test while full spot objects download | No TNBC94 program score or survival association viewed; only endpoint event counts and covariate distributions were audited |
| 2026-07-21 | Extend the frozen six-neighbor spatial test to the one pathologist-annotated array per TNBC94 patient; use dominant pixel annotation for filtering and permutation strata | Full non-batch-corrected spot counts, coordinates, and fractional pathology annotations became available after the patient-level outcome test | Patient-level pseudobulk survival result was viewed, but no TNBC94 spot-level P4, fibrovascular, or neighborhood result was viewed |
