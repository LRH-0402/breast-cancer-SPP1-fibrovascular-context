# Supplementary figure legends

## Figure S1. Discovery-cohort quality control and NMF rank audit

(A) Cophenetic correlation and consensus silhouette across candidate NMF ranks;
the vertical line marks the selected six-program solution. (B) Number of
author-annotated macrophages per GSE176078 patient, colored by breast cancer
subtype. (C) Distribution of patient-level P4 exposure in the discovery
cohort. NMF was performed after patient balancing so that patients with larger
cell yields did not dominate factor discovery.

## Figure S2. Independent factor matching and patient-level confounder audit

(A) Complete background-corrected best-factor matches for discovery programs
P1-P6 in three independently factorized validation cohorts. Color denotes
-log10(FDR). (B) Fraction of standardized patient-level P4 variation remaining
after joint adjustment for macrophage identity, hypoxia, and library depth.
(C) Patient-level correlations of raw P4 with the three adjustment variables.
(D) Number of discovery or validation datasets in which each locked P4 gene
occurred among the matched factor's top-loading genes. Recurrence establishes a
shared expression axis but does not make raw P4 biologically specific.

## Figure S3. Tumor-normal boundary-test robustness

(A) Paired patient-level P4 scores in normal and tumor immune cells from four
patients in GSE114725. (B) Median paired tumor-minus-normal difference for the
complete P4 score and after omitting each gene in turn. Every score retained
the same two-of-four directional split. (C) Coefficients and 95% confidence
intervals from the exploratory patient-tissue aggregate model. The tumor term
was not significant after adjustment for hypoxia, macrophage identity, and
library depth.

## Figure S4. Complete Wu discovery spatial effects and quality controls

(A) Fibrovascular neighbor effect minus the pathology-stratified null median
for all six patients. (B) Complete target-specific excess-effect matrix, with
patients ordered by the primary effect. (C) Retained and top-quartile P4 spot
counts by patient. (D) Primary excess effect versus the within-patient P4
residual model R-squared. Spatial effects compare the mean six-neighbor target
of top-quartile P4-residual spots with all remaining spots. (E) H&E image for
the representative 1160920F section used in the main spatial map. (F)
Author-provided pathology classifications overlaid on the same tissue section.

## Figure S5. Complete independent spatial-validation effects and quality controls

(A) Fibrovascular neighbor effect minus the composition-stratified null median
for all eight sections. (B) Complete target-specific excess-effect matrix,
with sections ordered by the primary effect. (C) Retained and top-quartile P4
spot counts by section. (D) Primary excess effect versus the within-section P4
residual model R-squared. Section identity is the inferential unit because
independent patient identity could not be verified.

## Figure S6. TNBC94 spatial quality control and complete patient effects

(A) Primary excess fibrovascular effect for all 94 patients. (B) Distribution
of retained spots per patient. (C) Median detected genes versus median UMI per
array; point area denotes retained spots. (D) Number of patients with a
positive excess effect for each spatial target; the dashed line marks half of
the cohort. The primary effect was positive in 82 of 94 patients.

## Figure S7. Neoadjuvant-response and longitudinal specificity controls

(A) P4 odds ratios for pathological response after one-at-a-time adjustment
for prespecified technical or biological covariates. These small models are
sensitivity analyses, not a jointly adjusted prediction model. (B) Univariable
response associations for frozen negative-control programs. (C) Median paired
change, scaled by the baseline standard deviation, for seven frozen programs
and three treatment contrasts. The decrease after radiotherapy plus
pembrolizumab was shared by P4, hypoxia, and interferon programs.

## Figure S8. Survival specificity controls and leave-one-gene-out robustness

(A) Adjusted continuous ecosystem-score hazard ratios for distant relapse-free,
invasive disease-free, and overall survival. (B) Adjusted distant-relapse
associations for the P4 residual, fibrovascular, macrophage-identity, and
hypoxia components. (C) Adjusted distant-relapse associations after omitting
each P4 gene in turn. All leave-one-gene-out estimates remained below one and
non-significant. Adding the frozen ecosystem score changed concordance by
0.018 and did not significantly improve the covariate-only model
(likelihood-ratio P = 0.173).

## Figure S9. Cell-compartment and gene-level expression boundaries

(A) Patient-level mean P4 score across six consolidated cell compartments in
GSE176078. Points are patient–compartment aggregates. (B) Mean normalized
expression of each locked P4 gene by compartment. (C) Fraction of cells with
detected expression of each P4 gene by compartment. (D) Empirical gene
classification based on the macrophage-to-strongest-other-compartment mean
expression ratio. A gene was macrophage-dominant only when macrophages had the
highest mean and the ratio was at least 1.5; remaining classes describe
myeloid-enriched/shared, multicellular/stromal, or lymphocyte-weighted
expression.

## Figure S10. Spatial robustness and post hoc T-cell-state boundary

(A) Fraction of spatial units with a positive fibrovascular excess effect after
varying neighborhood size, P4-high threshold, residual covariates, and
permutation stratification. The fibroblast-conditioned model is intentionally
shown because it conditions on a component of the target and attenuates the
effect. (B) Positive fraction after omitting each P4 gene in turn. Sensitivity
analyses used 250 permutations per unit and setting and assess directional
consistency; confirmatory P values derive only from the frozen
2,000-permutation analysis. (C) Fraction of spatial units with positive excess
expression for CD8/cytotoxic, exhaustion-associated, IFNG-response, and
TCF7/progenitor-associated modules. Panel C is post hoc and mixed spatial spots
do not identify single-cell T-cell phenotypes.
