# Limited orthogonal validation plan

## Purpose and claim boundary

The computational evidence supports a recurrent macrophage-associated P4 program and a fibrovascular, immune-infiltrated spatial context. It does not establish that every P4 transcript is macrophage-derived, that the implicated cells directly contact one another, or that SPP1 causes the neighborhood. The most efficient experimental addition should therefore validate cellular source and spatial co-occurrence first. Functional perturbation should be presented as a bounded test of one candidate mediator, not as proof that the entire 13-gene program is an SPP1-driven pathway.

## Tier 1 tissue validation

Use an independent retrospective cohort of formalin-fixed paraffin-embedded breast cancers, preferably enriched for TNBC while retaining a small non-TNBC comparison group. A practical pilot is 10–12 tumors, followed by a confirmatory set of at least 20 additional tumors if staining performance and between-patient variance are acceptable. Patient, not image field or cell, is the independent unit.

RNAscope or equivalent multiplex RNA in situ hybridization should combine a macrophage identity marker such as *C1QA* or *CD68* with *SPP1* and one second P4 feature such as *FN1*, *FABP5*, or *HMOX1*. Multiplex immunofluorescence should identify macrophages, fibroblasts, vessels, and cytotoxic lymphocytes with a compact panel such as CD68, FAP or COL1A1, CD31, CD8, and granzyme B. Panel design must be optimized on control tonsil, normal breast, and known positive tumor tissue before analysis of the study cohort.

The primary tissue endpoint is the patient-level excess proximity of P4-positive macrophages to fibrovascular structures relative to label-permuted macrophages within the same pathology compartment. A prespecified 50 µm radius or a nearest-neighbor distance may be used, but the radius and segmentation rules must be frozen before outcome analysis. Secondary endpoints are proximity to CD8-positive or granzyme-B-positive cells, the fraction of P4-positive macrophages at tumor–stroma boundaries, and concordance between RNA and protein-level P4 markers.

Segmentation and phenotype calls should be performed blinded to clinical variables. Each patient should contribute the same maximum number or area of quality-controlled regions. The primary analysis should summarize fields within patient and use either a patient-level paired test or a mixed model with patient as a random intercept. Pathology compartment, tumor subtype, tissue area, and macrophage density should be treated as prespecified covariates or stratification factors. Exact cell counts, excluded regions, failed stains, and all image-processing settings should be reported.

Negative controls should include bacterial negative-control probes, positive housekeeping probes, single-stain controls, isotype or secondary-only controls where applicable, and a spatial target not predicted to reproduce, such as malignant-epithelial enrichment. The locked P3 program can be assessed as a biological specificity comparator if sufficient markers can be measured without compromising the primary panel.

## Tier 2 bounded functional validation

If a functional experiment is feasible, use primary human monocyte-derived macrophages from at least three independent donors. THP-1-derived macrophages may be used for assay development but should not be the only biological system. Expose macrophages to normoxia versus hypoxia and to a defined lipid condition, then quantify the complete P4 panel by targeted RNA measurement rather than relying on SPP1 alone. Conditions should be selected before measuring fibroblast or endothelial outcomes.

Conditioned medium or transwell co-culture can test whether the induced macrophage condition changes primary breast fibroblast matrix expression or endothelial organization. Suggested readouts are fibroblast *COL1A1*, *POSTN*, and *FAP* expression; collagen deposition; endothelial migration or tube-network metrics; and CD8 T-cell chemotaxis where donor material permits. Independent macrophage donors or independently repeated differentiations are the biological replicates. Wells, images, and segmented cells are technical replicates and must not inflate sample size.

SPP1 knockdown or a validated neutralizing reagent can test whether one component contributes to the phenotype. A rescue arm with recombinant SPP1 is required for a strong mediator claim. Failure of SPP1 perturbation to collapse the full P4 state should be retained because P4 is a multi-gene co-expression program rather than an assumed linear pathway.

## Statistical and reporting rules

The tissue proximity endpoint and one functional endpoint should be designated primary before data collection. All other readouts are secondary and should use false-discovery-rate control. Effect sizes and confidence intervals should accompany exact P values. Formal sample-size calculation should use variance observed in the blinded pilot; no retrospective power calculation should be used to reinterpret a negative result. Raw microscopy fields, segmentation masks, analysis code, and uncropped assay images should be archived.

## Stopping and interpretation rules

If fewer than 80% of tumors pass staining and segmentation quality control, optimize the assay before expanding the cohort. If P4-positive macrophage identity cannot be distinguished from extracellular or non-myeloid SPP1 signal, restrict the conclusion to a macrophage-associated tissue program. If fibrovascular proximity reproduces but CD8 proximity does not, retain the fibrovascular result and treat immune coexistence as cohort-dependent. Only a perturbation with a prespecified phenotype and successful rescue can support causal language about SPP1.

Tier 1 can be incorporated as a new main validation figure or a focused supplementary figure. Tier 2 should only be added when the donor-level design and rescue are complete; a small unreplicated cell-line assay would weaken rather than strengthen the current evidence-bounded manuscript.
