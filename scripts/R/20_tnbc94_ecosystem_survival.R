#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(yaml)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grep("^--file=", args)])
root <- normalizePath(file.path(dirname(script_arg), "..", ".."))
raw_dir <- file.path(root, "data", "raw", "spatial", "tnbc_zenodo_14204217")
clinical_dir <- file.path(root, "data", "interim", "spatial_audit", "tnbc_clinical", "Clinical")
table_dir <- file.path(root, "results", "tables")
figure_dir <- file.path(root, "results", "figures")
log_dir <- file.path(root, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

counts <- readRDS(file.path(raw_dir, "PB_count.RDS"))
clinical <- readRDS(file.path(clinical_dir, "Clinical.RDS"))
lock <- read_yaml(file.path(root, "config", "locked_macrophage_programs_v1.yml"))
spec <- read_yaml(file.path(root, "config", "spatial_analysis_spec_v1.yml"))
stopifnot(setequal(colnames(counts), rownames(clinical)))
clinical <- clinical[colnames(counts), , drop = FALSE]

library_size <- colSums(counts)
log_cpm <- log1p(t(t(counts) / library_size * 1e6))
modules <- c(
  list(
    P4_primary = unlist(lock$primary_program$genes),
    macrophage_identity = c("C1QA", "C1QB", "C1QC", "CD68", "MSR1"),
    hypoxia = c(
      "HIF1A", "EGLN1", "CA9", "VEGFA", "LDHA", "PDK1", "SLC2A1",
      "BNIP3", "NDRG1", "ENO1", "ALDOA", "PGK1", "HK2", "PFKP"
    )
  ),
  lapply(spec$spatial_signatures[c("fibroblast_ecm", "endothelial_angiogenic")], unlist)
)

module_score <- function(matrix, genes) {
  present <- intersect(genes, rownames(matrix))
  if (length(present) == 0) stop("No module genes present")
  colMeans(matrix[present, , drop = FALSE])
}
z <- function(x) as.numeric(scale(x))

scores <- data.table(patient = colnames(counts), library_size = library_size)
for (module in names(modules)) scores[[module]] <- module_score(log_cpm, modules[[module]])
p4_model <- lm(
  P4_primary ~ macrophage_identity + hypoxia + log1p(library_size),
  data = scores
)
scores[, P4_residual := resid(p4_model)]
scores[, fibrovascular := (z(fibroblast_ecm) + z(endothelial_angiogenic)) / 2]
scores[, ecosystem_score := z((z(P4_residual) + z(fibrovascular)) / 2)]

for (gene in unlist(lock$primary_program$genes)) {
  remaining <- setdiff(unlist(lock$primary_program$genes), gene)
  p4_loo <- module_score(log_cpm, remaining)
  loo_model <- lm(p4_loo ~ scores$macrophage_identity + scores$hypoxia + log1p(library_size))
  scores[[paste0("ecosystem_without_", gene)]] <- z(
    (z(resid(loo_model)) + z(scores$fibrovascular)) / 2
  )
}

score_output <- copy(scores)
score_output[, patient := as.character(patient)]
fwrite(score_output, file.path(table_dir, "TNBC94_frozen_ecosystem_scores.tsv"), sep = "\t")

clinical_key <- data.table(
  patient = rownames(clinical),
  TIME_pathologist = clinical$TIME_classes.by.pathologist,
  TIME_pseudobulk = clinical$TIME_classes_expression_global_pseudobulk,
  spatial_archetype = clinical[["Spatial archetypes_defined_on_ST_global_pseudobulk"]],
  DRFS_time = clinical$DRFS[, "time"], DRFS_status = clinical$DRFS[, "status"],
  iDFS_time = clinical$iDFS[, "time"], iDFS_status = clinical$iDFS[, "status"],
  OS_time = clinical$OS[, "time"], OS_status = clinical$OS[, "status"]
)
fwrite(clinical_key, file.path(table_dir, "TNBC94_clinical_ecology_key.tsv"), sep = "\t")

analysis <- data.frame(
  patient = rownames(clinical),
  age10 = clinical$Age_at_diagnosis / 10,
  T_stage = clinical$T_TNM,
  N_stage = clinical$N_TNM,
  DRFS_time = clinical$DRFS[, "time"],
  DRFS_status = clinical$DRFS[, "status"],
  iDFS_time = clinical$iDFS[, "time"],
  iDFS_status = clinical$iDFS[, "status"],
  OS_time = clinical$OS[, "time"],
  OS_status = clinical$OS[, "status"],
  scores[, setdiff(names(scores), c("patient", "library_size")), with = FALSE],
  check.names = FALSE
)

fit_endpoint <- function(prefix, exposure, adjusted = TRUE) {
  survival_term <- sprintf("Surv(%s_time, %s_status)", prefix, prefix)
  rhs <- if (adjusted) {
    paste(exposure, "+ age10 + T_stage + N_stage")
  } else exposure
  model <- coxph(as.formula(paste(survival_term, "~", rhs)), data = analysis, x = TRUE)
  coefficient <- summary(model)$coefficients[exposure, ]
  interval <- summary(model)$conf.int[exposure, ]
  data.table(
    endpoint = prefix,
    exposure = exposure,
    model = ifelse(adjusted, "adjusted", "unadjusted"),
    n = model$n,
    events = model$nevent,
    hazard_ratio = interval["exp(coef)"],
    ci_low = interval["lower .95"],
    ci_high = interval["upper .95"],
    wald_p = coefficient["Pr(>|z|)"],
    concordance = summary(model)$concordance[1],
    ph_p = cox.zph(model)$table[exposure, "p"]
  )
}

main_results <- rbindlist(lapply(c("DRFS", "iDFS", "OS"), function(endpoint) {
  rbind(fit_endpoint(endpoint, "ecosystem_score", FALSE),
        fit_endpoint(endpoint, "ecosystem_score", TRUE))
}))
secondary <- main_results$model == "adjusted" & main_results$endpoint %in% c("iDFS", "OS")
main_results[secondary, secondary_fdr := p.adjust(wald_p, method = "BH")]

base <- coxph(Surv(DRFS_time, DRFS_status) ~ age10 + T_stage + N_stage, data = analysis)
full <- coxph(
  Surv(DRFS_time, DRFS_status) ~ ecosystem_score + age10 + T_stage + N_stage,
  data = analysis
)
comparison <- anova(base, full, test = "LRT")
primary_increment <- data.table(
  base_concordance = summary(base)$concordance[1],
  full_concordance = summary(full)$concordance[1],
  delta_concordance = summary(full)$concordance[1] - summary(base)$concordance[1],
  likelihood_ratio_chisq = comparison[2, "Chisq"],
  likelihood_ratio_p = comparison[2, "Pr(>|Chi|)"]
)

controls <- rbindlist(lapply(
  c("P4_residual", "fibrovascular", "macrophage_identity", "hypoxia"),
  function(exposure) fit_endpoint("DRFS", exposure, TRUE)
))
controls[, control_fdr := p.adjust(wald_p, method = "BH")]

loo_results <- rbindlist(lapply(
  grep("^ecosystem_without_", names(analysis), value = TRUE),
  function(exposure) fit_endpoint("DRFS", exposure, TRUE)
))

fwrite(main_results, file.path(table_dir, "TNBC94_ecosystem_survival_models.tsv"), sep = "\t")
fwrite(primary_increment, file.path(table_dir, "TNBC94_ecosystem_model_increment.tsv"), sep = "\t")
fwrite(controls, file.path(table_dir, "TNBC94_ecosystem_specificity_controls.tsv"), sep = "\t")
fwrite(loo_results, file.path(table_dir, "TNBC94_ecosystem_leave_one_P4_gene_out.tsv"), sep = "\t")

pdf(file.path(figure_dir, "TNBC94_ecosystem_DRFS_forest_and_KM.pdf"), width = 9, height = 4.2)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
adjusted_row <- main_results[endpoint == "DRFS" & model == "adjusted"]
plot(
  1, adjusted_row$hazard_ratio, log = "y", ylim = range(c(adjusted_row$ci_low, adjusted_row$ci_high, 1)),
  xlim = c(0.5, 1.5), xaxt = "n", xlab = "", ylab = "Hazard ratio (log scale)", pch = 19,
  main = "Adjusted DRFS association"
)
segments(1, adjusted_row$ci_low, 1, adjusted_row$ci_high, lwd = 2)
abline(h = 1, lty = 2, col = "grey40")
axis(1, at = 1, labels = "Ecosystem score\n(per SD)")
analysis$ecosystem_group <- ifelse(
  analysis$ecosystem_score >= median(analysis$ecosystem_score), "Higher", "Lower"
)
km <- survfit(Surv(DRFS_time, DRFS_status) ~ ecosystem_group, data = analysis)
plot(km, col = c("#2C7BB6", "#D7191C"), lwd = 2, xlab = "Years", ylab = "Distant relapse-free survival", main = "Descriptive median split")
legend("bottomleft", legend = c("Higher", "Lower"), col = c("#D7191C", "#2C7BB6"), lwd = 2, bty = "n")
dev.off()

writeLines(capture.output(sessionInfo()), file.path(log_dir, "20_tnbc94_ecosystem_survival_sessionInfo.log"))
cat("Primary and secondary models\n")
print(main_results)
cat("\nPrimary model increment\n")
print(primary_increment)
cat("\nSpecificity controls\n")
print(controls)
