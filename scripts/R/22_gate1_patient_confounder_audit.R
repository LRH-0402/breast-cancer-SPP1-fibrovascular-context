#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(data.table)
  library(yaml)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grep("^--file=", args)])
root <- normalizePath(file.path(dirname(script_arg), "..", ".."))
table_dir <- file.path(root, "results", "tables")
log_dir <- file.path(root, "logs")
lock <- read_yaml(file.path(root, "config", "locked_macrophage_programs_v1.yml"))

modules <- list(
  P4_primary = unlist(lock$primary_program$genes),
  macrophage_identity = c("C1QA", "C1QB", "C1QC", "CD68", "MSR1"),
  hypoxia = c(
    "HIF1A", "EGLN1", "CA9", "VEGFA", "LDHA", "PDK1", "SLC2A1",
    "BNIP3", "NDRG1", "ENO1", "ALDOA", "PGK1", "HK2", "PFKP"
  )
)

score_cells <- function(object, cohort, patient_column, keep = rep(TRUE, ncol(object))) {
  counts <- GetAssayData(object, assay = "RNA", layer = "counts")[, keep, drop = FALSE]
  metadata <- object@meta.data[keep, , drop = FALSE]
  library_size <- Matrix::colSums(counts)
  result <- data.table(
    patient = as.character(metadata[[patient_column]]),
    total_counts = as.numeric(library_size)
  )
  for (module in names(modules)) {
    present <- intersect(modules[[module]], rownames(counts))
    values <- counts[present, , drop = FALSE] %*% Diagonal(x = 10000 / library_size)
    values@x <- log1p(values@x)
    result[[module]] <- Matrix::colMeans(values)
  }
  result[, .(
    cohort = cohort,
    macrophage_cells = .N,
    P4_primary = median(P4_primary),
    macrophage_identity = median(macrophage_identity),
    hypoxia = median(hypoxia),
    median_total_counts = median(total_counts)
  ), by = patient]
}

message("Scoring GSE176078")
gse176078 <- readRDS(file.path(root, "data", "interim", "GSE176078_macrophages_seurat.rds"))
patient_tables <- list(score_cells(gse176078, "GSE176078", "orig.ident"))
rm(gse176078); gc()

message("Scoring GSE161529 macrophage-rich cluster")
gse161529 <- readRDS(file.path(root, "data", "interim", "GSE161529_TNBCSub_updated.rds"))
patient_tables[[2]] <- score_cells(
  gse161529, "GSE161529", "group", as.character(gse161529$seurat_clusters) == "1"
)
rm(gse161529); gc()

message("Loading blinded aggregate scores for GSE246613 and GSE114725")
gse246613 <- fread(file.path(table_dir, "GSE246613_locked_program_by_patient_timepoint_blinded.tsv"))[
  treatment == "Base",
  .(
    cohort = "GSE246613", patient = cohort, macrophage_cells = macrophage_like_cells,
    P4_primary, macrophage_identity, hypoxia, median_total_counts = median_n_counts
  )
]
gse114725 <- fread(file.path(table_dir, "GSE114725_locked_program_scores_by_patient_tissue.tsv"))[
  tissue == "TUMOR",
  .(
    cohort = "GSE114725", patient, macrophage_cells, P4_primary,
    macrophage_identity, hypoxia, median_total_counts
  )
]
patient_tables[[3]] <- gse246613
patient_tables[[4]] <- gse114725
patients <- rbindlist(patient_tables, use.names = TRUE)
fwrite(patients, file.path(table_dir, "gate1_patient_level_locked_P4_scores.tsv"), sep = "\t")

audit <- rbindlist(lapply(split(patients, patients$cohort), function(d) {
  scaled <- copy(d)
  for (column in c("P4_primary", "macrophage_identity", "hypoxia", "median_total_counts")) {
    scaled[[column]] <- as.numeric(scale(scaled[[column]]))
  }
  model <- lm(
    P4_primary ~ macrophage_identity + hypoxia + median_total_counts,
    data = scaled
  )
  data.table(
    cohort = unique(d$cohort),
    patients = nrow(d),
    P4_sd = sd(d$P4_primary),
    adjusted_residual_sd = sd(resid(model)),
    residual_sd_fraction = sd(resid(model)) / sd(scaled$P4_primary),
    model_r_squared = summary(model)$r.squared,
    cor_P4_macrophage_identity = cor(d$P4_primary, d$macrophage_identity),
    cor_P4_hypoxia = cor(d$P4_primary, d$hypoxia),
    cor_P4_depth = cor(d$P4_primary, d$median_total_counts)
  )
}))
fwrite(audit, file.path(table_dir, "gate1_patient_level_confounder_audit.tsv"), sep = "\t")

writeLines(capture.output(sessionInfo()), file.path(log_dir, "22_gate1_patient_confounder_audit_sessionInfo.log"))
print(audit)

