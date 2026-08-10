#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grep("^--file=", args)])
root <- normalizePath(file.path(dirname(script_arg), "..", ".."))
object_root <- file.path(root, "data", "interim", "tnbc94_spatial_objects", "Robjects")
clinical_root <- file.path(root, "data", "interim", "spatial_audit", "tnbc_clinical", "Clinical")
output_dir <- file.path(root, "data", "derived", "tnbc94_annotated_spot_scores")
table_dir <- file.path(root, "results", "tables")
log_dir <- file.path(root, "logs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

lock <- read_yaml(file.path(root, "config", "locked_macrophage_programs_v1.yml"))
spec <- read_yaml(file.path(root, "config", "spatial_analysis_spec_v1.yml"))
tcell_spec <- read_yaml(file.path(root, "config", "tcell_state_sensitivity_v1.yml"))
ids <- readRDS(file.path(clinical_root, "ids.RDS"))
annotated <- ids[ids$hasAnnot %in% TRUE, , drop = FALSE]
stopifnot(nrow(annotated) == 94, length(unique(annotated$id)) == 94)

modules <- c(
  list(
    P4_primary = unlist(lock$primary_program$genes),
    P3_secondary = unlist(lock$secondary_program$genes),
    macrophage_identity = c("C1QA", "C1QB", "C1QC", "CD68", "MSR1"),
    hypoxia = c(
      "HIF1A", "EGLN1", "CA9", "VEGFA", "LDHA", "PDK1", "SLC2A1",
      "BNIP3", "NDRG1", "ENO1", "ALDOA", "PGK1", "HK2", "PFKP"
    )
  ),
  lapply(spec$spatial_signatures, unlist),
  lapply(tcell_spec$tcell_state_signatures, unlist)
)

coverage <- list()
audit <- list()
for (row_index in seq_len(nrow(annotated))) {
  patient <- as.character(annotated$id[row_index])
  array <- rownames(annotated)[row_index]
  message("Scoring TNBC", patient, " array ", array)
  counts_object <- readRDS(file.path(object_root, "countsNonCorrected", paste0("TNBC", patient, ".RDS")))
  annotation_object <- readRDS(file.path(object_root, "annotsBySpot", paste0("TNBC", patient, ".RDS")))
  selected <- which(counts_object$spots$slide == array)
  stopifnot(
    length(selected) == nrow(annotation_object$annots),
    identical(rownames(counts_object$cnts)[selected], rownames(annotation_object$annots))
  )
  counts <- counts_object$cnts[selected, , drop = FALSE]
  spots <- counts_object$spots[selected, , drop = FALSE]
  annotations <- annotation_object$annots
  library_size <- rowSums(counts)
  detected_genes <- rowSums(counts > 0)
  dominant <- colnames(annotations)[max.col(annotations, ties.method = "first")]

  result <- data.table(
    patient = patient,
    array = array,
    barcode = rownames(counts),
    x = spots$new_x,
    y = spots$new_y,
    pixel_x = spots$pixel_x,
    pixel_y = spots$pixel_y,
    nCount_RNA = library_size,
    nFeature_RNA = detected_genes,
    Classification = dominant,
    artifact_fraction = annotations[, "Artefacts"] / pmax(rowSums(annotations), 1)
  )
  for (module in names(modules)) {
    present <- intersect(modules[[module]], colnames(counts))
    if (length(present) == 0) stop("No genes present for module ", module)
    normalized <- log1p(counts[, present, drop = FALSE] / library_size * 10000)
    result[[module]] <- rowMeans(normalized)
    coverage[[length(coverage) + 1]] <- data.table(
      patient = patient, module = module, requested_genes = length(modules[[module]]),
      present_genes = length(present), coverage = length(present) / length(modules[[module]]),
      genes_present = paste(present, collapse = ";")
    )
  }
  for (gene in unlist(lock$primary_program$genes)) {
    if (gene %in% colnames(counts)) {
      result[[paste0("P4gene_", gene)]] <- log1p(counts[, gene] / library_size * 10000)
    }
  }
  fwrite(result, file.path(output_dir, paste0("TNBC", patient, "_annotated_scores.tsv")), sep = "\t")
  audit[[length(audit) + 1]] <- data.table(
    patient = patient, array = array, spots = nrow(result),
    retained_minimum_genes = sum(result$nFeature_RNA >= spec$spot_filter$minimum_detected_genes & result$Classification != "Artefacts"),
    pathology_classes = uniqueN(result$Classification),
    median_umi = median(result$nCount_RNA), median_features = median(result$nFeature_RNA)
  )
  rm(counts_object, annotation_object, counts, spots, annotations, result)
  gc()
}

fwrite(rbindlist(coverage), file.path(table_dir, "TNBC94_spatial_signature_coverage.tsv"), sep = "\t")
fwrite(rbindlist(audit), file.path(table_dir, "TNBC94_annotated_spatial_array_audit.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(log_dir, "23_score_tnbc94_annotated_spatial_arrays_sessionInfo.log"))
