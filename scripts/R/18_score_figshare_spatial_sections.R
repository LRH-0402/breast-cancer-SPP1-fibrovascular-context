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
input_dir <- file.path(root, "data", "raw", "spatial", "figshare_21591429")
output_dir <- file.path(root, "data", "derived", "figshare_21591429_scores")
table_dir <- file.path(root, "results", "tables")
log_dir <- file.path(root, "logs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

lock <- read_yaml(file.path(root, "config", "locked_macrophage_programs_v1.yml"))
spec <- read_yaml(file.path(root, "config", "spatial_analysis_spec_v1.yml"))
tcell_spec <- read_yaml(file.path(root, "config", "tcell_state_sensitivity_v1.yml"))
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

score_module <- function(counts, library_size, genes) {
  present <- intersect(genes, rownames(counts))
  if (length(present) == 0) stop("No requested genes are present")
  values <- counts[present, , drop = FALSE]
  values <- values %*% Diagonal(x = 10000 / library_size)
  values@x <- log1p(values@x)
  list(score = Matrix::colMeans(values), present = present)
}

coverage <- list()
for (section in paste0("P", 1:8)) {
  message("Scoring ", section)
  object <- readRDS(file.path(input_dir, paste0(section, ".rds")))
  counts <- GetAssayData(object, assay = "Spatial", layer = "counts")
  cells <- colnames(counts)
  metadata <- object@meta.data[cells, , drop = FALSE]
  coordinates <- GetTissueCoordinates(object, image = names(object@images)[1])
  coordinates <- coordinates[cells, , drop = FALSE]
  keep <- metadata$nFeature_Spatial >= spec$spot_filter$minimum_detected_genes
  counts <- counts[, keep, drop = FALSE]
  metadata <- metadata[keep, , drop = FALSE]
  coordinates <- coordinates[keep, , drop = FALSE]
  library_size <- Matrix::colSums(counts)

  result <- data.table(
    barcode = colnames(counts),
    section = section,
    x = coordinates$imagerow,
    y = coordinates$imagecol,
    nCount_RNA = as.numeric(library_size),
    nFeature_RNA = metadata$nFeature_Spatial
  )
  for (module in names(modules)) {
    scored <- score_module(counts, library_size, modules[[module]])
    result[[module]] <- scored$score
    coverage[[length(coverage) + 1]] <- data.table(
      section = section,
      module = module,
      requested_genes = length(modules[[module]]),
      present_genes = length(scored$present),
      coverage = length(scored$present) / length(modules[[module]]),
      genes_present = paste(scored$present, collapse = ";")
    )
  }
  for (gene in unlist(lock$primary_program$genes)) {
    if (gene %in% rownames(counts)) {
      gene_values <- counts[gene, , drop = FALSE] %*% Diagonal(x = 10000 / library_size)
      result[[paste0("P4gene_", gene)]] <- as.numeric(log1p(gene_values))
    }
  }
  fwrite(result, file.path(output_dir, paste0(section, "_scores.tsv")), sep = "\t")
  rm(object, counts, metadata, coordinates, result)
  gc()
}

fwrite(rbindlist(coverage), file.path(table_dir, "Figshare_spatial_signature_coverage.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(log_dir, "18_score_figshare_spatial_sections_sessionInfo.log"))
