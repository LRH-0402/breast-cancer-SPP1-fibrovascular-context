#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(yaml)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grep("^--file=", args)])
root <- normalizePath(file.path(dirname(script_arg), "..", ".."))
input_dir <- file.path(root, "data", "raw", "Wu_etal_2021_BRCA_scRNASeq")
table_dir <- file.path(root, "results", "tables")
log_dir <- file.path(root, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

lock <- read_yaml(file.path(root, "config", "locked_macrophage_programs_v1.yml"))
p4 <- unlist(lock$primary_program$genes)
counts <- readMM(file.path(input_dir, "count_matrix_sparse.mtx"))
genes <- fread(file.path(input_dir, "count_matrix_genes.tsv"), header = FALSE)[[1]]
barcodes <- fread(file.path(input_dir, "count_matrix_barcodes.tsv"), header = FALSE)[[1]]
meta <- fread(file.path(input_dir, "metadata.csv"))
barcode_col <- names(meta)[1]
setnames(meta, barcode_col, "barcode")
stopifnot(nrow(counts) == length(genes), ncol(counts) == length(barcodes))
stopifnot(identical(meta$barcode, barcodes))

gene_index <- match(p4, genes)
if (anyNA(gene_index)) stop("Missing P4 genes: ", paste(p4[is.na(gene_index)], collapse = ", "))
library_size <- Matrix::colSums(counts)
selected <- counts[gene_index, , drop = FALSE] %*% Diagonal(x = 10000 / pmax(library_size, 1))
selected@x <- log1p(selected@x)
rownames(selected) <- p4
meta[, P4_primary := Matrix::colMeans(selected)]

major_map <- c(
  "Macrophage" = "Macrophage",
  "Myeloid" = "Other myeloid",
  "T-cells" = "Lymphocyte",
  "B-cells" = "Lymphocyte",
  "Endothelial" = "Endothelial",
  "CAFs" = "CAF",
  "Cancer Epithelial" = "Malignant epithelial"
)
meta[, analysis_celltype := unname(major_map[as.character(celltype_major)])]
meta[celltype_minor == "Macrophage", analysis_celltype := "Macrophage"]
meta[is.na(analysis_celltype) & grepl("Mono|DC|Neutro|Mast", celltype_minor, ignore.case = TRUE),
     analysis_celltype := "Other myeloid"]
meta[is.na(analysis_celltype) & grepl("T|B|NK", celltype_major, ignore.case = TRUE),
     analysis_celltype := "Lymphocyte"]
meta[is.na(analysis_celltype) & grepl("Fibro|CAF", celltype_major, ignore.case = TRUE),
     analysis_celltype := "CAF"]
meta[is.na(analysis_celltype) & grepl("Cancer|Epithelial|Malignant", celltype_major, ignore.case = TRUE),
     analysis_celltype := "Malignant epithelial"]

keep_types <- c("Macrophage", "Other myeloid", "CAF", "Endothelial", "Malignant epithelial", "Lymphocyte")
meta <- meta[analysis_celltype %in% keep_types]
patient_summary <- meta[, .(
  cells = .N,
  median_P4 = median(P4_primary),
  mean_P4 = mean(P4_primary),
  fraction_P4_detected = mean(P4_primary > 0)
), by = .(orig.ident, subtype, analysis_celltype)]
celltype_summary <- patient_summary[, .(
  patients = uniqueN(orig.ident),
  cells = sum(cells),
  median_patient_P4 = median(mean_P4),
  q1_patient_P4 = quantile(mean_P4, 0.25),
  q3_patient_P4 = quantile(mean_P4, 0.75)
), by = analysis_celltype]

gene_rows <- list()
for (gene in p4) {
  meta[, value := as.numeric(selected[gene, match(barcode, barcodes)])]
  gene_rows[[gene]] <- meta[, .(
    cells = .N,
    mean_expression = mean(value),
    fraction_detected = mean(value > 0)
  ), by = analysis_celltype][, gene := gene]
}
gene_summary <- rbindlist(gene_rows)
setcolorder(gene_summary, c("gene", "analysis_celltype", "cells", "mean_expression", "fraction_detected"))

fwrite(patient_summary, file.path(table_dir, "GSE176078_P4_celltype_specificity_by_patient.tsv"), sep = "\t")
fwrite(celltype_summary, file.path(table_dir, "GSE176078_P4_celltype_specificity_summary.tsv"), sep = "\t")
fwrite(gene_summary, file.path(table_dir, "GSE176078_P4_gene_expression_by_celltype.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(log_dir, "40_gse176078_celltype_specificity_sessionInfo.log"))
