#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(data.table)
  library(ggplot2)
})

set.seed(20260721)

root <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[1]), "..", ".."), mustWork = FALSE)
if (!file.exists(file.path(root, "README.md"))) root <- normalizePath(".")

input_dir <- file.path(root, "data", "raw", "Wu_etal_2021_BRCA_scRNASeq")
output_dir <- file.path(root, "data", "interim")
table_dir <- file.path(root, "results", "tables")
figure_dir <- file.path(root, "results", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading sparse count matrix")
counts <- readMM(file.path(input_dir, "count_matrix_sparse.mtx"))
genes <- fread(file.path(input_dir, "count_matrix_genes.tsv"), header = FALSE)[[1]]
barcodes <- fread(file.path(input_dir, "count_matrix_barcodes.tsv"), header = FALSE)[[1]]
metadata <- fread(file.path(input_dir, "metadata.csv"), data.table = FALSE)
rownames(metadata) <- metadata[[1]]
metadata[[1]] <- NULL

stopifnot(nrow(counts) == length(genes), ncol(counts) == length(barcodes))
stopifnot(setequal(barcodes, rownames(metadata)))
rownames(counts) <- make.unique(genes)
colnames(counts) <- barcodes
metadata <- metadata[barcodes, , drop = FALSE]

keep <- metadata$celltype_minor == "Macrophage"
mac_counts <- counts[, keep, drop = FALSE]
mac_meta <- metadata[keep, , drop = FALSE]
rm(counts)
invisible(gc())

message("Creating macrophage Seurat object")
obj <- CreateSeuratObject(counts = mac_counts, meta.data = mac_meta, project = "GSE176078_macrophages")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 1500, verbose = FALSE)

excluded_pattern <- "^(MT-|RPS|RPL|IG[HKL]|TR[ABDG])"
clean_variable_features <- VariableFeatures(obj)[!grepl(excluded_pattern, VariableFeatures(obj), ignore.case = FALSE)]
VariableFeatures(obj) <- clean_variable_features

qc_patient <- as.data.table(obj@meta.data)[, .(
  cells = .N,
  median_umis = as.numeric(median(nCount_RNA)),
  median_features = as.numeric(median(nFeature_RNA)),
  median_percent_mito = as.numeric(median(percent.mito)),
  subtype = as.character(unique(subtype)[1])
), by = orig.ident]
setorder(qc_patient, subtype, -cells)
fwrite(qc_patient, file.path(table_dir, "GSE176078_macrophage_patient_qc.tsv"), sep = "\t")

p1 <- ggplot(qc_patient, aes(x = reorder(orig.ident, cells), y = cells, fill = subtype)) +
  geom_col(width = 0.8) +
  coord_flip() +
  labs(x = NULL, y = "Macrophage cells", title = "GSE176078 macrophage representation by patient") +
  theme_classic(base_size = 10)
ggsave(file.path(figure_dir, "GSE176078_macrophage_cells_by_patient.pdf"), p1, width = 6.5, height = 6)

saveRDS(obj, file.path(output_dir, "GSE176078_macrophages_seurat.rds"), compress = FALSE)

provenance <- data.frame(
  dataset = "GSE176078",
  source_cells = nrow(metadata),
  macrophage_cells = ncol(obj),
  patients = length(unique(obj$orig.ident)),
  variable_features = length(VariableFeatures(obj)),
  stringsAsFactors = FALSE
)
fwrite(provenance, file.path(table_dir, "GSE176078_preparation_summary.tsv"), sep = "\t")

writeLines(capture.output(sessionInfo()), file.path(root, "logs", "01_prepare_gse176078_sessionInfo.log"))
message("Prepared object: ", ncol(obj), " macrophages from ", length(unique(obj$orig.ident)), " patients")
