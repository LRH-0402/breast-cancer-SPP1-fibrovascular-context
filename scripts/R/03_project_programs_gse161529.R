#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260721)
root <- normalizePath(".")
object_file <- file.path(root, "data", "interim", "GSE161529_TNBCSub_updated.rds")
program_file <- file.path(root, "results", "tables", "GSE176078_macrophage_nmf_program_genes_pilot.tsv")
table_dir <- file.path(root, "results", "tables")
figure_dir <- file.path(root, "results", "figures")

obj <- readRDS(object_file)
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

programs <- fread(program_file)
programs <- programs[rank <= 30]
gene_sets <- split(programs$gene, programs$program)

get_data <- function(object) {
  if (packageVersion("SeuratObject") >= "5.0.0") {
    return(GetAssayData(object, assay = "RNA", layer = "data"))
  }
  GetAssayData(object, assay = "RNA", slot = "data")
}

expr <- get_data(obj)
for (program in names(gene_sets)) {
  genes <- intersect(gene_sets[[program]], rownames(expr))
  obj[[paste0(program, "_projection")]] <- Matrix::colMeans(expr[genes, , drop = FALSE])
}

score_columns <- paste0(names(gene_sets), "_projection")
cell_table <- as.data.table(obj@meta.data, keep.rownames = "cell_id")
cluster_scores <- cell_table[, c(
  list(cells = .N),
  lapply(.SD, median)
), by = seurat_clusters, .SDcols = score_columns]
fwrite(cluster_scores, file.path(table_dir, "GSE161529_TNBCSub_program_projection_by_cluster.tsv"), sep = "\t")

# The author code identifies clusters 0 and 4 as T-cell populations. Marker
# expression independently identifies cluster 1 as the macrophage-rich cluster.
macrophage_cluster <- "1"
mac <- cell_table[as.character(seurat_clusters) == macrophage_cluster]
patient_scores <- mac[, c(
  list(macrophage_cells = .N),
  lapply(.SD, median)
), by = group, .SDcols = score_columns]
fwrite(patient_scores, file.path(table_dir, "GSE161529_TNBCSub_macrophage_program_projection_by_patient.tsv"), sep = "\t")

long_cluster <- melt(
  cluster_scores,
  id.vars = c("seurat_clusters", "cells"),
  measure.vars = score_columns,
  variable.name = "program",
  value.name = "median_score"
)
long_cluster[, program := sub("_projection$", "", program)]
p_cluster <- ggplot(long_cluster, aes(x = factor(seurat_clusters), y = median_score, fill = program)) +
  geom_col(position = "dodge") +
  theme_classic(base_size = 9) +
  labs(x = "Author-defined microenvironment cluster", y = "Median projected score",
       title = "Independent TNBC cohort: discovery-program projection")

long_patient <- melt(
  patient_scores,
  id.vars = c("group", "macrophage_cells"),
  measure.vars = score_columns,
  variable.name = "program",
  value.name = "median_score"
)
long_patient[, program := sub("_projection$", "", program)]
p_patient <- ggplot(long_patient, aes(x = reorder(group, median_score), y = median_score, color = program)) +
  geom_point(size = 2) +
  facet_wrap(~program, scales = "free_y") +
  coord_flip() +
  theme_classic(base_size = 9) +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Median macrophage score", title = "Patient-level program heterogeneity")

ggsave(
  file.path(figure_dir, "GSE161529_TNBCSub_program_projection_pilot.pdf"),
  p_cluster / p_patient,
  width = 9,
  height = 9
)

saveRDS(obj, file.path(root, "data", "interim", "GSE161529_TNBCSub_scored_pilot.rds"), compress = FALSE)
writeLines(capture.output(sessionInfo()), file.path(root, "logs", "03_project_programs_gse161529_sessionInfo.log"))

message("Projected ", length(gene_sets), " discovery programs into ", ncol(obj), " validation cells")
message("Macrophage cluster contains ", nrow(mac), " cells from ", uniqueN(mac$group), " samples")

