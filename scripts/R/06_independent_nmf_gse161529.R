#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(NMF)
  library(data.table)
  library(ggplot2)
})

set.seed(20260721)
root <- normalizePath(".")
object_file <- file.path(root, "data", "interim", "GSE161529_TNBCSub_updated.rds")
discovery_file <- file.path(root, "results", "tables", "GSE176078_macrophage_nmf_program_genes_pilot.tsv")
table_dir <- file.path(root, "results", "tables")
figure_dir <- file.path(root, "results", "figures")
log_dir <- file.path(root, "logs")

obj <- readRDS(object_file)
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

# Author metadata plus marker review identifies zero-based cluster 1 as the
# macrophage-rich cluster. Sampling is balanced by tumor before factorization.
mac <- subset(obj, subset = seurat_clusters == 1)
mac <- FindVariableFeatures(mac, selection.method = "vst", nfeatures = 1000, verbose = FALSE)
cells_by_patient <- split(colnames(mac), mac$group)
balanced_cells <- unlist(lapply(cells_by_patient, function(ids) {
  sample(ids, size = min(length(ids), 250), replace = FALSE)
}), use.names = FALSE)
mac <- subset(mac, cells = balanced_cells)

get_data <- function(object, genes) {
  if (packageVersion("SeuratObject") >= "5.0.0") {
    return(GetAssayData(object, assay = "RNA", layer = "data")[genes, , drop = FALSE])
  }
  GetAssayData(object, assay = "RNA", slot = "data")[genes, , drop = FALSE]
}

features <- head(VariableFeatures(mac), 500)
x <- as.matrix(get_data(mac, features))
x <- x[rowSums(x > 0) >= ceiling(0.05 * ncol(x)), , drop = FALSE]

ranks <- 5:7
message("Independent TNBC NMF: ", nrow(x), " genes x ", ncol(x), " balanced macrophages")
rank_estimate <- nmfEstimateRank(
  x,
  range = ranks,
  method = "brunet",
  nrun = 3,
  seed = 20260721,
  .options = "v-p"
)
metrics <- as.data.frame(summary(rank_estimate))
metrics$rank <- as.integer(rownames(metrics))
fwrite(metrics, file.path(table_dir, "GSE161529_macrophage_independent_nmf_rank_metrics.tsv"), sep = "\t")

valid <- metrics[is.finite(metrics$cophenetic), , drop = FALSE]
if (nrow(valid) == 0) stop("No valid validation NMF ranks")
best_rank <- valid$rank[which.max(valid$cophenetic)]
fit <- rank_estimate$fit[[as.character(best_rank)]]
w <- basis(fit)

top_n <- 50
validation_genes <- rbindlist(lapply(seq_len(ncol(w)), function(k) {
  ordering <- order(w[, k], decreasing = TRUE)
  data.table(
    validation_program = paste0("V", k),
    rank = seq_len(min(top_n, length(ordering))),
    gene = rownames(w)[ordering][seq_len(min(top_n, length(ordering)))],
    loading = w[ordering, k][seq_len(min(top_n, length(ordering)))]
  )
}))
fwrite(validation_genes, file.path(table_dir, "GSE161529_macrophage_independent_nmf_program_genes.tsv"), sep = "\t")

discovery <- fread(discovery_file)[rank <= top_n]
comparison <- rbindlist(lapply(unique(discovery$program), function(dp) {
  dg <- discovery[program == dp, unique(gene)]
  rbindlist(lapply(unique(validation_genes$validation_program), function(vp) {
    vg <- validation_genes[validation_program == vp, unique(gene)]
    overlap <- intersect(dg, vg)
    data.table(
      discovery_program = dp,
      validation_program = vp,
      overlap_n = length(overlap),
      jaccard = length(overlap) / length(union(dg, vg)),
      overlap_coefficient = length(overlap) / min(length(dg), length(vg)),
      overlapping_genes = paste(overlap, collapse = ";")
    )
  }))
}))
fwrite(comparison, file.path(table_dir, "GSE161529_vs_GSE176078_nmf_program_overlap.tsv"), sep = "\t")

best_matches <- comparison[order(discovery_program, -overlap_coefficient), .SD[1], by = discovery_program]
fwrite(best_matches, file.path(table_dir, "GSE161529_vs_GSE176078_nmf_best_matches.tsv"), sep = "\t")

p <- ggplot(comparison, aes(x = validation_program, y = discovery_program, fill = overlap_coefficient)) +
  geom_tile(color = "white") +
  geom_text(aes(label = overlap_n), size = 3) +
  scale_fill_viridis_c(limits = c(0, max(comparison$overlap_coefficient))) +
  theme_classic(base_size = 9) +
  labs(x = "Independent TNBC factor", y = "Discovery factor", fill = "Overlap\ncoefficient",
       title = "Independent factorization: top-50 gene overlap")
ggsave(file.path(figure_dir, "GSE161529_vs_GSE176078_nmf_program_overlap.pdf"), p, width = 6.5, height = 4.8)

saveRDS(
  list(rank_estimate = rank_estimate, best_rank = best_rank, fit = fit, features = rownames(x)),
  file.path(root, "data", "interim", "GSE161529_macrophage_independent_nmf.rds"),
  compress = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "06_independent_nmf_gse161529_sessionInfo.log"))
message("Selected rank ", best_rank)
print(best_matches)
