#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(NMF)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260721)

root <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[1]), "..", ".."), mustWork = FALSE)
if (!file.exists(file.path(root, "README.md"))) root <- normalizePath(".")
input_file <- file.path(root, "data", "interim", "GSE176078_macrophages_seurat.rds")
output_dir <- file.path(root, "data", "interim")
table_dir <- file.path(root, "results", "tables")
figure_dir <- file.path(root, "results", "figures")

obj <- readRDS(input_file)
features <- head(VariableFeatures(obj), 500)

# The pilot deliberately balances patients so that one macrophage-rich tumor
# cannot dominate the feasibility decomposition. The formal analysis will use
# all cells with substantially more random starts and independent replication.
cells_by_patient <- split(colnames(obj), obj$orig.ident)
pilot_cells <- unlist(lapply(cells_by_patient, function(ids) {
  sample(ids, size = min(length(ids), 100), replace = FALSE)
}), use.names = FALSE)
obj <- subset(obj, cells = pilot_cells)

get_data <- function(object, genes) {
  assay <- DefaultAssay(object)
  if (packageVersion("SeuratObject") >= "5.0.0") {
    return(GetAssayData(object, assay = assay, layer = "data")[genes, , drop = FALSE])
  }
  GetAssayData(object, assay = assay, slot = "data")[genes, , drop = FALSE]
}

message("Constructing non-negative matrix: ", length(features), " genes x ", ncol(obj), " cells")
x <- as.matrix(get_data(obj, features))
x <- x[rowSums(x > 0) >= ceiling(0.05 * ncol(x)), , drop = FALSE]

ranks <- 5:7
message("Estimating pilot NMF ranks: ", paste(ranks, collapse = ", "))
rank_estimate <- nmfEstimateRank(
  x,
  range = ranks,
  method = "brunet",
  nrun = 2,
  seed = 20260721,
  .options = "v-p"
)

measures <- as.data.frame(summary(rank_estimate))
measures$rank <- as.integer(rownames(measures))
fwrite(measures, file.path(table_dir, "GSE176078_macrophage_nmf_rank_metrics.tsv"), sep = "\t")

valid <- measures[is.finite(measures$cophenetic), , drop = FALSE]
if (nrow(valid) == 0) stop("No valid NMF ranks were produced")
best_rank <- valid$rank[which.max(valid$cophenetic)]
message("Pilot-selected rank by maximum cophenetic coefficient: ", best_rank)

fit <- rank_estimate$fit[[as.character(best_rank)]]
w <- basis(fit)
h <- coef(fit)

top_n <- 50
program_genes <- rbindlist(lapply(seq_len(ncol(w)), function(k) {
  ordering <- order(w[, k], decreasing = TRUE)
  data.table(
    program = paste0("P", k),
    rank = seq_len(min(top_n, length(ordering))),
    gene = rownames(w)[ordering][seq_len(min(top_n, length(ordering)))],
    loading = w[ordering, k][seq_len(min(top_n, length(ordering)))]
  )
}))
fwrite(program_genes, file.path(table_dir, "GSE176078_macrophage_nmf_program_genes_pilot.tsv"), sep = "\t")

cell_scores <- as.data.table(t(h), keep.rownames = "cell_id")
setnames(cell_scores, old = setdiff(names(cell_scores), "cell_id"), new = paste0("P", seq_len(nrow(h))))
meta <- as.data.table(obj@meta.data, keep.rownames = "cell_id")
cell_scores <- merge(cell_scores, meta[, .(cell_id, patient = orig.ident, subtype)], by = "cell_id", all.x = TRUE)

program_columns <- grep("^P[0-9]+$", names(cell_scores), value = TRUE)
patient_scores <- cell_scores[, lapply(.SD, median), by = .(patient, subtype), .SDcols = program_columns]
fwrite(patient_scores, file.path(table_dir, "GSE176078_macrophage_nmf_patient_scores_pilot.tsv"), sep = "\t")

long_scores <- melt(patient_scores, id.vars = c("patient", "subtype"), variable.name = "program", value.name = "score")
p_scores <- ggplot(long_scores, aes(x = subtype, y = score, color = subtype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  geom_jitter(width = 0.15, height = 0, size = 1.4) +
  facet_wrap(~program, scales = "free_y") +
  theme_classic(base_size = 9) +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Median NMF exposure", title = "Pilot patient-level macrophage program exposures")
ggsave(file.path(figure_dir, "GSE176078_macrophage_nmf_patient_scores_pilot.pdf"), p_scores, width = 8, height = 5.5)

saveRDS(
  list(rank_estimate = rank_estimate, best_rank = best_rank, fit = fit, features = rownames(x)),
  file.path(output_dir, "GSE176078_macrophage_nmf_pilot.rds"),
  compress = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(root, "logs", "02_nmf_macrophage_pilot_sessionInfo.log"))
message("Pilot NMF complete")
