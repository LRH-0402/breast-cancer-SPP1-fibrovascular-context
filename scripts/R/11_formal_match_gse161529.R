#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
root <- normalizePath(".")
fit <- readRDS(file.path(root, "data", "interim", "GSE161529_macrophage_independent_nmf.rds"))
universe <- fit$features
discovery <- fread(file.path(root, "results", "tables", "GSE176078_macrophage_nmf_program_genes_pilot.tsv"))[rank <= 50]
validation <- fread(file.path(root, "results", "tables", "GSE161529_macrophage_independent_nmf_program_genes.tsv"))

comparison <- rbindlist(lapply(unique(discovery$program), function(dp) {
  rbindlist(lapply(unique(validation$validation_program), function(vp) {
    dg <- intersect(discovery[program == dp]$gene, universe)
    vg <- validation[validation_program == vp]$gene
    overlap <- intersect(dg, vg)
    data.table(
      discovery_program = dp,
      validation_program = vp,
      overlap_n = length(overlap),
      universe_genes = length(universe),
      discovery_genes_in_universe = length(dg),
      hypergeom_p = phyper(
        length(overlap) - 1,
        length(dg),
        length(universe) - length(dg),
        length(vg),
        lower.tail = FALSE
      ),
      overlapping_genes = paste(overlap, collapse = ";")
    )
  }))
}))
comparison[, hypergeom_fdr := p.adjust(hypergeom_p, method = "BH")]
setorder(comparison, discovery_program, hypergeom_fdr, -overlap_n)
best <- comparison[, .SD[1], by = discovery_program]

fwrite(comparison, file.path(root, "results", "tables", "GSE161529_formal_program_matching.tsv"), sep = "\t")
fwrite(best, file.path(root, "results", "tables", "GSE161529_formal_program_best_matches.tsv"), sep = "\t")
print(best)
