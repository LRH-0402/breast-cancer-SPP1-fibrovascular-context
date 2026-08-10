#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

root <- normalizePath(".")
table_dir <- file.path(root, "results", "tables")
log_dir <- file.path(root, "logs")

get_counts <- function(object) {
  assay <- if ("Spatial" %in% Assays(object)) "Spatial" else DefaultAssay(object)
  if (packageVersion("SeuratObject") >= "5.0.0") {
    return(GetAssayData(object, assay = assay, layer = "counts"))
  }
  GetAssayData(object, assay = assay, slot = "counts")
}

figshare_dir <- file.path(root, "data", "raw", "spatial", "figshare_21591429")
figshare_rows <- rbindlist(lapply(sprintf("P%d", 1:8), function(section) {
  path <- file.path(figshare_dir, paste0(section, ".rds"))
  object <- readRDS(path)
  counts <- get_counts(object)
  coordinates <- tryCatch(
    GetTissueCoordinates(object, image = names(object@images)[1]),
    error = function(e) object@images[[1]]@coordinates
  )
  result <- data.table(
    cohort = "Figshare_21591429",
    section = section,
    genes = nrow(counts),
    spots = ncol(counts),
    nonzero_counts = length(counts@x),
    median_umi = median(object$nCount_Spatial),
    median_features = median(object$nFeature_Spatial),
    images = length(object@images),
    coordinate_rows = nrow(coordinates),
    clusters = uniqueN(object$seurat_clusters),
    file_md5 = unname(tools::md5sum(path))
  )
  rm(object, counts, coordinates)
  gc()
  result
}))
fwrite(figshare_rows, file.path(table_dir, "Figshare_21591429_spatial_audit.tsv"), sep = "\t")

wu_matrix_root <- file.path(
  root, "data", "interim", "spatial_audit", "wu_visium", "filtered_count_matrices"
)
wu_metadata_root <- file.path(root, "data", "interim", "spatial_audit", "wu_metadata", "metadata")
wu_rows <- rbindlist(lapply(c("1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"), function(patient) {
  matrix_dir <- file.path(wu_matrix_root, paste0(patient, "_filtered_count_matrix"))
  counts <- Matrix::readMM(file.path(matrix_dir, "matrix.mtx.gz"))
  genes <- readLines(file.path(matrix_dir, "features.tsv.gz"))
  barcodes <- readLines(file.path(matrix_dir, "barcodes.tsv.gz"))
  rownames(counts) <- genes
  colnames(counts) <- barcodes
  metadata <- fread(file.path(wu_metadata_root, paste0(patient, "_metadata.csv")))
  barcode_column <- names(metadata)[1]
  matched <- sum(colnames(counts) %in% metadata[[barcode_column]])
  data.table(
    cohort = "Zenodo_4739739",
    patient = patient,
    subtype = unique(metadata$subtype),
    genes = nrow(counts),
    spots = ncol(counts),
    metadata_spots = nrow(metadata),
    matched_barcodes = matched,
    median_umi = median(metadata$nCount_RNA),
    median_features = median(metadata$nFeature_RNA),
    pathology_classes = uniqueN(metadata$Classification)
  )
}))
fwrite(wu_rows, file.path(table_dir, "Wu_Visium_spatial_audit.tsv"), sep = "\t")

clinical <- readRDS(file.path(
  root, "data", "interim", "spatial_audit", "tnbc_clinical", "Clinical", "Clinical.RDS"
))
ids <- readRDS(file.path(
  root, "data", "interim", "spatial_audit", "tnbc_clinical", "Clinical", "ids.RDS"
))
clinical_summary <- data.table(
  metric = c(
    "patients", "arrays_or_subarrays", "patients_with_arrays", "annotated_arrays",
    "distant_relapse_events", "idfs_events", "os_events", "median_followup_years"
  ),
  value = c(
    nrow(clinical), nrow(ids), uniqueN(ids$id), sum(ids$hasAnnot, na.rm = TRUE),
    sum(clinical$DRFS[, "status"], na.rm = TRUE),
    sum(clinical$iDFS[, "status"], na.rm = TRUE),
    sum(clinical$OS[, "status"], na.rm = TRUE),
    median(clinical$OS[, "time"], na.rm = TRUE)
  )
)
fwrite(clinical_summary, file.path(table_dir, "TNBC_Zenodo_14204217_spatial_clinical_audit.tsv"), sep = "\t")

writeLines(capture.output(sessionInfo()), file.path(log_dir, "16_audit_spatial_cohorts_sessionInfo.log"))
print(figshare_rows)
print(wu_rows)
print(clinical_summary)
