#!/usr/bin/env Rscript
# =============================================================================
# Step 7b — RCTD deconvolution of the Visium spots with the single-cell
#           reference (GSE147287), used to localise the ZEB2-high population
# Packages: spacexr 2.2.1
#
# Parameters / dimensions as used in the reported analysis:
#   * reference: 11,562 cells, 20,955 genes
#   * doublet mode: the spacexr default ("doublet")
#   * platform effect normalisation enabled (the step reports the platform
#     effect normalisation DE genes being computed)
#
# INPUT : spatial_io_visium_qc.rds (step 7a) and the annotated single-cell
#         object from step 6b (its celltype column is the reference label).
# OUTPUT: spatial_rctd.rds, rctd_celltype_proportions.csv (spot x cell type)
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(spacexr)      # 2.2.1
})

SPATIAL_RDS  <- Sys.getenv("SPATIAL_RDS",  unset = "spatial_io_visium_qc.rds")
REFERENCE_RDS <- Sys.getenv("SC_ANNOTATED_RDS", unset = "05_annotated.rds")
OUT_DIR      <- Sys.getenv("SPATIAL_OUT_DIR", unset = ".")
REF_LABEL_COL <- "celltype"
DOUBLET_MODE <- Sys.getenv("RCTD_DOUBLET_MODE", unset = "doublet")

spat <- readRDS(SPATIAL_RDS)
ref  <- readRDS(REFERENCE_RDS)

counts <- GetAssayData(spat, assay = "Spatial", layer = "counts")
coords <- GetTissueCoordinates(spat)
coords <- data.frame(x = coords[, 1], y = coords[, 2],
                     row.names = colnames(spat)[match(rownames(coords), colnames(spat))])
coords <- coords[colnames(counts), , drop = FALSE]

## ---- build the spatial object ---------------------------------------------
puck <- SpatialRNA(coords, counts, colSums(counts))

## ---- build the reference --------------------------------------------------
ref_counts <- GetAssayData(ref, layer = "counts")
cell_types <- ref[[REF_LABEL_COL]][, 1]
names(cell_types) <- colnames(ref)
reference <- Reference(ref_counts, cell_types)
message("Reference: ", ncol(ref_counts), " cells, ", nrow(ref_counts), " genes, ",
        length(unique(cell_types)), " cell types")

## ---- run RCTD -------------------------------------------------------------
rctd <- create.RCTD(puck, reference, max_cores = 1)
rctd <- run.RCTD(rctd, doublet_mode = DOUBLET_MODE)

saveRDS(rctd, file.path(OUT_DIR, "spatial_rctd.rds"))

props <- as.data.frame(as.matrix(rctd@results$weights))
write.csv(cbind(spot = rownames(props), props),
          file.path(OUT_DIR, "rctd_celltype_proportions.csv"), row.names = FALSE)

## ---- ZEB2-high population: mean proportion per spot -----------------------
if (any(grepl("ZEB2", colnames(props)))) {
  zc <- colnames(props)[grepl("ZEB2", colnames(props))]
  write.csv(data.frame(spot = rownames(props),
                       ZEB2_high = rowSums(props[, zc, drop = FALSE])),
            file.path(OUT_DIR, "rctd_ZEB2high_per_spot.csv"), row.names = FALSE)
  message("Wrote per-spot ZEB2-high proportion (", paste(zc, collapse = ", "), ")")
}
message("Done. RCTD output in ", normalizePath(OUT_DIR))
