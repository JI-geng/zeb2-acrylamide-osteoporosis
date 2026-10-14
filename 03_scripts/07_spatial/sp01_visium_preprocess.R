#!/usr/bin/env Rscript
# =============================================================================
# Step 7a — spatial (Visium) loading, QC and normalisation
# Dataset : GSE284089
# Packages: Seurat 5.5.1 (Load10X_Spatial / VisiumV1 reader)
#
# Parameters actually used in the reported analysis (from the module parameter
# record of the original run):
#   platform            = Visium (not an imaging-based platform)
#   min_counts          = 0
#   max_counts          = 10,000,000
#   min_features        = 200
#   max_features        = 1,000,000
#   max_pct_mt          = 30
#   remove_empty_spots  = TRUE   (empty_umi_threshold = 100)
#   normalize           = TRUE   (log-normalisation)
#   sctransform         = FALSE
#   seed                = 2024
#
# INPUT : Visium output folder (filtered_feature_bc_matrix + spatial/, or
#         the RDS written by this script). NOT deposited: download from GEO.
# OUTPUT: spatial_io_visium_qc.rds, spatial_empty_spots.csv
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
})

SPATIAL_DIR <- Sys.getenv("SPATIAL_DIR", unset = "GSE284089_visium")
OUT_DIR     <- Sys.getenv("SPATIAL_OUT_DIR", unset = ".")
SEED        <- 2024
set.seed(SEED)

MIN_COUNTS   <- 0
MAX_COUNTS   <- 1e7
MIN_FEATURES <- 200
MAX_FEATURES <- 1e6
MAX_PCT_MT   <- 30
EMPTY_UMI_THRESHOLD <- 100

if (!dir.exists(SPATIAL_DIR)) {
  stop("Visium folder not found: ", SPATIAL_DIR,
       ". Download GSE284089 from GEO and point SPATIAL_DIR at it.")
}

sobj <- Load10X_Spatial(data.dir = SPATIAL_DIR,
                        filename = "filtered_feature_bc_matrix.h5",
                        assay = "Spatial", slice = "sample1")
message("Loaded: ", ncol(sobj), " spots x ", nrow(sobj), " genes")

sobj <- PercentageFeatureSet(sobj, pattern = "^MT-", col.name = "percent.mt")

## Empty spots: UMI count below the threshold
umi <- sobj$nCount_Spatial
empty <- names(umi)[umi < EMPTY_UMI_THRESHOLD]
write.csv(data.frame(spot = empty, nCount = umi[empty]),
          file.path(OUT_DIR, "spatial_empty_spots.csv"), row.names = FALSE)
message("Empty spots (UMI < ", EMPTY_UMI_THRESHOLD, "): ", length(empty))

keep <- sobj$nCount_Spatial  >= MIN_COUNTS   & sobj$nCount_Spatial  <= MAX_COUNTS &
        sobj$nFeature_Spatial >= MIN_FEATURES & sobj$nFeature_Spatial <= MAX_FEATURES &
        sobj$percent.mt <= MAX_PCT_MT
if (TRUE) keep <- keep & !(colnames(sobj) %in% empty)   # remove_empty_spots = TRUE

n_before <- ncol(sobj)
sobj <- subset(sobj, cells = colnames(sobj)[keep])
message("Spatial QC: ", n_before, " -> ", ncol(sobj), " spots")

if (TRUE) {                       # normalize = TRUE, sctransform = FALSE
  sobj <- NormalizeData(sobj, assay = "Spatial", normalization.method = "LogNormalize")
}

saveRDS(sobj, file.path(OUT_DIR, "spatial_io_visium_qc.rds"))
write.csv(sobj[[]], file.path(OUT_DIR, "spatial_spot_metadata.csv"))
message("Done. Object written to ", normalizePath(OUT_DIR))
