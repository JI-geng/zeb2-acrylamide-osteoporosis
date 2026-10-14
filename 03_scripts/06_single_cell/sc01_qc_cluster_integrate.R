#!/usr/bin/env Rscript
# =============================================================================
# Step 6a — single-cell QC, batch integration, clustering and embedding
# Dataset : GSE147287 (bone-marrow / skeletal single-cell RNA-seq)
# Packages: Seurat 5.5.1, scDblFinder, harmony
#
# Every numeric value below is the value actually used in the analysis reported
# in the manuscript (taken from the parameter records of the original analysis).
# Nothing here is tuned post hoc: if you change a value, you are no longer
# reproducing the published clustering.
#
# INPUT  (not deposited — see NOTE ON INPUTS below):
#   SEURAT_IN_RDS  a Seurat object built from the GSE147287 matrix
#                  (or the CellRanger filtered_feature_bc_matrix folders)
# OUTPUT:
#   03b_integrated.rds, 04_clustered.rds, umap/tsne coordinates (CSV)
#
# NOTE ON INPUTS: this script needs the GSE147287 expression matrix, which is a
# public GEO download and is NOT distributed with this repository (see README
# section 1). Point SEURAT_IN_RDS at your local copy before running.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)      # 5.5.1
  library(scDblFinder)
  library(harmony)
})

SEURAT_IN_RDS  <- Sys.getenv("SEURAT_IN_RDS",  unset = "GSE147287_seurat_raw.rds")
OUT_DIR        <- Sys.getenv("SC_OUT_DIR",     unset = ".")
BATCH_COL      <- "sample_id"      # 2 batches in GSE147287
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ---------------------------------------------------------------------------
## 0. Load
## ---------------------------------------------------------------------------
if (file.exists(SEURAT_IN_RDS)) {
  sobj <- readRDS(SEURAT_IN_RDS)
} else {
  # Alternative: build the object from CellRanger output.
  # GSE147287 supplies filtered_feature_bc_matrix per sample; read them with
  # Read10X() and CreateSeuratObject(min.cells = 3, min.features = 200).
  stop("SEURAT_IN_RDS not found (", SEURAT_IN_RDS,
       "). Download GSE147287 from GEO and build a Seurat object first.")
}
message("Loaded: ", ncol(sobj), " cells x ", nrow(sobj), " genes")

## ---------------------------------------------------------------------------
## 1. QC metrics and doublet scoring
##    scDblFinder is run BEFORE filtering and scored only; removal happens
##    together with the QC thresholds in step 2 (18,187 -> 15,481 cells).
## ---------------------------------------------------------------------------
sobj <- PercentageFeatureSet(sobj, pattern = "^MT-",  col.name = "percent.mt")
sobj <- PercentageFeatureSet(sobj, pattern = "^RP[SL]", col.name = "percent.ribo")
sobj <- PercentageFeatureSet(sobj, pattern = "^HB[^(P)]", col.name = "percent.hb")

sobj <- scDblFinder(sobj)          # scores only at this stage
message("scDblFinder done; doublet calls: ",
        sum(sobj$scDblFinder.class == "doublet"))

## ---------------------------------------------------------------------------
## 2. Threshold filtering — values used in the reported analysis
## ---------------------------------------------------------------------------
MIN_FEATURES <- 200
MAX_FEATURES <- 1e6
MIN_COUNTS   <- 0
MAX_COUNTS   <- 1e7
MAX_PCT_MT   <- 20      # percent.mt <= 20
MAX_PCT_RIBO <- 100     # effectively no filter
MAX_PCT_HB   <- 100     # effectively no filter

n_before <- ncol(sobj)
keep <- sobj$nFeature_RNA >= MIN_FEATURES & sobj$nFeature_RNA <= MAX_FEATURES &
        sobj$nCount_RNA  >= MIN_COUNTS   & sobj$nCount_RNA  <= MAX_COUNTS   &
        sobj$percent.mt  <= MAX_PCT_MT   &
        sobj$percent.ribo <= MAX_PCT_RIBO &
        sobj$percent.hb   <= MAX_PCT_HB
sobj <- subset(sobj, cells = colnames(sobj)[keep])
message("QC filter: ", n_before, " -> ", ncol(sobj), " cells")

saveRDS(sobj, file.path(OUT_DIR, "02_qc_filtered.rds"))

## ---------------------------------------------------------------------------
## 3. Normalisation, HVG selection, scaling, PCA
##    2000 HVG | ScaleData with no regression | npcs = 30
## ---------------------------------------------------------------------------
sobj <- NormalizeData(sobj, normalization.method = "LogNormalize",
                      scale.factor = 1e4)
sobj <- FindVariableFeatures(sobj, nfeatures = 2000)
sobj <- ScaleData(sobj, vars.to.regress = NULL)     # regress = none
sobj <- RunPCA(sobj, npcs = 30)
saveRDS(sobj, file.path(OUT_DIR, "03_pca_reduced.rds"))

## ---------------------------------------------------------------------------
## 4. Batch integration — Harmony over sample_id, dims 1:18
##    ✅ This step is described in the Methods ("batch effects between donors
##    were corrected using Harmony"). GSE147287 has two donors = two batches,
##    so the integration is required to reproduce the published clustering.
## ---------------------------------------------------------------------------
sobj <- RunHarmony(sobj, group.by.vars = BATCH_COL, reduction.use = "pca",
                   dims.use = 1:18)
saveRDS(sobj, file.path(OUT_DIR, "03b_integrated.rds"))

## ---------------------------------------------------------------------------
## 5. Neighbour graph, clustering, embeddings
##    FindNeighbors dims = 1:18 on the harmony reduction
##    FindClusters  resolution = 0.5  ->  17 clusters
##    UMAP / tSNE   dims = 1:18 on harmony (Seurat 5 uses uwot, metric cosine)
## ---------------------------------------------------------------------------
sobj <- FindNeighbors(sobj, reduction = "harmony", dims = 1:18)
sobj <- FindClusters(sobj, resolution = 0.5)
message("Clusters: ", length(unique(sobj$seurat_clusters)))
sobj <- RunUMAP(sobj, reduction = "harmony", dims = 1:18,
                umap.method = "uwot", metric = "cosine")
sobj <- RunTSNE(sobj, reduction = "harmony", dims = 1:18)

saveRDS(sobj, file.path(OUT_DIR, "04_clustered.rds"))
write.csv(Embeddings(sobj, "umap"), file.path(OUT_DIR, "umap_coordinates.csv"))
write.csv(Embeddings(sobj, "tsne"), file.path(OUT_DIR, "tsne_coordinates.csv"))
write.csv(sobj[[]], file.path(OUT_DIR, "cell_metadata.csv"))
message("Done. Objects written to ", normalizePath(OUT_DIR))
