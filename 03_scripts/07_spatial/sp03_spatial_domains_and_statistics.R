#!/usr/bin/env Rscript
# =============================================================================
# Step 7c — spatially variable genes, spatial domains and spatial statistics
#           (Methods 2.8)
#
# As described in the manuscript:
#   * spatially variable genes: SPARK-X
#   * spatial domains: spatially aware Louvain clustering on the first 15 PCs,
#     a six-nearest-neighbour graph, resolution 0.5
#   * domain markers: minimum expression fraction 25%, |avg log2FC| > 0.25
#   * spatial graph for the autocorrelation statistics: six nearest neighbours
#   * global autocorrelation: Moran's I and Geary's C, analytical P values from
#     the normal approximation, Benjamini-Hochberg adjusted
#   * local indicators of spatial association (LISA): 100 permutations
#   * local clustering: Getis-Ord Gi*
#   * genes examined: ZEB2, CXCL12, CXCR4, COL1A1, SPARC, BGLAP
#
# Reported values to check against (Results 3.7):
#   QC:  1,616 of 4,992 spots retained
#   Moran's I: CXCL12 0.063, CXCR4 0.083, COL1A1 0.367, SPARC 0.222,
#              BGLAP 0.314 (adjusted P < 0.001)
#   ZEB2: Moran's I = -0.007, adjusted P = 0.598  (NOT significant)
#
# INPUT : spatial_io_visium_qc.rds from sp01_visium_preprocess.R (kept local)
# OUTPUT: spatially variable genes, domain labels + markers, global and local
#         spatial statistics
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SPARK)        # SPARK-X
  library(spdep)        # Moran's I, Geary's C, LISA, Getis-Ord Gi*
})

SPATIAL_RDS <- Sys.getenv("SPATIAL_RDS",  unset = "spatial_io_visium_qc.rds")
OUT_DIR     <- Sys.getenv("SPATIAL_OUT_DIR", unset = ".")
N_PCS       <- 15
K_NN        <- 6
RESOLUTION  <- 0.5
N_PERM      <- 100
SEED        <- 2024
GENES_OF_INTEREST <- c("ZEB2", "CXCL12", "CXCR4", "COL1A1", "SPARC", "BGLAP")
set.seed(SEED)

spat <- readRDS(SPATIAL_RDS)
message("Spots: ", ncol(spat))

coords <- GetTissueCoordinates(spat)
coords <- as.matrix(coords[, 1:2])
rownames(coords) <- rownames(coords)   # spot barcodes

## ---------------------------------------------------------------------------
## 1. Spatially variable genes — SPARK-X
## ---------------------------------------------------------------------------
cnt <- GetAssayData(spat, assay = "Spatial", layer = "counts")
# SPARK-X expects genes x cells; subset as needed for runtime
spark_res <- sparkx(cnt, coords, numCores = 1, option = "mixture")
svg <- data.frame(gene = rownames(spark_res$res_mtest),
                  spark_res$res_mtest[, c("combinedPvalue", "adjustedPvalue")])
write.csv(svg[order(svg$adjustedPvalue), ],
          file.path(OUT_DIR, "spatial_svg_sparkx.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 2. Spatial domains — Louvain on the first 15 PCs with a 6-NN graph
## ---------------------------------------------------------------------------
spat <- RunPCA(spat, assay = "Spatial", npcs = N_PCS, reduction.name = "pca")

knn  <- knearneigh(coords, k = K_NN)
nb   <- knn2nb(knn)
lw   <- nb2listw(nb, style = "W", zero.policy = TRUE)

# Integrate the neighbour graph with the PCA reduction, then cluster
spat <- FindNeighbors(spat, reduction = "pca", dims = 1:N_PCS,
                      k.param = K_NN, graph.name = "spatial_nn")
spat <- FindClusters(spat, graph.name = "spatial_nn", resolution = RESOLUTION,
                     algorithm = 1)          # algorithm 1 = Louvain
spat$domain <- spat$seurat_clusters
message("Domains: ", length(unique(spat$domain)))

Idents(spat) <- "domain"
markers <- FindAllMarkers(spat, assay = "Spatial", only.pos = TRUE,
                          min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, file.path(OUT_DIR, "spatial_domain_markers.csv"), row.names = FALSE)
write.csv(data.frame(spot = colnames(spat), domain = spat$domain),
          file.path(OUT_DIR, "spatial_domains.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 3. Global autocorrelation — Moran's I and Geary's C
## ---------------------------------------------------------------------------
expr <- as.matrix(GetAssayData(spat, assay = "Spatial", layer = "data"))
genes <- intersect(GENES_OF_INTEREST, rownames(expr))

glob <- do.call(rbind, lapply(genes, function(g) {
  x <- expr[g, ]
  mi <- moran.test(x, lw, zero.policy = TRUE, alternative = "greater")
  gc <- geary.test(x, lw, zero.policy = TRUE, alternative = "less")
  data.frame(gene = g,
             morans_I = mi$estimate[1],
             morans_p = mi$p.value,
             gearys_C = gc$estimate[1],
             gearys_p = gc$p.value)
}))
glob$morans_p_adj  <- p.adjust(glob$morans_p,  "BH")
glob$gearys_p_adj  <- p.adjust(glob$gearys_p,  "BH")
write.csv(glob, file.path(OUT_DIR, "spatial_global_autocorrelation.csv"),
          row.names = FALSE)
print(glob)

## ---------------------------------------------------------------------------
## 4. Local statistics — LISA (100 permutations) and Getis-Ord Gi*
## ---------------------------------------------------------------------------
lisa <- do.call(rbind, lapply(genes, function(g) {
  x <- expr[g, ]
  li <- localmoran_perm(x, lw, nsim = N_PERM, zero.policy = TRUE)
  data.frame(gene = g, spot = names(x),
             Ii = li[, 1], p_perm = li[, 5])
}))
write.csv(lisa, file.path(OUT_DIR, "spatial_lisa.csv"), row.names = FALSE)

gi <- do.call(rbind, lapply(genes, function(g) {
  x <- expr[g, ]
  gstar <- localG(x, lw, zero.policy = TRUE)
  data.frame(gene = g, spot = names(x), Gi = as.numeric(gstar))
}))
write.csv(gi, file.path(OUT_DIR, "spatial_getis_ord_gistar.csv"), row.names = FALSE)

message("Done. Spatial outputs in ", normalizePath(OUT_DIR))
