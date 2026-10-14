#!/usr/bin/env Rscript
# =============================================================================
# Step 6b — cluster annotation and merging into the ten populations of Fig. 4
#
# How the labels were actually assigned in the reported analysis:
#   * the annotation step was run with method = "manual" (marker-based
#     assignment; see the parameter record of the original analysis), and
#   * a reference-based transfer against the Human Primary Cell Atlas ("hpca")
#     was used as supporting evidence (PCA reduction, dims 1:30, at most
#     1,000 cells per cluster for the anchor step).
#   The final labels are the manual ones; this script reproduces the marker
#   table and applies the deposited cluster -> cell type mapping.
#
# INPUT : 04_clustered.rds (output of sc01_qc_cluster_integrate.R; kept local)
#         cluster_to_celltype_10populations.csv  (deposited with this repo)
# OUTPUT: 05_annotated.rds, celltype_annotation.csv
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)      # 5.5.1
})

CLUSTERED_RDS <- Sys.getenv("SC_CLUSTERED_RDS", unset = "04_clustered.rds")
MAPPING_CSV   <- Sys.getenv(
  "SC_MAPPING_CSV",
  unset = file.path(dirname(sub("^--file=", "",
              grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))),
              "cluster_to_celltype_10populations.csv"))
OUT_DIR       <- Sys.getenv("SC_OUT_DIR", unset = ".")
CLUSTER_COL   <- "seurat_clusters"

sobj <- readRDS(CLUSTERED_RDS)
mapping <- read.csv(MAPPING_CSV, stringsAsFactors = FALSE)

## ---------------------------------------------------------------------------
## 1. Cluster-level markers (Wilcoxon) — the evidence used for manual labelling
## ---------------------------------------------------------------------------
Idents(sobj) <- CLUSTER_COL
markers <- FindAllMarkers(sobj, only.pos = TRUE,
                          min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, file.path(OUT_DIR, "all_cluster_markers.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 2. Apply the deposited mapping
##    NOTE ON THE CLUSTER COUNT: the annotation table deposited here covers
##    clusters 0-15 (16 clusters, the object that was actually annotated).
##    FindClusters(resolution = 0.5) on a freshly processed object returns 17
##    clusters; if your run does the same, re-derive the mapping from
##    all_cluster_markers.csv before continuing — do not assume the deposited
##    mapping applies unchanged.
## ---------------------------------------------------------------------------
obs <- as.character(unique(sobj[[CLUSTER_COL]][, 1]))
if (!all(obs %in% mapping$cluster)) {
  stop("Clusters ", paste(setdiff(obs, mapping$cluster), collapse = ", "),
       " have no entry in ", MAPPING_CSV,
       ". Re-derive the mapping from all_cluster_markers.csv.")
}

sobj$celltype <- mapping$celltype[match(as.character(sobj[[CLUSTER_COL]][, 1]),
                                        mapping$cluster)]
message("Populations: ", paste(sort(unique(sobj$celltype)), collapse = ", "))

## ---------------------------------------------------------------------------
## 3. ZEB2-high / ZEB2-low states
##    As described in the Methods: cells are ranked by normalised ZEB2
##    expression and the UPPER and LOWER QUARTILES are designated ZEB2-high and
##    ZEB2-low; the two intermediate quartiles are EXCLUDED from the comparison.
##    Applied within the Tissue stem cell population.
##    This is also what produces the 12-label version used by the communication
##    step (Tissue_stem_cells -> _ZEB2_high / _ZEB2_low); see sc04.
## ---------------------------------------------------------------------------
zeb2 <- as.numeric(GetAssayData(sobj, layer = "data")["ZEB2", ])
q25 <- quantile(zeb2, 0.25, na.rm = TRUE)
q75 <- quantile(zeb2, 0.75, na.rm = TRUE)
sobj$ZEB2_group <- ifelse(zeb2 >= q75, "ZEB2_high",
                   ifelse(zeb2 <= q25, "ZEB2_low", NA_character_))
message("ZEB2 states (NA = intermediate quartiles, excluded):")
print(table(sobj$ZEB2_group, useNA = "ifany"))

saveRDS(sobj, file.path(OUT_DIR, "05_annotated.rds"))
write.csv(sobj[[]], file.path(OUT_DIR, "celltype_annotation.csv"))
message("Done: ", length(unique(sobj$celltype)), " populations")
