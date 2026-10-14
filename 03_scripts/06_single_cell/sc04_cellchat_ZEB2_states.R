#!/usr/bin/env Rscript
# =============================================================================
# Step 6d — CellChat inference of cell-cell communication (Methods 2.6)
#
# As described in the manuscript:
#   * CellChat 2.2.0.9001 with CellChatDB.human
#   * populations containing fewer than 10 cells are excluded
#   * for the ZEB2-state comparison, ZEB2-high and ZEB2-low Tissue stem cells
#     are represented as SEPARATE populations while the remaining cell types
#     keep their original annotations  ->  12 populations
#   * the analysis covers secreted signalling and extracellular-matrix-receptor
#     interactions
#   * random seed = 2024
#   * compared: overall interaction strength, incoming/outgoing signalling,
#     pathway-level and ligand-receptor-level patterns
#
# Reported values to check against (Results 3.5):
#   outgoing strength  ZEB2-high 5.28 vs ZEB2-low 3.94
#   incoming strength  ZEB2-high 2.21 vs ZEB2-low 1.23
#
# INPUT : the annotated object from sc02; the 12-label column must contain
#         Tissue_stem_cells_ZEB2_high / Tissue_stem_cells_ZEB2_low
# OUTPUT: CellChat object, per-population strength table, pathway table
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)      # 2.2.0.9001
})

ANNOTATED_RDS <- Sys.getenv("SC_ANNOTATED_RDS", unset = "05_annotated.rds")
OUT_DIR       <- Sys.getenv("SC_OUT_DIR", unset = ".")
LABEL_COL     <- Sys.getenv("CELLCHAT_LABEL_COL", unset = "cellchat_label")
MIN_CELLS     <- 10
SEED          <- 2024
set.seed(SEED)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

sobj <- readRDS(ANNOTATED_RDS)
if (!LABEL_COL %in% colnames(sobj[[]])) {
  stop("Column '", LABEL_COL, "' not found. Build the 12-label version first ",
       "(Tissue_stem_cells split into _ZEB2_high / _ZEB2_low).")
}

## ---------------------------------------------------------------------------
## 1. Drop populations with fewer than MIN_CELLS cells
## ---------------------------------------------------------------------------
lab  <- sobj[[LABEL_COL]][, 1]
tab  <- table(lab)
drop <- names(tab)[tab < MIN_CELLS]
if (length(drop)) {
  message("Dropping (n < ", MIN_CELLS, "): ", paste(drop, collapse = ", "))
  sobj <- subset(sobj, cells = colnames(sobj)[!(lab %in% drop)])
  lab <- sobj[[LABEL_COL]][, 1]
}

## ---------------------------------------------------------------------------
## 2. Build the CellChat object
## ---------------------------------------------------------------------------
data.input <- GetAssayData(sobj, layer = "data")
meta <- data.frame(labels = lab, row.names = colnames(sobj))

cc <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cc@DB <- subsetDB(CellChatDB.human,
                  search = c("Secreted Signaling", "ECM-Receptor"))
cc <- subsetData(cc)
cc <- identifyOverExpressedGenes(cc)
cc <- identifyOverExpressedInteractions(cc)
cc <- computeCommunProb(cc, type = "triMean", raw.use = FALSE)
cc <- filterCommunication(cc, min.cells = MIN_CELLS)
cc <- computeCommunProbPathway(cc)
cc <- aggregateNet(cc)

saveRDS(cc, file.path(OUT_DIR, "cellchat_object.rds"))

## ---------------------------------------------------------------------------
## 3. Interaction strength per population (Fig. 5D)
## ---------------------------------------------------------------------------
strength <- data.frame(
  population    = levels(factor(lab)),
  outgoing      = cc$net$weight %*% rep(1, ncol(cc$net$weight)) |> as.numeric(),
  incoming      = rep(1, nrow(cc$net$weight)) %*% cc$net$weight |> as.numeric())
write.csv(strength, file.path(OUT_DIR, "cellchat_strength_per_population.csv"),
          row.names = FALSE)
print(strength[order(-strength$outgoing), ])

write.csv(cc$netP$pathways, file.path(OUT_DIR, "cellchat_pathways.csv"),
          row.names = FALSE)
message("Done. CellChat object and tables in ", normalizePath(OUT_DIR))
