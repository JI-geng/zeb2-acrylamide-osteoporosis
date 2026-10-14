#!/usr/bin/env Rscript
# =============================================================================
# Step 6c — CytoTRACE 2 differentiation-potential scoring
#
# CytoTRACE 2 (v1.0.0) ranks cells on a differentiation-potential score
# (0 = differentiated, 1 = stem-like) and was used to show that the
# ZEB2-high population sits at the less-differentiated end of the trajectory.
#
# PARAMETER NOTE: the run that produced the manuscript figure did not record
# the CytoTRACE 2 options in its log. This script therefore uses the package
# defaults on the same input (annotated object, cell type labels). If you need
# bit-identical numbers, set smooth / ncores explicitly and record them.
#
# INPUT : 05_annotated.rds (kept local)
# OUTPUT: cytotrace2_scores.csv (one row per cell: score + predicted class)
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(CytoTRACE2)
})

ANNOTATED_RDS <- Sys.getenv("SC_ANNOTATED_RDS", unset = "05_annotated.rds")
OUT_DIR       <- Sys.getenv("SC_OUT_DIR", unset = ".")
CELLTYPE_COL  <- "celltype"

sobj <- readRDS(ANNOTATED_RDS)

expr <- as.matrix(GetAssayData(sobj, layer = "counts"))
pheno <- data.frame(cell_name = colnames(sobj),
                    cell_type = sobj[[CELLTYPE_COL]][, 1],
                    stringsAsFactors = FALSE)

## CytoTRACE 2: defaults, on raw counts
res <- cytotrace2(expr, annotation = pheno, species = "human")

## res$CytoTRACE2_Score   : continuous potency score
## res$CytoTRACE2_Potency : discretised class (Differentiated / Unipotent /
##                          Oligopotent / Multipotent / Pluripotent / Totipotent)
## res$CytoTRACE2_Relative: score rescaled 0-1 within the sample
out <- data.frame(cell = colnames(sobj),
                  celltype = sobj[[CELLTYPE_COL]][, 1],
                  ZEB2_group = sobj$ZEB2_group,
                  score = res$CytoTRACE2_Score,
                  potency = res$CytoTRACE2_Potency,
                  relative = res$CytoTRACE2_Relative)

write.csv(out, file.path(OUT_DIR, "cytotrace2_scores.csv"), row.names = FALSE)

## Summary used for the manuscript comparison (ZEB2-high vs ZEB2-low)
agg <- aggregate(score ~ ZEB2_group, data = out, FUN = function(x)
                 c(n = length(x), mean = mean(x), sd = sd(x)))
print(do.call(rbind, lapply(seq_len(nrow(agg)), function(i)
        data.frame(group = agg$ZEB2_group[i], as.list(agg$score[, i])))))
write.csv(do.call(rbind, lapply(seq_len(nrow(agg)), function(i)
        data.frame(group = agg$ZEB2_group[i], as.list(agg$score[, i])))),
        file.path(OUT_DIR, "cytotrace2_ZEB2_group_summary.csv"), row.names = FALSE)
message("Done. Scores written to ", normalizePath(OUT_DIR))
