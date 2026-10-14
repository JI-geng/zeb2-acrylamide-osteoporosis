#!/usr/bin/env Rscript
# =============================================================================
# Step 6e — virtual knockout of ZEB2 with scTenifoldKnk
#
# Parameters as used in the reported analysis (values taken from the
# module parameter record of the run that produced the reported results;
# see 00_parameters/RUN_PROVENANCE.md):
#   * engine = scTenifoldKnk (bagging_n = 1)
#   * target gene = ZEB2
#   * subset: cells whose `celltype` label is "Tissue_stem_cells"
#   * input gene set: the top 2,000 highly variable genes, after removing
#     genes detected in fewer than 5% of cells (expression_ratio_cutoff = 0.05)
#     -> 965 genes enter the gene-regulatory-network step
#   * number of networks (nc_nNet): 5
#   * cells sampled per network (nc_nCells): 200
#   * QC: minLibSize = 1000, maxMitoRatio = 0.1
#   * significance: p.adjust < 0.05
#   * random seed: 2024; remaining parameters left at their defaults
#
# scTenifoldKnk builds a gene regulatory network from the wild-type data,
# deletes the target gene's edges, and compares the resulting manifold to the
# original; genes with a large displacement are the knock-out responsive genes.
#
# INPUT : raw counts of the single-cell object (kept local)
# OUTPUT: knockout_ZEB2_sigDiff.txt (significant shifted genes), full table
# =============================================================================

suppressPackageStartupMessages({
  library(scTenifoldKnk)
})

COUNTS_RDS   <- Sys.getenv("SC_COUNTS_RDS", unset = "05_annotated.rds")
OUT_DIR      <- Sys.getenv("SC_OUT_DIR", unset = ".")
TARGET       <- "ZEB2"
SUBSET_COL   <- "celltype"          # subset_col
SUBSET_VALUE <- "Tissue_stem_cells" # subset_value
N_HVG        <- 2000   # use_hvg = TRUE, n_hvg = 2000
MIN_DETECT   <- 0.05   # expression_ratio_cutoff
N_NET        <- 5      # nc_nNet
N_BAGGING    <- 1      # bagging_n
N_GENES      <- 965    # genes entering the network step in the reported run
N_CELLS      <- 200    # cells sampled per network (nc_nCells)
MITO_MAX     <- 0.1    # qc_mt_threshold
SEED         <- 2024

obj <- readRDS(COUNTS_RDS)
if (inherits(obj, "Seurat")) {
  counts <- GetAssayData(obj, layer = "counts")
  meta   <- obj[[]]
} else {
  counts <- as.matrix(obj)
  meta   <- NULL
}

## Restrict to the cell population used in the reported run.
if (!is.null(meta) && SUBSET_COL %in% colnames(meta)) {
  keep <- meta[[SUBSET_COL]] == SUBSET_VALUE
  keep[is.na(keep)] <- FALSE
  counts <- counts[, keep, drop = FALSE]
  message("Subset ", SUBSET_COL, " == ", SUBSET_VALUE, ": ", ncol(counts), " cells")
}

counts <- counts[rowSums(counts) > 0, , drop = FALSE]
stopifnot(TARGET %in% rownames(counts))

## QC mirroring the reported run (makeCounts / scQC)
qc <- scQC(counts, minLibSize = 1000, minCells = 3, maxMitoRatio = MITO_MAX)
message("QC: ", ncol(qc$counts), " cells x ", nrow(qc$counts), " genes")

## Highly variable genes, then the 5% detection-rate filter. The reported run
## took the top 2,000 HVG and removed genes detected in < 5% of cells, which
## left 965 genes; N_GENES is applied afterwards only as a guard so that the
## network is built on the same number of genes as in the reported run.
detected <- rowMeans(qc$counts > 0)
qc$counts <- qc$counts[detected >= MIN_DETECT, , drop = FALSE]
message("Genes passing the ", MIN_DETECT, " detection-rate filter: ", nrow(qc$counts))

n_hvg <- min(N_HVG, nrow(qc$counts))
hv <- names(sort(apply(log1p(qc$counts), 1, var), decreasing = TRUE)[seq_len(n_hvg)])
counts_net <- qc$counts[hv, , drop = FALSE]
if (nrow(counts_net) > N_GENES) counts_net <- counts_net[seq_len(N_GENES), , drop = FALSE]
message("Genes entering the network step: ", nrow(counts_net))

set.seed(SEED)
res <- scTenifoldKnk(countMat = counts_net, gKO = TARGET,
                     qc_nGenes = nrow(counts_net), nc_nNet = N_NET,
                     nc_nCells = N_CELLS, nc_bagging = N_BAGGING)

write.csv(res$diffRegulation,
          file.path(OUT_DIR, "knockout_ZEB2_diffRegulation_full.csv"))
sig <- res$diffRegulation[res$diffRegulation$p.adjust < 0.05, ]
write.csv(sig, file.path(OUT_DIR, "knockout_ZEB2_diffRegulation_sig_p050.csv"))
write.table(sig$gene, file.path(OUT_DIR, "knockout_ZEB2_sigDiff.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

message("Target: ", TARGET, " | networks: ", N_NET, " | genes: ", N_GENES)
message("Significant shifted genes (p.adjust < 0.05): ", nrow(sig))
