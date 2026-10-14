# ============================================================================
# DEG analysis: GSE7158 + GSE56815 -> ComBat -> limma -> DEGs -> CTD intersection
#
# Reproduction script for the transcriptomic layer of the manuscript.
#
# ----------------------------------------------------------------------------
# INPUT — THE PUBLIC FILES MUST BE DOWNLOADED BY HAND FIRST.
# This script downloads nothing. Put three files in ../../02_derived_data/
# (or point GEO_DATA_DIR at another directory):
#
#   1. GSE7158_series_matrix.txt.gz
#      https://ftp.ncbi.nlm.nih.gov/geo/series/GSE7nnn/GSE7158/matrix/GSE7158_series_matrix.txt.gz
#   2. GSE56815_series_matrix.txt.gz
#      https://ftp.ncbi.nlm.nih.gov/geo/series/GSE56nnn/GSE56815/matrix/GSE56815_series_matrix.txt.gz
#      (or use the "Series Matrix File(s)" link on the GEO record page)
#   3. ctd_acrylamide_genes.txt
#      CTD -> Downloads -> Chemical–gene interactions
#      (https://ctdbase.org/downloads/); filter the chemical to acrylamide
#      (CAS 79-06-1) and save one gene symbol per line.
#
#   Optional: sample_group_labels.csv with columns `sample` and `group`
#   (`high` / `low` BMD). Supply it if the GEO sample titles do not say which
#   group a sample belongs to — the script stops and asks for it in that case.
#
# ----------------------------------------------------------------------------
# Datasets
#   GSE7158  (GPL570, hgu133plus2, n = 26)  peripheral blood monocytes
#   GSE56815 (GPL96,  hgu133a,     n = 80)  peripheral blood monocytes
#   total after merge: 106 (54 high BMD / 52 low BMD)
#
# Parameters (00_parameters/WORKFLOW_PARAMETERS.md, Step 2)
#   batch correction   ComBat (sva)
#   DEG method         limma (lmFit + eBayes)
#   threshold          nominal P < 0.05
#
# Expected (manuscript): 446 up / 695 down = 1,141 DEGs; CTD ∩ DEG = 169 genes.
# ============================================================================

library(limma)
library(sva)

## --- configuration ----------------------------------------------------------
DERIVED    <- Sys.getenv("GEO_DATA_DIR", unset = "../../02_derived_data")
F_7158     <- file.path(DERIVED, "GSE7158_series_matrix.txt.gz")
F_56815    <- file.path(DERIVED, "GSE56815_series_matrix.txt.gz")
F_CTD      <- file.path(DERIVED, "ctd_acrylamide_genes.txt")
F_GROUPS   <- file.path(DERIVED, "sample_group_labels.csv")
PLATFORM   <- c(GSE7158 = "GPL570", GSE56815 = "GPL96")
DEG_P      <- 0.05
LOG2FC     <- 0

need <- c(F_7158, F_56815, F_CTD)
miss <- need[!file.exists(need)]
if (length(miss))
  stop("\nMissing input file(s):\n  ", paste(miss, collapse = "\n  "),
       "\nDownload them by hand — see the INPUT block at the top of this script.\n")

## --- helpers -----------------------------------------------------------------
read_series_matrix <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  lines <- readLines(con, warn = FALSE); close(con)
  b <- grep("^!series_matrix_table_begin", lines)
  e <- grep("^!series_matrix_table_end", lines)
  if (!length(b) || !length(e)) stop("not a GEO series matrix file: ", path)
  tbl <- read.delim(text = paste(lines[(b + 1):(e - 1)], collapse = "\n"),
                    header = TRUE, row.names = 1, check.names = FALSE,
                    quote = "", comment.char = "")
  meta <- function(key) {
    ln <- grep(paste0("^", key, "\t"), lines, value = TRUE)
    if (!length(ln)) return(NULL)
    strsplit(sub(paste0("^", key, "\t"), "", ln[1]), "\t")[[1]]
  }
  list(expr = as.matrix(tbl),
       titles = meta("!Sample_title"),
       chars  = meta("!Sample_characteristics_ch1"),
       geo    = meta("!Sample_geo_accession"))
}

map_probes_to_genes <- function(expr, platform) {
  # prefer the Bioconductor annotation package for the platform
  db <- switch(platform,
               GPL570 = "hgu133plus2.db",
               GPL96  = "hgu133a.db",
               stop("unknown platform: ", platform))
  if (!requireNamespace(db, quietly = TRUE))
    stop("annotation package '", db, "' is not installed.\n",
         "Install it with: BiocManager::install(\"", db, "\")\n",
         "or supply a local probe->gene table (see the README of this directory).")
  sym <- AnnotationDbi::mapIds(get(db), keys = rownames(expr), column = "SYMBOL",
                               keytype = "PROBEID", multiVals = "first")
  keep <- !is.na(sym) & sym != ""
  expr <- expr[keep, , drop = FALSE]; sym <- sym[keep]
  # collapse probes for the same gene by mean
  ex <- t(apply(expr, 1, function(x) x))   # probes x samples
  sp <- split(seq_len(nrow(ex)), sym)
  m  <- vapply(sp, function(i) if (length(i) == 1) ex[i, ] else colMeans(ex[i, , drop = FALSE]),
               numeric(ncol(ex)))
  t(m)   # genes x samples
}

## --- 1. read the two locally downloaded series matrices ----------------------
cat("=== 1. Reading local GEO series matrices ===\n")
s7158  <- read_series_matrix(F_7158)
s56815 <- read_series_matrix(F_56815)
cat("GSE7158 :", ncol(s7158$expr), "samples,", nrow(s7158$expr), "probes\n")
cat("GSE56815:", ncol(s56815$expr), "samples,", nrow(s56815$expr), "probes\n")

## --- 2. group labels --------------------------------------------------------
cat("\n=== 2. Group labels (high vs low BMD) ===\n")
guess_group <- function(s) {
  txt <- paste(s$titles, s$chars, sep = " ")
  txt <- tolower(paste(txt, collapse = " | "))
  txt <- strsplit(txt, "\\|")[[1]]
  vapply(strsplit(txt, " ")[seq_along(txt)], function(w) {
    if (any(grepl("low",  w))) "low"  else
    if (any(grepl("high", w))) "high" else NA_character_
  }, character(1))
}
g1 <- guess_group(s7158)
g2 <- guess_group(s56815)
if (any(is.na(g1)) || any(is.na(g2))) {
  if (!file.exists(F_GROUPS))
    stop("cannot read high/low BMD labels from the GEO sample titles.\n",
         "Write ", F_GROUPS, " with columns `sample,group` (high/low) and re-run.")
  lab <- read.csv(F_GROUPS, stringsAsFactors = FALSE)
  id  <- c(s7158$geo, s56815$geo)
  m   <- setNames(tolower(lab$group), lab$sample)
  g1  <- m[s7158$geo]; g2 <- m[s56815$geo]
  if (any(is.na(c(g1, g2)))) stop("sample_group_labels.csv does not cover every sample")
}
group <- factor(c(g1, g2), levels = c("high", "low"))
cat("high:", sum(group == "high"), " low:", sum(group == "low"), "\n")
if (sum(group == "high") != 54 || sum(group == "low") != 52)
  cat("⚠️ expected 54 high / 52 low — check the labels before trusting the DEG counts\n")

## --- 3. probe -> gene -------------------------------------------------------
cat("\n=== 3. Probe -> gene mapping ===\n")
e7158  <- map_probes_to_genes(s7158$expr,  PLATFORM[["GSE7158"]])
e56815 <- map_probes_to_genes(s56815$expr, PLATFORM[["GSE56815"]])
cat("genes: GSE7158", nrow(e7158), " GSE56815", nrow(e56815), "\n")

common <- intersect(rownames(e7158), rownames(e56815))
cat("common genes:", length(common), "\n")
expr_merged <- cbind(e7158[common, ], e56815[common, ])
stopifnot(ncol(expr_merged) == length(group))

## --- 4. ComBat --------------------------------------------------------------
cat("\n=== 4. ComBat batch correction ===\n")
batch <- c(rep("GSE7158", ncol(e7158)), rep("GSE56815", ncol(e56815)))
mod   <- model.matrix(~group)
expr_cb <- ComBat(dat = expr_merged, batch = batch, mod = mod, par.prior = TRUE)
cat("corrected matrix:", nrow(expr_cb), "genes x", ncol(expr_cb), "samples\n")

## --- 5. limma ---------------------------------------------------------------
cat("\n=== 5. limma differential expression ===\n")
design <- model.matrix(~group)
fit <- lmFit(expr_cb, design)
fit <- eBayes(fit)
deg <- topTable(fit, coef = "grouplow", number = Inf, adjust = "none")
deg <- deg[order(deg$P.Value), ]
deg_sig <- deg[deg$P.Value < DEG_P & abs(deg$logFC) > LOG2FC, ]

n_up   <- sum(deg_sig$logFC > 0)
n_down <- sum(deg_sig$logFC < 0)
cat("DEGs (nominal P <", DEG_P, "):", nrow(deg_sig), "\n")
cat("  up (high vs low):", n_up,   "  [manuscript: 446]\n")
cat("  down            :", n_down, "  [manuscript: 695]\n")

## --- 6. CTD intersection ----------------------------------------------------
cat("\n=== 6. CTD intersection ===\n")
ctd <- readLines(F_CTD, warn = FALSE)
ctd <- trimws(ctd); ctd <- ctd[nzchar(ctd)]
inter <- intersect(ctd, rownames(deg_sig))
cat("CTD acrylamide genes:", length(ctd), "\n")
cat("DEGs:", nrow(deg_sig), "\n")
cat("intersection:", length(inter), " [manuscript: 169]\n")

## --- 7. outputs -------------------------------------------------------------
cat("\n=== 7. Writing outputs ===\n")
write.csv(deg,     file.path(DERIVED, "deg_table.csv"),       row.names = TRUE)
write.csv(deg_sig, file.path(DERIVED, "deg_significant.csv"), row.names = TRUE)
write.csv(data.frame(gene = sort(inter)), file.path(DERIVED, "ctd_deg_intersection.csv"),
          row.names = FALSE)
# one gene per line, the format 03_scripts/01_enrichment/_enrich169.py reads
writeLines(sort(inter), file.path(DERIVED, "169_genes.csv"))
saveRDS(expr_cb, file.path(DERIVED, "merged_expr_combat.rds"))
cat("wrote deg_table.csv, deg_significant.csv, ctd_deg_intersection.csv,",
    "169_genes.csv, merged_expr_combat.rds to", DERIVED, "\n")
cat("\nDone. Compare the counts above with 446 / 695 / 1,141 / 169 before using the output.\n")
