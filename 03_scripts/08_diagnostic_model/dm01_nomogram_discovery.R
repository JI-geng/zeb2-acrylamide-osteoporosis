# ---------------------------------------------------------------------------
# Diagnostic model — discovery cohort (n = 106, 5 genes)
#     rms   (lrm / nomogram / calibrate)   -> logistic model, nomogram, calibration
#     pROC  (roc / auc)                    -> AUC
#     rmda  (decision_curve)               -> decision-curve analysis (net benefit)
#
# ---------------------------------------------------------------------------

EXP_FILE <- "../../02_derived_data/expression_matrix_all_12487genes_106samples.csv"
GRP_FILE <- "../../02_derived_data/group_labels_106samples.csv"
GENES    <- c("NAP1L2", "ZEB2", "MST1", "PLCB3", "ASNS")  # order recorded in genes.csv of the original run
CASE_LV  <- "low BMD"          # event coded as 1 (the group label used by the ML scripts)
B_BOOT   <- 200                # bootstrap replicates (use 1000 for the final run)
SAVE_PLOTS <- FALSE            # figure rendering is not part of this deposit (README §1)

# locate this script and run from its own directory, so the paths above work
fa0 <- commandArgs(trailingOnly = FALSE)
f0  <- sub("^--file=", "", fa0[grep("^--file=", fa0)])
if (length(f0)) setwd(dirname(normalizePath(f0[1])))

suppressMessages({
  library(rms); library(pROC); library(rmda)
})

expr <- read.csv(EXP_FILE, row.names = 1, check.names = FALSE)
grp  <- read.csv(GRP_FILE, stringsAsFactors = FALSE)

stopifnot(all(GENES %in% rownames(expr)))
stopifnot(ncol(expr) == nrow(grp))

df <- as.data.frame(t(expr[GENES, , drop = FALSE]), check.names = FALSE)
df$y <- as.integer(grp$Group == CASE_LV)
cat("samples:", nrow(df), " events (low BMD):", sum(df$y), " controls:", sum(df$y == 0), "\n")

dd <- datadist(df); options(datadist = "dd")
fit <- lrm(y ~ NAP1L2 + ZEB2 + MST1 + PLCB3 + ASNS, data = df, x = TRUE, y = TRUE)
print(fit)

## --- apparent AUC (the reported 0.856) -----------------------------------------
p  <- as.numeric(predict(fit, type = "fitted"))
r  <- pROC::roc(df$y, p, quiet = TRUE)
auc_apparent <- as.numeric(pROC::auc(r))
cat(sprintf("\nApparent AUC (in-sample, manuscript value 0.856): %.3f\n", auc_apparent))

## --- per-gene AUC (manuscript: NAP1L2 .748 / ZEB2 .760 / MST1 .645 / PLCB3 .710 / ASNS .707)
per_gene <- vapply(GENES, function(g) as.numeric(pROC::auc(pROC::roc(df$y, df[[g]], quiet = TRUE))), numeric(1))
per_gene <- data.frame(gene = GENES, AUC = round(per_gene, 3))
cat("\nSingle-gene AUC:\n"); print(per_gene, row.names = FALSE)

## --- bootstrap optimism-corrected AUC (transparency only) ----------------------
set.seed(411)
opt <- numeric(B_BOOT)
for (b in seq_len(B_BOOT)) {
  idx <- sample(seq_len(nrow(df)), replace = TRUE)
  db  <- df[idx, ]
  fb  <- try(lrm(y ~ NAP1L2 + ZEB2 + MST1 + PLCB3 + ASNS, data = db, x = TRUE, y = TRUE), silent = TRUE)
  if (inherits(fb, "try-error")) next
  pb  <- as.numeric(predict(fb, type = "fitted"))
  if (length(unique(db$y)) < 2) next
  a_boot <- as.numeric(pROC::auc(pROC::roc(db$y, pb, quiet = TRUE)))
  p_orig <- as.numeric(predict(fb, newdata = df, type = "fitted"))
  a_orig <- as.numeric(pROC::auc(pROC::roc(df$y, p_orig, quiet = TRUE)))
  opt[b] <- a_boot - a_orig
}
opt <- opt[is.finite(opt)]
cat(sprintf("Bootstrap optimism (%d replicates): %.3f  ->  optimism-corrected AUC: %.3f\n",
            length(opt), mean(opt), auc_apparent - mean(opt)))

## --- calibration (rms, bootstrap) ---------------------------------------------
cal <- calibrate(fit, method = "boot", B = B_BOOT)
cat("\nCalibration (bootstrap, B =", B_BOOT, "):\n"); print(cal)

## --- decision curve (rmda) -----------------------------------------------------
dca <- decision_curve(y ~ NAP1L2 + ZEB2 + MST1 + PLCB3 + ASNS,
                      data = df, thresholds = seq(0.01, 0.99, by = 0.01),
                      bootstraps = 0, study.design = "case-control")
nb <- as.data.frame(dca$derived.data)
cat("\nDecision curve: max net benefit =", round(max(nb$NB, na.rm = TRUE), 3), "\n")

if (SAVE_PLOTS) {
  pdf("diagnostic_model_discovery.pdf", width = 7, height = 6)
  plot(nomogram(fit, fun = plogis, funlabel = "P(low BMD)"))
  plot(cal, main = "Calibration (bootstrap)")
  plot_decision_curve(dca)
  dev.off()
}

write.csv(per_gene, "../../02_derived_data/dm01_single_gene_AUC.csv", row.names = FALSE)
cat("\nWrote ../../02_derived_data/dm01_single_gene_AUC.csv\n")
