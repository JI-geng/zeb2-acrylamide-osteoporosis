# Inputs are obtained by running the upstream step.
#
# External validation of the 5-gene model (GSE13850, n = 40).
#
# Packages: rms (lrm), pROC (AUC), rmda (decision curve).
# ---------------------------------------------------------------------------

EXT_EXP   <- "../../02_derived_data/external_validation_GSE13850_expression.csv"
EXT_GRP   <- "../../02_derived_data/external_validation_GSE13850_group_labels.csv"
GENES     <- c("NAP1L2", "ZEB2", "MST1", "PLCB3", "ASNS")
CASE_LV   <- "low BMD"
SAVE_PLOTS <- FALSE

fa0 <- commandArgs(trailingOnly = FALSE)
f0  <- sub("^--file=", "", fa0[grep("^--file=", fa0)])
if (length(f0)) setwd(dirname(normalizePath(f0[1])))

suppressMessages({ library(rms); library(pROC); library(rmda) })

read_cohort <- function(exp_file, grp_file) {
  e <- read.csv(exp_file, row.names = 1, check.names = FALSE)
  g <- read.csv(grp_file, stringsAsFactors = FALSE)
  stopifnot(all(GENES %in% rownames(e)), ncol(e) == nrow(g))
  d <- as.data.frame(t(e[GENES, , drop = FALSE]), check.names = FALSE)
  d$y <- as.integer(g$Group == CASE_LV)
  d
}

ext <- read_cohort(EXT_EXP, EXT_GRP)
cat("external n =", nrow(ext), "(events", sum(ext$y), ")\n")

ddx <- datadist(ext); options(datadist = "ddx")
fit_ext <- lrm(y ~ NAP1L2 + ZEB2 + MST1 + PLCB3 + ASNS, data = ext, x = TRUE, y = TRUE)
p_ext <- as.numeric(predict(fit_ext, type = "fitted"))
auc_ext <- as.numeric(pROC::auc(pROC::roc(ext$y, p_ext, quiet = TRUE)))
cat(sprintf("\nExternal cohort AUC: %.3f\n", auc_ext))
print(fit_ext)

cal <- calibrate(fit_ext, method = "boot", B = 200)
cat("\nExternal-cohort calibration (bootstrap, B = 200):\n"); print(cal)

dca <- decision_curve(y ~ NAP1L2 + ZEB2 + MST1 + PLCB3 + ASNS, data = ext,
                      thresholds = seq(0.01, 0.99, by = 0.01), bootstraps = 0,
                      study.design = "case-control")
cat("Decision curve (external): max net benefit =",
    round(max(as.data.frame(dca$derived.data)$NB, na.rm = TRUE), 3), "\n")

if (SAVE_PLOTS) {
  pdf("diagnostic_model_external.pdf", width = 7, height = 6)
  plot(nomogram(fit_ext, fun = plogis, funlabel = "P(low BMD)"))
  plot(cal, main = "Calibration — external cohort")
  plot_decision_curve(dca)
  dev.off()
}

write.csv(data.frame(estimand = "external_validation",
                     n = nrow(ext),
                     AUC = round(auc_ext, 3)),
          "../../02_derived_data/dm02_external_AUC.csv", row.names = FALSE)
cat("\nWrote ../../02_derived_data/dm02_external_AUC.csv\n")
