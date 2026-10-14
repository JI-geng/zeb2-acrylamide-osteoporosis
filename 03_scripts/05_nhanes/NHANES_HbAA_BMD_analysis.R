# ============================================================================
# NHANES 2005-2006 / 2013-2014
# Haemoglobin adduct of acrylamide (HbAA) vs proximal femur bone mineral density
#
# Reproduction script for the population (exposure-outcome) layer of the
# manuscript. It reads the public CDC files **from a directory you create** —
# the script downloads nothing, see the INPUT block below.
#
# ----------------------------------------------------------------------------
# Summary of the analysis
#   Cycles      2005-2006 and 2013-2014 (two-year MEC-examined samples)
#   Exposure    LBXACR  acrylamide adduct, pmol/g haemoglobin  (AMDGYD_D / AMDGYD_H)
#   Outcome     DXXOFBMD  total femur BMD, g/cm2              (DXXFEM_D / DXXFEM_H)
#   Population  age >= 50 years with a non-missing exposure, outcome, BMI,
#               poverty-income ratio and sample weight
#   Subgroup    pre-specified: never smokers aged >= 65 years
#   Models      survey-weighted linear regression (svyglm), three adjustment
#               levels; exposure entered as a z-score (per 1 SD)
#   Weights     2005-2006 -> WTMEC2YR (examined sample)
#               2013-2014 -> WTSA2YR  (Subsample A; acrylamide was measured in a
#                                      one-third subsample in this cycle)
#               each divided by 2 for the two-cycle combination
#
# ----------------------------------------------------------------------------
# Corrections relative to the earlier draft of this script (all verified against
# the live CDC files on 2026-10-12):
#   1. Download URL   /Nchs/Nhanes/<cycle>/<file> is dead (returns a "Page Not
#                     Found" HTML page).  Current path is
#                     /Nchs/Data/Nhanes/Public/<release-year>/DataFiles/<file>.XPT
#   2. Exposure file  ACD_APP_D / ACD_HAA_H do not exist.
#                     Real names: AMDGYD_D (2005-2006), AMDGYD_H (2013-2014).
#   3. Exposure name  LBDHBAA does not exist. Real variable: LBXACR.
#                     (LBDACRLC flags results below the LOD; NHANES already
#                      substituted LLOD/sqrt(2) for those, 0.6% of this sample.)
#   4. Outcome file   BMX_FEM_D does not exist (2005-2006) and BMX_H is body
#                     measures, not bone densitometry (2013-2014).
#                     Real names: DXXFEM_D, DXXFEM_H.
#   5. Outcome name   BMXFMD does not exist in these files.
#                     Real variable: DXXOFBMD ("Total femur BMD", g/cm2).
#   6. BMI            was never merged (models using it could not run); taken
#                     from BMX_D / BMX_H.
#   7. PIR            formula used `PIR`, the variable is named `pir`.
#   8. Weight         2013-2014 must use the subsample weight WTSA2YR, not
#                     WTMEC2YR. Using WTMEC2YR for that cycle changes the
#                     never-smoker estimate from beta = -0.031 (P = 0.015) to
#                     beta = -0.021 (P = 0.17) - i.e. it removes the finding.
#   9. Smoking        SMQ020 alone mislabels former smokers as "current";
#                     SMQ020 (lifetime >= 100 cigarettes) is combined with
#                     SMQ040 (current smoking frequency).
#  10. Confidence     intervals are taken from confint() (t, survey df) rather
#                     than a normal 1.96 x SE approximation.
#                     intervals
# ============================================================================

suppressMessages({
  library(survey)   # svydesign / svyglm / confint
  library(dplyr)
  library(haven)    # read_xpt
})

# ----------------------------------------------------------------------------
# INPUT — THE CDC FILES MUST BE DOWNLOADED BY HAND BEFORE THIS SCRIPT IS RUN.
# Nothing in this repository downloads anything.
#
#   Download pages (select the cycle, then the "Data Files" table):
#     2005-2006   https://wwwn.cdc.gov/Nchs/Nhanes/Search/DataPage.aspx?Component=Demographics&Cycle=2005-2006
#     2013-2014   https://wwwn.cdc.gov/Nchs/Nhanes/Search/DataPage.aspx?Component=Demographics&Cycle=2013-2014
#   Direct URL pattern for an individual file (verified 2026-10-12):
#     https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/<release-year>/DataFiles/<file>.XPT
#     e.g. https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/DEMO_D.XPT
#
#   Put all 16 unrenamed files in one flat directory and point NHANES_DATA_DIR at it:
#     export NHANES_DATA_DIR=/path/to/nhanes_xpt        # Linux / macOS / Git Bash
#     set    NHANES_DATA_DIR=C:\path\to\nhanes_xpt      # Windows cmd.exe
#     Sys.setenv(NHANES_DATA_DIR = "/path/to/nhanes_xpt")   # inside R
#
#   2005-2006: DEMO_D  AMDGYD_D  DXXFEM_D  BMX_D  SMQ_D  ALQ_D  BPX_D  DIQ_D
#   2013-2014: DEMO_H  AMDGYD_H  DXXFEM_H  BMX_H  SMQ_H  ALQ_H  BPX_H  DIQ_H
# ----------------------------------------------------------------------------
DATA_DIR <- Sys.getenv("NHANES_DATA_DIR", unset = "")

REQUIRED <- c("DEMO_D","AMDGYD_D","DXXFEM_D","BMX_D","SMQ_D","ALQ_D","BPX_D","DIQ_D",
              "DEMO_H","AMDGYD_H","DXXFEM_H","BMX_H","SMQ_H","ALQ_H","BPX_H","DIQ_H")

if (!nzchar(DATA_DIR) || !dir.exists(DATA_DIR)) {
  stop("\nNHANES_DATA_DIR is not set or does not exist.\n",
       "Download the 16 CDC .XPT files listed in the header of this script and set\n",
       "NHANES_DATA_DIR to the directory that contains them. Nothing is downloaded\n",
       "automatically.\n")
}

find_xpt <- function(file) {
  for (f in file.path(DATA_DIR, c(paste0(file, ".XPT"), paste0(file, ".xpt")))) {
    if (file.exists(f)) return(f)
  }
  NA_character_
}
missing_files <- REQUIRED[vapply(REQUIRED, function(f) is.na(find_xpt(f)), logical(1))]
if (length(missing_files))
  stop("missing .XPT file(s) in ", DATA_DIR, ": ", paste(missing_files, collapse = ", "), "\n")

fetch <- function(release_year, file) {
  # reads a file the user downloaded; never downloads
  dest <- find_xpt(file)
  if (is.na(dest)) stop("file not found in NHANES_DATA_DIR: ", file)
  if (file.size(dest) < 1000) stop("file too small to be a SAS XPORT file: ", dest)
  if (!identical(rawToChar(readBin(dest, "raw", n = 13)), "HEADER RECORD"))
    stop("not a SAS XPORT file (the CDC download pages sometimes return an HTML page): ", dest)
  read_xpt(dest)
}

# ----------------------------------------------------------------------------
# 1. Assemble one cycle
# ----------------------------------------------------------------------------
build_cycle <- function(release_year, demo_f, hb_f, dxx_f, bmx_f,
                        smq_f, alq_f, bpx_f, diq_f, cycle) {
  demo <- fetch(release_year, demo_f)
  hb   <- fetch(release_year, hb_f)
  dxx  <- fetch(release_year, dxx_f)
  bmx  <- fetch(release_year, bmx_f)
  smq  <- fetch(release_year, smq_f)
  alq  <- fetch(release_year, alq_f)
  bpx  <- fetch(release_year, bpx_f)
  diq  <- fetch(release_year, diq_f)

  # Keep both candidate weights. The laboratory file ships a subsample weight
  # when the analyte was measured on a subsample; otherwise the examined-sample
  # weight applies.
  if ("WTSA2YR" %in% names(hb)) {
    weight_source <- "WTSA2YR"
    wts <- hb %>% transmute(SEQN, wt_sub = WTSA2YR)
  } else {
    weight_source <- "WTMEC2YR"
    wts <- hb %>% transmute(SEQN, wt_sub = NA_real_)
  }
  mec <- demo %>% transmute(SEQN, wt_mec = WTMEC2YR)

  demo %>%
    transmute(SEQN, age = RIDAGEYR, sex = RIAGENDR, race = RIDRETH1,
              pir = INDFMPIR, sdmvpsu = SDMVPSU, sdmvstra = SDMVSTRA) %>%
    left_join(mec, by = "SEQN") %>%
    left_join(wts, by = "SEQN") %>%
    left_join(hb  %>% transmute(SEQN, HbAA = LBXACR, below_lod = LBDACRLC), by = "SEQN") %>%
    left_join(dxx %>% transmute(SEQN, bmd = DXXOFBMD, femur_status = DXAFMRST,
                                femur_invalidity = DXXFMBCC), by = "SEQN") %>%
    left_join(bmx %>% transmute(SEQN, bmi = BMXBMI), by = "SEQN") %>%
    left_join(smq %>% transmute(SEQN, SMQ020, SMQ040), by = "SEQN") %>%
    left_join(alq %>% transmute(SEQN, ALQ101), by = "SEQN") %>%
    left_join(bpx %>% transmute(SEQN, BPXSY1, BPXSY2, BPXSY3), by = "SEQN") %>%
    left_join(diq %>% transmute(SEQN, DIQ010), by = "SEQN") %>%
    mutate(cycle = cycle, weight_source = weight_source)
}

d05 <- build_cycle("2005", "DEMO_D", "AMDGYD_D", "DXXFEM_D", "BMX_D",
                   "SMQ_D", "ALQ_D", "BPX_D", "DIQ_D", "2005-2006")
d13 <- build_cycle("2013", "DEMO_H", "AMDGYD_H", "DXXFEM_H", "BMX_H",
                   "SMQ_H", "ALQ_H", "BPX_H", "DIQ_H", "2013-2014")

cat("Sample weight used per cycle:\n")
print(distinct(bind_rows(d05, d13), cycle, weight_source), n = 2)

# ----------------------------------------------------------------------------
# 2. Eligibility, recoding, pooling
# ----------------------------------------------------------------------------
df <- bind_rows(d05, d13) %>%
  # the weight actually used = subsample weight where one exists, else MEC weight
  mutate(wt_raw = ifelse(!is.na(wt_sub), wt_sub, wt_mec)) %>%
  # NHANES sets total femur BMD to missing whenever any femur sub-region is
  # invalid, so non-missing DXXOFBMD already implies a valid scan.
  filter(age >= 50,
         !is.na(HbAA), !is.na(bmd), !is.na(bmi), !is.na(pir),
         !is.na(wt_raw), wt_raw > 0,
         !is.na(SMQ020),            # smoking is a model covariate and defines
         !is.na(ALQ101)) %>%        # the pre-specified subgroup
  mutate(
    sex = factor(sex, 1:2, c("Male", "Female")),
    race = factor(race, 1:5, c("Mexican American", "Other Hispanic",
                               "Non-Hispanic White", "Non-Hispanic Black",
                               "Other Race")),
    # SMQ020 = ever smoked >= 100 cigarettes (1 yes / 2 no)
    # SMQ040 = current frequency (1 every day / 2 some days / 3 not at all)
    smoke = factor(
      case_when(SMQ020 == 2                   ~ "Never",
                SMQ020 == 1 & SMQ040 %in% 1:2 ~ "Current",
                SMQ020 == 1 & SMQ040 == 3     ~ "Former"),
      levels = c("Never", "Former", "Current")),
    alcohol = factor(ifelse(ALQ101 == 1, "Drinker",
                            ifelse(ALQ101 == 2, "Non-drinker", NA)),
                     levels = c("Non-drinker", "Drinker")),
    diabetes = factor(ifelse(DIQ010 == 1, "Diabetes",
                             ifelse(DIQ010 %in% 2:3, "No diabetes", NA)),
                      levels = c("No diabetes", "Diabetes")),
    sbp = rowMeans(cbind(BPXSY1, BPXSY2, BPXSY3), na.rm = TRUE),
    hypertension = factor(ifelse(!is.na(sbp) & sbp >= 140, "Yes",
                                 ifelse(!is.na(sbp), "No", NA)),
                          levels = c("No", "Yes")),
    # two cycles -> halve each cycle's own two-year weight
    wt = wt_raw / 2,
    # PSU and stratum identifiers are cycle-specific; make them unique so that
    # pooling the two cycles cannot alias two different PSUs into one.
    psu = as.numeric(interaction(cycle, sdmvpsu,  drop = TRUE)),
    stratum = as.numeric(interaction(cycle, sdmvstra, drop = TRUE))
  )

hb_sd <- sd(df$HbAA, na.rm = TRUE)
df <- df %>% mutate(HbAA_sd = (HbAA - mean(HbAA, na.rm = TRUE)) / hb_sd)

cat(sprintf("\nN pooled = %d  (2005-2006 = %d, 2013-2014 = %d)\n",
            nrow(df), sum(df$cycle == "2005-2006"), sum(df$cycle == "2013-2014")))
cat(sprintf("HbAA  mean = %.2f pmol/g Hb, SD = %.2f, range %.1f-%.1f; below LOD = %.1f%%\n",
            mean(df$HbAA, na.rm = TRUE), hb_sd,
            min(df$HbAA, na.rm = TRUE), max(df$HbAA, na.rm = TRUE),
            100 * mean(df$below_lod == 1, na.rm = TRUE)))
cat(sprintf("BMD   mean = %.4f g/cm2 (SD %.4f)\n",
            mean(df$bmd, na.rm = TRUE), sd(df$bmd, na.rm = TRUE)))

design <- svydesign(ids = ~ psu, strata = ~ stratum, weights = ~ wt,
                    data = df, nest = TRUE)

results <- list()
cur_sample <- "overall_age_ge50"

report <- function(fit, label, n_unweighted) {
  ci <- confint(fit, parm = "HbAA_sd", method = "Wald", df = degf(fit))
  b <- unname(coef(fit)["HbAA_sd"])
  p <- summary(fit)$coefficients["HbAA_sd", "Pr(>|t|)"]
  cat(sprintf("  %-34s beta = %+.4f g/cm2 per 1 SD  (95%% CI %+.4f to %+.4f)  P = %.3g\n",
              label, b, unname(ci[1]), unname(ci[2]), p))
  results[[length(results) + 1]] <<- data.frame(
    sample             = cur_sample,
    model              = label,
    n_unweighted       = n_unweighted,
    beta_per_1SD_g_cm2 = b,
    ci_low             = unname(ci[1]),
    ci_high            = unname(ci[2]),
    p_value            = p,
    stringsAsFactors   = FALSE
  )
  invisible(NULL)
}

# ----------------------------------------------------------------------------
# 3. Overall models (age >= 50)
# ----------------------------------------------------------------------------
cat("\n--- Overall, age >= 50 ---\n")
report(svyglm(bmd ~ HbAA_sd + age + sex + race + bmi, design = design),
       "Model 2 (age+sex+race+BMI)", nrow(df))
report(svyglm(bmd ~ HbAA_sd + age + sex + race + bmi + smoke + alcohol +
                hypertension + diabetes + pir, design = design),
       "Model 3 (fully adjusted)", nrow(df))

# ----------------------------------------------------------------------------
# 4. Pre-specified subgroup: never smokers aged >= 65
# ----------------------------------------------------------------------------
sub <- subset(design, smoke == "Never" & age >= 65)
cur_sample <- "never_smokers_age_ge65"
cat(sprintf("\n--- Never smokers, age >= 65 (unweighted n = %d) ---\n",
            nrow(sub$variables)))
report(svyglm(bmd ~ HbAA_sd + age + sex + race + bmi, design = sub),
       "Model 2 (age+sex+race+BMI)", nrow(sub$variables))
report(svyglm(bmd ~ HbAA_sd + age + sex + race + bmi + alcohol + hypertension +
                diabetes + pir, design = sub),
       "Model 3 (fully adjusted)", nrow(sub$variables))

# ----------------------------------------------------------------------------
# 5. Write the effect estimates out as a machine-readable table
#    (set NHANES_OUT_CSV to a path; skipped when the variable is unset)
# ----------------------------------------------------------------------------
res <- do.call(rbind, results)
print(res, row.names = FALSE)

out_csv <- Sys.getenv("NHANES_OUT_CSV", unset = "")
if (nzchar(out_csv)) {
  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  write.csv(res, out_csv, row.names = FALSE)
  write.csv(data.frame(
    variable = c("HbAA_LBXACR_pmol_per_g_Hb", "BMD_DXXOFBMD_g_per_cm2"),
    n        = c(sum(!is.na(df$HbAA)), sum(!is.na(df$bmd))),
    mean     = c(mean(df$HbAA, na.rm = TRUE), mean(df$bmd, na.rm = TRUE)),
    sd       = c(sd(df$HbAA, na.rm = TRUE),   sd(df$bmd, na.rm = TRUE)),
    min      = c(min(df$HbAA, na.rm = TRUE),  min(df$bmd, na.rm = TRUE)),
    max      = c(max(df$HbAA, na.rm = TRUE),  max(df$bmd, na.rm = TRUE)),
    stringsAsFactors = FALSE
  ), sub("\\.csv$", "_descriptives.csv", out_csv), row.names = FALSE)
  cat("\nwrote ", out_csv, "\n", sep = "")
}
