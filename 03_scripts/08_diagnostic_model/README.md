# Diagnostic model (nomogram, calibration, DCA, external validation)

| Script | Reproduces |
|---|---|
| `dm01_nomogram_discovery.R` | discovery cohort, n = 106: 5-gene logistic model (`rms::lrm`), nomogram, apparent AUC (0.856), single-gene AUCs, bootstrap optimism-corrected AUC, bootstrap calibration, decision curve |
| `dm02_external_validation.R` | external cohort, n = 40: AUC, calibration, decision curve |

## Inputs (not distributed — `02_derived_data/`)

| File | Used by |
|---|---|
| `expression_matrix_all_12487genes_106samples.csv` | dm01, dm02 |
| `group_labels_106samples.csv` (event label = `low BMD`) | dm01, dm02 |
| `external_validation_GSE13850_expression.csv` | dm02 |
| `external_validation_GSE13850_group_labels.csv` | dm02 |
