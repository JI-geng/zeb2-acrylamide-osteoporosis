# NHANES: HbAA vs femur BMD

## 1. Data files

| Cycle | Exposure | Bone | Demography | Smoking | Alcohol | BP | Diabetes |
|---|---|---|---|---|---|---|---|
| 2005–2006 | `AMDGYD_D.XPT` | `DXXFEM_D.XPT` | `DEMO_D` | `SMQ_D` | `ALQ_D` | `BPX_D` | `DIQ_D` |
| 2013–2014 | `AMDGYD_H.XPT` | `DXXFEM_H.XPT` | `DEMO_H` | `SMQ_H` | `ALQ_H` | `BPX_H` | `DIQ_H` |

## 2. Run

```bash
export NHANES_DATA_DIR=/path/to/nhanes_xpt
Rscript NHANES_HbAA_BMD_analysis.R
```

Requirements: R ≥ 4.1, `survey`, `haven`, `dplyr`.

## 3. Outputs

`nhanes_effect_estimates.csv` — one row per model.
