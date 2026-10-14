# Workflow parameters

---

## Step 1 — NHANES exposure–outcome association

| Item | Value |
|------|-------|
| Cycles | 2005–2006, 2013–2014 |
| Exposure | HbAA (acrylamide haemoglobin adduct), `LBXACR`, pmol/g Hb |
| Outcome | Total femur BMD, `DXXOFBMD`, g/cm² |
| Population | age ≥50, n = 2,015 |
| Subgroup | never smokers aged ≥65, n = 433 |
| Survey design | `svydesign(ids=~psu, strata=~stratum, weights=~wt, nest=TRUE)` |
| Weights | 2005–2006: `WTMEC2YR/2`; 2013–2014: `WTSA2YR/2` |
| Models | `svyglm`, Model 2: age+sex+race+BMI; Model 3: +smoking+alcohol+hypertension+diabetes+PIR |
| Software | R 4.5.3, survey 4.4.8 |
| Script | `03_scripts/05_nhanes/NHANES_HbAA_BMD_analysis.R` |

---

## Step 2 — DEG and 169-gene intersection

| Item | Value |
|------|-------|
| Transcriptomic data | GSE7158 (n=26) + GSE56815 (n=80), merged n=106 |
| Tissue | peripheral blood monocytes |
| Groups | 54 high BMD / 52 low BMD |
| Method | limma + ComBat batch correction |
| DEG threshold | exploratory screen, nominal P < 0.05 |
| DEG counts | 446 up / 695 down (total 1,141) |
| Chemical–gene source | CTD, acrylamide query |
| Intersection | CTD ∩ DEG = 169 genes |
| Script | `03_scripts/01_deg_analysis/00_geo_deg_combat_limma.R` |

---

## Step 3 — Functional enrichment

| Item | Value |
|------|-------|
| Input | 169 genes |
| Libraries | GO BP/CC/MF, KEGG, MSigDB Hallmark |
| Tool | clusterProfiler 4.20.0, org.Hs.eg.db 3.20.0 |
| Correction | Benjamini–Hochberg |
| R version | 4.6.1 |
| Script | `03_scripts/02_enrichment/_enrich169.py` |

---

## Step 4 — Machine learning feature selection

| Item | Value |
|------|-------|
| Input | 169 genes × 106 samples |
| Split | stratified 70/30, train=75, test=31, seed=410 |
| Cross-validation | 5-fold, 1 repeat |
| Learners (10) | LASSO, RF, SVM-RFE, XGBoost, DT, GBM, GLM, NNET, KNN, C5.0 — one per file in `02_derived_data/per_learner_rankings/`; **the denominator of every consensus count below is exactly 10** |
| Feature shortlist | top 10 per learner |
| Consensus panel | intersection of top-10 from RF, XGBoost, GBM |
| Final genes | NAP1L2, ZEB2, MST1, PLCB3, ASNS |
| Software | R 4.5.3, caret 7.0.1, glmnet 4.1.10, randomForest 4.7.1.2, e1071 1.7.17, pROC 1.19.0.1 |
| Scripts | `03_scripts/03_machine_learning/` |

---

## Step 5 — Docking and molecular dynamics

| Item | Value |
|------|-------|
| Receptor | ZEB2 CZF construct (AlphaFold AF-O60315-F1) |
| Ligands | 5 natural products from ZINC22 |
| Docking tool | AutoDock Vina 1.2.7, exhaustiveness 32 |
| MD software | GROMACS 2026.2 |
| Protein force field | amber14sb |
| Ligand force field | gaff2 |
| Water model | tip3pfb |
| Box | truncated octahedron, 1.5 nm padding |
| Ion concentration | 0.15 M NaCl |
| Temperature | 310 K |
| Pressure | 1.0 bar |
| Time step | 2.0 fs |
| Production time | 100 ns |
| Equilibration | 200 ps NVT + 500 ps NPT |
| MM-GBSA | gmx_MMPBSA v1.5.5, igb=8, 500 frames |
| Input files | `03_scripts/04_docking_md/inputs/` |

---

## Step 6 — Single-cell analysis

| Item | Value |
|------|-------|
| Dataset | GSE147287 |
| Software | Seurat 5.5.1, harmony, SingleR/celldex 1.22.0 |
| Annotation method | SingleR, reference HPCA |
| Steps | QC → integration → clustering → annotation |
| Trajectory | CytoTRACE 2 (default parameters) |
| Cell–cell communication | CellChat 2.2.0.9001, mode=single, seed=2024 |
| Perturbation | scTenifoldKnk, target ZEB2, subset Tissue_stem_cells, 2000 HVG, 5 networks × 200 cells, seed=2024 |
| Scripts | `03_scripts/06_single_cell/` |

---

## Step 7 — Spatial transcriptomics

| Item | Value |
|------|-------|
| Dataset | GSE284089, 1 donor, 1 Visium section |
| Platform | 10x Visium |
| Min features | 200 |
| Max mitochondrial | 30% |
| Clustering | 15 PCs, 6 neighbors, resolution 0.5 |
| Spatial statistics | Moran's I, Geary's C, Getis-Ord |
| Deconvolution | spacexr 2.2.1 (RCTD) |
| Scripts | `03_scripts/07_spatial/` |

---

## Step 8 — Diagnostic model

| Item | Value |
|------|-------|
| Discovery cohort | n = 106 |
| Model | 5-gene multivariable logistic regression |
| Genes | NAP1L2, ZEB2, MST1, PLCB3, ASNS |
| Apparent AUC | 0.856 (95% CI 0.780–0.931) |
| External validation | GSE13850, n = 40, AUC = 0.708 |
| Software | R 4.6.1, rms 8.1.1, pROC 1.19.0.1, rmda 1.6 |
| Scripts | `03_scripts/08_diagnostic_model/` |
