# Analysis workflow and code — acrylamide exposure, ZEB2 and bone mineral density

---

## Contents

All analyses are secondary analyses of publicly available datasets. This repository provides:

1. **Analysis code** (`03_scripts/`)
2. **Full parameter record** (`00_parameters/`)
3. **Data provenance** (`01_data_provenance/`)

No data are deposited — all public source data must be downloaded by the user. No manuscript figures or figure-drawing code are included.

---

## Software

**R** — original runs: R 4.5.3; reproduction runs: R 4.6.1.

Key packages: `survey` 4.5, `caret` 7.0.1, `Seurat` 5.5.1, `limma` 3.68.5, `clusterProfiler` 4.20.0, `CellChat` 2.2.0.9001, `spacexr` 2.2.1.

**External tools:** AutoDock Vina 1.2.7, GROMACS 2026.2, gmx_MMPBSA v1.5.5, Uni-Dock (GPU portable build).

Full version table: `00_parameters/WORKFLOW_PARAMETERS.md`.

---

## Workflow

| Step | Description | Script |
|------|-------------|--------|
| 1 | GEO DEGs ∩ CTD → 169 genes | `03_scripts/01_deg_analysis/` |
| 2 | GO / KEGG / Hallmark enrichment | `03_scripts/02_enrichment/` |
| 3 | 10-algorithm ML → 5-gene signature | `03_scripts/03_machine_learning/` |
| 4 | Docking + MD / MM-GBSA | `03_scripts/04_docking_md/inputs/` |
| 5 | NHANES exposure–outcome association | `03_scripts/05_nhanes/` |
| 6 | Single-cell analysis (GSE147287) | `03_scripts/06_single_cell/` |
| 7 | Spatial deconvolution (GSE284089) | `03_scripts/07_spatial/` |
| 8 | Diagnostic nomogram + external validation | `03_scripts/08_diagnostic_model/` |

---

## Running

```bash
git clone <this-repo> && cd zeb2-acrylamide-osteoporosis
```

Each script's header lists the required input files and environment variables. All public datasets are listed in `01_data_provenance/DATA_SOURCES.md`.

---

## Layout

```
00_parameters/        workflow parameters + software versions
01_data_provenance/  public data sources, versions, access dates
03_scripts/          analysis code
  01_deg_analysis/    GEO → ComBat → limma → CTD intersection
  02_enrichment/      Enrichr re-analysis
  03_machine_learning/ 10-learners feature selection
  04_docking_md/      docking + MD input files
  05_nhanes/          NHANES survey-weighted analysis
  06_single_cell/     QC → clustering → annotation → CytoTRACE → CellChat → scTenifoldKnk
  07_spatial/         Visium QC → RCTD → spatial domains
  08_diagnostic_model/ nomogram + calibration + DCA + external validation
```

---

## Licence

MIT — see `LICENSE`. Upstream public datasets are not covered by this licence.
