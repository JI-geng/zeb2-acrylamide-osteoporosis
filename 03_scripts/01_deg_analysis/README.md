# DEG analysis pipeline

This directory contains the differential expression analysis that produces the 169-gene
intersection of the manuscript.

```
GSE7158 (n=26)  ──┐
                   ├──> merge → ComBat → limma DEG → CTD intersection → 169 genes
GSE56815 (n=80) ──┘
```

## 1. Input files

| File | Source |
|---|---|
| `GSE7158_series_matrix.txt.gz` | GEO |
| `GSE56815_series_matrix.txt.gz` | GEO |
| `ctd_acrylamide_genes.txt` | CTD, acrylamide (CAS 79-06-1) |

## 2. Datasets and parameters

| Dataset | Platform | Samples | Tissue |
|---------|----------|---------|--------|
| GSE7158 | GPL570 (hgu133plus2) | 26 | peripheral blood monocytes |
| GSE56815 | GPL96 (hgu133a) | 80 | peripheral blood monocytes |

Total after merge: **106 samples** (54 high-BMD / 52 low-BMD).

| Parameter | Value |
|-----------|-------|
| Batch correction | ComBat |
| DEG method | limma |
| Significance | nominal P < 0.05 |

## 3. Results

| Metric | Value |
|--------|-------|
| Up-regulated | 446 genes |
| Down-regulated | 695 genes |
| Total DEGs | 1,141 genes |
| CTD ∩ DEG | 169 genes |

## 4. How to run

```bash
Rscript 00_geo_deg_combat_limma.R
```
