# 🧬 TCGA RNA-seq Batch Effect Correction Pipeline (Nextflow)

A reproducible **Nextflow (DSL2)** pipeline for **TCGA RNA-seq analysis**, including:

- TCGA barcode parsing  
- Batch effect detection (PCA + heatmap)  
- Conditional batch correction  
- Multi-region data merging  
- Automated QC plot generation  
- Docker & GitHub Codespaces compatible  

---

## 📌 Overview

This pipeline processes **TCGA STAR count RNA-seq data** from multiple cancer regions (e.g., **lung, cervical, head & neck**) and performs:

1. **Batch effect checking**
   - PCA (before correction)
   - Sample clustering heatmap
2. **Batch effect correction**
   - Uses `limma::removeBatchEffect`
   - Automatically skipped if batch effect is insignificant
3. **Region-wise outputs**
4. **Merged dataset analysis**
   - Combined PCA
   - Combined heatmap

---

## 📁 Project Structure

```
tcga-nf-pipeline/
├── input/
│   ├── cervical/
│   ├── headneck/
│   └── lung/
├── scripts/
│   ├── tcga_barcode_batch_pipeline_auto.R
│   ├── batch_correction.R
│   └── merge_all_files.R
├── envs/
│   └── r_tcga.yml
├── results/
├── main.nf
├── nextflow.config
├── Dockerfile
├── .dockerignore
└── README.md
```

---

## 📥 Input Requirements

Each region folder inside `input/` must contain:

| File | Description |
|---|---|
| `*.unstranded.txt` | STAR raw gene counts |
| `*.xlsx` | Sample metadata (cBioPortal / GDC) |

---

## ⚙️ Software Requirements

- Nextflow ≥ 25.04
- Docker
- Git

R packages (inside container):
edgeR, limma, sva, ggplot2, pheatmap, readxl, dplyr

---

## 🚀 Running the Pipeline

```bash
nextflow run main.nf
```

Resume:
```bash
nextflow run main.nf -resume
```

---

## 📤 Outputs

### Region-wise (`results/<region>/`)

- PCA_before_batch_correction.png  
- heatmap_before_batch_correction.png  
- PCA_after_batch_correction.png  
- heatmap_after_batch_correction.png  
- logCPM_normalized.txt  
- logCPM_batch_corrected.txt  
- metadata_with_batch_info.txt  

### Merged (`results/merged/`)

- PCA_merged_regions.png  
- heatmap_merged_regions.png  
- logCPM_merged_all_regions.txt  
- metadata_merged_all_regions.txt  

---

## 🐳 Docker

Build:
```bash
docker build -t tcga-nf-pipeline .
```

Run:
```bash
docker run -it -v $PWD:/workspace tcga-nf-pipeline nextflow run main.nf
```

---

## 👤 Author

**Ranjith Gowda (Ranji7204)**  
Bioinformatics | Cancer Genomics  

