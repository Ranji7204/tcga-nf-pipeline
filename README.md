🧬 TCGA RNA-seq Batch Effect Correction Pipeline (Nextflow)

A reproducible Nextflow (DSL2) pipeline for TCGA RNA-seq analysis, including:

TCGA barcode parsing

Batch effect detection (PCA + heatmap)

Conditional batch correction

Multi-region data merging

Automated QC plot generation

Docker & GitHub Codespaces compatible

📌 Overview

This pipeline processes TCGA STAR count RNA-seq data from multiple cancer regions (e.g., lung, cervical, head & neck) and performs:

Batch effect checking

PCA (before correction)

Sample clustering heatmap

Batch effect correction

Uses limma::removeBatchEffect

Automatically skipped if batch effect is insignificant

Region-wise outputs

Merged dataset analysis

Combined PCA

Combined heatmap

All results are saved in a structured results/ directory.

📁 Project Structure
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

📥 Input Requirements

Each region folder inside input/ must contain:

File	Description
*.unstranded.txt	STAR raw gene counts
*.xlsx	Sample metadata (cBioPortal / GDC)

Example:

input/lung/
├── lusc.unstranded.txt
└── normalized_metadata.xlsx

⚙️ Software Requirements
Local execution

Nextflow ≥ 25.04

Docker

Git

R packages (inside Docker)

edgeR

limma

sva

ggplot2

pheatmap

readxl

dplyr

🚀 Running the Pipeline
▶️ Local (Docker)
nextflow run main.nf

▶️ Resume after error
nextflow run main.nf -resume

📤 Outputs
Region-wise (results/<region>/)
File	Description
logCPM_normalized.txt	Normalized expression
batch_median.txt	Batch effect metric
metadata_with_batch_info.txt	Metadata + batch
PCA_before_batch_correction.png	QC PCA
heatmap_before_batch_correction.png	QC heatmap
logCPM_batch_corrected.txt	Corrected expression
PCA_after_batch_correction.png	Post-correction PCA
heatmap_after_batch_correction.png	Post-correction heatmap
Merged outputs (results/merged/)
File	Description
logCPM_merged_all_regions.txt	Combined expression
metadata_merged_all_regions.txt	Combined metadata
PCA_merged_regions.png	Combined PCA
heatmap_merged_regions.png	Combined heatmap
🧠 Pipeline Logic

Batch effect is assessed using TSS / barcode-derived batch proxy

If median R² < threshold, correction is skipped

PCA & heatmaps are always generated

Merging only happens when ≥2 regions are available

🐳 Docker Support
Build image
docker build -t tcga-nf-pipeline .

Run inside Docker
docker run -it \
  -v $PWD:/workspace \
  tcga-nf-pipeline \
  nextflow run main.nf

💻 GitHub Codespaces

Open repository

Click Code → Codespaces → Create

Run:

nextflow run main.nf


No local setup required.

🔐 Reproducibility

Fully containerized

Deterministic outputs

Resume-safe (-resume)

Portable across systems

📚 References

TCGA: https://portal.gdc.cancer.gov

edgeR: Robinson et al., Bioinformatics (2010)

limma: Ritchie et al., Nucleic Acids Research (2015)

👤 Author

Ranjith Gowda (Ranji7204)
Bioinformatics | RNA-seq | Cancer Genomics

GitHub: https://github.com/Ranji7204
