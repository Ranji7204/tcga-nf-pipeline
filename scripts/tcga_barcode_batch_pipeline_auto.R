############################################
# TCGA AUTO-DETECT BARCODE + BATCH CHECK
# FINAL, NEXTFLOW-SAFE VERSION
############################################

#############################
# 0. Load libraries
#############################
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(pheatmap)
})

set.seed(123)

#############################
# 1. Auto-detect input files
#############################

meta_file <- list.files(pattern = "\\.xlsx$", ignore.case = TRUE)
expr_file <- list.files(pattern = "\\.txt$",  ignore.case = TRUE)

if (length(meta_file) != 1)
  stop("❌ Exactly ONE metadata (.xlsx) file must be present")

if (length(expr_file) != 1)
  stop("❌ Exactly ONE expression (.txt) file must be present")

message("📄 Metadata file: ", meta_file)
message("📊 Expression file: ", expr_file)

#############################
# 2. Read input data
#############################

meta <- read_excel(meta_file)

expr <- read.table(
  expr_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

#############################
# 3. Separate gene annotation & counts
#############################

gene_annot <- expr[, 1:3]
count_mat  <- expr[, -(1:3)]

#############################
# 4. Create 16-character TCGA IDs
#############################

if (!"SAMPLE_ID" %in% colnames(meta))
  stop("❌ Metadata must contain column: SAMPLE_ID")

meta <- meta %>%
  mutate(SAMPLE_ID_16 = substr(SAMPLE_ID, 1, 16))

expr_sample_16 <- substr(colnames(count_mat), 1, 16)

#############################
# 5. Match samples
#############################

common_samples <- intersect(meta$SAMPLE_ID_16, expr_sample_16)

if (length(common_samples) < 5)
  stop("❌ Too few common samples after matching")

meta_filt <- meta %>% filter(SAMPLE_ID_16 %in% common_samples)

count_mat_filt <- count_mat[, expr_sample_16 %in% common_samples]

count_mat_filt <- count_mat_filt[, match(
  meta_filt$SAMPLE_ID_16,
  substr(colnames(count_mat_filt), 1, 16)
)]

#############################
# 6. Restore full TCGA barcode
#############################

sample_map <- data.frame(
  SAMPLE_ID_16 = substr(colnames(count_mat_filt), 1, 16),
  COMPLETE_SAMPLE_ID = colnames(count_mat_filt),
  stringsAsFactors = FALSE
)

meta_final <- meta_filt %>%
  left_join(sample_map, by = "SAMPLE_ID_16")

stopifnot(all(meta_final$COMPLETE_SAMPLE_ID == colnames(count_mat_filt)))

#############################
# 7. Normalization (edgeR)
#############################

rownames(count_mat_filt) <- gene_annot[, 1]

dge <- DGEList(counts = count_mat_filt)
dge <- dge[filterByExpr(dge), , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

logCPM <- cpm(dge, log = TRUE, prior.count = 1)

#############################
# 8. Extract batch variables
#############################

barcode_parts <- strsplit(meta_final$COMPLETE_SAMPLE_ID, "-")

meta_final$TSS     <- sapply(barcode_parts, `[`, 2)
meta_final$PORTION <- sapply(barcode_parts, `[`, 5)
meta_final$PLATE   <- sapply(barcode_parts, `[`, 6)
meta_final$CENTER  <- sapply(barcode_parts, `[`, 7)

#############################
# 9. PCA BEFORE batch correction (QC)
#############################

pca <- prcomp(t(logCPM), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  PLATE = meta_final$PLATE
)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = PLATE)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(title = "PCA before batch correction")

ggsave(
  "PCA_before_batch_correction.png",
  p_pca,
  width = 6,
  height = 5
)

#############################
# 10. Heatmap BEFORE batch correction
#############################

annotation_col <- data.frame(PLATE = meta_final$PLATE)
rownames(annotation_col) <- meta_final$COMPLETE_SAMPLE_ID

var_genes <- apply(logCPM, 1, var)
top_genes <- names(sort(var_genes, decreasing = TRUE))[1:50]

pheatmap(
  logCPM[top_genes, ],
  scale = "row",
  annotation_col = annotation_col,
  show_colnames = FALSE,
  filename = "heatmap_before_batch_correction.png"
)

#############################
# 11. Quantify batch effect (median R²)
#############################

batch <- factor(meta_final$PLATE)

if (length(levels(batch)) < 2) {
  median_r2 <- 0
} else {
  r2 <- apply(logCPM, 1, function(g) {
    summary(lm(g ~ batch))$r.squared
  })
  median_r2 <- median(r2, na.rm = TRUE)
}

write.table(
  data.frame(median_r2 = median_r2),
  file = "batch_median.txt",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

#############################
# 12. Save outputs for downstream steps
#############################

write.table(
  logCPM,
  "logCPM_normalized.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

write.table(
  meta_final,
  "metadata_with_batch_info.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("✅ Batch check + QC completed successfully")
