############################################
# Batch correction script (NEXTFLOW SAFE)
############################################

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(pheatmap)
})

message("▶ Starting batch correction")

#############################
# 1. Read median R2
#############################

threshold <- 0.1
median_file <- "batch_median.txt"

if (!file.exists(median_file)) {
  stop("❌ batch_median.txt not found")
}

median_r2 <- as.numeric(tail(readLines(median_file), 1))
if (is.na(median_r2)) stop("❌ Cannot parse median R2")

message("Median R2 = ", median_r2)

#############################
# 2. Read inputs
#############################

meta <- read.table(
  "metadata_with_batch_info.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

expr_norm <- read.table(
  "logCPM_normalized.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

expr_norm <- as.matrix(expr_norm)

stopifnot(all(colnames(expr_norm) == meta$COMPLETE_SAMPLE_ID))

#############################
# 3. Decide whether to correct
#############################

do_correction <- TRUE

if (median_r2 < threshold) {
  message("✔ Batch effect below threshold — skipping correction")
  do_correction <- FALSE
}

if (!"PLATE" %in% colnames(meta)) {
  message("✔ PLATE column missing — skipping correction")
  do_correction <- FALSE
}

batch <- factor(meta$PLATE)
if (length(levels(batch)) < 2) {
  message("✔ Only one batch — skipping correction")
  do_correction <- FALSE
}

#############################
# 4. Apply correction OR pass-through
#############################

if (do_correction) {
  
  message("▶ Applying limma batch correction")
  
  design <- model.matrix(~1, data = meta)
  
  expr_final <- removeBatchEffect(
    expr_norm,
    batch  = batch,
    design = design
  )
  
  meta_final <- meta
  
} else {
  
  expr_final <- expr_norm
  meta_final <- meta
}

#############################
# 5. Save expression & metadata (ALWAYS)
#############################

write.table(
  expr_final,
  "logCPM_batch_corrected.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

write.table(
  meta_final,
  "metadata_batch_corrected.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

#############################
# 6. PCA after batch correction (ALWAYS)
#############################

pca <- prcomp(t(expr_final), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Batch = batch
)

p <- ggplot(pca_df, aes(PC1, PC2, color = Batch)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(title = "PCA after batch correction")

ggsave(
  "PCA_after_batch_correction.png",
  p,
  width = 6,
  height = 5
)

#############################
# 7. Heatmap QC (ALWAYS)
#############################

annotation_col <- data.frame(Batch = batch)
rownames(annotation_col) <- meta$COMPLETE_SAMPLE_ID

pheatmap(
  expr_final[1:min(500, nrow(expr_final)), ],
  scale = "row",
  annotation_col = annotation_col,
  show_colnames = FALSE,
  filename = "heatmap_after_batch_correction.png"
)

message("✅ Batch correction step completed successfully")
