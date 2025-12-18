############################################
# Merge all regions (FINAL, FIXED)
############################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) > 0) args[1] else "results"

message("▶ Merging all regions from: ", results_dir)

# -----------------------------------------
# 1. Detect valid region directories
# -----------------------------------------
region_dirs <- list.dirs(results_dir, recursive = FALSE, full.names = TRUE)

valid_regions <- region_dirs[
  file.exists(file.path(region_dirs, "logCPM_batch_corrected.txt")) &
    file.exists(file.path(region_dirs, "metadata_batch_corrected.txt"))
]

message("✔ Regions detected: ", length(valid_regions))
message("✔ Regions: ", paste(basename(valid_regions), collapse = ", "))

if (length(valid_regions) < 2) {
  stop("❌ Need at least two regions to merge")
}

# -----------------------------------------
# 2. Read expression matrices
# -----------------------------------------
expr_list <- lapply(valid_regions, function(rdir) {
  x <- read.table(
    file.path(rdir, "logCPM_batch_corrected.txt"),
    header = TRUE,
    sep = "\t",
    check.names = FALSE
  )
  as.matrix(x)
})

# Intersect genes
common_genes <- Reduce(intersect, lapply(expr_list, rownames))
expr_list <- lapply(expr_list, function(x) x[common_genes, ])

expr_all <- do.call(cbind, expr_list)

# -----------------------------------------
# 3. Read metadata + ADD REGION
# -----------------------------------------
meta_list <- lapply(valid_regions, function(rdir) {
  
  m <- read.table(
    file.path(rdir, "metadata_batch_corrected.txt"),
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  # ADD REGION FROM DIRECTORY NAME
  m$REGION <- basename(rdir)
  m
})

meta_all <- bind_rows(meta_list)

# Ensure sample order matches expression
meta_all <- meta_all[
  match(colnames(expr_all), meta_all$COMPLETE_SAMPLE_ID),
]

stopifnot(all(meta_all$COMPLETE_SAMPLE_ID == colnames(expr_all)))

# -----------------------------------------
# 4. PCA
# -----------------------------------------
pca <- prcomp(t(expr_all), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  REGION = meta_all$REGION
)

p <- ggplot(pca_df, aes(PC1, PC2, color = REGION)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(title = "PCA of merged regions")

ggsave("PCA_merged_regions.png", p, width = 6, height = 5)

# -----------------------------------------
# 5. Heatmap
# -----------------------------------------
var_genes <- apply(expr_all, 1, var)
top_genes <- names(sort(var_genes, decreasing = TRUE))[1:50]

annotation_col <- data.frame(REGION = meta_all$REGION)
rownames(annotation_col) <- meta_all$COMPLETE_SAMPLE_ID

pheatmap(
  expr_all[top_genes, ],
  scale = "row",
  annotation_col = annotation_col,
  show_colnames = FALSE,
  filename = "heatmap_merged_regions.png"
)

# -----------------------------------------
# 6. Save outputs
# -----------------------------------------
write.table(
  expr_all,
  "logCPM_merged_all_regions.txt",
  sep = "\t",
  quote = FALSE
)

write.table(
  meta_all,
  "metadata_merged_all_regions.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("✅ Merge completed successfully")
