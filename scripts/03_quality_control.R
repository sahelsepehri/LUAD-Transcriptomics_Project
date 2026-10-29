
# Check expression distributions and sample patterns.

required_packages <- c("Biobase", "limma")

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Missing packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

data_file <- file.path("data", "raw", "GSE43458.rds")
groups_file <- file.path("data", "processed", "sample_groups.csv")

if (!file.exists(data_file)) {
  stop("Dataset not found. Run 01_download_data.R first.")
}

if (!file.exists(groups_file)) {
  stop("Sample groups not found. Run 02_sample_groups.R first.")
}

gset <- readRDS(data_file)
expr <- Biobase::exprs(gset)

samples <- read.csv(
  groups_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Match the sample information to the expression matrix.
if (!all(c("sample_id", "group") %in% names(samples))) {
  stop("The sample group file is missing required columns.")
}

sample_index <- match(colnames(expr), samples$sample_id)

if (anyNA(sample_index)) {
  stop("Some expression samples are missing from the group table.")
}

samples <- samples[sample_index, , drop = FALSE]

group <- factor(
  samples$group,
  levels = c("Normal", "Tumor")
)

if (anyNA(group)) {
  stop("Some samples have missing or unexpected group labels.")
}

if (sum(group == "Normal") != 30 ||
    sum(group == "Tumor") != 80) {
  stop("Sample counts do not match the expected group sizes.")
}

cat("Sample counts:\n")
print(table(group))

# Remove features with non-finite values from the plots.
finite_features <- apply(is.finite(expr), 1, all)
expr_qc <- expr[finite_features, , drop = FALSE]

if (nrow(expr_qc) < 2) {
  stop("Not enough valid features for quality control.")
}

cat("\nFeatures in original matrix:", nrow(expr), "\n")
cat("Features used for QC plots:", nrow(expr_qc), "\n")
cat("Non-finite expression values:", sum(!is.finite(expr)), "\n")

# Use simple, consistent colors for both groups.
group_colors <- c(
  Normal = "#3B82F6",
  Tumor = "#E76F51"
)

sample_colors <- unname(group_colors[as.character(group)])

# Create output folders.
plot_dir <- file.path("Images", "QualityControl")
output_dir <- file.path("data", "processed")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Boxplot: compare expression distributions across samples.
png(
  file.path(plot_dir, "expression_boxplot.png"),
  width = 2400,
  height = 1200,
  res = 180
)

boxplot(
  expr_qc,
  outline = FALSE,
  xaxt = "n",
  col = sample_colors,
  border = "grey35",
  main = "Expression Distribution Across Samples",
  xlab = "Samples (sample IDs omitted)",
  ylab = "Expression value"
)

legend(
  "topright",
  legend = names(group_colors),
  fill = group_colors,
  border = NA,
  bty = "n"
)

dev.off()

# 2. Density plot: compare the shape of expression distributions.
png(
  file.path(plot_dir, "expression_density.png"),
  width = 2400,
  height = 1200,
  res = 180
)

limma::plotDensities(
  expr_qc,
  group = group,
  col = unname(group_colors),
  legend = "topright",
  main = "Density of Expression Values",
  xlab = "Expression value"
)

dev.off()

# 3. PCA: use the 1,000 most variable valid features.
feature_variance <- apply(expr_qc, 1, var)

variable_features <- is.finite(feature_variance) &
  feature_variance > 0

expr_variable <- expr_qc[variable_features, , drop = FALSE]
feature_variance <- feature_variance[variable_features]

if (nrow(expr_variable) < 2) {
  stop("Not enough variable features for PCA.")
}

top_n <- min(1000, nrow(expr_variable))

top_index <- order(
  feature_variance,
  decreasing = TRUE
)[seq_len(top_n)]

expr_top <- expr_variable[top_index, , drop = FALSE]

# Samples are rows for PCA.
pca <- prcomp(
  t(expr_top),
  center = TRUE,
  scale. = FALSE
)

variance_percent <- 100 * pca$sdev^2 / sum(pca$sdev^2)

png(
  file.path(plot_dir, "pca_samples.png"),
  width = 1600,
  height = 1200,
  res = 180
)

plot(
  pca$x[, 1],
  pca$x[, 2],
  col = sample_colors,
  pch = ifelse(group == "Normal", 16, 17),
  cex = 1.1,
  xlab = sprintf("PC1 (%.1f%%)", variance_percent[1]),
  ylab = sprintf("PC2 (%.1f%%)", variance_percent[2]),
  main = "PCA of Tumor and Normal Samples"
)

legend(
  "topright",
  legend = levels(group),
  col = unname(group_colors[levels(group)]),
  pch = c(16, 17),
  bty = "n"
)

dev.off()

# Save sample-level QC measurements.
sample_summary <- data.frame(
  sample_id = colnames(expr),
  group = as.character(group),
  median_expression = apply(expr, 2, median, na.rm = TRUE),
  iqr_expression = apply(expr, 2, IQR, na.rm = TRUE),
  non_finite_values = colSums(!is.finite(expr))
)

write.csv(
  sample_summary,
  file.path(output_dir, "qc_sample_summary.csv"),
  row.names = FALSE
)

# Save a short overview of the QC step.
qc_overview <- data.frame(
  sample_count = ncol(expr),
  feature_count = nrow(expr),
  features_used_for_qc = nrow(expr_qc),
  features_used_for_pca = top_n,
  non_finite_values = sum(!is.finite(expr)),
  pc1_variance_percent = variance_percent[1],
  pc2_variance_percent = variance_percent[2]
)

write.csv(
  qc_overview,
  file.path(output_dir, "qc_overview.csv"),
  row.names = FALSE
)

cat("\nPCA variance explained:\n")
cat("PC1:", round(variance_percent[1], 2), "%\n")
cat("PC2:", round(variance_percent[2], 2), "%\n")

cat("\nPlots saved in:", plot_dir, "\n")
cat("QC summaries saved in:", output_dir, "\n")
