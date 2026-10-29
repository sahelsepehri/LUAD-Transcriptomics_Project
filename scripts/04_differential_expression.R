
# ============================================================
# 04_differential_expression.R
# LUAD Gene Expression Analysis
# Tumor vs Normal
# ============================================================

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

# 1. Load the existing dataset and sample groups
data_file <- file.path("data", "raw", "GSE43458.rds")
groups_file <- file.path(
  "data", "processed", "sample_groups.csv"
)

if (!file.exists(data_file)) {
  stop("Dataset not found: ", data_file)
}

if (!file.exists(groups_file)) {
  stop("Sample groups not found: ", groups_file)
}

gset <- readRDS(data_file)
expr <- as.matrix(Biobase::exprs(gset))
storage.mode(expr) <- "numeric"

samples <- read.csv(
  groups_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 2. Validate sample metadata
required_columns <- c("sample_id", "group")

if (!all(required_columns %in% names(samples))) {
  stop("sample_groups.csv must contain sample_id and group.")
}

samples$sample_id <- trimws(as.character(samples$sample_id))
samples$group <- trimws(as.character(samples$group))

sample_index <- match(colnames(expr), samples$sample_id)

if (anyNA(sample_index)) {
  stop("Some expression samples are missing from sample_groups.csv.")
}

samples <- samples[sample_index, , drop = FALSE]

group <- factor(
  samples$group,
  levels = c("Normal", "Tumor")
)

if (anyNA(group)) {
  stop("Unexpected or missing sample group labels.")
}

if (sum(group == "Normal") != 30 ||
    sum(group == "Tumor") != 80) {
  stop("Expected 30 Normal and 80 Tumor samples. Check sample groups.")
}

cat("Sample counts:\n")
print(table(group))

# 3. Remove features containing non-finite values
finite_features <- apply(is.finite(expr), 1, all)
expr_use <- expr[finite_features, , drop = FALSE]

if (nrow(expr_use) < 2) {
  stop("Not enough valid expression features.")
}

cat("Original features:", nrow(expr), "\n")
cat("Features retained:", nrow(expr_use), "\n")

# 4. Differential expression using limma
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

contrast_matrix <- limma::makeContrasts(
  TumorVsNormal = Tumor - Normal,
  levels = design
)

fit <- limma::lmFit(expr_use, design)
fit <- limma::contrasts.fit(fit, contrast_matrix)
fit <- limma::eBayes(fit)

results <- limma::topTable(
  fit,
  coef = "TumorVsNormal",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# 5. Add gene symbols if the existing dataset contains them
feature_ids <- rownames(results)
gene_symbols <- rep(NA_character_, length(feature_ids))

feature_data <- tryCatch(
  Biobase::fData(gset),
  error = function(e) NULL
)

if (!is.null(feature_data) &&
    nrow(feature_data) > 0 &&
    !is.null(rownames(feature_data))) {

  symbol_columns <- grep(
    "gene[._ ]?symbol|symbol",
    names(feature_data),
    ignore.case = TRUE,
    value = TRUE
  )

  if (length(symbol_columns) > 0) {
    feature_index <- match(
      feature_ids,
      rownames(feature_data)
    )

    gene_symbols <- as.character(
      feature_data[[symbol_columns[1]]][feature_index]
    )

    gene_symbols[
      is.na(gene_symbols) |
      trimws(gene_symbols) %in% c("", "---", "NA")
    ] <- NA_character_
  }
}

results$FeatureID <- feature_ids
results$GeneSymbol <- gene_symbols

# 6. Classify upregulated and downregulated features
# logFC = Tumor minus Normal
results$Status <- "Not significant"

significant <- !is.na(results$adj.P.Val) &
  results$adj.P.Val < 0.05

results$Status[
  significant & results$logFC >= 1
] <- "Upregulated"

results$Status[
  significant & results$logFC <= -1
] <- "Downregulated"

results <- results[
  ,
  c(
    "FeatureID", "GeneSymbol", "logFC", "AveExpr",
    "t", "P.Value", "adj.P.Val", "B", "Status"
  )
]

# 7. Save result tables
results_dir <- file.path(
  "results", "differential_expression"
)

image_dir <- "Images"

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  image_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

upregulated <- results[
  results$Status == "Upregulated",
  ,
  drop = FALSE
]

downregulated <- results[
  results$Status == "Downregulated",
  ,
  drop = FALSE
]

upregulated <- upregulated[
  order(-upregulated$logFC, upregulated$adj.P.Val),
  ,
  drop = FALSE
]

downregulated <- downregulated[
  order(downregulated$logFC, downregulated$adj.P.Val),
  ,
  drop = FALSE
]

write.csv(
  results,
  file.path(results_dir, "DEG_all.csv"),
  row.names = FALSE
)

write.csv(
  upregulated,
  file.path(results_dir, "DEG_upregulated.csv"),
  row.names = FALSE
)

write.csv(
  downregulated,
  file.path(results_dir, "DEG_downregulated.csv"),
  row.names = FALSE
)

write.csv(
  head(upregulated, 20),
  file.path(results_dir, "top20_upregulated.csv"),
  row.names = FALSE
)

write.csv(
  head(downregulated, 20),
  file.path(results_dir, "top20_downregulated.csv"),
  row.names = FALSE
)

# 8. Define plot colors
point_colors <- ifelse(
  results$Status == "Upregulated",
  "#D55E00",
  ifelse(
    results$Status == "Downregulated",
    "#0072B2",
    "grey70"
  )
)

# 9. Volcano plot
padj_safe <- pmax(
  results$adj.P.Val,
  .Machine$double.xmin
)

neg_log10_padj <- -log10(padj_safe)
neg_log10_padj[is.na(results$adj.P.Val)] <- NA_real_

png(
  file.path(image_dir, "volcano_plot.png"),
  width = 1800,
  height = 1400,
  res = 180
)

par(mar = c(5, 5, 4, 2))

plot(
  results$logFC,
  neg_log10_padj,
  pch = 20,
  cex = 0.65,
  col = point_colors,
  xlab = "log2 Fold Change (Tumor - Normal)",
  ylab = "-log10 Adjusted P-value",
  main = "Volcano Plot: Tumor vs Normal"
)

abline(v = c(-1, 1), h = -log10(0.05), lty = 2)

legend(
  "topright",
  legend = c(
    "Downregulated",
    "Not significant",
    "Upregulated"
  ),
  col = c("#0072B2", "grey70", "#D55E00"),
  pch = 20,
  bty = "n"
)

dev.off()

# 10. MA plot
png(
  file.path(image_dir, "MA_plot.png"),
  width = 1800,
  height = 1400,
  res = 180
)

par(mar = c(5, 5, 4, 2))

plot(
  results$AveExpr,
  results$logFC,
  pch = 20,
  cex = 0.65,
  col = point_colors,
  xlab = "Average Expression",
  ylab = "log2 Fold Change (Tumor - Normal)",
  main = "MA Plot: Tumor vs Normal"
)

abline(h = c(-1, 1), lty = 2)
abline(h = 0, lty = 1)

dev.off()

# 11. Heatmap of the top-ranked features
ranked_ids <- results$FeatureID[
  order(results$adj.P.Val, results$P.Value, na.last = NA)
]

top_ids <- head(ranked_ids, 50)

heat_matrix <- expr_use[
  top_ids,
  ,
  drop = FALSE
]

if (nrow(heat_matrix) > 0) {
  row_sd <- apply(heat_matrix, 1, sd)
  variable_rows <- is.finite(row_sd) & row_sd > 0

  heat_matrix <- heat_matrix[
    variable_rows,
    ,
    drop = FALSE
  ]
}

png(
  file.path(image_dir, "heatmap_top_DEGs.png"),
  width = 1800,
  height = 1600,
  res = 180
)

if (nrow(heat_matrix) >= 2) {

  sample_colors <- ifelse(
    group == "Normal",
    "#377EB8",
    "#E76F51"
  )

  heatmap(
    heat_matrix,
    scale = "row",
    ColSideColors = sample_colors,
    labRow = FALSE,
    labCol = FALSE,
    margins = c(5, 2),
    col = colorRampPalette(
      c("navy", "white", "firebrick3")
    )(100),
    main = "Top 50 Ranked Features: Tumor vs Normal"
  )

} else {

  plot.new()
  title("Heatmap: Tumor vs Normal")
  text(
    0.5, 0.5,
    "Insufficient variable features for heatmap"
  )
}

dev.off()

# 12. Report summary
cat("\nDifferential expression analysis completed.\n")
cat("Features tested:", nrow(results), "\n")
cat("Upregulated:", nrow(upregulated), "\n")
cat("Downregulated:", nrow(downregulated), "\n")
cat(
  "Rows with gene symbols:",
  sum(!is.na(results$GeneSymbol)),
  "\n"
)
cat("Results directory:", results_dir, "\n")
cat("Plot directory:", image_dir, "\n")
