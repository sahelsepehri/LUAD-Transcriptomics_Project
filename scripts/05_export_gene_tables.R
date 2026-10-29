
# 05_export_gene_tables.R
# Create clean, reproducible gene tables from limma results.
# Run from the project root after 04_differential_expression.R.

input_file <- file.path(
  "results",
  "differential_expression",
  "DEG_all.csv"
)

out_dir <- dirname(input_file)

if (!file.exists(input_file)) {
  stop("DEG_all.csv not found. Run 04_differential_expression.R first.")
}

deg <- read.csv(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required <- c(
  "FeatureID", "GeneSymbol", "logFC",
  "P.Value", "adj.P.Val", "Status"
)

missing_cols <- setdiff(required, names(deg))

if (length(missing_cols) > 0) {
  stop(
    "Missing columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# Prepare columns for consistent filtering.
deg$GeneSymbol <- trimws(as.character(deg$GeneSymbol))
deg$Status_clean <- tolower(trimws(as.character(deg$Status)))

deg$logFC <- as.numeric(deg$logFC)
deg$P.Value <- as.numeric(deg$P.Value)
deg$adj.P.Val <- as.numeric(deg$adj.P.Val)

# 1. Check the four genes against the COMPLETE results.
genes_check <- c("BUB1", "MELK", "CCNA2", "SDK1")

symbol_upper <- toupper(
  ifelse(is.na(deg$GeneSymbol), "", deg$GeneSymbol)
)

gene_check <- deg[
  symbol_upper %in% genes_check,
  required,
  drop = FALSE
]

write.csv(
  gene_check,
  file.path(out_dir, "gene_check_BUB1_MELK_CCNA2_SDK1.csv"),
  row.names = FALSE
)

cat("\nCHECK AGAINST COMPLETE RESULTS\n")
cat("Matching rows:", nrow(gene_check), "\n")

if (nrow(gene_check) > 0) {
  print(gene_check, row.names = FALSE)
} else {
  cat("None of the four symbols were found in DEG_all.csv.\n")
}

# 2. Keep significant DEGs using the thresholds shown
# in the differential-expression plots.
eligible <- deg[
  !is.na(deg$GeneSymbol) &
    nzchar(deg$GeneSymbol) &
    is.finite(deg$logFC) &
    is.finite(deg$P.Value) &
    is.finite(deg$adj.P.Val) &
    deg$adj.P.Val < 0.05 &
    abs(deg$logFC) >= 1 &
    deg$Status_clean %in% c("upregulated", "downregulated"),
  ,
  drop = FALSE
]

eligible$symbol_key <- toupper(eligible$GeneSymbol)

# 3. Identify symbols that have conflicting directions.
up_symbols <- unique(
  eligible$symbol_key[
    eligible$Status_clean == "upregulated"
  ]
)

down_symbols <- unique(
  eligible$symbol_key[
    eligible$Status_clean == "downregulated"
  ]
)

conflict_symbols <- intersect(up_symbols, down_symbols)

conflict_rows <- eligible[
  eligible$symbol_key %in% conflict_symbols,
  required,
  drop = FALSE
]

write.csv(
  conflict_rows,
  file.path(out_dir, "gene_symbol_conflicts.csv"),
  row.names = FALSE
)

# Exclude conflicting symbols from the final gene-level tables.
eligible <- eligible[
  !eligible$symbol_key %in% conflict_symbols,
  ,
  drop = FALSE
]

# 4. Select one representative feature per gene symbol.
# Rank by adjusted P-value, then by absolute logFC.
make_top <- function(data, direction, n = 10) {
  
  x <- data[
    data$Status_clean == direction,
    ,
    drop = FALSE
  ]
  
  x <- x[
    order(x$adj.P.Val, -abs(x$logFC), x$GeneSymbol),
    ,
    drop = FALSE
  ]
  
  x <- x[!duplicated(x$symbol_key), , drop = FALSE]
  
  x <- head(x, n)
  
  x$Status <- if (direction == "upregulated") {
    "Upregulated"
  } else {
    "Downregulated"
  }
  
  x <- x[
    ,
    c("FeatureID", "GeneSymbol", "logFC",
      "P.Value", "adj.P.Val", "Status"),
    drop = FALSE
  ]
  
  rownames(x) <- NULL
  x
}

top_up <- make_top(eligible, "upregulated", 10)
top_down <- make_top(eligible, "downregulated", 10)

# 5. Export separate tables and a combined table.
write.csv(
  top_up,
  file.path(out_dir, "top10_upregulated_genes.csv"),
  row.names = FALSE
)

write.csv(
  top_down,
  file.path(out_dir, "top10_downregulated_genes.csv"),
  row.names = FALSE
)

top_combined <- rbind(top_up, top_down)

write.csv(
  top_combined,
  file.path(out_dir, "top10_gene_tables.csv"),
  row.names = FALSE
)

# 6. Generate Markdown tables for the scientific report.
format_p <- function(x) {
  vapply(
    x,
    function(v) formatC(v, format = "e", digits = 2),
    character(1)
  )
}

markdown_table <- function(x) {
  
  if (nrow(x) == 0) {
    return("_No genes met the selection criteria._")
  }
  
  display <- data.frame(
    Gene = x$GeneSymbol,
    log2FC = sprintf("%.3f", x$logFC),
    `P-value` = format_p(x$P.Value),
    `Adjusted P-value` = format_p(x$adj.P.Val),
    check.names = FALSE
  )
  
  header <- paste0(
    "| ", paste(names(display), collapse = " | "), " |"
  )
  
  separator <- paste0(
    "| ", paste(rep("---", ncol(display)), collapse = " | "), " |"
  )
  
  rows <- apply(
    as.matrix(display),
    1,
    function(r) paste0("| ", paste(r, collapse = " | "), " |")
  )
  
  c(header, separator, rows)
}

report_lines <- c(
  "# Differential Expression: Top Gene Tables",
  "",
  paste(
    "Selection criteria: adjusted P-value < 0.05 and",
    "|log2FC| >= 1, consistent with the thresholds used",
    "for differential-expression classification."
  ),
  "",
  paste(
    "Genes are ranked by adjusted P-value, with absolute",
    "log2FC used to break ties. When multiple features map",
    "to the same gene symbol, only the feature with the",
    "smallest adjusted P-value is retained. Symbols found",
    "in both directions are excluded pending review."
  ),
  "",
  "## Top 10 Upregulated Genes",
  "",
  markdown_table(top_up),
  "",
  "## Top 10 Downregulated Genes",
  "",
  markdown_table(top_down),
  "",
  "## Interpretation",
  "",
  paste(
    "Positive log2FC indicates higher expression in Tumor",
    "than Normal samples; negative log2FC indicates lower",
    "expression in Tumor samples. These genes are candidates",
    "for downstream enrichment and protein interaction analyses."
  )
)

writeLines(
  report_lines,
  file.path(out_dir, "top10_gene_tables.md"),
  useBytes = TRUE
)

# 7. Report the outputs.
cat("\nTABLE GENERATION COMPLETED\n")
cat("Unique upregulated genes:", nrow(top_up), "\n")
cat("Unique downregulated genes:", nrow(top_down), "\n")
cat("Conflicting gene symbols:", length(conflict_symbols), "\n")
cat("Output directory:", out_dir, "\n")

cat("\nTOP UPREGULATED GENES\n")
print(top_up, row.names = FALSE)

cat("\nTOP DOWNREGULATED GENES\n")
print(top_down, row.names = FALSE)
