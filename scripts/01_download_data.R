
# Download GSE43458 from GEO and save a local copy.

required_packages <- c("GEOquery", "Biobase")

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

# Keep downloaded files outside Git.
out_dir <- file.path("data", "raw")

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

message("Downloading GSE43458 from GEO...")

gse_list <- GEOquery::getGEO(
  "GSE43458",
  GSEMatrix = TRUE,
  AnnotGPL = TRUE
)

if (length(gse_list) == 0) {
  stop("No expression data were returned.")
}

# Select GPL6244 if it is available.
platform_names <- names(gse_list)

platform_match <- if (!is.null(platform_names)) {
  grep(
    "GPL6244",
    platform_names,
    ignore.case = TRUE
  )
} else {
  integer(0)
}

if (length(platform_match) > 0) {
  gset <- gse_list[[platform_match[1]]]
} else if (length(gse_list) == 1) {
  gset <- gse_list[[1]]
} else {
  stop(
    "Could not select GPL6244. Available platforms: ",
    paste(platform_names, collapse = ", ")
  )
}

# Extract expression values and sample information.
expr <- Biobase::exprs(gset)
metadata <- as.data.frame(Biobase::pData(gset))

if (!is.numeric(expr)) {
  stop("The expression matrix is not numeric.")
}

if (nrow(expr) == 0 || ncol(expr) == 0) {
  stop("The expression matrix is empty.")
}

if (is.null(rownames(expr)) || is.null(colnames(expr))) {
  stop("Expression feature or sample IDs are missing.")
}

if (
  is.null(rownames(metadata)) ||
  !all(colnames(expr) %in% rownames(metadata))
) {
  stop("Expression samples do not match the metadata.")
}

# Keep the metadata in the same order as the expression matrix.
metadata <- metadata[
  colnames(expr),
  ,
  drop = FALSE
]

if (ncol(expr) != nrow(metadata)) {
  stop("Expression data and metadata do not match.")
}

if (anyNA(expr)) {
  warning("The expression matrix contains missing values.")
}

# Add the GEO sample ID as an explicit column.
metadata$sample_id <- rownames(metadata)

metadata <- metadata[
  ,
  c("sample_id", setdiff(names(metadata), "sample_id")),
  drop = FALSE
]

# Save the dataset and sample metadata.
saveRDS(
  gset,
  file = file.path(out_dir, "GSE43458.rds")
)

write.csv(
  metadata,
  file = file.path(out_dir, "sample_metadata.csv"),
  row.names = FALSE,
  na = ""
)

# Save a short summary of the downloaded data.
value_range <- range(expr, na.rm = TRUE)

dataset_summary <- data.frame(
  accession = "GSE43458",
  platform = if (
    !is.null(platform_names) &&
    length(platform_match) > 0
  ) {
    platform_names[platform_match[1]]
  } else {
    "Single returned platform"
  },
  feature_count = nrow(expr),
  sample_count = ncol(expr),
  min_expression = value_range[1],
  max_expression = value_range[2],
  retrieved_on = as.character(Sys.Date())
)

write.csv(
  dataset_summary,
  file = file.path(out_dir, "dataset_summary.csv"),
  row.names = FALSE
)

# Print the main checks.
cat("\nDownload completed.\n")
cat("Features:", nrow(expr), "\n")
cat("Samples:", ncol(expr), "\n")
cat("Metadata rows:", nrow(metadata), "\n")
cat(
  "Expression range:",
  paste(round(value_range, 3), collapse = " to "),
  "\n"
)

cat("\nMetadata columns:\n")
print(names(metadata))

if ("title" %in% names(metadata)) {
  cat("\nFirst five sample titles:\n")
  print(head(metadata$title, 5))
}
