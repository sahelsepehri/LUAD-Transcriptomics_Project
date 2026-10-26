
# Check sample metadata and assign tissue groups.

if (!requireNamespace("Biobase", quietly = TRUE)) {
  stop("Biobase is not installed.")
}

input_file <- file.path("data", "raw", "GSE43458.rds")
output_dir <- file.path("data", "processed")

if (!file.exists(input_file)) {
  stop("Dataset file not found. Run 01_download_data.R first.")
}

gset <- readRDS(input_file)

# Extract expression data and sample metadata.
expr <- Biobase::exprs(gset)
metadata <- as.data.frame(Biobase::pData(gset))

# Match sample IDs between the expression matrix and metadata.
if (
  is.null(rownames(metadata)) ||
  !all(colnames(expr) %in% rownames(metadata))
) {
  stop("Expression sample IDs do not match the metadata.")
}

metadata <- metadata[colnames(expr), , drop = FALSE]

# Use the sample title to identify the tissue type.
if (!"title" %in% names(metadata)) {
  stop("The title column was not found in the sample metadata.")
}

sample_title <- trimws(as.character(metadata[["title"]]))

# GEO titles distinguish tumor tissue from normal lung tissue.

is_tumor <- grepl(
  "^lung adenocarcinoma tissue([_[:space:]]|$)",
  sample_title,
  ignore.case = TRUE
)

is_normal <- grepl(
  "^normal lung tissue([_[:space:]]|$)",
  sample_title,
  ignore.case = TRUE
)


# Stop if any sample cannot be classified clearly.
unknown <- (is_tumor & is_normal) |
  (!is_tumor & !is_normal) |
  is.na(sample_title) |
  sample_title == ""

if (any(unknown)) {
  unclear_samples <- data.frame(
    sample_id = colnames(expr)[unknown],
    title = sample_title[unknown]
  )

  print(unclear_samples)

  stop(
    "Some samples could not be classified. ",
    "Check their titles before continuing."
  )
}

# Create the group variable.
group <- ifelse(is_tumor, "Tumor", "Normal")
group <- factor(group, levels = c("Normal", "Tumor"))

# Count samples in each group.
group_counts <- table(group)

cat("\nSample counts:\n")
print(group_counts)

# Check the expected dataset composition.
normal_count <- sum(group == "Normal")
tumor_count <- sum(group == "Tumor")

if (normal_count != 30 || tumor_count != 80) {
  stop(
    "Unexpected sample counts. Expected 30 Normal and 80 Tumor. ",
    "Check the sample titles before continuing."
  )
}

# Prepare a clean table for later analysis.
sample_groups <- data.frame(
  sample_id = colnames(expr),
  title = sample_title,
  group = as.character(group),
  stringsAsFactors = FALSE
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  sample_groups,
  file = file.path(output_dir, "sample_groups.csv"),
  row.names = FALSE,
  na = ""
)

cat("\nFirst six samples:\n")
print(head(sample_groups))

cat("\nSample group table saved successfully.\n")
