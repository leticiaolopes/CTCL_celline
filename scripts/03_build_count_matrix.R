# combine featurecounts files into count matrix 
# export sample annotations

# setup
source("scripts/00_config.R")

library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tibble)

counts_dir <- paths$featurecounts_raw

count_files <- list.files(counts_dir, pattern = "\\.tabular$", full.names = TRUE)

# read & merge featurecounts files
read_featurecount <- function(file) {
  sample_id <- basename(file) |>
    str_remove("\\.tabular$") |>
    str_extract("^25060a\\d{3}_01")

  sample_name <- sample_metadata |>
    filter(sample_id == !!sample_id) |>
    pull(sample_name)

  if (length(sample_name) != 1) {
    stop("Could not uniquely match sample ID: ", sample_id)}

  x <- read_tsv(file, show_col_types = FALSE) |>
    select(Geneid, 2)

  names(x)[2] <- sample_name

  x}

counts_tbl <- count_files |>
  purrr::map(read_featurecount) |>
  purrr::reduce(dplyr::full_join, by = "Geneid")

counts_tbl <- counts_tbl |>
  select(Geneid, all_of(sample_metadata$sample_name))

dim(counts_tbl)
names(counts_tbl)

assigned_counts <- colSums(counts_tbl[, -1])

assigned_counts

sum(duplicated(counts_tbl$Geneid))
sum(is.na(counts_tbl))

count_matrix <- counts_tbl |>
  tibble::column_to_rownames("Geneid") |>
  as.matrix()

storage.mode(count_matrix) <- "integer"

coldata <- sample_metadata |>
  dplyr::select(sample_name, analysis_sample_name, original_cell_line, inferred_cell_line, identity_status,
    analysis_cell_line) |>
  tibble::column_to_rownames("sample_name")

# export count data & sample annotations
write.csv(counts_tbl, file.path(paths$featurecounts, "count_matrix_raw.csv"), row.names = FALSE)

write.csv(sample_metadata, file.path(paths$metadata, "sample_metadata_analysis.csv"), row.names = FALSE)

saveRDS(count_matrix, file.path(paths$featurecounts, "count_matrix_raw.rds"))

saveRDS(coldata, file.path(paths$featurecounts, "coldata.rds"))

dim(count_matrix)
coldata
