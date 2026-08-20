# multiqc output

# setup
source("scripts/00_config.R")

library(readr)
library(dplyr)
library(stringr)
library(tidyr)

# load multiqc 
multiqc <- read_tsv(file.path(paths$qc, "multiqc_data.tabular"), show_col_types = FALSE)

glimpse(multiqc)
names(multiqc)
head(multiqc)

multiqc_parsed <- multiqc |>
  mutate(sample_id = str_extract(Sample, "^25060a[0-9]{3}_01"), sequencing_id = str_extract(Sample,
    "S[0-9]+"), lane = str_extract(Sample, "L[0-9]{3}"), read = str_extract(Sample, "R[12]"))

sample_metadata <- read_csv(file.path(paths$metadata, "sample_metadata.csv"), show_col_types = FALSE)

# attach sample metadata 
qc_metadata <- multiqc_parsed |>
  left_join(sample_metadata, by = "sample_id")

qc_summary <- qc_metadata |>
  select(sample_id, sample_name, cell_line, sequencing_id, lane, read, `falco-total_sequences`, `falco-percent_gc`,
    `falco-percent_duplicates`, `falco-avg_sequence_length`, `falco-percent_fails`) |>
  arrange(sample_id, lane, read)

qc_summary |>
  print(n = Inf, width = Inf)

write_csv(qc_summary, file.path(paths$qc, "qc_summary.csv"))

# summarize observed sequencing depth
sample_depth <- qc_summary |>
  filter(read == "R1") |>
  group_by(sample_id, sample_name, cell_line) |>
  summarise(read_pairs_million = sum(`falco-total_sequences`), .groups = "drop")

sample_depth |>
  print(n = Inf)

write_csv(sample_depth, file.path(paths$qc, "sample_depth.csv"))

# compare observed & reported read counts
depth_check <- qc_summary |>
  filter(read == "R1") |>
  group_by(sample_id, sample_name, cell_line) |>
  summarise(read_pairs_million = sum(`falco-total_sequences`), .groups = "drop") |>
  left_join(sample_metadata |>
    select(sample_id, facility_read_count), by = "sample_id") |>
  mutate(observed_reads = read_pairs_million * 2 * 1e+06, difference_reads = observed_reads - facility_read_count,
    difference_percent = 100 * difference_reads/facility_read_count)

depth_check |>
  print(n = Inf, width = Inf)

write_csv(depth_check, file.path(paths$qc, "depth_check.csv"))
