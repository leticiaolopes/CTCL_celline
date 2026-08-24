# Generate the top-50 heatmap from protein-coding one-vs-rest DE results and
# the protein-coding VST matrix.

local({
  previous_options <- options(ctcl.gene_scope = "protein_coding")
  on.exit(options(previous_options), add = TRUE)
  source("scripts/09_heatmap.R", local = TRUE)
})
