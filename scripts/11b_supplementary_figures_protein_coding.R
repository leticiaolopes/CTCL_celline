# Assemble the five supplementary figures from protein-coding panel PDFs.
# Outputs are written to figures/supplementary_protein_coding.

local({
  previous_options <- options(ctcl.gene_scope = "protein_coding")
  on.exit(options(previous_options), add = TRUE)
  source("scripts/11_supplementary_figures.R", local = TRUE)
})
