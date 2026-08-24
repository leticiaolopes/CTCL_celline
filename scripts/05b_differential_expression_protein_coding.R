# run all pairwise DE analyses using the protein coding DESeq2 object

local({
  previous_options <- options(ctcl.gene_scope = "protein_coding")
  on.exit(options(previous_options), add = TRUE)
  source("scripts/05_differential_expression.R", local = TRUE)
})
