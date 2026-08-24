# Replot supplementary panels from cached protein-coding RDS/CSV objects.

local({
  previous_options <- options(ctcl.gene_scope = "protein_coding")
  on.exit(options(previous_options), add = TRUE)
  source("scripts/10b_replot_supplementary_panels.R", local = TRUE)
})
