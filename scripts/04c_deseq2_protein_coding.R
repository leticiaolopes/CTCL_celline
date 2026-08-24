# build a DESeq2 object restricted to GENCODE v50 protein-coding genes

local({
  previous_options <- options(ctcl.gene_scope = "protein_coding")
  on.exit(options(previous_options), add = TRUE)
  source("scripts/04_deseq2.R", local = TRUE)
})
