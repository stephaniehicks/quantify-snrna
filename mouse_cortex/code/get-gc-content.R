# get-gc-content.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Sep 22, 2021
#
# Run cQN to normalize and then run DE analysis.

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(SingleCellExperiment)
  library(EDASeq)
  library(tictoc)
})

# Read in SingleCellExperiment objects
normalize = FALSE
run_number = "all" # give run_number or "all" for all of them together
sce_ls = list()
sce_ls[["transcripts"]] = readRDS(here("mouse_cortex", "salmon_quants", "transcripts_pipeline", paste0("sce_", run_number, ".rds")))
sce_ls[["preandmrna"]] = readRDS(here("mouse_cortex", "salmon_quants", "preandmrna_pipeline", paste0("sce_", run_number, ".rds")))
sce_ls[["introncollapse"]] = readRDS(here("mouse_cortex", "salmon_quants", "introncollapse_pipeline", paste0("sce_", run_number, ".rds")))
sce_ls[["intronseparate"]] = readRDS(here("mouse_cortex", "salmon_quants", "intronseparate_pipeline", paste0("sce_", run_number, ".rds")))

pipeline = "preandmrna"

id = rownames(sce_ls[[pipeline]])
id = sub("\\..*", "", id)
org = "mmusculus_gene_ensembl"
tic()
gc = getGeneLengthAndGCContent(id, org, mode=c("biomart"))
toc()

# Save results
saveRDS(gc, here(paste0("./mouse_cortex/output/gc_content_", pipeline, ".rds")))