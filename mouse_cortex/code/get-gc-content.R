# get-gc-content.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Feb 1, 2022
#
# Run cQN to normalize and then run DE analysis.

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(SingleCellExperiment)
  library(EDASeq)
  library(tictoc)
})

org = "mmusculus_gene_ensembl"

if(org == "mmusculus_gene_ensembl"){
  # Read in SingleCellExperiment objects
  normalize = FALSE
  run_number = "all" # give run_number or "all" for all of them together
  sce_ls = list()
  sce_ls[["transcripts"]] = readRDS(here("mouse_cortex", "salmon_quants", "transcripts_pipeline", paste0("sce_", run_number, ".rds")))
  sce_ls[["preandmrna"]] = readRDS(here("mouse_cortex", "salmon_quants", "preandmrna_pipeline", paste0("sce_", run_number, ".rds")))
  sce_ls[["introncollapse"]] = readRDS(here("mouse_cortex", "salmon_quants", "introncollapse_pipeline", paste0("sce_", run_number, ".rds")))
  sce_ls[["intronseparate"]] = readRDS(here("mouse_cortex", "salmon_quants", "intronseparate_pipeline", paste0("sce_", run_number, ".rds")))
  
  id = Reduce(union, list(rownames(sce_ls[["transcripts"]]), rownames(sce_ls[["preandmrna"]]), rownames(sce_ls[["introncollapse"]]), rownames(sce_ls[["intronseparate"]])))
  id = sub("\\..*", "", id)
  
  tic()
  gc = getGeneLengthAndGCContent(id, org, mode=c("biomart"))
  toc()
  
  # Save results
  saveRDS(gc, here("./mouse_cortex/output/gc_content.rds"))
} else if(org == "hsa"){
  # Load SCE objects
  load(here("mouse_cortex", "data", "SCE_AMY-n5_tran-etal.rda"))
  load(here("mouse_cortex", "data", "SCE_DLPFC-n3_tran-etal.rda"))
  load(here("mouse_cortex", "data", "SCE_HPC-n3_tran-etal.rda"))
  load(here("mouse_cortex", "data", "SCE_NAc-n8_tran-etal.rda"))
  load(here("mouse_cortex", "data", "SCE_sACC-n5_tran-etal.rda"))
  
  id = Reduce(union, list(rownames(sce.amy.tran), rownames(sce.dlpfc.tran),
                          rownames(sce.hpc.tran), rownames(sce.nac.tran),
                          rownames(sce.sacc.tran)))
  
  tic()
  gc = getGeneLengthAndGCContent(id, org, mode=c("biomart"))
  toc()
  
  # Save results
  saveRDS(gc, here("./mouse_cortex/output/gc_content_human.rds"))
}