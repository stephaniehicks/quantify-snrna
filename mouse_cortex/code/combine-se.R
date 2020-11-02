# combine-se.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Nov 2, 2020
#
# Combine SummarizedExperiment objects into one 

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(SummarizedExperiment)
})

# pipeline = "transcripts" # "transcripts" or "preandmrna" or "introncollapse/separate"
for(pipeline in c("transcripts", "preandmrna", "introncollapse", "intronseparate")){
  cortex_names = c("cortex1", "cortex2")
  
  # Get t2g metadata for each gene
  tx2gene = readRDS(here("mouse_cortex", "salmon_files", "gencode.vM25.transcripts.tx2gene.mouse.rds"))
  
  # Read in individual rds files
  se_ls = list()
  for(cortex_name in cortex_names){
    se = readRDS(here("mouse_cortex", "salmon_quants", paste0(pipeline, "_pipeline"), paste0("se_", cortex_name, ".rds")))
    
    # Add colData
    column_names = colnames(se)
    colData(se) = DataFrame(cortex = rep(cortex_name, length(column_names)))
    colnames(se) = column_names
    
    # Add rowData
    match_rows = match(rownames(se), tx2gene$gene_id)
    t2g = tx2gene[match_rows, ]
    rowData(se) = t2g
    
    se_ls[[cortex_name]] = se
  }
  
  # Use cbind to combine SummarizedExperiments
  se_all = Reduce(cbind, se_ls)
  saveRDS(se_all, here("mouse_cortex", "salmon_quants", paste0(pipeline, "_pipeline"), "se_all.rds"))
}


# Get mapping rate
# run_number = "all" # give run_number or "all" for all of them together
# se_ls = list()
# se_ls[["transcripts"]] = readRDS(here("mouse_cortex", "salmon_quants", "transcripts_pipeline", paste0("se_", run_number, ".rds")))
# se_ls[["preandmrna"]] = readRDS(here("mouse_cortex", "salmon_quants", "preandmrna_pipeline", paste0("se_", run_number, ".rds")))
# se_ls[["introncollapse"]] = readRDS(here("mouse_cortex", "salmon_quants", "introncollapse_pipeline", paste0("se_", run_number, ".rds")))
# se_ls[["intronseparate"]] = readRDS(here("mouse_cortex", "salmon_quants", "intronseparate_pipeline", paste0("se_", run_number, ".rds")))
# 
# 
# # Extract mapping rate
# for(se in se_ls){
#   p = c()
#   for(i in -4+6*1:16){
#     p = c(p, metadata(se)[i]$quantInfo$percent_mapped)
#   }
#   print(mean(p))
# }