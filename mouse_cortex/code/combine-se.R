# combine-se.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Oct 8, 2020
#
# Combine SummarizedExperiment objects into one 

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(SummarizedExperiment)
})

# pipeline = "transcripts" # "transcripts" or "preandmrna" or "introncollapse/separate"
for(pipeline in c("transcripts", "preandmrna", "introncollapse", "intronseparate")){
  # Get run names
  run_names = gsub("_quant", "", basename(list.dirs(here("mouse_cortex", "salmon_quants", paste0(pipeline, "_pipeline")), recursive = F)))
  run_names = run_names[!grepl("nodecoys", run_names)]
  
  # Get cortex/flow cell metadata for each run
  sra_meta = read_csv(here("mouse_cortex", "files", "SRA_geo_metadata.csv"))
  sra_meta = sra_meta %>%
    filter(Run %in% run_names) %>%
    select(Run, title) %>%
    separate(title, c("cortex", "flow_cell", "seq", "lane"))
  
  # Get t2g metadata for each gene
  tx2gene = readRDS(here("mouse_cortex", "salmon_files", "gencode.vM25.transcripts.tx2gene.mouse.rds"))
  
  # Read in individual rds files
  se_ls = list()
  for(run_name in run_names[1:16]){
    se = readRDS(here("mouse_cortex", "salmon_quants", paste0(pipeline, "_pipeline"), paste0("se_", run_name, ".rds")))
    # Add colData
    column_names = colnames(se)
    colData(se) = DataFrame(sra_meta %>% 
                              dplyr::filter(Run == run_name) %>%
                              dplyr::slice(rep(1, ncol(se))))
    colnames(se) = column_names
    
    # Add rowData
    match_rows = match(rownames(se), tx2gene$gene_id)
    t2g = tx2gene[match_rows, ]
    rowData(se) = t2g
    
    se_ls[[run_name]] = se
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