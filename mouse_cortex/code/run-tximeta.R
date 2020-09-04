# run-tximeta.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Sep 4, 2020
#
# Use tximeta to read in quant files (output from alevin) into SummarizedExperiment files

library(here)
suppressPackageStartupMessages({
  library(here)
  library(tximeta)
  library(fishpond)
  library(SummarizedExperiment)
})

pipeline = "transcripts" # "transcripts" or "preandmrna" or "introncollapse/separate"

# load linkedTxome json file if not already in cache (it should be in cache when you create it in quantify-salmon.Rmd)
# if(pipeline == "transcripts"){
#   json_file = here("mouse_cortex", "salmon_files", "gencode.vM25_salmon-index-v1.0.0-transcripts-mouse.json")
#   json_file = "/fastscratch/myscratch/akuo/alsf-filbin/mouse_cortex/salmon_files/gencode.vM25_salmon-index-v1.0.0-transcripts-mouse.json"
#   loadLinkedTxome(json_file)
# } else if(pipeline == "preandmrna"){
#   json_file = here("mouse_cortex", "salmon_files", "gencode.vM25_salmon-index-v1.0.0-preandmrna-mouse.json")
#   loadLinkedTxome(json_file)
# } else if(pipeline == "introncollapse"){
#   json_file = here("mouse_cortex", "salmon_files", "gencode.vM25_salmon-index-v1.0.0-introncollapse-mouse.json")
#   loadLinkedTxome(json_file)
# } else if(pipeline == "intronseparate"){
#   json_file = here("mouse_cortex", "salmon_files", "gencode.vM25_salmon-index-v1.0.0-intronseparate-mouse.json")
#   loadLinkedTxome(json_file)
# }

# Import with tximeta
# Note: alevin import currently only supports a single experiment at a time
run_names = gsub("_quant", "", basename(list.dirs(here("mouse_cortex", "salmon_quants", paste0(pipeline, "_pipeline")), recursive = F)))
run_names = run_names[!grepl("nodecoys", run_names)]

file_paths = here("mouse_cortex", "salmon_quants", paste0(pipeline, "_pipeline"), paste0(run_names, "_quant"), "alevin", "quants_mat.gz") 

for(i in 1:length(run_names)){
  print(file_paths[i])
  se = tximeta(coldata = data.frame(names = run_names[i],
                                    files = file_paths[i],
                                    stringsAsFactors = FALSE),
               type = "alevin")
  
  # Check SummarizedExperiment object
  colData(se)
  assayNames(se)
  rowRanges(se)
  
  # Save SummarizedExperiment object
  saveRDS(se, here("mouse_cortex", "salmon_quants", paste0(pipeline, "_pipeline"), paste0("se_", run_names[i], ".rds")))
}