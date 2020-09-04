# run-tximeta.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Sep 4, 2020
#
# Use tximeta to read in quant files (output from alevin) into SummarizedExperiment files

suppressPackageStartupMessages({
  library(here)
  library(tximeta)
})

pipeline = "transcripts" # "transcripts" or "preandmrna" or "introncollapse/separate"

# Create linkedTranscriptome
if(pipeline == "transcripts"){
  index_dir = here("mouse_cortex", "salmon_files", "gencode.vM25_salmon-index-v1.0.0-transcripts-mouse")
  fasta_path = here("mouse_cortex", "salmon_files", "gencode.vM25.transcripts.mouse.fa.gz")
  gtf_path = here("mouse_cortex", "salmon_files", "gencode.vM25.annotation.gtf.gz") 
  json_file = here("mouse_cortex", "salmon_files", paste0(basename(index_dir), ".json"))
  makeLinkedTxome(indexDir=index_dir,
                  source="other", organism="mouse",
                  release="other", genome="GRCm38",
                  fasta=fasta_path, gtf=gtf_path,
                  jsonFile=json_file) # this command will add the index to the cache automatically
  
} else if(pipeline == "preandmrna"){
  index_dir = here("mouse_cortex", "salmon_files", "gencode.vM25_salmon-index-v1.0.0-preandmrna-mouse")
  fasta_path = here("mouse_cortex", "salmon_files", "gencode.vM25.preandmrna.mouse.fa.gz")
  gtf_path = here("mouse_cortex", "salmon_files", "gencode.vM25.annotation.gtf.gz") 
  json_file = here("mouse_cortex", "salmon_files", paste0(basename(index_dir), ".json"))
  makeLinkedTxome(indexDir=index_dir,
                  source="other", organism="mouse",
                  release="other", genome="GRCm38",
                  fasta=fasta_path, gtf=gtf_path,
                  jsonFile=json_file) # this command will add the index to the cache automatically
} else if(pipeline == "introncollapse"){
  index_dir = here("mouse_cortex", "salmon_files", "gencode.vM25_salmon-index-v1.0.0-introncollapse-mouse")
  fasta_path = here("mouse_cortex", "salmon_files", "gencode.vM25.introncollapse.mouse.fa.gz")
  gtf_path = here("mouse_cortex", "salmon_files", "gencode.vM25.annotation.gtf.gz") 
  json_file = here("mouse_cortex", "salmon_files", paste0(basename(index_dir), ".json"))
  makeLinkedTxome(indexDir=index_dir,
                  source="other", organism="mouse",
                  release="other", genome="GRCm38",
                  fasta=fasta_path, gtf=gtf_path,
                  jsonFile=json_file) # this command will add the index to the cache automatically
} else if(pipeline == "intronseparate"){
  index_dir = here("mouse_cortex", "salmon_files", "gencode.vM25_salmon-index-v1.0.0-intronseparate-mouse")
  fasta_path = here("mouse_cortex", "salmon_files", "gencode.vM25.intronseparate.mouse.fa.gz")
  gtf_path = here("mouse_cortex", "salmon_files", "gencode.vM25.annotation.gtf.gz") 
  json_file = here("mouse_cortex", "salmon_files", paste0(basename(index_dir), ".json"))
  makeLinkedTxome(indexDir=index_dir,
                  source="other", organism="mouse",
                  release="other", genome="GRCm38",
                  fasta=fasta_path, gtf=gtf_path,
                  jsonFile=json_file) # this command will add the index to the cache automatically
}

# Import with tximeta
# Note: alevin import currently only supports a single experiment at a time
run_names = gsub("_quant", "", basename(list.dirs(here("mouse_cortex", "salmon_quants", paste0(pipeline, "_pipeline")), recursive = F)))
run_names = run_names[!grepl("nodecoys", run_names)]
file_paths = here("mouse_cortex", "salmon_quants", paste0(pipeline, "_pipeline"), paste0(run_names, "_quant"), "alevin", "quants_mat.gz") 

for(i in 1:length(run_names)){
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