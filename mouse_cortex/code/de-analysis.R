# de-analysis.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Sep 13, 2021
#
# Differential expression analysis

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(SingleCellExperiment)
  library(DESeq2)
  library(EDASeq)
  library(apeglm)
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
select_cells = colData(sce_ls[[pipeline]]) %>%
  as.data.frame() %>%
  filter(ding_labels %in% c("Excitatory neuron", "Astrocyte")) %>%
  filter(cortex == "cortex1") %>%
  row.names()
sce_sub = sce_ls[[pipeline]][, select_cells]

# Get gene lengths
genes_length_tb = rowData(sce_sub) %>%
  as_tibble() %>%
  select(start, end) %>%
  mutate(length = end - start,
         gene = rownames(rowData(sce_ls[[pipeline]]))) %>%
  select(gene, length)
genes_length_tb = as.data.frame(genes_length_tb)
rownames(genes_length_tb) = genes_length_tb$gene

# Normalize with EDASeq
seq_data = newSeqExpressionSet(counts = as.matrix(counts(sce_sub)),
                               featureData = as.data.frame(genes_length_tb),
                               phenoData = as.data.frame(colData(sce_sub)))
if(normalize){
  tic("normalize with EDASeq")
  seq_data_within = withinLaneNormalization(seq_data, "length", which = "full", offset = TRUE)
  seq_data_norm = betweenLaneNormalization(seq_data_within, which = "full")
  toc()
} else {
  seq_data_norm = seq_data
}

# Store counts and intermediate quantities (dds is a SummarizedExperiment)
# dds = DESeqDataSetFromMatrix(countData = ceiling(counts(sce_sub)[, ]), # take integers 
#                              colData = colData(sce_sub),
#                              design = ~ cortex + ding_labels)

dds = DESeqDataSetFromMatrix(countData = ceiling(counts(seq_data_norm)[, ]), # take integers 
                             colData = pData(seq_data_norm),
                             design = ~ ding_labels)

if(normalize){
  normFactors <- exp(-1 * offst(seq_data_norm))
  normFactors <- normFactors / exp(rowMeans(log(normFactors)))
  normalizationFactors(dds) <- normFactors
}

# Compute log2fold change and p-values
tic("DEseq")
dds = DESeq(dds)
toc()

res = results(dds)
res

# Shrinkage of LFC when count values are too low
resultsNames(dds)
resLFC = lfcShrink(dds, coef="ding_labels_Excitatory.neuron_vs_Astrocyte", type="apeglm")

# MA plot for shrunken log2 fold change
plotMA(resLFC)

# Add gene length to table
genes_length_tb = rowData(sce_ls[[pipeline]]) %>%
  as_tibble() %>%
  select(start, end) %>%
  mutate(length = end - start,
         gene = rownames(rowData(sce_ls[[pipeline]]))) %>%
  select(gene, length)

resLFC = resLFC %>%
  as_tibble() %>%
  mutate(gene = rownames(res)) %>%
  left_join(., genes_length_tb, by = "gene")

if(normalize){
  saveRDS(resLFC, here(paste0("./mouse_cortex/output/de_ea_", pipeline, "_lfc_norm.rds")))
} else {
  saveRDS(resLFC, here(paste0("./mouse_cortex/output/de_ea_", pipeline, "_lfc.rds")))
}