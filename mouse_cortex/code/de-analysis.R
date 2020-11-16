# de-analysis.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Nov 16, 2020
#
# Differential expression analysis

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(SingleCellExperiment)
  library(DESeq2)
  library(apeglm)
  library(tictoc)
})


# Read in SingleCellExperiment objects
run_number = "all" # give run_number or "all" for all of them together
sce_ls = list()
sce_ls[["transcripts"]] = readRDS(here("mouse_cortex", "salmon_quants", "transcripts_pipeline", paste0("sce_", run_number, ".rds")))
sce_ls[["preandmrna"]] = readRDS(here("mouse_cortex", "salmon_quants", "preandmrna_pipeline", paste0("sce_", run_number, ".rds")))
sce_ls[["introncollapse"]] = readRDS(here("mouse_cortex", "salmon_quants", "introncollapse_pipeline", paste0("sce_", run_number, ".rds")))
sce_ls[["intronseparate"]] = readRDS(here("mouse_cortex", "salmon_quants", "intronseparate_pipeline", paste0("sce_", run_number, ".rds")))



pipeline = "preandmrna"
sce_sub = sce_ls[[pipeline]][, colData(sce_ls[[pipeline]])$ding_labels %in% c("Excitatory neuron", "Inhibitory neuron")]

# Store counts and intermediate quantities (dds is a SummarizedExperiment)
dds = DESeqDataSetFromMatrix(countData = ceiling(counts(sce_sub)[, ]), # take integers 
                             colData = colData(sce_sub),
                             design = ~ cortex + ding_labels)

# Compute log2fold change and p-values
tic()
dds = DESeq(dds)
toc()

res = results(dds)
res

# Shrinkage of LFC
resultsNames(dds)
resLFC = lfcShrink(dds, coef="ding_labels_Inhibitory.neuron_vs_Excitatory.neuron", type="apeglm")

# MA plot for shrunken log2 fold change
plotMA(resLFC)

# Add gene length to table
genes_length_tb = rowData(sce_ls[[pipeline]]) %>%
  as_tibble() %>%
  select(start, end) %>%
  mutate(length = end - start,
         gene = rownames(rowData(sce_ls[[pipeline]]))) %>%
  select(gene, length)

res = res %>%
  as_tibble() %>%
  mutate(gene = rownames(res)) %>%
  left_join(., genes_length_tb, by = "gene")

saveRDS(res, here("./mouse_cortex/output/res.rds"))