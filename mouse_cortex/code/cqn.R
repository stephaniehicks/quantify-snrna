# cqn.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Oct 20, 2021
#
# Run cQN to normalize and then run DE analysis.

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(SingleCellExperiment)
  library(DESeq2)
  library(EDASeq)
  library(apeglm)
  library(tictoc)
  library(cqn)
  library(scales)
})

# Read in SingleCellExperiment objects
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
  filter(cortex == "cortex2") %>%
  row.names()
sce_sub = sce_ls[[pipeline]][, select_cells]

# Get gene lengths
genes_length_tb = rowData(sce_sub) %>%
  as_tibble() %>%
  select(start, end) %>%
  mutate(length = end - start,
         gene = rownames(rowData(sce_sub))) %>%
  select(gene, length)
genes_length_tb = as.data.frame(genes_length_tb)

# Get GC content
gc_content = readRDS(here("./mouse_cortex/output/gc_content.rds"))
gc_content = gc_content %>% 
  as_tibble(rownames = "gene") %>% 
  dplyr::rename(length_biomart = length)
genes_length_tb = genes_length_tb %>%
  mutate(gene_noversion = gsub("\\..*", "", gene)) %>%
  left_join(., gc_content, by = c("gene_noversion" = "gene"))
rownames(genes_length_tb) = genes_length_tb$gene

# Remove genes with missing GC content
keep_genes = which(!is.na(genes_length_tb$gc))
sce_sub = sce_sub[keep_genes, ]
counts_sub = as.matrix(round(counts(sce_sub)))

genes_length_tb = genes_length_tb[keep_genes, ]
size_factors = colData(sce_sub)$sizeFactor
# nonzero_sums = which(rowSums(counts_sub) != 0)
# counts_sub = counts_sub[nonzero_sums, ]

# Aggregate over cell type pseudo-groups (5 per cell type)
counts_sub = as.matrix(counts(sce_sub))
colnames(counts_sub) = colData(sce_sub)$ding_labels
counts_sub = t(rowsum(t(counts_sub), paste(colnames(counts_sub), sample(1:5, ncol(counts_sub), replace = TRUE))))
counts_sub = round(counts_sub)

pdata = as.data.frame(tibble(group = colnames(counts_sub),
                             ding_labels = sapply(colnames(counts_sub), function(x) gsub(".{2}$", "", x))))

rownames(pdata) = colnames(counts_sub)
seq_data = newSeqExpressionSet(counts = counts_sub,
                               featureData = as.data.frame(genes_length_tb),
                               phenoData = pdata)

dds = DESeqDataSetFromMatrix(countData = counts(seq_data), # take integers 
                             colData = pData(seq_data),
                             design = ~ ding_labels)

# Run cQN to get normalizationFactors
tic()
cqn_res = cqn(counts = counts_sub, 
              lengths = genes_length_tb$length, # length
              x = genes_length_tb$gc, # GC content
              # subindex = which(rowMeans(counts_sub) > 15), # Default is rowMeans > 50
              # sizeFactors = size_factors,
              verbose = FALSE)
toc()

# Get offset from cqn to use as normalizationfactors in DESeq function
cqnOffset <- cqn_res$glm.offset
cqnNormFactors <- exp(cqnOffset)
normFactors <- cqnNormFactors / exp(rowMeans(log(cqnNormFactors)))

# Compute log2fold change and p-values
tic("DESeq")
normalizationFactors(dds) <- normFactors
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

# Save results
saveRDS(seq_data, here(paste0("./mouse_cortex/output/counts_ea_", pipeline, "_cqn.rds")))
saveRDS(resLFC, here(paste0("./mouse_cortex/output/de_ea_", pipeline, "_lfc_cqn.rds")))
