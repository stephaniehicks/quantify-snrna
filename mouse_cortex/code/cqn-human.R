# cqn-human.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Feb 8, 2022
#
# Copy of cqn.R to run on human samples

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
  source(here("./mouse_cortex/code/distribution-plots-helpers.R"))
})

cqn = T        # boolean toggle for cqn step

# Load SCE objects
# load(here("mouse_cortex", "data", "SCE_AMY-n5_tran-etal.rda"))
load(here("mouse_cortex", "data", "SCE_DLPFC-n3_tran-etal.rda"))
# load(here("mouse_cortex", "data", "SCE_HPC-n3_tran-etal.rda"))
# load(here("mouse_cortex", "data", "SCE_NAc-n8_tran-etal.rda"))
# load(here("mouse_cortex", "data", "SCE_sACC-n5_tran-etal.rda"))

# Comparison
abb = "donor1-2_Oligo"
cell_type = c("Oligo")
samples = c("donor1", "donor2")
sce = sce.dlpfc.tran

select_cells = colData(sce) %>%
  as.data.frame() %>%
  filter(cellType %in% cell_type) %>%
  filter(donor %in% samples) %>%
  row.names()
sce_sub = sce[, select_cells]

# Get gene lengths
gc_content_human = readRDS(here("./mouse_cortex/output/gc_content_human.rds"))
gc_content_human = gc_content_human %>% 
  as_tibble(rownames = "gene") %>% 
  dplyr::rename(length_biomart = length) %>%
  filter(gene %in% rowData(sce)$gene_id)

rowData(sce) = rowData(sce) %>%
  as.data.frame() %>%
  left_join(., gc_content_human, by = c("gene_id" = "gene"))

# Remove genes with missing GC content or length
keep_genes = which(!is.na(gc_content_human$gc) & !is.na(gc_content_human$length_biomart))
sce_sub = sce_sub[keep_genes, ]
counts_sub = as.matrix(round(counts(sce_sub)))

gc_content_human = gc_content_human[keep_genes, ]
gc_content_human = as.data.frame(gc_content_human)
rownames(gc_content_human) = gc_content_human$gene
size_factors = colData(sce_sub)$sizeFactor
# nonzero_sums = which(rowSums(counts_sub) != 0)
# counts_sub = counts_sub[nonzero_sums, ]

# Aggregate over sample pseudo-groups (2 per donor)
counts_sub = as.matrix(counts(sce_sub))
colnames(counts_sub) = colData(sce_sub)$donor
counts_sub = t(rowsum(t(counts_sub), paste(colnames(counts_sub), sample(1:2, ncol(counts_sub), replace = TRUE))))
counts_sub = round(counts_sub)
rownames(counts_sub) = rownames(gc_content_human)

pdata = as.data.frame(tibble(group = colnames(counts_sub),
                             sample_labels = sapply(colnames(counts_sub), function(x) gsub(".{2}$", "", x))))

rownames(pdata) = colnames(counts_sub)
seq_data = newSeqExpressionSet(counts = counts_sub,
                               featureData = gc_content_human,
                               phenoData = pdata)

dds = DESeqDataSetFromMatrix(countData = counts(seq_data),
                             colData = pData(seq_data),
                             design = ~ sample_labels)

# Run cQN to get normalizationFactors
if(cqn){
  counts_leftover = counts_sub
  tic()
  cqn_res = cqn(counts = counts_leftover, 
                lengths = gc_content_human$length_biomart, # length
                x = gc_content_human$gc, # GC content
                # subindex = which(rowMeans(counts_sub) > 15), # Default is rowMeans > 50
                # sizeFactors = size_factors,
                verbose = FALSE)
  toc()
  
  # Get offset from cqn to use as normalizationfactors in DESeq function
  cqnOffset <- cqn_res$glm.offset
  cqnNormFactors <- exp(cqnOffset)
  normFactors <- cqnNormFactors / exp(rowMeans(log(cqnNormFactors)))
  normalizationFactors(dds) <- normFactors
}

# Compute log2fold change and p-values
tic("DESeq")
dds = DESeq(dds)
toc()

res = results(dds)
res

# Shrinkage of LFC when count values are too low
resultsNames(dds)
resLFC = lfcShrink(dds, coef=resultsNames(dds)[resultsNames(dds) != "Intercept"], 
                   type="apeglm") 

# MA plot for shrunken log2 fold change
plotMA(resLFC)

# Add gene length and GC content to table
genes_length_tb = rowData(sce) %>%
  as_tibble() %>%
  dplyr::select(gene_id, length_biomart, gc)

res = res %>%
  as_tibble() %>%
  mutate(gene = rownames(res)) %>%
  left_join(., genes_length_tb, by = c("gene" = "gene_id"))

resLFC = resLFC %>%
  as_tibble() %>%
  mutate(gene = rownames(resLFC)) %>%
  left_join(., genes_length_tb, by = c("gene" = "gene_id"))

# Save results
if(cqn){
  saveRDS(seq_data, here(paste0("./mouse_cortex/output/counts_", abb, "_", "human", "_cqn.rds")))
  saveRDS(res, here(paste0("./mouse_cortex/output/de_", abb, "_", "human", "_cqn.rds")))
  saveRDS(resLFC, here(paste0("./mouse_cortex/output/de_", abb, "_", "human", "_lfc_cqn.rds")))
} else {
  saveRDS(seq_data, here(paste0("./mouse_cortex/output/counts_", abb, "_", "human", ".rds")))
  saveRDS(res, here(paste0("./mouse_cortex/output/de_", abb, "_", "human", ".rds")))
  saveRDS(resLFC, here(paste0("./mouse_cortex/output/de_", abb, "_", "human", "_lfc.rds")))
}
