# cqn.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Feb 16, 2022
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
  source(here("./mouse_cortex/code/distribution-plots-helpers.R"))
})

downsample = F # boolean toggle for downsampling step
cqn = F        # boolean toggle for cqn step

abb = "ea"
cell_type_1 = "Neurons" # Neurons, Astrocytes, Endothelial cells
cell_type_2 = "Astrocytes"
cell_type_labels = c("Excitatory neuron", "Astrocyte") # Inhibitory neuron, Excitatory neuron, Astrocyte, Endothelial

# Read in SingleCellExperiment objects
run_number = "all" # give run_number or "all" for all of them together
sce_ls = list()
sce_ls[["transcripts"]] = readRDS(here("mouse_cortex", "salmon_quants", "transcripts_pipeline", paste0("sce_", run_number, ".rds")))
sce_ls[["preandmrna"]] = readRDS(here("mouse_cortex", "salmon_quants", "preandmrna_pipeline", paste0("sce_", run_number, ".rds")))
sce_ls[["introncollapse"]] = readRDS(here("mouse_cortex", "salmon_quants", "introncollapse_pipeline", paste0("sce_", run_number, ".rds")))
sce_ls[["intronseparate"]] = readRDS(here("mouse_cortex", "salmon_quants", "intronseparate_pipeline", paste0("sce_", run_number, ".rds")))

pipeline = "intronseparate"
lengths_force_preandmrna = F # if TRUE, always use preandmrna length regardless of pipeline

select_cells = colData(sce_ls[[pipeline]]) %>%
  as.data.frame() %>%
  filter(ding_labels %in% cell_type_labels) %>% 
  filter(cortex == "cortex2") %>%
  row.names()
sce_sub = sce_ls[[pipeline]][, select_cells]

# Estimate p_BC (p(L, C))
pred_celltypes = readRDS(here("mouse_cortex", "salmon_quants", "transcripts_pipeline", "singler_results.rds"))
all.markers <- metadata(pred_celltypes)$de.genes

celltype_markers = all.markers[unique(colData(sce_ls[[pipeline]])$singleR_labels_pruned)]
celltype_markers = sapply(celltype_markers, unlist)
celltype_markers_tb_ls = list()
for(i in seq_along(celltype_markers)){
  if(is.null(unlist(celltype_markers[i]))) next
  celltype_markers_tb_ls[[i]] = tibble(cell_type = names(celltype_markers)[i],
                                       gene_name = unique(unlist(celltype_markers[i])))
}
celltype_markers_tb = bind_rows(celltype_markers_tb_ls)
celltype_markers_length_tb = rowData(sce_ls[[pipeline]]) %>%
  as_tibble() %>%
  dplyr::select(gene_name, start, end, transcript_length) %>%
  mutate(length = end - start,
         gene = rownames(rowData(sce_ls[[pipeline]]))) %>%
  dplyr::select(gene, gene_name, length, transcript_length) %>%
  right_join(., celltype_markers_tb, by = "gene_name")

length_name = ifelse(pipeline == "transcripts" && !lengths_force_preandmrna, "transcript_length", "length")
density_1 = density(log10(celltype_markers_length_tb %>% filter(cell_type == cell_type_1) %>% pull(!!length_name)))
density_2 = density(log10(celltype_markers_length_tb %>% filter(cell_type == cell_type_2) %>% pull(!!length_name))) 
density_est = list(approxfun(density_1), approxfun(density_2))
names(density_est) = cell_type_labels

# Get gene lengths
genes_length_tb = rowData(sce_sub) %>%
  as_tibble() %>%
  dplyr::select(start, end, transcript_length) %>%
  mutate(length = end - start,
         gene = rownames(rowData(sce_sub))) %>%
  dplyr::select(gene, length, transcript_length)
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
colnames(rowData(sce_sub)) = paste(colnames(rowData(sce_sub)), "_temp_fix")
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

dds = DESeqDataSetFromMatrix(countData = counts(seq_data),
                             colData = pData(seq_data),
                             design = ~ ding_labels)

# Downsample according to density p_BC
if(downsample){
  if(pipeline == "transcripts" && !lengths_force_preandmrna){
    lengths_use = genes_length_tb$transcript_length
  } else {
    lengths_use = genes_length_tb$length
  }
  counts_leftover = downsample_by_density(m = counts_sub, density_est = density_est, 
                                          lengths = lengths_use)
} else {
  counts_leftover = counts_sub
}

# Run cQN on leftover counts to get normalizationFactors
if(cqn){
  tic()
  if(pipeline == "transcripts"){
    lengths_use = genes_length_tb$transcript_length
  } else {
    lengths_use = genes_length_tb$length
  }
  cqn_res = cqn(counts = counts_leftover, 
                lengths = lengths_use, # length
                x = genes_length_tb$gc, # GC content
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
resLFC = lfcShrink(dds, coef=resultsNames(dds)[resultsNames(dds) != "Intercept"], type="apeglm")

# MA plot for shrunken log2 fold change
# plotMA(resLFC)

# Add gene length to table
resLFC = resLFC %>%
  as_tibble() %>%
  mutate(gene = rownames(res)) %>%
  left_join(., genes_length_tb, by = "gene")

# Save results
if(downsample & cqn){
  saveRDS(seq_data, here(paste0("./mouse_cortex/output/counts_", abb, "_", pipeline, 
                                "_lengthsforcepreandmrna_", lengths_force_preandmrna, "_cqn_bc.rds"))) # ea, ie
  saveRDS(resLFC, here(paste0("./mouse_cortex/output/de_", abb, "_", pipeline, 
                              "_lengthsforcepreandmrna_", lengths_force_preandmrna, "_lfc_cqn_bc.rds")))
} else if(cqn){
  saveRDS(seq_data, here(paste0("./mouse_cortex/output/counts_", abb, "_", pipeline, 
                                "_lengthsforcepreandmrna_", lengths_force_preandmrna, "_cqn.rds")))
  saveRDS(resLFC, here(paste0("./mouse_cortex/output/de_", abb, "_", pipeline, 
                              "_lengthsforcepreandmrna_", lengths_force_preandmrna, "_lfc_cqn.rds")))
} else {
  saveRDS(seq_data, here(paste0("./mouse_cortex/output/counts_", abb, "_", pipeline, 
                                "_lengthsforcepreandmrna_", lengths_force_preandmrna, ".rds")))
  saveRDS(resLFC, here(paste0("./mouse_cortex/output/de_", abb, "_", pipeline, 
                              "_lengthsforcepreandmrna_", lengths_force_preandmrna, "_lfc.rds")))
}
