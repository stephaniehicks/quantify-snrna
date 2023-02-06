# gene-length-bias-plots-external.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Feb 5, 2023
#
# Make length bias plot for data coming from external sources
# Should be identical to `snuc-analysis-bioc.Rmd`, 
# section gene length bias - overall length bias section

# Packages
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(tictoc)
  library(scran)
  library(scater)
  library(BiocSingular)
  library(SingleCellExperiment)
  library(snakecase)
})

##############
# Read data #
##############
data_source = 2

if(data_source == 1){
  # 10x adult mouse brain nuclei (dataset provided by 10x Genomics)
  library(DropletUtils)
  source_name = "10x_mouse_brain"
  matrix_dir = here("./mouse_cortex/data/filtered_feature_bc_matrix/")
  sce = read10xCounts(matrix_dir)
} else if(data_source == 2){
  # Combined (sNuc-DropSeq, DroNc-seq and 10X Chromium) mouse kidney (dataset from Wu 2019)
  library(data.table)
  source_name = "Dropseq_mouse_kidney"
  matrix_dir = here("./mouse_cortex/data/Healthy.combined.dge.txt.gz")
  counts = fread(matrix_dir)
  gene_names = counts$V1
  counts_sub = counts %>% dplyr::select(starts_with("sNucDropseq"))
  sce = SingleCellExperiment(list(counts = counts_sub))
  
  # Map gene symbols to ENSEMBL IDs (eta: 10 sec)
  library('biomaRt')
  org = "mmusculus_gene_ensembl"
  mart = useDataset(org, useMart("ensembl"))
  g_list = getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
                 filters = "mgi_symbol",
                 values = gene_names,
                 mart = mart)
  g_list = distinct(g_list, mgi_symbol, .keep_all = TRUE)
  gene_match_tb = tibble(mgi_symbol = gene_names)
  gene_match_tb = gene_match_tb %>% left_join(g_list, by = "mgi_symbol")
  rownames(sce) = gene_match_tb$ensembl_gene_id
} else if(data_source == 3){
  # Combined (sNuc-DropSeq, DroNc-seq and 10X Chromium) mouse kidney (dataset from Wu 2019)
  library(data.table)
  source_name = "10x_mouse_kidney"
  matrix_dir = here("./mouse_cortex/data/Healthy.combined.dge.txt.gz")
  counts = fread(matrix_dir)
  gene_names = counts$V1
  counts_sub = counts %>% dplyr::select(starts_with("sNuc.10x"))
  sce = SingleCellExperiment(list(counts = counts_sub))
  
  # Map gene symbols to ENSEMBL IDs (eta: 10 sec)
  library('biomaRt')
  org = "mmusculus_gene_ensembl"
  mart = useDataset(org, useMart("ensembl"))
  g_list = getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
                 filters = "mgi_symbol",
                 values = gene_names,
                 mart = mart)
  g_list = distinct(g_list, mgi_symbol, .keep_all = TRUE)
  gene_match_tb = tibble(mgi_symbol = gene_names)
  gene_match_tb = gene_match_tb %>% left_join(g_list, by = "mgi_symbol")
  rownames(sce) = gene_match_tb$ensembl_gene_id
}

# Get counts
counts_sub = counts(sce)
print(dim(counts_sub))
summary(colSums(counts_sub))

####################
# Gene length bias #
####################
# Get rowsums
genes_rowsums = rowSums(counts_sub)
genes_rowsums_tb = tibble(gene = rownames(sce),
                          count = genes_rowsums)

# Get length
gc_content = readRDS(here("./mouse_cortex/output/gc_content.rds"))
gc_content = gc_content %>% 
  as_tibble(rownames = "gene") %>% 
  dplyr::rename(length_biomart = length)

# Bin by length
genes_rowsums_tb = genes_rowsums_tb %>%
  left_join(gc_content, by = "gene")

genes_rowsums_tb = genes_rowsums_tb %>%
  mutate(bin = ntile(length_biomart, 10))

# Plot lines with ribbon
p = genes_rowsums_tb %>%
  group_by(bin) %>%
  filter(count > 0) %>%
  summarize(boxplot = list(setNames(boxplot.stats(count)$stats,
                                    c('lower_whisker','lower_hinge','median','upper_hinge','upper_whisker')))) %>%
  unnest_wider(boxplot) %>%
  ungroup() %>%
  ggplot(aes(x = bin, y = median)) + 
  geom_ribbon(aes(ymin = lower_hinge, ymax = upper_hinge), alpha = 0.1, linetype = "dashed") +
  geom_point() +
  geom_line() + 
  scale_y_log10() +
  scale_x_continuous(breaks = 1:10) +
  labs(x = "Bin (by transcript length)",
       y = "count",
       title = source_name) +
  theme_bw() +
  theme(legend.position = "top")
saveRDS(p, here(paste0("./mouse_cortex/plots/", source_name, "_transcript_length_bias.rds")))
