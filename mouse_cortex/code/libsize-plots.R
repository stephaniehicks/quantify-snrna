# libsize-plots.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Sep 24, 2020
#
# Make all library size plots


suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(scran)
  library(scater)
  library(BiocSingular)
  library(SingleCellExperiment)
})


# Read in SingleCellExperiment objects
run_number = "all" # give run_number or "all" for all of them together
sce_ls = list()
sce_ls[["transcripts"]] = readRDS(here("mouse_cortex", "salmon_quants", "transcripts_pipeline", paste0("sce_", run_number, ".rds")))
sce_ls[["preandmrna"]] = readRDS(here("mouse_cortex", "salmon_quants", "preandmrna_pipeline", paste0("sce_", run_number, ".rds")))
sce_ls[["introncollapse"]] = readRDS(here("mouse_cortex", "salmon_quants", "introncollapse_pipeline", paste0("sce_", run_number, ".rds")))
sce_ls[["intronseparate"]] = readRDS(here("mouse_cortex", "salmon_quants", "intronseparate_pipeline", paste0("sce_", run_number, ".rds")))


#########################################
# Plot overall library size comparisons #
#########################################
# Get library sizes
lib_sizes_ls = list()
for(i in seq_along(sce_ls)){
  sce = sce_ls[[i]]
  pipeline = names(sce_ls)[i]
  
  lib_sizes_ls[[pipeline]] = colData(sce) %>%
    as_tibble() %>%
    mutate(cell_barcode = rownames(colData(sce)),
           pipeline = pipeline)
}

lib_sizes = bind_rows(lib_sizes_ls)

# Compute median label
summ = lib_sizes %>% 
  group_by(pipeline, cortex) %>% 
  summarize(median = median(lib_size))

# Plot library size comparison between pipelines
lib_sizes %>%
  ggplot(aes(x = pipeline, y = lib_size)) +
  geom_boxplot() + 
  # geom_text(data = summ, aes(x = pipeline, y = median, 
  # label = round(median, 0)), 
  # vjust = -0.5, size = 2.5) +
  facet_wrap(~ cortex) +
  labs(x = "Pipeline",
       y = "Library size") +
  theme_bw() +
  theme(axis.text = element_text(size = 8))

ggsave(file = here(paste0("./mouse_cortex/plots/lib_size_comparison.png")),
       width = 8, height = 5)

########################################################################
# Plot library size comparison grouped by the 4 most common gene types #
########################################################################
# Calculate sums for different gene types
gene_sum_tb_ls = list()
for(i in seq_along(sce_ls)){
  pipeline_name = names(sce_ls)[i]
  sce = sce_ls[[i]]
  
  lib_size_tb = colData(sce) %>%
    as_tibble() %>%
    mutate(cell_barcode = colnames(sce))
  
  for(gene_type_name in unique(rowData(sce)$gene_type)){ 
    genes_in_type = rowData(sce)[rowData(sce)$gene_type == gene_type_name, "gene_id"]
    m = counts(sce)[genes_in_type, ]
    if(length(genes_in_type) > 1){
      lib_size_gene_type = colSums(m)
    } else {
      lib_size_gene_type = m
    }
    lib_size_tb = lib_size_tb %>% 
      mutate(!!gene_type_name := lib_size_gene_type)
  }
  gene_sum_tb_ls[[pipeline_name]] = lib_size_tb %>%
    mutate(pipeline = pipeline_name)
}
gene_sum_tb = bind_rows(gene_sum_tb_ls)

# Plot
for(gene_type in c("protein_coding", "lincRNA", "antisense", "processed_pseudogene")){ # 4 with highest counts across all pipelines
  # Compute median label
  summ = gene_sum_tb %>%
    group_by(pipeline, cortex) %>%
    summarize(median = median(!!sym(gene_type)))
  
  p = gene_sum_tb %>%
    ggplot(aes(x = pipeline, y = !!sym(gene_type))) + 
    geom_boxplot() +
    # geom_label(data = summ, aes(x = pipeline, y = median, 
    #                             label = round(median, 0))) +
    labs(x = "Pipeline",
         y = "Library size") +
    facet_wrap(~ cortex) +
    theme_bw()
  
  ggsave(file = here(paste0("./mouse_cortex/plots/lib_size_", gene_type, "_comparison.png")),
         width = 8, height = 5)
  saveRDS(p, here(paste0("./mouse_cortex/plots/lib_size_", gene_type, "_comparison.rds")))
}