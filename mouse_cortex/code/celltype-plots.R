# celltype-plots.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Sep 23, 2020
#
# Make all cell type classification-related plots

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(SingleCellExperiment)
})


# Read in SingleCellExperiment objects
run_number = "all" # give run_number or "all" for all of them together
sce_ls = list()
sce_ls[["transcripts"]] = readRDS(here("mouse_cortex", "salmon_quants", "transcripts_pipeline", paste0("sce_", run_number, ".rds")))
sce_ls[["preandmrna"]] = readRDS(here("mouse_cortex", "salmon_quants", "preandmrna_pipeline", paste0("sce_", run_number, ".rds")))
sce_ls[["introncollapse"]] = readRDS(here("mouse_cortex", "salmon_quants", "introncollapse_pipeline", paste0("sce_", run_number, ".rds")))
sce_ls[["intronseparate"]] = readRDS(here("mouse_cortex", "salmon_quants", "intronseparate_pipeline", paste0("sce_", run_number, ".rds")))

##########################
# Plot misclassification #
##########################
tb_misclassify_ls = list()
for(i in seq_along(sce_ls)){
  pipeline = names(sce_ls)[i]
  # Get overlap with singleR labels
  tb = as_tibble(colData(sce_ls[[i]])) %>%
    mutate(cell = rownames(colData(sce_ls[[i]])))
  
  # Label cells where singleR label != ding label
  tb = tb %>%
    mutate(ding_labels_match = case_when(ding_labels == "Astrocyte" ~ "Astrocytes",
                                         ding_labels == "Endothelial" ~ "Endothelial cells",
                                         ding_labels == "Excitatory neuron" ~ "Neurons",
                                         ding_labels == "Inhibitory neuron" ~ "Neurons",
                                         ding_labels == "Oligodendrocyte" ~ "Oligodendrocytes",
                                         T ~ ding_labels))
  
  tb = tb %>%
    mutate(singleR_labels_match = case_when(singleR_labels_pruned == "Astrocytes activated" ~ "Astrocytes",
                                            singleR_labels_pruned == "Fibroblasts activated" ~ "Fibroblasts",
                                            singleR_labels_pruned == "Fibroblasts senescent" ~ "Fibroblasts",
                                            singleR_labels_pruned == "Macrophages activated" ~ "Macrophages",
                                            singleR_labels_pruned == "Microglia activated" ~ "Microglia",
                                            singleR_labels_pruned == "Neurons activated" ~ "Neurons",
                                            T ~ singleR_labels_pruned),
           non_matching_cells = (singleR_labels_match != ding_labels_match) | (is.na(singleR_labels_match) & !is.na(ding_labels_match)))
  
  tb_misclassify_ls[[i]] = tb %>%
    filter(!is.na(ding_labels_match)) %>%
    mutate(pipeline = pipeline)
}

tb_misclassify = bind_rows(tb_misclassify_ls)

p = tb_misclassify %>%
  ggplot(aes(x = pipeline, fill = singleR_labels_match)) +
  geom_bar(position = "stack") +
  facet_wrap(~ ding_labels_match, scales = "free") +
  theme_bw()

ggsave(file = here(paste0("./mouse_cortex/plots/celltypemis_plot.png")), plot = p)

##############################
# Plot ratio of marker genes #
##############################
# Read in SingleR results
pred_celltypes = list()
for(i in seq_along(sce_ls)){
  pipeline = names(sce_ls)[i]
  pred_celltypes[[pipeline]] = readRDS(here(paste0("./mouse_cortex/salmon_quants/"), pipeline, "_pipeline/singler_results.rds"))
}

# Plot ratio of astrocyte:qNSC marker genes for astrocyte cells between pipelines
astrocyte_mg = unique(unlist(all.markers$Astrocytes))
astrocyte_genes = t2g %>% filter(gene_name %in% astrocyte_mg) %>% pull(gene_id)
qnsc_mg = unique(unlist(all.markers$qNSCs))
qnsc_genes = t2g %>% filter(gene_name %in% qnsc_mg) %>% pull(gene_id)

gene_type_name = "Astrocyte qNSC ratio"
gene_sum_tb_ls = list()
for(i in seq_along(sce_ls)){
  pipeline_name = names(sce_ls)[i]
  sce_i = sce_ls[[i]]
  
  lib_size_tb = tibble(cell_barcode = astrocyte_cells)
  lib_size_gene_type_astrocyte = colSums(counts(sce_i)[intersect(rownames(counts(sce_i)), astrocyte_genes), astrocyte_cells, drop = F])
  lib_size_gene_type_qnsc = colSums(counts(sce_i)[intersect(rownames(counts(sce_i)), qnsc_genes), astrocyte_cells, drop = F])
  lib_size_gene_type_tb = tibble(cell_barcode = names(lib_size_gene_type_astrocyte),
                                 !!gene_type_name := lib_size_gene_type_astrocyte/lib_size_gene_type_qnsc)
  lib_size_tb = lib_size_tb %>% 
    left_join(., lib_size_gene_type_tb, by = "cell_barcode")
  gene_sum_tb_ls[[pipeline_name]] = lib_size_tb %>%
    pivot_longer(cols = -cell_barcode) %>%
    mutate(pipeline = pipeline_name)
}
gene_sum_tb = bind_rows(gene_sum_tb_ls)
gene_sum_tb = gene_sum_tb %>%
  pivot_wider(names_from = pipeline, values_from = value)

gene_sum_tb_2 = gene_sum_tb %>%
  pivot_longer(-c("cell_barcode", "name"), names_to = "pipeline", values_to = "ratio")
summ <- gene_sum_tb_2 %>% 
  left_join(., cd, by = c("cell_barcode" = "cell")) %>%
  group_by(pipeline) %>% 
  summarize(median = median(ratio))
gene_sum_tb_2 %>%
  left_join(., cd, by = c("cell_barcode" = "cell")) %>%
  ggplot(aes(x = pipeline, y = ratio)) +
  geom_boxplot() +
  geom_label(data = summ, aes(x = pipeline, y = median, 
                              label = round(median, 2))) +
  labs() +
  theme_bw() 
