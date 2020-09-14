# distribution-plots.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Sep 14, 2020
#
# Make all distribution-related plots 

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
source(here("./mouse_cortex/code/distribution-plots-helpers.R"))

# Parameters
task_id = as.integer(Sys.getenv("SGE_TASK_ID"))
pipeline_ls = c("transcripts", "preandmrna", "intronseparate", "introncollapse")
cell_type_ls = c("Excitatory neuron", "Inhibitory neuron", "Astrocyte", "Oligodendrocyte", "OPC", "Endothelial", "Microglia")
cortex_ls = c("Cortex1", "Cortex2")

pipeline = pipeline_ls[ceiling(task_id/(length(cortex_ls)*length(cell_type_ls)))]
cell_type = cell_type_ls[ceiling(task_id/length(cortex_ls)) %% length(cell_type_ls) + 1]
cortex = cortex_ls[task_id %% length(cortex_ls) + 1]

# pipeline = "transcripts" # Options are "transcripts", "preandmrna", "intronseparate", and "introncollapse"
# cell_type = "Inhibitory neuron" # Options are "Excitatory neuron", "Inhibitory neuron", "Astrocyte", "Oligodendrocyte", "OPC", "Endothelial", and "Microglia"
# cortex = "Cortex1" # Options are "Cortex1" or "Cortex2"

# Read in data
sce = readRDS(here("mouse_cortex", "salmon_quants", 
                   paste0(pipeline, "_pipeline"), 
                   "sce_all.rds"))

# Add Ding cell type labels
meta_ding = read_tsv(here("mouse_cortex", "files", "meta_combined.txt"))
meta_ding_10x = meta_ding %>%
  filter(Method == "10x Chromium") %>%
  mutate(cell_barcode = gsub("Cortex.*10xChromium", "", NAME))

match_cols = match(colnames(sce), meta_ding_10x$cell_barcode)
colData(sce)$ding_labels = meta_ding_10x$CellType[match_cols]

# Subset to cell type and cortex
sce_sub = sce[, colData(sce)$ding_labels == cell_type & 
                !is.na(colData(sce)$ding_labels) & # Subset to cell type
                colData(sce)$cortex == cortex]     # Subset to cortex
counts_sub = counts(sce_sub)
print(dim(counts_sub))
summary(colSums(counts_sub))

# Downsample cell counts
counts_sub_scaled = Down_Sample_Matrix(ceiling(counts_sub))
counts_sub_scaled = counts_sub_scaled[rowSums(counts_sub_scaled) != 0, ]
summary(colSums(counts_sub_scaled))

# Plot P(X_i = 0) against average expression level mu_i
prob_out = plot_prob(counts_sub_scaled)
p = prob_out$plot +
  labs(title = pipeline,
       subtitle = paste(cell_type, cortex, sep = " - "))
ggsave(file = here(paste0("./mouse_cortex/plots/prob0_plot_",
                          paste(pipeline, to_snake_case(cell_type), 
                                to_snake_case(cortex), sep = "_"), 
                          ".png")), plot = p)

# Plot BIC values
m = counts_sub_scaled
summary(colSums(m))
summary(rowSums(m))
dim(m)

bic_tb = tibble(multinomial = mult_bic(m),     
                # dmn = dmn_bic(m),     # Unknown time
                poisson = poi_bic(m),       
                negative_binomial_1 = nb_bic_1(m),    
                negative_binomial_2 = nb_bic_2(m),
                cell_type_name = cell_type)

p = bic_tb %>% 
  pivot_longer(cols = -cell_type_name) %>%
  ggplot(aes(x = reorder(name, value), y = value)) +
  geom_point(size = 3) +
  labs(title = pipeline,
       subtitle = paste(cell_type, cortex, sep = " - "),
       x = "Distribution",
       y = "BIC") +
  theme_bw()

ggsave(file = here(paste0("./mouse_cortex/plots/bic_plot_",
                          paste(pipeline, to_snake_case(cell_type), 
                                to_snake_case(cortex), sep = "_"), 
                          ".png")), plot = p)


