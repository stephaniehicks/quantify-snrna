# qq-plots-subclusters.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Mar 5, 2023
#
# Make QQ plots for subclusters of inhibitory neurons cortex 1

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
  library(cowplot)
})
source(here("./mouse_cortex/code/distribution-plots-helpers.R"))

##############
# Read data #
##############
tic("entire process")
run_number = "all" 
sce = readRDS(here("mouse_cortex", "salmon_quants", "preandmrna_pipeline", paste0("sce_", run_number, ".rds")))

# Subset to cell type
cell_type = "Inhibitory neuron"
sce_sub = sce[, colData(sce)$ding_labels == cell_type & !is.na(colData(sce)$ding_labels)]
dim(sce_sub)

# Subset to cortex
sce_sub = sce_sub[, colData(sce_sub)$cortex == "cortex1"]
dim(sce_sub)

# Use quickcluster (hierarchical clustering) to get subclusters 
set.seed(1)
tic("quickCluster")
clust = quickCluster(sce_sub, min.size = 0)
toc()

p_ls = list()
for(clust_num in 1:8){
  # Take chosen cluster
  if(clust_num < 8){ # there are 7 subclusters
    counts_sub = counts(sce_sub[, which(clust == clust_num)])
  } else if (clust_num == 8) {
    counts_sub = counts(sce_sub) # run all cells for cell type in last iteration
  }
  print(dim(counts_sub))
  summary(colSums(counts_sub))
  
  ###################
  # Plot QQ plots #
  ###################
  chi_obs_pois = p_chisq_test_2_grouped(counts_sub, distribution = "poisson")
  plot_dt = tibble(chi_obs = chi_obs_pois[[1]]) %>%
    # expression = rowSums(counts_sub_scaled)) %>%
    arrange(chi_obs) %>%
    drop_na() %>%
    mutate(chi_quantile = qchisq(p = (1:length(chi_obs))/length(chi_obs), df = chi_obs_pois[[2]] - 1),
           percentile_color = case_when(chi_obs > quantile(chi_obs, 0.995) ~ "> 99.5",
                                        chi_obs > quantile(chi_obs, 0.95) ~ "> 95",
                                        #chi_obs > quantile(chi_obs, 0.90) ~ "> 90",
                                        T ~ "rest"))
  
  
  p_ls[[clust_num]] = plot_dt %>%
    ggplot(aes(x = chi_quantile,
               y = chi_obs,
               color = percentile_color)) + 
    geom_point(size = 0.5) + 
    geom_abline(slope = 1, intercept = 0) +
    # coord_cartesian(ylim = c(0, 20)) +
    labs(x = "Chi-squared quantiles",
         y = "Goodness of fit statistic") +
    theme_bw()
  toc()
}


# Construct plot
# Apply plot aesthetics to QQ plots
apply_fig1_aesthetics_qq = function(plt){
  plt = plt +
    scale_x_continuous(limits = c(0, 15)) + 
    scale_y_continuous(limits = c(0, 30)) +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 15),
          #plot.title = element_text(size = 20),
          legend.position = "none")
}

# Get legend
legend_qq <- get_legend(
  p_ls[[8]] + 
    scale_color_manual(name = "Percentile",
                       labels = c("<95", ">95", ">99.5"),
                       values = c("#619cff", "#f8766d", "#00ba38")) +
    guides(color = guide_legend(nrow = 1, override.aes = list(size = 2))) +
    theme(legend.position = "bottom", 
          legend.title = element_text(size = 22),
          legend.text = element_text(size = 18))
)

p_ls = lapply(p_ls, FUN = apply_fig1_aesthetics_qq)
p = plot_grid(plot_grid(p_ls[[8]], p_ls[[2]], p_ls[[3]], p_ls[[4]],
              p_ls[[5]], p_ls[[6]], p_ls[[7]], ncol = 4, labels = "auto"),
              legend_qq,
              ncol = 1, rel_heights = c(3, .1))
p
ggsave(here("./mouse_cortex/plots/qq_plots_subclusters.png"), plot = p, width = 16, height = 8)

