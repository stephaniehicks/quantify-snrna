# distribution-plots.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Sep 22, 2020
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

# Plot fraction of zero droplets P(X_i = 0) against mean expression level mu_i
prob_out = plot_prob(counts_sub_scaled)
p = prob_out$plot +
  labs(title = pipeline,
       subtitle = paste(cell_type, cortex, sep = " - "))
ggsave(file = here(paste0("./mouse_cortex/plots/prob0_plot_",
                          paste(pipeline, to_snake_case(cell_type), 
                                to_snake_case(cortex), sep = "_"), 
                          ".png")), plot = p)



# Plot mean-variance
# Empirical mean and variance
mean_emp = rowMeans(counts_sub_scaled)
var_emp = genefilter::rowVars(counts_sub_scaled)

# Binomial
n = round(median(colSums(counts_sub_scaled)))
emp_props = rowSums(counts_sub_scaled)/sum(colSums(counts_sub_scaled))
var_binom = n*emp_props*(1-emp_props)

# Poisson
fit_pois = glmGamPoi::glm_gp(counts_sub_scaled, design = ~ 1, size_factors = FALSE, 
                             overdispersion = FALSE)
# Negative binomial
# Estimate overall size/dispersion parameter
model = lm(var_emp ~ 1*mean_emp + I(mean_emp^2) + 0, tibble(mean_emp, var_emp))
phi = 1/coef(model)["I(mean_emp^2)"]

# Estimate size/dispersion parameter for every gene
# library(glmGamPoi)
# fit_nb =  glmGamPoi::glm_gp(counts_sub_scaled, design = ~ 1, size_factors = FALSE, 
#                          overdispersion = TRUE)

# Note: rowMean(fit_pois$Mu) = rowMeans(fit_nb$Mu) = rowMeans(counts)
mean_var_tb = tibble(mean_emp = mean_emp,
                     var_emp = var_emp,
                     binomial = var_binom,
                     poisson = rowMeans(counts_sub_scaled),
                     nbinomial = mean_emp + mean_emp^2 * 1/phi) %>%
  # var_nb_v2 = mean_emp + mean_emp^2 * fit_nb$overdispersions) %>% 
  tidyr::pivot_longer(cols = -mean_emp, names_to = "model", values_to = "var_value")

# Plot
p = mean_var_tb %>%
  filter(model %in% c("var_emp")) %>%
  ggplot(aes(x = mean_emp, y = var_value)) + 
  geom_point(alpha = 0.3) + 
  geom_line(data = mean_var_tb %>% filter(model %in% c("binomial", "poisson", "nbinomial")),
            aes(x = mean_emp, y = var_value, color = model)) +
  scale_x_log10() + scale_y_log10() +
  labs(title = pipeline,
       subtitle = paste(cell_type, cortex, sep = " - "),
       x = "Log of mean expression",
       y = "Log of variance") +
  theme_bw()
ggsave(file = here(paste0("./mouse_cortex/plots/meanvar_plot_",
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
                nbinomial = nb_bic_1(m),    
                nbinomial_2 = nb_bic_2(m),
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


