# distribution-plots-external.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Jan 29, 2023
#
# Make all distribution-related plots for data coming from external sources
# Should be identical to distribution-plots.R other than the data source

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

##############
# Read data #
##############
tic()
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
  rownames(sce) = gene_names
}

# Subset (may need to subset by cell type if available?)
counts_sub = counts(sce)
print(dim(counts_sub))
summary(colSums(counts_sub))

# Downsample cell counts
counts_sub_scaled = Down_Sample_Matrix(ceiling(counts_sub[1:dim(sce)[1], 1:100]))
counts_sub_scaled = counts_sub_scaled[rowSums(counts_sub_scaled) != 0, ]
summary(colSums(counts_sub_scaled))

################################################################################
# Plot fraction of zero droplets P(X_i = 0) against mean expression level mu_i #
################################################################################
prob_out = plot_prob(counts_sub_scaled)
p = prob_out$plot +
  labs(title = source_name)
ggsave(file = here(paste0("./mouse_cortex/plots/prob0_plot_", source_name, 
                          ".png")), plot = p)
saveRDS(p, file = here(paste0("./mouse_cortex/plots/prob0_plot_", source_name, 
                              ".rds")))

######################
# Plot mean-variance #
######################
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
  labs(title = source_name,
       x = "Log of mean expression",
       y = "Log of variance") +
  theme_bw()
ggsave(file = here(paste0("./mouse_cortex/plots/meanvar_plot_",
                          source_name, 
                          ".png")), plot = p)
saveRDS(p, file = here(paste0("./mouse_cortex/plots/meanvar_plot_",
                              source_name, 
                              ".rds")))

###################
# Plot BIC values #
###################
m = counts_sub_scaled
summary(colSums(m))
summary(rowSums(m))
dim(m)

bic_tb = tibble(binomial = mult_bic(m),     
                # dmn = dmn_bic(m),     # Unknown time
                poisson = poi_bic(m),       
                nbinomial = nb_bic_1(m),    
                nbinomial_2 = nb_bic_2(m),
                cell_type_name = "placeholder")

p = bic_tb %>% 
  pivot_longer(cols = -cell_type_name) %>%
  ggplot(aes(x = reorder(name, value), y = value)) +
  geom_point(size = 3) +
  labs(title = source_name,
       x = "Distribution",
       y = "BIC") +
  theme_bw()

ggsave(file = here(paste0("./mouse_cortex/plots/bic_plot_",
                          source_name, 
                          ".png")), plot = p)
saveRDS(p, file = here(paste0("./mouse_cortex/plots/bic_plot_",
                              source_name, 
                              ".rds")))
saveRDS(bic_tb, here(paste0("./mouse_cortex/plots/bic_data_",
                            source_name, 
                            ".rds")))

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


p = plot_dt %>%
  ggplot(aes(x = chi_quantile,
             y = chi_obs,
             color = percentile_color)) + 
  geom_point(size = 0.5) + 
  geom_abline(slope = 1, intercept = 0) +
  # coord_cartesian(ylim = c(0, 20)) +
  labs(x = "Chi-squared quantiles",
       y = "Goodness of fit statistic") +
  # guides(color = FALSE) +
  theme_bw()
print(p)
toc()
saveRDS(p, here(paste0("./mouse_cortex/plots/qq_plot_poisson_", source_name, ".rds")))
