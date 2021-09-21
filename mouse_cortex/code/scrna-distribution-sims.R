# scrna-distribution-sims.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Apr 23, 2021
#
# Run chi-squared tests (by sample) for simulations/scRNA data in parallel

library(tidyverse)
library(here)
source(here("./mouse_cortex/code/distribution-plots-helpers.R"))

# Get task ID
task_id = as.integer(Sys.getenv("SGE_TASK_ID")) 
i = task_id
print(i)

# Run chi-squared tests
counts_ls = readRDS(here("./scrna/counts_more_ls.rds"))
counts = counts_ls[[i]]

# Poisson
chi_obs_pois = p_chisq_test_2(counts, distribution = "poisson")
plot_dt_pois = tibble(chi_obs = chi_obs_pois) %>%
  # expression = rowSums(counts),
  # var = apply(counts, 1, var),
  # prop_0 = apply(counts, 1, function(x) sum(x == 0))) %>%
  mutate(id = 1:n()) %>%
  arrange(chi_obs) %>%
  drop_na() %>%
  mutate(chi_quantile = qchisq(p = (1:length(chi_obs))/length(chi_obs), df = ncol(counts) - 1),
         percentile_color = case_when(chi_obs > quantile(chi_obs, 0.995, na.rm = TRUE) ~ "> 99.5",
                                      chi_obs > quantile(chi_obs, 0.95, na.rm = TRUE) ~ "> 95",
                                      chi_obs > quantile(chi_obs, 0.90, na.rm = TRUE) ~ "> 90",
                                      T ~ "rest"),
         i = i)

saveRDS(plot_dt_pois, here(paste0("./scrna/output/plot_dt_pois_ls_", i, ".rds")))

# Poisson (grouped)
chi_obs_pois = p_chisq_test_2_grouped(counts, distribution = "poisson")
plot_dt_pois = tibble(chi_obs = chi_obs_pois[[1]]) %>%
                         # expression = rowSums(counts),
                         # var = apply(counts, 1, var),
                         # prop_0 = apply(counts, 1, function(x) sum(x == 0))) %>%
  mutate(id = 1:n()) %>%
  arrange(chi_obs) %>%
  drop_na() %>%
  mutate(chi_quantile = qchisq(p = (1:length(chi_obs))/length(chi_obs), df = chi_obs_pois[[2]] - 1),
         percentile_color = case_when(chi_obs > quantile(chi_obs, 0.995, na.rm = TRUE) ~ "> 99.5",
                                      chi_obs > quantile(chi_obs, 0.95, na.rm = TRUE) ~ "> 95",
                                      chi_obs > quantile(chi_obs, 0.90, na.rm = TRUE) ~ "> 90",
                                      T ~ "rest"),
         i = i)

saveRDS(plot_dt_pois, here(paste0("./scrna/output/plot_dt_pois_grouped_ls_", i, ".rds")))

# Negative binomial (single overdispersion parameter)
chi_obs_nb = p_chisq_test_2(counts, distribution = "nb 1")
plot_dt_nb_1 = tibble(chi_obs = chi_obs_nb) %>%
                         # expression = rowSums(counts),
                         # var = apply(counts, 1, var),
                         # prop_0 = apply(counts, 1, function(x) sum(x == 0))) %>%
  mutate(id = 1:n()) %>%
  arrange(chi_obs) %>%
  drop_na() %>%
  mutate(chi_quantile = qchisq(p = (1:length(chi_obs))/length(chi_obs), df = ncol(counts) - 2),
         percentile_color = case_when(chi_obs > quantile(chi_obs, 0.995, na.rm = TRUE) ~ "> 99.5",
                                      chi_obs > quantile(chi_obs, 0.95, na.rm = TRUE) ~ "> 95",
                                      chi_obs > quantile(chi_obs, 0.90, na.rm = TRUE) ~ "> 90",
                                      T ~ "rest"),
         i = i)

saveRDS(plot_dt_nb_1, here(paste0("./scrna/output/plot_dt_nb_1_ls_", i, ".rds")))

# Negative binomial (single overdispersion parameter - grouped)
chi_obs_nb = p_chisq_test_2_grouped(counts, distribution = "nb 1")
plot_dt_nb = tibble(chi_obs = chi_obs_nb[[1]]) %>%
  # expression = rowSums(counts),
  # var = apply(counts, 1, var),
  # prop_0 = apply(counts, 1, function(x) sum(x == 0))) %>%
  mutate(id = 1:n()) %>%
  arrange(chi_obs) %>%
  drop_na() %>%
  mutate(chi_quantile = qchisq(p = (1:length(chi_obs))/length(chi_obs), df = chi_obs_nb[[2]] - 2),
         percentile_color = case_when(chi_obs > quantile(chi_obs, 0.995, na.rm = TRUE) ~ "> 99.5",
                                      chi_obs > quantile(chi_obs, 0.95, na.rm = TRUE) ~ "> 95",
                                      chi_obs > quantile(chi_obs, 0.90, na.rm = TRUE) ~ "> 90",
                                      T ~ "rest"),
         i = i)

saveRDS(plot_dt_nb, here(paste0("./scrna/output/plot_dt_nb_1_grouped_ls_", i, ".rds")))

# Negative binomial (gene-specific overdispersion parameter)
chi_obs_nb = p_chisq_test_2(counts, distribution = "nb 2")
plot_dt_nb_2 = tibble(chi_obs = chi_obs_nb[[1]],
                         f_phi = chi_obs_nb[[2]]) %>%
                         # expression = rowSums(counts),
                         # var = apply(counts, 1, var),
                         # prop_0 = apply(counts, 1, function(x) sum(x == 0))) %>%
  mutate(id = 1:n()) %>%
  arrange(chi_obs) %>%
  drop_na() %>%
  mutate(chi_quantile = qchisq(p = (1:length(chi_obs))/length(chi_obs), df = ncol(counts) - 2),
         percentile_color = case_when(chi_obs > quantile(chi_obs, 0.995, na.rm = TRUE) ~ "> 99.5",
                                      chi_obs > quantile(chi_obs, 0.95, na.rm = TRUE) ~ "> 95",
                                      chi_obs > quantile(chi_obs, 0.90, na.rm = TRUE) ~ "> 90",
                                      T ~ "rest"),
         i = i)

saveRDS(plot_dt_nb_2, here(paste0("./scrna/output/plot_dt_nb_2_ls_", i, ".rds")))

# Negative binomial (gene-specific overdispersion parameter - grouped)
chi_obs_nb = p_chisq_test_2_grouped(counts, distribution = "nb 2")
plot_dt_nb = tibble(chi_obs = chi_obs_nb[[1]]) %>%
  # expression = rowSums(counts),
  # var = apply(counts, 1, var),
  # prop_0 = apply(counts, 1, function(x) sum(x == 0))) %>%
  mutate(id = 1:n()) %>%
  arrange(chi_obs) %>%
  drop_na() %>%
  mutate(chi_quantile = qchisq(p = (1:length(chi_obs))/length(chi_obs), df = chi_obs_nb[[3]] - 2),
         percentile_color = case_when(chi_obs > quantile(chi_obs, 0.995, na.rm = TRUE) ~ "> 99.5",
                                      chi_obs > quantile(chi_obs, 0.95, na.rm = TRUE) ~ "> 95",
                                      chi_obs > quantile(chi_obs, 0.90, na.rm = TRUE) ~ "> 90",
                                      T ~ "rest"),
         i = i)

saveRDS(plot_dt_nb, here(paste0("./scrna/output/plot_dt_nb_2_grouped_ls_", i, ".rds")))