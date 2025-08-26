## Analysis File - METAGENOMICS
## Author: Kristine Gierz
## Date: 02/09/2023

## This file contains EDA, distance and clustering methods, and Bayesian
## analysis for the metagenomics data. 

## Check the directory to make sure you are reading in everything right
# getwd()

## THERE STILL NEEDS TO BE SOME EDITS IN THIS FOR FILE PATHWAYS - NOT
## IMMEDIATELY FULLY FUNCTIONAL

# ---- Source Helper Files ----

source("utils.R")
source("makedata.R")

# ---- Load Needed Libraries ----

using("tidyverse",
      "janitor",
      "readxl",
      "vegan",
      'rstan',
      'pbapply',
      'QuantileNPCI',
      'patchwork') #observe start up messages

options(mc.cores = parallel::detectCores())

rstan_options(auto_write = TRUE)

# ---- Bray-Curtis & Hellinger Distances ----

## All together - need to put data in wide format
hwide <- hseq %>%
  select(!c(sequence, new_id)) %>% 
  pivot_wider(names_from = asv_id,
              values_from = rel_abund) %>%
  mutate(box = ifelse(group == "water", 100, 
                      ifelse(group == "zymo", 200, box))) 

# Take out the pooled samples - we won't be looking at those. For distance
# matrix, we need to take out group, box, replicate info for now
huse <- hwide %>%
  # filter(group %in% c("Omnivore", "Vegetarian", "water", "zymo")) %>%
  dplyr::select(!c(group, box, replicate))

# We are using Bray-Curtis distance, but this data is also applicable to using
# a Hellinger distance. The Bray-Curtis distance is preferred by Stephanie and
# Sam due to the fact that is more commonly used in the field. 

# We don't need the first column that gives sample info
distances <- vegdist(as.matrix(huse[,-1]), method = 'bray')
distances_hellinger <- vegdist(as.matrix(huse[,-1]), method = 'hellinger')

# We are using multi-dimensional scaling as dimension reduction
mds <- cmdscale(distances) %>% as_tibble()
colnames(mds) <- c("Dim1", "Dim2")

mds_hellinger <- cmdscale(distances_hellinger) %>% as_tibble()
colnames(mds_hellinger) <- c("Dim1", "Dim2")

# Combine back with the full data set for ease of plotting
mdsfull <- mds %>%
  mutate(sample = huse$sample) %>%
  left_join(hwide)

mdsfull_hellinger <- mds_hellinger %>%
  mutate(sample = huse$sample) %>%
  left_join(hwide)
## ---- Distance Plot ----
plot <- ggplot(data = mdsfull %>%
                 filter(group %in% c("Omnivore",
                                      "Vegetarian")) %>%
                 mutate(box = factor(box, levels = c(1, 6, 11, 15, 20, 25, 30, 36, 44, 50))), 
               aes(x = Dim1, 
                   y = Dim2, 
                   # col = box,
                   col = group)) +
  geom_point(size = 3,
             position = position_jitter(),
             alpha = 0.75,
             show.legend = F) +
  labs(col = "Box",
       shape = "Group",
       title = "Homogeneity MDS Plot for Bray-Curtis Distance") 

plot

plot_h <- ggplot(data = mdsfull_hellinger %>%
                 filter(group %in% c("Omnivore",
                                     "Vegetarian")) %>%
                 mutate(box = factor(box, levels = c(1, 6, 11, 15, 20, 25, 30, 36, 44, 50))), 
               aes(x = Dim1, 
                   y = Dim2, 
                   # col = box,
                   col = group)) +
  geom_point(size = 3,
             position = position_jitter(),
             alpha = 0.75) +
  labs(col = "Group",
       shape = "Group",
       title = "Homogeneity MDS Plot for Hellinger Distance")

plot_h

combined_plot <- plot + 
  plot_h +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'right') 
combined_plot

# ---- Bayesian Models ----

## I WOULD RECOMMEND RUNNING THIS IN AN HPC OR ONLY RUNNING A PORTION IF YOU
## WANT TO VERIFY RESULTS. IT WILL TAKE QUITE SOME TIME AS THE RELATIVE 
## ABUNDANCE DECREASES IN THE ORDERED IDS. YOU SHOULD ALSO MAKE SURE YOUR
## DIRECTORY STRUCTURE IS SET UP CORRECTLY. 

## I ran this in a shell script over rhea.cam.nist.gov and saved results. For
## this many models, effectively storing all model information is not as easily
## tenable as it could be. If you want to verify results, the seed is set in the 
## fit_and_save function. You should be able to verify any subset of random 
## models that you would like to. Results are all saved. 

## This is the STAN code - you could also read this in as a .stan file

## ---- STAN Model ----
stan_model_code <- "
data {
  int<lower=0> N; // total number of observations
  int<lower=0> J; // number of vials
  int<lower=0> K; // number of replicates per vial
  int<lower=1,upper=J> vial[N]; // vial indicator
  vector<lower=0>[N] y; // outcome variable (relative abundance)
}

parameters {
  real<lower=0> mu; // overall mean
  real<lower=0> sigma_v; // standard deviation of vial effects
  real<lower=0> sigma_e; // standard deviation of replicate errors
  vector[J] theta; // vial effects
}

model {
  // Priors
  sigma_v ~ student_t(4, 0, 1); // Half Student's t with 4 df
  sigma_e ~ student_t(4, 0, 1); // Half Student's t with 4 df
  theta ~ normal(0, sigma_v);

  // Likelihood
  for (n in 1:N) {
    y[n] ~ normal(mu + theta[vial[n]], sigma_e);
  }
}

generated quantities {
  real<lower=0> rho;
  real<lower=0> cov;
  rho = sigma_v / sigma_e;
  cov = sigma_v / mu;
}

"

# We need to separate this into omnivores and vegetarians

### ---- Omnivores ----

# List of datasets, each corresponding to a different ASV ID
# This is the omnivore data
# uncomment to run
# data <- df$data[[1]]

# We need to nest this by the ordered IDs of sequences
# uncomment to run
# dataset_list <- data %>%
#   nest(.by = orderID)

# Iterate over each nested dataset in the list
# uncomment to run
# for (i in seq_along(dataset_list$data)) {
#   fit_and_save(data = dataset_list$data[[i]], 
#                asv_id = dataset_list$orderID[i],
#                group = "omnivore",
#                model_code = stan_model_code)
# }

### ---- Vegetarians ----

# List of datasets, each corresponding to a different ASV ID
# This is the vegetarian data

# uncomment to run
# data <- df$data[[5]]

# We need to nest this by the ordered IDs of sequences
# uncomment to run 
# dataset_list <- data %>%
#   nest(.by = orderID)

# Iterate over each nested dataset in the list
# uncomment to run
# for (i in seq_along(dataset_list$data)) {
#   fit_and_save(data = dataset_list$data[[i]], 
#                asv_id = dataset_list$orderID[i],
#                group = "vegetarian",
#                model_code = stan_model_code)
# }

## ---- Plots ----

## We need to create the plots for the analysis - for how we are defining
## homogeneity, we need the ratio of variances with their credible intervals
## and we need the coefficient of variation with its credible interval. 

### ---- Omnivores ----

file_path <- "../output/results/omnivore"

# List all RData files in the folder
files <- list.files(path = file_path, 
                    pattern = "\\.Rds$", 
                    full.names = TRUE)

#### ---- Rho plot ----

# Extract stats for each file for a specific element (e.g., "rho")
results <- pblapply(files, extract_element_stats, element_name = "rho")

# Create a dataframe
results_df <- do.call(rbind, results) %>%
  as.data.frame() %>%
  mutate(file_name = gsub(file_path, "", files)) %>%
  mutate(s_number = as.numeric(gsub(".*_(\\d+)\\.Rds", "\\1", file_name))) %>%
  arrange(s_number) %>%
  rename(lower = 'lower.2.5%', upper = 'upper.97.5%')

# Plot using ggplot2
plot <- ggplot(results_df %>% filter(s_number <= 300), aes(x = s_number, y = mean)) +
  geom_point(aes(color = ifelse(lower > 1 & upper > 1, 
                                "above_1", 
                                ifelse(lower <= 1 & upper > 1, 
                                       "above_1_below_1", 
                                       "below_1")))) + # Color dots based on condition
  geom_errorbar(aes(ymin = lower, 
                    ymax = upper, 
                    color = ifelse(lower > 1 & upper > 1,
                                   "above_1",
                                   ifelse(lower <= 1 & upper > 1, 
                                          "above_1_below_1", 
                                          "below_1"))), 
                width = 0.2) + # Set color for error bars
  geom_hline(yintercept = 1, 
             linetype = "dashed", 
             color = "red") +
  labs(x = "Sequences (order of relative abundance)", 
       y = "Variance Ratio",
       title = "Ratios and Intervals for 300\nMost Abundant Seqs - Omn") + 
  scale_color_manual(values = c(above_1 = "red", 
                                above_1_below_1 = "blue",
                                below_1 = "black"),
                     labels = c("Above 1", 
                                "Straddle 1", 
                                "Below 1"),
                     name = NULL) + # Define colors and hide legend
  theme_minimal() +
  theme(axis.text.x = element_blank()) # Remove x-axis labels

plot

#### ---- COV plot ----

covresults <- pblapply(files, extract_element_stats, element_name = "cov")

# Create a dataframe

results_c_df <- do.call(rbind, covresults) %>%
  as.data.frame() %>%
  mutate(file_name = gsub(file_path, "", files)) %>%
  mutate(s_number = as.numeric(gsub(".*_(\\d+)\\.Rds", "\\1", file_name))) %>%
  arrange(s_number) %>%
  rename(lower = 'lower.2.5%', upper = 'upper.97.5%')


plot <- ggplot(results_c_df %>% filter(s_number <= 300), aes(x = s_number, y = mean)) +
  geom_point() + 
  geom_errorbar(aes(ymin = lower, 
                    ymax = upper), 
                width = 0.2) +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "red") +
  labs(x = "Sequences (order of relative abundance)", 
       y = "Coefficient of Variation",
       title = "CV for 300\nMost Abundant Seqs - Omn") +
  theme_minimal() +
  ylim(0, 5) +
  theme(axis.text.x = element_blank()) # Remove x-axis labels

plot

### ---- Vegetarians ----

file_path <- "../output/results/vegetarian"

# List all RData files in the folder
files <- list.files(path = file_path, 
                    pattern = "\\.Rds$", 
                    full.names = TRUE)

#### ---- Rho plot ----

# Extract stats for each file for a specific element (e.g., "rho")
results <- pblapply(files, extract_element_stats, element_name = "rho")

# Create a dataframe
results_df <- do.call(rbind, results) %>%
  as.data.frame() %>%
  mutate(file_name = gsub(file_path, "", files)) %>%
  mutate(s_number = as.numeric(gsub(".*_(\\d+)\\.Rds", "\\1", file_name))) %>%
  arrange(s_number) %>%
  rename(lower = 'lower.2.5%', upper = 'upper.97.5%')

# Plot using ggplot2
plotv <- ggplot(results_df %>% filter(s_number <= 300), aes(x = s_number, y = mean)) +
  geom_point(aes(color = ifelse(lower > 1 & upper > 1, 
                                "above_1", 
                                ifelse(lower <= 1 & upper > 1, 
                                       "above_1_below_1", 
                                       "below_1")))) + # Color dots based on condition
  geom_errorbar(aes(ymin = lower, 
                    ymax = upper, 
                    color = ifelse(lower > 1 & upper > 1,
                                   "above_1",
                                   ifelse(lower <= 1 & upper > 1, 
                                          "above_1_below_1", 
                                          "below_1"))), 
                width = 0.2) + # Set color for error bars
  geom_hline(yintercept = 1, 
             linetype = "dashed", 
             color = "red") +
  labs(x = "Sequences (order of relative abundance)", 
       y = "Variance Ratio",
       title = "Ratios and Intervals for 300\nMost Abundant Seqs - Veg") + 
  scale_color_manual(values = c(above_1 = "red", 
                                above_1_below_1 = "blue",
                                below_1 = "black"),
                     labels = c("Above 1", 
                                "Straddle 1", 
                                "Below 1"),
                     name = NULL) + # Define colors and hide legend
  theme_minimal() +
  theme(axis.text.x = element_blank()) # Remove x-axis labels

plotv

#### ---- COV plot ----

covresults <- pblapply(files, extract_element_stats, element_name = "cov")

# Create a dataframe

results_c_df <- do.call(rbind, covresults) %>%
  as.data.frame() %>%
  mutate(file_name = gsub(file_path, "", files)) %>%
  mutate(s_number = as.numeric(gsub(".*_(\\d+)\\.Rds", "\\1", file_name))) %>%
  arrange(s_number) %>%
  rename(lower = 'lower.2.5%', upper = 'upper.97.5%')


plotv <- ggplot(results_c_df %>% filter(s_number <= 300), aes(x = s_number, y = mean)) +
  geom_point() + 
  geom_errorbar(aes(ymin = lower, 
                    ymax = upper), 
                width = 0.2) +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "red") +
  labs(x = "Sequences (order of relative abundance)", 
       y = "Coefficient of Variation",
       title = "CV for 300\nMost Abundant Seqs - Veg") +
  theme_minimal() +
  ylim(0, 5) +
  theme(axis.text.x = element_blank()) # Remove x-axis labels

plotv

# ---- Pairwise Comparisons ----

## An alternate way to guage the data variability - we look at
## all pairwise percent differences between the reps for the 
## 300 most abundant sequences.

## ---- Omnivores ----

omnivore <- df$data[[1]]
omnivore <- omnivore %>%
  filter(orderID < 301)

### ---- Pairwise differences ----

threshold <- omnivore %>%
  group_by(orderID) %>%
  summarise(avg_abund = mean(rel_abund)) %>%
  ungroup() %>%
  mutate(rel_abund = avg_abund) %>%
  select(!avg_abund) %>%
  filter(rel_abund > 0.0008)

means <- omnivore %>%
  filter(orderID < 125) %>%
  group_by(orderID, vial, replicate) %>%
  summarise(avg_abund = mean(rel_abund)) %>%
  ungroup() %>%
  mutate(rel_abund = avg_abund) %>%
  select(!avg_abund)

pairwise_results <- means %>%
  group_by(orderID) %>%
  do(pairwise_diff_with_pairs(., epsilon = 1e-8)) %>%
  ungroup()

# We will get INF for log(0) if a sequence was present at a low level
# and during a measurement became "absent"
nacheck <- pairwise_results %>%
  filter(is.na(percent_diff))

# To get the percent of values that won't be considered in the 
# 99th percentile
nrow(nacheck) / nrow(pairwise_results)
# [1] 0.002915633

coef_dis_results <- pairwise_results %>%
  mutate(across(percent_diff, ~ ifelse(is.finite(.), ., NA))) %>%
  group_by(pair1, pair2) %>%
  summarise(coef_dis = quantile(percent_diff, 0.99, na.rm = T),
            .groups = 'drop') 

### ---- Plots ----

# Duplicate each row, swapping pair1 and pair2
duplicated_data <-  coef_dis_results %>%
  mutate(
    pair1_temp = pair1,
    pair1 = pair2,
    pair2 = pair1_temp
  ) %>%
  select(-pair1_temp)

# Create the data frame from the matrix
lookup <- cbind(unique(omnivore$box), unique(omnivore$vial)) %>%
  as.data.frame() %>%
  rename(box = V1, vial = V2) %>%
  mutate(across(everything(), as.numeric)) %>%
  arrange(vial)

# Bind the original and duplicated data together
coef_dis_plot <- bind_rows(coef_dis_results, duplicated_data) %>%
  separate_wider_delim(pair1, delim = "_", names = c("vial", "rep")) %>%
  separate_wider_delim(pair2, delim = "_", names = c("vial2", "rep2")) %>%
  mutate(vial = as.numeric(vial)) %>%
  left_join(lookup, by = "vial") %>%
  rename_pair() %>%
  mutate(vial = as.numeric(vial2)) %>%
  left_join(lookup, by = 'vial') %>%
  mutate(box = factor(box, levels = c(1, 6, 11, 15, 20, 25, 30, 36, 44, 50))) 

plot <- ggplot(coef_dis_plot, aes(x = group,
                          y = coef_dis, 
                          col = as.factor(vial))) +
  geom_point() +
  ylim(c(0, 0.10)) +
  # geom_hline(yintercept = 0.07445026, color = 'red') +
  ylab("Coefficient of Disagreement") +
  xlab("Sample") +
  labs(col = "Vial") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot
### ---- Confidence Intervals ----

# We can estimate an approximate interval using the quantCI function. Since
# the values are actually dependent and this method assumes independence, this
# is truly only an approximation. 

ci <- quantCI(coef_dis_results$coef_dis, 0.99, 0.05, method = 'approximate')
print(ci$upper.ci)
# 99.5720557487675th percentile  
# 0.07267177 

# Alternately, we can do a bootstrap confidence interval
# However, this is currently commented out due to the fact that this will,
# essentially, be exactly the interval above since we do not have enough 
# samples to resample from. 

# Perform bootstrap
# set.seed(123)
# n_bootstrap <- 1000
# 
# original_stat <- quantile(pairwise_results$percent_diff, 0.99, na.rm = TRUE)
# 
# # Calculate z0
# bootstrap_samples <- bootstrap_cod(pairwise_results, n_bootstrap = 1000)
# z0 <- qnorm(mean(bootstrap_samples < original_stat))
# 
# # Jackknife method to calculate acceleration factor (a)
# jackknife_samples <- sapply(1:nrow(pairwise_results), function(i) {
#   jackknife_data <- pairwise_results[-i, ]
#   quantile(jackknife_data$percent_diff, 0.99, na.rm = TRUE)
# })
# 
# mean_jackknife <- mean(jackknife_samples)
# acceleration <- sum((mean_jackknife - jackknife_samples)^3) /
#   (6 * sum((mean_jackknife - jackknife_samples)^2)^(3/2))
# 
# # Adjusted percentiles
# alpha <- 0.95
# z_alpha <- qnorm(alpha)
# p1 <- pnorm(z0 + (z0 + z_alpha) / (1 - acceleration * (z0 + z_alpha)))
# p2 <- pnorm(z0 + (z0 - z_alpha) / (1 - acceleration * (z0 - z_alpha)))
# 
# upper_percentile <- quantile(bootstrap_samples, p1, na.rm = TRUE)
# lower_percentile <- quantile(bootstrap_samples, p2, na.rm = TRUE)
# 
# bc_ci <- c(lower_percentile, upper_percentile)
# 
# print(bc_ci)
# 34.84617%  99.82628% 
#   0.03454683 0.07265898 

## ---- Vegetarians ----

vegetarian <- df$data[[5]]
vegetarian <- vegetarian %>%
  filter(orderID < 301) 

### ---- Pairwise differences ----

threshold <- vegetarian %>%
  group_by(orderID) %>%
  summarise(avg_abund = mean(rel_abund)) %>%
  ungroup() %>%
  mutate(rel_abund = avg_abund) %>%
  select(!avg_abund) %>%
  filter(rel_abund > 0.0008)

means <- vegetarian %>%
  filter(orderID < 105) %>%
  group_by(orderID, vial, replicate) %>%
  summarise(avg_abund = mean(rel_abund)) %>%
  ungroup() %>%
  mutate(rel_abund = avg_abund) %>%
  select(!avg_abund)

pairwise_results <- means %>%
  group_by(orderID) %>%
  do(pairwise_diff_with_pairs(., epsilon = 1e-8)) %>%
  ungroup()

# We will get INF for log(0) if a sequence was present at a low level
# and during a measurement became "absent"
nacheck <- pairwise_results %>%
  filter(is.na(percent_diff))

# To get the percent of values that won't be considered in the 
# 99th percentile
nrow(nacheck) / nrow(pairwise_results)
# [1] 0.001012146

coef_dis_results <- pairwise_results %>%
  group_by(pair1, pair2) %>%
  summarise(coef_dis = quantile(percent_diff, 0.99, na.rm = T),
            .groups = 'drop') 

### ---- Plots ----

# Duplicate each row, swapping pair1 and pair2
duplicated_data <-  coef_dis_results %>%
  mutate(
    pair1_temp = pair1,
    pair1 = pair2,
    pair2 = pair1_temp
  ) %>%
  select(-pair1_temp)

# Create the data frame from the matrix
lookup <- cbind(unique(vegetarian$box), unique(vegetarian$vial)) %>%
  as.data.frame() %>%
  rename(box = V1, vial = V2) %>%
  mutate(across(everything(), as.numeric)) %>%
  arrange(vial)

# Bind the original and duplicated data together
coef_dis_plot <- bind_rows(coef_dis_results, duplicated_data) %>%
  separate_wider_delim(pair1, delim = "_", names = c("vial", "rep")) %>%
  separate_wider_delim(pair2, delim = "_", names = c("vial2", "rep2")) %>%
  mutate(vial = as.numeric(vial)) %>%
  left_join(lookup, by = "vial") %>%
  rename_pair() %>%
  mutate(vial = as.numeric(vial2)) %>%
  left_join(lookup, by = 'vial') %>%
  mutate(box = factor(box, levels = c(1, 6, 11, 15, 20, 25, 30, 36, 44, 50)))

ggplot(coef_dis_plot, aes(x = group,
                          y = coef_dis,
                          col = box)) +
  geom_point() +
  ylim(c(0, 0.25)) +
  ylab("Coefficient of Disagreement") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### ---- Confidence Intervals ----

ci <- quantCI(coef_dis_results$coef_dis, 0.99, 0.05, method = 'approximate')
print(ci$upper.ci)
# 99.5951703995335th percentile  
# 0.1057341 

# Alternately, we can do a bootstrap confidence interval

# Perform bootstrap
# set.seed(123)
# n_bootstrap <- 1000
# 
# original_stat <- quantile(pairwise_results$percent_diff, 0.99, na.rm = TRUE)
# 
# # Calculate z0
# bootstrap_samples <- bootstrap_cod(pairwise_results, n_bootstrap = 1000)
# z0 <- qnorm(mean(bootstrap_samples < original_stat))
# 
# # Jackknife method to calculate acceleration factor (a)
# jackknife_samples <- sapply(1:nrow(pairwise_results), function(i) {
#   jackknife_data <- pairwise_results[-i, ]
#   quantile(jackknife_data$percent_diff, 0.99, na.rm = TRUE)
# })
# 
# mean_jackknife <- mean(jackknife_samples)
# acceleration <- sum((mean_jackknife - jackknife_samples)^3) /
#   (6 * sum((mean_jackknife - jackknife_samples)^2)^(3/2))
# 
# # Adjusted percentiles
# alpha <- 0.95
# z_alpha <- qnorm(alpha)
# p1 <- pnorm(z0 + (z0 + z_alpha) / (1 - acceleration * (z0 + z_alpha)))
# p2 <- pnorm(z0 + (z0 - z_alpha) / (1 - acceleration * (z0 - z_alpha)))
# 
# upper_percentile <- quantile(bootstrap_samples, p1, na.rm = TRUE)
# lower_percentile <- quantile(bootstrap_samples, p2, na.rm = TRUE)
# 
# bc_ci <- c(lower_percentile, upper_percentile)
# 
# print(bc_ci)
# 46.18118%  99.93706% 
#  0.04750851 0.10995658  
