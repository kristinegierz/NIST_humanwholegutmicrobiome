# Stability Analysis - METAGENOMICS
# Author: Kristine Gierz
# Date: August 20, 2024

# This is the code file for the stability analysis for metagenomics.

# --- Source in helper files ----

source('utils.R')
source('makedata.R')

# ---- Load needed libraries ----

using('tidyverse',
      'janitor',
      'readxl',
      'vegan',
      'QuantileNPCI',
      'patchwork')

# ---- Bray-Curtis & Hellinger Distances ----

stability_wide <- stability_seq %>%
  select(!asv_id) %>%
  pivot_wider(names_from = new_id,
              values_from = rel_abund)

stability_use <- stability_wide %>%
  dplyr::select(!c(group, box, replicate, timepoint))

distances <- vegdist(as.matrix(stability_use[,-1]), method = 'bray')
distances_hellinger <- vegdist(as.matrix(stability_use[,-1]), method = "hellinger")

mds <- cmdscale(distances) %>% as_tibble()
colnames(mds) <- c("Dim1", "Dim2")

mds_hellinger <- cmdscale(distances_hellinger) %>% as_tibble()
colnames(mds_hellinger) <- c("Dim1", "Dim2")

mdsfull <- mds %>%
  mutate(sample = stability_use$sample) %>%
  left_join(stability_wide)

mdsfull_hellinger <- mds_hellinger %>%
  mutate(sample = stability_use$sample) %>%
  left_join(stability_wide)

plot <- ggplot(data = mdsfull %>% 
                 filter(group %in% c("Omnivore", "Vegetarian")) %>%
                 filter(!is.na(group)), aes(x = Dim1, y = Dim2, 
                                            col = as.factor(timepoint),
                                            shape = group)) +
  geom_point(size = 3,
             position = position_jitter(),
             alpha = 0.7,
             show.legend = F) +
  labs(color = "Date",
       title = "Stability MDS Plot for Bray-Curtis Distance")

plot

plot_h <- ggplot(data = mdsfull_hellinger %>% 
                 filter(group %in% c("Omnivore", "Vegetarian")) %>%
                 filter(!is.na(group)), aes(x = Dim1, y = Dim2, 
                                            col = as.factor(timepoint),
                                            shape = group)) +
  geom_point(size = 3,
             position = position_jitter(),
             alpha = 0.7) +
  labs(color = "Date",
       title = "Stability MDS Plot for Hellinger Distance")

plot_h

combined_plot <- plot + 
  plot_h +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'right')
combined_plot

# ---- Pairwise Differences ----

## ---- Omnivore ----

omnivore <- df$data[[1]]
test <- omnivore %>%
  filter(orderID < 301)

### ---- Differences ----

## For each pairwise difference of days/vials/reps we need to get a 
## set of percent differences for each orderID. So for a comparison between
## omnivore day 1 vial 1 rep 1 and omnivore day 1 vial 2 rep 1 we should have
## a collection of 300 pairwise percent differences. From this, we will find the
## 99th percentile to get to the Coefficient of Disagreement. quantCI can get
## us an approximate 95% CI on the COD. 

## We will combine omnivore/vegetarian together after this

threshold <- test %>%
  group_by(orderID) %>%
  summarise(avg_abund = mean(rel_abund)) %>%
  ungroup() %>%
  mutate(rel_abund = avg_abund) %>%
  dplyr::select(!avg_abund) %>%
  filter(rel_abund > 0.0008)

means <- test %>%
  filter(orderID < 122) %>%
  group_by(orderID, day, vial) %>%
  summarise(avg_abund = mean(rel_abund)) %>%
  ungroup() %>%
  mutate(rel_abund = avg_abund) %>%
  dplyr::select(!avg_abund)

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
# [1] 0.002089864

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
lookup <- cbind(unique(omnivore$box), unique(omnivore$vial)) %>%
  as.data.frame() %>%
  rename(box = V1, vial = V2) %>%
  mutate(across(everything(), as.numeric)) %>%
  arrange(vial)

# Bind the original and duplicated data together
coef_dis_plot <- bind_rows(coef_dis_results, duplicated_data) %>%
  separate_wider_delim(pair1, delim = "_", names = c("day", "vial")) %>%
  separate_wider_delim(pair2, delim = "_", names = c("day2", "vial2")) %>%
  mutate(vial = as.numeric(vial)) %>%
  left_join(lookup, by = "vial") %>%
  rename_pair() %>%
  mutate(day = factor(day2, levels = c(1, 10, 16, 31, 59, 87, 112, 140, 168)))

plot2 <- ggplot(coef_dis_plot, aes(x = group,
                          y = coef_dis, 
                          col = factor(day))) +
  geom_point() +
  ylim(c(0, 0.15)) +
  # geom_hline(yintercept = 0.07445026, color = 'red') +
  ylab("Coefficient of Disagreement") +
  xlab("Sample") +
  labs(col = "Time Point") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot2
### ---- Confidence Intervals ----

# First with quantCI approximation

upper <- quantCI(coef_dis_results$coef_dis, 0.99, 0.05, method = 
                   'approximate')$upper.ci
print(upper)
# 99.705983302224th percentile  
# 0.1140143 

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
# 36.73906%  99.85689% 
#   0.03746716 0.11367601 

## ---- Vegetarians ----

vegetarian <- df$data[[2]]
test <- vegetarian %>%
  filter(orderID < 301)

### ---- Differences ----

threshold <- test %>%
  group_by(orderID) %>%
  summarise(avg_abund = mean(rel_abund)) %>%
  ungroup() %>%
  mutate(rel_abund = avg_abund) %>%
  dplyr::select(!avg_abund) %>%
  filter(rel_abund > 0.0008)

means <- test %>%
  filter(orderID < 100) %>%
  group_by(orderID, day, vial) %>%
  summarise(avg_abund = mean(rel_abund)) %>%
  ungroup() %>%
  mutate(rel_abund = avg_abund) %>%
  dplyr::select(!avg_abund)

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
# [1] 0.002983683

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
  separate_wider_delim(pair1, delim = "_", names = c("day", "vial")) %>%
  separate_wider_delim(pair2, delim = "_", names = c("day2", "vial2")) %>%
  mutate(vial = as.numeric(vial)) %>%
  left_join(lookup, by = "vial") %>%
  rename_pair() %>%
  mutate(day = factor(day2, levels = c(1, 10, 16, 31, 59, 87, 112, 140, 168)))

plot2 <- ggplot(coef_dis_plot, aes(x = group,
                          y = coef_dis, 
                          col = day)) +
  geom_point() +
  ylim(c(0, 0.20)) +
  # geom_hline(yintercept = 0.07445026, color = 'red') +
  ylab("Coefficient of Disagreement") +
  xlab("Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot2

### ---- Confidence Intervals ----

# First with quantCI approximation

upper <- quantCI(coef_dis_results$coef_dis, 0.99, 0.05, method = 
                   'approximate')
quantile(coef_dis_results$coef_dis, upper$u2)
# 99.7737% 
# 0.1504442

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
# 30.10079%  99.74662% 
#   0.04904131 0.16027750 
