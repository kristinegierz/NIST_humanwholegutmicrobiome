## Helper Functions Required - METAGENOMICS
## Author: Kristine Gierz
## Date: Aug 20, 2024

## This file has all helper functions needed for the analysis.

# ---- Loading/Installing New Libraries ----

using <- function(...) {
  libs <- unlist(list(...))
  req <- unlist(lapply(libs, require, character.only = TRUE))
  need <- libs[req == FALSE]
  if (length(need) > 0) { 
    install.packages(need)
    lapply(need, require, character.only = TRUE)
  }
}

# ---- Data Restructuring Helpers ----

# Assuming your nested data frame is named nested_df
# Example of nesting, assuming you start with a data frame named df
# nested_df <- df %>% nest(data = -some_grouping_columns)

# Define a function to calculate orderID for a single nested data frame
bayes_restructure <- function(data) {
  data %>%
    group_by(sequence) %>%                          # Group by the DNA sequence column
    summarise(avg_rel_abund = mean(rel_abund)) %>%  # Calculate the average rel_abund
    arrange(desc(avg_rel_abund)) %>%                # Arrange in descending order
    mutate(orderID = row_number()) %>%              # Create the orderID column
    select(sequence, orderID) %>%                   # Select only the necessary columns
    right_join(data, by = "sequence") %>%           # Join back with the original data
    mutate(vial = dense_rank(box))                  # Renumber the boxes to vials
}

rename_pair <- function(df) {
  # Split the pair into components and create a new column with renamed pairs
  df <- df %>%
    mutate(
      new_pair = paste0("D", day, "B", box)
    )
  
  # Sort the unique pairs
  sorted_unique_pairs <- df %>%
    mutate(day = as.numeric(day)) %>%
    distinct(day, box, new_pair) %>%
    arrange(day, box) %>%
    pull(new_pair)
  
  # Convert new_pair to a factor with the levels in the sorted order
  df <- df %>%
    mutate(group = factor(new_pair,
                          levels = sorted_unique_pairs)) %>%
    select(-day, -box, -new_pair, -vial) # Remove the temporary columns
  
  return(df)
}

# ---- Pairwise Differences Helpers ----

# Function to calculate pairwise differences and include day/vial information
pairwise_diff_with_pairs <- function(df, epsilon = 1e-6) {
  # Ensure there are at least two observations to compare
  if(nrow(df) < 2) {
    return(data.frame(
      pair1_day = numeric(0),
      pair1_vial = numeric(0),
      pair2_day = numeric(0),
      pair2_vial = numeric(0),
      pairwise_diff = numeric(0)
    ))
  }
  
  # Get all combinations of rows without repetition
  combn(seq_len(nrow(df)), 2, function(indices) {
    i1 <- indices[1]
    i2 <- indices[2]
    value1 <- df$rel_abund[i1] 
    value2 <- df$rel_abund[i2] 
    
    percent_difference <- ifelse(value1 == 0 & value2 == 0, 0,
                                 abs((log(value1) - log(value2))/
                                       (log(value1) + log(value2)) ))
    
    data.frame(
      pair1 = paste(df$day[i1], df$vial[i1],
                    sep = "_"),
      pair2 = paste(df$day[i2], df$vial[i2],
                    sep = "_"),
      percent_diff = percent_difference
    )
  }, simplify = FALSE) %>%
    bind_rows()
}

# ---- Bootstrap Helpers ----

# Function to calculate the CoD for a single bootstrap sample
bootstrap_single_sample_cod <- function(data) {
  sampled_pair <- data %>%
    distinct(pair1, pair2) %>%
    sample_n(size = 1, replace = TRUE)
  
  resampled_data <- sampled_pair %>%
    inner_join(data, by = c("pair1", "pair2"))
  
  cod_value <- quantile(resampled_data$percent_diff, 0.99, na.rm = TRUE)
  
  return(cod_value)
}

bootstrap_cod <- function(data, n_bootstrap = 1000) {
  bootstrap_samples <- replicate(n_bootstrap, {
    bootstrap_single_sample_cod(data)
  }, simplify = TRUE)
  
  return(bootstrap_samples)
}
