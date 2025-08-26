## Make Data - METAGENOMICS
## Author: Kristine Gierz
## Date: Aug 20, 2024

## This file creates the data for analysis for metagenomics stability.

# ---- Source in Helpers ----

source('utils.R')

# ---- Load Additional Libraries ----

using('tidyverse',
      'janitor',
      'readxl')

# ---- Load the Data ----

# The homogeneity timepoint
homogenlong <- read.delim("../../../02 Processed Data/Homogeneity/Taxa.tab.meta.txt")
stabilitylong <- read.csv("../data/Taxa.tab3t12345678.meta.csv")

# ---- Cleaning Etc ----

homogen <- homogenlong %>% 
  clean_names() %>%
  filter(sample != "V11",
         sample != "V37") %>%
  dplyr::select(!c(operator, udi)) %>%
  filter(group %in% c("Omnivore", "Vegetarian")) %>%
  mutate(timepoint = "T0",
         box = as.numeric(box),
         day = 1) 

stability <- stabilitylong %>%
  clean_names() %>%
  dplyr::select(!c(sample_id, 
                   sequencing_run_date, 
                   sequencing_run, 
                   sequencing_code)) %>%
  filter(group %in% c("Omnivore", "Vegetarian")) %>%
  mutate(box = as.numeric(box))

fulldf <- rbind(homogen, stability) %>%
  as_tibble()

# Create the sequences
sequences <- data.frame(sequence = unique(stabilitylong$Sequence),
                        new_id = 1:2353)

# Match up sequences and new id - for the distance plots
stability_seq <- stabilitylong %>%
  clean_names() %>%
  left_join(sequences) %>%
  mutate(timepoint = 
           as_date(ifelse(timepoint == 'T1', 
                          mdy('6/1/2023'),
                          ifelse(timepoint == 'T2',
                                 mdy('6/7/2023'),
                                 ifelse(timepoint == 'T3', 
                                        mdy('6/22/2023'),
                                        ifelse(timepoint == 'T4',
                                               mdy('7/20/2023'),
                                               ifelse(timepoint == 'T5',
                                                      mdy('8/17/2023'),
                                                      ifelse(timepoint == 'T6',
                                                             mdy('9/11/2023'),
                                                             ifelse(timepoint == 'T7',
                                                                    mdy('10/9/2023'),
                                                                    mdy('11/6/2023')
                                                             )
                                                      )
                                               )
                                        )
                                 )
                          )
           )
           )) %>%
  dplyr::select(asv_id, sample, group, rel_abund, 
                box, replicate, new_id, timepoint)

# Nest your data - easier to restructure for the pairwise comparisons
df <- fulldf %>%
  nest(.by = group)

# Assuming your nested data frame is named nested_df
# Example of nesting, assuming you start with a data frame named df
# nested_df <- df %>% nest(data = -some_grouping_columns)

# Apply the function to each element of the nested data frame
df <- df %>%
  mutate(data = map(data, bayes_restructure))

