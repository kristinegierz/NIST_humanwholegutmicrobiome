## Making the Data - METAGENOMICS
## Author: Kristine Gierz
## Date: 02/09/2023

## This file makes the data for the analysis.

#### ---- Source in the helper functions, etc ----

source("utils.R")

#### ---- Load Libraries ----

using('tidyverse',
      'janitor',
      'readxl') 

#### ---- Load in the data ----

# Read in the data
# change pathway as appropriate
homogenlong <- read.delim("../../../02 Processed Data/Homogeneity/Taxa.tab.meta.txt")

# Clean variable names and take out vegetarian sample 11 (experimental error)
homogenlong <- homogenlong %>%
  clean_names() %>%
  filter(sample != 'V11',
         sample != 'V37') 

# Create the sequences
sequences <- data.frame(sequence = unique(homogenlong$sequence),
                        new_id = 1:1066)

# Match up sequences and new id
hseq <- homogenlong %>%
  left_join(sequences) %>%
  dplyr::select(asv_id, sample, group, rel_abund, 
                box, replicate, sequence, new_id)

# Nested dataframe by experimental group
df <- hseq %>%
  nest(.by = group)

# Apply the function to each element of the nested data frame
df <- df %>%
  mutate(data = map(data, bayes_restructure))

