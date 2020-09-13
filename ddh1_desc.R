library(tidyverse)
library(readr)
library(dplyr)

smiles <- read_csv("smiles.csv") %>%
  mutate(pIC50 = ifelse(pIC50 == "BLINDED", NA, pIC50)) %>%
  mutate(pIC50 = as.numeric(pIC50))

desc_files <- c("padel_descriptors_1d_2d.csv", "blue_desc.csv", "chemopy_desc.csv", "pybel_desc.csv", "rdkit_desc.csv")
descs <- lapply(desc_files, FUN = function(desc_file) {
  df <- read_csv(paste0("descriptors/", desc_file)) %>%
    select(where(is.numeric)) %>%
    mutate(id = 1:nrow(.)) %>%
    relocate(id) %>%
    mutate(pIC50 = smiles$pIC50[id]) %>%
    relocate(pIC50, .after = last_col()) %>%
    filter(!is.na(pIC50))
  write_csv(df, path = paste0("output/", desc_file))
  df
})




