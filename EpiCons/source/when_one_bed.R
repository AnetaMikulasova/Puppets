args <- commandArgs(trailing = TRUE)

library(dplyr)
library(tidyverse)
library(data.table)

INPUT       = args[1] %>% as.character()
OUTPUT      = args[2] %>% as.character()

# INPUT = "/Volumes/Aneta_10TB/BCK/Aneta_data_bck/190121/projects/IGH_project/blueprint_chromatin-state_data-segmentation/results/x_multiintersect/multiintersect_BcellPC_E9_chr14_106000000_106400000_n=15.bed"
# FRACTION = "0.05" %>% as.numeric()
# NofSAMPLES = "15" %>% as.numeric()
# OUTPUT = "/Volumes/Aneta_10TB/BCK/Aneta_data_bck/190121/projects/IGH_project/blueprint_chromatin-state_data-segmentation/results/x_multiintersect/multiintersect_BcellPC_E9_chr14_106000000_106400000_n=15_processed.bed"

TABLE = read.table(INPUT) %>%
  select(V1, V2, V3) %>%
  mutate(V4 = 1)

write_tsv(TABLE, path = paste0(OUTPUT), col_names = F)
