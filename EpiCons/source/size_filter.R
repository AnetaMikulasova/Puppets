args <- commandArgs(trailing = TRUE)

library(dplyr)
library(tidyverse)
library(data.table)

INPUT       = args[1] %>% as.character()
MINSIZE     = args[2] %>% as.numeric()
MAXSIZE     =  args[3] %>% as.numeric()
OUTPUT      = args[4] %>% as.character()

# INPUT = "/Volumes/Aneta_Data/projects/IGH_project/blueprint_ChIPseq_chromatin-state/results/working_bedfiles_BD-SE_genome-wide/BrD/H_plasma-cell_T_S00Y8QH1_11_BrD_sort_merge.bed" %>% as.character()
# MINSIZE = 15000 %>% as.numeric()
# OUTPUT = "/Volumes/Aneta_Data/projects/IGH_project/blueprint_ChIPseq_chromatin-state/results/working_bedfiles_BD-SE_genome-wide/BrD/H_plasma-cell_T_S00Y8QH1_11_BrD_sort_merge_sizefilter.bed" %>% as.character()

TABLE = read.table(INPUT) %>%
  # mutate(V4 = V3 - (V2+1)) %>% #+1 because of bed file file type; this can't be
  mutate(V4 = V3 - V2) %>%
  filter(V4 >= MINSIZE) %>%
  filter(V4 <= MAXSIZE) %>%
  select(V1, V2, V3, V4)

write_tsv(TABLE, path = paste0(OUTPUT), col_names = F)
