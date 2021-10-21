args <- commandArgs(trailing = TRUE)

library(dplyr)
library(tidyverse)
library(data.table)

INPUT       = args[1] %>% as.character()
OUTPUT      = args[2] %>% as.character()

# INPUT = "/Users/aneta/Documents/projects/IGH_project/blueprint_ChIPseq_chromatin-state/BD-SE_genome-wide/results_setting2/working_dir/BrD/AML_BrD_E10_E11_n=34_final.bed" %>% as.character()
# OUTPUT = "/Users/aneta/Documents/projects/IGH_project/blueprint_ChIPseq_chromatin-state/BD-SE_genome-wide/results_setting2/working_dir/BrD/AML_BrD_E10_E11_n=34_final.bed_temp" %>% as.character()

TABLE = read.table(INPUT) %>% mutate(V4 = V3 - V2) %>% arrange(V4)

NofTOT = nrow(TABLE)
Perc = 0.05 * NofTOT
Perc = round(Perc)

TABLE_FILTER = tail(TABLE,Perc)
MINSIZE = min(TABLE_FILTER$V4)

TABLE_FILTER2 = TABLE %>% filter(V4 > MINSIZE) %>%
  arrange(V1, V2, V3)

write_tsv(TABLE_FILTER2, path = paste0(OUTPUT), col_names = F)
