args <- commandArgs(trailing = TRUE)

library(tidyverse)

GET_BRD_INPUT = args[1] %>% as.character()
GET_ENH_INPUT = args[2] %>% as.character()
GET_BCGR_INPUT = args[3] %>% as.character()
GET_PROM_INPUT = args[4] %>% as.character()
GET_SETTING = args[5] %>% as.character()
GET_CELLTYPE = args[6] %>% as.character()
GET_CELLTYPE_ID = args[7] %>% as.character()
OUTPUT_MERGED_CALLS = args[8] %>% as.character()

# GET_BRD_INPUT = "/Users/aneta/Documents/projects/IGH_project/blueprint_ChIPseq_chromatin-state/BD-SE_genome-wide/results_setting18/working_dir/BrD/Bcell_BrD_E10_E11_n=15_final.bed"
# GET_ENH_INPUT = "/Users/aneta/Documents/projects/IGH_project/blueprint_ChIPseq_chromatin-state/BD-SE_genome-wide/results_setting18/working_dir/ENH/Bcell_ENH_E9_n=15_final.bed"
# GET_CELLTYPE = "^H_T-cell_CD8"
# GET_CELLTYPE_ID = "TcellCD8"
# OUTPUT_MERGED_CALLS = "/Users/aneta/Documents/projects/IGH_project/blueprint_ChIPseq_chromatin-state/BD-SE_genome-wide/results_setting18/BD_ENH_calls.txt"

# GET_BRD_INPUT="/media/aneta/Data/projects/IGH_project/blueprint_ChIPseq_chromatin-state/BD-SE_genome-wide/results_setting18/working_dir/BrD/Bcell_BrD_E10_E11_n=15_final.bed"
# GET_ENH_INPUT="/media/aneta/Data/projects/IGH_project/blueprint_ChIPseq_chromatin-state/BD-SE_genome-wide/results_setting18/working_dir/ENH/Bcell_ENH_E9_n=15_final.bed"
# GET_CELLTYPE="^H_B-cell_|^H_plasma-cell_"
# GET_CELLTYPE_ID="Bcell"
# OUTPUT_MERGED_CALLS="/media/aneta/Data/projects/IGH_project/blueprint_ChIPseq_chromatin-state/BD-SE_genome-wide/results_setting18/BD_ENH_calls.txt"


#get data
if(file.exists(OUTPUT_MERGED_CALLS)==F) {RESULTS_MERGED_CALLS = data.frame()}
if(file.exists(OUTPUT_MERGED_CALLS)==T) {
  RESULTS_MERGED_CALLS = read.delim(paste0(OUTPUT_MERGED_CALLS), stringsAsFactors = F)
}

# message(GET_CELLTYPE_ID)
#get data
#-------------------------------------------------------------------------------------------
#inputs
BRD_INPUT = read.delim(GET_BRD_INPUT, stringsAsFactors = F, header = F) %>%
  select(V1,V2,V3) %>%
  mutate(ELEMENT = "BRD") %>%
  mutate(SETTING = GET_SETTING) %>%
  mutate(CELLTYPE = GET_CELLTYPE) %>%
  mutate(CELLTYPE_ID = GET_CELLTYPE_ID)
names(BRD_INPUT) = c("CONTIG", "START", "END", "ELEMENT", "SETTING" ,"CELLTYPE", "CELLTYPE_ID")

ENH_INPUT = read.delim(GET_ENH_INPUT, stringsAsFactors = F, header = F) %>%
  select(V1,V2,V3) %>%
  mutate(ELEMENT = "ENH") %>%
  mutate(SETTING = GET_SETTING) %>%
  mutate(CELLTYPE = GET_CELLTYPE) %>%
  mutate(CELLTYPE_ID = GET_CELLTYPE_ID)
names(ENH_INPUT) = c("CONTIG", "START", "END", "ELEMENT", "SETTING" , "CELLTYPE", "CELLTYPE_ID")

BCGR_INPUT = read.delim(GET_BCGR_INPUT, stringsAsFactors = F, header = F) %>%
  select(V1,V2,V3) %>%
  mutate(ELEMENT = "BCGR") %>%
  mutate(SETTING = GET_SETTING) %>%
  mutate(CELLTYPE = GET_CELLTYPE) %>%
  mutate(CELLTYPE_ID = GET_CELLTYPE_ID)
names(BCGR_INPUT) = c("CONTIG", "START", "END", "ELEMENT", "SETTING" , "CELLTYPE", "CELLTYPE_ID")

PROM_INPUT = read.delim(GET_PROM_INPUT, stringsAsFactors = F, header = F) %>%
  select(V1,V2,V3) %>%
  mutate(ELEMENT = "PROM") %>%
  mutate(SETTING = GET_SETTING) %>%
  mutate(CELLTYPE = GET_CELLTYPE) %>%
  mutate(CELLTYPE_ID = GET_CELLTYPE_ID)
names(PROM_INPUT) = c("CONTIG", "START", "END", "ELEMENT", "SETTING" , "CELLTYPE", "CELLTYPE_ID")

#-------------------------------------------------------------------------------------------

#merge to output
#-------------------------------------------------------------------------------------------
RESULTS_MERGED_CALLS = rbind(RESULTS_MERGED_CALLS, BRD_INPUT, ENH_INPUT, BCGR_INPUT, PROM_INPUT)
write_tsv(RESULTS_MERGED_CALLS, path = OUTPUT_MERGED_CALLS, col_names = T)
#-------------------------------------------------------------------------------------------
