args <- commandArgs(trailing = TRUE)

library(dplyr)
library(tidyverse)
library(data.table)

INDIR       = args[1] %>% as.character()
CELLTYPE    = args[2] %>% as.character()
ELEMENT     = args[3] %>% as.character()
CHRSTATE    = args[4] %>% strsplit("_", fixed=T) %>% as_vector() %>% as.list()
OUTDIR      = args[5] %>% as.character()

# INDIR = "/Volumes/Aneta_Data/projects/IGH_project/blueprint_ChIPseq_chromatin-state/results/IGV_bedfiles/"
# ELEMENT = "BrD"
# CHRSTATE = "E10_E11" %>% strsplit("_", fixed=T) %>% as_vector() %>% as.list()
# OUTDIR = "/Volumes/Aneta_Data/projects/IGH_project/blueprint_ChIPseq_chromatin-state/results/working_bedfiles_BD-SE_genome-wide/BrD/"

#get data
files = list.files(INDIR, recursive = T, full.names = T, pattern = CELLTYPE)

#read file, subset for chrstate, position and export files
# sample = paste0(INDIR, "H_plasma-cell_T_S00Y8QH1_11.bed")
for (sample in files){
  NAME          <- basename(sample) %>% gsub(".bed", paste0("_", ELEMENT, ".bed"),.)
  FILE          <- read.table(sample)
  names(FILE)   <- c("CONTIG", "START", "END", "CHR_STATE", "X", "Y", "START2", "END2", "CHR_STATE_COLOR")
  FILE          <- subset(FILE, CHR_STATE %in% CHRSTATE)
  write_tsv(FILE, path = paste0(OUTDIR, NAME), col_names = F)
}
