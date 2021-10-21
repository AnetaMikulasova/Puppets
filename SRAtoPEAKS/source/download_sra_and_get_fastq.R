library(dplyr)
library(expss)  #vlookup

# http://www.sthda.com/english/wiki/from-sra-to-fastq-file


#----------------------------------------------------------------------------------------------------
#1) KOPT-T1 and DND-41 - SRP035662
#Read SRA file infos
# info<-read.csv("/Users/aneta/Documents/projects/IGH_project/T-ALL_ChIPseq/raw_data/SRP035662/doc/sra_result.csv", stringsAsFactors=FALSE) 
info<-read.csv("/media/aneta/Data/projects/IGH_project/T-ALL_ChIPseq/from_raw/raw_files/SRP035662/doc/sra_result.csv", stringsAsFactors=FALSE) 

# sri<-read.csv("/Users/aneta/Documents/projects/IGH_project/T-ALL_ChIPseq/raw_data/SRP035662/doc/SraRunInfo.csv", stringsAsFactors=FALSE) %>%
sri<-read.csv("/media/aneta/Data/projects/IGH_project/T-ALL_ChIPseq/from_raw/raw_files/SRP035662/doc/SraRunInfo.csv", stringsAsFactors=FALSE) %>%
  mutate(Experiment.Title = vlookup(Experiment, info, lookup_column = 1, result_column = 2)) %>% 
  # filter(grepl("KOPT-K1", Experiment.Title)) %>%
  filter(!grepl("persisters", Experiment.Title)) %>%
  filter(!grepl("BRD4", Experiment.Title))
  # filter(Run == "SRR1143136")
  

files<-basename(sri$download_path)
# setwd("/Users/aneta/Documents/projects/IGH_project/T-ALL_ChIPseq/raw_data/SRP035662/")
setwd("/media/aneta/Data/projects/IGH_project/T-ALL_ChIPseq/from_raw/raw_files/SRP035662/")

for(i in 1:length(files)) download.file(sri$download_path[i], files[i])

stopifnot(all(file.exists(files)))
# f = "SRR1143134.1"
for(f in files) {
  # cmd = paste("/Users/aneta/Documents/applications/sratoolkit.2.9.6-1-mac64/bin/fasterq-dump --split-3", f)
  # cmd = paste("/home/aneta/Applications/sratoolkit.2.9.6-1-ubuntu64/bin/fasterq-dump --split-3", f)
  cmd = paste("/home/aneta/Applications/sratoolkit.2.9.6-1-ubuntu64/bin/fasterq-dump", f) #fasterq-dump is equal to fast-dump --split-3 --skip-technical, see github for fasterq-dump
  cat(cmd,"\n")#print the current command
  system(cmd) # invoke command
}

#----------------------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------------
#2) JURKAT - SRP053266 and SRP103118
#Read SRA file infos

#SRP053266
info<-read.csv("/media/aneta/Data/projects/IGH_project/T-ALL_ChIPseq/from_raw/raw_files/SRP053266/doc/sra_result.csv", stringsAsFactors=FALSE) 
sri<-read.csv("/media/aneta/Data/projects/IGH_project/T-ALL_ChIPseq/from_raw/raw_files/SRP053266/doc/SraRunInfo.csv", stringsAsFactors=FALSE) %>%
  mutate(Experiment.Title = vlookup(Experiment, info, lookup_column = 1, result_column = 2)) %>%
  filter(grepl("GFP_Input|H3K4me1_GFP|H3K4me3_GFP|H3K27Ac_GFP|H3K36me3_GFP", Experiment.Title))
files<-basename(sri$download_path)
setwd("/media/aneta/Data/projects/IGH_project/T-ALL_ChIPseq/from_raw/raw_files/SRP053266/")
for(i in 1:length(files)) download.file(sri$download_path[i], files[i])
stopifnot(all(file.exists(files)))
for(f in files) {
  # cmd = paste("/Users/aneta/Documents/applications/sratoolkit.2.9.6-1-mac64/bin/fasterq-dump --split-3", f)
  # cmd = paste("/home/aneta/Applications/sratoolkit.2.9.6-1-ubuntu64/bin/fasterq-dump --split-3", f)
  cmd = paste("/home/aneta/Applications/sratoolkit.2.9.6-1-ubuntu64/bin/fasterq-dump", f) #fasterq-dump is equal to fast-dump --split-3 --skip-technical, see github for fasterq-dump
  cat(cmd,"\n")#print the current command
  system(cmd) # invoke command
}


#SRP103118
info<-read.csv("/media/aneta/Data/projects/IGH_project/T-ALL_ChIPseq/from_raw/raw_files/SRP103118/doc/sra_result.csv", stringsAsFactors=FALSE) 
sri<-read.csv("/media/aneta/Data/projects/IGH_project/T-ALL_ChIPseq/from_raw/raw_files/SRP103118/doc/SraRunInfo.csv", stringsAsFactors=FALSE) %>%
  mutate(Experiment.Title = vlookup(Experiment, info, lookup_column = 1, result_column = 2)) %>%
  filter(grepl("H3K27me3", Experiment.Title)) %>%
  filter(grepl("DMSO", Experiment.Title))
files<-basename(sri$download_path)
setwd("/media/aneta/Data/projects/IGH_project/T-ALL_ChIPseq/from_raw/raw_files/SRP103118/")
for(i in 1:length(files)) download.file(sri$download_path[i], files[i])
stopifnot(all(file.exists(files)))
for(f in files) {
  # cmd = paste("/Users/aneta/Documents/applications/sratoolkit.2.9.6-1-mac64/bin/fasterq-dump --split-3", f)
  # cmd = paste("/home/aneta/Applications/sratoolkit.2.9.6-1-ubuntu64/bin/fasterq-dump --split-3", f)
  cmd = paste("/home/aneta/Applications/sratoolkit.2.9.6-1-ubuntu64/bin/fasterq-dump", f) #fasterq-dump is equal to fast-dump --split-3 --skip-technical, see github for fasterq-dump
  cat(cmd,"\n")#print the current command
  system(cmd) # invoke command
}

#----------------------------------------------------------------------------------------------------


