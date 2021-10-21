#!/bin/bash

FILES_SHEET="./files_sheet.txt"
dos2unix ${FILES_SHEET}
PICARD="/home/aneta/Applications/picard_2.21.1/picard.jar"

BAM_DIR="/path/to/folder/with/bam/"
FASTQ_OUT_DIR="/path/to/folder/for/fastq/output/"
mkdir -p ${FASTQ_OUT_DIR}

#dealign and split lanes for each input bam file
#--------------------------------------------------------------------
cat ${FILES_SHEET} | awk -F"\t" '$1=="yes"' | while read -r sample || [[ -n "$sample" ]]; do
    
    SAMPLE_ID=$(echo "${sample}" | cut -f 2)
    BAM_ID=$(echo "${sample}" | cut -f 3)

    BAM_IN=${BAM_DIR}${BAM_ID}
    FASTQ_R1_OUT=${FASTQ_OUT_DIR}${SAMPLE_ID}"_R1.fastq"
    FASTQ_R2_OUT=${FASTQ_OUT_DIR}${SAMPLE_ID}"_R2.fastq"
    FASTQ_SPLIT_BEG=${FASTQ_OUT_DIR}${SAMPLE_ID}

    #dealign
    java -Xmx2g -jar ${PICARD} SamToFastq I=${BAM_IN} F=${FASTQ_R1_OUT} F2=${FASTQ_R2_OUT}

    #split lanes
    #R1
    awk -v FASTQ_SPLIT_BEG="$FASTQ_SPLIT_BEG" 'BEGIN {FS = ":"} {lane=$4 ; print > FASTQ_SPLIT_BEG"_L00"lane"_R1.fastq" ; for (i = 1; i <= 3; i++) {getline ; print > FASTQ_SPLIT_BEG"_L00"lane"_R1.fastq"}}' < ${FASTQ_R1_OUT}
    rm ${FASTQ_R1_OUT}  
    #R2
    awk -v FASTQ_SPLIT_BEG="$FASTQ_SPLIT_BEG" 'BEGIN {FS = ":"} {lane=$4 ; print > FASTQ_SPLIT_BEG"_L00"lane"_R2.fastq" ; for (i = 1; i <= 3; i++) {getline ; print > FASTQ_SPLIT_BEG"_L00"lane"_R2.fastq"}}' < ${FASTQ_R2_OUT}
    rm ${FASTQ_R2_OUT}

    gzip ${FASTQ_OUT_DIR}${SAMPLE_ID}*".fastq"

    ##compress
    ## gzip ${FASTQ_R1_OUT}
    ## gzip ${FASTQ_R2_OUT}
    ## gzip -c ${FASTQ_R1_OUT} > ${FASTQ_R1_OUT}".gz"
    ## gzip -c ${FASTQ_R2_OUT} > ${FASTQ_R2_OUT}".gz"

done
#--------------------------------------------------------------------

