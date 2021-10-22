#!/bin/bash

#VARIABLES
#reference prep
PICARD="/home/aneta/Applications/picard_2.21.1/picard.jar"
BWA="/home/aneta/Applications/bwa-0.7.17/bwa"
FASTA_ORIG="/media/aneta/Data/databases/UCSC/fasta/hg38.fa"
FASTA_22XY="/media/aneta/Data/databases/UCSC/fasta/hg38_1-22XY.fa"
FASTA_22X="/media/aneta/Data/databases/UCSC/fasta/hg38_1-22X.fa"
FASTA_22="/media/aneta/Data/databases/UCSC/fasta/hg38_1-22.fa"
#simuG
SIMUG="/home/aneta/Applications/simuG/simuG.pl"
OUTDIR_REF="/media/aneta/Data/projects/ClinVar_in_silico/ref_edit/"
mkdir -p ${OUTDIR_REF}
#ART
ART="/home/aneta/Applications/artbinmountrainier2016.06.05linux64/art_bin_MountRainier/art_illumina"
OUTDIR_FASTQ="/media/aneta/Data/projects_raw_files/ClinVar_in_silico/fastq/raw/"
mkdir -p ${OUTDIR_FASTQ}
MASTER_IN_SILICO_FASTQ="./source/master_in_silico_fastq.txt"
dos2unix ${MASTER_IN_SILICO_FASTQ}



# 1) prepare reference (needed for simuG and ART)
#----------------------------------------------------------------------------------------------------

#functions
function filter_reference {
     
     #$1: list of contigs with separator "_"
     #$2: input fasta
     #$3: output fasta

     contigs=${1//_/ }
     for contig in $contigs ; do
          samtools faidx $2 $contig  >> $3
     done

     #dict index
     REF_BASE="${3%.*}"
     REF_DICT=${REF_BASE}".dict"
     java -jar ${PICARD} CreateSequenceDictionary \
     R=$2 \
     O=${REF_DICT}
     #fai index
     samtools faidx $3
     #bwa indes
     ${BWA} index -a bwtsw $3
}

#a) hg38 -> hg38_1-22XY
chr="chr1_chr2_chr3_chr4_chr5_chr6_chr7_chr8_chr9_chr10_chr11_chr12_chr13_chr14_chr15_chr16_chr17_chr18_chr19_chr20_chr21_chr22_chrX_chrY"
filter_reference $chr $FASTA_ORIG $FASTA_22XY
#b) hg38 -> hg38_1-22X
chr="chr1_chr2_chr3_chr4_chr5_chr6_chr7_chr8_chr9_chr10_chr11_chr12_chr13_chr14_chr15_chr16_chr17_chr18_chr19_chr20_chr21_chr22_chrX"
filter_reference $chr $FASTA_ORIG $FASTA_22X
#c) hg38 -> hg38_1-22
chr="chr1_chr2_chr3_chr4_chr5_chr6_chr7_chr8_chr9_chr10_chr11_chr12_chr13_chr14_chr15_chr16_chr17_chr18_chr19_chr20_chr21_chr22"
filter_reference $chr $FASTA_ORIG $FASTA_22

#----------------------------------------------------------------------------------------------------



#2) prepare aberrant fastq with given abnormalities, see https://github.com/yjx1217/simuG for more examples and VCF formats
#----------------------------------------------------------------------------------------------------

for EACH in "./source/clinvar_selection_A.vcf" "./source/clinvar_selection_B.vcf"; do

     PART=$(basename ${EACH})
     PART=`echo "${PART//.vcf/}"`
     PART=`echo "${PART//clinvar_selection_/}"`
     echo $PART

     perl ${SIMUG} \
          -refseq ${FASTA_22XY} \
          -cnv_vcf ${EACH} \
          -prefix ${OUTDIR_REF}"hg38_"${PART}

done

#----------------------------------------------------------------------------------------------------



# 3) generate fastq files
#----------------------------------------------------------------------------------------------------

# HS20 = HiSeq 2000 (100bp)
# -l = read length
# -f = depth
# -p = paired
# -m = the mean size of DNA/RNA fragments for paired-end simulations
# -s = the standard deviation of DNA/RNA fragment size for paired-end simulations

cat ${MASTER_IN_SILICO_FASTQ} | awk -F"\t" '$1=="yes"' | while read -r line || [[ -n "$line" ]]; do

     ID=$(echo "${line}" | cut -f 2)
     LANE=$(echo "${line}" | cut -f 3)
     REF_FOR_FASTQ=$(echo "${line}" | cut -f 4)

     echo ${ID}"_"${LANE}" with "$REF_FOR_FASTQ
     $ART -ss HS20 -i ${REF_FOR_FASTQ} -o ${OUTDIR_FASTQ}${ID}"_"${LANE}"_R" -l 100 -f 15 -p -m 500 -s 10
     #remove aln files
     rm ${OUTDIR_FASTQ}${ID}"_"${LANE}*aln
     #compress
     gzip ${OUTDIR_FASTQ}${ID}"_"${LANE}"_R"*".fq"

done 

#----------------------------------------------------------------------------------------------------


