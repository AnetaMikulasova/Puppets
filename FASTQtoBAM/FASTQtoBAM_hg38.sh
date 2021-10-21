#!/bin/bash

#VARIABLES TO BE REVIEWED
#----------------------------------------------------------------------------------------------------
#folders
RAW_FASTQ_DIR="/media/aneta/Data/projects_raw_files/ClinVar_in_silico/fastq/raw/"
RESULTS_DIR="/media/aneta/Data/projects_raw_files/ClinVar_in_silico/"
FASTQ_DIR=${RESULTS_DIR}"fastq/"
FASTQC_DIR=${RESULTS_DIR}"qc/fastqc/"
MULTIQC_DIR=${RESULTS_DIR}"qc/multiqc/"
BAM_DIR=${RESULTS_DIR}"/bam/"
BAM_BWA_DIR=${RESULTS_DIR}"/bam/temp/bwa/"
BAM_SORTED_DIR=${RESULTS_DIR}"/bam/temp/sorted/"
BAM_DUPMARK_DIR=${RESULTS_DIR}"/bam/temp/markdup/"
BAM_GATK_DIR=${RESULTS_DIR}"/bam/temp/bqsr_indelrealign/"

#folder for Picard temp files
TEMP_DIR=${RESULTS_DIR}"/temp/"
mkdir -p $TEMP_DIR


#master table
MASTER_SAMPLE="./source/master_preprocess.txt"
dos2unix ${MASTER_SAMPLE}

#reference, databases and softwares
REF="/media/aneta/Data/databases/UCSC/fasta/hg38_1-22XY.fa"
PICARD="/home/aneta/Applications/picard_2.21.1/picard.jar"
GATK="/home/aneta/Applications/gatk-4.1.4.0/gatk"
# https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk
GATK3="/home/aneta/Applications/gatk-3.8-1-0/GenomeAnalysisTK.jar"
FASTQC="/home/aneta/Applications/fastqc_v0.11.8/FastQC/fastqc"
BWA="/home/aneta/Applications/bwa-0.7.17/bwa"
#change following based on genome version, download from ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/
GATK_BUNDLE_DIR="/media/aneta/Data/databases/GATK/gatk-boundle/hg38/UCSC_hg38_1-22XY/"
MILLS=${GATK_BUNDLE_DIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf
KNOWNINDELS=${GATK_BUNDLE_DIR}/Homo_sapiens_assembly38.known_indels.vcf
DBSNP=${GATK_BUNDLE_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf
#----------------------------------------------------------------------------------------------------


#index reference
#----------------------------------------------------------------------------------------------------
#1) for bwa mem
if ! ls ${REF}".bwt" 1> /dev/null 2>&1; then
	${BWA} index -a bwtsw ${REF}
fi
#2) for GATK
REF_BASE="${REF%.*}"
REF_DICT=${REF_BASE}".dict"
if ! ls ${REF_DICT} 1> /dev/null 2>&1; then
	java -jar ${PICARD} CreateSequenceDictionary \
		REFERENCE=${REF} \
		OUTPUT=${REF_DICT}
fi
#3) for GATK
if ! ls ${REF}".fai" 1> /dev/null 2>&1; then
	samtools faidx ${REF}
fi
#----------------------------------------------------------------------------------------------------



#loop master table
#----------------------------------------------------------------------------------------------------
cat ${MASTER_SAMPLE} | awk -F"\t" '$1=="yes"' | while read -r sample || [[ -n "$sample" ]]; do
	SAMPLE_ID=$(echo "${sample}" | cut -f 2)
	BATCH_ID=$(echo "${sample}" | cut -f 4)
	PLATFORM=$(echo "${sample}" | cut -f 5)
	CAPTURE=$(echo "${sample}" | cut -f 6)
	SEQUENCER=$(echo "${sample}" | cut -f 7)

	FASTQ_R1=${FASTQ_DIR}${SAMPLE_ID}"_R1.fastq.gz"
	FASTQ_R2=${FASTQ_DIR}${SAMPLE_ID}"_R2.fastq.gz"
	BAM_BWA=${BAM_BWA_DIR}${SAMPLE_ID}"_bwa.bam"
	BAM_SORTED=${BAM_SORTED_DIR}${SAMPLE_ID}"_sorted.bam"
	BAM_DUPMARK=${BAM_DUPMARK_DIR}${SAMPLE_ID}"_dupmarked.bam"
	BAM_DUPMARK_METRICS=${BAM_DUPMARK_DIR}${SAMPLE_ID}"_dupmarked_metrics.txt"
	BAM_REALIGN_TARGETS=${BAM_GATK_DIR}${SAMPLE_ID}"_realignment_targets.list"
	BAM_REALIGN=${BAM_GATK_DIR}${SAMPLE_ID}"_realigned.bam"
	BAM_BQSR_RECALTABLE=${BAM_GATK_DIR}${SAMPLE_ID}"_recal_data.table"
	BAM_FINAL=${BAM_DIR}${SAMPLE_ID}"_final.bam"
	BAM_RMDUP=${BAM_DIR}${SAMPLE_ID}"_final_remdup.bam"

	#merge fastq files and run fastqc
    mkdir -p ${FASTQ_DIR}
    mkdir -p ${FASTQC_DIR}
    for R in "R1" "R2"; do
	    #prepare, compress and run fastqc
	    if ! ls ${FASTQ_DIR}${SAMPLE_ID}"_"${R}".fastq.gz" 1> /dev/null 2>&1; then
		    echo ${SAMPLE_ID}" - "${R}" - merging fastq files:"
		    touch ${FASTQ_DIR}${SAMPLE_ID}"_"${R}".fastq.gz"
		    for rawfile in ${RAW_FASTQ_DIR}${SAMPLE_ID}*${R}*; do
		        echo ${rawfile}
		        cat ${rawfile} >> ${FASTQ_DIR}${SAMPLE_ID}"_"${R}".fastq.gz"
		    done
		   	echo ${SAMPLE_ID}" - "${R}" - fastqc"
		    ${FASTQC} ${FASTQ_DIR}${SAMPLE_ID}"_"${R}".fastq.gz" -o ${FASTQC_DIR}
		fi
	done

	#align to human genome reference (BWA MEM)
	mkdir -p ${BAM_BWA_DIR}
	if ! ls ${BAM_FINAL} 1> /dev/null 2>&1; then
	if ! ls ${BAM_BWA} 1> /dev/null 2>&1; then
		echo ${SAMPLE_ID}" - mapping to human genome reference"
		RG='@RG\tID:'${BATCH_ID}'\tSM:'${SAMPLE_ID}'\tPL:${PLATFORM}\tLB:'${CAPTURE}'\tPU:'${SEQUENCER}
		${BWA} mem \
	        -t 25 \
	        -M \
	        -R ${RG} \
	        ${REF} \
	        ${FASTQ_R1} \
	        ${FASTQ_R2} | samtools view -bS - > ${BAM_BWA}
	fi
	fi

	#sort alignments (PICARD)
	mkdir -p ${BAM_SORTED_DIR}
	if ! ls ${BAM_FINAL} 1> /dev/null 2>&1; then
	if ! ls ${BAM_SORTED} 1> /dev/null 2>&1; then
		echo ${SAMPLE_ID}" - sorting bam"
		java -Xmx20g -Djava.io.tmpdir=${TEMP_DIR} -jar ${PICARD} SortSam \
			INPUT=${BAM_BWA} \
			OUTPUT=${BAM_SORTED} \
			SORT_ORDER=coordinate \
			TMP_DIR=${TEMP_DIR}
	fi
	fi

	#mark duplicates (PICARD)
	mkdir -p ${BAM_DUPMARK_DIR}
	if ! ls ${BAM_FINAL} 1> /dev/null 2>&1; then
	if ! ls ${BAM_DUPMARK} 1> /dev/null 2>&1; then
		echo ${SAMPLE_ID}" - marking duplicates in bam"
		java -Xmx20g -Djava.io.tmpdir=${TEMP_DIR} -jar ${PICARD} MarkDuplicates \
			INPUT=${BAM_SORTED} \
			OUTPUT=${BAM_DUPMARK} \
			METRICS_FILE=${BAM_DUPMARK_METRICS} \
			TMP_DIR=${TEMP_DIR}
		#index bam (PICARD)
		java -Xmx20g -Djava.io.tmpdir=${TEMP_DIR} -jar ${PICARD} BuildBamIndex \
			INPUT=${BAM_DUPMARK} \
			TMP_DIR=${TEMP_DIR}
	fi
	fi

	#indel realignment
	mkdir -p ${BAM_GATK_DIR}
	if ! ls ${BAM_FINAL} 1> /dev/null 2>&1; then
	if ! ls ${BAM_REALIGN} 1> /dev/null 2>&1; then
		#find targets
		echo ${SAMPLE_ID}" - finding targets for indel realignment"
		java -Xmx25g -jar \
		    ${GATK3} \
		        -T RealignerTargetCreator \
		        -R ${REF} \
		        -known ${MILLS} \
		        -known ${KNOWNINDELS} \
		        -I ${BAM_DUPMARK} \
		        -o ${BAM_REALIGN_TARGETS}
		#apply realignment
		echo ${SAMPLE_ID}" - indel realignment application"
		java -Xmx25g -jar \
		    ${GATK3} \
		        -T IndelRealigner \
		        -R ${REF} \
		        -known ${MILLS} \
		        -known ${KNOWNINDELS} \
		        -LOD 0.4 \
		        --maxReadsForRealignment 10000000 \
		        --maxConsensuses 300 \
		        --maxReadsForConsensuses 1200 \
		        -targetIntervals ${BAM_REALIGN_TARGETS} \
		        -I ${BAM_DUPMARK} \
		        -o ${BAM_REALIGN}
	fi
	fi

	#base quality scores recalibration (BQSR, GATK)
	if ! ls ${BAM_FINAL} 1> /dev/null 2>&1; then
		#build BQSR model
		echo ${SAMPLE_ID}" - bqsr model"
	    ${GATK} BaseRecalibrator \
	    	--java-options "-Xmx25g" \
	        -R ${REF} \
	        --known-sites ${MILLS} \
	        --known-sites ${KNOWNINDELS} \
	        --known-sites ${DBSNP} \
	        -I ${BAM_REALIGN} \
	        -O ${BAM_BQSR_RECALTABLE}
		#apply BQSR
		echo ${SAMPLE_ID}" - bqsr application"
	    ${GATK} ApplyBQSR \
	    	--java-options "-Xmx25g" \
	        -R ${REF} \
	        -bqsr ${BAM_BQSR_RECALTABLE} \
	        -I ${BAM_REALIGN} \
	        -O ${BAM_FINAL}
	fi

	# #remove duplicates (SAMTOOLS)
	# if ! ls ${BAM_RMDUP} 1> /dev/null 2>&1; then
	# 	echo ${SAMPLE_ID}" - removing duplicates in bam"
	# 	samtools view -b -F 0x0400 ${BAM_FINAL} > ${BAM_RMDUP}
	# 	samtools index ${BAM_RMDUP}	
	# fi

	rm -r ${BAM_DIR}"/temp"
done
#----------------------------------------------------------------------------------------------------

rm -r ${TEMP_DIR}

#multiqc
#----------------------------------------------------------------------------------------------------
if ! ls ${MULTIQC_DIR}multiqc_report.html 1> /dev/null 2>&1; then
	echo "running multiqc..."
	mkdir -p ${MULTIQC_DIR}
	multiqc --outdir ${MULTIQC_DIR} ${FASTQC_DIR}
fi
#----------------------------------------------------------------------------------------------------





