#!/bin/bash

#notes:
# pipeline done by Bluprint http://dcc.blueprint-epigenome.eu/#/md/chip_seq_grch38
# PhantomPeakQualTools: not working on Mac (some awk problem); follow the instruction by https://github.com/crazyhottommy/phantompeakqualtools
# PhantomPeakQualTools: then download the R script "run_spp.R" + I made changes as runmean function was not recognized (replaced by caTools::runmean())
# MACS2: follow instruction INSTALL.md in downloaded dir https://pypi.org/project/MACS2/#files; 
# MACS2: create Python Virtual Environmnet (python3 -m venv ~/Applications/PythonVirtEnv) and activate (source ~/Applications/PythonVirEnv/bin/activate)
# MACS2: install numpy (pip install numpy) and MACS2 (pip install MACS2) within the activated virtual environment


#variables
SCRIPTS="./source/"

MASTER_SAMPLE="./master_TALL_ChIPseq_raw_files.txt"
dos2unix $MASTER_SAMPLE
FASTQ_DIR_BASE="/media/aneta/Data/projects/IGH_project/T-ALL_ChIPseq/from_raw/raw_files/"
FASTQC="/home/aneta/Applications/fastqc_v0.11.8/FastQC/fastqc"
BWA="/home/aneta/Applications/bwa-0.7.17/bwa"
REF="/media/aneta/Data/databases/Blueprint/Homo_sapiens.GRCh37.chromosomes.fa"

WORKING_DIR_BASE="/media/aneta/Data/projects/IGH_project/T-ALL_ChIPseq/from_raw/working_dir/"

PICARD="/home/aneta/Applications/picard_2.21.1/picard.jar"
PHANTOM="/home/aneta/Applications/PhantomPeakQualTools/run_spp_edit.R"
UCSC_CHROMOSOME_SIZE="./hg19.chrom.sizes.txt"



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


#download SRA files and convert to FASTQ
#----------------------------------------------------------------------------------------------------
#note internal variables in following script, not automatic
Rscript --vanilla ${SCRIPTS}"download_sra_and_get_fastq.R"
#----------------------------------------------------------------------------------------------------


#merge fastq files
#----------------------------------------------------------------------------------------------------
cat ${MASTER_SAMPLE} | awk -F"\t" '$1=="yes"' | while read -r sample || [[ -n "$sample" ]]; do
# cat ${MASTER_SAMPLE} | awk -F"\t" '$1!="sample_id"' | while read -r sample || [[ -n "$sample" ]]; do

	SAMPLE_ID=$(echo "${sample}" | cut -f 2)
	HISTONE_MARK=$(echo "${sample}" | cut -f 3)
	RUN=$(echo "${sample}" | cut -f 7)
	SRP_ID=$(echo "${sample}" | cut -f 9)

	FASTQ_DIR=${FASTQ_DIR_BASE}${SRP_ID}"/"
	FASTQ_FILE_IN=${FASTQ_DIR}${RUN}".1.fastq"

	WORKING_DIR=${WORKING_DIR_BASE}${SRP_ID}"/"
	FASTQC_DIR=${WORKING_DIR}"qc/fastqc/"

	TEMP_DIR=${WORKING_DIR}"temp/"${SAMPLE_ID}"/"
	mkdir -p ${TEMP_DIR}
	FASTQ_FILE_OUT=${TEMP_DIR}${SAMPLE_ID}"_"${HISTONE_MARK}".fastq"
	FASTQ_FILE_OUT_GZ=${TEMP_DIR}${SAMPLE_ID}"_"${HISTONE_MARK}".fastq.gz"

	if ! ls ${FASTQC_DIR}${SAMPLE_ID}"_"${HISTONE_MARK}* 1> /dev/null 2>&1; then

		#merge fastq files
		if ! ls ${FASTQ_FILE_OUT} 1> /dev/null 2>&1; then
	    	touch ${FASTQ_FILE_OUT}
		fi
		echo -e "\033[1;31m$(date): "${FASTQ_FILE_IN}"   --->   "${FASTQ_FILE_OUT}"\033[0m"
		cat ${FASTQ_FILE_IN} >> ${FASTQ_FILE_OUT}
	fi
done
#----------------------------------------------------------------------------------------------------	

#compress fastq run fastqc
#----------------------------------------------------------------------------------------------------
cat ${MASTER_SAMPLE} | awk -F"\t" '$1=="yes"' | while read -r sample || [[ -n "$sample" ]]; do
# cat ${MASTER_SAMPLE} | awk -F"\t" '$1!="sample_id"' | while read -r sample || [[ -n "$sample" ]]; do

	SAMPLE_ID=$(echo "${sample}" | cut -f 2)
	HISTONE_MARK=$(echo "${sample}" | cut -f 3)
	RUN=$(echo "${sample}" | cut -f 7)
	SRP_ID=$(echo "${sample}" | cut -f 9)

	# FASTQ_DIR=${FASTQ_DIR_BASE}${SRP_ID}"/"
	# FASTQ_FILE_IN=${FASTQ_DIR}${RUN}".1.fastq"

	WORKING_DIR=${WORKING_DIR_BASE}${SRP_ID}"/"
	FASTQC_DIR=${WORKING_DIR}"qc/fastqc/"

	TEMP_DIR=${WORKING_DIR}"temp/"${SAMPLE_ID}"/"
	mkdir -p ${TEMP_DIR}
	FASTQ_FILE_OUT=${TEMP_DIR}${SAMPLE_ID}"_"${HISTONE_MARK}".fastq"
	FASTQ_FILE_OUT_GZ=${TEMP_DIR}${SAMPLE_ID}"_"${HISTONE_MARK}".fastq.gz"

	if ! ls ${FASTQC_DIR}${SAMPLE_ID}"_"${HISTONE_MARK}* 1> /dev/null 2>&1; then

		#compress fastq to fastq.gz
		if ! ls ${FASTQ_FILE_OUT_GZ} 1> /dev/null 2>&1; then
			echo -e "\033[1;31m$(date): "${SAMPLE_ID}" - "${HISTONE_MARK}" - compressing fastq""\033[0m"
			gzip ${FASTQ_FILE_OUT}
		fi

		#run fastqc
		mkdir -p ${FASTQC_DIR}
		echo -e "\033[1;31m$(date): "${SAMPLE_ID}" - "${HISTONE_MARK}" - running fastq""\033[0m"
		${FASTQC} ${FASTQ_FILE_OUT_GZ} -o ${FASTQC_DIR}
	fi

done
#----------------------------------------------------------------------------------------------------


#run multiqc on fastqc
#----------------------------------------------------------------------------------------------------
cat ${MASTER_SAMPLE} | awk -F"\t" '$1=="yes"' | while read -r sample || [[ -n "$sample" ]]; do

	SRP_ID=$(echo "${sample}" | cut -f 9)

	WORKING_DIR=${WORKING_DIR_BASE}${SRP_ID}"/"
	FASTQC_DIR=${WORKING_DIR}"qc/fastqc/"
	MULTIQC_DIR=${WORKING_DIR}"qc/multiqc/"

	if ! ls ${MULTIQC_DIR}multiqc_report.html 1> /dev/null 2>&1; then
		echo -e "\033[1;31m$(date): running multiqc...""\033[0m"
		mkdir -p ${MULTIQC_DIR}
		multiqc --outdir ${MULTIQC_DIR} ${FASTQC_DIR}
	fi

done
#----------------------------------------------------------------------------------------------------




#align, process bams and get fragment size
#----------------------------------------------------------------------------------------------------
cat ${MASTER_SAMPLE} | awk -F"\t" '$1=="yes"' | while read -r sample || [[ -n "$sample" ]]; do

	SAMPLE_ID=$(echo "${sample}" | cut -f 2)
	HISTONE_MARK=$(echo "${sample}" | cut -f 3)
	SRP_ID=$(echo "${sample}" | cut -f 9)
	PLATFORM=$(echo "${sample}" | cut -f 10)
	SEQUENCER=$(echo "${sample}" | cut -f 11)

	WORKING_DIR=${WORKING_DIR_BASE}${SRP_ID}"/"

	#temp work
	TEMP_DIR=${WORKING_DIR}"temp/"${SAMPLE_ID}"/"
	OUTPUT_PATTERN_TEMP=${TEMP_DIR}${SAMPLE_ID}"_"${HISTONE_MARK}

	FASTQ_INPUT=${OUTPUT_PATTERN_TEMP}.fastq.gz


	SAIFILE=${OUTPUT_PATTERN_TEMP}".sai"

	BAMFILE_ALIGN=${OUTPUT_PATTERN_TEMP}"_align.bam"
	BAMFILE_SORTED=${OUTPUT_PATTERN_TEMP}"_sorted.bam"
	BAMFILE_MARKED=${OUTPUT_PATTERN_TEMP}"_dupmarked.bam"
	METRICS_MARKED=${OUTPUT_PATTERN_TEMP}"_metrics.txt"
	BAMFILE_FILTER=${OUTPUT_PATTERN_TEMP}"_filter.bam"

	#final files
	mkdir -p ${WORKING_DIR}"BAM/"
	mkdir -p ${WORKING_DIR}"PhantomPeakQualTools/"
	BAMFILE_RMDUP=${WORKING_DIR}"BAM/"${SAMPLE_ID}"_"${HISTONE_MARK}"_final.bam"
	PAROUT=${WORKING_DIR}"PhantomPeakQualTools/"${SAMPLE_ID}"_"${HISTONE_MARK}"_params.out"


	#align fastq - single end sequencing
	if ! ls ${BAMFILE_ALIGN} 1> /dev/null 2>&1; then
		echo -e "\033[1;31m$(date): "${SAMPLE_ID}" - "${HISTONE_MARK}" - bwa alignment""\033[0m"
		RG='@RG\tID:'${SRP_ID}'\tSM:'${SAMPLE_ID}'_'${HISTONE_MARK}'\tPL:'${PLATFORM}'\tPU:'${SEQUENCER}

		${BWA} aln ${REF} ${FASTQ_INPUT} > ${SAIFILE} ; ${BWA} samse -r ${RG} ${REF} ${SAIFILE} ${FASTQ_INPUT} | samtools view -bS - > ${BAMFILE_ALIGN}

	fi

	# #align fastq - pair-end read sequencing; prededing code needs to be edited
	# if ! ls ${BAMFILE_ALIGN} 1> /dev/null 2>&1; then
	# 	echo -e "\033[1;31m$(date): "${SAMPLE_ID}" - "${HISTONE_MARK}" - bwa alignment""\033[0m"
	# 	RG='@RG\tID:'${SAMPLE_ID}'\tSM:'${SAMPLE_ID}'_'${HISTONE_MARK}'\tPL:'${PLATFORM}'\tPU:'${SEQUENCER}		
	# 	if ! ls ${SAIFILE_1} 1> /dev/null 2>&1; then
	# 		${BWA} aln ${REF} ${FASTQ_R1} > ${SAIFILE_1}
	# 	fi
	# 	if ! ls ${SAIFILE_2} 1> /dev/null 2>&1; then
	# 		${BWA} aln ${REF} ${FASTQ_R2} > ${SAIFILE_2}
	# 	fi
	# 	${BWA} sampe -r ${RG} ${REF} ${SAIFILE_1} ${SAIFILE_2} ${FASTQ_R1} ${FASTQ_R2} | samtools view -bS - > ${BAMFILE_ALIGN}
	# fi

	#sort bam file
	if ! ls ${BAMFILE_SORTED} 1> /dev/null 2>&1; then
		echo -e "\033[1;31m$(date): "${SAMPLE_ID}" - "${HISTONE_MARK}" - sorting bam file""\033[0m"
		java -Xmx20g -jar ${PICARD} \
							SortSam \
							INPUT=${BAMFILE_ALIGN} \
							OUTPUT=${BAMFILE_SORTED} \
							SORT_ORDER=coordinate \
							VALIDATION_STRINGENCY=SILENT
	
	fi

	#mark bam file
	if ! ls ${BAMFILE_MARKED} 1> /dev/null 2>&1; then
		echo -e "\033[1;31m$(date): "${SAMPLE_ID}" - "${HISTONE_MARK}" - marking duplicates in bam""\033[0m"
		java -Xmx20g -jar ${PICARD} \
							MarkDuplicates \
							INPUT=${BAMFILE_SORTED} \
							OUTPUT=${BAMFILE_MARKED} \
							METRICS_FILE=${METRICS_MARKED} \
							REMOVE_DUPLICATES=false \
							ASSUME_SORTED=true \
							VALIDATION_STRINGENCY=SILENT
	fi



	# filter to remove unmapped reads and reads with Mapping Quality less than 5
	if ! ls ${BAMFILE_FILTER} 1> /dev/null 2>&1; then
		echo -e "\033[1;31m$(date): "${SAMPLE_ID}" - "${HISTONE_MARK}" - filtering bam file""\033[0m"
		samtools view -b -F 4 -q 5 ${BAMFILE_MARKED} > ${BAMFILE_FILTER}
	fi

	# remove duplicate reads
	if ! ls ${BAMFILE_RMDUP} 1> /dev/null 2>&1; then
		echo -e "\033[1;31m$(date): "${SAMPLE_ID}" - "${HISTONE_MARK}" - removind duplicates in bam""\033[0m"
		samtools view -b -F 1024 ${BAMFILE_FILTER} > ${BAMFILE_RMDUP}
		samtools index ${BAMFILE_RMDUP}
	fi

	# # # remove intermed files
	# # rm ${BAMFILE_SORTED}
	# # rm ${BAMFILE_MARKED}
	# # rm ${METRICS_MARKED}
	# # rm ${BAMFILE_FILTER}

	# model fragment size
	if ! ls ${PAROUT} 1> /dev/null 2>&1; then
		echo -e "\033[1;31m$(date): "${SAMPLE_ID}" - "${HISTONE_MARK}" - modeling fragment size""\033[0m"
		Rscript ${PHANTOM} -c=${BAMFILE_RMDUP} -rf -out=${PAROUT}
	fi

done
#----------------------------------------------------------------------------------------------------




#peak calling
#----------------------------------------------------------------------------------------------------
cat ${MASTER_SAMPLE} | awk -F"\t" '$1=="yes"' | awk -F"\t" '$3!="input"' | while read -r sample || [[ -n "$sample" ]]; do

	SAMPLE_ID=$(echo "${sample}" | cut -f 2)
	HISTONE_MARK=$(echo "${sample}" | cut -f 3)
	SRP_ID=$(echo "${sample}" | cut -f 9)
	INPUT_ID=$(echo "${sample}" | cut -f 12)
	# INPUT_ID="WholeCellExtracts"

	WORKING_DIR=${WORKING_DIR_BASE}${SRP_ID}"/"

	SAMPLE_BAM=${WORKING_DIR}"BAM/"${SAMPLE_ID}"_"${HISTONE_MARK}"_final.bam"
	INPUT_BAM=${WORKING_DIR}"BAM/"${SAMPLE_ID}"_"${INPUT_ID}"_final.bam"
	PHANTOM_FILE=${WORKING_DIR}"PhantomPeakQualTools/"${SAMPLE_ID}"_"${HISTONE_MARK}"_params.out"

	OUT_DIR=${WORKING_DIR}"MACS2/"
	mkdir -p ${OUT_DIR}

	source /home/aneta/Applications/PythonVirtEnv/bin/activate

	#get half fragments size
	FRAG_SIZE=`sed -r 's/,[^\t]+//g' ${PHANTOM_FILE} | awk -F"\t" '{print $3}'`
	HALF_FRAG_SIZE=`echo $((FRAG_SIZE/2))`

	if [ ${HISTONE_MARK} != ${INPUT_ID} ]; then

		if ! ls ${OUT_DIR}${SAMPLE_ID}"_"${HISTONE_MARK}"_peaks.xls" 1> /dev/null 2>&1; then

			echo -e "\033[1;31m$(date): "${SAMPLE_ID}" - "${HISTONE_MARK}" - running MACS2 peak calling (half-frag-size = "${HALF_FRAG_SIZE}")""\033[0m"
			
			if [ ${HISTONE_MARK} == "H3K27me3" ] || [ ${HISTONE_MARK} == "H3K36me3" ] || [ ${HISTONE_MARK} == "H3K9me3" ] || [ ${HISTONE_MARK} == "H3K4me1" ]; then
				macs2 callpeak \
					-t ${SAMPLE_BAM} \
					-c ${INPUT_BAM} \
					-n ${SAMPLE_ID}"_"${HISTONE_MARK} \
					--outdir ${OUT_DIR} \
					--gsize hs \
					--nomodel \
					--extsize ${HALF_FRAG_SIZE} \
					--bdg \
					--broad
			fi
			
			if [ ${HISTONE_MARK} == "H3K27ac" ] || [ ${HISTONE_MARK} == "H3K4me3" ] || [ ${HISTONE_MARK} == "H3K9/14ac" ] || [ ${HISTONE_MARK} == "H2A.Zac" ]; then
				macs2 callpeak \
					-t ${SAMPLE_BAM} \
					-c ${INPUT_BAM} \
					-n ${SAMPLE_ID}"_"${HISTONE_MARK} \
					--outdir ${OUT_DIR} \
					--gsize hs \
					--nomodel \
					--extsize ${HALF_FRAG_SIZE} \
					--bdg
			fi
		fi

		#convert bedGraph to BigWig
		if ! ls ${OUT_DIR}${SAMPLE_ID}"_"${HISTONE_MARK}"_control_lambda.bw" 1> /dev/null 2>&1; then
			echo -e "\033[1;31m$(date): "${SAMPLE_ID}" - "${HISTONE_MARK}" - converting control wig""\033[0m"
			bedGraphToBigWig ${OUT_DIR}${SAMPLE_ID}"_"${HISTONE_MARK}"_control_lambda.bdg" ${UCSC_CHROMOSOME_SIZE} ${OUT_DIR}${SAMPLE_ID}"_"${HISTONE_MARK}"_control_lambda.bw"
		fi

		if ! ls ${OUT_DIR}${SAMPLE_ID}"_"${HISTONE_MARK}"_treat_pileup.bw" 1> /dev/null 2>&1; then
			echo -e "\033[1;31m$(date): "${SAMPLE_ID}" - "${HISTONE_MARK}" - converting sample wig""\033[0m"
			bedGraphToBigWig ${OUT_DIR}${SAMPLE_ID}"_"${HISTONE_MARK}"_treat_pileup.bdg" ${UCSC_CHROMOSOME_SIZE} ${OUT_DIR}${SAMPLE_ID}"_"${HISTONE_MARK}"_treat_pileup.bw"
		fi
	
	fi

done
#----------------------------------------------------------------------------------------------------































