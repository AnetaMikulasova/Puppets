#!/bin/bash

#variables
BASE="/path/to/outpu/folder/"
CTRLS=${BASE}"controls/vcf/"
mkdir $CTRLS
PON=${BASE}"controls/pon.vcf.gz"

GATK="/path/to/gatk-4.0.6.0/"
REF="/path/to/gatk-boundle/b37/human_g1k_v37.fasta"
CAPTURE="/path/to/Capture.bed"
GNOMAD1="/path/to/gatk-boundle/mutect2/af-only-gnomad.raw.sites.b37.vcf.gz"
GNOMAD2="/path/to/gatk-boundle/mutect2/small_exac_common_3_b37.vcf.gz"

SAMPLES="./samples.txt"
dos2unix $SAMPLES
CONTROLS="/controls.txt"
dos2unix $CONTROLS

#These databases were in GATK3 - PoN and germline resource replaced cosmic input
#COSMIC="/Users/anetamikulasova/Documents/applications/gatk-boundle/cosmic/GRCh37_v86/COSMIC_b37_v86.vcf"
#DBSNP="/Users/anetamikulasova/Documents/applications/gatk-boundle/b37/dbsnp_138.b37.vcf"

#PoN does not exist in GATK3
#Contamination filter does not exist in GAKT3



#==============[1]=PREPARE-CONTROLS========================================

#apply Mutect2 to controls
cat ${CONTROLS} | awk -F"\t" '$1=="yes"' | while read -r control || [[ -n "$control" ]]; do

	PATID=$(echo "${samples}" | cut -f 2)
	CTRLID=$(echo "${samples}" | cut -f 3)
	CTRLBAM=$(echo "${samples}" | cut -f 4)

	VCF_OUT=${CTRLS}${PAT_ID}"_"${CTRLID}".vcf.gz"

	if ! ls ${VCF_OUT} 1> /dev/null 2>&1; then
	${GATK}gatk Mutect2 \
		-R ${REF} \
		-I ${CTRLBAM} \
		-tumor ${PATID}"_"${CTRLID} \
		-L ${CAPTURE} \
		-ip 76 \
		-O ${VCF_OUT}
	fi

done

#create PON
vcflist=`ls ${CTRLS}*.vcf.gz | sed 's/^/-vcfs /'`
${GATK}gatk CreateSomaticPanelOfNormals \
	${vcflist} \
	-O ${PON}
#============================================================================



#Process samples by master table to get IDs

cat ${SAMPLES} | awk -F"\t" '$1=="yes"' | while read -r samples || [[ -n "$samples" ]]; do

	PATIENT=$(echo "${samples}" | cut -f 2)
	TUMOR=$(echo "${samples}" | cut -f 3)
	TUMOR_BAM=$(echo "${samples}" | cut -f 4)
	GERMLINE=$(echo "${samples}" | cut -f 5)
	GERMLINE_BAM=$(echo "${samples}" | cut -f 6)


	echo "Processing *** patient ${PATIENT} *** tumor ${TUMOR} *** germline ${GERMLINE}"

	OUTDIR=${BASE}"Results_snv/"${PATIENT}"/"
	mkdir ${OUTDIR}
	OUTPATTERN=${OUTDIR}${PATIENT}"_"${TUMOR}
	OUTPATTERN_GERM=${OUTDIR}${PATIENT}"_"${GERMLINE}


#==============[1]=RUN-MUTECT2-ON-TUMOR-SAMPLE===============================

	#Mutect2-to-get-raw-set
	${GATK}gatk Mutect2 \
		-R ${REF} \
		-I ${TUMOR_BAM} \
		-I ${GERMLINE_BAM} \
		-tumor ${TUMsOR} \
		-normal ${GERMLINE} \
		-pon ${PON} \
		--germline-resource ${GNOMAD1} \
		-L ${CAPTURE} \
		-ip 76 \
		-O ${OUTPATTERN}"_somatic.vcf.gz" \
		-bamout ${OUTPATTERN}".bam"

#============================================================================



#==============[2]=PREPARE-CONTAMINTAION-FOR-FILTERING=======================

	#note -L and -ip were added to speed it up

	#Get-summaries-for-tumor-sample
	${GATK}gatk GetPileupSummaries \
		-I ${TUMOR_BAM} \
		-V ${GNOMAD2} \
		-L ${CAPTURE} \
		-ip 76 \
		-O ${OUTPATTERN}"_getpileupsummaries.table"

	#Get-summaries-for-germline-sample
	${GATK}gatk GetPileupSummaries \
		-I ${GERMLINE_BAM} \
		-V ${GNOMAD2} \
		-L ${CAPTURE} \
		-ip 76 \
		-O ${OUTPATTERN_GERM}"_getpileupsummaries.table"

	#Calculate-contamination-in-tumor-sample
	${GATK}gatk CalculateContamination \
		-I ${OUTPATTERN}"_getpileupsummaries.table" \
		-matched ${OUTPATTERN_GERM}"_getpileupsummaries.table" \
		-O ${OUTPATTERN}"_calculatecontamination.table" \
		-segments ${OUTPATTERN}"_tumorsegments.table"

#============================================================================



#==============[3]=FILTER-RAWSET-BY-STANDARD-FILTERS=========================

	#FILTERING-14-filters-including-contamination
	${GATK}gatk FilterMutectCalls \
		-V ${OUTPATTERN}"_somatic.vcf.gz" \
		--contamination-table ${OUTPATTERN}"_calculatecontamination.table" \
		-O ${OUTPATTERN}"_somatic_filtered_level1.vcf.gz"

#============================================================================



#==============[4]=ADDITIONAL-FILTERS-FOR-SEQUENCING-ARTIFACT================

	#Collect-sequencing-artifact-metrics
	${GATK}gatk CollectSequencingArtifactMetrics \
		-R ${REF} \
		-I ${TUMOR_BAM} \
		-O ${OUTPATTERN}"_artifact" \
		--FILE_EXTENSION ".txt"

	#Filter-by-orientation-bias; G->T transversion = OxoG, C->T transition = FFPE 
	${GATK}gatk FilterByOrientationBias \
		-AM G/T \
		-AM C/T \
		-V ${OUTPATTERN}"_somatic_filtered_level1.vcf.gz" \
		-P ${OUTPATTERN}"_artifact.pre_adapter_detail_metrics.txt" \
		-O ${OUTPATTERN}"_somatic_filtered_level2.vcf.gz"

#============================================================================

done