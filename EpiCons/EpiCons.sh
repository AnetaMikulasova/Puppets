#!/bin/bash

#results folder
DIR="/media/aneta/Data/projects/IGH_project/blueprint_ChIPseq_chromatin-state/BD-SE_genome-wide/"

#folder with chromatin states bed files
BEDFILES_DIR="/media/aneta/Data/projects_raw_files/IGH_project/blueprint_ChIPseq_chromatin-state/results/IGV_bedfiles/"

SETTING="./source/setting.txt"
dos2unix ${SETTING}
CELLTYPES_TO_RUN="./source/celltypes.txt"
dos2unix ${CELLTYPES_TO_RUN}

SCRIPTS="./source/"

cat ${SETTING} | awk -F"\t" '$1=="yes"' | while read -r setting || [[ -n "$setting" ]]; do

	#load setting
	SETTING_ID=$(echo "${setting}" | cut -f 2)
	SETTING_ID_BASE_DIR="${SETTING_ID%_*}"

	BRD_ID=$(echo "${setting}" | cut -f 3)
	BRD_CHR_STATE=$(echo "${setting}" | cut -f 4)
	BRD_GAP_SKIP=$(echo "${setting}" | cut -f 5)
	BRD_MINSIZE=$(echo "${setting}" | cut -f 6)
	BRD_MAXSIZE=$(echo "${setting}" | cut -f 7)

	ENH_ID=$(echo "${setting}" | cut -f 8)
	ENH_CHR_STATE=$(echo "${setting}" | cut -f 9)
	ENH_GAP_SKIP=$(echo "${setting}" | cut -f 10)
	ENH_MINSIZE=$(echo "${setting}" | cut -f 11)
	ENH_MAXSIZE=$(echo "${setting}" | cut -f 12)

	BCGR_ID=$(echo "${setting}" | cut -f 13)
	BCGR_CHR_STATE=$(echo "${setting}" | cut -f 14)
	BCGR_GAP_SKIP=$(echo "${setting}" | cut -f 15)
	BCGR_MINSIZE=$(echo "${setting}" | cut -f 16)
	BCGR_MAXSIZE=$(echo "${setting}" | cut -f 17)

	PROM_ID=$(echo "${setting}" | cut -f 18)
	PROM_CHR_STATE=$(echo "${setting}" | cut -f 19)
	PROM_GAP_SKIP=$(echo "${setting}" | cut -f 20)
	PROM_MINSIZE=$(echo "${setting}" | cut -f 21)
	PROM_MAXSIZE=$(echo "${setting}" | cut -f 22)

	DISTANCES=$(echo "${setting}" | cut -f 23)


	#other variables
	RESDIR=${DIR}${SETTING_ID_BASE_DIR}"/results_"${SETTING_ID}"/"
	WORKDIR=${RESDIR}"working_dir/"
	RESULTS_MERGED=${RESDIR}"BD_ENH_calls.txt"

	mkdir -p ${WORKDIR}

	cat ${CELLTYPES_TO_RUN} | awk -F"\t" '$1=="yes"' | while read -r region || [[ -n "$region" ]]; do

		CELLTYPE=$(echo "${region}" | cut -f 2)
		CELLTYPE_ID=$(echo "${region}" | cut -f 3)
		FRACTION=$(echo "${region}" | cut -f 4)

		for ELEMENT in ${BRD_ID} ${ENH_ID} ${BCGR_ID} ${PROM_ID}; do
		# for ELEMENT in ${PROM_ID}; do
		# for ELEMENT in "ENH"; do

			[ ${ELEMENT} == ${BRD_ID} ] && CHR_STATE=${BRD_CHR_STATE} && GAP_SKIP=${BRD_GAP_SKIP} && MINSIZE=${BRD_MINSIZE} && MAXSIZE=${BRD_MAXSIZE}
			[ ${ELEMENT} == ${ENH_ID} ] && CHR_STATE=${ENH_CHR_STATE} && GAP_SKIP=${ENH_GAP_SKIP} && MINSIZE=${ENH_MINSIZE} && MAXSIZE=${ENH_MAXSIZE}
			[ ${ELEMENT} == ${BCGR_ID} ] && CHR_STATE=${BCGR_CHR_STATE} && GAP_SKIP=${BCGR_GAP_SKIP} && MINSIZE=${BCGR_MINSIZE} && MAXSIZE=${BCGR_MAXSIZE}
			[ ${ELEMENT} == ${PROM_ID} ] && CHR_STATE=${PROM_CHR_STATE} && GAP_SKIP=${PROM_GAP_SKIP} && MINSIZE=${PROM_MINSIZE} && MAXSIZE=${PROM_MAXSIZE}
			
			[ ${MAXSIZE} == "max" ] && MAXSIZE=250000000

			WORKDIR_ELEMENT=${WORKDIR}${ELEMENT}"/"
			mkdir -p ${WORKDIR_ELEMENT}"temp/"

			# clean the temp
			rm ${WORKDIR_ELEMENT}"temp/"*

			# subset each file for chromatin states of interest
			Rscript --vanilla ${SCRIPTS}subset_bed_files.R \
				${BEDFILES_DIR} \
				${CELLTYPE} \
				${ELEMENT} \
				${CHR_STATE} \
				${WORKDIR_ELEMENT}"temp/"

			# multiintersect; merge files together and get freq of each fragment
			Nx=`ls ${WORKDIR_ELEMENT}"temp/"* | wc -l`
			N=$(echo $Nx | tr -d ' ')
			OUTFILE=${WORKDIR_ELEMENT}${CELLTYPE_ID}"_"${ELEMENT}"_"${CHR_STATE}"_n="$N
			list=`ls ${WORKDIR_ELEMENT}"temp/"*`
			
			if [ ${N} != 1 ]; then
				multiIntersectBed -i ${list} > ${OUTFILE}"_temp1.bed"
			fi

			#because EC and LC only one sample; multiintersect can't work with one file
			if [ ${N} == 1 ]; then
			Rscript --vanilla ${SCRIPTS}when_one_bed.R \
			${list} \
			${OUTFILE}"_temp1.bed"
			fi

			# clean the temp
			rm ${WORKDIR_ELEMENT}"temp/"*

			# filter fragments for those that are present in minimal number of samples defined by $FRACTION
			Rscript --vanilla ${SCRIPTS}process_subset.R \
			${OUTFILE}"_temp1.bed" \
			${FRACTION} \
			${N} \
			${OUTFILE}"_temp2.bed"

			# sort and delete gaps defined by $GAP_SKIP
			bedtools sort -i ${OUTFILE}"_temp2.bed" > ${OUTFILE}"_temp3.bed"
			bedtools merge -i ${OUTFILE}"_temp3.bed" -d ${GAP_SKIP} > ${OUTFILE}"_temp4.bed"

			# filter fragments for minimal size defined by $MINSIZE
			if [ ${MINSIZE} != "perc" ]; then
				Rscript --vanilla ${SCRIPTS}size_filter.R \
				${OUTFILE}"_temp4.bed" \
				${MINSIZE} \
				${MAXSIZE} \
				${OUTFILE}"_final.bed"
			fi

			# or, filter fragments for top 5% of biggest
			if [ ${MINSIZE} == "perc" ]; then 
				Rscript --vanilla ${SCRIPTS}top_5perc_filter.R \
				${OUTFILE}"_temp4.bed" \
				${OUTFILE}"_final.bed"
			fi

			rm ${WORKDIR_ELEMENT}*_temp*

		done

		# find pairs SE-BrD
		BRD_FILE=${WORKDIR}${BRD_ID}"/"${CELLTYPE_ID}"_"${BRD_ID}"_"${BRD_CHR_STATE}"_n="*"final.bed"
		ENH_FILE=${WORKDIR}${ENH_ID}"/"${CELLTYPE_ID}"_"${ENH_ID}"_"${ENH_CHR_STATE}"_n="*"final.bed"
		BCGR_FILE=${WORKDIR}${BCGR_ID}"/"${CELLTYPE_ID}"_"${BCGR_ID}"_"${BCGR_CHR_STATE}"_n="*"final.bed"
		PROM_FILE=${WORKDIR}${PROM_ID}"/"${CELLTYPE_ID}"_"${PROM_ID}"_"${PROM_CHR_STATE}"_n="*"final.bed"

		Rscript --vanilla ${SCRIPTS}merge_calls.R \
			${BRD_FILE} \
			${ENH_FILE} \
			${BCGR_FILE} \
			${PROM_FILE} \
			${N}"_"${BRD_ID}"-"${FRACTION}"-"${BRD_CHR_STATE}"-"${BRD_GAP_SKIP}"-"${BRD_MINSIZE}"_"${ENH_ID}"-"${FRACTION}"-"${ENH_CHR_STATE}"-"${ENH_GAP_SKIP}"-"${ENH_MINSIZE} \
			${CELLTYPE} \
			${CELLTYPE_ID} \
			${RESULTS_MERGED}

	done

done

