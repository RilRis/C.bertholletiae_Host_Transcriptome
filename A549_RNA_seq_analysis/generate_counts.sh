#!/bin/bash


PYTHON_BIN_DIR=/usr/local/packages/python-3.8.2/bin
PROJECT_DIR=/local/projects-t3/EFUNG/riley.risteen/Cb_A549_inf
	BAM_DIR=${PROJECT_DIR}/bam
	COUNTS_DIR=${PROJECT_DIR}/counts
SAMPLES=/home/riley.risteen/Cb_A549_inf_rnaseq/list_of_samples_trimmed
THREADS=10
GFF3=/local/projects-t3/EFUNG/riley.risteen/GRCh38_latest_genomic.gff


#first, make sure the specified files/folders exist
if [[ ! -f ${GFF3} || ! -f ${SAMPLES} || ! -d ${BAM_DIR} ]]; then
	echo "at least one of the files/folders needed could not be found"
	echo "aborting process"
	exit
fi

#now generate the counts
mkdir -p ${COUNTS_DIR}

for SAMPLE in $(cat ${SAMPLES}); do
   
    SAMPLE_NAME=$(echo ${SAMPLE} | cut -d',' -f2)       
	BAM=${BAM_DIR}/${SAMPLE_NAME}.sorted.bam
	COUNTS_FILE=${COUNTS_DIR}/${SAMPLE_NAME}.counts

         ## Must keep HISAT2 alignment and HTSeq quantitation seperate during annotation pipeline because HISAT2 needs the sorted/indexed bam file but HTSeq needs the fully annotated GFF3 file
 echo -e " ${PYTHON_BIN_DIR}/htseq-count -n ${THREADS} -s reverse --max-reads-in-buffer 3000000000 -r pos --nonunique none -f bam -m union -t gene --idattr ID ${BAM} ${GFF3} \
        | sed -e 's/ /\t/g' > ${COUNTS_FILE} " | qsub -V -q threaded.q -pe thread "$THREADS" -P jhotopp-gcid-proj4b-filariasis -N htseq_counts."$SAMPLE_NAME" -wd "$COUNTS_DIR" -l mem_free=70G -hold_jid hisat2.align.${SAMPLE_NAME}

done

