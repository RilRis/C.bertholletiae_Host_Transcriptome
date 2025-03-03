#!/bin/bash

HISAT2_DIR=/usr/local/packages/hisat2-2.2.1
SAMTOOLS_BIN_DIR=/usr/local/packages/samtools-1.9/bin
BAM_DIR=/local/projects-t3/EFUNG/riley.risteen/Cb_A549_inf/bam
INDEXED_REFERENCE_DIR=/local/projects-t3/EFUNG/riley.risteen/indexed_reference
SAMPLES=/home/riley.risteen/Cb_A549_inf_rnaseq/list_of_samples_trimmed
THREADS=8


#first, make sure the specified files exist
if [[ ! -d ${INDEXED_REFERENCE_DIR} || ! -f ${SAMPLES} ]]; then
		echo "indexed reference genome and/or sample list not found"
		echo "aborting process"
		exit
fi



#now do stuff
mkdir -p ${BAM_DIR} 


for SAMPLE in $(cat ${SAMPLES}); do
    
    FASTQ_NAME=$(echo ${SAMPLE} | cut -d',' -f1)
	SAMPLE_NAME=$(echo ${SAMPLE} | cut -d',' -f2)


	FASTQ1=${FASTQ_NAME}_R1.fastq.gz
	FASTQ2=${FASTQ_NAME}_R2.fastq.gz

	BAM=${BAM_DIR}/${SAMPLE_NAME}.sorted.bam

        ## Must keep HISAT2 alignment and HTSeq quantitation seperate during annotation pipeline because HISAT2 needs the sorted/indexed bam file but HTSeq needs the fully annotated GFF3 file

    echo -e " "$HISAT2_DIR"/hisat2 -p ${THREADS} --rna-strandness RF --max-intronlen 50000 -x ${INDEXED_REFERENCE_DIR}/GRCh38_latest_genomic -1 "$FASTQ1" -2 "$FASTQ2" | \
            "$SAMTOOLS_BIN_DIR"/samtools sort -@ ${THREADS} -o ${BAM} && \
            "$SAMTOOLS_BIN_DIR"/samtools index -@ ${THREADS} ${BAM}  
    " | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread "$THREADS" -q threaded.q -l mem_free=70G -wd "$BAM_DIR" -N hisat2.align."$SAMPLE_NAME"


done

#NOTE: each sample will have its own output and error file
	#output file name:	hisat2.align.[$SAMPLE_NAME].o[Job ID]
	#error file name: 	hisat2.align.[$SAMPLE_NAME].e[Job ID]
		#(Job ID is an arbitrary number assigned by the qsub command)

