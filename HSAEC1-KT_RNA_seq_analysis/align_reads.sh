#!/bin/bash

HISAT2_DIR=/usr/local/packages/hisat2-2.2.1
SAMTOOLS_BIN_DIR=/usr/local/packages/samtools-1.9/bin
BAM_DIR=/local/projects-t3/EFUNG/riley.risteen/Cb_HSAEC1-KT_inf/bam
INDEXED_REFERENCE_DIR=/local/projects-t3/EFUNG/riley.risteen/indexed_reference
SAMPLES=/home/riley.risteen/Cb_hsaec_inf_rnaseq/list_of_samples_trimmed
THREADS=8

mkdir -p ${BAM_DIR} 


for SAMPLE in $(cat ${SAMPLES}); do
    
    FASTQ_NAME=$(echo ${SAMPLE} | cut -d',' -f1)
	SAMPLE_NAME=$(echo ${SAMPLE} | cut -d',' -f2)


	FASTQ1=${FASTQ_NAME}_R1.fastq.gz
	FASTQ2=${FASTQ_NAME}_R2.fastq.gz

	BAM=${BAM_DIR}/${SAMPLE_NAME}.sorted.bam

        ## Need to split this as annotation pipeline needs the sorted/indexed bam file but htseq needs the fully annotated GFF3 file

    echo -e " "$HISAT2_DIR"/hisat2 -p ${THREADS} --rna-strandness RF --max-intronlen 50000 -x ${INDEXED_REFERENCE_DIR}/GRCh38_latest_genomic -1 "$FASTQ1" -2 "$FASTQ2" | \
            "$SAMTOOLS_BIN_DIR"/samtools sort -@ ${THREADS} -o ${BAM} && \
            "$SAMTOOLS_BIN_DIR"/samtools index -@ ${THREADS} ${BAM}  
    " | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread "$THREADS" -q threaded.q -l mem_free=70G -wd "$BAM_DIR" -N hisat2.align."$SAMPLE_NAME"


done



