#!/bin/bash

FASTQ_TYPE=original_fastq

PROJECT_DIR=/local/projects-t3/EFUNG/riley.risteen/Cb_HSAEC1-KT_inf

FASTQC_DIR=/usr/local/packages/fastqc-0.11.9
THREADS=8

#State if your reads are trimmed or untrimmed. Write "trimmed" or "untrimmed"
READ_FORMAT=trimmed
#Identify the list of your samples -- untrimmed or trimmed
SAMPLES=/home/riley.risteen/Cb_hsaec_inf_rnaseq/list_of_samples_${READ_FORMAT}

##############################################################################################################################

echo "$SAMPLES"
for SAMPLE in $(cat "$SAMPLES"); do

        SAMPLE_PATH=$(echo ${SAMPLE} | cut -d',' -f1)
        SAMPLE_NAME=$(echo ${SAMPLE} | cut -d',' -f2)
        OUTPUT_DIR="$PROJECT_DIR"/fastqc/"$READ_FORMAT"_fastqc/"$SAMPLE_NAME"
	#Makes a directory to put the fastqc results in. The -p argument means that if the directory already exists then it wont make a new one.
        mkdir -p "$OUTPUT_DIR"


        FASTQ1=${SAMPLE_PATH}_R1.fastq.gz
        FASTQ2=${SAMPLE_PATH}_R2.fastq.gz


   echo -e " \
    "$FASTQC_DIR"/fastqc -o "$OUTPUT_DIR" -f fastq "$FASTQ1" -extract " | \
        qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread "$THREADS" -q threaded.q -l mem_free=30G -wd "$OUTPUT_DIR" -N fastqc.R1."$SAMPLE_NAME"

    echo -e "\
    "$FASTQC_DIR"/fastqc -o "$OUTPUT_DIR" -f fastq "$FASTQ2" -extract " | \
        qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread "$THREADS" -q threaded.q -l mem_free=30G -wd "$OUTPUT_DIR" -N fastqc.R2."$SAMPLE_NAME"

done
