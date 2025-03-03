#!/bin/bash

#NOTE: When you run this code, run it with "trimmed" or "untrimmed" as an argument to specify if your reads are trimmed or untrimmed.
#alternatively, you can write if your reads are trimmed or untrimmed here by replacing "$1" with "trimmed" or "untrimmed"
READ_FORMAT=$1

#Identify the list of your samples -- untrimmed or trimmed
SAMPLES=/home/riley.risteen/Cb_A549_inf_rnaseq/list_of_samples_${READ_FORMAT}

#Specify where you want the fastqc output to go
PROJECT_DIR=/local/projects-t3/EFUNG/riley.risteen/Cb_A549_inf


#do not change anything below here!                                        

FASTQ_TYPE=original_fastq
FASTQC_DIR=/usr/local/packages/fastqc-0.11.9
THREADS=8

#Check that you correctly input the file name and read format for the list of samples
if [[ ${READ_FORMAT} != untrimmed && ${READ_FORMAT} != trimmed ]]; then
	echo "please specify if your sequences are <trimmed> or <untrimmed>"
	exit
else
	if [ -f $SAMPLES ];then
		:
	else
		FILENAME=$(echo ${SAMPLES} | rev | cut -d/ -f1 | rev)
		echo "${FILENAME} does not exist"
		exit
	fi
fi


echo "$SAMPLES"

#Use the sample list to retrieve each fastq file
for SAMPLE in $(cat "$SAMPLES"); do

        SAMPLE_PATH=$(echo ${SAMPLE} | cut -d',' -f1)
        SAMPLE_NAME=$(echo ${SAMPLE} | cut -d',' -f2)
        OUTPUT_DIR="$PROJECT_DIR"/fastqc/"$READ_FORMAT"_fastqc/"$SAMPLE_NAME"
	
	#Make a directory to put the fastqc results in. The -p argument means that if the directory already exists then it wont make a new one.
        mkdir -p "$OUTPUT_DIR"

	#Define the file paths for a given sample's fastq files
        FASTQ1=${SAMPLE_PATH}_R1.fastq.gz
        FASTQ2=${SAMPLE_PATH}_R2.fastq.gz

	#Check that those fastq files exist
		if [[ -f ${FASTQ1} && -f ${FASTQ2} ]];then :
		else echo "at least one fastq file not found for ${SAMPLE}"
		fi

	#Run fastQC on those fastq files 
    echo -e " \
    "$FASTQC_DIR"/fastqc -o "$OUTPUT_DIR" -f fastq "$FASTQ1" -extract " | \
        qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread "$THREADS" -q threaded.q -l mem_free=30G -wd "$OUTPUT_DIR" -N fastqc.R1."$SAMPLE_NAME"

    echo -e "\
    "$FASTQC_DIR"/fastqc -o "$OUTPUT_DIR" -f fastq "$FASTQ2" -extract " | \
        qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread "$THREADS" -q threaded.q -l mem_free=30G -wd "$OUTPUT_DIR" -N fastqc.R2."$SAMPLE_NAME"

done
