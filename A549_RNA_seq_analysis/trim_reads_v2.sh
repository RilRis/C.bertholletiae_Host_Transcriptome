#!/bin/bash
#modified from code provided by Chris Holt of the Dunning Hotopp Lab
#this code has two main functions: 
	#1. use trimmer to trim the ends of the RNA seq reads
	#2. create a list of the file paths for all of the trimmed reads (which is needed for later steps)

##specify where you want the trimmed sequences to be stored
RNA_TRIMMED=/local/projects-t3/EFUNG/riley.risteen/Cb_A549_inf/trimmed_sequences
THREADS=8

##list of fastq.gz files (with full filepaths for each) along with sample name in a CSV. The filenames have the ends removed from them (_R1.fastq.gz)
	##EXAMPLE: The filename is EFUNG_N5UD-B01_L002_R1.fastq.gz and it is located in /local/projects-t3/EFUNG/UI_3h_B_RNA/ILLUMINA_DATA/
	##/local/projects-t3/EFUNG/UI_3h_B_RNA/ILLUMINA_DATA/EFUNG_N5UD-B01_L002,Cb175_UI_3h_B
SAMPLES=/home/riley.risteen/Cb_A549_inf_rnaseq/list_of_samples_untrimmed



#ALL VARIABLES ARE ABOVE THIS LINE
#################################################################################

#the following paths are specific to trimmomatic on our servers. DO NOT CHANGE
JDK_DIR=/usr/local/packages/jdk/bin/
TRIMMOMATIC_DIR=/usr/local/packages/trimmomatic-0.38/

#create a file in which your filepaths for the trimmed reads will be recorded (take sample list filepath and replace "untrimmed" with "trimmed")
touch ${SAMPLES/_untrimmed/_trimmed}

#create a directory in which your trimmed reads for all samples will be stored
mkdir -p "$RNA_TRIMMED"

for SAMPLE in $(cat "$SAMPLES"); do

        FASTQ_NAME=$(echo ${SAMPLE} | cut -d',' -f1)
        SAMPLE_NAME=$(echo ${SAMPLE} | cut -d',' -f2)

        FASTQ1=${FASTQ_NAME}_R1.fastq.gz
        FASTQ2=${FASTQ_NAME}_R2.fastq.gz
	
	OUTPUT_DIR="$RNA_TRIMMED"/"$SAMPLE_NAME"
	
	#adds the sample to the list of trimmed read filepaths
	echo "${OUTPUT_DIR}/${SAMPLE_NAME}_paired,$SAMPLE_NAME" >> ${SAMPLES/_untrimmed/_trimmed}

	#create a subdirectory in which your trimmed reads for a particular sample will be stored
    mkdir -p "$OUTPUT_DIR"
	
		
##The following lines are just to make sure your variables are correct. Comment out when not troubleshooting.
#echo "$FASTQ1" "$FASTQ2" "$SAMPLE_NAME"	
#exit


    ## TruSeq3-PE-2.fa downloaded from https://github.com/timflutre/trimmomatic/blob/master/adapters/TruSeq3-PE-2.fa

     echo -e " \
         "$JDK_DIR"/java -Xmx20g -jar "$TRIMMOMATIC_DIR"/trimmomatic-0.38.jar PE -phred33 -trimlog "$OUTPUT_DIR"/trim.log \
         "$FASTQ1" "$FASTQ2" "$OUTPUT_DIR"/"$SAMPLE_NAME"_paired_R1.fastq.gz "$OUTPUT_DIR"/"$SAMPLE_NAME"_unpaired_R1.fastq.gz "$OUTPUT_DIR"/"$SAMPLE_NAME"_paired_R2.fastq.gz  \
         "$OUTPUT_DIR"/"$SAMPLE_NAME"_unpaired_R2.fastq.gz \
         ILLUMINACLIP:"$TRIMMOMATIC_DIR"/adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 \
         " | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread "$THREADS" -q threaded.q -l mem_free=30G -wd "$OUTPUT_DIR" -N trimmomatic."$SAMPLE_NAME"

done

#remove the duplicate lines from your sample list
sort -u -o ${SAMPLES/_untrimmed/_trimmed}{,}
