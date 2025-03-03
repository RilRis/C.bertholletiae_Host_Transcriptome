# Analyzing the A549 Data

## Generating the count table
###1) Download and prepare reference files if you have not already done so
#### Retrieve and unzip the reference genome 
```
cd /local/projects-t3/EFUNG/riley.risteen
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz
gunzip /local/projects-t3/EFUNG/riley.risteen/GRCh38_latest_genomic.fna.gz
```

#### Index the reference genome
Run **build_index.sh**
```
#!/bin/bash

HISAT2_DIR=/usr/local/packages/hisat2-2.2.1
REFERENCE_GENOME=/local/projects-t3/EFUNG/riley.risteen/GRCh38_latest_genomic.fna
INDEXED_REFERENCE_DIR=/local/projects-t3/EFUNG/riley.risteen/indexed_reference
THREADS=4

mkdir -p ${SAM_DIR} ${INDEXED_REFERENCE_DIR} 

echo -e "
${HISAT2_DIR}/hisat2-build ${REFERENCE_GENOME} ${INDEXED_REFERENCE_DIR}/GRCh38_latest_genomic
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread "$THREADS" -q threaded.q -l mem_free=20G -wd ${INDEXED_REFERENCE_DIR} -N hisat2.build
```

#### Retrieve and unzip the annotation file associated with the reference genome
```
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gtf.gz
gunzip /local/projects-t3/EFUNG/riley.risteen/GRCh38_latest_genomic.gtf.gz
```

###2) Create a list of file paths to sequences and name it list_of_samples_untrimmed
Each sample will have several sequence files within its folder. Identify just the part of the file names that are shared by all files for that sample. For each sample, include the absolute file path and just the part of the file name that is common between all seq files for that sample. Then, write a comma and add a sample name for that sequence file. 
![image](https://github.com/user-attachments/assets/02bd0e05-8dba-4686-b995-0c72d8823fa4)



###3) Use trimmomatic to trim the reads for each sample
Run **trim_reads_v2.sh**

Note: this will automatically create a list of of file paths to your trimmed reads
```
#!/bin/bash
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
	
     echo -e " \
         "$JDK_DIR"/java -Xmx20g -jar "$TRIMMOMATIC_DIR"/trimmomatic-0.38.jar PE -phred33 -trimlog "$OUTPUT_DIR"/trim.log \
         "$FASTQ1" "$FASTQ2" "$OUTPUT_DIR"/"$SAMPLE_NAME"_paired_R1.fastq.gz "$OUTPUT_DIR"/"$SAMPLE_NAME"_unpaired_R1.fastq.gz "$OUTPUT_DIR"/"$SAMPLE_NAME"_paired_R2.fastq.gz  \
         "$OUTPUT_DIR"/"$SAMPLE_NAME"_unpaired_R2.fastq.gz \
         ILLUMINACLIP:"$TRIMMOMATIC_DIR"/adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 \
         " | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread "$THREADS" -q threaded.q -l mem_free=30G -wd "$OUTPUT_DIR" -N trimmomatic."$SAMPLE_NAME"

done

#remove the duplicate lines from your sample list
sort -u -o ${SAMPLES/_untrimmed/_trimmed}{,}
```

###4) Use FastQC to assess trimmed reads
Run **fastqc_reads_v2.sh** with the "trimmed" argument
```
#!/bin/bash

#NOTE: When you run this code, run it with "trimmed" or "untrimmed" as an argument to specify if your reads are trimmed or untrimmed.
#alternatively, you can write if your reads are trimmed or untrimmed here by replacing "$1" with "trimmed" or "untrimmed"
READ_FORMAT=$1

#Identify the list of your samples -- untrimmed or trimmed
SAMPLES=/home/riley.risteen/Cb_A549_inf_rnaseq/list_of_samples_${READ_FORMAT}

#Specify where you want the fastqc output to go
PROJECT_DIR=/local/projects-t3/EFUNG/riley.risteen/Cb_A549_inf

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
Done
```
###5) Align reads to reference genome using HISAT2
Run **align_reads.sh**
```
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

#now you can use HISAT2 to align the reads
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
```

###6) Use HTSeq to create a counts file for each sample
Run **generate_counts.sh**
```
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
```

###7) Create a single file containing the counts from all samples
Run **generate_ master_count_table.sh**
```
#!/bin/bash

#Define path of the directory containing each sample's count file. This is also where your master count table will be output.
FILE_PATH=/local/projects-t3/EFUNG/riley.risteen/Cb_A549_inf/counts

#IMPORTANT: all of your counts files should end with ".counts" If it doesnt, you will have to write in a a new pattern to match the files you want on line 17

#########################################################################

#list the files you want to work with. In order, the below chunk...
#   1) fetches all the file names in FILE_PATH
#   2) filters to only keep the names that end in ".counts"
#   3) if you have already run this script, excludes any output files to only get the count tables for each sample
#   4) keep only the file names (remove paths)

SAMPLE_LIST=$(find ${FILE_PATH} \
  | grep "\.counts$" \
  | grep -v "master_count_table" \
  | xargs -L 1 basename)   


#check if the output file ("master_count_table.counts") exists. If it does, we will change the name of the new output
if [[ -f ${FILE_PATH}/master_count_table.counts ]]; then
		TIMESTAMP=$(date +%Y%m%d_%H%M%S)
		MASTER_TABLE=NEW_master_count_table_${TIMESTAMP}.counts
		echo "'master_count_table.counts' already exists"
		echo "new master count table will be saved as '${MASTER_TABLE}'"
	else
		MASTER_TABLE=master_count_table.counts
fi


#before we start the loop, create a starting amount (0) so we can add up the number of lines in each file as we go
CALC_TOTAL=0


#for each file in the list of samples we created....
for f in $SAMPLE_LIST; do 
	
	#copy the contents of the count file and add it to a master file
	cat ${FILE_PATH}/${f} >> ${FILE_PATH}/${MASTER_TABLE}
		
	#QUALITY CONTROL
	##count the number of lines in this sample's file
	INDIV_COUNT=$(wc -l ${FILE_PATH}/${f} | cut -d " " -f1)
	
	##add this file's line count to the running total
	CALC_TOTAL=$((${CALC_TOTAL} + ${INDIV_COUNT})) 

	#EXTRA OUTPUT (if you want more details of what is happening while this runs, uncomment the below lines)
	#echo "Processing ${f}..."
	#echo "Total lines: ${CALC_TOTAL}"

done


#now compare the calculated line total with the actual line total
##count the number of lines in the MASTER_TABLE file
COUNT_TOTAL=$(wc -l ${FILE_PATH}/${MASTER_TABLE} | cut -d " " -f1)

##compare, throwing an error if COUNT_TOTAL does not equal CALC_TOTAL
if (($COUNT_TOTAL != $CALC_TOTAL)); then
	echo "ERROR!!! INCORRECT NUMBER OF LINES IN ${MASTER_TABLE}!"
fi
```

## Differential expression analysis

### Create a metadata csv file (one row per sample) with the following details for each sample:
- sample name (with replicate identifier)
- fungal strain (for uninfected samples, write the strain that they acted as a control for)
  * in this case, our uninfected samples acted as controls for both Cb strains so I just left it as "Cb" 
- condition (infected vs uninfected)
- time point
- experimental group

Use the same column names as provided in this example:
| sample       | strain | infec_status | time_point | group        |
|-------------|--------|----------------|------------|-------------|
| Cb175_I_3h_A | Cb175  | I              | 3h         | Cb175_I_3h  |
| Cb175_I_3h_B | Cb175  | I              | 3h         | Cb175_I_3h  |
| Cb175_I_3h_C | Cb175  | I              | 3h         | Cb175_I_3h  |
| Cb175_I_6h_A | Cb175  | I              | 6h         | Cb175_I_6h  |
| Cb175_I_6h_B | Cb175  | I              | 6h         | Cb175_I_6h  |
| Cb175_I_6h_C | Cb175  | I              | 6h         | Cb175_I_6h  |
| Cb182_I_3h_A | Cb182  | I              | 3h         | Cb182_I_3h  |
| Cb182_I_3h_B | Cb182  | I              | 3h         | Cb182_I_3h  |
| Cb182_I_3h_C | Cb182  | I              | 3h         | Cb182_I_3h  |
| Cb182_I_6h_A | Cb182  | I              | 6h         | Cb182_I_6h  |
| Cb182_I_6h_B | Cb182  | I              | 6h         | Cb182_I_6h  |
| Cb182_I_6h_C | Cb182  | I              | 6h         | Cb182_I_6h  |
| Cb_UI_3h_A   | Cb     | UI             | 3h         | Cb_UI_3h    |
| Cb_UI_3h_B   | Cb     | UI             | 3h         | Cb_UI_3h    |
| Cb_UI_3h_C   | Cb     | UI             | 3h         | Cb_UI_3h    |
| Cb_UI_6h_A   | Cb     | UI             | 6h         | Cb_UI_6h    |
| Cb_UI_6h_B   | Cb     | UI             | 6h         | Cb_UI_6h    |
| Cb_UI_6h_C   | Cb     | UI             | 6h         | Cb_UI_6h    |

Make sure the control samples (in this case, the uninfected [UI] samples) are at the bottom of the table to ensure the following differential expression analysis calculates difference in experimental group expression relative to the control group.

### Calculate differential expression in R using edgeR
1) Open **A549_DE_analysis_R.Rmd** and modify the filepaths within the third chunk
2) Run **A549_DE_analysis_R.Rmd**
