# Transcriptional analysis

## Overview
The process outlined here (described by Holt & Dunning Hotopp, 2024) is composed of two major parts:
1) Processing the raw RNA seq reads into a count table
2) Using the count table to carry out differential expression analysis

The code used for each cell line (A549 and HSAEC1-KT) is divided into two seperate folders. Before you start, you will need to download and index the reference human genome. This reference will be used for aligning the reads from both cell lines (for step 1 of the process outlined above).


## Download and prepare reference files
### Retrieve and unzip the reference genome 
```
cd /local/projects-t3/EFUNG/riley.risteen
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz
gunzip /local/projects-t3/EFUNG/riley.risteen/GRCh38_latest_genomic.fna.gz
```

### Index the reference genome
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




## References:
Holt, C. I., & Dunning Hotopp, J. C. (2024). Updated annotation and meta-analysis of Brugia malayi transcriptomics data reveals consistent transcriptional profiles across time and space with some study-specific differences in adult female worm transcriptional profiles. *PLOS Neglected Tropical Diseases*, 18(9), e0012511.
