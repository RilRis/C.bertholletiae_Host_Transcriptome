#!/bin/bash


HISAT2_DIR=/usr/local/packages/hisat2-2.2.1
REFERENCE_GENOME=/local/projects-t3/EFUNG/riley.risteen/GRCh38_latest_genomic.fna
INDEXED_REFERENCE_DIR=/local/projects-t3/EFUNG/riley.risteen/indexed_reference
THREADS=4


mkdir -p ${SAM_DIR} ${INDEXED_REFERENCE_DIR} 


echo -e "
${HISAT2_DIR}/hisat2-build ${REFERENCE_GENOME} ${INDEXED_REFERENCE_DIR}/GRCh38_latest_genomic
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread "$THREADS" -q threaded.q -l mem_free=20G -wd ${INDEXED_REFERENCE_DIR} -N hisat2.build



