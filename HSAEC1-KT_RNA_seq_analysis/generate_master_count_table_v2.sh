#!/bin/bash

#Define path of the directory containing each sample's count file. This is also where your master count table will be output.
FILE_PATH=/local/projects-t3/EFUNG/riley.risteen/Cb_HSAEC1-KT_inf/counts

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