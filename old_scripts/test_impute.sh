#!/bin/bash

# Test the imputation workflow
INPUT=/media/data/Daniel/Scripts/cell.counts.txt
IMPUTEDIR="/media/data/Daniel/test_data/imputation_output/"
IMPOUTPUT=$IMPUTEDIR/imputed.$(basename $INPUT)
THREADS=4
IMPUTE="yes"
CELLCOUNTS="/media/data/Daniel/Scripts/cell.counts.txt"

echo "transform file into csv file"
cat $INPUT | tr "\\t" "," | tr " " "," > ./cell.counts.temp.csv
head ./cell.counts.temp.csv

if [ ${IMPUTE} == "yes" ]; then
    if [ -z ${IMPUTOUT+x} ]; then
        DIR=$(dirname "${CELLCOUNTS}")
        BASE=$(basename "${CELLCOUNTS}")
        IMPOUTPUT=$DIR/imputed.$BASE
    fi
    echo "Perform imputation on raw count table"
    Rscript ./impute_counts.R "./cell.counts.temp.csv" $IMPOUTPUT $IMPUTEDIR $THREADS
fi
rm ./cell.counts.temp.csv

R
