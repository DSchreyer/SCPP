#!/bin/bash
#
# Daniel Schreyer
# Generate count table from RSEM output files
# Input: Directory path with RSEM output files
# Output: Cell count table

while [[ $# -gt 0 ]]; do
  option="$1"
  case $option in
    --data)
      DATA="$2"
      shift
      ;;
    --rsemOutputDir)
      RSEMOUTPUTDIR="$2"
      shift
      ;;
    --rsemResult)
      RSEMRESULT="$2" # "genes" or "isoforms"
      shift
      ;;
    --project)
      PROJECT="$2"
      shift
      ;;
    --output)
      OUTPUT="$2"
      shift
      ;;
    --animal)
      ANIMAL="$2"
      shift
      ;;
    --sex)
      SEX="$2"
      shift
      ;;
    --treatment)
      TREATMENT="$2"
      shift
      ;;
    --condition)
      CONDITION="$2"
      shift
      ;;
  esac
  shift
done

if [[ ! $RSEMRESULT =~ ^genes$|^isoforms$ ]]; then
  exit "Error: --rsemResult 'genes' or 'isoforms'!"
fi

RSEMOUTPUTFILES=()
INPUTFILES=$(ls $DATA)

for file in $(ls $RSEMOUTPUTDIR); do
  if [[ $file =~ $RSEMRESULT.results$ ]]; then
    RSEMOUTPUTFILES+=(RSEMOUTPUTDIR/$file)
  fi
done

echo RSEM output files = "${RSEMOUTPUTFILES[@]}"

if [[ ! ${#RSEMOUTPUTFILES[@]} == ${#INPUTFILES[@]} ]]; then
  exit "Didn't generate an RSEM output for each sequencing file!"
fi


METADATA="${OUTPUT}/${PROJECT}.cells.metadata.txt"
CELLCOUNTS="${OUTPUT}/${PROJECT}.cells.counts.txt"
GENEINF="${OUTPUT}/${PROJECT}.gene.information.txt"
TEMP="${OUTPUT}/${PROJECT}.cells.counts.temp.txt"

i=0
for file in "${RSEMOUTPUTFILES[@]}"; do
  echo "Generate Metadata file for $file"
  CELL=$(basename $file)
  RSEMOUTPUT="${file}.${RSEMRESULT}.results"
  if [[ $i == 0 ]]; then
    if [[ ! -f ${CELLCOUNTS} ]]; then
      awk -v cell="$CELL" 'NR==1 {print $1,cell} NR>1 {print $1,$5}' $RSEMOUTPUT \
        > "${CELLCOUNTS}"
    else
      paste $CELLCOUNTS \
        <(awk -v cell="$CELL" 'NR==1 {print cell} NR>1 {print $5}' $RSEMOUTPUT) \
        > ${TEMP}
      cat ${TEMP} > ${CELLCOUNTS}
      rm ${TEMP}
    fi
    if [[ ! -f $GENEINF ]]; then
      awk '{print $1,$3}' $RSEMOUTPUT > ${GENEINF}
    fi
    if [[ ! -f $METADATA ]]; then
      echo -e "cell\tanimal\tsex\tcondition\ttreatment" > ${METADATA}
    fi
    echo -e "$(basename $file)\t${ANIMAL}\t${SEX}\t${CONDITION}\t${TREATMENT}" \
      >> "${METADATA}"
    i=1
  else
    echo -e "$(basename $file)\t${ANIMAL}\t${SEX}\t${CONDITION}\t${TREATMENT}" \
      >> ${METADATA}
    paste $CELLCOUNTS \
      <(awk -v cell="$CELL" 'NR==1 {print cell} NR>1 {print $5}' $RSEMOUTPUT) \
      > $TEMP
    cat $TEMP > $CELLCOUNTS
    rm $TEMP
  fi
done
echo "DONE: Generated count table. Stored in ${CELLCOUNTS}."
echo "DONE: Generated metadata file. Stored in ${METADATA}."
exit

# If imputation is requested: Perform imputation on the raw count table
if [[ $IMPUTE =~ yes|Yes|YES|Y|y ]]; then
  echo "Perform imputation with scImpute to reduce dropouts."
  if [[ -z ${IMPUTOUT+x} ]]; then
    DIR=$(dirname "${CELLCOUNTS}")
    BASE=$(basename "${CELLCOUNTS}")
    IMPOUTPUT=$DIR/imputed.$BASE
  fi
  echo "Perform imputation on count table ${CELLCOUNTS}: Store in $IMPOUTPUT."
  Rscript impute_counts.R $CELLCOUNTS $IMPOUTPUT
fi

# shell script finished -- start R Pipeline
# Generating Heatmaps, PCA, tSNE, 
# Find Marker genes
# Do Differential expression analysis 
echo "Start Quality Control and Normalization with scater and scran in R!"
Rscript R_pipeline.R $CELLCOUNTS $METADATA
