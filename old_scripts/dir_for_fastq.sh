#!/bin/bash
# Daniel Schreyer
# Create Directory for each scRNA-seq fastq file
# Input: directory with all fastq files

while [[ $# -gt 0 ]]
do
  option="$1"
  case $option in
    --data)
      DATA="$2"
      shift
      ;;
  esac
  shift
done

FILES=($(ls ${DATA}))
for file in ${FILES[@]}; do
  if [[ $file =~ ^(.*)+(\.fq.*|\.fastq.*)$ ]]; then
    echo file = $file
    base=${BASH_REMATCH[1]}
    format=${BASH_REMATCH[2]}
    echo "base:" $base
    echo "format:" $format
    echo $DATA/$file
    mv $DATA/$file $DATA/temp
    mkdir $DATA/$base
    mv $DATA/temp $DATA/$base/$file
  else
    echo "$file is not a fastq file"
  fi
done

echo "Generate one directory for each fastq file!"
echo "Directories stored in: $DATA"
