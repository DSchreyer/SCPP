#!/bin/bash
#
# Daniel Schreyer
# Single Cell Analysis Pipeline
# Input: 10x single cell data
# Output:   1. Cell count table
#           2. Metadata file
#           3. Differentially expressed genes
#           4. PCA, tSNE
#           5. Marker genes for sequenced cell types

# TODO: ADD HTSeq-count
# TODO: Change Variable names for feature counts and rsem --> umi-tools count
# TODO: Start R-script from pipeline

##########################
# print out help message, if an error occured
# Arguments: 
#   None
# Returns:
#   exit program
##########################
help_message(){
  echo "Define Arguments in start.pipeline.sh"
  exit
}

while [[ $# -gt 0 ]]
do
  option="$1"
  case $option in
    --qualityControl) # Perform quality control?
      QC="$2"
      shift # shift arguments to the left = $2 -> $1
      ;;
    --projectName) # Project name to generate specific file names
      PROJECT="$2"
      shift # shift arguments to the left = $2 -> $1
      ;;
    --genome) # reference genome file
      GENOME="$2"
      shift 
      ;;
    --annotation) # annotation file
      ANNOTATION="$2"
      shift 
      ;;
    --indicesDir)
      INDICESDIR="$2" # path to directory where genome indices are stored
      shift
      ;;
    --index) # "yes": indices are available, "no": generate indices 
      GENOMEINDEX="$2"
      shift
      ;;
    --data)  # path to data directory with read files
      DATA="$2"
      shift
      ;;
    --read)
      READ="$2"
      shift
      ;;
    --barcode) # Specify "R1" or "R2", which is the Barcode + UMI read
      BARCODE="$2"
      shift
      ;;
    --num-cells) # Number of cells sequenced per sample
      NUMCELLS="$2"
      shift
      ;;
    --threads) # how many threads are available
      THREADS="$2"
      shift
      ;;
    --umi-tools) # how many threads are available
      UMITOOLS="$2"
      shift
      ;;
    --samtools) # executable trimmomatic path
      SAMTOOLS="$2"
      shift
      ;;
    --trimmomatic) # executable trimmomatic path
      TRIMMOMATIC="$2"
      shift
      ;;
    --HTSeq) # executable HTSeq
      HTSEQ="$2"
      shift
      ;;
    --fastqc) # executable fastqc
      FASTQC="$2"
      shift
      ;;
    --multiqc)
      MULTIQC="$2"
      shift
      ;;
    --STAR)
      STAR="$2"
      shift
      ;;
    --featureCounts)
      FEATURECOUNTS="$2"
      shift
      ;;
    --impute)
      IMPUTE="$2" # "no": no imputation, "yes": imputation 
      shift
      ;;
    --output) # Output directory
      OUTPUT="$2"
      shift
      ;;
    --trimOptions) # Output directory
      TRIMOPTIONS="$2"
      shift
      ;;
    --starOptions) # Output directory
      STAROPTIONS="$2"
      shift
      ;;
    *)
      echo -e "ERROR: \"${option}\" is an unknown option!"
      help_message
      ;;
  esac
  shift
done

# creates directory, if it does not exist
##########################
# creates missing directory 
# Arguments: 
#   directory path
# Returns:
#   None 
##########################
make_dir (){
  local dir="${1}"
  if [[ ! -d $dir ]]; then
    mkdir -p $dir
    echo "Created $dir!"
  fi
}

# Controls existance of file
##########################
# Arguments: 
#   file path
# Returns:
#   exits shell script with error message 
##########################
file_exists (){
  local file="${1}"
  if [[ ! -f $file ]]; then
    echo "Error: $file does not exist!"
    exit
  fi
}

    
# if [[ ${#EXPECTEDARGS[@]} != ${#ARGARRAY[@]} ]]; then
#   echo "One or multiple argument/s is/are missing!"
#   help_message
# fi

# create directories for essential outputs  
TRIMDIR=$OUTPUT/trimmomatic_output 
FASTQCDIR=$OUTPUT/FastQC_output    
STARDIR=$OUTPUT/STAR_output    
MULTIQCDIR=$OUTPUT/MultiQC_output  
UMITOOLSDIR=${OUTPUT}/Umi-Tools
COUNTS=${OUTPUT}/counts
MERGED=${OUTPUT}/Files

echo "Create multiple directories in $OUTPUT"
make_dir $OUTPUT 
make_dir $TRIMDIR 
make_dir $FASTQCDIR   
make_dir $OUTPUT  
make_dir $INDICESDIR  
make_dir $STARDIR 
make_dir $UMITOOLSDIR 
make_dir $COUNTS
make_dir $MERGED

# Concatenate fastq files together
FILES=$(ls -d $DATA/*)

# READ1 stores each fastq file with barcodes and umis, READ2 contains cDNA read
READ1_ARRAY=()
READ2_ARRAY=()
GZCOMPRESSED=()
BZ2COMPRESSED=()
echo "Start: Concatenating all Fastq files of a sample"
date

# Lanes to use: "all"
# Lanes to use: 1,2,3,4,5

if [[ $LANES == "all" ]]; then
  use_lanes="[0-9]*"
else
  OIFS=$IFS
  IFS=","
  use_lanes=($LANES)
  IFS=$OIFS
fi


for file in ${FILES[@]}; do
  for lane in ${use_lanes[@]}; do
  if [[ $file =~ ^(.*)/(.*)(_L[0-9]{3}_)${BARCODE}_001\.(fastq|fq)(\.gz|\.bz2)* ]]; then
    R1_SAMPLE=${BASH_REMATCH[2]}
    R1_LANES=${BASH_REMATCH[3]}
    R1_FORMAT=${BASH_REMATCH[4]}
    R1_ZIP=${BASH_REMATCH[5]}
    
    NEWFILE_R1=${DATA}/${R1_SAMPLE}${R1_LANES}${BARCODE}_001.${R1_FORMAT}${R1_ZIP}
    if [[ $R1_ZIP == ".gz" ]]; then
      echo "Decompress gz compressed File: $file"
      gunzip $file 
      GZCOMPRESSED+=($NEWFILE_R1)
    elif [[ $R1_ZIP == ".bz2" ]]; then
      echo "Decompress bzip2 compressed File: $file"
      bzip2 -d $file
      BZ2COMPRESSED+=($NEWFILE_R1)
    fi
    R1_ARRAY+=($NEWFILE_R1)
    R1=${MERGED}/${R1_SAMPLE}_ALL_${BARCODE}_001.${R1_FORMAT}
  elif [[ $file =~ ^(.*)/(.*)(_L[0-9]{3}_)${READ}_001\.(fastq|fq)(\.gz|\.bz2)* ]]; then
    R2_SAMPLE=${BASH_REMATCH[2]}
    R2_LANES=${BASH_REMATCH[3]}
    R2_FORMAT=${BASH_REMATCH[4]}
    R2_ZIP=${BASH_REMATCH[5]}
    NEWFILE_R2=${DATA}/${R2_SAMPLE}${R2_LANES}${READ}_001.${R2_FORMAT}${R2_ZIP}
  exit
fi


TRIMMING="no"
if [[ $TRIMMING = "no" ]]; then
  echo "Don't perform Trimming with Trimmomatic"
elif [[ $TRIMMING == "yes" ]]; then
  echo "Perform Trimming with Trimmomatic"
else
  echo "Perform Trimming with Trimmomatic?"
  echo "'yes' | 'no'"
  exit
fi
# Quality Control with all files in the data Directory --- FastQC
if [[ $QC == "yes" ]]; then
  echo "Perform quality control"
  echo "Quality Control of ${R2}"
  $FASTQC -o $FASTQCDIR -t $THREADS ${R2}
  echo "Performed quality control"
fi

# Umi-Tools
# Identify correct cell barcodes

EXTRACTED_UMIS=()
DEMUX_FILES=()
echo "Umi-Tools: $R2"
echo "Using Umi-Tools to extract Barcode"
date 
EXT="_extracted"
R1_EXT=$UMITOOLSDIR/${R1_SAMPLE}${R1_LANES}${BARCODE}${EXT}_001.${R1_FORMAT}${R1_ZIP}
WHITELIST=$UMITOOLSDIR/${R1_SAMPLE}.whitelist.txt
R2_EXT=$UMITOOLSDIR/${R2_SAMPLE}${R2_LANES}${READ}${EXT}_001.${R2_FORMAT}${R2_ZIP}

# if w != yes -> no umi_tools whitelist
w="no"
if [[ $w == "yes" ]]; then
  # Identify correct cell barcodes
  echo "Identify correct cell barcodes with Umi-Tools!"
  echo "START Umi-Tools whitelist: $R1"
  date
  $UMITOOLS whitelist \
    --stdin $R1 \
    --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
    --method "reads" \
    --log2stderr > $WHITELIST
  echo "END Umi-Tools whitelist: $R1" 
  date
  echo "Extract barcodes and UMIs and add to read names"
  echo "START Umi-Tools extract: $R1 and $R2"
  date
  $UMITOOLS extract \
    --stdin $R1 \
    --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
    --stdout $R1_EXT \
    --read2-in $R2 \
    --read2-out $R2_EXT \
    --quality-encoding="phred33" \
    --quality-filter-threshold 15 \
    --filter-cell-barcode \
    --error-correct-cell \
    --whitelist=$WHITELIST

  echo "END Umi-Tools extract: $R1 and $R2"
  echo "Stored: $R1_EXT and $R2_EXT"
  date
else
  echo "Extract barcodes and UMIs and add to read names"
  echo "START Umi-Tools extract: $R1 and $R2"
  date
  $UMITOOLS extract \
    --stdin $R1 \
    --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
    --stdout $R1_EXT \
    --read2-in $R2 \
    --read2-out $R2_EXT \
    --quality-filter-threshold 15 
  echo "END Umi-Tools extract: $R1 and $R2"
  echo "Stored: $R1_EXT and $R2_EXT"
  date
fi
echo "Finished Umi-Tools Whitelist and extract. Next Steps are Alignment and counting"

# Trimming 10x sequencing read2
if [[ $TRIMMING == "yes" ]]; then
  echo "Start trimming reads with Trimmomatic!"
  echo "Input: $R2_EXT"
  echo "OUTPUT: R2_TRIM"
  date
  TRIM="_trim"
  R2_TRIM=$TRIMDIR/${R2_SAMPLE}${R2_LANES}${READ}${EXT}${TRIM}_001.${R2_FORMAT}${R2_ZIP}
  java -jar $TRIMMOMATIC \
    SE -phred33 $R2_EXT $R2_TRIM \
    -threads ${THREADS} \
    LEADING:20 TRAILING:20 MINLEN:75 
  echo "Finished trimming!"
  echo "Stored trimmed reads in $R2_TRIM"
date
else
  R2_TRIM="${R2_EXT}"
fi

# Quality Control of trimmed single end files with FastQC
if [[ $QC == "yes" && $TRIMMING == "yes" ]]; then
  echo "Start Quality control: $R2_TRIM"
  date
  $FASTQC -o $FASTQCDIR -t $THREADS ${R2_TRIM}
  echo "Finished quality control: $R2_TRIM"
  date
fi


# ANNOZIP=0: annotation file is unzipped, ANNOZIP=1: annotation file was zipped
ANNOZIP=0
# gunzip gzipped annotation file --- expects a .gtf annotation file
if [[ $ANNOTATION =~ .*gtf.gz$|.*gff.gz$ ]]; then
  file_exists $ANNOTATION
  ANNOTATION_unzipped=$(echo "${ANNOTATION}" | sed 's/.gz$//')
  if [[ -f $ANNOTATION_unzipped ]]; then
    ANNOTATION=$ANNOTATION_unzipped
  else
    echo  "gunzip annotation file ${ANNOTATION}"
    gunzip --keep "${ANNOTATION}"
    ANNOZIP=1
    ANNOTATION=$(echo "${ANNOTATION}" | sed 's/.gz$//')
    echo "Finished unzipping annotation file"
    date
  fi
fi
echo Annotation File = "${ANNOTATION}"

# Trimming 10x sequencing read2
if [[ $TRIMMING == "yes" ]]; then
  echo "Start trimming reads with Trimmomatic!"
  echo "Input: $R2_EXT"
  echo "OUTPUT: R2_TRIM"
  date
  TRIM="_trim"
  R2_TRIM=$TRIMDIR/${R2_SAMPLE}${R2_LANES}${READ}${EXT}${TRIM}_001.${R2_FORMAT}${R2_ZIP}
  java -jar $TRIMMOMATIC \
    SE -phred33 $R2_EXT $R2_TRIM \
    -threads ${THREADS} \
    $TRIM_OPTIONS
  echo "Finished trimming!"
  echo "Stored trimmed reads in $R2_TRIM"
  date
else
  R2_TRIM="${R2_EXT}"
fi

# Quality Control of trimmed single end files with FastQC
if [[ $QC == "yes" && $TRIMMING == "yes" ]]; then
  echo "Start Quality control: $R2_TRIM"
  date
  $FASTQC -o $FASTQCDIR -t $THREADS ${R2_TRIM}
  echo "Finished quality control: $R2_TRIM"
  date
fi

# gunzip gzipped reference genome fasta file
if [[ $GENOME =~ .*fa.gz ]]
then
  echo "unzip genome fasta file"
  gunzip --keep "${GENOME}"
  GENOME=$(echo "${GENOME}" | sed 's/.gz$//')
fi

# ANNOZIP=0: annotation file is unzipped, ANNOZIP=1: annotation file was zipped
ANNOZIP=0
# gunzip gzipped annotation file --- expects a .gtf annotation file
if [[ $ANNOTATION =~ .*gtf.gz$|.*gff.gz$ ]]; then
  file_exists $ANNOTATION
  ANNOTATION_unzipped=$(echo "${ANNOTATION}" | sed 's/.gz$//')
  if [[ -f $ANNOTATION_unzipped ]]; then
    ANNOTATION=$ANNOTATION_unzipped
  else
    echo  "gunzip annotation file ${ANNOTATION}"
    gunzip --keep "${ANNOTATION}"
    ANNOZIP=1
    ANNOTATION=$(echo "${ANNOTATION}" | sed 's/.gz$//')
    echo "Finished unzipping annotation file"
    date
  fi
fi


echo Annotation File = "${ANNOTATION}"
echo "Generate Genome Index"

# GENOMEINDEX stores user input, if genome indices were already generated
echo "${GENOMEINDEX}"
if [[ ! ${GENOMEINDEX} =~ yes|y|n|no|Yes|No|Y|N|NO|YES ]]; then
  echo "Please specify, if the STAR indices are available or not!"
  help_message
fi

# with anno file: generate indices with it | without: generate indices without
if [ -z ${GENOMEINDEX+x} ] || [[ ${GENOMEINDEX} =~ no|n|No|N|NO ]]; then
  if [[ -z ${ANNOTATION+x} ]]; then
    echo "No annotation file. Indexing with annotation file is recommended!"
    $STAR --runMode genomeGenerate \
      --genomeDir "${INDICESDIR}" \
      --genomeFastaFiles "${GENOME}" \
      --runThreadN ${THREADS}
  else
    if [[ $ANNOTATION =~ .*\.gtf ]]; then
      echo "Run genomeGenerate with gtf annotation file"
      $STAR --runMode genomeGenerate \
        --genomeDir "${INDICESDIR}" \
        --sjdbGTFfile "${ANNOTATION}" \
        --genomeFastaFiles "${GENOME}" \
        --runThreadN ${THREADS}
    elif [[ $ANNOTATION =~ .*\.gff ]]; then
      echo "Run genomeGenerate with gff annotation file"
      $STAR --runmode genomeGenerate \
        --genomeDir "${INDICESDIR}" \
        --sjdbGTFtagExonParentTranscript Parent \
        --genomeFastaFiles "${GENOME}" \
        --runThreadN "${THREADS}"
    else
      echo "Annotation file ${ANNOTATION} has the wrong format"
      echo "Enter .gtf or .gtf annotation file"
      help_message
    fi
  fi
  echo "Stored genome indices in ${INDICESDIR}"
  echo "Finished generating STAR indices!"
  date
fi


echo "Start STAR alignment: $R2_TRIM"
if [[ $TRIMMING == "yes" ]]; then
  R2_STAROUT=$STARDIR/${R2_SAMPLE}${R2_LANES}${READ}${EXT}${TRIM}.
  echo "STAR Output = $R2_STAROUT"
else
  R2_STAROUT=$STARDIR/${R2_SAMPLE}${R2_LANES}${READ}${EXT}.
fi

if [[ $R2_ZIP == ".bz2" ]]; then
  echo "bzipped"
  $STAR --runThreadN $THREADS \
    --genomeDir "${INDICESDIR}" \
    --readFilesIn $R2_TRIM \
    --readFilesCommand bunzip2 -c \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix $R2_STAROUT 
  echo "Saved STAR output of $R2_TRIM in $R2_STAROUT"
elif [[ $R2_ZIP == ".gz" ]]; then
  echo "gzipped"
  $STAR --runThreadN $THREADS \
    --genomeDir "${INDICESDIR}" \
    --readFilesIn $R2_TRIM \
echo "Finished: Samtools sort!"
date
echo "Start: Samtools index"
echo "Input: $SAMOUTPUT"
$SAMTOOLS index $SAMOUTPUT
file_exists ${SAMOUTPUT}
echo "Finished: Samtools index"


COUNTFILE=${COUNTS}/$(basename $SAMOUTPUT .bam).tsv.gz
echo "Start: Demultiplexing counts with Umi-tools count!"
echo "Input: $SAMOUTPUT"
echo "Output: $COUNTFILE"
date

$UMITOOLS count --per-gene \
  --gene-tag=XT --assigned-status-tag=XS \
  --wide-format-cell-counts \
  --per-cell -I ${SAMOUTPUT} -S ${COUNTFILE}

echo "Finished: Umi-tools count"
echo "Stored count table in $COUNTFILE"
date

echo "Compress previous compressed files"
for file in ${GZCOMPRESSED[@]}; do
  echo "Gzip File: $file"
  gzip $file
done

for file in ${BZ2COMPRESSED[@]}; do
  echo "Bz2 File: $file"
  bzip2 $file
done

echo "Finished: Compress Files"
date
# remove unzipped annotation file - restore original file structure
if [[ $ANNOZIP == 1 ]]; then
  echo "Remove unzipped annotation file"
  rm $ANNOTATION
fi
echo "End of Pipeline"
