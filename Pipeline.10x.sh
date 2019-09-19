#!/bin/bash
#
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
# TODO: Add STAR date format: "Jul 03 11:33:36 ..... started STAR run"

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
    --trimming) # Perform quality control?
      TRIMMING="$2"
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
    --threads) # how many threads are available
      THREADS="$2"
      shift
      ;;
    --umi-tools) # 
      UMITOOLS="$2"
      shift
      ;;
    --useCellranger) # 
      USECELLRANGER="$2"
      shift
      ;;
    --useSTARsolo) # 
      USESTARSOLO="$2"
      shift
      ;;
    --useUMItools) # 
      USEUMITOOLS="$2"
      shift
      ;;
    --cellranger) # 
      CELLRANGER="$2"
      shift
      ;;
    --cellrangerTranscriptome) # 
      CR_TRANSCRIPTOME="$2"
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
    --htseq) # executable HTSeq
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
    --star)
      STAR="$2"
      shift
      ;;
    --featureCounts)
      FEATURECOUNTS="$2"
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
    --whitelist)
      WHITELIST="$2"
      shift
      ;;
    --useLanes)
      LANES="$2"
      shift
      ;;
    --nGenes)
      NGENES="$2"
      shift
      ;;
    --nUMIs)
      NUMIS="$2"
      shift
      ;;
    --MAD)
      MAD="$2"
      shift
      ;;
    --abundantMT)
      ABUNDANTMT="$2"
      shift
      ;;
    --filterGenes)
      FILTERGENES="$2"
      shift
      ;;
    --normalize)
      NORMALIZE="$2"
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

# Set default parameters
if [[ -z ${GENOMEINDEX+x} || $GENOMEINDEX == "" ]]; then
  GENOMEINDEX="yes"
fi
if [[ -z ${QC+x} || $QC == "" ]]; then
  QC="yes"
fi
if [[ -z ${READ+x} || $READ == "" ]]; then
  READ="R2"
fi
if [[ -z ${BARCODE+x} || $BARCODE == "" ]]; then
  BARCODE="R1"
fi
if [[ -z ${LANES+x} || $LANES == "" ]]; then
  LANES="all"
fi
if [[ -z ${THREADS+x} || $THREADS == "" ]]; then
  THREADS=1
fi
if [[ -z ${TRIMOPTIONS+x} || $TRIMOPTIONS == "" ]]; then
  TRIMOPTIONS="TRAILING: 20 HEADING:20 MINLEN:75"
fi 
if [[ -z ${USEUMITOOLS+x} || $USEUMITOOLS == "" ]]; then
  USEUMITOOLS="no"
fi
if [[ -z ${USESTARSOLO+x} || $USESTARSOLO == "" ]]; then
  USESTARSOLO="no"
fi
if [[ -z ${USECELLRANGER+x} || $USECELLRANGER == "" ]]; then
  USECELLRANGER="yes"
fi
if [[ -z ${USESTARSOLO+x} || $USESTARSOLO == "" ]]; then
  USESTARSOLO="no"
fi
if [[ -z ${NGENES+x} || $NGENES == "" ]]; then
  NGENES="100"
fi
if [[ -z ${NUMIS+x} || $NUMIS == "" ]]; then
  NUMIS="125"
fi
if [[ -z ${MAD+x} || $MAD == "" ]]; then
  MAD="5"
fi
if [[ -z ${ABUNDANTMT+x} || $ABUNDANTMT == "" ]]; then
  ABUNDANTMT="0.5"
fi
if [[ -z ${FILTERGENES+x} || $FILTERGENES == "" ]]; then
  FILTERGENES="0.001"
fi
if [[ -z ${NORMALIZE+x} || $NORMALIZE == "" ]]; then
  NORMALIZE="yes"
fi


# create directories for essential outputs  
STARDIR=$OUTPUT/STAR_output    
MULTIQCDIR=$OUTPUT/MultiQC_output  
MERGED=${OUTPUT}/Files

echo "Create multiple directories in $OUTPUT"
make_dir $OUTPUT 
make_dir $STARDIR 
make_dir $MERGED


# Sequencing files
FILES=$(ls -d $DATA/*)

echo "${FILES[@]}"
# READ1 stores each fastq file with barcodes and umis, READ2 contains cDNA read
READ1_ARRAY=()
READ2_ARRAY=()
GZCOMPRESSED=()
BZ2COMPRESSED=()
echo "Start: Concatenating all Fastq files of a sample"
date


if [[ $LANES == "all" ]]; then
  use_lanes="L[0-9]*"
else
  OIFS=$IFS
  IFS=","
  use_lanes=($LANES)
  IFS=$OIFS
  LANES=$(echo $LANES | sed 's/,/_/g')
fi

for file in ${FILES[@]}; do
  for lane in ${use_lanes[@]}; do
    if [[ $lane =~ "[0-9]{3}" ]]; then
      lane=$(printf "L%03d" $lane)
    fi
    if [[ $file =~ ^(.*)/(.*)(_${lane}_)${BARCODE}_001\.(fastq|fq)(\.gz|\.bz2)* ]]; then
      R1_SAMPLE=${BASH_REMATCH[2]}
      R1_LANES=${BASH_REMATCH[3]}
      R1_FORMAT=${BASH_REMATCH[4]}
      R1_ZIP=${BASH_REMATCH[5]}
      NEWFILE_R1=${DATA}/${R1_SAMPLE}${R1_LANES}${BARCODE}_001.${R1_FORMAT}
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
      break
    elif [[ $file =~ ^(.*)/(.*)(_${lane}_)${READ}_001\.(fastq|fq)(\.gz|\.bz2)* ]]; then
      R2_SAMPLE=${BASH_REMATCH[2]}
      R2_LANES=${BASH_REMATCH[3]}
      R2_FORMAT=${BASH_REMATCH[4]}
      R2_ZIP=${BASH_REMATCH[5]}
      NEWFILE_R2=${DATA}/${R2_SAMPLE}${R2_LANES}${READ}_001.${R2_FORMAT}
      if [[ $R2_ZIP == ".gz" ]]; then
        echo "Decompress gz compressed File: $file"
        gunzip $file 
        GZCOMPRESSED+=($NEWFILE_R2)
      elif [[ $R2_ZIP == ".bz2" ]]; then
        echo "Decompress bzip2 compressed File: $file"
        bzip2 -d $file
        BZ2COMPRESSED+=($NEWFILE_R2)
      fi
      R2_ARRAY+=($NEWFILE_R2)
      break
    fi
  done
done

R1=${MERGED}/${R1_SAMPLE}_${LANES}_${BARCODE}_001.${R1_FORMAT}
R2=${MERGED}/${R2_SAMPLE}_${LANES}_${READ}_001.${R2_FORMAT}

if [[ ${#R2_ARRAY[@]} -eq 1 ]]; then
  R2=${R2_ARRAY}
  R1=${R1_ARRAY}
  echo "Only one sequencing file"
elif [[ ${#R2_ARRAY[@]} -gt 1 ]]; then
  echo "Concatenate Read 1 Files"
  cat ${R1_ARRAY[@]} > $R1
  echo "Concatenate Read 2 Files"
  cat ${R2_ARRAY[@]} > $R2
  echo "Concatenate ${R2_ARRAY[@]}"
else
  echo "No read file is stored in $DATA!"
  help_message
fi
echo "Finished: Stored Fastq file in $MERGED"

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
  FASTQCDIR=$OUTPUT/FastQC_output    
  make_dir $FASTQCDIR
  echo "Perform quality control"
  echo "Quality Control of ${R2}"
  $FASTQC -o $FASTQCDIR -t $THREADS ${R2}
  echo "Performed quality control"
fi

if [[ $USECELLRANGER == "yes" ]]; then
  $CELLRANGER count \
    --id=$(basename $DATA) \
    --sample=$(basename $DATA) \
    --transcriptome=$CR_TRANSCRIPTOME \
    --fastqs=$DATA
fi


WHITEZIP=0
if [[ $WHITELIST =~ .*.gz ]]
then
  echo "unzip whitelist: $WHITELIST"
  WHITELIST_UNZIP=$(echo "${WHITELIST}" | sed 's/.gz$//')
  if [[ ! -f $WHITELIST_UNZIP ]]; then
    echo "Unzip whitelist file: $WHITELIST"
    gunzip --keep "${WHITELIST}"
    WHITEZIP=1
  fi
  WHITELIST=$WHITELIST_UNZIP
fi

# Umi-Tools
# Identify correct cell barcodes
umi_whitelist="no"
if [[ umi_whitelist == "yes" ]]; then
  WHITELIST=$UMITOOLSDIR/${R1_SAMPLE}.whitelist.txt
  echo "Identify correct cell barcodes with Umi-Tools!"
  echo "START Umi-Tools whitelist: $R1"
  date
  $UMITOOLS whitelist \
    --stdin $R1 \
    --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
    --method=umis \
    --log2stderr > $WHITELIST
  echo "END Umi-Tools whitelist: $R1" 
  date
fi

if [[ $use_umi == "yes" ]]; then
  UMITOOLSDIR=${OUTPUT}/Umi-Tools
  make_dir $UMITOOLSDIR
  EXTRACTED_UMIS=()
  DEMUX_FILES=()
  echo "Umi-Tools: $R2"
  echo "Using Umi-Tools to extract Barcode"
  date 
  EXT="_extracted"
  R1_EXT=$UMITOOLSDIR/${R1_SAMPLE}${R1_LANES}${BARCODE}${EXT}_001.${R1_FORMAT}${R1_ZIP}
  R2_EXT=$UMITOOLSDIR/${R2_SAMPLE}${R2_LANES}${READ}${EXT}_001.${R2_FORMAT}${R2_ZIP}

  # if w != yes -> no umi_tools whitelist
  if [[ $umi_whitelist == "yes" ]]; then
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
  R1=$R1_EXT
  R2=$R2_EXT
fi
# Trimming 10x sequencing read2
if [[ $TRIMMING == "yes" ]]; then
  TRIMDIR=$OUTPUT/trimmomatic_output 
  make_dir $TRIMDIR
  echo "Start trimming reads with Trimmomatic!"
  echo "Input: $R2"
  echo "OUTPUT: R2_TRIM"
  date
  TRIM="_trim"
  R2_TRIM=$TRIMDIR/${R2_SAMPLE}${R2_LANES}${READ}${EXT}${TRIM}_001.${R2_FORMAT}${R2_ZIP}
  java -jar $TRIMMOMATIC \
    SE -phred33 $R2 $R2_TRIM \
    -threads ${THREADS} \
    $TRIM_OPTIONS
  echo "Finished trimming!"
  echo "Stored trimmed reads in $R2_TRIM"
  date
  R2=${R2_TRIM}
fi

# Quality Control of trimmed single end files with FastQC
if [[ $QC == "yes" && $TRIMMING == "yes" ]]; then
  echo "Start Quality control: $R2"
  date
  $FASTQC -o $FASTQCDIR -t $THREADS ${R2}
  echo "Finished quality control: $R2"
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

# with anno file: generate indices with it | without: generate indices without
if [[ ${GENOMEINDEX} == "no" ]]; then
  make_dir $INDICESDIR  
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

R2_STAROUT=$STARDIR/$(basename $R2 fastq)

if [[ $use_umi == "yes" ]]; then
  echo "Start STAR alignment: $R2"
  echo "STAR Output = $R2_STAROUT"
  $STAR --runThreadN $THREADS \
    --genomeDir "${INDICESDIR}" \
    --readFilesIn $R2 \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix $R2_STAROUT 
  echo "Saved STAR output of $R2 in $R2_STAROUT"
  echo "Finished STAR Alignment!"
  date
fi
if [[ $USESTARSOLO == "yes" ]]; then
  echo "Start STARsolo: Mapping, Demultiplexing and gene quantification"
  echo "Input: $R1 and $R2"
  $STAR --runThreadN $THREADS \
    --genomeDir "${INDICESDIR}" \
    --readFilesIn $R2 $R1 \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix $R2_STAROUT \
    --soloType Droplet \
    --soloCBwhitelist $WHITELIST
  echo "Finished STARsolo run with $R1 and $R2"
fi

# Count reads per gene with featureCounts
if [[ $use_umi == "yes" ]]; then
  COUNTS=${OUTPUT}/counts
  make_dir $COUNTS
  INPUT=${R2_STAROUT}Aligned.out.bam
  COUNTOUT=${COUNTS}/$(basename $INPUT out.bam)counts.txt
  echo "Start counting features with featureCounts!"
  echo "Input: $INPUT"
  echo "Output: $COUNTOUT"
  $FEATURECOUNTS \
    -a $ANNOTATION \
    -o $COUNTOUT \
    -R BAM $INPUT \
    -T ${THREADS}
  SAMINPUT=${INPUT}.featureCounts.bam
  SAMOUTPUT=$COUNTS/$(basename $COUNTOUT counts.txt)Aligned.out.featureCounts.sorted.bam  
  echo "Start: Samtools sort!"
  echo "Input: $SAMINPUT"
  echo "Output: $SAMOUTPUT"
  date
  $SAMTOOLS sort --threads ${THREADS} ${SAMINPUT} \
    -o ${SAMOUTPUT}
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
fi

echo "Gzip File: $R1 and $R2"
date
gzip -f $R1 $R2
echo "Finished: Gzip File $R1 and $R2"
date

# remove unzipped annotation file - restore original file structure
if [[ $ANNOZIP == 1 ]]; then
  echo "Remove unzipped annotation file"
  rm $ANNOTATION
fi
# Remove unzipped whitelist file
if [[ $WHITEZIP == 1 ]]; then
  echo "Remove unzipped whitelist file: $WHITELIST"
  rm $WHITELIST
fi

echo "Compress previous compressed files"
for file in ${GZCOMPRESSED[@]}; do
  echo "Gzip File: $file"
  gzip -f $file
done

for file in ${BZ2COMPRESSED[@]}; do
  echo "Bz2 File: $file"
  bzip2 -f $file
done

echo "Finished: Compress Files"
echo "End of Pipeline"

