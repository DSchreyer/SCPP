#!/bin/bash
#
#
# Daniel Schreyer
# Single Cell Analysis Pipeline
# Input: 10x single cell data

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
    --pipeline) # Perform quality control?
      PIPELINE="$2"
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
    --CRoptions) # Output directory
      CROPTIONS="$2"
      shift
      ;;
    --STARoptions) # Output directory
      STAROPTIONS="$2"
      shift
      ;;
    --STARwhitelist)
      STARWHITELIST="$2"
      shift
      ;;
    --UMITOOLSwhitelist)
      UMITOOLSWHITELIST="$2"
      shift
      ;;
    --genWhitelist)
      gen_whitelist="$2"
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
    --thresholdMT)
      THRESHOLDMT="$2"
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
if [[ -z ${PIPELINE+x} || $PIPELINE == "" ]]; then
  PIPELINE=$(pwd)/Pipeline.10x.sh
fi
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
  TRIMOPTIONS="TRAILING:20 HEADCROP:20 MINLEN:75"
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
if [[ -z ${TRIMMING+x} || $TRIMMING == "" ]]; then
  TRIMMING="no"
fi
if [[ -z ${UMITOOLSWHITELIST+x} ]]; then
  UMITOOLSWHITELIST=""
fi
if [[ -z ${STARWHITELIST+x} ]]; then
  STARWHITELIST=""
fi
if [[ -z ${gen_whitelist+x} ]]; then
  gen_whitelist="yes"
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
if [[ -z ${THRESHOLDMT+x} || $THRESHOLDMT == "" ]]; then
  THRESHOLDMT="0.5"
fi
if [[ -z ${FILTERGENES+x} || $FILTERGENES == "" ]]; then
  FILTERGENES="0.001"
fi
if [[ -z ${NORMALIZE+x} || $NORMALIZE == "" ]]; then
  NORMALIZE="yes"
fi
# Create Output Directory
make_dir $OUTPUT 

# create time variable
now=$(date +%d-%m-%Y" "%X" |")

# Sequencing files
FILES=$(ls -d $DATA/*)
# Sample Name
sample=$(basename $DATA)
PIPELINE_DIR=$(dirname $PIPELINE)
# READ1 stores each fastq file with barcodes and umis, READ2 contains cDNA read
READ1_ARRAY=()
READ2_ARRAY=()
GZCOMPRESSED=()
BZ2COMPRESSED=()

echo "$now START Single-Cell RNA-seq Processing Pipeline"
echo "$now Input: $sample"
echo "$now Output: $OUTPUT"
if [[ $LANES == "all" ]]; then
  use_lanes="L[0-9]*"
  LANES="_all_"
else
  OIFS=$IFS
  IFS=","
  use_lanes=($LANES)
  IFS=$OIFS
  CRLANES=$LANES
  LANES=$(echo $LANES | sed 's/,/_/g')
fi
if [[ $USESTARSOLO == "yes" || $USEUMITOOLS == "yes" ]]; then
  for file in ${FILES[@]}; do
    MERGED=$OUTPUT/Files
    make_dir $MERGED
    for lane in ${use_lanes[@]}; do
      if [[ $lane =~ "[0-9]{3}" ]]; then
        lane=$(printf "L%03d" $lane)
      fi
      if [[ $file =~ ^(.*)/(.*)(_${lane}_)${BARCODE}_001\.(fastq|fq)(\.gz|\.bz2)* ]]; then
        R1_SAMPLE=${BASH_REMATCH[2]}
        R1_LANES=${BASH_REMATCH[3]}
        R1_FORMAT=${BASH_REMATCH[4]}
        R1_ZIP=${BASH_REMATCH[5]}
        NEWFILE_R1=${MERGED}/${R1_SAMPLE}${R1_LANES}${BARCODE}_001.${R1_FORMAT}
        if [[ $R1_ZIP == ".gz" ]]; then
          echo "Decompress gz compressed File: $file"
          gunzip --force --keep -c $file > $NEWFILE_R1
          GZCOMPRESSED+=($NEWFILE_R1)
        elif [[ $R1_ZIP == ".bz2" ]]; then
          echo "Decompress bzip2 compressed File: $file"
          bzip2 --decompress --keep -c $file > $NEWFILE_R1
          BZ2COMPRESSED+=($NEWFILE_R1)
        fi
        R1_ARRAY+=($NEWFILE_R1)
        break
      elif [[ $file =~ ^(.*)/(.*)(_${lane}_)${READ}_001\.(fastq|fq)(\.gz|\.bz2)* ]]; then
        R2_SAMPLE=${BASH_REMATCH[2]}
        R2_LANES=${BASH_REMATCH[3]}
        R2_FORMAT=${BASH_REMATCH[4]}
        R2_ZIP=${BASH_REMATCH[5]}
        NEWFILE_R2=${MERGED}/${R2_SAMPLE}${R2_LANES}${READ}_001.${R2_FORMAT}
        if [[ $R2_ZIP == ".gz" ]]; then
          echo "Decompress gz compressed File: $file"
          gunzip --force --keep -c $file > $NEWFILE_R2
          GZCOMPRESSED+=($NEWFILE_R2)
        elif [[ $R2_ZIP == ".bz2" ]]; then
          echo "Decompress bzip2 compressed File: $file"
          bzip2 -d --force --keep -c $file > $NEWFILE_R2
          BZ2COMPRESSED+=($NEWFILE_R2)
        fi
        R2_ARRAY+=($NEWFILE_R2)
      fi
    done
  done
fi

# Quality Control with all files in the data Directory --- FastQC
if [[ $QC == "yes" ]]; then
  FASTQCDIR=$OUTPUT/FastQC_output    
  make_dir $FASTQCDIR
  echo "$now Perform quality control"
  for file in ${R2_ARRAY[@]}; do
    echo "Start quality control of $file"
    dir=${FASTQCDIR}/$(basename $file .fastq)
    make_dir $dir
    $FASTQC -o $dir -t $THREADS $file
    echo "$now Performed quality control of $file."
  done
fi

# Trimming 10x sequencing read2
if [[ $TRIMMING == "yes" ]]; then
  TRIMDIR=$OUTPUT/trimmomatic_output 
  make_dir $TRIMDIR
  TRIM_ARRAY=()
  echo "$now Start trimming reads with Trimmomatic!"
  for file in ${R2_ARRAY[@]}; do
    file_trim=$TRIMDIR/$(basename $file .fastq)_trim.fastq
    echo "Perform trimming with $file"
    echo "Options: $TRIM_OPTIONS"
    echo "OUTPUT: $file_trim"
    TRIM="_trim"
    java -jar $TRIMMOMATIC \
      SE -phred33 $file $file_trim \
      -threads ${THREADS} \
      $TRIMOPTIONS
    echo "$now Finished trimming. Stored in $file_trim"
    TRIM_ARRAY+=($file_trim)
  done
  R2_ARRAY=${TRIM_ARRAY[@]}
  # Quality Control of trimmed single end files with FastQC
  if [[ $QC == "yes" ]]; then
    for file in ${R2_ARRAY[@]}; do
      echo "$now Start quality control of $file"
      dir_trim=${FASTQCDIR}/$(basename $file .fastq)
      make_dir $dir_trim
      $FASTQC -o $dir_trim -t $THREADS $file
      echo "$now Performed quality control of $file."
    done
  fi
fi



# Merge fastq files together. Only used in UMI tools or STARsolo branch
if [[ $USEUMITOOLS == "yes" || $USESTARSOLO == "yes" ]]; then
  make_dir $MERGED
  R1=${MERGED}/${R1_SAMPLE}${LANES}${BARCODE}_001.${R1_FORMAT}
  R2=${MERGED}/${R2_SAMPLE}${LANES}${READ}_001.${R2_FORMAT}
  STARDIR=$OUTPUT/STAR_output
  make_dir $STARDIR
  if [[ ${#R2_ARRAY[@]} -eq 1 ]]; then
    R2=${R2_ARRAY}
    R1=${R1_ARRAY}
  elif [[ ${#R2_ARRAY[@]} -gt 1 ]]; then
    echo "$now Concatenate Read 1 Files"
    cat ${R1_ARRAY[@]} > $R1
    R1_STAR=$R1
    echo "$now Concatenate Read 2 Files"
    cat ${R2_ARRAY[@]} > $R2
    R2_STAR=$R2
  else
    echo "$now No read file is stored in $DATA!"
    help_message
  fi
fi
echo "$now Stored Fastq file in $MERGED"

# Use CellRanger branch by specifying useCellranger "yes"
if [[ $USECELLRANGER == "yes" ]]; then
  if [[ $LANES == "_all_" ]]; then
    CRLANES=""
  else
    CRLANES="--lanes $CRLANES"
  fi
  CROUTPUT=$OUTPUT/CellRanger/
  make_dir $CROUTPUT
  old_pwd=$(pwd)
  cd $CROUTPUT
  echo "$now Start CellRanger"
  $CELLRANGER count \
    --id=$sample \
    --sample=$sample \
    --transcriptome=$CR_TRANSCRIPTOME \
    --fastqs=$DATA \
    --nosecondary \
    $CRLANES \
    $CROPTIONS
  echo "$now Finished CellRanger"
fi
# Unzip given whitelist for UMI-tools
UMIWHITEZIP=0
if [[ $USEUMITOOLS && $UMITOOLSWHITELIST =~ .*.gz ]]; then
  WHITELIST_UNZIP=$(echo "${UMITOOLSWHITELIST}" | sed 's/.gz$//')
  if [[ ! -f $WHITELIST_UNZIP ]]; then
    echo "$now Unzip whitelist file: $UMITOOLSWHITELIST"
    gunzip --keep --force "${UMITOOLSWHITELIST}"
    UMIWHITEZIP=1
    echo "$now Finished."
  fi
  UMIWHITELIST=$WHITELIST_UNZIP
fi
if [[ $USESTARSOLO && $STARWHITELIST =~ .*.gz ]]; then
  STAR_WHITELIST_UNZIP=$(echo "${STARWHITELIST}" | sed 's/.gz$//')
  if [[ ! -f $STAR_WHITELIST_UNZIP ]]; then
    echo "$now Unzip whitelist file: $STARWHITELIST"
    gunzip --keep --force "${STARWHITELIST}"
    STARWHITEZIP=1
    echo "$now Finished."
  fi
  STARWHITELIST=$STAR_WHITELIST_UNZIP
fi

# Umi-Tools
# Identify correct cell barcodes
if [[ $gen_whitelist == "yes" ]]; then
  if [[ ($USEUMITOOLS == "yes" && $UMITOOLSWHITELIST == "") || ($USESTARSOLO == "yes" && $STARWHITELIST == "") ]]; then
    UMITOOLSDIR=${OUTPUT}/Umi-Tools
    make_dir $UMITOOLSDIR
    WHITELIST=$UMITOOLSDIR/${sample}.whitelist.txt
    echo "$now START Umi-Tools whitelist: $R1"
    $UMITOOLS whitelist \
      --stdin $R1 \
      --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
      --method=umis \
      --log2stderr > $WHITELIST
    echo "$now END Umi-Tools whitelist: $R1" 
    if [[ $UMITOOLSWHITELIST == "" && $USEUMITOOLS == "yes" ]]; then
      UMITOOLSWHITELIST=$WHITELIST
    fi
    if [[ $STARWHITELIST == "" && $USESTARSOLO == "yes" ]]; then
      STARWHITELIST=$UMITOOLSDIR/${sample}.STAR.whitelist.txt
      cat $WHITELIST | awk '{print $1}' > $STARWHITELIST
    fi
  fi
fi


if [[ $USEUMITOOLS == "yes" ]]; then
  UMITOOLSDIR=${OUTPUT}/Umi-Tools
  make_dir $UMITOOLSDIR
  EXTRACTED_UMIS=()
  DEMUX_FILES=()
  EXT="_extracted"
  R1_EXT=$UMITOOLSDIR/${R1_SAMPLE}${LANES}${BARCODE}${EXT}_001.${R1_FORMAT}${R1_ZIP}
  R2_EXT=$UMITOOLSDIR/${R2_SAMPLE}${LANES}${READ}${EXT}_001.${R2_FORMAT}${R2_ZIP}
  echo "$now START Umi-Tools extract: $R1 and $R2"
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
    --whitelist=$UMITOOLSWHITELIST

  echo "$now END Umi-Tools extract: $R1 and $R2"
  echo "$now Stored: $R1_EXT and $R2_EXT"
  R1=$R1_EXT
fi

# gunzip gzipped reference genome fasta file
if [[ $GENOME =~ .*fa.gz ]]
then
  echo "$now START Unzip genome fasta file"
  gunzip --keep "${GENOME}"
  GENOME=$(echo "${GENOME}" | sed 's/.gz$//')
  echo "$now Finished."
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
    echo  "$now Gunzip annotation file ${ANNOTATION}"
    gunzip --keep "${ANNOTATION}"
    ANNOZIP=1
    ANNOTATION=$(echo "${ANNOTATION}" | sed 's/.gz$//')
    echo "$now Finished unzipping annotation file"
  fi
fi
echo "$now" Annotation File = "${ANNOTATION}"
# GENOMEINDEX stores user input, if genome indices were already generated

# with anno file: generate indices with it | without: generate indices without
if [[ ${GENOMEINDEX} == "no" ]]; then
  echo "$now START generate genome index"
  make_dir $INDICESDIR  
  if [[ -z ${ANNOTATION+x} ]]; then
    echo "$now No annotation file. Indexing with annotation file is recommended!"
    $STAR --runMode genomeGenerate \
      --genomeDir "${INDICESDIR}" \
      --genomeFastaFiles "${GENOME}" \
      --runThreadN ${THREADS}
  else
    if [[ $ANNOTATION =~ .*\.gtf ]]; then
      echo "$now Run genomeGenerate with gtf annotation file"
      $STAR --runMode genomeGenerate \
        --genomeDir "${INDICESDIR}" \
        --sjdbGTFfile "${ANNOTATION}" \
        --genomeFastaFiles "${GENOME}" \
        --runThreadN ${THREADS}
    elif [[ $ANNOTATION =~ .*\.gff ]]; then
      echo "$now Run genomeGenerate with gff annotation file"
      $STAR --runmode genomeGenerate \
        --genomeDir "${INDICESDIR}" \
        --sjdbGTFtagExonParentTranscript Parent \
        --genomeFastaFiles "${GENOME}" \
        --runThreadN "${THREADS}"
    else
      echo "$now Annotation file ${ANNOTATION} has the wrong format"
      echo "$now Enter .gtf or .gtf annotation file"
      help_message
    fi
  fi
  echo "$now Stored genome indices in ${INDICESDIR}"
  echo "$now Finished generating STAR indices!"
fi

if [[ $USEUMITOOLS == "yes" ]]; then
  R2_STAROUT=$STARDIR/$(basename $R2 fastq)
  echo "$now START STAR alignment: $R2"
  echo "$now STAR Output = $R2_STAROUT"
  $STAR --runThreadN $THREADS \
    --genomeDir "${INDICESDIR}" \
    --readFilesIn $R2 \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix $R2_STAROUT \
    $STAROPTIONS
  echo "$now Finished STAR Alignment!"
fi
if [[ $USESTARSOLO == "yes" ]]; then
  STARSOLO_DIR=$OUTPUT/STARsolo
  make_dir $STARSOLO_DIR
  R2_SOLO_OUT=$STARSOLO_DIR/${sample}${LANES}
  echo "$now START STARsolo: Mapping, Demultiplexing and gene quantification"
  echo "$now Input: R2: $R2_STAR & R1: $R1_STAR"
  echo "$now Whitelist: $STARWHITELIST"
  $STAR --runThreadN $THREADS \
    --genomeDir "${INDICESDIR}" \
    --readFilesIn $R2_STAR $R1_STAR \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix $R2_SOLO_OUT \
    --soloType Droplet \
    --soloCBwhitelist $STARWHITELIST \
    $STAROPTIONS
  echo "$now Finished STARsolo run with $R1_STAR and $R2_STAR"
fi

# Count reads per gene with featureCounts
if [[ $USEUMITOOLS == "yes" ]]; then
  COUNTS=${UMITOOLSDIR}/counts
  make_dir $COUNTS
  INPUT=${R2_STAROUT}Aligned.out.bam
  COUNTOUT=${COUNTS}/$(basename $INPUT out.bam)counts.txt
  echo "$now START counting features with featureCounts!"
  echo "$now Input: $INPUT"
  echo "$now Output: $COUNTOUT"
  $FEATURECOUNTS \
    -a $ANNOTATION \
    -o $COUNTOUT \
    -R BAM $INPUT \
    -T ${THREADS}
  SAMINPUT=$COUNTS/$(basename $INPUT).featureCounts.bam
  SAMOUTPUT=$COUNTS/$(basename $SAMINPUT .bam).sorted.bam  
  echo "$now Start: Samtools sort!"
  echo "$now Input: $SAMINPUT"
  echo "$now Output: $SAMOUTPUT"
  $SAMTOOLS sort --threads ${THREADS} ${SAMINPUT} \
    -o ${SAMOUTPUT}
  echo "$now Finished: Samtools sort!"
  echo "$now START Samtools index"
  echo "$now Input: $SAMOUTPUT"
  $SAMTOOLS index $SAMOUTPUT
  file_exists ${SAMOUTPUT}
  echo "$now Finished: Samtools index"

  COUNTFILE=${COUNTS}/$(basename $COUNTOUT .txt).matrix.tsv
  echo "$now Start: Demultiplexing counts with Umi-tools count!"
  echo "$now Input: $SAMOUTPUT"
  echo "$now Output: $COUNTFILE"

  $UMITOOLS count --per-gene \
    --gene-tag=XT --assigned-status-tag=XS \
    --per-cell -I ${SAMOUTPUT} -S ${COUNTFILE}

  echo "$now Finished: Umi-tools count"
  echo "$now Stored count table in $COUNTFILE"
fi
echo $(pwd)
cd $old_pwd
echo $(pwd)
echo "old pwd $old_pwd"
QC_SCRIPT=$PIPELINE_DIR/sc_analysis_qc.R
if [[ $USEUMITOOLS == "yes" ]]; then
  echo "$now START processing of UMI-tools count table in R"
  echo "$now Input: $COUNTFILE"
  output_dir=$COUNTS/R_processed/
  make_dir $output_dir
  Rscript $QC_SCRIPT "UMItools" $COUNTFILE $NGENES $NUMIS \
    $MAD $THRESHOLDMT $NORMALIZE $FILTERGENES $output_dir 
  echo "$now Finished."
  echo "$now Output: $output_dir"
fi
if [[ $USESTARSOLO == "yes" ]]; then
  count_dir=${R2_SOLO_OUT}Solo.out/
  output_dir=$STARSOLO_DIR/R_processed/
  make_dir $output_dir
  echo "$now START processing of STARsolo count table in R"
  echo "$now Input: $count_dir"
  Rscript $QC_SCRIPT "STARsolo" $count_dir $NGENES $NUMIS \
    $MAD $THRESHOLDMT $NORMALIZE $FILTERGENES $output_dir
  echo "$now Finished."
  echo "$now Output: $output_dir"
fi
if [[ $USECELLRANGER == "yes" ]]; then
  echo "$now START processing of CellRanger count table in R"
  count_dir=$CROUTPUT/$sample/outs/raw_feature_bc_matrix/
  output_dir=$CROUTPUT/$sample/R_processed/
  make_dir $output_dir
  Rscript $QC_SCRIPT "CellRanger" $count_dir $NGENES $NUMIS \
    $MAD $THRESHOLDMT $NORMALIZE $FILTERGENES $output_dir
  echo "$now Finished."
  echo "$now Output: $output_dir"
fi

echo "$now Remove generated FastQ files:"
if [[ $USEUMITOOLS == "yes" || $USESTARSOLO == "yes" ]]; then
  rm -f $R1 $R2
  rm -f $R1_STAR $R2_STAR
  rmdir -f $MERGED
fi

# remove unzipped annotation file - restore original file structure
if [[ $ANNOZIP == 1 ]]; then
  echo "$now Remove unzipped annotation file"
  rm $ANNOTATION
  echo "$now Finished."
fi
# Remove unzipped whitelist file
if [[ $UMIWHITEZIP == 1 ]]; then
  echo "Remove unzipped whitelist file: $UMITOOLSWHITELIST"
  rm $UMITOOLSWHITELIST
  echo "$now Finished."
fi
if [[ $STARWHITEZIP == 1 ]]; then
  echo "Remove unzipped whitelist file: $STARWHITELIST"
  rm $STARWHITELIST
  echo "$now Finished."
fi

echo "Compress previous compressed files"
for file in ${GZCOMPRESSED[@]}; do
  echo "Remove unipped file: $file"
  rm $file
  echo "$now Finished."
done

for file in ${BZ2COMPRESSED[@]}; do
  echo "$now Remove unzipped file: $file"
  rm $file
  echo "$now Finished."
done

echo "$now End of Pipeline"

