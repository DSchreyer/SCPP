#!/bin/bash
#
# Daniel Schreyer
# Single Cell Analysis Pipeline
# Input: Single cell data
# Output:   1. Cell count table
#           2. Metadata file
#           3. Differentially expressed genes
#           4. PCA, tSNE
#           5. Marker genes for sequenced cell types


# TODO: Imputation R code
# TODO: Demultiplexing integration
# TODO: change sequencing file input -> run pipeline only one time
# TODO: Start R-script from pipeline

while [[ $# -gt 0 ]]
do
  option="$1"
  case $option in
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
    --threads) # how many threads are available
      THREADS="$2"
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
    --STAR)
      STAR="$2"
      shift
      ;;
    --RSEM)
      RSEM="$2"
      shift
      ;;
    --rsemResult) # use "gene" or "isoform" counts for downstream analyses
      RSEMRESULT="$2" # gene counts = "genes" --- isoform counts ="isoforms"
      shift
      ;;
    --rsemRef)
      RSEMREF="$2" # "no": create rsemref in RSEMREFDIR,
                   # "yes": specify ref path and prefix in --rsemRefDir
      shift
      ;;
    --rsemRefDir)
      RSEMREFDIR="$2" # directory for Rsem Reference [Default: "$RSEMOUT/ref"]
      shift
      ;;
    --bamFiles)
      BAMFILES="$2" # "no": Alignment, "yes": specify bam file directory
      shift
      ;;
    --bamSuffix)
      BAMSUFFIX="$2" # use this bam suffix ["Aligned.toTranscriptome.out.bam"]
      shift
      ;;
    --impute)
      IMPUTE="$2" # "no": no imputation, "yes": imputation 
      shift
      ;;
    --animal) # Metadata information: used animal for scRNA-seq
      ANIMAL="$2"
      shift
      ;;
    --sex) # Metadata information: sex of the animal
      SEX="$2" 
      shift
      ;;
    --condition) # Metadata information: genetic or physical condition
      CONDITION="$2"
      shift
      ;;
    --treatment) # Metadata information: control or treatment groups?
      TREATMENT="$2"
      shift
      ;;
    --output) # Output directory
      OUTPUT="$2"
      shift
      ;;
    *)
      echo -e "ERROR: \"${option}\" is an unknown option!"
      help_message
      ;;
  esac
  shift
done

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

# check, if all arguments were passed to the Pipeline
EXPECTEDARGS=(PROJECT GENOME ANNOTATION INDICESDIR GENOMEINDEX DATA THREADS \
    TRIMMOMATIC FASTQC STAR RSEM RSEMRESULT RSEMREF BAMFILES IMPUTE ANIMAL SEX \
    CONDITION TREATMENT OUTPUT)
ARGARRAY=($PROJECT $GENOME $ANNOTATION $INDICESDIR $GENOMEINDEX $DATA $THREADS \
    $TRIMMOMATIC $FASTQC $STAR $RSEM $RSEMRESULT $RSEMREF $BAMFILES $IMPUTE \
    $ANIMAL $SEX $CONDITION $TREATMENT $OUTPUT)
    
if [[ ${#EXPECTEDARGS[@]} != ${#ARGARRAY[@]} ]]; then
  echo "One or multiple argument/s is/are missing!"
  help_message
fi

# print directory paths
echo STAR Directory = "${STARDIR}"
echo FASTQC PATH = "${FASTQC}"
echo Trimmomatic Path = "${TRIMMOMATIC}"
echo STAR Path = "${STAR}"
echo RSEM Path = "${RSEM}"
echo RSEM REF Prepaired = "${RSEMREF}"
echo Prepared RSEM reference in "${RSEMREFDIR}"
echo Use calculated expression of "${RSEMRESULT}"
echo BAM Files = "$BAMFILES"
echo Directory for paired end BAM files = "$BAMPAIRED"
echo Directory for single end BAM files = "$BAMSINGLE"

# Quality Control with all files in the data Directory --- FastQC
FILES=($(ls -d $DATA/*))
echo "Quality Control of all files in ${DATA}"
$FASTQC -o $FASTQCDIR -t $THREADS ${FILES[@]}

# skip trimming if trim == 1
trim=1
if [[ $trim ==  1 ]]; then
  for file in $(ls $DATA); do
    FILE="$DATA/$file"
    if [[ $FILE =~ .*\.fq.*|.*\.fastq.* ]]; then
      if [[ $FILE =~ ^.*/(.*)(1)(\.f[a-z]*.*)$|^.*/(.*)([^2])(\.f[a-z]*.*)$ ]]
      then
        SAMPLENAME=${BASH_REMATCH[1]}
        READ=${BASH_REMATCH[2]}
        FORMAT=${BASH_REMATCH[3]}
        echo Sample Name: $SAMPLENAME
        if [[ ! ${READ} == "1" ]]; then # if there is no read2
          echo "${FILE} has no pair --> Single End Trimming"
          echo "Single End Trimming: ${FILE}"
          FILTERED="$TRIMSEDIR/${SAMPLENAME}Filtered${FORMAT}"

          java -jar $TRIMMOMATIC \
            SE -phred33 $FILE $FILTERED \
            -threads ${THREADS} \
            LEADING:20 TRAILING:20 MINLEN:60
          SE+=($FILTERED)
        else
          FILE2="${DATA}/${SAMPLENAME}2${FORMAT}"
          if [[ -f $FILE2 ]]; then
            echo "PairedEnd Trimming: ${FILE} + ${FILE2}"

            java -jar $TRIMMOMATIC \
              PE -phred33 $FILE $FILE2 \
              -baseout "$TRIMPEDIR/${SAMPLENAME}${FORMAT}" \
              -threads ${THREADS} \
              LEADING:20 TRAILING:20 MINLEN:60

            PE+=($TRIMPEDIR/${SAMPLENAME}${FORMAT})
            # Quality Control with FastQC
            $FASTQC -o $FASTQCDIR \
                -t $THREADS \
                ${TRIMPEDIR}/${SAMPLENAME}_[12]P${FORMAT}
          else
            echo "${FILE2} does not exist or has a different compression state."
            echo "${FILE} is Single End"
            echo "Single End Trimming: ${FILE}"

            FILTERED="$TRIMSEDIR/${SAMPLENAME}Filtered${FORMAT}"
            java -jar $TRIMMOMATIC \
              SE -phred33 $FILE $FILTERED \
              -threads ${THREADS} \
              LEADING:20 TRAILING:20 MINLEN:60

            SE+=(${FILTERED})
          fi
        fi
      fi
    fi
  done
else
  echo "Skipped trimming with Trimmomatic"
fi

# Quality Control of trimmed single end files with FastQC
echo "Quality Control of trimmed single end files"
$FASTQC -o $FASTQCDIR -t $THREADS ${SE[@]}

echo "SINGLE END: ${SE[@]}"
echo "PAIRED END: ${PE[@]}"

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
  echo  "gunzip annotation file ${ANNOTATION}"
  gunzip --keep "${ANNOTATION}"
  ANNOTATION=$(echo "${ANNOTATION}" | sed 's/.gz$//')
  ANNOZIP=1
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
      --runThreadN ${THREADS} \
      --quantMode TranscriptomeSAM
  else
    if [[ $ANNOTATION =~ .*\.gtf ]]; then
      echo "Run genomeGenerate with gtf annotation file"
      $STAR --runMode genomeGenerate \
        --genomeDir "${INDICESDIR}" \
        --sjdbGTFfile "${ANNOTATION}" \
        --genomeFastaFiles "${GENOME}" \
        --runThreadN ${THREADS} \
        --quantMode TranscriptomeSAM
    elif [[ $ANNOTATION =~ .*\.gff ]]; then
      echo "Run genomeGenerate with gff annotation file"
      $STAR --runmode genomeGenerate \
        --genomeDir "${INDICESDIR}" \
        --sjdbGTFtagExonParentTranscript Parent \
        --genomeFastaFiles "${GENOME}" \
        --runThreadN "${THREADS}" \
        --quantMode TransriptomeSAM 
    else
      echo "Annotation file ${ANNOTATION} has the wrong format"
      echo "Enter .gtf or .gtf annotation file"
      help_message
    fi
  fi
  echo "Stored genome indices in ${INDICESDIR}"
fi

# skip star alignment
star=0
if [[ $star == 0 ]]; then
  # align only paired reads with STAR -- 1P and 2P extension
  PE_STAR=()
  for file in "${PE[@]}";do
    echo "Start aligning $FILE with STAR"
    if [[ $file =~ ^(.*)(\.f[a-z]*)(\.gz|\.bz2)?$ ]]; then
      SAMPLENAME=${BASH_REMATCH[1]}
      FORMAT=${BASH_REMATCH[2]}
      ZIP=${BASH_REMATCH[3]}
      STAROUT=$STARPAIREDDIR/$(basename $SAMPLENAME)_
      PE_STAR+=($STAROUT)
      echo SAMPLE NAME = $SAMPLENAME
      echo FORMAT = $FORMAT
      echo ZIP = $ZIP
    fi
    READ1=${SAMPLENAME}_1P${FORMAT}${ZIP}
    READ2=${SAMPLENAME}_2P${FORMAT}${ZIP}

    if [[ $ZIP == ".bz2" ]]; then
      echo "bzipped"
      $STAR --runThreadN $THREADS \
        --genomeDir "${INDICESDIR}" \
        --readFilesIn $READ1 $READ2 \
        --readFilesCommand bunzip2 -c \
        --outFileNamePrefix $STAROUT \
        --quantMode TranscriptomeSAM

    elif [[ $ZIP == ".gz" ]]; then
      echo "gzipped"
      $STAR --runThreadN $THREADS \
        --genomeDir "${INDICESDIR}" \
        --readFilesIn $READ1 $READ2 \
        --readFilesCommand gunzip -c \
        --outFileNamePrefix $STAROUT \
        --quantMode TranscriptomeSAM
    elif [[ -z "$ZIP" ]]; then
      echo "unzip"
      $STAR --runThreadN $THREADS \
        --genomeDir "${INDICESDIR}" \
        --readFilesIn $READ1 $READ2 \
        --outFileNamePrefix $STAROUT \
        --quantMode TranscriptomeSAM
    else
      echo "$READ1 and/or $READ2 has the wrong format"
    fi
  done

  SE_STAR=()
  for file in "${SE[@]}"; do
    echo $file
    if [[ $file =~ ^(.*)Filtered(\.f[a-z]*)(\.gz|\.bz2)?$ ]]; then
      echo "Start Aligning $file to $(basename $GENOMEINDEX)"
      SAMPLENAME=${BASH_REMATCH[1]}
      FORMAT=${BASH_REMATCH[2]}
      ZIP=${BASH_REMATCH[3]}
      STAROUT=$STARSINGLEDIR/$(basename $SAMPLENAME)_
      SE_STAR+=($STAROUT)
      echo SAMPLE NAME = $SAMPLENAME
      echo FORMAT = $FORMAT
      echo ZIP = $ZIP
    else
      echo "$FILE has the wrong format or does not exist"
      exit
    fi
    if [[ $ZIP == ".bz2" ]]; then
      echo "bzipped"
      $STAR --runThreadN $THREADS \
        --genomeDir "${INDICESDIR}" \
        --readFilesIn $file \
        --readFilesCommand bunzip2 -c \
        --outFileNamePrefix $STAROUT \
        --quantMode TranscriptomeSAM
    elif [[ $ZIP == ".gz" ]]; then
      echo "gzipped"
      $STAR --runThreadN $THREADS \
        --genomeDir "${INDICESDIR}" \
        --readFilesIn $file \
        --readFilesCommand gunzip -c \
        --outFileNamePrefix $STAROUT \
        --quantMode TranscriptomeSAM
    elif [[ -z "$ZIP" ]]; then
      echo "unzip"
      $STAR --runThreadN $THREADS \
        --genomeDir "${INDICESDIR}" \
        --readFilesIn $file \
        --outFileNamePrefix $STAROUT \
        --quantMode TranscriptomeSAM
    else
      echo "$file has the wrong format"
      help_message
    fi
  done
else
  echo "Skipped STAR Alignment"
fi
echo "Quantification with RSEM"

# create RSEM reference
if [[ $RSEMREF == "no" ]]; then

  echo "Prepare Reference with RSEM"
  echo $ANNOTATION
  if [[ $ANNOTATION =~ ^.*/(.*)(\..*)$ ]]; then
    ANNONAME=${BASH_REMATCH[1]}
    FORMAT=${BASH_REMATCH[2]}
    PREP_REF=$RSEMREFDIR/$ANNONAME

    if [[ $FORMAT == ".gtf" ]]; then
      ${RSEM}rsem-prepare-reference --gtf ${ANNOTATION} \
        ${GENOME} \
        ${PREP_REF}
    elif [[ $FORMAT == ".gff" ]]; then
      ${RSEM}rsem-prepare-reference --gff ${ANNOTATION} \
        ${GENOME} \
        ${PREP_REF}
    else
      echo "$ANNOTATION has the wrong format. GFF and GTF Format accepted"
      help_message
    fi
    echo "Prepare Reference with RSEM: DONE!"
    echo "Reference Files are stored in ${PREP_REF}"
  fi
else
  echo "Skipped preparing reference with RSEM"
  PREP_REF=${RSEMREFDIR}
fi

# remove unzipped annotation file - restore original file structure
if [[ $ANNOZIP == 1 ]]; then
  echo "Remove unzipped annotation file"
  rm $ANNOTATION
fi

echo "PAIRED END FILES: ${PE_STAR[@]}"
echo "SINGLE END FILES: ${SE_STAR[@]}"

# Quantification step with RSEM function: rsem-calculate-expression
RSEMOUTPUTFILES=()
for file in "${PE_STAR[@]}"; do
  echo $file
  RSEMINPUT=${file}Aligned.toTranscriptome.out.bam
  if [[ -f ${RSEMINPUT} ]]; then
    RSEMOUTPUTFILES+=(${RSEMOUT}/$(basename $file))
    echo "Quantification with RSEM: ${RSEMINPUT}"
    ${RSEM}/rsem-calculate-expression \
      --no-bam-output \
      --paired-end \
      -p $THREADS \
      --alignments ${RSEMINPUT} \
      ${PREP_REF} \
      ${RSEMOUT}/$(basename $file)
  else
    echo "${RSEMINPUT} does not exist. Please check STAR INPUT and STAR OUTPUT"
    exit
  fi
done

for file in "${SE_STAR[@]}"; do
  RSEMINPUT=${file}Aligned.toTranscriptome.out.bam
  if [[ -f ${RSEMINPUT} ]]; then
    echo "Quantification with RSEM: ${RSEMINPUT}"
    RSEMOUTPUTFILES+=(${RSEMOUT}/$(basename $file))
    ${RSEM}/rsem-calculate-expression \
      --no-bam-output \
      -p $THREADS \
      --alignments ${RSEMINPUT} \
      ${PREP_REF} \
      ${RSEMOUT}/$(basename $file)
  else
    echo "${RSEMINPUT} does not exist. Please check STAR INPUT and STAR OUTPUT"
    exit
  fi
done

# cells.counts.txt, cells.metadata.txt, gene.information.txt

echo RSEM output files = "${RSEMOUTPUTFILES[@]}"

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
echo "Generated count table. Stored in ${CELLCOUNTS}."
echo "Generated metadata file. Stored in ${METADATA}."
echo "Start Quality Control and Normalization with scater and scran in R!"
Rscript R_pipeline.R $CELLCOUNTS $METADATA


