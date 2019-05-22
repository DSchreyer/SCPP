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

# TODO: Demultiplexing after Bam file generation
# TODO: change sequencing file input -> run pipeline only one time
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

# create directories for essential outputs  # print directory paths
TRIMDIR=$OUTPUT/trimmomatic_output 
FASTQCDIR=$OUTPUT/FastQC_output    
STARDIR=$OUTPUT/STAR_output    
RSEMOUT=$OUTPUT/RSEM_output    
BAMSINGLE=$STARDIR/SingleEnd   
BAMPAIRED=$STARDIR/PairedEnd   
MULTIQCDIR=$OUTPUT/MultiQC_output  
TRIMSEDIR="$TRIMDIR/SingleEnd" 
TRIMPEDIR="$TRIMDIR/PairedEnd" 
STARPAIREDDIR=${STARDIR}/PairedEnd 
STARSINGLEDIR=${STARDIR}/SingleEnd 
UMITOOLSDIR=${OUTPUT}/Umi-Tools

make_dir $OUTPUT 
make_dir $TRIMDIR 
make_dir $FASTQCDIR   
make_dir "${TRIMDIR}/PairedEnd"   
make_dir "${TRIMDIR}/SingleEnd"   
make_dir $RSEMOUT 
make_dir $OUTPUT  
make_dir $BAMPAIRED   
make_dir $BAMSINGLE   
make_dir $INDICESDIR  
make_dir $STARDIR 
make_dir "${STARDIR}/SingleEnd"   
make_dir "${STARDIR}/PairedEnd"   
make_dir $RSEMREFDIR
make_dir $UMITOOLSDIR 


# print directory paths
echo Umi-Tools path = ${UMITOOLS}
echo Umi-Tools dir
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
echo Sequencing Read = $READ
echo Barcode and Umi Read = $BARCODE

FILES=$(ls -d $DATA/*)
date
# Quality Control with all files in the data Directory --- FastQC
echo "Perform quality control"
echo "Quality Control of ${FILES[@]}"
# $FASTQC -o $FASTQCDIR -t $THREADS ${FILES[@]}
echo "Performed quality control"

# Umi-Tools
# Identify correct cell barcodes

EXTRACTED_UMIS=()
DEMUX_FILES=()
for file in ${FILES[@]}; do
  echo "Umi-Tools: $file"
  if [[ $file =~ ^(.*)/(.*)_${BARCODE}_001\.(fastq|fq)(\.gz|\.bz2)* ]]; then
    echo "Using Umi-Tools to extract Barcode"
    date 
    dir=${BASH_REMATCH[1]}
    SAMPLENAME=${BASH_REMATCH[2]}
    FORMAT=${BASH_REMATCH[3]}
    ZIP=${BASH_REMATCH[4]}
    EXTRACTED=$UMITOOLSDIR/${SAMPLENAME}_${BARCODE}_extracted_001.${FORMAT}${ZIP}
    WHITELIST="${SAMPLENAME}.whitelist.txt"
    READ2=$(ls $dir/${SAMPLENAME}_${READ}_001.*)
    READ2EXTRACTED=$UMITOOLSDIR/${SAMPLENAME}_${READ}_extracted_001.fastq.gz
    if [[ ! -f $READ2 ]]; then
      echo "Barcoding file: $file has no sequencing file $READ2!"
      echo "Error: $READ2 does not exist!"
      date
      exit
    fi
    # Identify correct cell barcodes
    echo "Identify correct cell barcodes with Umi-Tools!"
    echo "START Umi-Tools whitelist: $file"
    date
    # $UMITOOLS whitelist \
    #   --stdin $file \
    #   --bc-pattern='(?P<cell_1>.{16})(?P<umi_1>.{10})' \
    #   --extract-method=regex \
    #   --log2stderr > $WHITELIST
    echo "END Umi-Tools whitelist: $file" 
    date

    echo "Extract barcodes and UMIs and add to read names"
    echo "START Umi-Tools extract: $file"
    date
    $UMITOOLS extract \
      --stdin $file \
      --bc-pattern='(?P<cell_1>.{16})(?P<umi_1>.{10})' \
      --stdout $EXTRACTED \
      --read2-in $READ2 \
      --read2-out $READ2EXTRACTED \
      --extract-method=regex \
      --quality-filter-mask=20 \
      --quality-encoding="phred33"

    # --filter-cell-barcode \ # Retain Barcodes in whitelist
    # --whitelist=$WHITELIST
    echo "END Umi-Tools extract: $file"
    date
    DEMUX_FILES+=($READ2EXTRACTED)
  fi
done

# Trimming 10x sequencing read2
SE=()
for file in ${DEMUX_FILES[@]}; do
  echo file = $file
  file_exists $file
  if [[ $file =~ ^.*/(.*)_${READ}_extracted_001\.(fastq|fq)(\.gz|\.bz2)* ]]; then
    echo "$file is sequencing read!"
    SAMPLENAME=${BASH_REMATCH[1]}
    FORMAT=${BASH_REMATCH[2]}
    ZIP=${BASH_REMATCH[3]}
    FILTERED="${TRIMSEDIR}/${SAMPLENAME}_extracted_trim_${READ}_001.${FORMAT}${ZIP}"
    echo "Perform trimming with $file"
    date
    java -jar $TRIMMOMATIC \
      SE -phred33 $file $FILTERED \
      -threads ${THREADS} \
      LEADING:20 TRAILING:20 MINLEN:50 SLIDINGWINDOW:4:20
    echo "Trimmed filed saved as $FILTERED"
    date
    SE+=($FILTERED) 
  else
    echo "$file was not trimmed!"
  fi
done

# Quality Control of trimmed single end files with FastQC
echo "Quality Control of trimmed single end files"
$FASTQC -o $FASTQCDIR -t $THREADS ${SE[@]}
echo "SINGLE END: ${SE[@]}"
echo "PAIRED END: ${PE[@]}"
echo "Finished Trimming"
date

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
  echo  "gunzip annotation file ${ANNOTATION}"
  gunzip --keep "${ANNOTATION}"
  ANNOTATION=$(echo "${ANNOTATION}" | sed 's/.gz$//')
  ANNOZIP=1
  echo "Finished unzipping annotation file"
  date
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
  echo "Finished generating STAR indices!"
  date
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
      STAROUT=$STARPAIREDDIR/$(basename $SAMPLENAME)
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
        --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM

    elif [[ $ZIP == ".gz" ]]; then
      echo "gzipped"
      $STAR --runThreadN $THREADS \
        --genomeDir "${INDICESDIR}" \
        --readFilesIn $READ1 $READ2 \
        --readFilesCommand gunzip -c \
        --outFileNamePrefix $STAROUT \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM
    elif [[ $ZIP == "" ]]; then
      echo "unzip"
      $STAR --runThreadN $THREADS \
        --genomeDir "${INDICESDIR}" \
        --readFilesIn $READ1 $READ2 \
        --outFileNamePrefix $STAROUT \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM
    else
      echo "$READ1 and/or $READ2 has the wrong format"
    fi
  done

  SE_STAR=()
  for file in "${SE[@]}"; do
    echo $file
    file_exists $file
    if [[ $file =~ ^.*/(.*)\.(fastq|fq)(\.gz|\.bz2)?$ ]]; then
      echo "Start STAR alignment with $file!"
      SAMPLENAME=${BASH_REMATCH[1]}
      FORMAT=${BASH_REMATCH[2]}
      ZIP=${BASH_REMATCH[3]}
      echo SampleName = $SAMPLENAME
      echo Format = $FORMAT
      echo ZIP = $ZIP
      STAROUT=$STARSINGLEDIR/${SAMPLENAME}
      echo "STAR Output = $STAROUT"
      SE_STAR+=($STAROUT)
    else
      echo "$file has the wrong format or does not exist"
      exit
    fi
    if [[ $ZIP == ".bz2" ]]; then
      echo "bzipped"
      $STAR --runThreadN $THREADS \
        --genomeDir "${INDICESDIR}" \
        --readFilesIn $file \
        --readFilesCommand bunzip2 -c \
        --outSAMtype BAM \
        --outFileNamePrefix $STAROUT \
        --quantMode TranscriptomeSAM 
      echo "Saved STAR output of $file in $STAROUT"
    elif [[ $ZIP == ".gz" ]]; then
      echo "gzipped"
      $STAR --runThreadN $THREADS \
        --genomeDir "${INDICESDIR}" \
        --readFilesIn $file \
        --readFilesCommand gunzip -c \
        --outSAMtype BAM \
        --outFileNamePrefix $STAROUT \
        --quantMode TranscriptomeSAM
      echo "Saved STAR output of $file in $STAROUT"
    elif [[ $ZIP == "" ]]; then
      echo "unzip"
      $STAR --runThreadN $THREADS \
        --genomeDir "${INDICESDIR}" \
        --readFilesIn $file \
        --outSAMtype BAM \
        --outFileNamePrefix $STAROUT \
        --quantMode TranscriptomeSAM
      echo "Saved STAR output of $file in $STAROUT"
    else
      echo "$file has the wrong format"
      help_message
    fi
  done
else
  echo "Skipped STAR Alignment"
fi
echo "Finished STAR Alignment!"
date
echo "Pipeline done!"
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
  echo "Finished generating RSEM reference!"
  date
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
  RSEMINPUT=${file}.Aligned.toTranscriptome.out.bam
  file_exists $RSEMINPUT
  if [[ -f ${RSEMINPUT} ]]; then
    RSEMOUTPUTFILE=${RSEMOUT}/$(basename $file)
    echo "Quantification with RSEM: ${RSEMINPUT}"
    ${RSEM}/rsem-calculate-expression \
      --paired-end \
      --quite \
      -p $THREADS \
      --alignments ${RSEMINPUT} \
      ${PREP_REF} \
      $RSEMOUTPUTFILE
    $SAMTOOLS
    file_exists $RSEMOUTPUTFILE
  else
    echo "${RSEMINPUT} does not exist. Please check STAR INPUT and STAR OUTPUT"
    exit
  fi
done

RSEMSORTED=()
for file in "${SE_STAR[@]}"; do
  RSEMINPUT=${file}.Aligned.toTranscriptome.out.bam
  file_exists $RSEMINPUT
  echo "Quantification with RSEM: ${RSEMINPUT}"
  RSEMOUTPUTFILE=${RSEMOUT}/$(basename $file)
  ${RSEM}/rsem-calculate-expression \
    -p $THREADS \
    --quite \
    --alignments ${RSEMINPUT} \
    ${PREP_REF} \
    $RSEMOUTPUTFILE
  echo "Finished: RSEM rsem-calculate-expression $RSEMINPUT"
  date
  SAMINPUT=${RSEMOUTPUTFILE}.transcript.bam
  SAMOUTPUT=${RSEMOUTPUTFILE}.transcript.sorted.bam
  echo "Start: Samtools sort ${SAMINPUT}"
  $SAMTOOLS sort --threads ${THREADS} ${SAMINPUT} \
    -o ${SAMOUTPUT}
  echo "Finished: Samtools sort ${SAMOUTPUT}"
  date
  RSEMSORTED+=(${SAMOUTPUT})
  file_exists ${SAMOUTPUT}
done

echo "Start counting genes with RSEM bam output!"
echo "Using Umi-tools count"
date

for file in ${RSEMSORTED[@]}; do
  FILE=$(basename $file .bam)
  echo "Umi-tools count: ${FILE}"
  umi_tools count --per-gene \
    --gene-tag=XT --assigned-status-tag=XS \
    --per-cell -I $FILE -S ${file}.tsv.gz
  echo "Umi-tools count DONE: ${FILE}"
done
# cells.counts.txt, cells.metadata.txt, gene.information.txt

echo "Finished generating count table with Umi-tools count"
date
exit

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


