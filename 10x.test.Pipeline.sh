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

# TODO: Add FeatureCounts in Pipeline and start.pipeline
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
    --RSEM)
      RSEM="$2"
      shift
      ;;
    --featureCounts)
      FEATURECOUNTS="$2"
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
MULTIQCDIR=$OUTPUT/MultiQC_output  
UMITOOLSDIR=${OUTPUT}/Umi-Tools
COUNTS=${OUTPUT}/counts
MERGED=${OUTPUT}/Files

make_dir $OUTPUT 
make_dir $TRIMDIR 
make_dir $FASTQCDIR   
make_dir $RSEMOUT 
make_dir $OUTPUT  
make_dir $INDICESDIR  
make_dir $STARDIR 
make_dir $RSEMREFDIR
make_dir $UMITOOLSDIR 
make_dir $COUNTS
make_dir $MERGED


# print directory paths
echo Umi-Tools path = ${UMITOOLS}
echo STAR Directory = "${STARDIR}"
echo FASTQC PATH = "${FASTQC}"
echo Trimmomatic Path = "${TRIMMOMATIC}"
echo STAR Path = "${STAR}"
echo RSEM Path = "${RSEM}"
echo RSEM REF Prepaired = "${RSEMREF}"
echo Prepared RSEM reference in "${RSEMREFDIR}"
echo Use calculated expression of "${RSEMRESULT}"
echo Sequencing Read = $READ
echo Barcode and Umi Read = $BARCODE


# Concatenate fastq files together
FILES=$(ls -d $DATA/*)

# READ1 stores each fastq file with barcodes and umis, READ2 contains cDNA read
READ1_ARRAY=()
READ2_ARRAY=()
echo "Start: Concatenating all Fastq files of a sample"
date
for file in ${FILES[@]}; do
  if [[ $file =~ ^(.*)/(.*)_L[0-9]{3}_${BARCODE}_001\.(fastq|fq)(\.gz|\.bz2)* ]]; then
    R1_ARRAY+=($file)
    R1_SAMPLE=${BASH_REMATCH[2]}
    R1_LANES=_ALL_
    R1_FORMAT=${BASH_REMATCH[3]}
    R1_ZIP=${BASH_REMATCH[4]}
    R1=${MERGED}/${R1_SAMPLE}${R2_LANES}${BARCODE}_001.${R1_FORMAT}${R1_ZIP}
  elif [[ $file =~ ^(.*)/(.*)_L[0-9]{3}_${READ}_001\.(fastq|fq)(\.gz|\.bz2)* ]]; then
    R2_ARRAY+=($file)
    R2_SAMPLE=${BASH_REMATCH[2]}
    R2_LANES=_ALL_
    R2_FORMAT=${BASH_REMATCH[3]}
    R2_ZIP=${BASH_REMATCH[4]}
    R2=${MERGED}/${R2_SAMPLE}${R2_LANES}${READ}_001.${R2_FORMAT}${R2_ZIP}
  fi 
done
if [[ ${#R2_ARRAY[@]} -eq 1 ]]; then
  R2=${R2_ARRAY}
  R1=${R1_ARRAY}
  echo "Only one sequencing file"
elif [[ ${#R2_ARRAY[@]} -gt 1 ]]; then
  cat ${R1_ARRAY[@]} > $R1
  cat ${R2_ARRAY[@]} > $R2
  echo "Concatenate ${R2_ARRAY[@]}"
else
  echo "No read file is stored in $DATA!"
  help_message
fi
echo "Finished: Stored Fastq file in $MERGED"


if [[ $QC = "no" ]]; then
  echo "Don't perform Quality Control with FastQC"
elif [[ $QC == "yes" ]]; then
  echo "Perform Quality Control with FastQC"
else
  echo "Perform Quality Control?"
  echo "'yes' | 'no'"
  exit
fi


TRIMMING="yes"
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

# if w != 1 -> no umi_tools whitelist
w=1
if [[ $w == 1 ]]; then
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
    # --quality-filter-mask=20 \
    --quality-encoding="phred33" \
    --filter-cell-barcode \
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
    --quality-filter-mask=20 \
    --quality-encoding="phred33" \
  echo "END Umi-Tools extract: $R1 and $R2"
  echo "Stored: $R1_EXT and $R2_EXT"
  date
fi

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
    LEADING:20 TRAILING:20 HEADCROP:10 MINLEN:75 
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
    --outFileNamePrefix $R2_STAROUT \
    --quantMode TranscriptomeSAM 
  echo "Saved STAR output of $R2_TRIM in $R2_STAROUT"
elif [[ $R2_ZIP == ".gz" ]]; then
  echo "gzipped"
  $STAR --runThreadN $THREADS \
    --genomeDir "${INDICESDIR}" \
    --readFilesIn $R2_TRIM \
    --readFilesCommand gunzip -c \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix $R2_STAROUT \
    --quantMode TranscriptomeSAM
  echo "Saved STAR output of $R2_TRIM in $R2_STAROUT"
elif [[ $R2_ZIP == "" ]]; then
  echo "unzip"
  $STAR --runThreadN $THREADS \
    --genomeDir "${INDICESDIR}" \
    --readFilesIn $R2_TRIM \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix $R2_STAROUT \
    --quantMode TranscriptomeSAM
  echo "Saved STAR output of $R2_TRIM in $R2_STAROUT"
else
  echo "$R2_TRIM has the wrong format"
  help_message
fi
echo "Finished STAR Alignment!"
date

echo "Start Quantification with RSEM"
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

# Quantification step with RSEM function: rsem-calculate-expression
r=1
if [[ $r == 0 ]]; then
  RSEMINPUT="${R2_STAROUT}Aligned.toTranscriptome.out.bam"
  file_exists $RSEMINPUT
  echo "Quantification with RSEM: ${RSEMINPUT}"
  RSEMOUTFILE=${RSEMOUT}/$(basename $RSEMINPUT .out.bam)
  echo "RSEM file prefix: $RSEMOUTFILE"
  ${RSEM}/rsem-calculate-expression \
    -p $THREADS \
    --quiet \
    --alignments ${RSEMINPUT} \
    ${PREP_REF} \
    $RSEMOUTFILE
  echo "Finished: RSEM rsem-calculate-expression $RSEMINPUT"
  date
  TRANSCRIPT_BAM=${RSEMOUTFILE}.transcript.bam
  GENOME_BAM=${RSEMOUTFILE}.genes.bam
  $RSEM/rsem-tbam2gbam $PREP_REF $TRANSCRIPT_BAM $GENOME_BAM
  SAMOUTPUT=${RSEMOUTFILE}.genes.sorted.bam
fi

f=0
# Count reads per gene with featureCounts
if [[ $f == 0 ]]; then
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
fi

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


