#!/bin/bash
while [[ $# -gt 0 ]]
do
  option="$1"

  case $option in
      --projectName)
      PROJECT="$2"
      shift
      ;;
      --genome)
      GENOME="$2"
      shift # past argument
      ;;
      --annotation)
      ANNOTATION="$2"
      shift # past argument
      ;;
      --indicesDir)
      INDICESDIR="$2" # path to directory where genome indices are stored
      shift
      ;;
      --index) # y/yes if indices is available n/no if indices need to be generated
      GENOMEINDEX="$2"
      shift
      ;;
      --data)  # path to data directory with read files
      DATA="$2"
      shift
      ;;
      --threads)
      THREADS="$2"
      shift
      ;;
      --trimmomatic)
      TRIMMOMATIC="$2"
      shift
      ;;
      --fastqc)
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
      --rsemResult)
      RSEMRESULT="$2" # use "genes"/"isoforms"
      shift
      ;;
      --rsemRef)
      RSEMREF="$2" # "no" -> prepare reference else skipped, if yes
      # specify rsemRefDir and Prefix at --rsemRefDir
      shift
      ;;
      --rsemRefDir)
      RSEMREFDIR="$2" # specifiy directory for Rsem Reference [Default: "$RSEMDIR/ref"]
      shift
      ;;
      --bamFiles)
      BAMFILES="$2" # "no" -> do alignment, yes -> specifiy bam dir for paired end reads and single end
      shift
      ;;
      --bamSuffix)
      BAMSUFFIX="$2" # Use Files with this suffix [e.g. "Aligned.toTranscriptome.out.bam"]
      shift
      ;;
      --impute)
      IMPUTE="$2" # Yes, if cell counts should be imputed to identify real dropouts
      shift
      ;;
      --animal) # used animal for scRNA-seq
      ANIMAL="$2"
      shift
      ;;
      --sex)
      SEX="$2" # sex of the animal
      shift
      ;;
      --condition) # genetic or physical condition of animal used for scRNA-seq
      CONDITION="$2"
      shift
      ;;
      --treatment) # control or treatment groups?
      TREATMENT="$2"
      shift
      ;;
      --output) # Store output files in this directory 
      OUTPUT="$2"
      shift
      ;;
      *)
      echo -e "ERROR: \"${option}\" is an unknown option!"
      exit
      ;;
  esac
  shift
done

make_dir (){
    local dir="${1}"
    if [[ ! -d $dir ]]; then
        mkdir -p $dir
    fi
}

if [[ $METADATA == "" ]]; then
  METADATA="./cells.metadata.txt"
fi

if [[ ! $RSEMRESULT =~ ^(genes|isoforms)$ ]]; then
  echo "Please enter the RSEM output, which should be used for further analyses!"
  echo "Enter 'genes' or 'isoforms'."
  exit
fi


if [[ ! -f $GENOME || $GENOME =~ [^.*fa(.gz)?]$ ]]
then
  echo -e "${GENOME} is not a file or has the wrong format.
  Please enter a reference genome in fasta format.\n"
fi

# create directories for essential outputs

TRIMDIR=$OUTPUT/trimmomatic_output
FASTQCDIR=$OUTPUT/FastQC_output
STARDIR=$OUTPUT/STAR_output
RSEMDIR=$OUTPUT/RSEM_output
BAMSINGLE=$STARDIR/SingleEnd
BAMPAIRED=$STARDIR/PairedEnd
MULTIQCDIR=$OUTPUT/MultiQC_output
SEDIR="$TRIMDIR/SingleEnd"
PEDIR="$TRIMDIR/PairedEnd"
STARPAIREDDIR=${STARDIR}/PairedEnd
STARSINGLEDIR=${STARDIR}/SingleEnd

make_dir $OUTPUT
make_dir $TRIMDIR
make_dir $FASTQCDIR
make_dir "${TRIMDIR}/PairedEnd"
make_dir "${TRIMDIR}/SingleEnd"
make_dir $RSEMDIR
make_dir $OUTPUT
make_dir $BAMPAIRED
make_dir $BAMSINGLE
make_dir $INDICESDIR
make_dir $STARDIR
make_dir "${STARDIR}/SingleEnd"
make_dir "${STARDIR}/PairedEnd"
make_dir $RSEMREFDIR

echo Project Name = "${PROJECT}"
echo Reference Genome  = "${GENOME}"
echo Annotation File = "${ANNOTATION}"
echo Data Directory = "${DATA}"
echo Trimmomatic Output Directory = "${TRIMDIR}"
echo Genome Index available = "${GENOMEINDEX}"
echo Indices Directory = "${INDICESDIR}"
echo Number of Threads = "${THREADS}"
echo RSEM Directory = "$RSEMDIR"
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

for FILE in $(ls $DATA); do
  FILE="$DATA/$FILE"
  if [[ $FILE =~ .*\.fq.*|.*\.fastq.* ]]; then
    if [[ $FILE =~ ^.*/(.*)(1)(\.f[a-z]*.*)$|^.*/(.*)([^2])(\.f[a-z]*.*)$ ]]
    then
      SAMPLENAME=${BASH_REMATCH[1]}
      READ=${BASH_REMATCH[2]}
      FORMAT=${BASH_REMATCH[3]}
      echo Sample Name: $SAMPLENAME
      if [ ! ${READ} == "1" ]; then # if there is no read2
        echo "${FILE} has no pair --> Single End Trimming"
        echo "Single End Trimming: ${FILE}"
        FILTERED="$SEDIR/${SAMPLENAME}Filtered${FORMAT}"

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
          PE -phred33 $FILE $FILE2 -baseout "$PEDIR/${SAMPLENAME}${FORMAT}" \
          -threads ${THREADS} \
          LEADING:20 TRAILING:20 MINLEN:60

          PE+=($PEDIR/${SAMPLENAME}${FORMAT})
          # Quality Control with FastQC
          $FASTQC -o $FASTQCDIR -t $THREADS ${PEDIR}/${SAMPLENAME}_[12]P${FORMAT}
        else
          echo "${FILE2} does not exist or has a different compression state"
          echo "${FILE} is Single End"
          echo "Single End Trimming: ${FILE}"

          FILTERED="$SEDIR/${SAMPLENAME}Filtered${FORMAT}"
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
  echo "Skip Trimming with Trimmomatic"
fi

# Quality Control of trimmed single end files with FastQC
echo "Quality Control of trimmed single end files"
$FASTQC -o $FASTQCDIR -t $THREADS ${SE[@]}

echo "SINGLE END: ${SE[@]}"
echo "PAIRED END: ${PE[@]}"


if [[ $GENOME =~ .*fa.gz ]]
then
	echo "unzip genome fasta file"
	gunzip --keep "${GENOME}"
	GENOME=$(echo "${GENOME}" | sed 's/.gz$//')
fi

ANNOZIP=0
# gunzip gzipped annotation file
if [[ $ANNOTATION =~ .*gtf.gz ]]; then
	echo  "unzip annotation file ${ANNOTATION}"
	gunzip --keep "${ANNOTATION}"
	ANNOTATION=$(echo "${ANNOTATION}" | sed 's/.gz$//')
    ANNOZIP=1
fi


echo Annotation File = "${ANNOTATION}"
echo "Generate Genome Index"

echo "${GENOMEINDEX}"
if [[ ! ${GENOMEINDEX} =~ yes|y|n|no|Yes|No|Y|N|NO|YES ]]; then
  echo "Please specify, if the STAR indices are available or need to be generated"
  exit
fi

if [ -z ${GENOMEINDEX+x} ] || [[ ${GENOMEINDEX} =~ no|n|No|N|NO ]]; then
  if [[ -z ${ANNOTATION+x} ]]; then
    echo "No annotation file. Indexing with an annotation file is higly recommended"
    $STAR --runMode genomeGenerate \
    --genomeDir "${INDICESDIR}" \
    --genomeFastaFiles "${GENOME}" \
    --runThreadN ${THREADS} \
    --quantMode TranscriptomeSAM
  else
    echo "Run genomeGenerate with annotation file"
    $STAR --runMode genomeGenerate \
    --genomeDir "${INDICESDIR}" \
    --sjdbGTFfile "${ANNOTATION}" \
    --genomeFastaFiles "${GENOME}" \
    --runThreadN ${THREADS} \
    --quantMode TranscriptomeSAM
  fi
  echo "Stored genome indices in ${INDICESDIR}"
fi


# skip star alignment
star=0
if [[ $star == 0 ]]; then
# align only paired reads with STAR -- 1P and 2P extension
PE_STAR=()
for FILE in "${PE[@]}";do
  echo "Start aligning $FILE with STAR"
  if [[ $FILE =~ ^(.*)(\.f[a-z]*)(\.gz|\.bz2)?$ ]]; then
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
  elif [[ $ZIP == "" ]]; then
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
for FILE in "${SE[@]}"; do
  echo $FILE
  if [[ $FILE =~ ^(.*)Filtered(\.f[a-z]*)(\.gz|\.bz2)?$ ]]; then
    echo "Start Aligning $FILE to $(basename $GENOMEINDEX)"
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
    --readFilesIn $FILE \
    --readFilesCommand bunzip2 -c \
    --outFileNamePrefix $STAROUT \
    --quantMode TranscriptomeSAM
  elif [[ $ZIP == ".gz" ]]; then
    echo "gzipped"
    $STAR --runThreadN $THREADS \
    --genomeDir "${INDICESDIR}" \
    --readFilesIn $FILE \
    --readFilesCommand gunzip -c \
    --outFileNamePrefix $STAROUT \
    --quantMode TranscriptomeSAM
  elif [[ $ZIP == "" ]]; then
    echo "unzip"
    $STAR --runThreadN $THREADS \
    --genomeDir "${INDICESDIR}" \
    --readFilesIn $FILE \
    --outFileNamePrefix $STAROUT \
    --quantMode TranscriptomeSAM
  else
    echo "$FILE has the wrong format"
    exit
  fi
done
else
  echo "Skipped STAR Alignment"
fi
echo "Quantification with RSEM"


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
        ${GENOME} \test_datest_da
        ${PREP_REF}
    else
        echo "$ANNOTATION has the wrong format. GFF and GTF Format accepted"
        exit
    fi
    echo "Prepare Reference with RSEM: DONE!"
    echo "Reference Files are stored in ${PREP_REF}"
fi
else
  echo "Skipped preparing reference with RSEM"
  PREP_REF=${RSEMREFDIR}
fi

if [[ $ANNOZIP == 1 ]]; then
    echo "Remove unzipped annotation file"
    rm $ANNOTATION
fi

RSEMOUTPUTDIR=${RSEMDIR/Output}
make_dir RSEMOUTPUTDIR


echo "PAIRED END FILES: ${PE_STAR[@]}"
echo "SINGLE END FILES: ${SE_STAR[@]}"

RSEMOUTPUTFILES=()
for FILE in "${PE_STAR[@]}"; do
    echo $FILE

    RSEMINPUT=${FILE}Aligned.toTranscriptome.out.bam
    if [[ -f ${RSEMINPUT} ]]; then
        RSEMOUTPUTFILES+=(${RSEMOUTPUTDIR}/$(basename $FILE))
        echo "Quantification with RSEM: ${RSEMINPUT}"
        ${RSEM}/rsem-calculate-expression \
        --no-bam-output \
        --paired-end \
        -p $THREADS \
        --alignments ${RSEMINPUT} \
        ${PREP_REF} \
        ${RSEMOUTPUTDIR}/$(basename $FILE)
    else
        echo "${RSEMINPUT} does not exist. Please check STAR INPUT and STAR OUTPUT"
        exit
    fi
done

for FILE in "${SE_STAR[@]}"; do
  RSEMINPUT=${FILE}Aligned.toTranscriptome.out.bam
  if [[ -f ${RSEMINPUT} ]]; then
    echo "Quantification with RSEM: ${RSEMINPUT}"
    RSEMOUTPUTFILES+=(${RSEMOUTPUTDIR}/$(basename $FILE))
    ${RSEM}/rsem-calculate-expression \
    --no-bam-output \
    -p $THREADS \
    --alignments ${RSEMINPUT} \
    ${PREP_REF} \
    ${RSEMOUTPUTDIR}/$(basename $FILE)
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
for FILE in "${RSEMOUTPUTFILES[@]}"; do
  echo "Generate Metadata file for $FILE"
  CELL=$(basename $FILE)
  RSEMOUTPUT="${FILE}.${RSEMRESULT}.results"
  if [[ $i == 0 ]]; then
    if [[ ! -f ${CELLCOUNTS} ]]; then
      awk -v cell="$CELL" 'NR==1 {print $1,cell} NR>1 {print $1,$5}' $RSEMOUTPUT > "${CELLCOUNTS}"
    else
      paste $CELLCOUNTS <(awk -v cell="$CELL" 'NR==1 {print cell} NR>1 {print $5}' $RSEMOUTPUT) > ${TEMP}
      cat ${TEMP} > ${CELLCOUNTS}
      rm ${TEMP}
    fi
    if [[ ! -f $GENEINF ]]; then
      awk '{print $1,$3}' $RSEMOUTPUT > ${GENEINF}
    fi
    if [[ ! -f $METADATA ]]; then
      echo -e "cell\tanimal\tsex\tcondition\ttreatment" > ${METADATA}
    fi
    echo -e "$(basename $FILE)\t${ANIMAL}\t${SEX}\t${CONDITION}\t${TREATMENT}" >> "${METADATA}"
    i=1
  else
    echo -e "$(basename $FILE)\t${ANIMAL}\t${SEX}\t${CONDITION}\t${TREATMENT}" >> ${METADATA}
    paste $CELLCOUNTS <(awk -v cell="$CELL" 'NR==1 {print cell} NR>1 {print $5}' $RSEMOUTPUT) > $TEMP
    cat $TEMP > $CELLCOUNTS
    rm $TEMP
  fi
done

# If imputation is requested: Perform imputation on the raw count table
if [[ $IMPUTE =~ yes|Yes|YES|Y|y ]]; then
    echo "Perform imputation with scImpute before normalization, to reduce dropouts"
    if [[ -z ${IMPUTOUT+x} ]]; then
        DIR=$(dirname "${CELLCOUNTS}")
        BASE=$(basename "${CELLCOUNTS}")
        IMPOUTPUT=$DIR/imputed.$BASE
    fi
    echo "Perform imputation on raw count table ${CELLCOUNTS} and save it in $IMPOUTPUT"
    ./impute_counts.R $CELLCOUNTS $IMPOUTPUT
fi

# shell script finished -- start R Pipeline
# Generating Heatmaps, PCA, tSNE, 
# Find Marker genes
# Do Differential expression analysis 
echo "Generated count table. Stored in ${CELLCOUNTS}."
echo "Generated metadata file. Stored in ${METADATA}."
echo "Start Quality Control and Normalization with scater and scran in R!"
./R_pipeline.R $CELLCOUNTS $METADATA

