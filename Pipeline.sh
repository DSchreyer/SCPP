#!/bin/bash
while [[ $# -gt 0 ]]
do
  option="$1"

  case $option in
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
      --trimDir) # specify trimmomatic output directory
      TRIMDIR="$2"
      shift
      ;;
      --starDir)
      STARDIR="$2"
      shift
      ;;
      --rsemDir)
      RSEMDIR="$2"
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
      --FastQC)
      FASTQC="$2"
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
      --bamPaired)
      BAMPAIRED="$2" # Bam files directory for paired end readss
      shift
      ;;
      --bamSingle)
      BAMSINGLE="$2"
      shift
      ;;
      --bamSuffix)
      BAMSUFFIX="$2" # Use Files with this suffix [e.g. "Aligned.toTranscriptome.out.bam"]
      shift
      ;;
      --metaData) # File Name and Path to store the metadata for each cell
      METADATA="$2" # [Default: cells.metadata.txt]
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
      *)
      echo -e "ERROR: \"${option}\" is an unknown option!"
      exit
      ;;
  esac
  shift
done

if [[ $METADATA == "" ]]; then
  METADATA="./cells.metadata.txt"
fi

if [[ ! $RSEMRESULT =~ ^(genes|isoforms)$ ]]; then
  echo "Please enter the RSEM output, which should be used for further analyses!"
  echo "Enter 'genes' or 'isoforms'."
  exit
fi


echo Reference Genome  = "${GENOME}"
echo Annotation File     = "${ANNOTATION}"
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

if [[ ! -f $GENOME || $GENOME =~ [^.*fa(.gz)?]$ ]]
then
  echo -e "${GENOME} is not a file or has the wrong format.
  Please enter a fasta file.\n"
fi

# create directories for trimmomatic output files
if [[ ! -d "${TRIMDIR}" ]]; then
  mkdir "${TRIMDIR}"
fi

if [[ ! -d "${TRIMDIR}/PairedEnd" ]]; then
  mkdir "${TRIMDIR}/PairedEnd"
fi
if [[ ! -d "${TRIMDIR}/SingleEnd" ]]; then
  mkdir "${TRIMDIR}/SingleEnd"
fi
if [[ ! -d "${RSEMDIR}" ]]; then
  mkdir "${RSEMDIR}"
fi

SEDIR="$TRIMDIR/SingleEnd"
PEDIR="$TRIMDIR/PairedEnd"

# skip trimming with this command
SE=(/media/data/Daniel/test_data/trimmomatic_output//SingleEnd/ERR522934_Filtered.fastq.gz)
PE=(/media/data/Daniel/test_data/trimmomatic_output//PairedEnd/ERR522959_.fastq /media/data/Daniel/test_data/trimmomatic_output//PairedEnd/ERR523111_.fastq.gz)

# skip trimming if trim == 1
trim=0
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

echo "SINGLE END: ${SE[@]}"
echo "PAIRED END: ${PE[@]}"


if [ ! -d "${INDICESDIR}" ]; then
  mkdir ${INDICESDIR}
fi


if [[ $GENOME =~ .*fa.gz ]]
then
	echo "unzip genome fasta file"
	gunzip --keep "${GENOME}"
	GENOME=$(echo "${GENOME}" | sed 's/.gz$//')
fi


# gunzip gzipped annotation file
if [[ $ANNOTATION =~ .*gtf.gz ]]; then
	echo -e "unzip annotation file ${ANNOTATION}\n"
	gunzip --keep "${ANNOTATION}"
	ANNOTATION=$(echo "${ANNOTATION}" | sed 's/.gz$//')
fi


echo -e Annotation File     = "${ANNOTATION}\n"
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

STARPAIREDDIR=${STARDIR}/PairedEnd
STARSINGLEDIR=${STARDIR}/SingleEnd
if [ ! -d "${STARDIR}" ]; then
  mkdir "${STARDIR}"
  mkdir "${STARDIR}/PairedEnd"
  mkdir "${STARDIR}/SingleEnd"
else
  if [ ! -d "${STARPAIREDDIR}" ]; then
    mkdir "${STARPAIREDDIR}"
  fi
  if [ ! -d "${STARSINGLEDIR}" ]; then
    mkdir "${STARSINGLEDIR}"
  fi
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
    OUTPUT=$STARPAIREDDIR/$(basename $SAMPLENAME)_
    PE_STAR+=($OUTPUT)
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
    --outFileNamePrefix $OUTPUT \
    --quantMode TranscriptomeSAM

  elif [[ $ZIP == ".gz" ]]; then
    echo "gzipped"
    $STAR --runThreadN $THREADS \
    --genomeDir "${INDICESDIR}" \
    --readFilesIn $READ1 $READ2 \
    --readFilesCommand gunzip -c \
    --outFileNamePrefix $OUTPUT \
    --quantMode TranscriptomeSAM
  elif [[ $ZIP == "" ]]; then
    echo "unzip"
    $STAR --runThreadN $THREADS \
    --genomeDir "${INDICESDIR}" \
    --readFilesIn $READ1 $READ2 \
    --outFileNamePrefix $OUTPUT \
    --quantMode TranscriptomeSAM
  else
    echo "$READ1 and/or $READ2 has the wrong format"
  fi
done

SE_STAR=()
for FILE in "${SE[@]}"; do
  echo $FILE
  if [[ $FILE =~ ^(.*)Filtered(\.f[a-z]*)(\.gz|\.bz2)?$ ]]; then
    SAMPLENAME=${BASH_REMATCH[1]}
    FORMAT=${BASH_REMATCH[2]}
    ZIP=${BASH_REMATCH[3]}
    OUTPUT=$STARSINGLEDIR/$(basename $SAMPLENAME)_
    SE_STAR+=($OUTPUT)
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
    --outFileNamePrefix $OUTPUT \
    --quantMode TranscriptomeSAM
  elif [[ $ZIP == ".gz" ]]; then
    echo "gzipped"
    $STAR --runThreadN $THREADS \
    --genomeDir "${INDICESDIR}" \
    --readFilesIn $FILE \
    --readFilesCommand gunzip -c \
    --outFileNamePrefix $OUTPUT \
    --quantMode TranscriptomeSAM
  elif [[ $ZIP == "" ]]; then
    echo "unzip"
    $STAR --runThreadN $THREADS \
    --genomeDir "${INDICESDIR}" \
    --readFilesIn $FILE \
    --outFileNamePrefix $OUTPUT \
    --quantMode TranscriptomeSAM
  else
    echo "$FILE has the wrong format"
    exit
  fi
done
else
  echo "Skipped STAR Alignment"
fi
echo -e "Quantification with RSEM\n"


if [[ $RSEMREF == "no" ]]; then

if [ ! -d "${RSEMREFDIR}" ]; then
  mkdir "${RSEMREFDIR}"
fi

echo -e "Prepare Reference with RSEM\n"
echo $ANNOTATION
if [[ $ANNOTATION =~ ^.*/(.*)(\..*)$ ]]; then
  ANNONAME=${BASH_REMATCH[1]}
  FORMAT=${BASH_REMATCH[2]}
  PREP_REF=$RSEMREF/$ANNONAME

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
  echo "Prepare Reference with RSEM: DONE!\n"
  echo "Reference Files are stored in ${PREP_REF}"
fi
else
  echo "Skipped preparing reference with RSEM"
  PREP_REF=${RSEMREFDIR}
fi

RSEMOUTPUTDIR=${RSEMDIR/Output}
if [[ ! -d "${RSEMOUTPUTDIR}" ]]; then
  mkdir "${RSEMOUTPUTDIR}"
fi


echo "PAIRED END FILES: ${PE_STAR[@]}"
echo "SINGLE END FILES: ${SE_STAR[@]}"

RSEMOUTPUTFILES=()
for FILE in "${PE_STAR[@]}"; do
  echo $FILE
  echo ""

  RSEMINPUT=${FILE}Aligned.toTranscriptome.out.bam
  if [[ -f ${RSEMINPUT} ]]; then
    RSEMOUTPUTFILES+=(${RSEMOUTPUTDIR}/$(basename $FILE))
    echo "Quantification with RSEM: ${RSEMINPUT}"
    ${RSEM}rsem-calculate-expression \
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
    ${RSEM}rsem-calculate-expression \
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
echo "The next step is quality control and visualization with scater"
#if [[ -f "./cell.counts.txt" ]]; then
#  i=1
#else
#  i=0
#fi
i=0
for FILE in "${RSEMOUTPUTFILES[@]}"; do
  CELL=$(basename $FILE)
  RSEMOUTPUT="${FILE}.${RSEMRESULT}.results"
  if [[ $i=0 ]]; then
    awk 'NR==1 {print $RSEMRESULT,$FILE} NR>1 {print $1,$5}' $RSEMOUTPUT > cell.counts.txt
    awk '{print $1,$3}' $RSEMOUTPUT > gene.information.txt

    echo -e "cell\tanimal\tsex\tcondition\ttreatment" > ${METADATA}
    echo -e "$(basename $FILE)\t${ANIMAL}\t${SEX}\t${CONDITION}\t${TREATMENT}" >> ${METADATA}
    i=1
  else
    echo -e "$(basename $FILE)\t${ANIMAL}\t${SEX}\t${CONDITION}\t${TREATMENT}" >> ${METADATA}
    paste cell.counts.txt <(awk 'NR==1 {print $FILE} NR>1 {print $5}' $RSEMOUTPUT) > cell.counts.txt
  fi
done
