#!/bin/bash

input="/media/data2/Daniel/test_pipeline_data//"
output="/media/data2/Daniel/Output/test/"

for dir in $(ls $input); do
  echo "Start pipeline with Sample $dir!"
  ./Pipeline.10x.sh \
    --output ${output}/$dir \
    --genome /media/data/Daniel/data/genome/Mus_musculus.GRCm38.dna.primary_assembly.fa \
    --annotation /media/data/Daniel/data/genome/Mus_musculus.GRCm38.96.gtf \
    --indicesDir /media/data/Daniel/data/indices \
    --qualityControl "no" \
    --trimming "no" \
    --data ${input}/${dir} \
    --read "R2" \
    --barcode "R1" \
    --useCellranger "yes" \
    --CRoptions "" \
    --useSTARsolo "no" \
    --useUMItools "no" \
    --STARoptions "" \
    --cellrangerTranscriptome /media/data2/Daniel/cellranger_test/refdata-cellranger-mm10-3.0.0 \
    --umi-tools /media/data/Daniel/src/anaconda3/bin/umi_tools \
    --samtools /media/data/Daniel/src/anaconda3/bin/samtools \
    --trimmomatic /media/data/tools/Trimmomatic-0.36/trimmomatic-0.36.jar \
    --fastqc /media/data/tools/FastQC/fastqc \
    --featureCounts /media/data/tools/subread-1.6.4-Linux-x86_64/bin/featureCounts \
    --star /media/data/tools/STAR-2.7.0e/bin/Linux_x86_64/STAR \
    --index "y" \
    --threads 4 \
    --STARwhitelist "/media/data2/Daniel/cellranger_test/cellranger-3.0.2/cellranger-cs/3.0.2/lib/python/cellranger/barcodes/737K-august-2016.txt" \
    --UMITOOLSwhitelist "" \
    --genWhitelist "no" \
    --trimOptions "TRAILING:20 HEADCROP:20 MINLEN:75" \
    --useLanes "all" \
    --nGenes "100" \
    --nUMIs "" \
    --MAD "" \
    --abundantMT 0.5 \
    --filterGenes 0.001 \
    --normalize "yes" \
    --cellranger "/media/data/Daniel/cellranger-3.0.2/cellranger" \
    | tee ${output}/${dir}.Pipeline.log.out
done

# For every single cell sample create one directory with the sequencing files
# Do not change project name, while running pipeline multiple times
# reference genome/transcriptome needs to be in .fa format
# annotation file needs .gtf format

# Software used:
# Trimmomatic: -preinstalled --- Version 0.36
# FastQC: -preinstalled --- FastQC v0.11.3
# STAR: https://github.com/alexdobin/STAR --- Version 2.7.0e
# scater: R-Package - Bioconductor --- Version 1.10.1
# scran: R-Package - Bioconductor --- Version 1.10.2
# Seurat: R-Package --- Version 3.0.0
# dplyr: R-Package --- Version 0.8.0.1
# ggplot2: R-Package --- Version 3.1.1

# reference genome file downloaded from ensembl - 30.04.19 - GRCm38.dna.primary_assembly.fa.gz
# annotation file downloaded from enselmbl - 30.04.19 - GRCm38.96.gtf.gz

# -- transcript to gene file 
# http://www.ensembl.org/biomart/martview/5327b67ede7aff19118cafed11db2417
# Mouse genes (GRCm38.p6) Download tsv.gz table 04.07.19
