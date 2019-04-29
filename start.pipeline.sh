#!/bin/bash
./Pipeline.sh \
--projectName "Bassoon_scRNA-seq" \
--output /media/data/Daniel/test_data/Output \
--genome /media/data/Daniel/test_data/genome/Mus_musculus.GRCm38.dna.primary_assembly.fa \
--annotation /media/data/Daniel/test_data/genome/Mus_musculus.GRCm38.95.gtf \
--indicesDir /media/data/Daniel/test_data/indices \
--index "no" \
--data /media/data/Daniel/test_data/data/ \
--trimmomatic /media/data/tools/Trimmomatic-0.36/trimmomatic-0.36.jar \
--fastqc /media/data/tools/FastQC/fastqc \
--RSEM /media/data/tools/RSEM-1.3.1/ \
--STAR /media/data/tools/STAR-2.7.0e/bin/Linux_x86_64/STAR \
--threads 8 \
--rsemResult "genes" \
--rsemRef "no" \
--rsemRefDir /media/data/Daniel/test_data/RSEM/ref/ \
--bamFiles "no" \
--impute "no" \
--metaData /media/data/Daniel/cells.metadata.txt \
--animal "mouse" \
--sex "F" \
--condition "Bassoon Knockout" \
--treatment "treatment" \
> /media/data/Daniel/Scripts/pipeline.log.out

# For every single cell sample create one directory with the sequencing files
# Do not change project name, while running pipeline multiple times
# reference genome/transcriptome needs to be in .fa format
# annotation file needs .gtf format

# Software used:
# Trimmomatic: -preinstalled --- Version 0.36
# FastQC: -preinstalled --- FastQC v0.11.3
# STAR: https://github.com/alexdobin/STAR --- Version 2.7.0e
# RSEM: https://deweylab.github.io/RSEM/ --- Version 1.3.1
# scater: R-Package - Bioconductor --- Version 1.10.1
# scran: R-Package - Bioconductor --- Version 1.10.2
# Seurat: R-Package --- Version 3.0.0
# dplyr: R-Package --- Version 0.8.0.1
# ggplot2: R-Package --- Version 3.1.1
