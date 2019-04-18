#!/bin/bash
./Pipeline.sh \
--projectName "Bassoon_scRNA-seq" \
--output /media/data/Daniel/test_data/Output \
--genome /media/data/Daniel/test_data/genome/Mus_musculus.GRCm38.dna.primary_assembly.fa \
--annotation /media/data/Daniel/test_data/genome/Mus_musculus.GRCm38.95.gtf \
--indicesDir /media/data/Daniel/test_data/indices \
--index "no" \
--data /media/data/Daniel/test_data/data/ \
--starDir /media/data/Daniel/test_data/STAR_output \
--trimDir /media/data/Daniel/test_data/trimmomatic_output/ \
--rsemDir /media/data/Daniel/test_data/RSEM \
--trimmomatic /media/data/tools/Trimmomatic-0.36/trimmomatic-0.36.jar \
--FastQC /media/data/tools/FastQC/fastqc \
--RSEM /media/data/tools/RSEM-1.3.1/ \
--STAR /media/data/tools/STAR-2.7.0e/bin/Linux_x86_64/STAR \
--threads 8 \
--rsemResult "genes" \
--rsemRef "no" \
--rsemRefDir /media/data/Daniel/test_data/RSEM/ref/ \
--bamFiles "no" \
--bamPaired /media/data/Daniel/test_data/STAR_output/PairedEnd \
--bamSingle /media/data/Daniel/test_data/STAR_output/SingleEnd \
--metaData /media/data/Daniel/cells.metadata.txt \
--animal "mouse" \
--sex "F" \
--condition "Bassoon Knockout" \
--treatment "no treatment" \
> /media/data/Daniel/Output/pipeline.log.out
