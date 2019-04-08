#!/bin/bash
./Pipeline.sh \
--genome test_data/genome/Mus_musculus.GRCm38.dna.primary_assembly.fa \
--annotation test_data/genome/Mus_musculus.GRCm38.95.gtf \
--indicesDir test_data/indices \
--index "yes" \
--data test_data/data/ \
--starDir test_data/STAR_output \
--trimDir test_data/trimmomatic_output/ \
--rsemDir test_data/RSEM \
--trimmomatic ../tools/Trimmomatic-0.36/trimmomatic-0.36.jar \
--FastQC ../tools/FastQC/fastqc \
--RSEM ../tools/RSEM-1.3.1/ \
--STAR ../tools/STAR-2.7.0e/bin/Linux_x86_64/STAR \
--threads 8 \
--rsemResult "genes" \
--rsemRef "yes" \
--rsemRefDir "test_data/RSEM/ref/Mus_musculus.GRCm38.95" \
--bamFiles "no" \
--bamPaired test_data/STAR_output/PairedEnd \
--bamSingle test_data/STAR_output/SingleEnd \
> pipeline.log.out
