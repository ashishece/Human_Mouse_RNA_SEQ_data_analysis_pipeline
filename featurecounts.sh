#!/bin/bash/
for i in *.bam

do

featureCounts -p -O -T 40 -a /scratch/ashi/Jan22a-fastq_10/human_analysis_scripts/Mus_musculus.GRCm39.108.gtf -o ${i}.counts.txt  ${i}

done
