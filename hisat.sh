#!/bin/bash

for read1 in *R1_trim.fastq.gz

do
read2=${read1/R1_trim.fastq.gz/R2_trim.fastq.gz}
echo here is R1 $read1
echo here is R2 $read2
hisat2 -p 30 -q --very-sensitive -x /scratch/ashi/Jan22a-fastq_10/human_analysis_scripts/genome -1 ${read1} -2 ${read2} -S ${read1/R1_trim.fastq.gz/.sam} --summary-file ${read1/R1_trim.fastq.gz/}.txt
done
