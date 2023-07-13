#!/bin/bash
for i in *_R1_001.fastq.gz
do
j=${i/_R1_001.fastq.gz/_R2_001.fastq.gz}
fastp -i ${i} -I ${j}  -o ${i/_R1_001.fastq.gz/_R1_trim.fastq.gz} -O ${j/_R2_001.fastq.gz/_R2_trim.fastq.gz} -l 36 \
     ${i} ${j}

done

