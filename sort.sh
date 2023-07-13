#!/bin/bash

for i in *.sam

do
samtools view -bS ${i}| samtools sort -@8 -o ${i/.sam/.bam}
samtools index ${i/.sam/.bam}
done
