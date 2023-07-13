#!/bin/bash/
for i in *.aligned.sortedByCoord.bam
 do

htseq-count -s no -r pos —t exon -i gene_name -f bam ${i}.sortedByCoord.out.bam /scratch/ashi/Jan22a-fastq_10/new_data/hisat/Homo_sapiens.GRCh38.108.gtf > ${i}.counts
