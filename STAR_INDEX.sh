#!/bin/bash/

STAR --runMode genomeGenerate --genomeDir /scratch/ashi/Jan22a-fastq_10/new_data/hisat  --genomeFastaFiles /scratch/ashi/Jan22a-fastq_10/new_data/hisat/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile /scratch/ashi/Jan22a-fastq_10/new_data/hisat/Homo_sapiens.GRCh38.108.gtf --sjdbOverhang 100 --runThreadN 20  
