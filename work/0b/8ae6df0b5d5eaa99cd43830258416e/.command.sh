#!/bin/bash -ue
mkdir genomeDir
cp ebola_ref.fasta genomeDir/genome.fa
samtools faidx genomeDir/genome.fa
gatk CreateSequenceDictionary --REFERENCE genomeDir/genome.fa --OUTPUT=genomeDir/genome.dict
STAR --runMode genomeGenerate --genomeDir genomeDir --genomeFastaFiles genomeDir/genome.fa --genomeSAindexNbases 6 --runThreadN 4
