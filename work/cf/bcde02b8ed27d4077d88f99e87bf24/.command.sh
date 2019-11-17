#!/bin/bash -ue
mkdir gatk_rundir
cd gatk_rundir
gatk SamFormatConverter -I ../star_rundir/Aligned.out.sam -O Aligned.out.bam
gatk SortSam -I Aligned.out.bam -O Aligned_sorted.bam -SO coordinate
gatk AddOrReplaceReadGroups -I Aligned_sorted.bam -O Aligned_sorted_RG.bam --SORT_ORDER=coordinate --RGID=445_RG --RGLB=445_LB_PAIRED --RGPL=ILLUMINA --RGSM=ebola --RGPU=NONE
gatk BuildBamIndex -I Aligned_sorted_RG.bam
gatk MarkDuplicates -I Aligned_sorted_RG.bam -O Aligned_sorted_RG_index.marked.bam --METRICS_FILE=Aligned_sorted_RG_dup_metrics --VALIDATION_STRINGENCY=LENIENT --CREATE_INDEX=true --REMOVE_DUPLICATES=true
gatk SplitNCigarReads -R ../genomeDir/genome.fa -I Aligned_sorted_RG_index.marked.bam -O Aligned_sorted_RG_index.marked_split.bam
gatk HaplotypeCaller -R ../genomeDir/genome.fa -I Aligned_sorted_RG_index.marked_split.bam -O variants.vcf
gatk SelectVariants -R ../genomeDir/genome.fa -V variants.vcf -O indels.vcf --select-type-to-include INDEL
mv Aligned_sorted_RG_index.marked_split.bam Aligned_sorted_RG_dup_metrics variants.vcf indels.vcf ../
