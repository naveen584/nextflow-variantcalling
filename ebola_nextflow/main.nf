#!/usr/bin/env nextflow

// Validate inputs
// Genome
if ( params.genome ) {
  g = Channel
    .fromPath(params.genome, checkIfExists: true)
    .ifEmpty { exit 1, "Genome FASTA file not found: ${params.genome}" }
} else {
  exit 1, "No genome FASTA file specified!"
}

// Reads 1
if ( params.reads1 ) {
  r1 = Channel
    .fromPath(params.reads1, checkIfExists: true)
    .ifEmpty { exit 1, "Reads 1 FASTQ file not found: ${params.reads1}" }
} else {
  exit 1, "No reads 1 FASTQ file specified!"
}

//Reads 2
if ( params.reads2 ) {
  r2 = Channel
    .fromPath(params.reads2, checkIfExists: true)
    .ifEmpty { exit 1, "Reads 2 FASTQ file not found: ${params.reads2}" }
} else {
  exit 1, "No reads 2 FASTQ file specified!"
}

log.info """========================================================================
Ebola - Nextflow pipeline v${params.version}
========================================================================"""

def summary = [:]
summary['Pipeline Name']  = 'Ebola - Nextflow pipeline'
summary['Pipeline Version'] = params.version
summary['Run Name'] = workflow.runName
summary['User'] = System.getProperty("user.name")
summary['Config Profile'] = workflow.profile
summary['Working Dir'] = workflow.workDir
summary['Genome'] = params.genome
summary['Reads 1'] = params.reads1
summary['Reads 2'] = params.reads2
summary['Output Dir'] = params.outdir
summary['Threads'] = params.threads
summary['SA Index Nbases'] = params.sa
summary['RGID'] = params.rgid
summary['RGLB'] = params.rglb
summary['RGPL'] = params.rgpl
summary['RGSM'] = params.rgsm
summary['RGPU'] = params.rgpu
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "========================================================================"

// Check Nextflow version
try {
  if( ! nextflow.version.matches(">= $params.nf_required_version") ){
    throw GroovyException('Nextflow version too old')
  }
} catch (all) {
  log.error "====================================================\n" +
            "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
            "  Pipeline execution will continue, but things may break.\n" +
            "  Please run `nextflow self-update` to update Nextflow.\n" +
            "============================================================"
}

/*
  STEP 1: Genome Preparation
*/
process genomePreparation {
  input:
  file genome from g

  output:
  file "genomeDir" into genome_dir

  script:
  genomeSAindexNbases = "--genomeSAindexNbases ${params.sa}"
  runThreadN = "--runThreadN ${params.threads}"
  """
  mkdir genomeDir
  cp $genome genomeDir/genome.fa
  samtools faidx genomeDir/genome.fa
  gatk CreateSequenceDictionary --REFERENCE genomeDir/genome.fa --OUTPUT=genomeDir/genome.dict
  STAR --runMode genomeGenerate --genomeDir genomeDir --genomeFastaFiles genomeDir/genome.fa $genomeSAindexNbases $runThreadN
  """
}

genome_dir.into { genome_dir_1; genome_dir_2 }

/*
  STEP 2: Run STAR
*/
process runStar {
  input:
  file genome from genome_dir_1.collect()
  file reads1 from r1
  file reads2 from r2

  output:
  file "star_rundir" into star_rundir

  script:
  runThreadN = "--runThreadN ${params.threads}"
  """
  mkdir star_rundir
  cd star_rundir
  STAR --genomeDir ../$genome/ --readFilesIn ../$reads1 ../$reads2 $runThreadN
  """
}

/*
  STEP 3: Run GATK4
*/
process runGatk {
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file genome from genome_dir_2.collect()
  file star_rundir from star_rundir.collect()

  output:
  file "Aligned_sorted_RG_index.marked_split.bam" into bam
  file "Aligned_sorted_RG_dup_metrics" into metrics
  file "variants.vcf" into variants
  file "indels.vcf" into indels

  script:
  rgid = "--RGID=${params.rgid}"
  rglb = "--RGLB=${params.rglb}"
  rgpl = "--RGPL=${params.rgpl}"
  rgsm = "--RGSM=${params.rgsm}"
  rgpu = "--RGPU=${params.rgpu}"
  """
  mkdir gatk_rundir
  cd gatk_rundir
  gatk SamFormatConverter -I ../$star_rundir/Aligned.out.sam -O Aligned.out.bam
  gatk SortSam -I Aligned.out.bam -O Aligned_sorted.bam -SO coordinate
  gatk AddOrReplaceReadGroups -I Aligned_sorted.bam -O Aligned_sorted_RG.bam --SORT_ORDER=coordinate $rgid $rglb $rgpl $rgsm $rgpu
  gatk BuildBamIndex -I Aligned_sorted_RG.bam
  gatk MarkDuplicates -I Aligned_sorted_RG.bam -O Aligned_sorted_RG_index.marked.bam --METRICS_FILE=Aligned_sorted_RG_dup_metrics --VALIDATION_STRINGENCY=LENIENT --CREATE_INDEX=true --REMOVE_DUPLICATES=true
  gatk SplitNCigarReads -R ../$genome/genome.fa -I Aligned_sorted_RG_index.marked.bam -O Aligned_sorted_RG_index.marked_split.bam
  gatk HaplotypeCaller -R ../$genome/genome.fa -I Aligned_sorted_RG_index.marked_split.bam -O variants.vcf
  gatk SelectVariants -R ../$genome/genome.fa -V variants.vcf -O indels.vcf --select-type-to-include INDEL
  mv Aligned_sorted_RG_index.marked_split.bam Aligned_sorted_RG_dup_metrics variants.vcf indels.vcf ../
  """
}
