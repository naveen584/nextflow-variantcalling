#!/bin/bash -ue
mkdir star_rundir
cd star_rundir
STAR --genomeDir ../genomeDir/ --readFilesIn ../ebola_mut_reads1.fq ../ebola_mut_reads2.fq --runThreadN 4
