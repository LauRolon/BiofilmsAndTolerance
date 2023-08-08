#!/bin/bash
#Usage: sh fastqc.sh <path to input files>
#Note: provide a full path to the program

cd $1

for f in *.fastq.gz
do
/storage/home/mlr355/work/Programs/FastQC/fastqc $f -o /storage/home/mlr355/work/USDA/0_reads/New/1_fastqc
done
