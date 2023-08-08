#!/bin/bash
#Usage: sh trimmomatic.sh <path to input files>

cd $1
echo | pwd
for f in *.fastq.gz
do
if [ -f "${f%_R1_001.fastq.gz}_R1_001.trimmedP.fastq.gz}" ]
then
echo 'skip'${f}
continue
fi
echo 'trim' ${f}
java -jar /storage/home/mlr355/work/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 10 -phred33 -trimlog log $f ${f%_R1_001.fastq.gz}_R2_001.fastq.gz ${f%_R1_001.fastq.gz}_R1_001.trimmedP.fastq.gz ${f%_R1_001.fastq.gz}_R1_001.trimmedS.fastq.gz ${f%_R1_001.fastq.gz}_R2_001.trimmedP.fastq.gz ${f%_R1_001.fastq.gz}_R2_001.trimmedS.fastq.gz ILLUMINACLIP:/storage/home/mlr355/work/Programs/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36;
done
