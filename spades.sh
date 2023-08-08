#!/bin/bash
#Usage: sh spades.sh <path to input files>

cd $1

for f in *_R1_001.trimmedP.fastq.gz
do
if [ -d "${f%_R1_001.trimmedP.fastq.gz}" ]
then
echo 'skip'${f}
continue
fi
echo 'assemble' ${f%_R1_001.trimmedP.fastq.gz}
python3 /storage/home/mlr355/work/Programs/SPAdes-3.15.3-Linux/bin/spades.py -k 99,127 --isolate -1 $f -2 ${f%_R1_001.trimmedP.fastq.gz}_R2_001.trimmedP.fastq.gz -o ${f%_R1_001.trimmedP.fastq.gz} -t 4;
done

mkdir contigs
for f in *_R1_001.trimmedP.fastq.gz
do
        cd ${f%_R1_001.trimmedP.fastq.gz}
        cat contigs.fasta > ${f%_R1_001.trimmedP.fastq.gz}_contigs.fasta
        cp ${f%_R1_001.trimmedP.fastq.gz}_contigs.fasta ../contigs
        cd ..;
        done

mkdir scaffolds
for f in *_R1_001.trimmedP.fastq.gz
        do
        cd ${f%_R1_001.trimmedP.fastq.gz}
        cat scaffolds.fasta > ${f%_R1_001.trimmedP.fastq.gz}_scaffolds.fasta
        cp ${f%_R1_001.trimmedP.fastq.gz}_scaffolds.fasta ../scaffolds
        cd ..;
done

