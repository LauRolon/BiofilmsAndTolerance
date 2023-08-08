#!/bin/bash
# Usage: sh average_coverage.sh <path to input files>

cd $1
for f in *_contigs.fasta
do

echo "Indexing $f..."
/storage/home/mlr355/work/Programs/bwa/bwa index $f

echo "Mapping reads to $f..."
#Make sure to change the sequence end tags as they sometimes change depending on the sequencing company. Either _1.trimmedP.fastq.gz or _R1_001.trimmedP.fastq.gz
/storage/home/mlr355/work/Programs/bwa/bwa mem -t 8 $f ${f%_contigs.fasta}_R1_001.trimmedP.fastq.gz ${f%_contigs.fasta}_R2_001.trimmedP.fastq.gz > ${f%_contigs.fasta}.sam
echo "SAM file created"

echo "Converting SAM to BAM with samtools..."
samtools view -Sb ${f%_contigs.fasta}.sam -o ${f%_contigs.fasta}.bam
echo "BAM file created."

echo "Removing sam file..."
rm *.sam

echo "Sorting BAM file with samtools..."
#try if this works: samtools sort ${f%_contigs.fasta}.bam ${f%_contigs.fasta}_sorted
samtools sort ${f%_contigs.fasta}.bam -o ${f%_contigs.fasta}_sorted.bam
echo "Finished sorting."

echo "Indexing sorted BAM file..."
samtools index ${f%_contigs.fasta}_sorted.bam
echo "Index complete."

echo "Removing unnecessary files"
rm *.sa
rm *.bai
rm *.bwt
rm *.pac
rm *.ann
rm *.amb

echo "Using samtools depth to obtain average genome coverage..."
X=$(samtools depth ${f%_contigs.fasta}_sorted.bam | awk '{sum+=$3} END { print sum/NR}');
echo "${f%_contigs.fasta}_sorted.bam";
echo "$X";
echo "${f%_contigs.fasta}_sorted.bam $X" > ${f%_contigs.fasta}_coverage.txt;
done

cat *_coverage.txt > Coverage_strains.txt
rm *.bam

