#!/bin/bash
#Usage: sh quast.sh <path to input files>

cd $1

mkdir quast_results

for f in *.fasta
do
python /storage/home/mlr355/work/Programs/quast-5.0.2/quast.py -o ./quast_results/quast_${f%_contigs.fasta} --min-contig 200 $f
done

mkdir quast_reports
for f in *.fasta
do
cd quast_results/quast_${f%_contigs.fasta}
cat transposed_report.tsv > ${f%_contigs.fasta}_report.tsv
cp ${f%_contigs.fasta}_report.tsv ../../quast_reports
cd ../../
done
