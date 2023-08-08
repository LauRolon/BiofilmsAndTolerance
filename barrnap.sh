#!/bin/bash
#Usage: sh barrnap.sh <path to input files>

cd $1

for f in *.fasta
do
barrnap -q -k bac $f --outseq ${f}_16s.fasta;
done

