##Analyses of WGS of environmental isolates from tree fruit packing houses
#Last updated: 05/28/22 MLR

cd /storage/home/mlr355/work/USDA

#Raw reads are stored in 0_reads folder within the directory

#1. Assess read quality with FasQC v0.11.9
mkdir 1_fastqc
cd ..
sh fastqc.sh 0_reads

#2 Trim the adapters and poor-quality bases with Trimmomatic v 0.39
mkdir 2_trimmed
sh trimmomatic.sh 0_reads

cd 0_reads
cp *.trimmedP.fastq.gz ../2_trimmed
cp *.trimmedS.fastq.gz ../2_trimmed
rm *.trimmedP.fastq.gz
rm *.trimmedS.fastq.gz

cd ..

#3 de novo assembly using SPAdes v 3.15.3 - long processing time 
mkdir 3_contigs
sh spades.sh 2_trimmed

#4 assembly quality control using Quast 5.0.2
cd contigs
cp *.fasta ../../3_contigs
cd ../..
sh quast.sh 3_contigs

#Combine all tsv reports in one file
cd quast_reports
cat *_reports.tsv > quast_summary.tsv

#The report will have all the headers, manually remove all odd rows except #1

#Calculate average coverage for the genomes
#Navigate to USDA folder
sh average_coverage.sh <input dir>


####SNP ANALYSIS
#Will do SNP analysis within each bacterial family
#Make a directory for every bacterial family and sort contig files to each expected taxonomy

cd ../..
mkdir 4_genomes
cd 4_genomes
mkdir Pseudo
mkdir Micro
mkdir Xantho
mkdir Flavo

#Sort files manually

module use /gpfs/group/RISE/sw7/modules

#We will then load a module
module load cfsan-snp

#Make infile for ksnp3
MakeKSNP3infile Micro Micro_infile A
MakeKSNP3infile Pseudo Pseudo_infile A
MakeKSNP3infile Xantho Xantho_infile A
MakeKSNP3infile Flavo Flavo_infile A

#Run kSNP3 for each bacterial family
kSNP3 -in Micro_infile -outdir ksnp_micro -k 21 -core -CPU 4 -ML
kSNP3 -in Pseudo_infile -outdir ksnp_pseudo -k 21 -core -CPU 4 -ML
kSNP3 -in Xantho_infile -outdir ksnp_xantho -k 21 -core -CPU 4 -ML
kSNP3 -in Flavo_infile -outdir ksnp_flavo -k 21 -core -CPU 4 -ML

#Make TREE based on kSNP results
#NAvigate to ksnp output directory and make directory for raxml output
#copy SNP file named core_SNPs_matrix.fasta
#Run tree with 500 bootstraps for publication quality

tmux
conda activate raxml
raxmlHPC -f a -x 165 -m ASC_GTRGAMMA --asc-corr=lewis -p 596 -N 500 -s core_SNPs_matrix.fasta -n ASC_core_SNPs.tre

#Open .tre file in FigTree to view results

## Use BARNAP to extract 16s sequence from WGS contigs
conda install -c bioconda -c conda-forge barrnap
sh barrnap.sh <path to contig files>  


#Move to directory containing barrnap results
rm *fai
#Add filename to headers in barrnap output
for f in *_16s.fasta; do sed -i "s/^>/>${f}_/" "$f"; done

#Extract 16s seq from barrnap output into one fasta file
cat *_16s.fasta | grep -A 1 '16S_rRNA' > Xantho_selected.fasta 

#Move 16ss files to indepedent folder
mkdir 16s
cp *_16.fasta 16s
rm *_16s.fasta




