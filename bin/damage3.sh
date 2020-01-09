#!/usr/bin/bash
##################


##This script analyzes ancient DNA damage patterns.

##Authors:
#Robert Harbert: rharbert@stonehill.edu
#Grace Moore



###################
module load samtools-1.5 
module load python-3.6.1 
module load bwa-0.7.15
#Use conda R 3.5 'conda install -c r r-base'
#jobdir=$1
taxsearch=$1
refdir=~/nas5/NeotomaSeq/refs

reffa=$refdir/genomes/plastids.fna
#reffa=$refdir/nt.fa

#python3 ../../bin/parse_reads.py Taxa.ResultsAfter *.centclass *.rmdup.fq $reffa $taxsearch

awk -v pattern="$taxsearch " '/^>/ { p = ($0 ~ pattern)} p' $reffa > $taxsearch.fa

 
bwa index $taxsearch.fa

bwa mem -I 0 -N 0.02 -L 1024 -E 7 $taxsearch.fa *.rmdup.fq > $taxsearch.sam
bwa aln -t 8 -i 0 -o 2 -n 0.02 -l 1024 -e 7 $taxsearch.fa *.rmdup.fq > $taxsearch.sai
bwa samse $taxsearch.fa $taxsearch.sai *.rmdup.fq

samtools view -Shb $taxsearch.sam > $taxsearch.bam

samtools view -b -F4 $taxsearch.bam > $taxsearch.final.bam
samtools sort $taxsearch.final.bam > $taxsearch.sort.bam
samtools index $taxsearch.sort.bam

mapDamage -i $taxsearch.sort.bam -r $taxsearch.fa



#use R3.5 and aRchaic for alternative
#bam to mff
#if needed run: wget https://raw.githubusercontent.com/kkdey/aRchaic/master/bin/generate_summary_bams.py ##Probably can't distribute this script in the GitHub repo


python ../../bin/generate_summary_bams.py -b $taxsearch.sort.bam -f $taxsearch.fa -o $taxsearch.mff.csv




