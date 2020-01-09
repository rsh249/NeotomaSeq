#!/usr/bin/bash
##################


##This script analyzes ancient DNA damage patterns.

##Authors:
#Robert Harbert: rharbert@stonehill.edu
#Grace Moore



###################

taxsearch=Artemisia
refdir=/usr/share/data/ncbi/genomes/plastids

reffa=$refdir/plastid_all.fna

#python3 ../../bin/parse_reads.py Taxa.ResultsAfter *.centclass *.rmdup.fq $reffa $taxsearch

awk -v pattern="$taxsearch " '/^>/ { p = ($0 ~ pattern)} p' $reffa > $taxsearch.fa

 
bwa index $reffa
bwa mem -I 0 -L 1024 -E 7 -t 24 $reffa *.rmdup.fq > $search.sam

bwa index $taxsearch.fa
bwa mem -I 0 -L 1024 -E 7 -t 24 $taxsearch.fa *.rmdup.fq > $search.sam
#bwa mem -I 0 -N 0.02 -L 1024 -E 7 -t 24 $reffa *.rmdup.fq > $search.sam
#bwa aln -t 8 -i 0 -o 2 -n 0.02 -l 1024 -e 7 $search.fa *.filtered > $search.sai
#bwa samse $search.fa $search.sai *.filtered

samtools view -Shb $search.sam > $search.bam

samtools view -b -F4 $search.bam > $search.final.bam
samtools sort $search.final.bam > $search.sort.bam
samtools index $search.sort.bam

mapDamage -i $search.sort.bam -r $reffa


