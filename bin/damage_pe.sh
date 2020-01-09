#!/usr/bin/bash
##################


##This script analyzes ancient DNA damage patterns.

##Authors:
#Robert Harbert: rharbert@stonehill.edu
#Grace Moore



################### Filter plastid genomes for expected plants
sample=$1
taxsearch=Artemisia
refdir=/usr/share/data/ncbi/genomes/plastids

reffa=$refdir/plastid_all.fna

#python3 ../../bin/parse_reads.py Taxa.ResultsAfter *.centclass *.rmdup.fq $reffa $taxsearch


awk -v pattern="$taxsearch " '/^>/ { p = ($0 ~ pattern)} p' $reffa > $taxsearch.fa
bwa index $taxsearch.fa

#### Generate bam file

bwa mem -Y -I 0 -L 1024 -E 7 -t 24 $taxsearch.fa ../../fastq/TS211A* > $taxsearch.sam
samtools view -S -b $taxsearch.sam > $taxsearch.bam
samtools view -b -F4 $taxsearch.bam > $taxsearch.final.bam
samtools sort $taxsearch.final.bam > $taxsearch.sort.bam
samtools index $taxsearch.sort.bam

### Filter with PMDtools 

samtools view -h $taxsearch.bam | python2 ~/GitHub/NeotomaSeq/bin/PMDtools/pmdtools.0.60.py --threshold 0 --header | samtools view -Sb - > pmdbam.pmds3filter.bam


### mapDamage
mapDamage -i pmdbam.pmds3filter.bam --seq-length=100 -r $taxsearch.fa --no-stats

### vis with PMDtools

samtools view $taxsearch.bam | python2 ~/GitHub/NeotomaSeq/bin/PMDtools/pmdtools.0.60.py -- platypus > PMD_temp.txt

R CMD BATCH ~/GitHub/NeotomaSeq/bin/plotPMD.v2.R
