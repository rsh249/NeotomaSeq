#!/usr/bin/bash
##################


##This script analyzes ancient DNA damage patterns.

##Authors:
#Robert Harbert: rharbert@stonehill.edu
#Grace Moore



################### Filter plastid genomes for expected plants
sample=$1
echo $sample

taxsearch=$2

reffa=$3



mkdir $sample
cd $sample


### Generate read length file
cat ../../results/$sample*/$sample*.dust.rmdup.fq | sed 's/N//g' | awk '{if(NR%4==2) print length($1)}' > read_length.txt


#Set up reference index
cp $reffa $taxsearch.fa
bwa index $taxsearch.fa


#### Generate bam file
bwa mem -Y -I 0 -L 1024 -E 7 -t 5 $taxsearch.fa ../../results/$sample*/$sample*.dust.rmdup.fq > $taxsearch.sam

samtools view -S -b $taxsearch.sam > $taxsearch.bam
samtools view -b -F4 $taxsearch.bam > $taxsearch.final.bam
samtools sort $taxsearch.final.bam > $taxsearch.sort.bam
samtools index $taxsearch.sort.bam

rm $taxsearch.bam
rm $taxsearch.sam
rm $taxsearch.fa

### Filter with PMDtools 

#samtools view -h $taxsearch.sort.bam | python2 ~/GitHub/NeotomaSeq/bin/PMDtools/pmdtools.0.60.py --threshold 0 --header | samtools view -Sb - > $taxsearch.pmds0filter.bam
#samtools view -h $taxsearch.sort.bam | python2 ~/GitHub/NeotomaSeq/bin/PMDtools/pmdtools.0.60.py --threshold 3 --header | samtools view -Sb - > $taxsearch.pmds3filter.bam

#PMD vis C->T 5' deamination

samtools view -h $taxsearch.sort.bam | python2 ~/GitHub/NeotomaSeq/bin/PMDtools/pmdtools.0.60.py --deamination --range 125 --CpG > $taxsearch.CpG
samtools view -h $taxsearch.sort.bam | python2 ~/GitHub/NeotomaSeq/bin/PMDtools/pmdtools.0.60.py --deamination --range 125 > $taxsearch.deam

#samtools view -h $taxsearch.pmds0filter.bam | python2 ~/GitHub/NeotomaSeq/bin/PMDtools/pmdtools.0.60.py --deamination --range 125 --CpG > $taxsearch.pmds0filter.CpG
#samtools view -h $taxsearch.pmds0filter.bam | python2 ~/GitHub/NeotomaSeq/bin/PMDtools/pmdtools.0.60.py --deamination --range 125 > $taxsearch.pmds0filter.deam


#samtools view -h $taxsearch.pmds3filter.bam | python2 ~/GitHub/NeotomaSeq/bin/PMDtools/pmdtools.0.60.py --deamination --range 125 --CpG > $taxsearch.pmds3filter.CpG
#samtools view -h $taxsearch.pmds3filter.bam | python2 ~/GitHub/NeotomaSeq/bin/PMDtools/pmdtools.0.60.py --deamination --range 125 > $taxsearch.pmds3filter.deam

