#!/usr/bin/bash



home=/home/user/path/to/repository
cd $home #Just in case

#put your sequence files in a folder called fastq in the repository:
targs=$(ls fastq | awk -F "_" '{print $1}' | uniq)

#put a Centrifuge database at refs/centr called nucind2
centdb=$home/refs/centr/nucind2
#put a fasta copy of the NCBI 'nt' database in refs as well.
#'wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz' ./refs/
reffa=$home/refs/nt.fa

IFS=' ' read -r -a array <<< $targs

#Run for 0 to n fastq read pair files:
./bin/NeotomaSeq.sh $centdb $reffa ${array[0]}
./bin/NeotomaSeq.sh $centdb $reffa ${array[1]}
./bin/NeotomaSeq.sh $centdb $reffa ${array[2]}
./bin/NeotomaSeq.sh $centdb $reffa ${array[3]}


