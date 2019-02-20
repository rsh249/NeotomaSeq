#!/usr/bin/bash
##################


###This is the main pipeline script for the NeotomaSeq ancient DNA metagenomics project.

##Authors:
#Robert Harbert: rharbert@stonehill.edu
#Grace Moore


###################

##NOTE: This script should be run in the main NeotomaSeq directory.

if [ $# -gt 3 ]; then
	echo "Accepted $# arguments"
else
	echo "Please enter arguments for input and output files"
fi

scriptname=$0
centdb=$1
reffa=$2
targ=$3
home=$(pwd)


echo "Setup:
	Main dir at $home
	Centrifuge DB at $centdb	
	FALCON DB at $reffa
	Target file name $targ
	"
cd results
mkdir $targ

cd $targ

#remove adapters and merge reads
AdapterRemoval --collapse  --file1 $home/fastq/$targ*.R1.fastq --file2 $home/fastq/$targ*.R2.fastq --basename $targ --threads 6

########### PREPROCESS ##############
$home/bin/sga_prep $targ

########### CENTRIFUGE ##############
centrifuge -p 32 -q -x $centdb -U $targ.dust.rmdup.fq -S $targ.centclass --report-file $targ.report.tab --out-fmt tab --ignore-quals -t --non-deterministic --min-hitlen 75 -k 500
Rscript --vanilla $home/bin/TableScript.R $targ.report.tab $home/taxdmp

############# FALCON ################
echo $reffa
FALCON -v -n 48 -t 300000 -F -Z -m 6:1:0:0/0 -m 12:20:0:0/0 -m 14:50:1:1/1 -m 20:200:1:5/10 -c 150 $targ.dust.rmdup.fq $reffa
Rscript --vanilla $home/bin/falconvis.R top.csv $home/taxdmp
