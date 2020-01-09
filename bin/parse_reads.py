#!/usr/bin/python3

#Import from biopython library (included in anaconda)
from Bio import SeqIO
import sys
#from multiprocessing import Pool
#from itertools import product
import os
    
#usage
#python3 parse_reads Taxa.ResultsAfter *.centclass *.dust.rmdup.fq Streptophyta

###Arguments to be passed
#sys.argv[1] is the path to the Taxa.ResultsAfter file
#sys.argv[2] is the path to the centclass file (usually *.centclass if run from results dir)
#sys.argv[3] is the path to the reads fq file (Usually *.rmdup.fq)
#sys.argv[4] is the genus of interest

def filter(records, readIDs):
    """
    Filter out parsed fastq file for read ID values and compare them to those found in the .centclass file
    @param records Parsed fastq file
    @param readIDS List of IDs from .centclass file
    """
    matches = []
    for rec in records:
        for readID in readIDs:
            #If the ID in the fastq matches an ID in the .centclass, add it to matches
            if rec.id == readID:
                matches.append(rec)
            else:
                continue
    return matches

def main():
    linesTaxID = []
    #Read in Tax.Results file
    with open (str(sys.argv[1]), 'rt') as orderfile:  
        for line in orderfile:
            linesTaxID.append(line)
    #Specific classification of interest
    substr = sys.argv[4]
    rodentID = []
    #Pull out taxonomic IDs of interest based on classification level
    for line in linesTaxID[1: ]:
        order = line.split("\t")[1]
        taxID = int(line.split("\t")[8])
        if order == substr:
            if taxID not in rodentID:
                rodentID.append(taxID)
            else:
                continue
        else:
            continue   
    linesReadID = []
    with open (str(sys.argv[2]), 'rt') as readfile:  
        for line in readfile:
            linesReadID.append(line)
    fileP = str(sys.argv[2]).split("/")
    name = fileP[-1].split(".")
    ccIDFile = name[0] + ".ReadIDs"
    accIDFile = name[0] + ".accIDs"
    ccIDs = open(ccIDFile, "w")
    accIDs = open(accIDFile, "w")
    #Pull out read IDs that align with the taxonomic ID
    for line in linesReadID[1: ]:
        try:
            readID = line.split('\t')[0]
            accID = line.split('\t')[1]
            taxRead = int(line.split("\t")[2])
        except:
            continue
        #Remove the 'M_' that is in front of the read ID in .centclass but not in fastq
        readID = readID[:0] + readID[2:]
        for ID in rodentID:
            if taxRead == ID:
                accIDs.write(accID)
                accIDs.write('\n')
                ccIDs.write(readID)
                ccIDs.write('\n')
            else:
                continue
    ccIDs.close()
    accIDs.close()
    sort = 'sort ' + ccIDFile + " | uniq > " + ccIDFile + "2"
    os.system(sort)
    
    os.system("LC_ALL=C")
    outfile = name[0] + ".filtered"
    fqfile = sys.argv[3]
    grep = 'grep -A3 -f ' + ccIDFile + '2' + ' ' + fqfile + " | sed -e 's/^--$//g' | sed '/^$/d'" + ' > ' + outfile
    os.system(grep) 

    sort = 'sort ' + accIDFile + " | uniq > " + accIDFile +"2"
    os.system(sort)

# THEN: Not necessary
# grep -A 3 TS211A-2-6-19.dust.rmdup.fq -f Streptophyta.ReadIDs2 > Streptophyta.fq




main()
