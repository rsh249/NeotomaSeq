#!/usr/bin/python3

#Import from biopython library (included in anaconda)
from Bio import SeqIO
import sys
#from multiprocessing import Pool
#from itertools import product
import os
    
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
    #Read in Tax.Results file './PPC-Pr.analysis/Tax.Results'
    with open (str(sys.argv[1]), 'rt') as orderfile:  
        for line in orderfile:
            linesTaxID.append(line)
    #Specific classification of interest
    substr = "Rodentia"
    rodentID = []
    #Pull out taxonomic IDs of interest based on classification level
    for line in linesTaxID[1: ]:
        order = line.split("\t")[2]
        taxID = int(line.split("\t")[6])
        if order == substr:
            if taxID not in rodentID:
                rodentID.append(taxID)
            else:
                continue
        else:
            continue   
    linesReadID = []
    #Read in .centclass file './PPC-Pr.analysis/PPCPr3kA.centsmall'
    with open (str(sys.argv[2]), 'rt') as readfile:  
        for line in readfile:
            linesReadID.append(line)
    fileP = str(sys.argv[2]).split("/")
    name = fileP[-1].split(".")
    ccIDFile = name[0] + ".ReadIDs"
    ccIDs = open(ccIDFile, "w")
    #Pull out read IDs that align with the taxonomic ID
    for line in linesReadID[1: ]:
        try:
            readID = line.split('\t')[0]
            taxRead = int(line.split("\t")[2])
        except:
            continue
        #Remove the 'M_' that is in front of the read ID in .centclass but not in fastq
        readID = readID[:0] + readID[2:]
        for ID in rodentID:
            if taxRead == ID:
                ccIDs.write(readID)
                ccIDs.write('\n')
            else:
                continue
    ccIDs.close()
    sort = 'sort ' + ccIDFile + " | uniq > " + ccIDFile + "2"
    os.system(sort)
    
    os.system("LC_ALL=C")
    outfile = name[0] + ".filtered"
    fqfile = sys.argv[3]
    grep = 'grep -A3 -f ' + ccIDFile + '2' + ' ' + fqfile + " | sed -e 's/^--$//g' | sed '/^$/d'" + ' > ' + outfile
    os.system(grep) 




    
    #print(readIDs)
##    #parse the fastq file and run the filter above to write a file of lines of interest './PPC-Pr.analysis/PPCPr3kA.smallfast'
##    with open( str(sys.argv[3]) ) as fqfile:
##        parser = SeqIO.parse(fqfile, "fastq")
##        pool = Pool(processes=4)
##        SeqIO.write(pool.starmap(filter, product(parser, readIDs)), 'Filtered_fastq', "fastq")

main()
