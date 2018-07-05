#!/usr/bin/python3
from Blastn import Blastn
from Bio import SeqIO
import sys

fileP = str(sys.argv[2]).split("/")
name = fileP[-1].split("_")
path =  '/'.join(fileP[0:len(fileP)-1]) + '/' + name[0] + ".fasta"
print(path)
SeqIO.convert(str(sys.argv[2]), "fastq", path, "fasta")
with Blastn(ref_path=str(sys.argv[1])) as blastn:
    hit_list = blastn(query_path=path)
    hits = open(name[0] + "hits.txt", "w")
    for hit in hit_list:
        #print (hit)
        hits.writelines(hit)
    hits.close()

