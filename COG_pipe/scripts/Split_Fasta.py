#!/usr/bin/env python
# -*- coding: latin-1 -*-
from Bio.SeqIO.FastaIO import *
from collections import defaultdict,Counter
import numpy as np
import argparse
import os 

def main(fasta_file,nb_chunks,output) :
	Sorted_Names=[]
	Dico_genome_seq={}
	Dico_genome_len={}
	for header,seq in SimpleFastaParser(open(fasta_file)) :
		Sorted_Names.append(header)
		Dico_genome_seq[header]=seq
		Dico_genome_len[header]=len(seq)
	Total_length=sum(Dico_genome_len.values())
	Chunk_size=Total_length/float(nb_chunks)
	os.system("mkdir -p "+output)
	# Start of loop
	num=0
	extension="."+fasta_file.split(".")[-1]
	fasta_path=output+"/"+fasta_file.split('/')[-1].split(".")[0]
	Current_filename=lambda x:fasta_path+"_"+str(x)+extension
	Handle=open(Current_filename(num),"w")
	Temp_length=0
	for header in Sorted_Names :
		if Temp_length>Chunk_size :
			Temp_length=0
			num+=1
			Handle.close()
			Handle=open(Current_filename(num),"w")
		Seq=Dico_genome_seq[header]
		Temp_length+=len(Seq)
		Handle.write(">"+header+"\n"+Seq+"\n")
	Handle.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("Fasta_file", help="fasta file you want to split in chuncks")
    parser.add_argument("nb_chunks", help="nb of Chunks")
    parser.add_argument("-o", help="Name of the folder where you want to store the chunks",default="Split_")
    args = parser.parse_args()
    fasta_file=args.Fasta_file
    n=int(args.nb_chunks)
    output=args.o
    if output[-1]=="/" :
    	output=output[:-1]
    if output=="Split_" :
    	output+=".".join(fasta_file.split(".")[:-1])
    main(fasta_file,n,output)
