#!/usr/bin/env python
# -*- coding: latin-1 -*-
from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp
from collections import defaultdict,Counter
import numpy as np
import argparse
import os 

def main(fasta_file,nb_chunks,output):
    header_len=[[header,len(seq)] for header,seq in sfp(open(fasta_file))]
    header_len = sorted(header_len,key=lambda x:-x[1] )
    # open nb_chunks handles 
    os.system("mkdir -p "+output)
    ext=fasta_file.split(".")[-1]
    fasta_path=output+"/"+fasta_file.split('/')[-1].split(".")[0]
    handles =[open("%s_%s.%s"%(fasta_path,nb,ext),"w") for nb in range(100)]
    # map header to handle, biggest to smallest seq, so it sort of even out
    index = 0
    header_to_handle={}
    for header,_ in header_len:
        header_to_handle[header] = handles[index]
        index = (index+1)%nb_chunks
    # read the file again, and this time write on all handles at the same time
    for header,seq in sfp(open(fasta_file)):
        header_to_handle[header].write(">%s\n%s\n"%(header,seq))
    # close all handles
    for handle in handles:
        handle.close()

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
