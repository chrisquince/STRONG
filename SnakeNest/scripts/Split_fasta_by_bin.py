#!/usr/bin/env python

from __future__ import division
import argparse,os
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import defaultdict
import resource



def  main(fasta_file,bin_composition,set_bins,output,folder):
    Dico_contigs_bin={line.rstrip().split(",")[0]:line.rstrip().split(",")[1] for line in open(bin_composition) if line.rstrip().split(",")[1] in set_bins}
    Dico_bin_Handle={}
    for bins in set(Dico_contigs_bin.values()) :
        if folder :
            os.system("mkdir -p "+output+"/Bin_"+bins)
        Dico_bin_Handle[bins]=open(output+("/Bin_"+bins+"/")*folder+"/Bin_"+bins+"."+fasta_file.split(".")[-1],"w")
    for contig_id,seq in SimpleFastaParser(open(fasta_file)) :
        contig_id2=contig_id.split()[0]
        if contig_id2 in Dico_contigs_bin :
            Dico_bin_Handle[Dico_contigs_bin[contig_id2]].write(">"+contig_id+"\n"+seq+"\n")
    for handle in Dico_bin_Handle.values() :
        handle.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_file", help="fasta file containing sequences binned")  
    parser.add_argument("bin_composition", help="csv file giving bin composition, output of concoct, first column=contig id, second is bin number")
    parser.add_argument("path_output", help="specify the place you want all these bin folder to be put")
    parser.add_argument("-l", nargs='+',help="restrict to a list of specifics bins ")
    parser.add_argument("-d",help=" actually create a folder for each bin", action='store_true')
    args = parser.parse_args()
    fasta_file=args.fasta_file
    bin_composition=args.bin_composition
    path_output=args.path_output
    if args.d :
        folder=True
    else :
        folder=False
    #---------- Batch strategy ----------------
    # there is a hard limit of number of handles which can be opened
    max_nb_handles=resource.getrlimit(resource.RLIMIT_NOFILE)[0]
    max_nb_handles=max_nb_handles-100
    List_all_bins=list({line.rstrip().split(",")[1] for line in open(bin_composition)})
    Nb_bins=len(List_all_bins)
    nb_batch=Nb_bins//max_nb_handles+1
    List_Batchs=[List_all_bins[batch*max_nb_handles:min((batch+1)*max_nb_handles,Nb_bins)] for batch in range(nb_batch)]
    for Batch in List_Batchs :
        if args.l :
            set_bins=set(Batch)&set(args.l)
        else :
            set_bins=set(Batch)
        main(fasta_file,bin_composition,set_bins,path_output,folder)
