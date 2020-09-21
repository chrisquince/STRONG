#!/usr/bin/env python
import os
import glob
import argparse
import numpy as np 
from collections import defaultdict


def format_row(header,line):
    row = header
    for el in line:
        row += "\t{:.4g}".format(el)
    return row


def main(coverage,bed):

    # translate contig to split contigs, so that we don't need to get cov for non split contigs

    contigs_to_splits = defaultdict(list)
    split_to_contigs = {}
    split_to_len = {}
    for line in open(bed):
        contig,start,end,split = line.rstrip().split("\t")
        split_to_len[split]=np.abs(int(end)-int(start))
        contigs_to_splits[contig].append(split)
        split_to_contigs[split] = contig
    # read coverage file and only look at selected splits
    # in the mean time, fill a matrix of number of nucleotide per mag

    contig_to_len = defaultdict(int)
    with open(coverage) as handle:
        sorted_samples = next(handle).rstrip().split('\t')[1:]
        
        nS = len(sorted_samples)
        
        contig_cov = defaultdict(lambda: np.zeros(nS))
        
        for line in handle:
            splitline = line.rstrip().split("\t")
            if splitline[0] in split_to_contigs:
                # get things
                contig = split_to_contigs[splitline[0]]

                split_len = split_to_len[splitline[0]]
                
                contig_cov[contig] += split_len*np.array([float(val) for val in splitline[1:]])
                
                contig_to_len[contig] += split_len


    sString = "contig\t%s"%"\t".join(sorted_samples)

    print(sString)
    
    for contig,clength in contig_to_len.items():
        covP = contig_cov[contig]/float(clength)
        
        cString = format_row(contig,covP)
     
        print(cString)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("coverage", help="split contigs coverage")
    parser.add_argument("bed", help="bed file defining split contigs relative to contigs")
  
    args = parser.parse_args()
    main(args.coverage,args.bed)
