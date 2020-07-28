#!/usr/bin/env python
import os
import glob
import numpy as np 
from collections import defaultdict


def format_row(header,line):
    row = header
    for el in line:
        row += "\t{:.2g}".format(el)
    return row+"\n"
def matrix_write(matrix,file_name,col_names,row_names):
    with open(file_name,"w") as handle:
        handle.write("/\t%s\n"%"\t".join(col_names))
        handle.writelines(format_row(row_names[index],line) for index,line in enumerate(matrix))

def main(cluster,coverage,bed,output):
    # get list of contigs corresponding ot selected mags 
    mags = {mag.replace("Bin_","") for mag in ALL_MAGS}
    mag_to_contigs = defaultdict(list)
    with open(cluster) as handle:
        _=next(handle)
        for line in handle:
            contig,mag = line.rstrip().split(",")
            if mag in mags:
                mag_to_contigs["Bin_%s"%mag].append(contig)
    # translate contig to split contigs, so that we don't need to get cov for non split contigs
    contigs = {contig for cntgs in mag_to_contigs.values() for contig in cntgs}
    contigs_to_splits = defaultdict(list)
    split_to_len = {}
    for line in open(bed):
        contig,start,end,split = line.rstrip().split("\t")
        if contig in contigs:
            split_to_len[split]=np.abs(int(end)-int(start))
            contigs_to_splits[contig].append(split)
    # read coverage file and only look at selected splits
    # in the mean time, fill a matrix of number of nucleotide per mag
    splits_to_mag = {split:mag for mag,contigs in mag_to_contigs.items() for contig in contigs for split in contigs_to_splits[contig]}
    mag_to_len = defaultdict(int)
    with open(coverage) as handle:
        sorted_samples = next(handle).rstrip().split('\t')[1:]
        sorted_mags = sorted(ALL_MAGS)
        nuc_matrix = np.zeros((len(sorted_mags),len(sorted_samples)))
        for line in handle:
            splitline = line.rstrip().split("\t")
            if splitline[0] in splits_to_mag:
                # get things
                mag = splits_to_mag[splitline[0]]
                row_index = sorted_mags.index(mag)
                split_len = split_to_len[splitline[0]]
                profile = split_len*np.array([float(val) for val in splitline[1:]])
                # store things
                nuc_matrix[row_index,:] += profile
                mag_to_len[mag] += split_len
    # get mag coverage, from mag mapped nucleotides
    sorted_mag_len = np.array([mag_to_len[mag] for mag in sorted_mags])
    cov_matrix = (nuc_matrix.T/sorted_mag_len).T
    # output result:
    matrix_write(cov_matrix,output,sorted_samples,sorted_mags)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("cluster", help="concoct cluster file")
    parser.add_argument("coverage", help="split contigs coverage")
    parser.add_argument("bed", help="bed file defining split contigs relative to contigs")
    parser.add_argument("output")    
    args = parser.parse_args()
    main(args.cluster,args.coverage,args.bed,args.output)
