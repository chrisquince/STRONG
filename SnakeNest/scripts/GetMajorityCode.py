#!/usr/bin/env python3
from Bio import SeqIO
import sys
import argparse
import os

def main(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument("table_file", help="prodigal gff file")

    parser.add_argument("fna_file", help="prodigal ORF sequences")

    args = parser.parse_args()

    #import ipdb; ipdb.set_trace() 

    map_code = {}

    with open(args.table_file) as fin:
        for line in fin:
            line = line.rstrip()
            toks = line.split('\t')
            map_code[toks[0]] = toks[1]

    freq11 = 0
    freq4 = 0

    handle = open(args.fna_file, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        contig = '_'.join(record.id.split('_')[:-1])

        if contig in map_code:

            if map_code[contig] == '4':
                freq4 += 1
            if map_code[contig] == '11':
                freq11 += 1

    if freq11 >= freq4:
        print('code=11',end='')
    else:
        print('code=4',end='')

if __name__ == "__main__":
    main(sys.argv[1:])
