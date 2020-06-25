#!/usr/bin/env python

from Bio import SeqIO
import sys, getopt
import os
import argparse

def main(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument("fasta_file", help="fasta input file")

    args = parser.parse_args()

    handle = open(args.fasta_file, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        seq = record.seq
    
        print(record.id + '\t0\t' + str(len(seq)))

if __name__ == "__main__":
    main(sys.argv[1:])
