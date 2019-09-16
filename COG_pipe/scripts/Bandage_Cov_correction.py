#!/usr/bin/env python
# -*- coding: latin-1 -*-
import argparse
import os
from subprocess import Popen, PIPE
from random import randint
from Bio.SeqIO.FastaIO import *


def gfa_correction(gfa_file, kmer, output):
    def correct_KC(split_line, kmer):
        KC = next(element for element in split_line if "KC:i:" in element)
        L = len(split_line[2])
        index = split_line.index(KC)
        corected_KC = 'KC:i:'+str(float(KC.split('KC:i:')[1])*L/(L-kmer))
        split_line[index] = corected_KC
    NewGfa = ""
    for line in open(gfa_file):
        line = line.rstrip().split('\t')
        if line[0] == "S":
            correct_KC(line, kmer)
        line.append('\n')
        NewGfa += "\t".join(line)
    with open(output, 'w') as H:
        H.write(NewGfa)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("gfa_file", help="gfa file you want to correct")
    parser.add_argument("kmer", help="kmer length")
    parser.add_argument("output", help="corrected gfa file")
    args = parser.parse_args()
    gfa_file = args.gfa_file
    kmer = int(args.kmer)
    output = args.output
    gfa_correction(gfa_file, kmer, output)
