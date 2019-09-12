#!/usr/bin/env python
# -*- coding: latin-1 -*-
import os
import argparse
import numpy as np
from collections import Counter, defaultdict
from subprocess import Popen, PIPE

# FIXME normalize spaces, local variable names, handles, scpacec, absolute indexing in paths


def prodigal_gff_parser(Handle):
    # specific to prodigal output, does not seem to be universal gff
    # directly adapted from  SimpleFastaParser in Bio.SeqIO.FastaIO
    while True:
        line = Handle.readline()
        if line == "":
            return
        if line[:16] == "# Sequence Data:":
            break
    while True:
        if not line:
            break
        if line[:16] != "# Sequence Data:":
            print(line)
            raise ValueError(
                "GFF Output from prodigal should start with '# Sequence Data:'")
        seq_data = line.rstrip()
        Model_data = Handle.readline().rstrip()
        line = Handle.readline().rstrip()
        ORF_list = []
        if not line:
            break
        while line[0] != '#':
            if not line:
                break
            ORF_list.append(line)
            line = Handle.readline().rstrip()
            if not line:
                break
        yield seq_data, Model_data, ORF_list
    if not line:
        return  # StopIteration


def gff_to_bed(gff_file, bed_file):
    to_write = []
    with open(gff_file, 'r') as gff:
        for seq_data, Model_data, ORF_list in prodigal_gff_parser(gff):
            for index_orf, ORF in enumerate(ORF_list):
                contig = ORF.split()[0]
                start = ORF.split()[3]
                end = ORF.split()[4]
                to_write.append("\t".join(
                    [contig, start, end, contig+"_"+ORF.split()[-1].split(";")[0].split("_")[1]]))

    with open(bed_file, 'w') as bed:
        # TODO premature pissimisation
        bed.write("\n".join(to_write))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("gff", help="Input .gff file (e.g. prodigal output)")
    parser.add_argument("bed", help="Output .bed file")
    args = parser.parse_args()

    gff_to_bed(args.gff, args.bed)
