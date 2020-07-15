#!/usr/bin/env python
# -*- coding: latin-1 -*-
import os
import argparse
import itertools
import numpy as np
from os.path import basename
from collections import Counter, defaultdict
from subprocess import Popen, PIPE


def collate_iterator(cov_files, output):
    get_cov = lambda line : str(float(line.rstrip().split("\t")[4])) # change "0.0000000" to "0.0"
    get_feature = lambda line : line.rstrip().split("\t")[3]
    sample_to_cov = {".".join(basename(path).split('.')[:-2]):path for path in cov_files}
    sorted_samples = sorted(sample_to_cov.keys())
    sorted_cov_files = [sample_to_cov[sample] for sample in sorted_samples]
    List_handle = [open(file) for file in sorted_cov_files]
    with open(output, "w") as handle:
        handle.write("\t".join(["cov"]+sorted_samples)+"\n")
        handle.writelines("%s\t%s\n"%(get_feature(line[0]),"\t".join(list(map(get_cov, line)))) for line in zip(*List_handle))
    for handle in List_handle:
        handle.close()


def collate_in_memory(cov_files, output):
    # we know that all coverage are sorted in the same fashion since they originiate from the same .bed file, then let's drop the dictionary, we still need to store everything in memory, mainly since I do have more than 1024 files.
    sample_to_cov = {".".join(basename(path).split('.')[:-2]):path for path in cov_files}
    sorted_samples = sorted(sample_to_cov.keys())
    sorted_feature = [line.rstrip().split("\t")[3] for line in open(cov_files[0])]
    sorted_cov_files = [sample_to_cov[sample] for sample in sorted_samples]
    for file in sorted_cov_files:
        for index_row, line in enumerate(open(file)):
            sorted_feature[index_row] += "\t"+str(float(line.rstrip().split("\t")[4]))
    sorted_feature = ["\t".join(["cov"]+sorted_samples)]+sorted_feature
    with open(output, "w") as handle:
        handle.writelines(line+"\n" for line in sorted_feature)


def collate_dictionary(cov_files, output):
    sorted_files = sorted(cov_files)
    samples_names = sorted([".".join(basename(path).split('.')[:-2]) for path in sorted_files])
    sorted_feature = [line.rstrip().split("\t")[3] for line in open(sorted_files[0])]
    cov_matrix = np.zeros((len(sorted_feature), len(samples_names)))
    # the same bed file is used to generate all these coverage files. In consequence the features are ordered in the same way in each file
    # TOFIX : maybe not optimal way of storgin/writting this data.
    for index_col, file in enumerate(sorted_files):
        print(file)
        for index_row, line in enumerate(open(file)):
            feature, cov = line.rstrip().split("\t")[3:]
            cov_matrix[index_row, index_col] = str(float(cov))
    # output :
    with open(output, "w") as handle:
        handle.write("\t".join(["feature"]+samples_names)+"\n")
        for index, line in enumerate(cov_matrix):
            handle.write("\t".join([sorted_feature[index]]+line)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", required=True, help="output file name")
    parser.add_argument("-l", nargs="+", required=True,
                        help="list of coverage files")
    args = parser.parse_args()
    output = args.o
    cov_files = args.l
    if len(cov_files) < 1020:
        # TOFIX : I could do a hierarchical collate : collate by batch of 1024 and collate the resulting collated files. 
        # it would be more disk intensive and less ram intensive than collate_in_memory. Not sure if faster in general.
        collate_iterator(cov_files, output)
    else:
        collate_in_memory(cov_files, output)
