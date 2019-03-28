#!/usr/bin/env python

import argparse
import sys
import glob
import os

from collections import defaultdict

from Bio import SeqIO

def main(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument("root_dir", help="directory with bin dirs in")

    parser.add_argument("cogs_list", help="list of COGs to collate from bins")

    parser.add_argument("output_folder", help="output folder name, where do we want to put the sequences")

    args = parser.parse_args()

#    import ipdb; ipdb.set_trace()

    with open(args.cogs_list, 'r') as f:
        core_cogs = set([x.rstrip() for x in f.readlines()])


    corecog_bin = defaultdict(dict)
    bins = set()
    for filename in glob.iglob(args.root_dir + "/Bin_*/SCG.fna"):
        filename.rstrip()

        bin_name = filename.split("/")[-2]

        # print(bin_name)
        bins.add(bin_name)
        for record in SeqIO.parse(filename, "fasta"):
            cogId = record.description.split(' ')[1]

            if cogId in core_cogs:
                # print("%s %i" % (record.id, len(record)))
            
                if bin_name in corecog_bin[cogId]:
                    if len(record) > len(corecog_bin[cogId][bin_name]):
                        corecog_bin[bin_name][cogId] = record
                else:
                    corecog_bin[cogId][bin_name] = record

    path_output=args.output_folder+"/SCGs/"
    if not os.path.isdir(path_output) :
        os.system("mkdir -p "+path_output)


    for corecog in core_cogs:
        with open(path_output+ corecog + ".fna", "w") as output_handle:
            for bin_name in bins:
                if bin_name in corecog_bin[corecog]:
                    record_rename = record
                    record_rename.id = bin_name + "_" + corecog
                    SeqIO.write(record_rename, output_handle, "fasta")


if __name__ == "__main__":
    main(sys.argv[1:])
