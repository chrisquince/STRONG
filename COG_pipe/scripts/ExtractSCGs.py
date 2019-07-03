#!/usr/bin/env python

import argparse
import sys
import glob
import os

from collections import defaultdict

from Bio.SeqIO.FastaIO import SimpleFastaParser as SPF

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
    for filename in glob.iglob(args.root_dir + "/*/SCG.fna"):


        bin_name = filename.split("/")[-2]

        # print(bin_name)
        bins.add(bin_name)
        for header,seq in SPF(open(filename)):
            cogId=header.split(" ")[1]

            if cogId in core_cogs:
                # print("%s %i" % (record.id, len(record)))
            
                if bin_name in corecog_bin[cogId]:
                    if len(seq) > len(corecog_bin[cogId][bin_name]):
                        corecog_bin[cogId][bin_name] = [header,seq]
                else:
                    corecog_bin[cogId][bin_name] = [header,seq]

    path_output=args.output_folder+"/SCGs/"
    if not os.path.isdir(path_output) :
        os.system("mkdir -p "+path_output)


    for corecog,dict_List_bins in corecog_bin.items():
        with open(path_output+ corecog + ".fna", "w") as output_handle:
            for bin_name,(header,seq) in dict_List_bins.items() :
                new_name=bin_name + "_" + corecog+" "+header
                output_handle.write(">"+new_name+"\n"+seq+"\n")


if __name__ == "__main__":
    main(sys.argv[1:])
