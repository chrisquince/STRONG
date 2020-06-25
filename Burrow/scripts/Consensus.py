#!/usr/bin/env python
from __future__ import print_function
import argparse
from collections import Counter, defaultdict
import sys

def main(bin_file):
    with open(bin_file) as handle:
        _ = next(handle)
        dict_contig_bins = {}
        for line in handle:
            contig, bin = line.rstrip().split(",")
            dict_contig_bins[contig] = bin
    number_of_dot = Counter(map(lambda x: len(x.split(".")), dict_contig_bins.keys()))
    if len(number_of_dot) > 3:
        print("Warning the number of dot in your contigs name is varying too much, it is not possible to distinguich cuts contigs from normal ones")
        exit()
    num_split = max(number_of_dot.keys())
    dict_normal_contigs = {contig: bins for (contig, bins) in dict_contig_bins.items() if len(contig.split(".")) != num_split}
    # get splits contigs and their bins
    dict_split = defaultdict(list)
    for contig, bins in dict_contig_bins.items():
        if len(contig.split(".")) == num_split:
            contig = ".".join(contig.split(".")[:-1])
            dict_split[contig].append(bins)
    dict_split = {contig: Counter(list_bins)
                  for contig, list_bins in dict_split.items()}
    # Assign each contig to an unique bin
    dict_splitcontig_bins = {contig: max(CounterBins.items(), key=lambda x: x[1])[
        0] for contig, CounterBins in dict_split.items()}
    # Output results
    print("contig_id,0")
    for contig, bins in list(dict_splitcontig_bins.items())+list(dict_normal_contigs.items()):
        print(contig+","+bins)
    # Add some warnings in case assignment is not so obvious :
    for contig, CouterBins in dict_split.items():
        if len(CouterBins) > 1:
            Sorted_Counter = sorted(CouterBins.items(), key=lambda x: x[1])
            if Sorted_Counter[-1][1] == Sorted_Counter[-2][1]:
                print('Warning contig %s is split equally between bins %s and %s' % (
                    contig, Sorted_Counter[-1][0], Sorted_Counter[-2][0]), file=sys.stderr)
            elif Sorted_Counter[-1][1] <= 0.5*sum([counter[1] for counter in Sorted_Counter]):
                print('Warning contig %s is oversplit : less than 50%% of it is present in bin %s' % (
                    contig, Sorted_Counter[-1][0]), file=sys.stderr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "bin_file", help="Binning result file, csv, first column is the contig, second column is the bin")
    args = parser.parse_args()
    bin_file = args.bin_file
    main(bin_file)
