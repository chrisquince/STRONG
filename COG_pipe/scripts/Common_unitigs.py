#!/usr/bin/env python
# FIXME this one needs refactoring and factoring out the hardcoded paths
# FIXME normalize names and spaces

from __future__ import print_function
import re
import sys
import glob
import argparse
from collections import defaultdict


def get_overlaping_bins(Dico_CogBin_unitigs, cog_threshold):
    Set_bins = {key for dico in Dico_CogBin_unitigs.values() for key in dico}
    Dico_bins_common_cogs = defaultdict(list)
    for Cog, dico_bin_unitig in Dico_CogBin_unitigs.items():
        for index, (bin1, set1) in enumerate(list(dico_bin_unitig.items())[:-1]):
            for bin2, set2 in list(dico_bin_unitig.items())[index+1:]:
                if set1 & set2:
                    if max(len(set1 & set2)/float(len(set1)), len(set1 & set2)/float(len(set2))) >= 0.1:
                        Dico_bins_common_cogs[tuple(
                            sorted([bin1, bin2]))].append(Cog)
    # Summarize for each bin how many cogs are shared
    Dico_bin_COGs = defaultdict(set)
    for (cog1, cog2), list_cogs in Dico_bins_common_cogs.items():
        Dico_bin_COGs[cog1] |= set(list_cogs)
        Dico_bin_COGs[cog2] |= set(list_cogs)

    # So which bins should be merged and which should just be flagged
    Dico_to_merge = {}
    Dico_to_flag = {}
    for bins in Set_bins:
        list_cog = Dico_bin_COGs[bins]
        if len(list_cog) < cog_threshold:
            Dico_to_flag[bins] = list_cog
        if len(list_cog) >= cog_threshold:
            Dico_to_merge[bins] = list_cog

    # list bins to merge
    List_sets_tomerge = [set([bin1, bin2]) for (bin1, bin2) in Dico_bins_common_cogs.keys(
    ) if (bin1 in Dico_to_merge) and (bin2 in Dico_to_merge)]
    Len = 0
    while len(List_sets_tomerge) != Len:
        New_List_sets_tomerge = []
        Len = len(List_sets_tomerge)
        for Set_bins in List_sets_tomerge:
            intersect = 0
            for index, element in enumerate(New_List_sets_tomerge):
                if Set_bins & element:
                    intersect = 1
                    New_List_sets_tomerge[index] |= Set_bins
            if not intersect:
                New_List_sets_tomerge.append(Set_bins)
        List_sets_tomerge = New_List_sets_tomerge
    Dico_merge_bins = {
        "Bin_merged_"+str(index+1): list_bins for index, list_bins in enumerate(List_sets_tomerge)}
    return Dico_to_flag, Dico_merge_bins


def main(gfa_regex, cog_threshold, output):
    # dangerous but needed
    List_gfa_files = [files for files in glob.glob(gfa_regex)]
    Dico_CogBin_unitigs = defaultdict(lambda: defaultdict(set))
    for file in List_gfa_files:
        COG = file.split("/")[-1].split(".gfa")[0]
        # more robust but still inherently not robust to naming changes
        Bin = re.findall(".*?Bin.*?(Bin.+?)/", file)[0]
        Dico_CogBin_unitigs[COG][Bin] = {line.split(
            "\t")[2] for line in open(file) if line[0] == "S"}

    Dico_to_flag, Dico_merge_bins = get_overlaping_bins(
        Dico_CogBin_unitigs, cog_threshold)

    # check for each of the merge list if any cog need to be deleted :
    Dico_CogBin_unitigs_merged = {COG: {Bin: unitigs for Bin, unitigs in dico_bin_unitig.items() if Bin not in {
        bin for Tuple in Dico_merge_bins.values() for bin in Tuple}} for COG, dico_bin_unitig in Dico_CogBin_unitigs.items()}
    for COG in Dico_CogBin_unitigs_merged:
        for name, list_bins in Dico_merge_bins.items():
            Dico_CogBin_unitigs_merged[COG][name] = set().union(
                *[Dico_CogBin_unitigs[COG][bins] for bins in list_bins])

    Dico_to_flag2, Dico_merge_bins2 = get_overlaping_bins(
        Dico_CogBin_unitigs_merged, cog_threshold)
    if Dico_merge_bins2:
        print("code contains errors : first pass merging was not enough to merge all bins which should be merged ", file=sys.stderr)
        exit(1)
    # output cog to be ignored
    with open(output+"/List_bin_cogs_to_ignore.tsv", "w") as Handle:
        Handle.write("\n".join(["\t".join([key]+list(List))
                                for key, List in Dico_to_flag2.items()]))
    # output bins to be merged
    with open(output+"/List_bin_to_merge.tsv", "w") as Handle:
        Handle.write("\n".join(["\t".join([key]+list(List))
                                for key, List in Dico_merge_bins.items()]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "gfa_regex", help='expect a regular expression matching all cog simplified subgraphs of all bins, be advised that bins folders must explicitely start by "Bin".  example : subgraphs/Merged_Bin/Bin*/simplif/COG*.gfa')
    parser.add_argument(
        "cog_threshold", help="number of cogs in common between bins to require merging", default="10")
    parser.add_argument("output", help="output file")
    args = parser.parse_args()
    gfa_regex = args.gfa_regex
    cog_threshold = int(args.cog_threshold)
    output = args.output
    main(gfa_regex, cog_threshold, output)
