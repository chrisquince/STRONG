#!/usr/bin/env python
# FIXME this one needs refactoring and factoring out the hardcoded paths
# FIXME normalize names and spaces

from __future__ import print_function
import re
import sys
import glob
import argparse
from collections import defaultdict


def get_overlaping_bins(dict_cogbin_unitigs, cog_threshold):
    set_bins = {key for dict_ in dict_cogbin_unitigs.values() for key in dict_}
    dict_bins_common_cogs = defaultdict(list)
    for Cog, dict__bin_unitig in dict_cogbin_unitigs.items():
        for index, (bin1, set1) in enumerate(list(dict__bin_unitig.items())[:-1]):
            for bin2, set2 in list(dict__bin_unitig.items())[index+1:]:
                if set1 & set2:
                    if max(len(set1 & set2)/float(len(set1)), len(set1 & set2)/float(len(set2))) >= 0.1:
                        dict_bins_common_cogs[tuple(
                            sorted([bin1, bin2]))].append(Cog)
    # Summarize for each bin how many cogs are shared
    dict_bin_cogs = defaultdict(set)
    for (cog1, cog2), list_cogs in dict_bins_common_cogs.items():
        dict_bin_cogs[cog1] |= set(list_cogs)
        dict_bin_cogs[cog2] |= set(list_cogs)

    # So which bins should be merged and which should just be flagged
    dict_to_merge = {}
    dict_to_flag = {}
    for bins in set_bins:
        list_cog = dict_bin_cogs[bins]
        if len(list_cog) < cog_threshold:
            dict_to_flag[bins] = list_cog
        if len(list_cog) >= cog_threshold:
            dict_to_merge[bins] = list_cog

    # list bins to merge
    #FIXME don't split lines on brackets!
    list_sets_tomerge = [set([bin1, bin2]) for (bin1, bin2) in dict_bins_common_cogs.keys(
    ) if (bin1 in dict_to_merge) and (bin2 in dict_to_merge)]
    #TODO is it to emulate do .. while?
    Len = 0
    while len(list_sets_tomerge) != Len:
        new_list_sets_tomerge = []
        Len = len(list_sets_tomerge)
        for set_bins in list_sets_tomerge:
            intersect = 0
            for index, element in enumerate(new_list_sets_tomerge):
                if set_bins & element:
                    intersect = 1
                    new_list_sets_tomerge[index] |= set_bins
            if not intersect:
                new_list_sets_tomerge.append(set_bins)
        list_sets_tomerge = new_list_sets_tomerge
    dict_merge_bins = {
        "bin__merged_"+str(index+1): list_bins for index, list_bins in enumerate(list_sets_tomerge)}
    return dict_to_flag, dict_merge_bins


def main(gfa_regex, cog_threshold, bins_to_merge, cogs_to_ignore):
    # dangerous but needed
    list_gfa_files = [files for files in glob.glob(gfa_regex)]
    dict_cogbin_unitigs = defaultdict(lambda: defaultdict(set))
    for file in list_gfa_files:
        cog = file.split("/")[-1].split(".gfa")[0]
        # more robust but still inherently not robust to naming changes
        bin_ = re.findall(".*?Bin.*?(Bin_.+?)/", file)[0]
        dict_cogbin_unitigs[cog][bin_] = {line.split(
            "\t")[2] for line in open(file) if line[0] == "S"}

    dict_to_flag, dict_merge_bins = get_overlaping_bins(
        dict_cogbin_unitigs, cog_threshold)

    #TODO improve naming
    # check for each of the merge list if any cog need to be deleted :
    dict_cogbin_unitigs_merged = {cog: {bin_: unitigs for bin_, unitigs in dict__bin_unitig.items() if bin_ not in {
        bin for Tuple in dict_merge_bins.values() for bin in Tuple}} for cog, dict__bin_unitig in dict_cogbin_unitigs.items()}
    for cog in dict_cogbin_unitigs_merged:
        for name, list_bins in dict_merge_bins.items():
            dict_cogbin_unitigs_merged[cog][name] = set().union(
                *[dict_cogbin_unitigs[cog][bins] for bins in list_bins])

    dict_to_flag2, dict_merge_bins2 = get_overlaping_bins(dict_cogbin_unitigs_merged, cog_threshold)
    if dict_merge_bins2:
        print("code contains errors : first pass merging was not enough to merge all bins which should be merged ", file=sys.stderr)
        exit(1)
    # output cog to be ignored
    with open(bins_to_merge, "w") as out:
      #TODO list(List) looks very weird!
        out.write("\n".join(["\t".join([key]+list(List))
                                for key, List in dict_to_flag2.items()]))
    # output bins to be merged
    with open(cogs_to_ignore, "w") as out:
      #TODO list(List) looks very weird!
        out.write("\n".join(["\t".join([key]+list(List))
                                for key, List in dict_merge_bins.items()]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "gfa_regex", help='expect a regular expression matching all cog simplified subgraphs of all bins, be advised that bins folders must explicitely start by "Bin_".  example : subgraphs/Merged_Bin/Bin*/simplif/COG*.gfa')
    parser.add_argument("cog_threshold", help="number of cogs in common between bins to require merging",default="10")
    parser.add_argument("bins_to_merge", help="Output merge plan (.tsv)")
    parser.add_argument("cogs_to_ignore", help="Output bin cogs to ignore (.tsv)")
    args = parser.parse_args()
    
    main(args.gfa_regex, int(args.cog_threshold), args.bins_to_merge, args.cogs_to_ignore)
