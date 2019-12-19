#!/usr/bin/env python
# FIXME this one needs refactoring and factoring out the hardcoded paths
# FIXME normalize names and spaces

from __future__ import print_function
import re
import sys
import glob
import argparse
from os.path import basename, join, dirname
from collections import defaultdict


def get_overlaping_bins(dict_cogbin_unitigs, cog_threshold, overlap_threshold):
    set_bins_global = {key for dict_ in dict_cogbin_unitigs.values() for key in dict_}
    dict_bins_common_cogs = defaultdict(list)
    for Cog, dict__bin_unitig in dict_cogbin_unitigs.items():
        for index, (bin1, set1) in enumerate(list(dict__bin_unitig.items())[:-1]):
            for bin2, set2 in list(dict__bin_unitig.items())[index+1:]:
                if set1 & set2:
                    if max(len(set1 & set2)/float(len(set1)), len(set1 & set2)/float(len(set2))) >= overlap_threshold:
                        dict_bins_common_cogs[tuple(sorted([bin1, bin2]))].append(Cog)

    # Summarize for each bin how many cogs are shared
    dict_bin_cogs = defaultdict(set)
    for (bin1, bin2), list_cogs in dict_bins_common_cogs.items():
        dict_bin_cogs[bin1] |= set(list_cogs)
        dict_bin_cogs[bin2] |= set(list_cogs)

    # So which bins should be merged and which should just be flagged
    candidate_to_merge = {}
    dict_to_flag = {}
    for bins in set_bins_global:
        list_cog = dict_bin_cogs[bins]
        if len(list_cog) < cog_threshold:
            dict_to_flag[bins] = list_cog
        if len(list_cog) >= cog_threshold:
            candidate_to_merge[bins] = list_cog

    # list bins to merge
    list_sets_tomerge = [{bin1, bin2} for (bin1, bin2) in dict_bins_common_cogs.keys() if (bin1 in candidate_to_merge) and (bin2 in candidate_to_merge)]

    # take into accounts bins with too many shared COGs but not going to be merged 
    bins_going_to_merge={bins for set_bin in list_sets_tomerge for bins in set_bin}
    for bins in candidate_to_merge.keys():
        if bins not in bins_going_to_merge:
            dict_to_flag[bins] = dict_bin_cogs[bins]

    # deal with the possibility that more than 2 bins must be merged together.
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

    # name bins from bin_to_merge, in a way that names don't collide if get_overlaping_bins is run multiple time
    dict_merge_bins={}
    for list_bins in list_sets_tomerge:
        index=1
        new_name="Bin_merged_"+str(index)
        while new_name in set_bins_global:
            index+=1
            new_name="Bin_merged_"+str(index)
        set_bins_global.add(new_name)
        dict_merge_bins[new_name] = list_bins
    return dict_to_flag, dict_merge_bins


def update_cogbin_unitigs(dict_merge_bins,dict_cogbin_unitigs):
    # remove merged bins
    Set_bins_merged={bin for Tuple in dict_merge_bins.values() for bin in Tuple}
    dict_cogbin_unitigs_merged = {cog: {bin_: unitigs for bin_, unitigs in dict_bin_unitig.items() if bin_ not in Set_bins_merged} for cog, dict_bin_unitig in dict_cogbin_unitigs.items()}

    # add merged bins in the datastructure
    for cog in dict_cogbin_unitigs_merged:
        for name, list_bins in dict_merge_bins.items():
            dict_cogbin_unitigs_merged[cog][name] = set.union(*[dict_cogbin_unitigs[cog][bins] for bins in list_bins])
    return dict_cogbin_unitigs_merged


def main(cog_threshold, bins_to_merge, cogs_to_ignore, bins_to_process, rel_path, overlap_threshold):
    dict_cogbin_unitigs = defaultdict(lambda: defaultdict(set))

    for bin_path in bins_to_process:
        bin_ = basename(bin_path)
        for cog_file in glob.glob(join(bin_path,rel_path)) :
            cog = basename(cog_file).replace(".gfa","")
            test = {line.split("\t")[2] for line in open(cog_file) if line[0] == "S"}
            dict_cogbin_unitigs[cog][bin_] = test

    # find out which bin needs to be merged or have cogs to flags
    dict_to_flag, dict_merge_bins = get_overlaping_bins(dict_cogbin_unitigs, cog_threshold, overlap_threshold)

    # add the merged bin in this datastructure before checking again if they share cogs or need to be merged
    flag=True
    dict_merg_bins_global={}
    dict_merg_bins_global.update(dict_merge_bins)
    while flag==True:
        dict_cogbin_unitigs_merged=update_cogbin_unitigs(dict_merge_bins,dict_cogbin_unitigs)
        dict_to_flag_new, dict_merge_bins_new = get_overlaping_bins(dict_cogbin_unitigs_merged, cog_threshold, overlap_threshold)
        if dict_merge_bins_new!={}:
            # merge the dict_merge_bins
            to_forget=[]
            for bins,composition in dict_merge_bins_new.items():
                real_composition=[]
                for comp_bin in composition :
                    if comp_bin in dict_merg_bins_global:
                        real_composition+=dict_merg_bins_global[comp_bin]
                        to_forget.append(comp_bin) 
                    else:
                        real_composition.append(comp_bin)
                dict_merg_bins_global[bins]=real_composition
            for bins in to_forget :
                del dict_merg_bins_global[bins]
        else :
            flag=False

    # output cog to be ignored
    with open(bins_to_merge, "w") as out:
      #TODO list(List) looks very weird!
        out.write("\n".join(["\t".join([key]+list(List))for key, List in dict_merg_bins_global.items()]))
    # output bins to be merged
    with open(cogs_to_ignore, "w") as out:
      #TODO list(List) looks very weird!
        out.write("\n".join(["\t".join([key]+list(List))for key, List in dict_to_flag_new.items()]))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", nargs='+', help="list of all bin considered",required=True)    
    parser.add_argument("-g", help="relative path of graph files from the Bin directory, must end in COG*.gfa",required=True)
    parser.add_argument("cog_threshold", help="number of cogs in common between bins to require merging",default="10")
    parser.add_argument("bins_to_merge", help="Output merge plan (.tsv)")
    parser.add_argument("cogs_to_ignore", help="Output bin cogs to ignore (.tsv)")
    parser.add_argument("-t", help="overlap treshold, percent of unitigs shared between graphs, to consider the graphs shared by multiple bins",default='0.1')
    args = parser.parse_args()
    
    main(int(args.cog_threshold), args.bins_to_merge, args.cogs_to_ignore, args.b, args.g, float(args.t))
