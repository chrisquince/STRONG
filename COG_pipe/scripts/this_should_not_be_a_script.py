#!/usr/bin/env python
from os.path import basename, dirname
import argparse
import glob

parser = argparse.ArgumentParser()
parser.add_argument("mags")
parser.add_argument("merge_recipe")
parser.add_argument("output")
parser.add_argument("path_to_bins")
parser.add_argument("nb_orfs")

args = parser.parse_args()

#import ipdb; ipdb.set_trace()
bin_files = configfiles = glob.glob(args.path_to_bins + '/*/selected_cogs.tsv')   

bins = [basename(dirname(cogs_fn)) for cogs_fn in bin_files if len(list(open(cogs_fn).readlines())) >= int(args.nb_orfs)]
# we have been considering indistinctly bins and mags previously, we will now only retain mags
mags = set(["Bin_%s" % line.rstrip() for line in open(args.mags)])
merged_recipe = {line.rstrip().split('\t')[0]:line.rstrip().split('\t')[1:] for line in open(args.merge_recipe)}

for merged_bin, bins_to_merge in merged_recipe.items():
    
    mags.add(merged_bin)

    for bin_to_merge in bins_to_merge:
        mags.discard(bin_to_merge)


bins_selected=list(mags &set(bins))
with open(args.output, 'w') as out:
  out.write("\n".join(bins_selected))
