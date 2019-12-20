#!/usr/bin/env python
from os.path import basename, dirname
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-s", nargs='+',required=True)    
parser.add_argument("mags")
parser.add_argument("merge_recipe")
parser.add_argument("output")
parser.add_argument("nb_orfs")

args = parser.parse_args()


bins = [basename(dirname(cogs_fn)) for cogs_fn in args.s if len(list(open(cogs_fn).readlines())) >= int(args.nb_orfs)]
# we have been considering indistinctly bins and mags previously, we will now only retain mags
mags = ["Bin_%s" % line.rstrip() for line in open(args.mags)]
merged_recipe = {line.rstrip().split('\t')[0]:line.rstrip().split('\t')[1:] for line in open(args.merge_recipe)}
mags += [merged_bin for mergerd_bin,list_bins in merged_recipe.items() if (set(mags)&set(merged_recipe))!=set()]
bins_selected=list(set(mags)&set(bins))
with open(args.output, 'w') as out:
  out.write("\n".join(bins_selected))
