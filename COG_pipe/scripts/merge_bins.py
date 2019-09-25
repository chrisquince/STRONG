#!/usr/bin/env python
import argparse 
import glob
import shutil
import os

def make_link(source_root, subfolder, target_root):
    relpath = os.path.relpath(source_root, target_root)
    os.symlink(os.path.join(relpath, subfolder), os.path.join(target_root, subfolder))

def make_copy(source_root, subfolder, target_root):
    shutil.copytree(os.path.join(source_root, subfolder), os.path.join(target_root, subfolder))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("merge_plan", help="File describing the bins to be merged")
    parser.add_argument("bins_dir", help="Folder with the original bins (NB: should be relative to the root!)")
    parser.add_argument("contig_assign", help="File with original contigs assignment")
    parser.add_argument("out_dir", help="Output folder (NB: should be relative to the root!)")
    args = parser.parse_args()

    os.system('rm -rf %s/Bin_*' % args.out_dir)

    #FIXME can we just consider all the subfolders?
    all_bins = [os.path.dirname(path) for path in glob.glob("%s/Bin_*/SCG.fna" % args.bins_dir)]
    all_bins = [os.path.basename(p) for p in all_bins]
    merged_bins = set()
    name_map = dict()
    with open(args.merge_plan) as merge_plan:
        for line in merge_plan:
            #TODO just split()?
            split_line=line.rstrip().split("\t")
            if len(split_line)>1 : # so that no empty folder is created, making os.symlink fail
                merged_name = split_line[0]
                to_merge = split_line[1:]
                merged_bins |= set(to_merge)
                merged_path = os.path.join(args.out_dir, merged_name)
                os.makedirs(merged_path)
                for b in to_merge:
                    bin_dir = os.path.join(args.bins_dir, b)
                    #os.system("cat %s/contigs.fa >> %s/contigs.fa" % (bin_dir, merged_path))
                    os.system("cat %s/SCG.fna >> %s/SCG.fna" % (bin_dir, merged_path))
                    name_map[b.replace("Bin_","")] = merged_name.replace("Bin_","")
  
    # link all non merged bin in the merged bin folder
    for b in set(all_bins) - merged_bins:
        #os.system("ln -s %s/%s %s" % (args.bins_dir, b, args.out_dir))
        # make_link(args.bins_dir, b, args.out_dir)
        make_copy(args.bins_dir, b, args.out_dir)

    # create a new contig assignment file 
    with open(args.out_dir + "/clustering.csv", "w") as out:
        with open(args.contig_assign) as contig_assign:
            for line in contig_assign:
                #TODO remove rstrip?
                contig, name = line.rstrip().split(",")
                new_name = name_map[name] if name in name_map else name
                out.write("%s,%s\n" % (contig, new_name))
