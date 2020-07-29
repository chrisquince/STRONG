#!/usr/bin/env python
import os
import argparse
import numpy as np
from os.path import basename
from collections import Counter, defaultdict
from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp

# 20 "distinguishable" colours
color_scheme = ["#0075DC","#FFCC99", "#2BCE48", "#F0A3FF", "#993F00", "#4C005C", "#808080", "#94FFB5", "#8F7C00", "#9DCC00","#C20088", "#003380", "#FFA405", "#FFA8BB", "#426600", "#FF0010", "#5EF1F2", "#00998F", "#740AFF", "#990000", "#FFFF00"]


def merge_color(colors):
    # black magic, look away (mostly padding issue though)
    def int_to_hex(x):
        hexagram = hex(int(x))[2:].upper()
        if hexagram == "0":
            return "00"
        else: 
            if len(hexagram) == 1:
                return "0"+hexagram
            else:
                return hexagram
    # convert hexa to int, get mean 
    mean_color = np.array([[int(color[1:3], 16), int(color[3:5], 16), int(color[5:], 16)] for color in colors]).mean(0)
    merged_color = "#"+"".join([int_to_hex(value) for value in mean_color])
    return merged_color

def color_gfa(gfa_file,contig_to_color):
    colored_gfa=""
    with open(gfa_file) as handle : 
        for line in handle : 
            splitline = line.rstrip().split('\t')
            # if line is vertex definition add color
            if line[0]=="S" :
                contig=splitline[1]
                if contig not in contig_to_color :
                    color="#d3d3d3"
                else :
                    color=contig_to_color[contig]
                colored_gfa+='%s\tCL:z:%s\tC2:z:%s\n'%(line.rstrip(),color,color)
            else : 
                colored_gfa+=line
    return colored_gfa

def main(gfa_file, contig_to_strains, output, split):
    sorted_strains=sorted({strain for strains in contig_to_strains.values() for strain in strains if "NA" not in strain})
    contig_to_color={contig:merge_color([color_scheme[sorted_strains.index(strain)] for strain in strains]) for contig,strains in contig_to_strains.items()}

    # next add this information to the gfa file 
    colored_gfa = color_gfa(gfa_file,contig_to_color)

    # and output 
    with open(output,"w") as handle:
        handle.write(colored_gfa)

    ### case we want 1 gfa per strain, everything else in grey 
    if split:
        for strain in sorted_strains:
            contig_to_color = {contig:color_scheme[sorted_strains.index(strain)] for contig,strains in contig_to_strains.items() if strain in strains }
            colored_gfa = color_gfa(gfa_file,contig_to_color)
            new_output = ".".join(output.split(".")[:-1])+"_%s.gfa"%strain
            with open(new_output,"w") as handle:
                handle.write(colored_gfa)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("gfa_file", help="Gfa file to color") 
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-p',help="haplotypes paths in all cogs, output of bayespaths")
    group.add_argument('-a',help="contig assignmeent file : bin/group/strain assigned to contig, first term is contig name following are list of bin/... contig it is in. Some contigs may not be assigned to anything, they are then colored grey, .tsv file")
    parser.add_argument("-s", action='store_true', help="flag, output 1 gfa per strain ")
    parser.add_argument("output",help="output name for colored gfa")
    args = parser.parse_args()
    # deal with multiple way of passing contig assignment
    if args.p : 
        # we assume the name of the gfa file contain the name of the cog of interest
        cog = basename(args.gfa_file.replace(".gfa",""))
        cog_lines = [[header]+[s.replace("-","").replace("+","") for s in seq.split(",")] for header,seq in sfp(open(args.p)) if cog in header]
        contig_to_strains = defaultdict(list)
        for line in cog_lines:
            strain = line[0].split("_")[1]
            for contig in line[1:]:
                contig_to_strains[contig].append(strain)
    else :
        contig_to_strains={line.split()[0]:line.split()[1:] for line in open(args.a)}

    main(args.gfa_file,contig_to_strains,args.output,args.s)


