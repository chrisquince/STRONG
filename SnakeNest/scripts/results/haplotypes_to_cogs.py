#!/usr/bin/env python
# -*- coding: latin-1 -*-
import os
import argparse
from os.path import basename, dirname
from collections import Counter, defaultdict
from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp


def map_strain_to_seq(fna_file,fa_file,strain_to_cog_to_orfs):
    orfs = set([orf for strain,cog_to_orf in strain_to_cog_to_orfs.items() for cog,orfs in cog_to_orf.items() for orf in orfs])
    orf_to_seq_full = {header:seq for header,seq in sfp(open(fna_file)) if header.split(" ")[0] in orfs}
    orf_to_coordinate = {header.split(" ")[0]:header.split(" # ")[1:3] for header in orf_to_seq_full}
    orf_to_seq = {header.split(" ")[0]:seq for header,seq in orf_to_seq_full.items()}
    def merge_orfs(contig,start,end):
        for header,seq in sfp(open(fa_file)):
            if contig in header:
                return seq[start:end+1]
    def same_contig(contig,orfs):
        # case where all orfs are on the same contig
        nbs = [int(orf.split("_")[-1]) for orf in orfs]
        if list(range(min(nbs),max(nbs)+1))==sorted(nbs):
            # case where they are all next to each other
            coordinates = [int(coords) for orf in orfs for coords in orf_to_coordinate[orf]]
            start = min(coordinates)
            end = max(coordinates)
            seq = merge_orfs(contig,start,end)
            assert(seq!=""),"no sequence found for %"%contig
            return [seq]
        else:
            # in this case, we assume they are distincts sequences and should be accounted for.
            return [orf_to_seq[orf] for orf in orfs]

    #### deal with case where multiple orfs are present
    # we merge cogs if are on the same contig next to each other (see same_contig)
    # we append the seq if they are on different contigs or not following each other
    strain_to_cog_to_seq=defaultdict(lambda:defaultdict(list))
    for strain,cog_to_orfs in strain_to_cog_to_orfs.items():
        for cog,orfs in cog_to_orfs.items():
            if len(orfs)==1:
                orf = orfs[0]
                strain_to_cog_to_seq[strain][cog] = [orf_to_seq[orf]]
            else :
                contigs_to_orfs = defaultdict(list)
                for orf in orfs :
                    contig = "_".join(orf.split("_")[:-1])
                    contigs_to_orfs[contig].append(orf)
                # we need orfs always on the same order for futur concatenation
                sorted_contigs = sorted(contigs_to_orfs.keys())
                for contig in sorted_contigs:
                    orfs = sorted(contigs_to_orfs[contig])
                    if len(orfs)==1:
                        strain_to_cog_to_seq[strain][cog]+=[orf_to_seq[orf]]
                    else :
                        strain_to_cog_to_seq[strain][cog]+=same_contig(contig,orfs)
    return strain_to_cog_to_seq


def main(fa_file,fna_file,cog_file,output):
    strain_to_cog_to_orfs=defaultdict(lambda:defaultdict(list))
    with open(cog_file) as handle:
        _ = next(handle)
        for line in handle:
            orf_name,annotation = line.rstrip().split("\t")[:2]
            # from previous concatenation, orf_name is a list of cog_strain separated by |
            # this also end by the orf_number
            names = "_".join(orf_name.split("_")[:-1]).split("|")
            # each name is cog_strain, there can be more than 1 strain: imagine 2 strains sharing the same unitig for 1 or more cogs
            strains = list({name.split("_")[1] for name in names})
            # get cogs
            cogs = {name.split("_")[0] for name in names}
            cog = {annotation}&cogs
            if cog:
                cog = list(cog)[0]
                for strain in strains:
                    strain_to_cog_to_orfs[strain][cog].append(orf_name) 

    # get sequences corresponding
    strain_to_cog_to_seq = map_strain_to_seq(fna_file,fa_file,strain_to_cog_to_orfs)

    # output
    sorted_strain = sorted(strain_to_cog_to_orfs.keys())
    sorted_cogs = sorted({cog for cogs_to_orfs in strain_to_cog_to_orfs.values() for cog in cogs_to_orfs})
    with open(output,"w") as handle:
        for strain in sorted_strain:
            for cog in sorted_cogs:
                seqs = strain_to_cog_to_seq[strain][cog]
                # in some cases it is not possible to find the corresponding orf, we're ignoring these cases
                if seqs!=[]:
                    # it is possible that more than one cog is found for each strain
                    # if len(seqs)>1:
                    #     mag = basename(dirname(output))
                    #     print(mag,strain,cog,len(seqs))
                    for index,seq in enumerate(seqs):
                        name = "%s_%s_nb%s"%(cog,strain,index)
                        handle.write(">%s\n%s\n"%(name,seq))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fa_file", help="contig sequences")
    parser.add_argument("fna_file", help="prodigal nucleotide orfs")
    parser.add_argument("cog_file", help="cog annotation")        
    parser.add_argument("output")
    args = parser.parse_args()
    main(args.fa_file,args.fna_file,args.cog_file,args.output)


