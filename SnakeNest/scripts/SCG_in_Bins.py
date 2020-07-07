#!/usr/bin/env python
import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import Counter, defaultdict
import os

def get_contig_split(contig_bed, contig, scg_bed, orf) : 
    start,end = scg_bed[orf]
    if end<start : 
        (end,start) = (start,end)
    best_overlap = []
    for split_contig,start_split,end_split in contig_bed[contig] :
        if end_split<start_split : 
            (end_split,start_split) = (start_split,end_split)
        # test if included in split 
        is_start_in = (start_split<=start<=end_split)
        is_end_in = (start_split<=end<=end_split)
        # simplest case : 
        if is_end_in&is_start_in :
            return split_contig 
        # shit overlap case
        elif is_end_in|is_start_in :
            best_overlap.append([split_contig,min(end_split,end)-max(start_split,start)])
    return max(best_overlap,key=lambda x:x[1])[0]

def get_mag_list(bins_to_scgs, scgs, threshold): 
    nb_scg = float(len(scgs))
    mags = []
    for bin_,scg_to_contigs in bins_to_scgs.items():
        scg_to_nb = Counter([scg for scg,fastas in scg_to_contigs.items() for header,seq in fastas])
        nb_unique_scg = sum([nb==1 for nb in scg_to_nb.values()])
        if nb_unique_scg>=threshold*nb_scg:
            mags.append(bin_)
    return mags


def main(Bin_file, Fasta_file, C10K_bed, orf_bed, path, Table, path_list, flag_all_bins,  scgs, threshold):
    # map contig + SCG to orf name and sequence
    Dico_contig_SCGSeq = defaultdict(lambda: defaultdict(list))
    for header, seq in SimpleFastaParser(open(Fasta_file)):
        contig = "_".join(header.split(" ")[0].split('_')[:-1])
        SCG = header.split(" ")[1]
        Dico_contig_SCGSeq[contig][SCG].append([header, seq])
    # ---- deal with case where we have split contigs  ------
    # get split contigs bed infos 
    contig_bed = defaultdict(list)
    for line in open(C10K_bed) :
        contig,start,end,split_contig = line.rstrip().split("\t")
        NB = len(contig.split("."))
        if (len(split_contig.split("."))==NB)|(contig not in Dico_contig_SCGSeq) :
            continue
        contig_bed[contig].append([split_contig,int(start),int(end)])
    # get SCG orf definitions 
    scg_bed = {line.rstrip().split("\t")[3]:list(map(int,line.rstrip().split("\t")[1:3])) for line in open(orf_bed) if line.rstrip().split("\t")[0] in Dico_contig_SCGSeq}
    # for each split contigs assign its SCG 
    Dico_contig_SCGSeq_splits = defaultdict(lambda: defaultdict(list))
    for contig, dict_scg in Dico_contig_SCGSeq.items() : 
        if contig not in contig_bed :
            continue
        for scg, list_header in dict_scg.items() :
            for header, seq in list_header :
                orf = header.split(" ")[0]
                split_contig = get_contig_split(contig_bed, contig, scg_bed, orf)
                Dico_contig_SCGSeq_splits[split_contig][scg].append([header, seq])
    Dico_contig_SCGSeq.update(Dico_contig_SCGSeq_splits)

    # map bins to SCG
    Dico_bins_SCG = defaultdict(lambda: defaultdict(list))
    Dico_bins_nbcontigs = defaultdict(int)
    for index, line in enumerate(open(Bin_file)):
        if index == 0:
            continue
        contig, Bin = line.rstrip().split(',')
        Dico_bins_nbcontigs[Bin] += 1
        if contig in Dico_contig_SCGSeq:
            for SCG, list_values in Dico_contig_SCGSeq[contig].items():
                Dico_bins_SCG[Bin][SCG] += list_values

# --------------- SCG output for concoct refine-----------------------------------
    if Table:
        List_SCG = sorted({key for dict in Dico_bins_SCG.values() for key in dict.keys()})
        Dico_bin_Array = {bin_nb: [len(Dico_bins_SCG[bin_nb][SCG]) if Dico_bins_SCG[bin_nb] else 0 for SCG in List_SCG] for bin_nb in Dico_bins_nbcontigs}
        List_bins = sorted(Dico_bins_nbcontigs.keys(), key=lambda x: int(x))
        SCG_table = [[Bin]+Dico_bin_Array[Bin] for Bin in List_bins]
        Handle = open(Table, "w")
        Handle.write(",".join(["Cluster"]+List_SCG)+"\n")
        Handle.write(
            "\n".join(map(lambda List: ','.join(map(str, List)), SCG_table)))
        Handle.close()

# -------------- create a folder by Mag with a folder by COG and their sequences inside -------------------
    if path:
        # fings mags which have at least threshold percent of their scg in a single copy
        List_Mags = get_mag_list(Dico_bins_SCG, scgs, threshold)
        # create a folder by Mag with a folder by COG and their sequences inside
        if flag_all_bins:
            bins=list(Dico_bins_SCG.keys())
        else :
            bins=List_Mags
        for Mag in bins:
            List_contigs_SCG = []
            Mag_path = path+"Bin_"+Mag
            os.system("mkdir -p "+Mag_path)
            Handle = open(Mag_path+"/SCG.fna", "w")
            Handle.write("".join(map(lambda x: ">"+x[0]+"\n"+x[1]+"\n", [
                         Fasta for List_fasta in Dico_bins_SCG[Mag].values() for Fasta in List_fasta])))
            Handle.close()

# -------------- create a list of all Mags ----------------------------------------------------------------
    if path_list:
        # fings mags which have at least threshold percent of their scg in a single copy
        List_Mags = get_mag_list(Dico_bins_SCG, scgs, threshold)
        # create a folder by Mag with a folder by COG and their sequences inside
        Handle = open(path_list, "w")
        Handle.write("\n".join(List_Mags))
        Handle.close()



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("Bin_file", help="Binning result file, csv, first column is the contig, second column is the bin")
    parser.add_argument("SCG_Fasta", help="fasta file of Orfs annotated as SCG")
    parser.add_argument("orf_bed", help="bed file of orfs definition : needed to handle concoct cut contigs")
    parser.add_argument("C10K_bed", help="bed file of cut contigs definition : needed to handle concoct cut contigs ")
    parser.add_argument("scg_file", help="file with one scg name by line")    
    parser.add_argument("-f", help="path to where you want to store the bins folders")
    parser.add_argument("-all", help="same as -f, but output for all bins, not just mags")
    parser.add_argument("-l", help="just output the list of MAG in this file")
    parser.add_argument("-t", help="name of a the SCG table")
    parser.add_argument("-T", help="percentage of unique copy mags",default="0.75")
    args = parser.parse_args()
    Bin_file = args.Bin_file
    Fasta_file = args.SCG_Fasta
    C10K_bed = args.C10K_bed
    orf_bed = args.orf_bed
    path = ""
    Table = ""
    path_list = ""
    flag_all_bins=0
    if args.f:
        path = args.f
    if args.all:
        path = args.all
        flag_all_bins=1
    if args.t:
        Table = args.t
    if args.l:
        path_list = args.l
    threshold = float(args.T)
    scgs = [line.rstrip() for line in open(args.scg_file)]
    main(Bin_file, Fasta_file, C10K_bed, orf_bed, path, Table, path_list, flag_all_bins, scgs, threshold)
