#!/usr/bin/env python
import argparse 
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import Counter,defaultdict
import os

def main(Bin_file,path,Fasta_file):
	# which SCG on which contig with which sequence ?
	Dico_contig_SCGSeq=defaultdict(lambda :defaultdict(list))
	for header,seq in SimpleFastaParser(open(Fasta_file)):
		contig="_".join(header.split(" ")[0].split('_')[:-1])
		SCG=header.split(" ")[1]
		Dico_contig_SCGSeq[contig][SCG].append([header,seq])
	# which SCG in which bins
	Dico_bins_SCG=defaultdict(lambda :defaultdict(list))
	for line in open(Bin_file) :
		contig,Bin=line.rstrip().split(',')
		if contig in Dico_contig_SCGSeq :
			for SCG,list_values in Dico_contig_SCGSeq[contig].items() :
				Dico_bins_SCG[Bin][SCG]+=list_values
	# which are 75% complete
	List_Mags=[Bin for Bin,List_contigs in Dico_bins_SCG.items() if sum(map(lambda x:x==1,Counter([SCG for SCG,list_fasta in List_contigs.items() for header,seq in list_fasta]).values()))>0.75*36]
	# create a folder by Mag with a folder by COG and their sequences inside
	for Mag in List_Mags :
		List_contigs_SCG=[]
		Mag_path=path+"Bin_"+Mag
		os.system("mkdir "+Mag_path)
		for COG,List_fasta in Dico_bins_SCG[Mag].items():
			COG_Path=Mag_path+"/"+COG+"/"
			os.system("mkdir "+COG_Path)
			Handle=open(COG_Path+"seq.fna","w")
			Handle.write("".join(map(lambda x:">"+x[0]+"\n"+x[1]+"\n",List_fasta)))
			Handle.close()

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("Bin_file",help="Binning result file, csv, first column is the contig, second column is the bin")
  parser.add_argument("SCG_Fasta",help="fasta file of Orfs annotated as SCG")
  parser.add_argument("folder",help="path to wher you want to store the bins folders",default=".")
  args = parser.parse_args()
  Bin_file=args.Bin_file
  path=args.folder
  Fasta_file=args.SCG_Fasta
  main(Bin_file,path,Fasta_file)
