#!/usr/bin/env python
import argparse 
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import Counter,defaultdict
import os

def main(Bin_file,Fasta_file,path,Table):
	# which SCG on which contig with which sequence ?
	Dico_contig_SCGSeq=defaultdict(lambda :defaultdict(list))
	for header,seq in SimpleFastaParser(open(Fasta_file)):
		contig="_".join(header.split(" ")[0].split('_')[:-1])
		SCG=header.split(" ")[1]
		Dico_contig_SCGSeq[contig][SCG].append([header,seq])
	# which SCG in which bins
	Dico_bins_SCG=defaultdict(lambda :defaultdict(list))
	Dico_bins_nbcontigs=defaultdict(int)
	for line in open(Bin_file) :
		contig,Bin=line.rstrip().split(',')
		if contig!="contig_id" :
			Dico_bins_nbcontigs[Bin]+=1
			if contig in Dico_contig_SCGSeq :
				for SCG,list_values in Dico_contig_SCGSeq[contig].items() :
					Dico_bins_SCG[Bin][SCG]+=list_values

#--------------- SCG output for concoct refine------------------------------------------------------------
	if Table : 
		List_SCG=sorted({key for dict in Dico_bins_SCG.values() for key in dict.keys()})
		Dico_bin_Array={bin_nb:[len(Dico_bins_SCG[bin_nb][SCG]) if Dico_bins_SCG[bin_nb] else 0 for SCG in List_SCG] for bin_nb in Dico_bins_nbcontigs}
		List_bins=sorted(Dico_bins_nbcontigs.keys(),key=lambda x:int(x))
		SCG_table=[[Bin]+Dico_bin_Array[Bin] for Bin in List_bins]
		Handle=open(Table,"w")
		Handle.write(",".join(["Cluster"]+List_SCG)+"\n")
		Handle.write("\n".join(map(lambda List:','.join(map(str,List)),SCG_table)))
		Handle.close()
#-------------- create a folder by Mag with a folder by COG and their sequences inside -------------------
	if path :
		# which are 75% complete
		List_Mags=[Bin for Bin,List_contigs in Dico_bins_SCG.items() if sum(map(lambda x:x==1,Counter([SCG for SCG,list_fasta in List_contigs.items() for header,seq in list_fasta]).values()))>0.75*36]
		# create a folder by Mag with a folder by COG and their sequences inside
		for Mag in List_Mags :
			List_contigs_SCG=[]
			Mag_path=path+"Bin_"+Mag
			os.system("mkdir -p "+Mag_path)
			Handle=open(Mag_path+"/SCG.fna","w")
			Handle.write("".join(map(lambda x:">"+x[0]+"\n"+x[1]+"\n",[Fasta for List_fasta in Dico_bins_SCG[Mag].values() for Fasta in List_fasta])))
			Handle.close()

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("Bin_file",help="Binning result file, csv, first column is the contig, second column is the bin")
  parser.add_argument("SCG_Fasta",help="fasta file of Orfs annotated as SCG")
  parser.add_argument("-f",help="path to where you want to store the bins folders")
  parser.add_argument("-t",help="name of a the SCG table")
  args = parser.parse_args()
  Bin_file=args.Bin_file
  Fasta_file=args.SCG_Fasta
  path=""
  Table=""
  if args.f :
  	path=args.f
  if args.t:
  	Table=args.t
  main(Bin_file,Fasta_file,path,Table)
