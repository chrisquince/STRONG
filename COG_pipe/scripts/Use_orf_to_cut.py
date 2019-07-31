#!/usr/bin/env python
# -*- coding: latin-1 -*-

import os 
from Bio.SeqIO.FastaIO import *
from Bio.Seq import Seq
import argparse
import numpy as np
from collections import Counter,defaultdict

def prodigal_gff_parser(Handle) :
	# specific to prodigal output, does not seem to be universal gff	
	# directly adapted from	SimpleFastaParser in Bio.SeqIO.FastaIO
	while True: 
		line = Handle.readline()
		if line == "": 
			return 
		if line[:16] == "# Sequence Data:" :
			break 
	while True: 
		if not line :
			break
		if line[:16] != "# Sequence Data:" : 
			print (line)
			raise ValueError("GFF Output from prodigal should start with '# Sequence Data:'") 
		seq_data = line.rstrip()
		Model_data = Handle.readline().rstrip() 
		line=Handle.readline().rstrip()
		ORF_list=[]
		if not line: 
			break 
		while line[0]!='#' :
			if not line :
				break
			ORF_list.append(line)
			line=Handle.readline().rstrip()
			if not line: 
				break 
		yield seq_data,Model_data,ORF_list 
	if not line: 
		return	# StopIteration



def get_gff_dico(Chunk_size,Gff_file):
	Limit_size=2*Chunk_size
	Dico_contigid_gff=defaultdict(list)
	Handle = open(Gff_file)
	Dico_contigs={}
	for Seq_data,Model,ORF_list in prodigal_gff_parser(Handle) :
		Seq_len=int(Seq_data.split(";")[1].split("=")[1])
		for ORF in ORF_list :
			(seqid,source,type,start,end,score,strand,phase,attributes)=ORF.split()
			Dico_contigs[seqid]=Seq_len
			start=int(start)-1
			end=int(end)-1
			if Seq_len>Limit_size :
				Dico_contigid_gff[seqid].append([start,end])
		if Seq_len>Limit_size :
			Dico_contigid_gff[seqid].append([end,Seq_len])
	return Dico_contigid_gff,Dico_contigs

def Cut_contigs(Chunk_size,Dico_contigid_gff):
	Contigid_Cut_location={}
	for key,List_start_end in Dico_contigid_gff.items() :
		List_cut_location=[0]
		List_end_orf=[end for start,end in List_start_end]
		Previous_cut=-1
		for end in List_end_orf :
			if (end - Previous_cut)>Chunk_size :
				if List_end_orf[-1]-end > Chunk_size : 
					List_cut_location.append(end+1)
					Previous_cut=end+1
		List_cut_location.append(List_end_orf[-1])
		Contigid_Cut_location[key]=[[List_cut_location[a],List_cut_location[a+1]] for a in range(len(List_cut_location)-1)]
	return Contigid_Cut_location

def Rename_ORF(ORF_list,Dico_contigid_Dico_ORfnb_index,Dico_Contigid_Cutlocation) :
	if ORF_list :
		Contig=ORF_list[0].split()[0]
		if Contig in Dico_contigid_Dico_ORfnb_index :
			Renamed_ORF_list=[]
			for index,ORF in enumerate(ORF_list) :
				#print ORF
				ORF_s=ORF.split("\t")
				Split_position=Dico_contigid_Dico_ORfnb_index[Contig][index]
				Starting_contig_position=Dico_Contigid_Cutlocation[Contig][Split_position][0]
				New_ORF_position=[int(string)-Starting_contig_position for string in ORF_s[3:5]]
				New_ORF_position[0]=(New_ORF_position[0]<0)*0+(New_ORF_position[0]>=0)*New_ORF_position[0] 
				New_ORF_name=Contig+"."+str(Split_position)
				New_ORF="\t".join([New_ORF_name]+ORF_s[1:3]+[str(number) for number in New_ORF_position]+ORF_s[5:])
				"\t".join([]+ORF_s[1:])
				Renamed_ORF_list+=[New_ORF]
			return Renamed_ORF_list
		else :
			return ORF_list
	else :
		return ORF_list

def delete_ending(sep,name) :
	return sep.join((name.split(sep)[:-1]))	

def Rename_GFF(Gff_file,Chunk_size,Dico_contigid_Dico_ORfnb_index,Dico_Contigid_Cutlocation) :
	Handle=open(delete_ending(".",Gff_file)+"_C"+str(Chunk_size//1000)+"K.gff","w")
	for (index,(Seq_data,Model,ORF_list)) in enumerate(prodigal_gff_parser(open(Gff_file))) :
		ORF_list=Rename_ORF(ORF_list,Dico_contigid_Dico_ORfnb_index,Dico_Contigid_Cutlocation)
		Chunk_gff="\n".join([Seq_data,Model]+ORF_list)+"\n"
		Handle.write(Chunk_gff)	
	Handle.close()

def Rename_FA(FA_file,Chunk_size,Dico_contigid_Dico_ORfnb_index) :
	# also take out ugly prodigal info on header, though they would still be availlable in NO_CUT directory
	Handle=open(delete_ending(".",FA_file)+"_C"+str(Chunk_size//1000)+"K."+FA_file.split(".")[-1],"w")
	Text_output=""
	for (index,(name,seq)) in enumerate(SimpleFastaParser(open(FA_file))) :
		ORF=name.split(" ")[0]
		Contig=delete_ending("_",ORF)
		if Contig in Dico_contigid_Dico_ORfnb_index :
			ORF_nb=int(ORF.split("_")[-1])
			contig_nb=Dico_contigid_Dico_ORfnb_index[Contig][ORF_nb]
			ORF=Contig+"."+str(contig_nb)+"_"+str(ORF_nb)
		Chunk_gff=">"+ORF+"\n"+seq+"\n"
		Text_output+=Chunk_gff
	Handle.write(Text_output)	
	Handle.close()

def main(Fasta_file,Gff_file,output_bed,Chunk_size,Replace) :
	Dico_contigid_gff,Dico_contigs=get_gff_dico(Chunk_size,Gff_file)
	Dico_Contigid_Cutlocation=Cut_contigs(Chunk_size,Dico_contigid_gff)	
	Dico_contigid_Dico_ORfnb_index={Contig:{index_orf:next(index for index,value in enumerate([i[1] for i in Dico_Contigid_Cutlocation[Contig]]) if ORF[1]<=value) for index_orf,ORF in enumerate(List_ORF)} for Contig,List_ORF in Dico_contigid_gff.items() if Contig in Dico_Contigid_Cutlocation}
	# Output cuts contigs		
	Lim=2*Chunk_size
	for title, sequence in SimpleFastaParser(open(Fasta_file)) :
		contig=title.split()[0]
		Dico_contigs[contig]=len(sequence)
		if (len(sequence)>Lim)&(contig in Dico_Contigid_Cutlocation) :
			for (index,(start,end)) in enumerate(Dico_Contigid_Cutlocation[contig]) :
				print (">"+contig+"."+str(index)+"\n"+sequence[start:end])
		else :
			print(">"+title+"\n"+sequence)
	# Output cuts contigs, as feature on intial contigs, in a bed file
	Handle=open(output_bed,"w")
	List_uncut_contigs=[[contig,"1",str(length),contig] for contig,length in Dico_contigs.items() if contig not in Dico_Contigid_Cutlocation]
	List_cut_contigs=[[contig,str(start+1),str(end),contig+"."+str(index)] for contig,list_coordinate in Dico_Contigid_Cutlocation.items() for index,(start,end) in enumerate(list_coordinate)]
	Handle.write("\n".join("\t".join(List) for List in List_uncut_contigs+List_cut_contigs))
	Handle.close()
	if Replace :
		# After cuting big contigs in chunks we need to update all prodigal annotations
		# Rename contigs in GFF file
		Rename_GFF(Gff_file,Chunk_size,Dico_contigid_Dico_ORfnb_index,Dico_Contigid_Cutlocation)
		# Rename contigs in FAA file
		FAA_file=delete_ending(".",Gff_file)+".faa"
		Rename_FA(FAA_file,Chunk_size,Dico_contigid_Dico_ORfnb_index)
		# Rename contigs in FNA file
		FNA_file=delete_ending(".",Gff_file)+".fna"
		Rename_FA(FNA_file,Chunk_size,Dico_contigid_Dico_ORfnb_index)
		# Move non cut files in a new folder
		CWD=os.getcwd()
		GFF_Directory="/".join(Gff_file.split("/")[:-1])
		if GFF_Directory :
			os.chdir(GFF_Directory)
		AGFFD=os.getcwd()
		os.chdir(CWD)
		os.system("mkdir -p " +AGFFD+"/No_Cut_Prodigal")
		os.system("mv " +AGFFD+"/"+delete_ending(".",Gff_file).split('/')[-1]+".* "+AGFFD+"/No_Cut_Prodigal/" )

	


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("Fasta_file", help="Fasta file from your assembly that you want to cut in bits")
	parser.add_argument("GFF", help="Gff file from, for instance prodigal output")
	parser.add_argument("output_bed", help="where to write the contig bed file")
	parser.add_argument("-c", help="chunk size",default="10000")
	parser.add_argument("-r", help="option to replace prodigal output with updated cut up entries",action='store_true',default=False)
	args = parser.parse_args()
	Chunk_size=int(args.c)
	Fasta_file=args.Fasta_file
	Gff_file=args.GFF
	Replace=args.r
	main(Fasta_file,Gff_file,output_bed,Chunk_size,Replace)



