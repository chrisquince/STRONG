#!/usr/bin/env python3.5
# -*- coding: latin-1 -*-

from Bio.SeqIO.FastaIO import SimpleFastaParser
from subprocess import Popen, PIPE
import argparse
import os

def main(fasta_file,Number_batch,Temp_dir) :
	extension=fasta_file.split(".")[-1]
	os.system("mkdir -p "+Temp_dir)
	process = Popen(['grep -c ">" ' +fasta_file], stdout=PIPE, stderr=PIPE,shell=True)
	nb_line = int(process.communicate()[0].split()[0])
	Nb_line_batch=nb_line/Number_batch
	Batch_num=0
	file_name=Temp_dir+"Batch_0."+extension
	Handle=open(file_name,"w")
	for (index,(header,seq)) in enumerate(SimpleFastaParser(open(fasta_file))) :
		if index > Nb_line_batch*(Batch_num+1) :
			if Batch_num==Number_batch-1 :
				Handle.write(">"+header+"\n"+seq+"\n")	
			else :
				Handle.close()
				Batch_num+=1
				file_name=Temp_dir+"Batch_"+str(Batch_num)+"."+extension
				Handle=open(file_name,"w")
		else :
			Handle.write(">"+header+"\n"+seq+"\n")
	Handle.close()
	

	
if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("fasta_file", help="fasta file you want to cut in batchs") 
	parser.add_argument("Number_batch", help="total number of batchs" )
	parser.add_argument("-t", help="temp directory where you want all your batchs",default='./temp_batchs/' )
	args = parser.parse_args()
	fasta_file=args.fasta_file
	Number_batch=float(args.Number_batch)
	Temp_dir=args.t
	main(fasta_file,Number_batch,Temp_dir)
