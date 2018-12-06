#!/usr/bin/env python3.5
import argparse

def main(bed_file,Annotation_file,SCG_file) :
	Set_scg={line.rstrip() for line in open(SCG_file)}
	Dico_orfs_cogs={"_".join(line.split("\t")[0].split("_")[:-1]).split('.')[0]+"_"+line.split("\t")[0].split("_")[-1]:line.split("\t")[1] for line in open(Annotation_file) if line.split("\t")[1] in Set_scg}
	for line in open(bed_file) :
		Splitline=line.rstrip().split("\t")
		orf=Splitline[0]
		if orf in Dico_orfs_cogs :
			print( "\t".join(Splitline[1:]+[Dico_orfs_cogs[orf]]) )

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("bed", help=('bed file where the first column correspond to the unitig'))  
  parser.add_argument('Annotation',help = ("tsv file, first column correspond unitig, second correspond to COG id"))
  parser.add_argument('SCG',help = ("list of cogs we want to extract"))
  args = parser.parse_args()
  bed_file=args.bed
  Annotation_file=args.Annotation
  SCG_file=args.SCG
  main(bed_file,Annotation_file,SCG_file)































































