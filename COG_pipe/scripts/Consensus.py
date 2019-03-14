#!/usr/bin/env python3
from __future__ import print_function
import argparse 
from collections import Counter,defaultdict
import sys

def main(Bin_file):
	Handle=open(Bin_file)
	_=next(Handle)
	Dico_contig_bins={}
	for line in Handle:
		contig,bin=line.rstrip().split(",")
		Dico_contig_bins[contig]=bin
	number_of_dot=Counter(map(lambda x:len(x.split(".")),Dico_contig_bins.keys()))
	if len(number_of_dot)>3:
		print("Warning the number of dot in your contigs name is varying too much, it is not possible to distinguich cuts contigs from normal ones")
		exit()
	num_split=max(number_of_dot.keys())
	Dico_normal_contigs={contig:bins for (contig,bins) in  Dico_contig_bins.items() if len(contig.split("."))!=num_split}
	# get splits contigs and their bins 
	Dico_split=defaultdict(list)
	for contig,bins in Dico_contig_bins.items() :
		if len(contig.split("."))==num_split :
			contig=".".join(contig.split(".")[:-1])
			Dico_split[contig].append(bins)
	Dico_split={contig:Counter(list_bins) for contig,list_bins in Dico_split.items()}
	# Assign each contig to an unique bin  
	Dico_splitcontig_bins={contig:max(CounterBins.items(),key=lambda x:x[1])[0] for contig,CounterBins in Dico_split.items()}
	# Output results
	print("contig_id,0")
	for contig,bins in list(Dico_splitcontig_bins.items())+list(Dico_normal_contigs.items()) :
		print(contig+","+bins)
	# Add some warnings in case assignment is not so obvious :
	for  contig,CouterBins in Dico_split.items() :
		if len(CouterBins)>1 : 
			Sorted_Counter=sorted(CouterBins.items(),key=lambda x:x[1])
			if Sorted_Counter[-1][1]==Sorted_Counter[-2][1] :
				print('Warning contig %s is split equally between bins %s and %s' % (contig,Sorted_Counter[-1][0],Sorted_Counter[-2][0]),file=sys.stderr)
			elif Sorted_Counter[-1][1]<=0.5*sum([counter[1] for counter in Sorted_Counter]) :
				print('Warning contig %s is oversplit : less than 50%% of it is present in bin %s' % (contig,Sorted_Counter[-1][0]),file=sys.stderr)


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("Bin_file",help="Binning result file, csv, first column is the contig, second column is the bin")
  args = parser.parse_args()
  Bin_file=args.Bin_file
  main(Bin_file)
