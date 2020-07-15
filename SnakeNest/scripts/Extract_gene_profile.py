#!/usr/bin/env python
from collections import defaultdict
from linecache import getline
import numpy as np
import argparse
import marshal
import time 
import sys 

version = "%s.%s"%sys.version_info[:2]

def Best_solution_yet(Coverage_profile,Indexed_coverage_profile,Set_genes) :
	# save line number of each contig with marshal and using linecache, pretty sure this is not optimal yet, because it does not exploit the ordered part of the list of line.
	Indexed_coverage_profile+="3"
	try :
		Dico_gene_index=marshal.load(open(Indexed_coverage_profile,"rb"))
		# TOFIX : sometimes encoding shenanigans : python 3.5, will need the .decode(),  seems like all string loaded are binary. So the dictionary throw key error since the genes from Set_genes are not binary string... So there is a need to decode binary to utf-8
		# TOFIX : sometimes encoding shenanigans : on python 3.7,  the .decode() make it fail. 
		if version == "3.5" :
			Dico_gene_index={key.decode():values for key,values in Dico_gene_index.items()}
		else :
			Dico_gene_index={key:values for key,values in Dico_gene_index.items()}
		# +2 is because linecache start at 1 (fucking assholes ) and because I got a line of headers.
		Sorted_index=sorted([Dico_gene_index[gene]+2 for gene in Set_genes])
		Nb_genes=len(Sorted_index)
		List_reduced=[getline(Coverage_profile, line_nb) for line_nb in Sorted_index]
		Dico_gene_profile={line.rstrip().split('\t')[0]:line for line in List_reduced}
	except IOError :
		Handle=open(Coverage_profile)
		Header=next(Handle)
		Dico_gene_index={}
		Dico_gene_profile={}
		for index,line in enumerate(Handle) :
			gene=line.split("\t")[0]
			Dico_gene_index[gene]=index
			if gene in Set_genes : 
				Dico_gene_profile[gene]=line
		marshal.dump(Dico_gene_index, open(Indexed_coverage_profile, "wb" ) )
	Reduced_Dico_gene_profile={key:np.array(Dico_gene_profile[key].rstrip().split('\t')[1:]).astype(float) for key in Set_genes}
	return Reduced_Dico_gene_profile

def main(Annotation_file,Coverage_profile,output):
	Dico_annotation_list_gene=defaultdict(list)
	for line in open(Annotation_file) :
		List_line=line.rstrip().split("\t")
		Annotation=List_line[0]
		if len(List_line)==1 :
			Dico_annotation_list_gene[Annotation]=[Annotation]
		else :
			Gene_list=List_line[1:]
			for gene in Gene_list :
				Dico_annotation_list_gene[Annotation].append(gene)
	Indexed_coverage_profile=Coverage_profile+".pickle"
	List_genes=[gene for list_gene in Dico_annotation_list_gene.values() for gene in list_gene]
	Set_genes=set(List_genes)
	Header=next(open(Coverage_profile))
	Sorted_sample=Header.rstrip().split("\t")[1:]
	Sorted_Annotation=sorted(list(Dico_annotation_list_gene.keys()))
	Reduced_Dico_gene_profile=Best_solution_yet(Coverage_profile,Indexed_coverage_profile,Set_genes)
	# better solutions (adress indexing/parrallisation) can be found at http://effbot.org/zone/wide-finder.htm
	Annotation_Matrix=np.array([sum([Reduced_Dico_gene_profile[gene] for gene in Dico_annotation_list_gene[Annotation]]) for Annotation in Sorted_Annotation])
	Handle=open(output,"w")
	Handle.write("Annotation\t"+"\t".join(Sorted_sample)+"\n")
	Handle.write("\n".join(["\t".join([Annotation]+[str(number) for number in Annotation_Matrix[index]]) for index,Annotation in enumerate(Sorted_Annotation) ])
	)
	Handle.close()


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("Annotation_file", help='Annotation file, should be a tsv with first field an annotation id, exemple KO/ARO/Species... and subsequent fields should be the list of all "gene/contig/orf" annotated by previous entry ')
	parser.add_argument("Coverage_profile", help="Coverage profile file of gene regarding samples, so that each line gives the coverage of a 'gene/orf/contig' in all samples, should be a tsv, first line should be a Header giving the name of all samples.")
	parser.add_argument("output", help="output file")
	args = parser.parse_args()
	Annotation_file=args.Annotation_file
	Coverage_profile=args.Coverage_profile
	output=args.output
	main(Annotation_file,Coverage_profile,output)

