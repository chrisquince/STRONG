#!/usr/bin/env python
from collections import defaultdict
from linecache import getline
import numpy as np
import argparse
import marshal
import time 
import sys 

def Best_solution_yet(Coverage_profile,Indexed_coverage_profile,Annotation_Matrix,gene_to_annotation_index) :
    # save line number of each contig with marshal and using linecache, pretty sure this is not optimal yet, because it does not exploit the ordered part of the list of line.
    try :
        index_to_gene = marshal.load(open(Indexed_coverage_profile,"rb"))
        # TOFIX : sometimes encoding shenanigans : python 3.5, will need the .decode(),  seems like all string loaded are binary. So the dictionary throw key error since the genes from Set_genes are not binary string... So there is a need to decode binary to utf-8
        # TOFIX : sometimes encoding shenanigans : on python 3.7,  the .decode() make it fail. 
        try:
            index_to_gene={key.decode():values for key,values in index_to_gene.items()}
        except:
            index_to_gene={key:values for key,values in index_to_gene.items()}
        # unclear if summing each time is better than summing at the end. At least it's not stored in memory
        Sorted_index = sorted([index for index in index_to_gene])
        for line_nb in Sorted_index:
	        # +1 because there is the header
            # +1 is because  linecache start at 1 
            line = getline(Coverage_profile, line_nb+2)
            Annotation_Matrix[gene_to_annotation_index[index_to_gene[line_nb]]] += np.array(line.rstrip().split("\t")[1:]).astype(float) 
        return Annotation_Matrix
    except IOError :
        with open(Coverage_profile) as handle:
            _ = next(handle)
            index_to_gene={}
            for index,line in enumerate(handle) :
                gene=line.split("\t")[0]
                index_to_gene[index] = gene
                if gene in gene_to_annotation_index :
                    Annotation_Matrix[gene_to_annotation_index[gene]] += np.array(line.rstrip().split("\t")[1:]).astype(float)
        marshal.dump(index_to_gene, open(Indexed_coverage_profile, "wb" ) )
        return Annotation_Matrix

def main(Annotation_file,Coverage_profile,output):
    annotation_to_genes=defaultdict(list)
    for line in open(Annotation_file) :
        List_line=line.rstrip().split("\t")
        Annotation=List_line[0]
        if len(List_line)==1 :
            annotation_to_genes[Annotation]=[Annotation]
        else :
            Gene_list=List_line[1:]
            for gene in Gene_list :
                annotation_to_genes[Annotation].append(gene)

    Header=next(open(Coverage_profile))
    Sorted_sample=Header.rstrip().split("\t")[1:]
    Sorted_Annotation=sorted(list(annotation_to_genes.keys()))

    gene_to_annotation_index = {gene:Sorted_Annotation.index(annotation) for annotation,genes in annotation_to_genes.items() for gene in genes}

    Annotation_Matrix = np.zeros((len(Sorted_Annotation),len(Sorted_sample)))
    Indexed_coverage_profile=Coverage_profile+".pickle"

    # sum while reading:
    Annotation_Matrix=Best_solution_yet(Coverage_profile,Indexed_coverage_profile,Annotation_Matrix,gene_to_annotation_index)

    # better solutions (adress indexing/parrallisation) can be found at http://effbot.org/zone/wide-finder.htm
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

