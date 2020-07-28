#!/usr/bin/env python
# -*- coding: latin-1 -*-
import os
import argparse
from os.path import basename
from collections import Counter, defaultdict
from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp



def generate_edges(paths,new_vertices):
	edges=[]
	for index,path in enumerate(paths[:-1]):
		next_path = paths[index+1]
		new_vertex = "%s_to_%s"%(path[0].split("_")[0],next_path[0].split("_")[0])
		new_vertices|={new_vertex}
		# edge from path to new_vertex
		edges.append("\t".join(["L",path[-1][:-1],path[-1][-1],new_vertex,"+","0M\n"]))
		# edge from path to new_vertex		
		edges.append("\t".join(["L",new_vertex,"+",next_path[0][:-1],next_path[0][-1],"0M\n"]))
	return edges

def main(paths,output,gfa_files):
	cog_to_gfa = {basename(gfa).split("_")[0]:gfa for gfa in gfa_files}
	cog_lines = [[header]+[s for s in seq.split(",")] for header,seq in sfp(open(paths))]
	strain_paths = defaultdict(list)
	# sort by strain
	cogs2 = set() # selected cogs for bayespath are not always outputed 
	for line in cog_lines:
		cog,strain = line[0].split("_")
		cogs2|={cog}
		strain_paths[strain].append(["%s_%s"%(cog,unitig) for unitig in line[1:]])
	# get new edges and buffer vertices
	new_vertices = set()
	new_strain_edges = "".join([new_edgs for strain,paths in strain_paths.items() for new_edgs in generate_edges(paths,new_vertices)])
	# add buffer vertices
	new_vertices = "".join(["S\t%s\t%s\tKC:i:10\tCL:z:#000000\tC2:z:#000000\t\n"%(name,100*'N') for name in new_vertices])
	# rename contigs and add edges,
	new_edges = ""
	for cog in cogs2:
		with open(cog_to_gfa[cog]) as handle:
			for line in handle:
				splitline = line.rstrip().split("\t")
				if line[0]=="S":
					new_contig = "%s_%s"%(cog,splitline[1])
					new_vertices += "S\t%s\t%s\n"%(new_contig,"\t".join(splitline[2:]))
				if line[0]=="L":
					new_contig1 = "%s_%s"%(cog,splitline[1])
					new_contig2 = "%s_%s"%(cog,splitline[3])
					new_edges += "\t".join(["L",new_contig1,splitline[2],new_contig2,splitline[4],splitline[5]])+"\n"
	new_edges += new_strain_edges
	# output joint gfa
	with open(output,"w") as handle:
		handle.write(new_vertices+new_edges)


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("haplo_paths", help="output from bayespath, path of haplotypes in each simplif .gfa subgraph")
	parser.add_argument("output")
	parser.add_argument('-l', nargs=argparse.REMAINDER,help="list of gfa files to joins")
	args = parser.parse_args()
	main(args.haplo_paths,args.output,args.l)
