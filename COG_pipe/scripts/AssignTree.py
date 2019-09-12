#!/usr/bin/env python2
import sys, getopt
import os
import numpy as np
import math
import argparse
import re
import operator

from ete3 import Tree
from collections import Counter
from collections import defaultdict
from operator import mul, div, eq, ne, add, ge, le, itemgetter

MIN_SUPPORT = 0.0

def read_lineage_file(lineage_file): 
    
    mapping = {}
    
    for line in open(lineage_file):
        line = line.rstrip()
    
        (taxaid, domainid, phylumid, classid, orderid, 
        familyid, genusid, speciesid) = line.split("\t")
    
        mapping[int(taxaid)]=[domainid, phylumid, classid, orderid, familyid, genusid, speciesid]

    return mapping

def read_mapping_file(mapping_file): 
    
    mapping = {}
    
    for line in open(mapping_file):
        line = line.rstrip()
    
        (taxaid,name) = line.split(",")
    
        mapping[name]= int(taxaid)
        

    return mapping


def main(argv):
    parser = argparse.ArgumentParser()
    
    parser.add_argument("tree_file", help="Newick tree file")

    parser.add_argument("map_file", help="Maps names to lineages")
    
    parser.add_argument("lineage_file", help="Lineages")

    args = parser.parse_args()

    #import ipdb; ipdb.set_trace()

    t = Tree(args.tree_file)

    mapping = read_mapping_file(args.map_file)
    
    lineages = read_lineage_file(args.lineage_file)

    leaves = t.get_leaves()
    clusters = []
    references = []
    for leaf in leaves:
        if re.search(r".*Bin.*", leaf.name):
            clusters.append(leaf)
        else:
            references.append(leaf)
    
    #get closest matching reference for each cluster
    
    clusterAssign = defaultdict(dict)
    
    for cluster in clusters:
        closest, closest_dist = cluster.get_farthest_node()
    #    print cluster.name
     #   sys.stdout.flush()
        for reference in references:
            dist = t.get_distance(cluster,reference, topology_only=False)
            #print dist
            if dist < closest_dist:
                closest = reference
                closest_dist = dist
        ancestor = t.get_common_ancestor(cluster,closest)
        assigns = Counter()
        matches = ancestor.get_leaves()
    
        collate_hits = list()
        for depth in range(7):
            collate_hits.append(Counter())
            
        for match in matches:
            if match.name in mapping:
                taxaid = mapping[match.name]
                if taxaid in lineages:
                    line_match = lineages[taxaid]
                
                    for depth in range(7):
                        collate_hits[depth][line_match[depth]]+=1.0
        
        for depth in range(7):
            collate = collate_hits[depth]
            dWeight = sum(collate.values())
        
            sortCollate = sorted(collate.items(), key=operator.itemgetter(1),reverse=True)

            nL = len(collate)
            if nL > 0:
                dP = 0.0
                if dWeight > 0.0:
                    dP = float(sortCollate[0][1])/dWeight
                    clusterAssign[cluster.name][depth] = (sortCollate[0][0],dP,dWeight)    
            else:    
                clusterAssign[cluster.name][depth] = ('None',1.0,1.0)

        

    for cluster in clusterAssign.keys():
        (assign,p,W) = clusterAssign[cluster][0]
        sys.stdout.write('%s'%cluster)
        if W >= MIN_SUPPORT:
    
            for depth in range(7):
                (assign,p,W) = clusterAssign[cluster][depth]
                sys.stdout.write('\t%d,%s,%f,%f'%(depth,assign,p,W))
            sys.stdout.write('\n')
            sys.stdout.flush()
        else:
            for depth in range(7):
                (assign,p,W) = clusterAssign[cluster][depth]
                sys.stdout.write('\t%d,NA,0.0,0.0'%(depth))
            sys.stdout.write('\n')
            sys.stdout.flush()
           
if __name__ == "__main__":
    main(sys.argv[1:])
