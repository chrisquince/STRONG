#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys
import operator
from collections import defaultdict
from collections import Counter
import numpy as np

def main(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument("blast_file", help="tab delimited blast file")

    parser.add_argument("map_file", help="csv map file")

    args = parser.parse_args()
#    import ipdb;ipdb.set_trace()
    mapSeq = {}
    with open(args.map_file,'r') as source:
        for line in source:
            line = line.rstrip()
            toks = line.split(",")

            species = toks.pop(0)
            strain = toks.pop(0)
            idx = toks.pop(0)

            for tok in toks:
                mapSeq[tok] = strain 

    genes = set()
    haplo_gene_unc = defaultdict(lambda: defaultdict(list)) 
    #import ipdb; ipdb.set_trace()
    #COG0201_248730281,0.0,1.0,2.87583e-123
    haplo_matches = defaultdict(Counter)
    haplo_match_id = defaultdict(lambda: defaultdict(list))    
    haplo_match_length = defaultdict(Counter) 
    
    haplo_match_gene_ref = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
    haplo_length_gene_ref = defaultdict(lambda: defaultdict(dict))

    with open(args.blast_file,'r') as source:
        for line in source:
            line = line.rstrip()
            fields = line.split('\t')
                        
            haplo = fields[0].split('_')[-1]
            gene = "_".join(fields[0].split('_')[0:-1])
            ref = mapSeq[fields[1]]
            pid = float(fields[2])
            alignlength = int(fields[3])
            mismatches = int(fields[4]) + int(fields[5])
            matches = int(fields[3]) - int(fields[4]) -  int(fields[5])

            if float(matches) > haplo_match_gene_ref[haplo][gene][ref]:
                haplo_match_gene_ref[haplo][gene][ref] = float(matches)
                haplo_length_gene_ref[haplo][gene][ref] = float(alignlength)
                haplo_match_id[haplo][ref].append(pid)
    
    

    for haplo, gene_ref in haplo_match_gene_ref.items():
        for gene, ref_vals in haplo_match_gene_ref[haplo].items():
            for ref, matches in ref_vals.items():

                haplo_matches[haplo][ref] += matches
                haplo_match_length[haplo][ref] += haplo_length_gene_ref[haplo][gene][ref]

    haplo_matches_pid = defaultdict(Counter)
    for haplo, matches in haplo_matches.items():
        for ref, match in matches.items():
            haplo_matches_pid[haplo][ref] = match/haplo_match_length[haplo][ref]

    
    haplo_match_flat = {}
    for haplo, matches in haplo_matches.items():
        for ref, val in matches.items():
            haplo_match_flat[haplo + "-" + ref ] = val 

    hit_ref = {}
    haplo_hit = {}
    for (haplo_ref,val) in sorted(haplo_match_flat.items(), key=operator.itemgetter(1),reverse=True):
        haplo = haplo_ref.split("-")[0]
        ref = haplo_ref.split("-")[1]

        if haplo not in haplo_hit:
            hit_ref[ref]     = (haplo,val)
            haplo_hit[haplo] = (ref,val)

    

    for haplo, matches in haplo_matches.items():
        if haplo in haplo_hit:    
            bestMatch = haplo_hit[haplo][0]
            nHits = len(haplo_match_id[haplo][bestMatch])
            meanPid = np.mean(np.asarray(haplo_match_id[haplo][bestMatch]))
            pid = haplo_matches[haplo][bestMatch]/haplo_match_length[haplo][bestMatch]    
            diff = haplo_match_length[haplo][bestMatch] - haplo_matches[haplo][bestMatch]
        
            print(haplo + "\t" + bestMatch + "\t" + str(nHits) + "\t" + "{:10.4f}".format(pid) + "\t" + str(diff) +  "\t" + str(haplo_match_length[haplo][bestMatch]))    
        #print(haplo + "\t" + bestMatch + "\t" + str(haplo_matches[haplo][bestMatch]) + "\t" + str(meanPid) + "\t" + str(diff) + "\t" + str(haplo_match_length[haplo][bestMatch]) + "\t" + str(pid) + "\t" + str(mean_unc))

if __name__ == "__main__":
    main(sys.argv[1:])


