#!/usr/bin/env python
from __future__ import print_function

import sys, getopt
import os
import argparse
from collections import defaultdict

def main(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument("assign_file", help="fastq input")

    parser.add_argument("bin_file", help="text input file")

    parser.add_argument("bin", help="bin to select")


    #import ipdb; ipdb.set_trace()

    args = parser.parse_args()

    strains = {}
    with open(args.bin_file) as f:
        for line in f:
            line = line.rstrip()
            toks = line.split(',')
            #Bin_59,2096,0.99805502037024,3,1159201_0,710128_0,708616_0,48.4334169187,184.024767788,179.788788771

            if toks[0] == args.bin:
                nStrains = int(toks[3])

                for i in range(nStrains):
                    straini = toks[4 + i]
                    
                    strains[straini] = float(toks[4 + nStrains + i])


    try:
        strainBest = defaultdict(list)
        with open(args.assign_file) as g:
            #1  1123863_0   12      0.9964  7.8383e-03  37.0    10198.0 0.0035955056179775282
            for line in g:
                line = line.rstrip()
                toks = line.split('\t')
            
                haploj = toks[2]
                strainj = toks[1]
                err = 1.0 - float(toks[3])
                marg = float(toks[4])
                diver = float(toks[7])
                covG = float(toks[8])            
                strainBest[strainj].append((haploj,err,marg,diver,covG))
    
    except FileNotFoundError:
    #    print('This file doesn\'t exist')
        strainBest = defaultdict(list)
    for strain, cov in strains.items():

        hits = strainBest[strain]

        if len(hits) > 0:
            sorted_by_err = sorted(hits, key=lambda tup: tup[1])

            print(strain + ",Found," + str(cov) + "," + str(sorted_by_err[0][1]) + "," + str(sorted_by_err[0][2]) +"," + str(sorted_by_err[0][3]) + "," + str(sorted_by_err[0][4]))
            
            for h in range(1,len(hits)):
                print(strain + ",Repeated," + str(cov) + "," + str(sorted_by_err[h][1]) + "," + str(sorted_by_err[h][2]) +"," + str(sorted_by_err[h][3]) + "," + str(sorted_by_err[h][4]))


        else:
            print(strain + ",NotFound," + str(cov) + ",na,na,na")

if __name__ == "__main__":
    main(sys.argv[1:])
