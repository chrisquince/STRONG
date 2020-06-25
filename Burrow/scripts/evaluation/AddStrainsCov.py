#!/usr/bin/env python
from __future__ import print_function
import sys
import re
import argparse
import numpy as np
from collections import defaultdict

# 0  663 1219076_0   0.212580070954  0.000583505311272

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("select_file", help="contig bin assignments")
    parser.add_argument("cov_file", help="long read maps")
    parser.add_argument("assign_file", help="read lengths")
    parser.add_argument("merge_file", help="read lengths")
    parser.add_argument('--nomerge', dest='merge', action='store_false')
    args = parser.parse_args()
#    import ipdb; ipdb.set_trace()
    species_strains = defaultdict(set)
    maxSampleNo = 96
    strain_covs = defaultdict(lambda: np.zeros(maxSampleNo))
    with open(args.cov_file) as sin:
        for line in sin:
            line = line.rstrip()
            (sample, species, strain, cov, frac) = line.split('\t')
            strain_covs[strain][int(sample)] = float(cov)
            species_strains[species].add(strain)
    

    mergeBins = defaultdict(list)
    mergeMap = {}
    if args.merge:
        with open(args.merge_file) as sin:
            for line in sin:
                line = line.rstrip()
                toks = line.split('\t')

                for mbin in toks[1:]:
                    mergeMap[mbin] = toks[0]
                    mergeBins[toks[0]].append(mbin)
    mergeSpecies = defaultdict(list)
    mergeStrains = defaultdict(list)
    mergeCovs = defaultdict(list)
    mergeProbs = defaultdict(list)
    b_first = True
    with open(args.assign_file) as sin:
        for line in sin:
            if b_first:
                b_first = False
            else:
                line = line.rstrip()
                (cluster, species, prob) = line.split(',')
                n_strains = len(list(species_strains[species]))
                strains = list(species_strains[species])
                total_covs = []
                for strain in strains:
                    total_covs.append(np.sum(strain_covs[strain][0:maxSampleNo]))
                
                binName = re.sub("^D", "Bin_", cluster)
                if binName not in mergeMap:
                    nameString = ",".join(strains)
                    covString = ",".join([str(x) for x in total_covs])
                    print(binName + "," + species + "," + prob + "," +
                        str(n_strains) + "," + nameString + "," + covString)
                else:
                    mbin = mergeMap[binName]
                    mergeStrains[mbin].extend(strains)
                    mergeCovs[mbin].extend(total_covs)
                    mergeSpecies[mbin].append(species)
                    mergeProbs[mbin].append(prob)

    for mbin,bins in mergeBins.items():
        nameString = ",".join(mergeStrains[mbin])
        covString = ",".join([str(x) for x in mergeCovs[mbin]])
        n_strains = len(mergeStrains[mbin])
        speciesM = "-".join(mergeSpecies[mbin])
        probs = "-".join(mergeProbs[mbin])
        print(mbin + "," + speciesM + "," + probs + "," +
                        str(n_strains) + "," + nameString + "," + covString)



if __name__ == "__main__":
    main(sys.argv[1:])
