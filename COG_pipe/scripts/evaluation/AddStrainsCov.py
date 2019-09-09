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

    args = parser.parse_args()

#    import ipdb; ipdb.set_trace()

    species_strains = defaultdict(set)
    strain_covs = defaultdict(lambda: np.zeros(96))

    with open(args.cov_file) as sin:

        for line in sin:

            line = line.rstrip()

            (sample, species, strain, cov, frac) = line.split('\t')

            strain_covs[strain][int(sample)] = float(cov)
            species_strains[species].add(strain)

    bFirst = True
    with open(args.assign_file) as sin:

        for line in sin:
            if bFirst:
                bFirst = False
            else:
                line = line.rstrip()

                (cluster, species, prob) = line.split(',')

                nStrains = len(list(species_strains[species]))

                strains = list(species_strains[species])
                total_covs = []

                for strain in strains:
                    total_covs.append(np.sum(strain_covs[strain][0:32]))
                nameString = ",".join(strains)
                covString = ",".join([str(x) for x in total_covs])
                print(re.sub("^D","Bin_",cluster) + "," + species + "," + prob + "," +
                      str(nStrains) + "," + nameString + "," + covString)


if __name__ == "__main__":
    main(sys.argv[1:])
