import sys, getopt
import os
import pandas as p
import numpy as np
import scipy.stats as ss
import scipy as sp
import scipy.misc as spm
import math
import argparse
import cPickle
import logging
import json
import glob
import re
import subprocess
import gzip
import atexit
import operator

from pprint import pprint
from Bio import SeqIO
from random import shuffle
from collections import defaultdict
from collections import Counter
from numpy.random import RandomState

def main(argv):

    parser = argparse.ArgumentParser()
    
    parser = argparse.ArgumentParser(description='Simulate a metagenomic data set')
    
    parser.add_argument("genomes_dir", help="directory with genomes")

    parser.add_argument("select_file", help="mapping strains to species")

    parser.add_argument("count_file", help="count file")

    parser.add_argument("--amb", help="use ambiguous reads",
                    action="store_true")

    args = parser.parse_args()

#    import ipdb; ipdb.set_trace()

    mapStrain = {}
    
    for fasta_file in glob.glob(args.genomes_dir + "/*seq.tmp"):
        seq_records = list(SeqIO.parse(fasta_file, "fasta"))

        m = re.search(r"(\d+_\d+)seq.tmp", fasta_file)

        strain_id = m.group(1)

        for record in seq_records:
            mapStrain[record.id] = strain_id
    
        print "Read fasta ids for:" + str(strain_id)


    mapSpecies = {}
    for line in open(args.select_file):
        line = line.rstrip()
            
        (species, strain, nSeqs, dir) = line.split("\t") 
        mapSpecies[strain] = species
    
    countContigStrainAmb  = defaultdict(Counter)
    countContigStrainUnamb  = defaultdict(Counter)
        

    countContigSpeciesAmb = defaultdict(Counter)
    countContigSpeciesUnamb = defaultdict(Counter)

    strain_file = open("Strain.csv", "w")
    species_file = open("Species.csv", "w")

    header = True
    mapIdx = {}
    for line in open(args.count_file):
        line = line.rstrip()
    
        if header:
            names = line.split("\t")
            header = False
            idx = 0
            names.pop(0)
            for name in names:
                m = re.search(r"(.*)_read_count_(.*)", name)
                
                amb = False
                if m.group(1) == 'amb':
                    amb = True
                mapIdx[idx] = (amb,m.group(2)) 
                idx = idx + 1 
        else:
            counts = line.split("\t")
            contig = counts.pop(0)
            idx = 0
            total = 0
            
            countSpecies = Counter()
            countStrain = Counter()
            for count in counts:
                (amb,seqId) = mapIdx[idx]
                strainId  = mapStrain[seqId]
                speciesId = mapSpecies[strainId]

                if count > 0:
                    if amb:
                        if args.amb:
                            countStrain[strainId] += int(count)
                            countSpecies[speciesId] += int(count)
                    else:
                        countStrain[strainId] += int(count)
                        countSpecies[speciesId] += int(count)

                    if amb:
                        countContigStrainAmb[contig][strainId]  += int(count)
                        countContigSpeciesAmb[contig][speciesId] += int(count)
                        if args.amb:
                            total += int(count)
                    else:
                        total += int(count)
                        countContigStrainUnamb[contig][strainId]  += int(count)
                        countContigSpeciesUnamb[contig][speciesId] += int(count)         
                
                idx = idx + 1
            
            sorted_strains = sorted(countStrain.items(), key=operator.itemgetter(1),reverse=True)
            sorted_species = sorted(countSpecies.items(), key=operator.itemgetter(1),reverse=True)            
            
            assignedStrains = []
            assignedString = ""
            if total > 0:
                for strain_tuple in sorted_strains:
                    (strain, count) = strain_tuple
        
                    if count > 0:
                        frac = float(count)/float(total)
                        assignedStrains.append((strain,count,frac))
                        assignedString += "," + ','.join(map(str, [strain,count,frac])) 
            nAssigned = len(assignedStrains)
            
            strain_file.write(contig + "," + str(total) + "," + str(nAssigned))            
            strain_file.write(assignedString)
            strain_file.write("\n")

            assignedSpecies = []
            assignedString  = ""
            if total > 0:
                for species_tuple in sorted_species:
                    (species, count) = species_tuple

                    if count > 0:
                        frac = float(count)/float(total)
                        assignedSpecies.append((species,count,frac))
                        assignedString += "," + ','.join(map(str, [species,count,frac]))
            nAssigned = len(assignedSpecies)

            species_file.write(contig + "," + str(total) + "," + str(nAssigned))
            species_file.write(assignedString)
            species_file.write("\n")

    strain_file.close()
    species_file.close()
    
if __name__ == "__main__":
    main(sys.argv[1:])
