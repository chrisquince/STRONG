#!/usr/bin/env python
import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import Counter, defaultdict
import os
import errno
import sys

def main(argv):

    parser = argparse.ArgumentParser()
    
    parser.add_argument("bin_file", help="Binning result file, csv, first column is the contig, second column is the bin")
    
    parser.add_argument("contigs_fasta", help="fasta file of binning contigs")
    
    parser.add_argument("bin_path", help="path to where you want to store the bins folders")
    
    args = parser.parse_args()
#    import ipdb; ipdb.set_trace()
    
    contigs_seq = {}
    
    with open(args.contigs_fasta, 'r') as f:
        for header, seq in SimpleFastaParser(f):
            contig = header.split(" ")[0]
            
            contigs_seq[contig] = seq

    bin_assign = defaultdict(list)
    
    with open(args.bin_file, 'r') as fp:
        for cnt, line in enumerate(fp):
            if cnt > 0:
                line = line.rstrip()
                toks = line.split(",")
                bin_assign[toks[1]].append(toks[0])
    
    for bin_idx, contigs in bin_assign.items():
        
        binPath = args.bin_path + '/' + "Bin_" + bin_idx + '/'
        
        if not os.path.exists(binPath):
            try:
                os.makedirs(binPath, 0o700)
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise

        FastaFileName = binPath + '/contigs.fasta'
        
        with open(FastaFileName, 'w') as ff:
            for contig in contigs:
                ff.write('>' + contig + '\n')
                ff.write(contigs_seq[contig] + '\n')
            

if __name__ == "__main__":
    main(sys.argv[1:])

