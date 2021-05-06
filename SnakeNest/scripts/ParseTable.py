#!/usr/bin/env python
from Bio import SeqIO
import sys
import argparse
import os

def main(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument("gff_file", help="prodigal gff file")

    parser.add_argument("fna_file", help="prodigal ORF sequences")

    args = parser.parse_args()

    #import ipdb; ipdb.set_trace() 

    contigs = set()

    handle = open(args.fna_file)
    for record in SeqIO.parse(handle, "fasta"):
        contig = '_'.join(record.id.split('_')[:-1])

        contigs.add(contig)
# Sequence Data: seqnum=7861;seqlen=78;seqhdr="NODE_7861_length_78_cov_10.000000"
# Model Data: version=Prodigal.v2.6.2;run_type=Metagenomic;model="34|Prochlorococcus_marinus_MIT_9313|B|50.7|11|0";gc_cont=50.70;transl_table=11;uses_sd=0    
    
    with open(args.gff_file) as fin:
        for line in fin:
            line = line.rstrip()
            if line.startswith('# Sequence Data:'):
                ttok= line.split(';')[-1]
                seqid=ttok[7:].strip('\"')

                if seqid in contigs:
                    line = next(fin, '')

                    line = line.rstrip()
                    if line.startswith('# Model Data:'):
                        tran_table = line.split(';')[4]
                        
                        assert (tran_table.startswith('transl_table='))

                        code=tran_table[len('transl_table='):]
                        
                        print (seqid + '\t' + code)

if __name__ == "__main__":
    main(sys.argv[1:])
