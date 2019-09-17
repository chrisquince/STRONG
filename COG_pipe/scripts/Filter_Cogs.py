#!/usr/bin/env python

import argparse
from collections import defaultdict, Counter
import numpy as np
import sys
np.warnings.filterwarnings('ignore')


def Print_Final_annotation(querry_annotation):
    # best evalue
    querry_final_annotation = max(querry_annotation, key=lambda x: float(x[2]))
    # retranslate evalue again and format for printing
    querry_final_annotation = querry_final_annotation[0:2]+list(map(
        lambda x: "%.4g" % x, [10**(-querry_final_annotation[2])]+querry_final_annotation[3:]))
    # print
    print("\t".join(querry_final_annotation))

def main(rpsblast_ouptut, database_file, min_evalue, min_pid, min_subject_pid, min_coverage, min_query_coverage):
    dict_ref_annotation = {line.rstrip().split("\t")[0]: line.rstrip().split("\t")[
        1] for line in open(database_file)}
    querry_annotation = []
    current_query = ''
    print("\t".join(["Query", "Subject", "Evalue", "PID",
                     "Subject_Pid", "Coverage", "Query_coverage"]))
    for line in open(rpsblast_ouptut):
        # Let's read the blast output with custom output format
        # qseqid sseqid evalue pident length slen qlen
        split_line = line.rstrip().split("\t")
        query = split_line[0]
        # if querry is new and it not the first line and previous query was annotated, lets write previous querry annotation
        if (querry_annotation != []) & (current_query != "") & (current_query != query):
            Print_Final_annotation(querry_annotation)
            querry_annotation = []
        current_query = query
        full_subject = split_line[1]

        refid = full_subject.split("|")[2]
        if refid in dict_ref_annotation:
            subject = dict_ref_annotation[refid]
        else:
            sys.stderr.write(
                "Possible incompatibilities in COG RPSBlast database and " + str(database_file) + "\n")
            continue
        subject = dict_ref_annotation[full_subject.split("|")[2]]
        # sometimes for really low evalue, float conversion and/or comparison does not work, lets just take the -log10 of the evalue
        if 'e-' in split_line[2]:
            evalue = float(split_line[2].split('e-')[1]) - \
                np.log10(float(split_line[2].split('e-')[0]))
        else:
            evalue = -np.log10(float(split_line[2]))
        pid = float(split_line[3])/100.
        len_alignemnt = float(split_line[4])
        subject_len = float(split_line[5])
        querry_len = float(split_line[6])
        # usefull quantities
        # compute coverage : size of Alignment over size of subject
        coverage = len_alignemnt/float(subject_len)
        # compute number of Match over length of refference, subject_Pid
        subject_Pid = coverage*pid
        # compute len of Alignment over size of the query
        query_coverage = len_alignemnt/querry_len
        # Filter things out
        if (evalue >= min_evalue) & (pid >= min_pid) & (subject_Pid >= min_subject_pid) & (coverage >= min_coverage) & (query_coverage >= min_query_coverage):
            querry_annotation.append(
                [query, subject, evalue, pid, subject_Pid, coverage, query_coverage])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("rpsblast_ouptut", help=(' Output of rpsblast run, assumed to be in tabular format whith. '
                                                 'columns: qseqid sseqid evalue pident score qstart qend sstart send length slen.'))
    parser.add_argument('--cdd_cog_file', help=('Supply a cdd to cog mapping file in a tsv format '
                                                'to take precedence over eutils fetching of name. '
                                                'Useful if running this script in parallel, since '
                                                'NCBI eutils has a limit on the number of requests per '
                                                'time unit you can make.'))
    parser.add_argument("-E", help="cutoff for evalue", default=10**-10)
    parser.add_argument(
        "-P", help="cutoff for Percentage IDentity (pid), between 0 and 1", default=0)
    parser.add_argument(
        "-Q", help="Querry coverage cutoff : Percentage of the Querry the Alignment does cover, between 0 and 1 ", default=0)
    parser.add_argument(
        "-C", help="subject coverage cutoff : Percentage of the subject the Alignment does cover, between 0 and 1 ", default=0.05)
    parser.add_argument(
        "-R", help="subject pid : subject coverage times the percentage identity, between 0 and 1", default=0)
    args = parser.parse_args()
    rpsblast_ouptut = args.rpsblast_ouptut
    database_file = args.cdd_cog_file
    if args.E != 0:
        min_evalue = -np.log10(args.E)
    else:
        min_evalue = args.E
    min_pid = float(args.P)
    min_ref_pid = float(args.R)
    min_coverage = float(args.C)
    min_query_coverage = float(args.Q)
    main(rpsblast_ouptut, database_file, min_evalue, min_pid,
         min_ref_pid, min_coverage, min_query_coverage)
