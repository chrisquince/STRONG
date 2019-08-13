#!/usr/bin/env python

import argparse
from collections import defaultdict,Counter 
import numpy as np
import sys
np.warnings.filterwarnings('ignore')

def Print_Final_annotation(Querry_Annotation) :
  # best evalue
  Querry_final_annotation=max(Querry_Annotation,key=lambda x:float(x[2]))
  # retranslate evalue again and format for printing
  Querry_final_annotation=Querry_final_annotation[0:2]+list(map(lambda x:"%.4g" %x,[10**(-Querry_final_annotation[2])]+Querry_final_annotation[3:]))
  # print
  print ("\t".join(Querry_final_annotation))

 
def main(Rpsblast_ouptut,Database_file,min_Evalue,min_Pid,min_Subject_Pid,min_coverage,min_Query_coverage) :
  Dico_Ref_annotation={line.rstrip().split("\t")[0]:line.rstrip().split("\t")[1] for line in open(Database_file)}
  Querry_Annotation=[]
  Current_query=''
  print( "\t".join(["Query","Subject","Evalue","PID","Subject_Pid","Coverage","Query_coverage"]))
  for line in open(Rpsblast_ouptut) :
    # Let's read the blast output with custom output format
    # qseqid sseqid evalue pident length slen qlen
    Splitline=line.rstrip().split("\t")
    Query=Splitline[0]
    # if querry is new and it not the first line and previous query was annotated, lets write previous querry annotation
    if (Querry_Annotation!=[])&(Current_query!="")&(Current_query!=Query) :
      Print_Final_annotation(Querry_Annotation)
      Querry_Annotation=[]
    Current_query=Query
    Full_Subject=Splitline[1]

    refid = Full_Subject.split("|")[2]
    if refid in Dico_Ref_annotation:
        Subject=Dico_Ref_annotation[refid]
    else:
        sys.stderr.write("Possible incompatibilities in COG RPSBlast database and " + str(Database_file) + "\n")
        continue
    Subject=Dico_Ref_annotation[Full_Subject.split("|")[2]]
    # sometimes for really low evalue, float conversion and/or comparison does not work, lets just take the -log10 of the evalue
    if 'e-' in Splitline[2] :
      Evalue=float(Splitline[2].split('e-')[1])-np.log10(float(Splitline[2].split('e-')[0]))
    else :
      Evalue=-np.log10(float(Splitline[2]))
    PID=float(Splitline[3])/100.
    Len_Alignemnt=float(Splitline[4])
    Subject_len=float(Splitline[5])
    Querry_len=float(Splitline[6])
    # usefull quantities 
    # compute coverage : size of Alignment over size of Subject
    Coverage=Len_Alignemnt/float(Subject_len)
    # compute number of Match over length of refference, Subject_Pid
    Subject_Pid=Coverage*PID
    # compute len of Alignment over size of the query
    Query_coverage=Len_Alignemnt/Querry_len
    # Filter things out
    if (Evalue>=min_Evalue)&(PID>=min_Pid)&(Subject_Pid>=min_Subject_Pid)&(Coverage>=min_coverage)&(Query_coverage>=min_Query_coverage) :
      Querry_Annotation.append([Query,Subject,Evalue,PID,Subject_Pid,Coverage,Query_coverage])

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("Rpsblast_ouptut", help=(' Output of rpsblast run, assumed to be in tabular format whith. '
         'columns: qseqid sseqid evalue pident score qstart qend sstart send length slen.'))  
  parser.add_argument('--cdd_cog_file',help = ('Supply a cdd to cog mapping file in a tsv format '
         'to take precedence over eutils fetching of name. '
         'Useful if running this script in parallel, since '
         'NCBI eutils has a limit on the number of requests per '
         'time unit you can make.'))
  parser.add_argument("-E",help="cutoff for Evalue",default=10**-10)
  parser.add_argument("-P",help="cutoff for Percentage IDentity (PID), between 0 and 1",default=0)
  parser.add_argument("-Q",help="Querry Coverage cutoff : Percentage of the Querry the Alignment does cover, between 0 and 1 ",default=0)
  parser.add_argument("-C",help="Subject Coverage cutoff : Percentage of the Subject the Alignment does cover, between 0 and 1 ",default=0.05)
  parser.add_argument("-R",help="Subject pid : Subject Coverage times the percentage identity, between 0 and 1",default=0)
  args = parser.parse_args()
  Rpsblast_ouptut=args.Rpsblast_ouptut
  Database_file=args.cdd_cog_file
  if args.E!=0 :
    min_Evalue=-np.log10(args.E)
  else :
    min_Evalue=args.E
  min_Pid=float(args.P)
  min_ref_pid=float(args.R)
  min_coverage=float(args.C)
  min_Query_coverage=float(args.Q)
  main(Rpsblast_ouptut,Database_file,min_Evalue,min_Pid,min_ref_pid,min_coverage,min_Query_coverage)




