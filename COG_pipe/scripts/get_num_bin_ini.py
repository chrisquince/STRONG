#!/usr/bin/env python
import argparse 
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import Counter

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("SCG_annotation",help="")
  args = parser.parse_args()
  # we take a look at all the SCG, take the median over the 36 of them, and multiply that by 2
  print(str(int(2*np.median(list(Counter([header.split(" ")[1] for header,seq in SimpleFastaParser(open(args.SCG_annotation))]).values())))))
