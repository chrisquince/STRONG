#!/usr/bin/env python
import argparse 
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import Counter

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("SCG_annotation",help="")
  parser.add_argument("Bin_Multiplier",help="Multiply the initial estimate of bins numbers by this")
  parser.add_argument("Bin_Maximum",help="The maximum number of bins allowed")
  args = parser.parse_args()
  binNumberInit = np.median(list(Counter([header.split(" ")[1] for header,seq in SimpleFastaParser(open(args.SCG_annotation))]).values()))
  binNumberFinal = np.min(binNumberInit*int(args.Bin_Multiplier), int(args.Bin_Maximum))
  # we take a look at all the SCG, take the median over the 36 of them, and multiply that by 2
  print(str(binNumberFinal))
