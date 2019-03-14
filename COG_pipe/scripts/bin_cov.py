#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

import sys
import re

def parse_len_cov(ctg_name):
    m = re.search("_length_(\d+)_cov_([\d\.]+)", ctg_name)
    assert m
    return int(m.group(1)), float(m.group(2))

if len(sys.argv) < 4:
    print("Usage: %s <binning.csv> <output> <K>" % sys.argv[0])
    sys.exit(1)

k = int(sys.argv[3])
first_line = True
total_len = dict()
total_cov = dict()

for l in open(sys.argv[1], "r"):
    if first_line:
        first_line = False
        continue

    ctg, b = l.strip().split(',')
    length, cov = parse_len_cov(ctg)

    length -= k
    if b in total_len:
        total_len[b] += length
        total_cov[b] += length * cov
    else:
        total_len[b] = length
        total_cov[b] = length * cov


with open(sys.argv[2], "w") as output:
    output.write("bin\tavg_coverage\n")
    for b in total_len:
        output.write("%s\t%f\n" % (b, total_cov[b] / total_len[b]))
