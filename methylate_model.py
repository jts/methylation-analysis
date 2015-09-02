#!/usr/bin/env python

from __future__ import print_function
import string
import sys
import os
import operator
import itertools
import argparse

alphabet = [ 'A', 'C', 'G', 'M', 'T' ]

def make_all_mers(K):
    l = list()
    l.append("")

    for i in range(0, K):
        nl = list()
        for s in l:
            for b in alphabet:
                nl.append(s + b)
        l = nl
    return l

description = """
Modify an existing model to contain methyl-cytosines
"""
parser = argparse.ArgumentParser(description=description, epilog='')
parser.add_argument('input', action='store', help='model file')
args = parser.parse_args()

fn = args.input
K = 5

# Make a hash table containing all strings of length K over the alphabet
model = dict()
for kmer in make_all_mers(K):
    model[kmer] = tuple([0.0, 0.0, 0.0, 0.0, 0.0])

# Read the initial model from a file and set the model
out = list()
f = open(fn)
for line in f:
    line = line.rstrip()

    # copy then skip header lines
    if line[0] == '#' or line.find("kmer") == 0:
        print(line)
        continue

    fields = line.split()
    model[fields[0]] = tuple(fields[1:6])

intab = "M"
outtab = "C"
transtab = string.maketrans(intab, outtab)

# Set the values for the methylated k-mers depending upon the non-methylated versions
for kmer in model:
    no_m_kmer = kmer.translate(transtab)
    model[kmer] = model[no_m_kmer]

for key in sorted(model):
    print(key + "\t" + "\t".join([str(f) for f in model[key]]))
