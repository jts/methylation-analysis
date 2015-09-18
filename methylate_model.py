#!/usr/bin/env python
from __future__ import print_function
import sys
import os
import operator
import itertools
import argparse

# see: http://stackoverflow.com/questions/6628306/attributeerror-module-object-has-no-attribute-maketrans
try:
    maketrans = ''.maketrans
except AttributeError:
    # fallback for Python 2
    from string import maketrans

alphabet = [ 'A', 'C', 'G', 'M', 'T' ]

def make_all_mers(K):
    return itertools.product(alphabet, repeat=K)

description = """
Modify an existing model to contain methyl-cytosines
"""
parser = argparse.ArgumentParser(description=description, epilog='')
parser.add_argument('input', action='store', help='model file')
args = parser.parse_args()

fn = args.input

# Make a hash table containing all strings of length K over the alphabet
model = dict()

# Read the initial model from a file and set the model
firstkmer = True
out = list()
f = open(fn)
for line in f:
    line = line.rstrip()

    # copy then skip header lines
    if line[0] == '#' or line.find("kmer") == 0:
        print(line)
        continue

    fields = line.split()
    if firstkmer:
        K = len(fields[0])
        firstkmer = False
    else:
        assert len(fields[0]) == K
    model[fields[0]] = tuple(fields[1:6])

# Add all missing kmers
for kmer in make_all_mers(K):
    kmer = "".join(kmer)
    if kmer not in model:
        model[kmer] = tuple([0.0, 0.0, 0.0, 0.0, 0.0])

intab = "M"
outtab = "C"
transtab = maketrans(intab, outtab)

# Set the values for the methylated k-mers depending upon the non-methylated versions
for kmer in model:
    no_m_kmer = kmer.translate(transtab)
    model[kmer] = model[no_m_kmer]

for key in sorted(model):
    print(key + "\t" + "\t".join([str(f) for f in model[key]]))
