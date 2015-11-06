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

METHYLATION_SYMBOL = 'M'

alphabet = [ 'A', 'C', 'G', METHYLATION_SYMBOL, 'T' ]
def make_all_mers(K):
    return itertools.product(alphabet, repeat=K)

description = """
Modify an existing model to contain methyl-cytosines
"""
parser = argparse.ArgumentParser(description=description, epilog='')
parser.add_argument('input', action='store', help='model file')
parser.add_argument('--alphabet', action='store', help='the alphabet to expand to')
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

from_char = METHYLATION_SYMBOL
to_char = "N"

# Set up a translation table which will make unmethylated
# versions of the k-mers containing methylation sites.
if args.alphabet == "cpg" or args.alphabet == "dcm":
    to_char = "C"
elif args.alphabet == "dam":
    to_char = "A"
else:
    sys.stderr.write("unknown alphabet: " + args.alphabet)
    sys.exit(1)

# print a header line containing the alphabet
print("#alphabet\t" + args.alphabet)

translation_table = maketrans(from_char, to_char)

IDX_LEVEL_STDV = 1

# Set the values for the methylated k-mers depending upon the non-methylated versions
for kmer in model:
    no_m_kmer = kmer.translate(translation_table)
    values = list(model[no_m_kmer])
    
    # we offset level_stdv when the kmer contains a methylated site
    # to provide a bit of separation between the gaussians
    if kmer.find('M') >= 0:
        values[IDX_LEVEL_STDV] = float(values[IDX_LEVEL_STDV]) + 1
    model[kmer] = tuple(values)

for key in sorted(model):
    print(key + "\t" + "\t".join([str(f) for f in model[key]]))
