#! /usr/bin/env python

from collections import namedtuple
import sys
import csv
import os
import math
import argparse
from methylation_utils import *
from methylation_parsers import *

parser = argparse.ArgumentParser( description='Calculate the difference between methyltrain means and ONT models')
parser.add_argument('--summary', type=str, required=True)
parser.add_argument('--ont-fofn', type=str, required=True)
args, files = parser.parse_known_args()

ont_model_set = dict()
load_ont_models_from_fofn("ont.alphabet_cpg.R7.fofn", ont_model_set)
load_ont_models_from_fofn("ont.alphabet_cpg.R9.fofn", ont_model_set)

summary = TrainingSummary(args.summary)

kmer_pattern = ['a', 'b', 'c', 'd', 'e', 'f']

print "\t".join(["model", "kmer", "m_position", "m_pattern",  "difference"])
for m in summary.models:
    ont_model = summary.get_ont_model_name(m)
    for k in summary.models[m]:
        sk = summary.models[m][k]
        mk = ont_model_set[ont_model].kmers[k]

        if sk.was_trained == 1 and k.count("M") == 1:
            m_pos = k.find("M")
            tmp = kmer_pattern[m_pos]
            kmer_pattern[m_pos] = 'M'
            print "\t".join([display_model(m), k, str(m_pos), ''.join(kmer_pattern), str(sk.trained_level_mean - mk.level_mean)])
            kmer_pattern[m_pos] = tmp
