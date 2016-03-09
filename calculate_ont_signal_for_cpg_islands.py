#! /usr/bin/env python

import math
import sys
from collections import namedtuple
from methylation_parsers import MethyltestRecord
from methylation_parsers import CpGIslandRecord
import argparse

parser = argparse.ArgumentParser( description='Calculate the percentage of methylated CpGs called by nanopore data')
parser.add_argument('-c', '--call-threshold', type=float, required=True)
args = parser.parse_args()

class IslandStats:
    def __init__(self):
        self.num_covered_sites = 0
        self.total_liklihood_ratio = 0
        self.sum_posterior = 0
        
        self.called_sites = 0
        self.called_sites_methylated = 0

island = dict()

for line in sys.stdin:
    fields = line.split()
    methyltest_record = MethyltestRecord(fields[0:4])
    cpg_island = CpGIslandRecord(fields[4:])

    key = cpg_island.key()

    #print e_chr, c_chr, e_m_percent, c_score, key
    if key not in island:
        island[key] = IslandStats()
        island[key].gene = cpg_island.gene

    island[key].sum_posterior += methyltest_record.posterior_num_methylated_sites()
    island[key].num_covered_sites += methyltest_record.num_cpgs

    # Stats on the sites with strong enough evidence to make a call
    if methyltest_record.is_region_callable(args.call_threshold):
        island[key].called_sites += methyltest_record.get_num_called_sites(args.call_threshold)
        island[key].called_sites_methylated += methyltest_record.get_num_called_methylated(args.call_threshold)

# header
print "\t".join(["key", "n_covered_sites", "posterior_methylated_sites", "called_sites", "called_sites_methylated", "gene"])

for key in island:
    print "\t".join([str(x) for x in [key, island[key].num_covered_sites, island[key].sum_posterior, island[key].called_sites, island[key].called_sites_methylated, island[key].gene]])
