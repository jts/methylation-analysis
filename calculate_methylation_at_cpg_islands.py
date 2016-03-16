#! /usr/bin/env python

import math
import sys
from collections import namedtuple
from methylation_parsers import BisulfiteRecord
from methylation_parsers import MethyltestRecord
from methylation_parsers import CpGIslandRecord
import argparse

class IslandStats:
    def __init__(self):
        self.num_covered_sites = 0
        self.total_liklihood_ratio = 0
        self.sum_posterior = 0
        
        self.called_sites = 0
        self.called_sites_methylated = 0

parser = argparse.ArgumentParser( description='Calculate the percent methylated at Cpg Islands')
parser.add_argument('-c', '--call-threshold', type=float, required=False, default=None)
parser.add_argument('-i', '--input', type=str, required=False)
parser.add_argument('-t', '--type', type=str)
args = parser.parse_args()

island = dict()

if args.input:
    in_file = open(args.input)
else:
    in_file = sys.stdin

for line in in_file:
    fields = line.split()
    
    num_covered_sites = 0
    num_methylated_sites = 0

    if args.type == "ont":
        assert(args.call_threshold is not None)

        # parse ONT call data
        record = MethyltestRecord(fields[0:4])

        # is the evidence strong enough at this site to make a call?
        if record.is_region_callable(args.call_threshold):

            num_covered_sites = record.get_num_called_sites(args.call_threshold)
            num_methylated_sites = record.get_num_called_methylated(args.call_threshold)

    elif args.type == "bisulfite":
        record = BisulfiteRecord(fields[0:11])
               
        # the bisulfite calls contain 1 site per entry but N reads
        num_covered_sites = record.get_num_reads()
        num_methylated_sites = record.get_num_methylated_reads()
    else:
        sys.stderr("Unknown input type" + args.type)
        sys.exit(1)

    cpg_island = CpGIslandRecord(fields[-4:])
    key = cpg_island.key()

    #print e_chr, c_chr, e_m_percent, c_score, key
    if key not in island:
        island[key] = IslandStats()
        island[key].gene = cpg_island.gene

    island[key].called_sites += num_covered_sites
    island[key].called_sites_methylated += num_methylated_sites

# header
print "\t".join(["key", "called_sites", "called_sites_methylated", "feature"])

for key in island:
    print "\t".join([str(x) for x in [key, island[key].called_sites, island[key].called_sites_methylated, island[key].gene]])
