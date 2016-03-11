#! /usr/bin/env python

import math
import sys
from collections import namedtuple
from methylation_parsers import MethyltestRecord
from methylation_parsers import CpGIslandRecord
from methylation_parsers import BisulfiteRecord
import argparse

parser = argparse.ArgumentParser( description='Calculate the percentage of methylated CpGs called by nanopore data as a function of distance from a gene')
parser.add_argument('-c', '--call-threshold', type=float, required=False, default=None)
parser.add_argument('-i', '--input', type=str, required=False)
parser.add_argument('-t', '--type', type=str)
args = parser.parse_args()

if args.input:
    in_file = open(args.input)
else:
    in_file = sys.stdin

distance_cuts = range(-5000, 5000, 50)
distance_sum_sites = [0] * len(distance_cuts)
distance_sum_methylated = [0] * len(distance_cuts)

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

    # the distance to the feature is the last field
    distance = int(fields[-1])
    
    for i in range(0, len(distance_cuts) - 1):
        if distance > distance_cuts[i] and distance < distance_cuts[i + 1]:
            distance_sum_sites[i] += num_covered_sites
            distance_sum_methylated[i] += num_methylated_sites

print "\t".join(["max_distance", "sites", "sites_methylated", "percent_methylated"])
for i in range(0, len(distance_cuts) - 1):
    p = float(distance_sum_methylated[i]) / distance_sum_sites[i] if distance_sum_sites[i] > 0 else 0
    print "%d\t%d\t%d\t%.3lf" % (distance_cuts[i], distance_sum_sites[i], distance_sum_methylated[i], p)
