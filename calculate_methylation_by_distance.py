#! /usr/bin/env python

import math
import sys
from collections import namedtuple
from methylation_parsers import MethyltestRecord
from methylation_parsers import CpGIslandRecord
from methylation_parsers import BisulfiteRecord
import argparse

def update_stats(chrom, distance, num_calls, num_methylated):

    # ignore sites outside of min/max distance
    if distance < distance_cuts[0] or distance > distance_cuts[-1]:
        return

    # calculate the distance bin this site falls in
    distance_bin = int( (distance - distance_cuts[0]) / bin_width )
    #print distance, distance_bin, distance_cuts[distance_bin], distance_cuts[distance_bin] + bin_width
    assert(distance >= distance_cuts[distance_bin] and distance < distance_cuts[distance_bin] + bin_width)

    if chrom not in total_sites:
        total_sites[chrom] = [0] * len(distance_cuts)
        methylated_sites[chrom] = [0] * len(distance_cuts)

    total_sites[chrom][distance_bin] += num_calls
    methylated_sites[chrom][distance_bin] += num_methylated

parser = argparse.ArgumentParser( description='Calculate the percentage of methylated CpGs called by nanopore data as a function of distance from a gene')
parser.add_argument('-c', '--call-threshold', type=float, required=False, default=None)
parser.add_argument('-i', '--input', type=str, required=False)
parser.add_argument('-t', '--type', type=str)
args = parser.parse_args()

if args.input:
    in_file = open(args.input)
else:
    in_file = sys.stdin

bin_width = 50
distance_cuts = range(-3000, 3000, bin_width)
autosome_names = dict([ ("chr" + str(i), 1) for i in range(1, 23)])

total_sites = dict()
methylated_sites = dict()

for line in in_file:
    fields = line.split()

    num_covered_sites = 0
    num_methylated_sites = 0

    if args.type == "ont":
        assert(args.call_threshold is not None)

        # parse ONT call data
        record = MethyltestRecord(fields[0:4])
        chrom = record.chromosome

        # is the evidence strong enough at this site to make a call?
        if record.is_region_callable(args.call_threshold):

            num_covered_sites = record.get_num_called_sites(args.call_threshold)
            num_methylated_sites = record.get_num_called_methylated(args.call_threshold)

    elif args.type == "bisulfite":
        record = BisulfiteRecord(fields[0:11])
        chrom = record.chromosome
               
        # the bisulfite calls contain 1 site per entry but N reads
        num_covered_sites = record.get_num_reads()
        num_methylated_sites = record.get_num_methylated_reads()
    else:
        sys.stderr("Unknown input type" + args.type)
        sys.exit(1)

    # ignore unplaced, epstein-barr virus sequence
    if chrom.find("chrUn") != -1 or chrom.find("random") != -1 or chrom == "chrEBV":
        continue

    # the strand of the feature is second last
    # this field is "." if a feature could not be found
    match_strand = fields[-2]
    if match_strand != "+" and match_strand != "-":
        continue

    # the distance to the feature is the last field
    distance = int(fields[-1])
    
    update_stats(chrom, distance, num_covered_sites, num_methylated_sites)
    update_stats("all", distance, num_covered_sites, num_methylated_sites)
    if chrom in autosome_names:
        update_stats("autosomes", distance, num_covered_sites, num_methylated_sites)

print "\t".join(["chromosome", "distance", "sites", "sites_methylated", "percent_methylated"])

for chrom in total_sites:
    for i in range(0, len(distance_cuts) - 1):
        f = float(methylated_sites[chrom][i]) / total_sites[chrom][i] if total_sites[chrom][i] > 0 else 0
        p = f * 100
        print "%s\t%d\t%d\t%d\t%.3lf" % (chrom, distance_cuts[i], total_sites[chrom][i], methylated_sites[chrom][i], p)

