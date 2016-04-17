#! /usr/bin/env python

import math
import sys
import csv
import re
from methylation_parsers import MethyltestRecord

import argparse


parser = argparse.ArgumentParser( description='Convert from bed file to a tsv for methylation phasing and average analysis in R')
parser.add_argument('-c', '--call-threshold', type=float, required=False, default=2.5)
parser.add_argument('-i', '--input', type=str, required=False)

args=parser.parse_args()

if args.input:
    in_file = open(args.input)
else:
    in_file = sys.stdin



##Generate tsv file for output of processed bed
outfile=args.input
outfile=outfile.replace('sites.bed', 'phase.tsv')
tab=open(outfile, mode='w')
tabwriter=csv.writer(tab, delimiter='\t')
tabwriter.writerow(['chr', 'start', 'readix', 'loglikratio', 'numcpgs', 'meth', 'trinuc'])


for line in in_file:
    fields = line.split()
    assert(args.call_threshold is not None)

    ## parse ONT call data
    record = MethyltestRecord(fields[0:4])

    

    cgs=re.finditer("CG", record.sequence)
        
    firstcgloc=re.search("CG", record.sequence).start()

    for cgloc in cgs:
        coord=str(record.start+cgloc.start()-firstcgloc)
        thresh=args.call_threshold*record.num_cpgs
        if record.loglik_ratio > thresh:
            meth="1"
        elif record.loglik_ratio < -thresh:
            meth="-1"
        else:
            meth="0"
        
        trinuc=record.sequence[cgloc.start():(cgloc.start()+3)]
        tabwriter.writerow([record.chromosome, coord, record.readidx, record.loglik_ratio, record.num_cpgs, meth, trinuc])


tab.close()


                  





