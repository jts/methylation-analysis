#! /usr/bin/env python

import math
import sys
import csv
import re
from methylation_parsers import MethyltestRecord

import argparse

def update_meths(loc_tuple,meth, trinuc):

    if loc_tuple not in bigtab:
        ##coverage first, then num methylated
        bigtab[loc_tuple]=[1,meth, trinuc]
    else:
        bigtab[loc_tuple][0] += 1
        bigtab[loc_tuple][1] += meth


parser = argparse.ArgumentParser( description='Take in per CG TSV file and summarize methylation information at each position as a bismark style cov file')
parser.add_argument('-i', '--input', type=str, required=False)

args=parser.parse_args()

if args.input:
    in_file = open(args.input)
else:
    in_file = sys.stdin

reader=csv.reader(in_file, delimiter='\t')
next(reader, None) ##skip header


##Generate bedMethyl file
outfile=args.input
outfile=outfile.replace('.phase.tsv', '.cyto.txt')
tab=open(outfile, mode='w')
tabwriter=csv.writer(tab, delimiter='\t')

##Coverage table and methylated table, indexed by tuples of chr and coord
bigtab = dict()

for line in reader:

    
    location=(line[0], int(line[1]))
    
    meth=int(line[5])

    if (meth != 0):
        meth=(meth+1)/2
        update_meths(location, meth, line[6])


for key, value in bigtab.iteritems():
    ##Bismark cytosine report style file 
    tabwriter.writerow([key[0], key[1]+1, "*", value[1], value[0]-value[1], "CG", value[2]])

tab.close()


                  





