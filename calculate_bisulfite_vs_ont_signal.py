#! /usr/bin/env python

import math
import sys
from collections import namedtuple

class SiteStats:
    def __init__(self):
        self.n = 0
        self.n_reads = 0
        self.n_methylated = 0

sites = dict()


for line in sys.stdin:
    (e_chr, e_start, e_end, _, _, e_strand, _, _, _, e_reads, e_m_percent, o_chr, o_start, o_end, o_kv) = line.split()
    key = o_chr + ":" + o_start + "-" + o_end
    kv_dict = dict( (k, v) for k, v in (pair.split("=") for pair in o_kv.split(";")) )
    ratio = kv_dict["LL_RATIO"]
    p_methylated = 1 / (1 + math.exp(-float(ratio)))
    
    if key not in sites:
        sites[key] = SiteStats()
        sites[key].ont_ratio = ratio
        sites[key].ont_p_methylated = p_methylated

    sites[key].n_reads += int(e_reads)
    sites[key].n_methylated += (int(e_reads) * float(e_m_percent) / 100)

# header
print "\t".join(["key", "bisulfite_depth", "bisulfite_percent_methylated", "ont_ratio", "ont_p_methylated"])

for key in sites:
    bisulfite_p = sites[key].n_methylated / sites[key].n_reads
    print "\t".join([key, str(sites[key].n_reads), str(bisulfite_p * 100), sites[key].ont_ratio, str(sites[key].ont_p_methylated)])

