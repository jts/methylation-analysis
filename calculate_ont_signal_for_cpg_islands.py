#! /usr/bin/env python

import math
import sys
from collections import namedtuple

call_threshold = 2.5

class IslandStats:
    def __init__(self):
        self.n = 0
        self.total = 0
        self.sum_posterior = 0
        self.estimated_methylated_sites = 0

island = dict()

for line in sys.stdin:
    (s_chr, s_start, s_end, s_kv, i_chr, i_start, i_end, i_kv) = line.split()
    
    key = i_chr + ":" + i_start + "-" + i_end
    s_kv_dict = dict( (k, v) for k, v in (pair.split("=") for pair in s_kv.split(";")) )
    i_kv_dict = dict( (k, v) for k, v in (pair.split("=") for pair in i_kv.split(";")) )

    #print e_chr, c_chr, e_m_percent, c_score, key
    if key not in island:
        island[key] = IslandStats()
        island[key].gene = i_kv_dict['Gene']

    ll = float(s_kv_dict['LL_RATIO'])
    
    if abs(ll) > call_threshold:
        island[key].n += 1
        island[key].total += ll
        island[key].sum_posterior += float(1 / (1 + math.exp(-ll)))
        island[key].estimated_methylated_sites += ll > 0

# header
print "\t".join(["key", "n_ont_sites", "sum_ll_ratio", "sum_posterior", "estimated_methylated_sites", "gene"])

for key in island:
    print "\t".join([str(x) for x in [key, island[key].n, island[key].total, island[key].sum_posterior, island[key].estimated_methylated_sites, island[key].gene]])
