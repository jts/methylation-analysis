#! /usr/bin/env python

import math
import sys
from collections import namedtuple

class IslandStats:
    def __init__(self):
        self.n = 0
        self.total = 0

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

    island[key].n += 1
    island[key].total += float(s_kv_dict['LL_RATIO'])

# header
print "\t".join(["key", "n_ont_sites", "sum_ll_ratio", "gene"])

for key in island:
    print "\t".join([str(x) for x in [key, island[key].n, island[key].total, island[key].gene]])
