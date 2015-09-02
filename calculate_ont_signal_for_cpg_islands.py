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
    (s_chr, s_start, s_end, s_kmer, s_score, i_chr, i_start, i_end, _, _, _, _, _, i_gene) = line.split()
    key = i_chr + ":" + i_start + "-" + i_end

    #print e_chr, c_chr, e_m_percent, c_score, key
    if key not in island:
        island[key] = IslandStats()
        island[key].gene = i_gene

    island[key].n += 1
    island[key].total += float(s_score)
for key in island:
    print "\t".join([str(x) for x in [key, island[key].n, island[key].total, island[key].gene]])
