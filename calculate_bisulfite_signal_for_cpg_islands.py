#! /usr/bin/env python

import sys
from collections import namedtuple

class IslandStats:
    def __init__(self):
        self.n = 0
        self.n_methylated = 0

island = dict()

for line in sys.stdin:
    (e_chr, e_start, e_end, _, _, e_strand, _, _, _, e_reads, e_m_percent, i_chr, i_start, i_end, i_kv) = line.split()
    key = i_chr + ":" + i_start + "-" + i_end
    
    kv_dict = dict( (k, v) for k, v in (pair.split("=") for pair in i_kv.split(";")) )

    if key not in island:
        island[key] = IslandStats()
        island[key].gene = kv_dict['Gene']

    island[key].n += int(e_reads)
    island[key].n_methylated += (int(e_reads) * float(e_m_percent))
    
# header
print "\t".join(["key", "bisulfite_depth", "bisulfite_percent_methylated", "gene"])

# data
for key in island:
    print "\t".join([str(x) for x in [key, island[key].n, island[key].n_methylated / island[key].n, island[key].gene]])

