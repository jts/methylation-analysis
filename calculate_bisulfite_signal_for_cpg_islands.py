#! /usr/bin/env python

import sys
from collections import namedtuple
from methylation_parsers import BisulfiteRecord
from methylation_parsers import CpGIslandRecord

class IslandStats:
    def __init__(self):
        self.n = 0
        self.n_methylated = 0

island = dict()

for line in sys.stdin:
    fields = line.split()
    bisulfite_record = BisulfiteRecord(fields[0:11])
    cpg_island = CpGIslandRecord(fields[11:])

    key = cpg_island.key()

    if key not in island:
        island[key] = IslandStats()
        island[key].gene = cpg_island.gene

    island[key].n += bisulfite_record.get_num_reads()
    island[key].n_methylated += bisulfite_record.get_num_methylated_reads()
    
# header
print "\t".join(["key", "bisulfite_depth", "bisulfite_methylated_reads", "bisulfite_percent_methylated", "gene"])

# data
for key in island:
    print "\t".join([str(x) for x in [key, island[key].n, island[key].n_methylated, (island[key].n_methylated / island[key].n) * 100, island[key].gene]])

