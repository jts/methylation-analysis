#! /usr/bin/env python

import math
import sys
from collections import namedtuple

class SiteStats:
    def __init__(self):
        self.n = 0
        self.n_reads = 0
        self.n_methylated = 0
        self.ont_ratio = 0
        self.ont_depth = 0
        self.ont_p_methylated = 0
        self.ont_template_ratio = 0
        self.ont_complement_ratio = 0

sites = dict()
used_records = dict()

for line in sys.stdin:
    (e_chr, e_start, e_end, _, _, e_strand, _, _, _, e_reads, e_m_percent, o_chr, o_start, o_end, o_kv) = line.split()
    key = o_chr + ":" + o_start + "-" + o_end
    kv_dict = dict( (k, v) for k, v in (pair.split("=") for pair in o_kv.split(";")) )
    
    if key not in sites:
        sites[key] = SiteStats()
        sites[key].context = kv_dict["SEQUENCE"]
        sites[key].num_sites = int(kv_dict["N_CPG"])
     
    # gather ont stats
    ont_key = kv_dict["READ_IDX"] + "." + key
    if ont_key not in used_records:
        used_records[ont_key] = 1
        
        sites[key].ont_depth += 1
        sites[key].ont_ratio += float(kv_dict["LL_RATIO"])
        p_methylated = 0
        
        # avoid overflow
        if sites[key].ont_ratio > -700:
            p_methylated = 1 / (1 + math.exp(-float(sites[key].ont_ratio)))
        sites[key].ont_p_methylated = p_methylated
    
        ratio_by_strand = kv_dict['LL_RATIO_BY_STRAND'].split(",")
        sites[key].ont_template_ratio += float(ratio_by_strand[0])
        sites[key].ont_complement_ratio += float(ratio_by_strand[1])
    
    bs_key = ".".join([e_chr, e_start, e_end, e_strand, key])
    if bs_key not in used_records:
        used_records[bs_key] = 1
        sites[key].n_reads += int(e_reads)
        sites[key].n_methylated += (int(e_reads) * float(e_m_percent) / 100)
        
# header
print "\t".join(["key", 
                 "bisulfite_depth", 
                 "bisulfite_percent_methylated", 
                 "ont_depth", 
                 "ont_ratio", 
                 "ont_template_ratio", 
                 "ont_complement_ratio", 
                 "context", 
                 "num_sites",
                 "ont_p_methylated"])

for key in sites:
    bisulfite_p = sites[key].n_methylated / sites[key].n_reads
    print "\t".join([key, 
                     str(sites[key].n_reads), 
                     str(bisulfite_p * 100), 
                     str(sites[key].ont_depth),
                     str(sites[key].ont_ratio),
                     str(sites[key].ont_template_ratio),
                     str(sites[key].ont_complement_ratio),
                     str(sites[key].context),
                     str(sites[key].num_sites),
                     str(sites[key].ont_p_methylated)])
