import math

#
# Utility functions
#

# Split a string formated like key1=a;key2=b into a dictionary
def str2dict(s):
    return dict( (k, v) for k, v in (pair.split("=") for pair in s.split(";")) )

#
# This class parses a line from a *.methyltest.sites.bed file
# 
class MethyltestRecord:
    def __init__(self, fields):
        assert(len(fields) == 4)
        self.chromosome = fields[0]
        self.start = int(fields[1])
        self.end = int(fields[2])
        
        kv_dict = str2dict(fields[3])
        self.loglik_ratio = float(kv_dict['LogLikRatio'])
        self.num_cpgs = int(kv_dict['NumCpGs'])
        
    def is_region_callable(self, call_threshold):
        return abs(self.loglik_ratio) > (self.num_cpgs * call_threshold)

    def get_num_called_sites(self, call_threshold):
        if self.is_region_callable(call_threshold):
            return self.num_cpgs
        else:
            return 0

    def get_num_called_methylated(self, call_threshold):
        if self.is_region_callable(call_threshold) and self.loglik_ratio > 0:
            return self.num_cpgs
        else:
            return 0

    # the posterior probability that all sites in this island
    # are methylated 
    def posterior_methylated(self):
        return float(1 / (1 + math.exp(-self.loglik_ratio)))
    
    # the expected number of methylated sites calculated from the posterior
    def posterior_num_methylated_sites(self):
        return self.num_cpgs * self.posterior_methylated()

class CpGIslandRecord:
    def __init__(self, fields):
        assert(len(fields) == 4)
        self.chromosome = fields[0]
        self.start = int(fields[1])
        self.end = int(fields[2])
        kv_dict = str2dict(fields[3])
        self.gene = kv_dict["Gene"]

    def key(self):
        return self.chromosome + ":" + str(self.start) + "-" + str(self.end)

class BisulfiteRecord:
    def __init__(self, fields):
        assert(len(fields) == 11)
        self.chromosome = fields[0]
        self.start = int(fields[1])
        self.end = int(fields[2])
        self.num_reads = float(fields[9])
        self.percent_methylated = float(fields[10])

    def get_num_reads(self):
        return self.num_reads

    def get_num_methylated_reads(self):
        return self.num_reads * (self.percent_methylated / 100.0)
