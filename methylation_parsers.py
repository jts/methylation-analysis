import math
import os
import csv
import methylation_utils
from collections import namedtuple
from methylation_utils import *

#
# Parse a record from a methyltest.sites.bed file
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
        self.sequence = kv_dict['Sequence']
        self.readidx = kv_dict['ReadIdx']

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

    # the posterior probability that all sites in this island are methylated
    def posterior_methylated(self):
        return float(1 / (1 + math.exp(-self.loglik_ratio)))

    # the expected number of methylated sites calculated from the posterior
    def posterior_num_methylated_sites(self):
        return self.num_cpgs * self.posterior_methylated()

#
# Parse a record from a CPG island file
#
class CpGIslandRecord:
    def __init__(self, fields):
        assert(len(fields) == 4)
        self.chromosome = fields[0]
        self.start = int(fields[1])
        self.end = int(fields[2])
        kv_dict = str2dict(fields[3])
        self.gene = kv_dict["Feature"]

    def key(self):
        return self.chromosome + ":" + str(self.start) + "-" + str(self.end)

#
# Parse a BED record from a bisulfite file
#
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

#
# Parse an ONT model file
#
ModelKmer = namedtuple('ModelKmer', ['kmer', 'level_mean', 'level_stdv', 'sd_mean', 'sd_stdv'])
class ONTModel:
    def __init__(self, filename):
        fh = open(filename)
        self.kmers = dict()

        for line in fh:
            line = line.rstrip()

            # copy then skip header lines
            if line[0] == '#' or line.find("kmer") == 0:
                continue

            fields = line.split()
            a = ModelKmer(fields[0],
                          float(fields[1]),
                          float(fields[2]),
                          float(fields[3]),
                          float(fields[4]))
            self.kmers[a.kmer] = a

    def get_num_kmers(self):
        return len(self.kmers)

def load_ont_models_from_fofn(ont_fofn, out_model_set):
    f = open(ont_fofn)
    for filename in f:
        filename = filename.rstrip()
        out_model_set[filename] = ONTModel(filename)

#
# Parse a .summary file output by methyltrain
#
SummaryRecord = namedtuple('SummaryRecord', ['model', 'kmer', 'was_trained', 'num_training_events', 'trained_level_mean', 'trained_level_stdv'])
class TrainingSummary:
    def __init__(self, filename):
        
        # parse the structured file name
        fn_fields = os.path.basename(filename).rstrip().split(".")
        
        assert(fn_fields[0] == "methyltrain")
        assert(fn_fields[-1] == "summary")
        self.sample = fn_fields[1]
        self.treatment = fn_fields[2]

        # Switch between R7/R9 name parsing (8 vs 9) fields
        is_r9_file = filename.find(".r9.") != -1
        if is_r9_file:
            assert(len(fn_fields) == 9)
            self.pore = fn_fields[3]

        else:
            assert(len(fn_fields) == 8)
            # r7 files don't have the pore field
            self.pore = "r7"

        # R9 files have these fields offset by 1 to account for the pore string
        self.lab = fn_fields[3 + is_r9_file]
        self.date = fn_fields[4 + is_r9_file]
        self.alphabet = fn_fields[5 + is_r9_file]
        self.short_alphabet = self.alphabet.split("_")[1]

        # read kmers
        fh = open(filename)
        self.models = dict()
        reader = csv.DictReader(fh, delimiter="\t")
        
        for record in reader:
            s = SummaryRecord(record["model_short_name"],
                              record["kmer"],
                              int(record["was_trained"]),
                              int(record["num_events_for_training"]),
                              float(record["trained_level_mean"]),
                              float(record["trained_level_stdv"]))

            if s.model not in self.models:
                self.models[s.model] = dict()
            self.models[s.model][s.kmer] = s

    def get_num_kmers(self, model_name):
        return len(self.models[model_name])

    def get_ont_model_name(self, model_short_name):
        return model_short_name + ".ont." + self.alphabet + ".model"

