#! /usr/bin/env python

from collections import namedtuple
import sys
import csv
import os

# define datatypes
KmerModel = namedtuple('KmerModel', ['kmer', 'level_mean', 'level_stdv', 'sd_mean', 'sd_stdv'])
SummaryRecord = namedtuple('KmerModel', ['model', 'kmer', 'was_trained', 'num_training_events', 'trained_level_mean', 'trained_level_stdv'])

def read_model(fn):
    fh = open(fn)
    model = dict()

    for line in fh:
        line = line.rstrip()

        # copy then skip header lines
        if line[0] == '#' or line.find("kmer") == 0:
            continue

        fields = line.split()
        a = KmerModel(fields[0],
                      float(fields[1]),
                      float(fields[2]),
                      float(fields[3]),
                      float(fields[4]))
        model[a.kmer] = a
    return model

def read_summary(fn):
    trained_model_set = dict()

    fh = open(fn)
    reader = csv.DictReader(fh, delimiter="\t")
    for record in reader:
        s = SummaryRecord(record["model_short_name"], 
                          record["kmer"], 
                          int(record["was_trained"]), 
                          int(record["num_events_for_training"]),
                          float(record["trained_level_mean"]),
                          float(record["trained_level_stdv"]))

        if s.model not in trained_model_set:
            trained_model_set[s.model] = dict()
        trained_model_set[s.model][s.kmer] = s
    return trained_model_set

def load_ont_models_from_fofn(ont_fofn, out_model_set):
    f = open(ont_fofn)
    for filename in f:
        filename = filename.rstrip()
        out_model_set[filename] = read_model(filename)

#
# main
#a
ont_models = dict()
load_ont_models_from_fofn("ont.alphabet_nucleotide.fofn", ont_models)
load_ont_models_from_fofn("ont.alphabet_cpg.fofn", ont_models)

# We calculate the number of kmers that trained more than 0.1pA, 0.5pA, etc
# different than the reference model
diff_cuts = [ 0.1, 0.5, 1.0, 2.0]

# Print a header
header_fields = ["sample", "treatment", "lab", "date", "model", "alphabet", "total_events", "total_kmers", "trained_kmers"]
header_fields += [str(x) for x in diff_cuts]
print "\t".join(header_fields)

for summary_file in sys.argv[1:]:
    summary = read_summary(summary_file)

    # Parse the summary filename
    fn_fields = os.path.basename(summary_file).rstrip().split(".")
    assert(len(fn_fields) == 8)
    assert(fn_fields[0] == "methyltrain")
    assert(fn_fields[-1] == "summary")

    sample = fn_fields[1]
    treatment = fn_fields[2]
    lab = fn_fields[3]
    date = fn_fields[4]
    alphabet = fn_fields[5]

    for model_short_name in summary:
        ont_model_name = model_short_name + ".ont." + alphabet + ".model"

        # Sanity check that the number of kmers is correct
        assert(len(ont_models[ont_model_name]) == len(summary[model_short_name]))

        total_kmers = len(ont_models[ont_model_name])
        total_events = 0
        total_trained = 0

        diff_cut_count = [0] * len(diff_cuts)

        for kmer in summary[model_short_name]:
            s = summary[model_short_name][kmer]
            m = ont_models[ont_model_name][kmer]
            diff = abs(s.trained_level_mean - m.level_mean)
            total_trained += s.was_trained
            total_events += s.num_training_events

            for (i,c) in enumerate(diff_cuts):
                if diff >= c:
                    diff_cut_count[i] += 1

        out = [sample, treatment, lab, date, model_short_name, alphabet, total_events, total_kmers, total_trained]
        out += diff_cut_count

        print "\t".join([str(x) for x in out])
