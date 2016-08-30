#! /usr/bin/env python

from collections import namedtuple
import sys
import csv
import os
import math
import argparse
from methylation_utils import *
from methylation_parsers import ONTModel, TrainingSummary

# define datatypes
TableRow = namedtuple('TableRow', ['sample', 'treatment', 'lab', 'date',
                                   'model', 'alphabet', 'total_events',
                                   'total_kmers', 'trained_kmers',
                                   'd0_1', 'd0_5', 'd1_0', 'd2_0', 'd4_0'])

def load_ont_models_from_fofn(ont_fofn, out_model_set):
    f = open(ont_fofn)
    for filename in f:
        filename = filename.rstrip()
        out_model_set[filename] = ONTModel(filename)

def print_table_latex(table, treatment, alphabet):

    field_sep = " & "
    # Print a header
    header_fields = ["model", "sample", "run", "training events", "trained kmers"]
    header_fields += [str(x) for x in diff_cuts]

    table_spec_fields = "c" * len(header_fields)
    table_spec_str = "|" + "|".join(table_spec_fields) + "|"

    print r'\begin{table}[h]'
    print r'\begin{adjustbox}{center}'
    print r'\begin{tabular}{' + table_spec_str + "}"
    print r'\hline'
    print field_sep.join(header_fields) + r'\\'

    for model in model_names:

        print r'\hline'
        # Count the number of rows in this block
        num_rows = 0
        for row in table:
            if row.model != model or row.treatment != treatment or row.alphabet != alphabet:
                continue
            num_rows += 1

        # Print the rows
        first_row_in_model = True
        for row in table:
            if row.model != model or row.treatment != treatment or row.alphabet != alphabet:
                continue

            assert(diff_cuts[0] == 0.1)
            assert(diff_cuts[1] == 0.5)
            assert(diff_cuts[2] == 1.0)
            assert(diff_cuts[3] == 2.0)
            assert(diff_cuts[4] == 4.0)

            model_out = ""

            if first_row_in_model:
                model_out = r'\multirow{' + str(num_rows) + '}{*}{' + display_model(model) + r'}'
                first_row_in_model = False

            out = [ model_out,
                    display_sample_name(row.sample) + " (" + display_treatment(row.treatment) + ")",
                    row.lab + " " + row.date,
                    display_number(row.total_events),
                    row.trained_kmers,
                    row.d0_1,
                    row.d0_5,
                    row.d1_0,
                    row.d2_0,
                    row.d4_0]
            print " & ".join([str(x) for x in out]) + r'\\'
        print r'\hline'
    print r'\end{tabular}'
    print r'\end{adjustbox}'

    caption_str = "Model training results for %s-treated DNA over the %s alphabet.\n" % (display_treatment(treatment), alphabet)
    caption_str += "The final four fields are the number of k-mers where the mean of the trained Gaussian differs from the ONT-trained mean by more than x pA."
    print r'\caption{' + caption_str + '}'
    print r'\end{table}'

#
# main
ont_models = dict()
pore_version = "007"
model_names = [ mn + "." + pore_version for mn in ["t", "c.p1", "c.p2"]]

load_ont_models_from_fofn("ont.alphabet_nucleotide.fofn", ont_models)
load_ont_models_from_fofn("ont.alphabet_cpg.fofn", ont_models)

parser = argparse.ArgumentParser( description='Generate latex-formatted tables after model training')
parser.add_argument('--treatment', type=str, required=True)
parser.add_argument('--alphabet', type=str, required=True)
args, files = parser.parse_known_args()

# We calculate the number of kmers that trained more than 0.1pA, 0.5pA, etc
# different than the reference model
diff_cuts = [ 0.1, 0.5, 1.0, 2.0, 4.0]

table_rows = list()

for summary_file in files:
    summary = TrainingSummary(summary_file)

    for model_short_name in summary.models:
        ont_model_name = summary.get_ont_model_name(model_short_name)

        # Sanity check that the number of kmers is correct
        assert(ont_models[ont_model_name].get_num_kmers() == summary.get_num_kmers(model_short_name))

        total_kmers = ont_models[ont_model_name].get_num_kmers()
        total_events = 0
        total_trained = 0

        diff_cut_count = [0] * len(diff_cuts)

        for kmer in summary.models[model_short_name]:
            s = summary.models[model_short_name][kmer]
            m = ont_models[ont_model_name].kmers[kmer]
            diff = abs(s.trained_level_mean - m.level_mean)
            total_trained += s.was_trained
            total_events += s.num_training_events

            for (i,c) in enumerate(diff_cuts):
                if diff >= c:
                    diff_cut_count[i] += 1

        result = [summary.sample, summary.treatment, summary.lab, summary.date, model_short_name, summary.short_alphabet, total_events, total_kmers, total_trained]
        result += diff_cut_count
        table_rows.append(TableRow(*result))

        #print "\t".join([str(x) for x in out])

print_table_latex(table_rows, args.treatment, args.alphabet)
