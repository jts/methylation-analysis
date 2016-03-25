#! /usr/bin/env python

from collections import namedtuple
import sys
import csv
import os
import math
import argparse

# define datatypes
KmerModel = namedtuple('KmerModel', ['kmer', 'level_mean', 'level_stdv', 'sd_mean', 'sd_stdv'])
SummaryRecord = namedtuple('KmerModel', ['model', 'kmer', 'was_trained', 'num_training_events', 'trained_level_mean', 'trained_level_stdv'])
TableRow = namedtuple('TableRow', ['sample', 'treatment', 'lab', 'date',
                                   'model', 'alphabet', 'total_events',
                                   'total_kmers', 'trained_kmers',
                                   'd0_1', 'd0_5', 'd1_0', 'd2_0'])

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

def display_sample_name(s):
    fields = s.split("_")
    assert(fields[0] == "ecoli")
    return "\\emph{E. coli}" + " " + fields[1].upper()

def display_model(s):
    model_dict = {'t.006':'template', 'c.p1.006':'comp.pop1', 'c.p2.006':'comp.pop2'}
    return model_dict[s]

# From: http://stackoverflow.com/questions/3154460/python-human-readable-large-numbers
def display_number(n):
    suffix= ['','K','M','B' ]
    n = float(n)
    millidx = max(0,min(len(suffix)-1,
                    int(math.floor(0 if n == 0 else math.log10(abs(n))/3))))

    return '{:.0f}{}'.format(n / 10**(3 * millidx), suffix[millidx])

def display_treatment(s):
    return s.replace("_", "+").replace("MSssI", "M.SssI").replace("pcr", "PCR")

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

    for model in ["t.006", "c.p1.006", "c.p2.006"]:

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
                    row.d2_0]
            print " & ".join([str(x) for x in out]) + r'\\'
        print r'\hline'
    print r'\end{tabular}'
    print r'\end{adjustbox}'

    caption_str = "Model training results for %s-treated DNA over the %s alphabet.\n" % (display_treatment(treatment), alphabet)
    caption_str += "The final four fields are the number of k-mers where the mean of the trained Gaussian differs from the ONT-trained mean by more than x pA"
    caption_str += "\nTODO finalize."
    print r'\caption{' + caption_str + '}'
    print r'\end{table}'

#
# main
ont_models = dict()
load_ont_models_from_fofn("ont.alphabet_nucleotide.fofn", ont_models)
load_ont_models_from_fofn("ont.alphabet_cpg.fofn", ont_models)

parser = argparse.ArgumentParser( description='Generate latex-formatted tables after model training')
parser.add_argument('--treatment', type=str, required=True)
parser.add_argument('--alphabet', type=str, required=True)
args, files = parser.parse_known_args()

# We calculate the number of kmers that trained more than 0.1pA, 0.5pA, etc
# different than the reference model
diff_cuts = [ 0.1, 0.5, 1.0, 2.0]

table_rows = list()

for summary_file in files:
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
    short_alphabet = alphabet.split("_")[1]

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

        result = [sample, treatment, lab, date, model_short_name, short_alphabet, total_events, total_kmers, total_trained]
        result += diff_cut_count
        table_rows.append(TableRow(*result))

        #print "\t".join([str(x) for x in out])

print_table_latex(table_rows, args.treatment, args.alphabet)
