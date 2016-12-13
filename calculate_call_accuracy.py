#! /usr/bin/env python
import math
import sys
import numpy
import random
from collections import namedtuple
from methylation_parsers import MethyltestRecord
from methylation_parsers import CpGIslandRecord
from methylation_parsers import BisulfiteRecord
import argparse

# tuples
CallRecord = namedtuple('CallRecord', ['loglik_ratio', 'is_true_methylated', 'context'])

# classes
class KmerStats:
    def __init__(self):
        self.called = 0
        self.correct = 0

# functions
def sample_sites(filename, is_methylated):
    all_crs = list()
    for line in open(filename):
        mt_record = MethyltestRecord(line.rstrip().split())

        # only use well-isolated singleton CpGs in this analysis
        if mt_record.num_cpgs != 1:
            continue

        assert(mt_record.sequence.count("CG") == 1)
        cg_position = mt_record.sequence.find("CG")
        context = mt_record.sequence[0:cg_position]
        all_crs.append( CallRecord(mt_record.loglik_ratio, is_methylated, context) )

    random.shuffle(all_crs)
    return all_crs[0:args.num_sites]

parser = argparse.ArgumentParser( description='Calculate call accuracy stats as a function of likelihood thresholds')
parser.add_argument('--unmethylated', type=str, required=True)
parser.add_argument('--methylated', type=str, required=True)
parser.add_argument('--pore', type=str, required=True)
parser.add_argument('--num-sites', type=str, required=False, default=100000)
args = parser.parse_args()

# Set the range of likelihood thresholds to use
likelihood_thresholds = numpy.arange(-20, 20.0, 0.25)

# Set a single threshold for the k-mer based analysis
kmer_analysis_threshold = 2.5
context_length = 3

# Read in the truth data
unmethylated_sites = sample_sites(args.unmethylated, False)
methylated_sites = sample_sites(args.methylated, True)
all_sites = methylated_sites + unmethylated_sites

#
# Accuracy contingency table
#

true_methylated = "TrueMethylated"
true_not_methylated = "TrueNotMethylated"
called_methylated = "CalledMethylated"
called_not_methylated = "CalledNotMethylated"

table_data = dict()
for s in all_sites:
    truth_tag = true_methylated if s.is_true_methylated else true_not_methylated
    call_tag = called_methylated if s.loglik_ratio > 0 else called_not_methylated
    if truth_tag not in table_data:
        table_data[truth_tag] = dict()
    if call_tag not in table_data[truth_tag]:
        table_data[truth_tag][call_tag] = 0
    table_data[truth_tag][call_tag] += 1

table_spec = """ \\begin{tabular}{ |c c|c|c| }
\hline
 & & \multicolumn{2}{|c|}{Called} \\\\ \cline{3-4}
 & & Methylated & Not Methylated  \\\\ \hline
 \multicolumn{1}{ |c  }{\multirow{2}{*}{Truth} } &
 \multicolumn{1}{ |c| }{Methylated} & %d & %d \\\\ \cline{2-4}
 \multicolumn{1}{ |c  }{}                        &
 \multicolumn{1}{ |c| }{Not Methylated} & %d & %d \\\\ \cline{1-4}
 \end{tabular} """

table_out = table_spec % (table_data[true_methylated][called_methylated],
                          table_data[true_methylated][called_not_methylated],
                          table_data[true_not_methylated][called_methylated],
                          table_data[true_not_methylated][called_not_methylated])

caption_str = "Contingency table for the positive and negative control accuracy assessment for %s data" % (args.pore)

table_writer = open("accuracy.table.%s.tex" % (args.pore), 'w')
table_writer.write(r'\begin{table}[h]' + "\n")
table_writer.write(r'\begin{adjustbox}{center}' + "\n")
table_writer.write(table_out + "\n")
table_writer.write(r'\end{adjustbox}' + "\n")
table_writer.write(r'\caption{' + caption_str + '}' + "\n")
table_writer.write(r'\end{table}' + "\n")

#
# Accuracy by kmer context
#

kmer_stats = dict()
for s in all_sites:

    # skip sites below the calling threshold
    if abs(s.loglik_ratio) < kmer_analysis_threshold:
        continue

    curr_context = s.context[len(s.context) - context_length:]
    if curr_context not in kmer_stats:
        kmer_stats[curr_context] = KmerStats()

    kmer_stats[curr_context].called += 1
    kmer_stats[curr_context].correct += ((s.loglik_ratio > 0) == s.is_true_methylated)

kmer_writer = open("accuracy.by_kmer.%s.tsv" % (args.pore), 'w')
kmer_writer.write("kmer\tcalled\tcorrect\taccuracy\n")
for kmer in kmer_stats:
    result = kmer_stats[kmer]
    kmer_writer.write("%s\t%d\t%d\t%.3f\n" % (kmer, result.called, result.correct, float(result.correct) / result.called))
kmer_writer.close()

#
# Accuracy by liklihood threshold
#

# open output files and write headers
pr_writer = open("accuracy.precision_recall.%s.tsv" % (args.pore), 'w')
pr_writer.write("threshold\ttrue_positive\tfalse_positive\ttrue_negative\tfalse_negative\tprecision\trecall\tspecificity\n")

acc_writer = open("accuracy.by_threshold.%s.tsv" % (args.pore), 'w')
acc_writer.write("threshold\tcalled\tcorrect\taccuracy\n")

for t in likelihood_thresholds:
    sys.stderr.write("t:" + str(t) + "\n")
    tp = 0
    fp = 0
    tn = 0
    fn = 0

    called = 0
    correct = 0

    for s in all_sites:

        # results set 1: precision/recall as a function of the likelihood ratio to call a site as methylated
        called_methylated = s.loglik_ratio > t
        tp += called_methylated and s.is_true_methylated
        fp += called_methylated and not s.is_true_methylated
        tn += not called_methylated and not s.is_true_methylated 
        fn += not called_methylated and s.is_true_methylated 
        
        # results set 2: accuracy as a function of the likelihood required for making a call
        called += abs(s.loglik_ratio) > t
        correct += abs(s.loglik_ratio) > t and ( (s.loglik_ratio > 0) == s.is_true_methylated)

    precision = float(tp) / (tp + fp)
    recall = float(tp) / (tp + fn)
    specificity = float(tn) / (tn + fp)
    accuracy = float(correct) / called
    pr_writer.write("%.2f\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\n" % (t, tp, fp, tn, fn, precision, recall, specificity))

    acc_writer.write("%.2f\t%d\t%d\t%.3f\n" % (t, called, correct, accuracy))

acc_writer.close()
pr_writer.close()
