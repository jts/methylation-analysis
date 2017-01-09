#! /usr/bin/env python

import sys

def parse_stderr_file(filename):
    cpu_time = ""
    cpu_time_key = "User time (seconds):"
    f = open(filename)
    for line in f:
        pos = line.find(cpu_time_key)
        if pos != -1:
            cpu_time = float(line[pos + len(cpu_time_key):])
    return cpu_time

print r'\begin{table}[h]'
print r'\begin{adjustbox}{center}'
print r'\begin{tabular}{|c|c|}'
print r'\hline'
print r'Dataset & CPU time (h) \\'
print r'\hline'

for filename in sys.argv[1:]:
    input_set = filename.replace(".sorted.bam.methlytest.stderr", "")
    print "%s & %.1f" % (input_set.replace("_", "\_"), parse_stderr_file(filename) / 3600) + r'\\'
print r'\hline'
print r'\end{tabular}'
print r'\end{adjustbox}'
print r'\end{table}'
