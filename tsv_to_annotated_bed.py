# /usr/bin/env python

import sys

header = sys.stdin.readline().rstrip().split()

for line in sys.stdin:
    fields = line.rstrip().split()

    kv_str = ";".join([ x + "=" + y for (x,y) in zip(header[3:], fields[3:]) ])
    print "\t".join(fields[0:3] + [kv_str])
