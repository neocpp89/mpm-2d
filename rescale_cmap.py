#!/usr/bin/env python

import sys
import csv

if len(sys.argv) < 2:
    print sys.argv[0], "COLORMAP"
    print "Rescales COLORMAP for the viz program (new map written to stdout)."
    exit(0)

with open(sys.argv[1], 'r') as f:
    reader = csv.reader(f)
    rows = list(line for line in reader)
    maxval = float(rows[-1][1])
    minval = float(rows[0][1])
    rangeval = maxval - minval
    s = lambda x: (x - minval) / rangeval
    for line in rows:
        line[1] = str(s(float(line[1])))
        print ','.join(line)
