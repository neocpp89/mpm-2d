#!/usr/bin/env python

import sys
import csv

with open(sys.argv[1], 'r') as f:
    reader = csv.reader(f)
    for line in reader:
        line[1] = str(float(line[1])*3)
        print ','.join(line)
