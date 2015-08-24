#!/usr/bin/env python

import csv
import sys

if len(sys.argv) != 2:
    print sys.argv[0]+": FILE"
    sys.exit(0)

filename = sys.argv[1]

with open(filename) as f:
    r = csv.reader(f)
    header = next(r)
    midx = header.index('m')
    xidx = header.index('x')
    yidx = header.index('y')
    data = map(lambda row:map(float, row), r)
    total_mass = sum([drow[midx] for drow in data])
    average_mass = total_mass / len(data)
    tmx = 0
    tmy = 0
    tm = 0
    for drow in data:
        m = drow[midx]
        if m > average_mass:
            x = drow[xidx]
            y = drow[yidx]
            tm = tm + m
            tmx = tmx + m*x
            tmy = tmy + m*y
    comx = 0
    comy = 0
    if (tm > 0):
        comx = tmx / tm
        comy = tmy / tm
    print comx, comy
