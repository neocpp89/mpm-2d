#!/usr/bin/env python
import sys
import csv
import math
import numpy

import framereader

def frame_callback(frame, pvars):
    # only care about rho, delete other keys
    frame['mu_avg'] = sum(p['mu']for p in frame['particles'] if p['active'] != 0)/len(frame['particles'])
    frame['particles'] = None
    return frame

def particle_callback(particle):
    # calculate these for particles while reading
    pressure = -0.5 * (particle['sxx']+particle['syy'])
    s0_xx = particle['sxx'] + pressure
    s0_yy = particle['syy'] + pressure
    s0_xy = particle['sxy']
    tau = math.sqrt(0.5 * (s0_xx ** 2 + 2*s0_xy ** 2 + s0_yy ** 2))
    particle['pressure'] = pressure
    particle['tau'] = tau
    if (abs(pressure) > 0):
        particle['mu'] = particle['tau'] / pressure
    else:
        particle['mu'] = numpy.inf
    return particle

if (len(sys.argv) <= 1):
    print 'usage: FRAME_DATA [outfile]'
    print '       outfile consists of [frameID, frametime, [particle densities]] rows in csv format.'
    exit(127)

infile = sys.argv[1]
if (len(sys.argv) > 2):
    outfile = sys.argv[2]
else:
    outfile = ".".join(infile.split('.')[:-1]) + '_mu_avg.csv'

print 'Using input file:', infile
print 'Creating output file:', outfile

i = 0
with open(infile, 'r') as f_in:
    with open(outfile, 'w') as f_csv:
        wcsv = csv.writer(f_csv)
        row = ['Frame', 't', 'mu_avg']
        wcsv.writerow(row)
        while True:
            frame = framereader.read_frame_cb(f_in,
                        lambda x: frame_callback(x, 'm'),
                        lambda x: particle_callback(x))
            if (frame is None):
                break
            else:
                print "Read Frame:", i
                row = [i, frame['time'], frame['mu_avg']]
                wcsv.writerow(row)
                i = i + 1

print "Bye."

