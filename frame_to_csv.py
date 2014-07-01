#!/usr/bin/env python
import sys
import numpy
import pylab
import math
import csv

import framereader

def frame_callback(frame, pvars):
    # since we want to plot things, we don't care about most keys
    # delete them so we can store more frames in memory
    for var in pvars:
        frame[var] = [0]*len(frame['particles'])
    for i,p in enumerate(frame['particles']):
        for var in pvars:
            frame[var][i] = p[var]
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




if (len(sys.argv) <= 3):
    print 'usage: FRAME_DATA plotvariables frame [outfile]'
    print '       plotvariables are of the form \'xdata,ydata,...\'.'
    exit(127)

infile = sys.argv[1]
plotvars = sys.argv[2]
framevar = int(sys.argv[3])
if (len(sys.argv) > 4):
    outfile = sys.argv[4]
else:
    outfile = ".".join(infile.split('.')[:-1]) + str(framevar) + '.csv'

print 'Using input file:', infile
print 'Creating output file:', outfile

pvars = plotvars.split(',')

i = 0
with open(infile, 'r') as f_in:
    while True:
        frame = framereader.read_frame_cb(f_in,
                    lambda x: frame_callback(x, pvars),
                    particle_callback)
        if (frame is None):
            break
        else:
            print "Read Frame:", i
            if ((framevar != -1 and i >= framevar)):
                break
            i = i + 1

if (len(pvars) > 1):
    print 'Writing', ', '.join(pvars[:-1]), 'and', pvars[-1], 'of frame', i, 'to csv file.'
else:
    print 'Writing', pvars[-1], 'of frame', i, 'to csv file.'

with open(outfile, 'w') as f_csv:
    wcsv = csv.writer(f_csv)
    row = ['Particle ID at ' + str(frame['time']) + 's']
    for var in pvars:
        row.append(var)
    wcsv.writerow(row)
    for i in xrange(0, len(frame[pvars[0]])):
        row = [i]
        for var in pvars:
            row.append(frame[var][i])
        wcsv.writerow(row)

print "Bye."

