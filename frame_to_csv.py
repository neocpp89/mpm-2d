#!/usr/bin/env python
import sys
import numpy
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
    print 'usage: FRAME_DATA plotvariables frame [outprefix]'
    print '       plotvariables are of the form \'xdata,ydata,...\'.'
    print '       frame can be either a single number \'#\', a range \'begin:end\','
    print '             or range with increment \'begin:step:end\'.'
    print '             end is not included in the count (e.g. 0:10 gives frames 0,1,...,9).'
    exit(127)

infile = sys.argv[1]
plotvars = sys.argv[2]

if (sys.argv[3] == '*'):
    framevars == None
else:
    toks = sys.argv[3].split(':')
    ftoks = map(int, toks)
    if (len(toks) == 1):
        framevars = [ftoks[0]]
    elif (len(toks) == 2):
        framevars = list(range(ftoks[0], ftoks[1]))
    elif (len(toks) == 3):
        framevars = list(range(ftoks[0], ftoks[2], ftoks[1]))

print framevars

print 'Using input file:', infile

pvars = plotvars.split(',')
try:
    pvars.index('active')
except:
    pvars.append('active')

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
            if (framevars == None or i in framevars):
                if (len(pvars) > 1):
                    print 'Writing', ', '.join(pvars[:-1]), 'and', pvars[-1], 'of frame', i, 'to csv file.'
                else:
                    print 'Writing', pvars[-1], 'of frame', i, 'to csv file.'

                if (len(sys.argv) > 4):
                    outfile = sys.argv[4] + '_' + str(i) + '.csv'
                else:
                    outfile = ".".join(infile.split('.')[:-1]) + '_' + str(i) + '.csv'

                print 'Creating output file:', outfile
                with open(outfile, 'w') as f_csv:
                    wcsv = csv.writer(f_csv)
                    row = ['Particle ID at ' + str(frame['time']) + 's']

                    for var in pvars:
                        row.append(var)
                    wcsv.writerow(row)
                    for j in xrange(0, len(frame[pvars[0]])):
                        row = [j]
                        for var in pvars:
                            row.append(frame[var][j])
                        wcsv.writerow(row)
            if i >= max(framevars):
                break
            i = i + 1

print "Bye."

