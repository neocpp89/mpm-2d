#!/usr/bin/env python
import sys
import numpy
import pylab
import math

import framereader

def frame_callback(frame, xvar, yvar, plothandle):
    # since we want to plot things, we don't care about most keys
    # delete them so we can store more frames in memory
    frame['xplot'] = [0]*len(frame['particles'])
    frame['yplot'] = [0]*len(frame['particles'])
    for i,p in enumerate(frame['particles']):
        frame['xplot'][i] = p[xvar]
        frame['yplot'][i] = p[yvar]
    frame['particles'] = None
    if (plothandle is not None):
        plothandle.set_xdata(frame['xplot'])
        plothandle.set_ydata(frame['yplot'])
        pylab.ylim([min(frame['yplot']),max(frame['yplot'])])
        pylab.xlim([min(frame['xplot']),max(frame['xplot'])])
        pylab.draw()
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




if (len(sys.argv) <= 2):
    print 'usage: FRAME_DATA plotvariables [frame]'
    print '       plotvariables are of the form \'xdata,ydata\'.'
    exit(127)

infile = sys.argv[1]
plotvars = sys.argv[2]
framevar = -1

if len(sys.argv) > 3:
    framevar = int(sys.argv[3])

print 'Using input file:', infile

pvars = plotvars.split(',')
xvar = pvars[0]
yvar = pvars[1]

pylab.ion()
phandle, = pylab.plot([], [], marker='o')

pylab.xlabel(xvar)
pylab.ylabel(yvar)

print 'Plotting', yvar, '(Y axis) vs', xvar, '(X axis).'

data = []
i = 0
with open(infile, 'r') as f_in:
    while True:
        frame = framereader.read_frame_cb(f_in,
                    lambda x: frame_callback(x, xvar, yvar, phandle),
                    particle_callback)
        if (frame is None):
            break
        else:
            print "Read Frame:", i
            data.append(frame)
            i = i + 1
            if ((framevar != -1 and i >= framevar)):
                break

print 'Got', len(data), 'frames of data. Enter for full movie.'
raw_input()

allxdata = []
allydata = []
for frame in data:
    allxdata.extend(frame['xplot'])
    allydata.extend(frame['yplot'])

xlims = [min(allxdata), max(allxdata)]
ylims = [min(allydata), max(allydata)]
i = 0
for frame in data:
    if framevar != -1 and i >= framevar:
        break
    pylab.title('Frame ' + str(i) +  ' Time ' + str(frame['time']))
    i = i + 1
    phandle.set_xdata(frame['xplot'])
    phandle.set_ydata(frame['yplot'])
    pylab.ylim(ylims)
    pylab.xlim(xlims)

    pylab.draw()

raw_input()
print "Bye."

