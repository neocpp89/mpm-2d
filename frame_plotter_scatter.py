#!/usr/bin/env python
import sys
import numpy
import pylab
import time

import matplotlib.pyplot as plt
import matplotlib.animation as animation

def update_scatter(i, data, xvar, yvar, zvar, scat):
    xdata = [p[xvar] for p in data[i]['particles']]
    ydata = [p[yvar] for p in data[i]['particles']]
    cdata = [p[zvar] for p in data[i]['particles']]
    # scat.set_offsets(zip(xdata, ydata))
    scat.set_array(numpy.array(cdata))
    return scat,

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

# if state is 0, next line is a frame descriptor
# if state is 1, next line is a particle descriptor
FRAMEDESC = 0
PARTICLEDESC = 1
state = FRAMEDESC
particles_left = 0;
framenum = 0

data = []

with open(infile, 'r') as f_in:
    for line in f_in:
        if (state == FRAMEDESC):
            tok = line.split(' ')
            particles_left = int(tok[2])
            state = PARTICLEDESC
            data.append({'time':float(tok[1]), 'particles':[]})
        elif (state == PARTICLEDESC):
            particles_left = particles_left - 1
            tok = line.split(' ')
            tokf = map(float, tok);
            particle = {'m':tokf[0], 'v':tokf[1], 'x':tokf[2], 'y':tokf[3],
                        'x_t':tokf[4], 'y_t':tokf[5], 'sxx':tokf[6],
                        'sxy':tokf[7], 'syy':tokf[8], 'ux':tokf[9],
                        'uy':tokf[10]}
            data[framenum]['particles'].append(particle)
            if (particles_left <= 0):
                state = FRAMEDESC
                framenum = framenum + 1

print 'Got', len(data), 'frames of data.'

pylab.ion()

pvars = plotvars.split(',')
xvar = pvars[0]
yvar = pvars[1]
if (len(pvars) >= 3):
    zvar = pvars[2]

print 'Plotting', yvar, '(Y axis) vs', xvar, '(X axis).'

if framevar == -1:
#    frame = data[0]
#    fig = plt.figure()
#    pylab.title(zvar + ' at Time ' + str(frame['time']))
#    pylab.xlabel(xvar)
#    pylab.ylabel(yvar)
#    xdata = [p[xvar] for p in frame['particles']]
#    ydata = [p[yvar] for p in frame['particles']]
#    cdata = [p[zvar] for p in frame['particles']]
#    scat = plt.scatter(xdata, ydata, c=cdata, s=100, cmap=pylab.cm.jet)
#    fig.colorbar(scat)

#    ani = animation.FuncAnimation(fig, update_scatter, frames=xrange(len(data)),
#                                  fargs=(data, xvar, yvar, zvar, scat))
#    plt.show()

    fig = plt.figure()
    for frame in data:
        pylab.title(zvar + ' at Time ' + str(frame['time']))
        pylab.xlabel(xvar)
        pylab.ylabel(yvar)
        pylab.ylim([0, 1])
        pylab.xlim([0, 1])
        xdata = [p[xvar] for p in frame['particles']]
        ydata = [p[yvar] for p in frame['particles']]
        if (zvar == 'p'):
            cdata = [-0.5*(p['sxx']+p['syy']) for p in frame['particles']]
        else:
            cdata = [p[zvar] for p in frame['particles']]
        scat = plt.scatter(xdata, ydata, c=cdata, s=10, cmap=pylab.cm.jet, edgecolors='none')
        fig.colorbar(scat)
        pylab.draw()
        time.sleep(0.05)
        # raw_input()
        pylab.clf()

else:
    frame = data[framevar]
    fig = plt.figure()
    pylab.title(zvar + ' at Time ' + str(frame['time']))
    pylab.xlabel(xvar)
    pylab.ylabel(yvar)
    xdata = [p[xvar] for p in frame['particles']]
    ydata = [p[yvar] for p in frame['particles']]
    if (zvar == 'p'):
        cdata = [-0.5*(p['sxx']+p['syy']) for p in frame['particles']]
    else:
        cdata = [p[zvar] for p in frame['particles']]
    scat = plt.scatter(xdata, ydata, c=cdata, s=100, cmap=pylab.cm.jet, edgecolors='none')
    fig.colorbar(scat)
raw_input()

