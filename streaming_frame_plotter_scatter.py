#!/usr/bin/env python
import sys
import numpy
import time
import os
import math

import matplotlib
#matplotlib.use("WXAgg")

import pylab
import matplotlib.pyplot as plt
#import matplotlib.animation as animation

particlesize = 10
xlims = [0, 1]
ylims = [0, 1]
scaling = 10
fsize = (scaling,scaling*((ylims[1] - ylims[0]) / (xlims[1] - xlims[0])))

def update_scatter(i, data, xvar, yvar, zvar, scat):
    xdata = [p[xvar] for p in data[i]['particles']]
    ydata = [p[yvar] for p in data[i]['particles']]
    cdata = [p[zvar] for p in data[i]['particles']]
    # scat.set_offsets(zip(xdata, ydata))
    scat.set_array(numpy.array(cdata))
    return scat,

def draw_scatter(fig, data, framenum, plotvars, crange):
    pylab.ion()

    pvars = plotvars.split(',')
    xvar = pvars[0]
    yvar = pvars[1]
    if (len(pvars) >= 3):
        zvar = pvars[2]

    if (zvar == 'gammap'):
        pylab.title(r'$\gamma^p$ at Time ' + str(data['time']))
    else:
        pylab.title(zvar + ' at Time ' + str(data['time']))
    pylab.xlabel(xvar)
    pylab.ylabel(yvar)
    pylab.ylim(ylims)
    pylab.xlim(xlims)
    xdata = [p[xvar] for p in data['particles']]
    ydata = [p[yvar] for p in data['particles']]
    if (zvar == 'p'):
        cdata = [-0.5*(p['sxx']+p['syy']) for p in data['particles']]
    elif (zvar == 'mp'):
        pr = [-0.5*(p['sxx']+p['syy']) for p in data['particles']]
        q0 = [numpy.sqrt((p['sxx']-pr[i])**2 + (p['syy']-pr[i])**2 + 2*p['sxy']**2) for i,p in enumerate(data['particles'])]
        mu = [(pr[i] / (q0[i] + 1e-10)) for i,p in enumerate(pr)]
        m = numpy.tan(30.0 * 3.14159 / 180.0)
        cdata = [q0[i] - m*pr[i] - 1e-2 for i,p in enumerate(pr)]
    elif (zvar == 'sxxmyy'):
        cdata = [p['sxx'] - p['syy'] for p in data['particles']]
    elif (zvar == 'umag'):
        cdata = [numpy.sqrt(p['ux']**2 + p['uy']**2) for p in data['particles']]
    else:
        cdata = [p[zvar] for p in data['particles']]
    if len(crange) > 1:
        scat = plt.scatter(xdata, ydata, c=cdata, s=particlesize, cmap=pylab.cm.RdBu, edgecolors='none', vmin=crange[0], vmax=crange[1])
    else:
        scat = plt.scatter(xdata, ydata, c=cdata, s=particlesize, cmap=pylab.cm.RdBu, edgecolors='none')
    plt.axes().set_aspect('equal')
    fig.colorbar(scat)
    pylab.draw()
#    time.sleep(0.05)
    return

if (len(sys.argv) <= 2):
    print 'usage: FRAME_DATA plotvariables [frame]'
    print '       plotvariables are of the form \'xdata,ydata\'.'
    exit(127)

infile = sys.argv[1]
plotvars = sys.argv[2]
framevar = -1
crange = [0]

if len(sys.argv) > 3:
    framevar = int(sys.argv[3])

if len(sys.argv) > 4:
    crange = map(float, sys.argv[4].split(','))
    
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
    fig = plt.figure(figsize=fsize, dpi=100)
    for line in f_in:
        if (state == FRAMEDESC):
            tok = line.split(' ')
            particles_left = int(tok[2])
            state = PARTICLEDESC
            data = {'time':float(tok[1]), 'particles':[]}
        elif (state == PARTICLEDESC):
            particles_left = particles_left - 1
            tok = line.split(' ')
            tokf = map(float, tok);
            if (len(tokf) == 14):
                particle = {'m':tokf[0], 'v':tokf[1], 'x':tokf[2], 'y':tokf[3],
                            'x_t':tokf[4], 'y_t':tokf[5], 'sxx':tokf[6],
                            'sxy':tokf[7], 'syy':tokf[8], 'ux':tokf[9],
                            'uy':tokf[10], 'gammap':tokf[11], 'color':tokf[12],
                            'magEf':tokf[13]}
            elif (len(tokf) == 13):
                particle = {'m':tokf[0], 'v':tokf[1], 'x':tokf[2], 'y':tokf[3],
                            'x_t':tokf[4], 'y_t':tokf[5], 'sxx':tokf[6],
                            'sxy':tokf[7], 'syy':tokf[8], 'ux':tokf[9],
                            'uy':tokf[10], 'gammap':tokf[11], 'color':tokf[12]}
            elif (len(tokf) == 12):
                particle = {'m':tokf[0], 'v':tokf[1], 'x':tokf[2], 'y':tokf[3],
                            'x_t':tokf[4], 'y_t':tokf[5], 'sxx':tokf[6],
                            'sxy':tokf[7], 'syy':tokf[8], 'ux':tokf[9],
                            'uy':tokf[10], 'gammap':tokf[11]}
            elif (len(tokf) > 12):
                particle = {'m':tokf[0], 'v':tokf[1], 'x':tokf[2], 'y':tokf[3],
                            'x_t':tokf[4], 'y_t':tokf[5], 'sxx':tokf[6],
                            'sxy':tokf[7], 'syy':tokf[8], 'ux':tokf[9],
                            'uy':tokf[10], 'gammap':tokf[11]}
            data['particles'].append(particle)
            if (particles_left <= 0):
                state = FRAMEDESC
                pylab.clf()
                if (framenum == framevar or framevar == -1):
                    draw_scatter(fig, data, framenum, plotvars, crange)
                if (framenum == framevar and framevar != -1):
                    pylab.savefig("figs/" + os.path.basename(infile)[:-4]+"_"+str(framenum)+".png", bbox_inches=0)
                    pylab.savefig("figs/" + os.path.basename(infile)[:-4]+"_"+str(framenum)+".svg", bbox_inches=0)
                    break; 
                framenum = framenum + 1

if (framevar == -1):
    print 'Got', framenum, 'total frames of data. Enter to exit.'
    raw_input()
