#!/usr/bin/env python
import sys
import numpy
import pylab
import gc

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
                gc.collect()
                if (framenum >= 100):
					break;
                if (len(data) * sys.getsizeof(data[0]['particles']) >= 128 * 1048576):
                    break;
                else:
                    print 'Data size is approx ', len(data) * sys.getsizeof(data[0]['particles']), ' bytes.'

print 'Got', len(data), 'frames of data.'

pylab.ion()
phandle, = pylab.plot([], [], marker='o')
p_analytical, = pylab.plot([], [], color='k', linestyle='--')
# pscatter, = pylab.scatter([], [], c=

pvars = plotvars.split(',')
xvar = pvars[0]
yvar = pvars[1]
if (len(pvars) >= 3):
    zvar = pvars[3]

pylab.xlabel(xvar)
pylab.ylabel(yvar)

print 'Plotting', yvar, '(Y axis) vs', xvar, '(X axis).'

if framevar == -1:
    allxdata = reduce(lambda x,y: x+y, map(lambda x: map(lambda y: y[xvar], x['particles']), data))
    allydata = reduce(lambda x,y: x+y, map(lambda x: map(lambda y: y[yvar], x['particles']), data))
    xlims = [min(allxdata), max(allxdata)]
    ylims = [min(allydata), max(allydata)]

    for frame in data:
        pylab.title('Time ' + str(frame['time']))
        xdata = map(lambda x: x[xvar], frame['particles'])
        ydata = map(lambda x: x[yvar], frame['particles'])
        phandle.set_xdata(xdata)
        phandle.set_ydata(ydata)
        p_analytical.set_xdata([min(xdata), max(xdata)])
        p_analytical.set_ydata([-1*1500*max(xdata), 0])
    #    print 'error =', ydata[0]+1*1500*max(xdata)

        pylab.ylim(ylims)
        pylab.xlim(xlims)

#        pylab.ylim([min(ydata), max(ydata)])
#        pylab.xlim([min(xdata), max(xdata)])

        pylab.draw()
else:
    frame = data[framevar]
    pylab.title('Time ' + str(frame['time']))
    xdata = map(lambda x: x[xvar], frame['particles'])
    ydata = map(lambda x: x[yvar], frame['particles'])
    phandle.set_xdata(xdata)
    phandle.set_ydata(ydata)
    pylab.ylim([min(ydata), max(ydata)])
    pylab.xlim([min(xdata), max(xdata)])
    pylab.draw()

raw_input()

