#!/usr/bin/env python
import sys
import numpy
import pylab
import gc
import math

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
i = 0
framenum = 0

data = []
pvars = plotvars.split(',')
xvar = pvars[0]
yvar = pvars[1]

pylab.ion()
phandle, = pylab.plot([], [], marker='o')

pylab.xlabel(xvar)
pylab.ylabel(yvar)

print 'Plotting', yvar, '(Y axis) vs', xvar, '(X axis).'

with open(infile, 'r') as f_in:
    for line in f_in:
        if (state == FRAMEDESC):
            tok = line.split(' ')
            particles_left = int(tok[2])
            state = PARTICLEDESC
            data.append({'time':float(tok[1]), 'xplot':[0]*particles_left, 'yplot':[0]*particles_left})
            i = 0
        elif (state == PARTICLEDESC):
            particles_left = particles_left - 1
            tok = line.split(' ')
            tokf = map(float, tok);
            particle = {'m':tokf[0], 'v':tokf[1], 'x':tokf[2], 'y':tokf[3],
                        'x_t':tokf[4], 'y_t':tokf[5], 'sxx':tokf[6],
                        'sxy':tokf[7], 'syy':tokf[8], 'ux':tokf[9],
                        'uy':tokf[10], 'gammap':tokf[11]}
            data[framenum]['xplot'][i] = particle[xvar]
            if (yvar == 'pressure'):
                data[framenum]['yplot'][i] = -0.5 * (particle['sxx']+particle['syy'])
            elif (yvar == 'tau'):
                pressure = -0.5 * (particle['sxx']+particle['syy'])
                s0_xx = particle['sxx'] + pressure
                s0_yy = particle['syy'] + pressure
                s0_xy = particle['sxy']
                data[framenum]['yplot'][i] = math.sqrt(0.5 * (s0_xx ** 2 + 2*s0_xy ** 2 + s0_yy ** 2))
            elif (yvar == 'mu'):
                pressure = -0.5 * (particle['sxx']+particle['syy'])
                s0_xx = particle['sxx'] + pressure
                s0_yy = particle['syy'] + pressure
                s0_xy = particle['sxy']
                tau = math.sqrt(0.5 * (s0_xx ** 2 + 2*s0_xy ** 2 + s0_yy ** 2))
                if (pressure == 0):
                    data[framenum]['yplot'][i] = 0
                else:
                    data[framenum]['yplot'][i] = float(tau) / float(pressure)
#                    print tau, pressure
            else:
                data[framenum]['yplot'][i] = particle[yvar]
            i = i + 1
            if (particles_left <= 0):
                state = FRAMEDESC
                print "Read Frame:", framenum
                framenum = framenum + 1
                if (framevar != -1 and framenum >= framevar):
                    break

                phandle.set_xdata(data[framenum-1]['xplot'])
                phandle.set_ydata(data[framenum-1]['yplot'])
                pylab.ylim([min(data[framenum-1]['yplot']),max(data[framenum-1]['yplot'])])
                pylab.xlim([min(data[framenum-1]['xplot']),max(data[framenum-1]['xplot'])])

                pylab.draw()
#                gc.collect()
#                if (framenum >= 100):
#                    break;
#                if (len(data) * sys.getsizeof(data[0]['particles']) >= 128 * 1048576):
#                    break;
#                else:
#                    print 'Data size is approx ', len(data) * sys.getsizeof(data[0]['particles']), ' bytes.'

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

