#!/usr/bin/env python
import sys
import numpy
import time
import os

from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg
import numpy as np

app = QtGui.QApplication([])
mw = QtGui.QMainWindow()
mw.resize(800,800)
view = pg.GraphicsLayoutWidget()  ## GraphicsView with GraphicsLayout inserted by default
mw.setCentralWidget(view)
mw.show()

## create area to add plots
w1 = view.addPlot()
print("Generating data, this takes a few seconds...")

## There are a few different ways we can draw scatter plots; each is optimized for different types of data:

## 1) All spots identical and transform-invariant (top-left plot). 
## In this case we can get a huge performance boost by pre-rendering the spot 
## image and just drawing that image repeatedly.

n = 1
s1 = pg.ScatterPlotItem(size=5, pen=pg.mkPen(None), brush=pg.mkBrush(255, 255, 255, 120))
pos = np.random.normal(size=(2,n), scale=1)
spots = [{'pos': pos[:,i], 'data': 1} for i in range(n)] + [{'pos': [0,0], 'data': 1}]
s1.addPoints(spots)
w1.setXRange(0, 1)
w1.setYRange(0, 1)
w1.addItem(s1)

QtGui.QApplication.processEvents()

## Start Qt event loop unless running in interactive mode.
#if __name__ == '__main__':
#    import sys
#    timer = QtGui.QTimer()
#    timer.timeout.connect(update_scatter)
#    timer.start(1000)
#    QtGui.QApplication.exec_()
#    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
#        QtGui.QApplication.instance().exec_()

def draw_scatter(s, data, framenum, plotvars, crange):
    pvars = plotvars.split(',')
    zvar = pvars[2]

    xdata = [p[pvars[0]] for p in data['particles']]
    ydata = [p[pvars[1]] for p in data['particles']]

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

    cmax = max(cdata)
    cmin = min(cdata)
    if ((cmax - cmin) != 0):
        brushes = [pg.hsvColor(0.7 * (c - cmin)/(cmax-cmin)) for c in cdata]
    else:
        brushes = [pg.hsvColor(0) for c in cdata]
    s.setData(x=xdata, y=ydata, brush=brushes)
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
            data['particles'].append(particle)
            if (particles_left <= 0):
                state = FRAMEDESC
                draw_scatter(s1, data, framenum, plotvars, crange)
                QtGui.QApplication.processEvents()
                print 'Read Frame', framenum, ': t =', data['time']
                framenum = framenum + 1

if (framevar == -1):
    print 'Got', framenum, 'total frames of data. Enter to exit.'
    raw_input()
