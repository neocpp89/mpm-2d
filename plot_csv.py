#!/usr/bin/env python
import sys
import csv
import pylab
import math
import matplotlib
import numpy

if (len(sys.argv) <= 2):
    print 'usage: CSV_FILE plotvariable'
    exit(127)

infile = sys.argv[1]
pvar = sys.argv[2]
framevar = -1

if len(sys.argv) > 3:
    framevar = int(sys.argv[3])

print 'Using input file:', infile

xvar = 'x'
yvar = 'y'

phandle, = pylab.plot([], [], marker='o')

pylab.xlabel(xvar)
pylab.ylabel(yvar)

print 'Plotting', yvar, '(Y axis) vs', xvar, '(X axis).'

with open(infile, 'r') as f_in:
    reader = csv.reader(f_in)
    # get header
    header = reader.next()
    print header

    # any of these will throw if we can't find the key
    idx = header.index(pvar)
    active = header.index('active')
    xidx = header.index(xvar)
    yidx = header.index(yvar)

    tups = map(lambda row: (float(row[xidx]), float(row[yidx]), float(row[idx]), bool((float(row[active])))), reader)

tups = filter(lambda t: t[3], tups)
x,y,c,a = zip(*tups)

#matplotlib.cm.get_cmap('summer')
cm = matplotlib.cm.get_cmap("YlOrBr", 100) #generate a jet map with 10 values
v = cm(numpy.arange(100)) #extract those values as an array
tcm = matplotlib.colors.LinearSegmentedColormap.from_list("truncPuBuGn", v[30:]) 

pylab.scatter(x,y,s=30,c=c, cmap=tcm, edgecolor='none', linewidth=1)
pylab.title('Framefile: ' + infile)
pylab.draw()
pylab.colorbar()
pylab.xlim(0,1)
pylab.ylim(0,1)
pylab.savefig('.'.join(infile.split('.')[:-1]) + '.png', dpi=100)

print "Bye."

