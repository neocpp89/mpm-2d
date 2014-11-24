#!/usr/bin/env python
import sys
import os
import numpy as np

mp_per_width = 25

width = 0.25
height = 0.15

x_center = 0.5
y_center = 0.25

x_ll = 0.0
y_ll = 0.0

# position the block with the x_ll and y_ll coordinates instead of the center
useLowerLeft = True

# material properties
rho = 1500

# cell spacing for material points (sampling in both directions is the same)
cs = width / mp_per_width
mp_per_height = int(height / cs)
sx = map(lambda k: (0.5 + k) * cs, range(0, mp_per_width))
sy = map(lambda k: (0.5 + k) * cs, range(0, mp_per_height))
xy_s = [(x,y) for x in sx for y in sy]

if (len(sys.argv) > 1):
    outfile = open(sys.argv[-1], "w")
else:
    print sys.argv[0], "outfile"
    exit(0)

if (useLowerLeft):
    disp = (x_ll, y_ll)
else:
    disp = (x_center - width/2.0, y_center - height/2.0)

xy_material_points = [(disp[0]+sp[0], disp[1]+sp[1]) for sp in xy_s]

#print xy_nodes
print xy_material_points
#sys.exit(0)

material_points = map(lambda xy: {
                        'm':cs*cs*rho, 'v':cs*cs,
                        'x':xy[0], 'y':xy[1],
                        'x_t':0, 'y_t':0,
                        'sxx':0, 'sxy':0, 'syy':0}, xy_material_points)

print "Have", len(material_points), "material points."
print "Writing to file:", sys.argv[-1]
outfile.write("%d\n" % len(material_points))
for point in material_points:
    outfile.write("%g %g %g %g %g %g %g %g %g\n" % (
        point['m'], point['v'],
        point['x'], point['y'],
        point['x_t'], point['y_t'],
        point['sxx'], point['sxy'], point['syy'])
    )
print "Done."
