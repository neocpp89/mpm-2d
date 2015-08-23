#!/usr/bin/env python
import sys
import os
import numpy as np

if (len(sys.argv) > 2):
    slope = float(sys.argv[-2])
    particle_filename = sys.argv[-1] + '.particles'
    grid_filename = sys.argv[-1] + '.grid'
else:
    print sys.argv[0], "slope outprefix"
    exit(0)

velmag = 1
rho_intruder = 8910.0

Lx = 0.91
Ly = Lx

const_area = 0.0107 # m^2
w_intruder = 0.0965 # m
h_tip = 0.5*slope*w_intruder
A_tip = 0.5*w_intruder*h_tip
h_added = 0
if (A_tip > const_area):
    print 'Warning, intruder tip alone exceeds area.'
else:
    h_added = (const_area - A_tip) / w_intruder

# linear material points per pixel (linear, so 2 -> 4 points in 1 pixel)
lmpp = 1
Ne = 200
rho_bulk = 1500.0

# cell spacing for material point
cs = Lx / lmpp
s = map(lambda k: (0.5 + k) * cs, range(0, lmpp))

# material point positions inside cell
xy_s = [(x,y) for x in s for y in s]

# pixel spacing
dx = float(Lx) / Ne
dy = float(Ly) / Ne

with open(grid_filename, 'w') as f:
    f.write('{}\n{}\n'.format(Ne+1, Lx))

ij_array = [(i, j) for i in xrange(0, Ne) for j in xrange(0, Ne)]

# i and j are measured from top left of image (i column, j row).
# mpm_2d measures from bottom left, so convert coordinates.
# xy nodes now has a tuple containing the coordinates of the bottom left node
# of the filled element
xy_nodes = map(lambda ij:
                (Lx*(float(ij[0])/Ne), Ly*(1.0-(float(ij[1])+1)/Ne)),
                ij_array)

xy_material_point_candidates = [(n[0]+dx*sp[0], n[1]+dy*sp[1]) for n in xy_nodes for sp in xy_s]

def intruder_intersection(x, y):
    return (   (y > 0) and
    (y < (h_added + h_tip)) and
    ((y - slope*x) > 0) and
    (x < 0.5 * w_intruder)  )

def bulk_intersection(x, y):
    return (y < 0)

def global_xform(xy):
    return (abs(xy[0] - 0.5*Lx), xy[1] - 0.5*Ly)

xy_material_points = filter(lambda xy:
                        intruder_intersection(*global_xform(xy)) or bulk_intersection(*global_xform(xy)),
                        xy_material_point_candidates)

material_points = map(lambda xy: {
                        'm':cs*cs*rho_bulk*dx*dy, 'v':cs*cs*dx*dy,
                        'x':xy[0], 'y':xy[1],
                        'x_t':0, 'y_t':0,
                        'sxx':0, 'sxy':0, 'syy':0}, xy_material_points)

for i, point in enumerate(material_points):
    if point['y'] > 0.5*Lx:
        material_points[i]['y_t'] = -velmag
        material_points[i]['m'] = (rho_intruder / rho_bulk) * material_points[i]['m']

print "Have", len(material_points), "material points."
print "Writing to file", particle_filename
with open(particle_filename, 'w') as f:
    f.write("%d\n" % len(material_points))
    for point in material_points:
        f.write("%g %g %g %g %g %g %g %g %g\n" % (point['m'], point['v'], point['x'], point['y'], point['x_t'], point['y_t'], point['sxx'], point['sxy'], point['syy']))
