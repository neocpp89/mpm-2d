#!/usr/bin/env python
import sys
import os
try:
    import Image
except:
    from PIL import Image
import numpy as np

# linear material points per pixel (linear, so 2 -> 4 points in 1 pixel)
lmpp = 2
rho = 1500
# cell spacing for material point
cs = 1.0 / lmpp
s = map(lambda k: (0.5 + k) * cs, range(0, lmpp))
xy_s = [(x,y) for x in s for y in s]

if (len(sys.argv) > 2):
    imgfile = Image.open(sys.argv[-2])
    outfile = open(sys.argv[-1], "w")
else:
    print sys.argv[0], "imagefile outfile"
    exit(0)

# make image black/white
im = imgfile.convert("1")

print "File", sys.argv[-2], "has size", im.size, "pixels."

# pixel spacing
dx = float(1) / im.size[0]
dy = float(1) / im.size[1]

ij_array = [(i, j) for i in xrange(0, im.size[0]) for j in xrange(0, im.size[1])]
pixel_array = map(im.getpixel, ij_array)

# if pixel is black, we have material points in that location
filled_array = [ij_array[idx] for idx, pixel in enumerate(pixel_array) if pixel == 0]

#print filled_array

# i and j are measured from top left of image (i column, j row).
# mpm_2d measures from bottom left, so convert coordinates.
# xy nodes now has a tuple containing the coordinates of the bottom left node
# of the filled element
xy_nodes = map(lambda ij:
                ((float(ij[0])/im.size[0]),(1.0-(float(ij[1])+1)/im.size[1])),
                filled_array)

xy_material_points = [(n[0]+dx*sp[0], n[1]+dy*sp[1]) for n in xy_nodes for sp in xy_s]

#print xy_nodes
#print xy_material_points
#sys.exit(0)

material_points = map(lambda xy: {
                        'm':cs*cs*rho*dx*dy, 'v':cs*cs*dx*dy,
                        'x':xy[0], 'y':xy[1],
                        'x_t':0, 'y_t':0,
                        'sxx':0, 'sxy':0, 'syy':0}, xy_material_points)

print "Have", len(material_points), "material points."
print "Writing to file", sys.argv[-1]
outfile.write("%d\n" % len(material_points))
for point in material_points:
    outfile.write("%g %g %g %g %g %g %g %g %g\n" % (point['m'], point['v'], point['x'], point['y'], point['x_t'], point['y_t'], point['sxx'], point['sxy'], point['syy']))
