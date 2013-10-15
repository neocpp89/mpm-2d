#!/usr/bin/env python

import sys

num_particles = int(sys.argv[1])
stride = int(sys.argv[2])


left_color = [228/255.0, 28/255.0, 28/255.0]
right_color = [153/255.0,153/255.0,153/255.0]

f = open("visualization.cfg", "w")
f.write("override-particle-colors = 1\n")

f.write("color-by-index = {\n")
for i in xrange(0, num_particles):
    c = (i / stride) % 2
    if (c == 0): 
        f.write("    "+str(left_color[0])+", "+str(left_color[1])+", "+str(left_color[2]))
    else:
        f.write("    "+str(right_color[0])+", "+str(right_color[1])+", "+str(right_color[2]))
    if (i < num_particles-1):
        f.write(",\n")

f.write("\n")
f.write("}")
