#!/usr/bin/env python

from pylab import *
from numpy import outer
rc('text', usetex=False)
a=outer(arange(0,1,0.01),ones(10))
figure(figsize=(10,5))
subplots_adjust(top=0.8,bottom=0.05,left=0.01,right=0.99)
maps=[m for m in cm.datad if not m.endswith("_r")]
maps.sort()
l=len(maps)+1
num_samples = 256
sample_range = linspace(0,1,num_samples)
for i, m in enumerate(maps):
    vals = get_cmap(m, num_samples)
    filename = "cmaps/"+m+".cm"
    print filename
    f = open(filename, 'w')
    for j in xrange(0,num_samples):
        s = "I "+str(float(sample_range[j]))+","+str(int(255*vals(j)[0]))+","+str(int(255*vals(j)[1]))+","+str(int(255*vals(j)[2]))+","+str(int(255*vals(j)[3]))
#        print int(255*vals(j)[0])
#        print int(255*vals(j)[1])
#        print int(255*vals(j)[2])
#        print int(255*vals(j)[3])
#        print sample_range[j]
        print s
        f.write(s+"\n")
#    subplot(1,l,i+1)
#    axis("off")
#    imshow(a,aspect='auto',cmap=get_cmap(m),origin="lower")
#    title(m,rotation=90,fontsize=10)
#savefig("colormaps.png",dpi=100,facecolor='gray')
