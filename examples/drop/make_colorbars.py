#!/usr/bin/env python
import matplotlib as mpl
from matplotlib import pyplot
import numpy

# load defaults for figures
mpl.rc_file(r'./mpl.rc')

lightgray = [0.95,0.95,0.95]

# gammadotp small and large colorbars
fig = pyplot.figure(figsize=(4,0.5))
ax1 = fig.add_axes([0.05, 0.75, 0.9, 0.20])
cm = mpl.cm.get_cmap("bone", 100)
cb1 = mpl.colorbar.ColorbarBase(ax1, extend='both', cmap=cm, norm=mpl.colors.Normalize(vmin=-2.0, vmax=1.0, clip = False), orientation='horizontal')
cb1.set_label(r'$\dot{\bar{\gamma}}^p$ [$\mathrm{s}^{-1}$]')
cb1.solids.set_rasterized(True)

from matplotlib import ticker

tick_locator = ticker.MaxNLocator(nbins=4)
cb1.locator = tick_locator
cb1.update_ticks()
cb1.ax.set_xticklabels(map(lambda x: '$10^{{{}}}$'.format(x), range(-2,2)))

pyplot.savefig('bone_colorbar.png')

