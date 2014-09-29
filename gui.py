#!/usr/bin/env python
import matplotlib
matplotlib.use('TkAgg')

from numpy import arange, sin, pi
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
# implement the default mpl key bindings
from matplotlib.backend_bases import key_press_handler

from matplotlib.figure import Figure

import matplotlib.patches as patches

import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk

class View:
    def __init__(self, master=None):
        self.master = master;
        if (master is not None):
            self.createWidgets()
        return

    def createWidgets(self):
        return

import tkFileDialog as Tkfd

root = Tk.Tk()
root.wm_title("mpm2d GUI")

Nn = 11
Ne = Nn - 1
dx = 1.0 / (Nn - 1)
dy = 1.0 / (Nn - 1)
nodes = [x*Nn + y for x in range(0, Nn) for y in range(0, Nn)]
nx = [dx*(n % Nn) for n in nodes]
ny = [dy*(n / Nn) for n in nodes]

elems = [x*Ne + y for x in range(0, Ne) for y in range(0, Ne)]
ex = [dx*(e % Ne + 0.5) for e in elems]
ey = [dy*(e / Ne + 0.5) for e in elems]

f = Figure(figsize=(5,4), dpi=100)
a = f.add_subplot(111)

#a.scatter(nx, ny, nodes)

qq = zip(nx, ny, nodes)
rr = zip(ex, ey, nodes)

for q in qq:
    a.text(*q, fontsize=8, horizontalalignment='center', verticalalignment='center')

for r in rr:
    a.text(*r, fontsize=8, color='red', horizontalalignment='center', verticalalignment='center')


# a tk.DrawingArea
canvas = FigureCanvasTkAgg(f, master=root)
canvas.show()
canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

toolbar = NavigationToolbar2TkAgg( canvas, root )
toolbar.update()
canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

def on_key_event(event):
    print('you pressed %s'%event.key)
    key_press_handler(event, canvas, toolbar)

canvas.mpl_connect('key_press_event', on_key_event)

def _quit():
    root.quit()     # stops mainloop
    root.destroy()  # this is necessary on Windows to prevent
                    # Fatal Python Error: PyEval_RestoreThread: NULL tstate

button = Tk.Button(master=root, text='Quit', command=_quit)
button.pack(side=Tk.BOTTOM)

Tk.mainloop()
