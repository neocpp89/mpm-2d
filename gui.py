#!/usr/bin/env python
#from Tkinter import *

#root = Tk();

#def callback(event):
#    print "clicked at", event.x, event.y 

#frame = Frame(root, width=100, height=100)
#frame.bind("<Button-1>", callback)
#frame.pack()

#root.mainloop()

##import matplotlib
##matplotlib.use('TkAgg')

##from numpy import arange, sin, pi
##from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
### implement the default mpl key bindings
##from matplotlib.backend_bases import key_press_handler


##from matplotlib.figure import Figure

##impor
##t sys
##if sys.version_info[0] < 3:
##    import Tkinter as Tk
##else:
##    import tkinter as Tk

##root = Tk.Tk()
##root.wm_title("Embedding in TK")


##f = Figure(figsize=(5,4), dpi=100)
##a = f.add_subplot(111)
##t = arange(0.0,3.0,0.01)
##s = sin(2*pi*t)

##a.plot(t,s)


### a tk.DrawingArea
##canvas = FigureCanvasTkAgg(f, master=root)
##canvas.show()
##canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

##toolbar = NavigationToolbar2TkAgg( canvas, root )
##toolbar.update()
##canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

##def on_key_event(event):
##    print('you pressed %s'%event.key)
##    key_press_handler(event, canvas, toolbar)

##canvas.mpl_connect('key_press_event', on_key_event)

##def _quit():
##    root.quit()     # stops mainloop
##    root.destroy()  # this is necessary on Windows to prevent
##                    # Fatal Python Error: PyEval_RestoreThread: NULL tstate

##button = Tk.Button(master=root, text='Quit', command=_quit)
##button.pack(side=Tk.BOTTOM)

##Tk.mainloop()
### If you put root.destroy() here, it will cause an error if
### the window is closed with the window manager.

#####import matplotlib
#####matplotlib.use('TkAgg')

#####from numpy import arange, sin, pi
#####from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
#####from matplotlib.figure import Figure

#####import sys
#####if sys.version_info[0] < 3:
#####    import Tkinter as Tk
#####else:
#####    import tkinter as Tk

#####def destroy(e): sys.exit()

#####def click_cb(event):
#####    print "clicked at", event.x, event.y 

#####root = Tk.Tk()
#####root.wm_title("Embedding in TK")
######root.bind("<Destroy>", destroy)


#####f = Figure(figsize=(5,4), dpi=100)
#####a = f.add_subplot(111)
#####t = arange(0.0,3.0,0.01)
#####s = sin(2*pi*t)

#####a.plot(t,s)
######a.set_title('Tk embedding')
######a.set_xlabel('X axis label')
######a.set_ylabel('Y label')


###### a tk.DrawingArea
#####canvas = FigureCanvasTkAgg(f, master=root)
#####canvas.show()
#####canvas.get_tk_widget().bind("<Button-1>", click_cb)
#####canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

######toolbar = NavigationToolbar2TkAgg( canvas, root )
######toolbar.update()
#####canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

#####button = Tk.Button(master=root, text='Quit', command=sys.exit)
#####button.pack(side=Tk.BOTTOM)

#####Tk.mainloop()

#######import numpy as np


#######class PointBrowser:
#######    """
#######    Click on a point to select and highlight it -- the data that
#######    generated the point will be shown in the lower axes.  Use the 'n'
#######    and 'p' keys to browse through the next and previous points
#######    """
#######    def __init__(self):
#######        self.lastind = 0

#######        self.text = ax.text(0.05, 0.95, 'selected: none',
#######                            transform=ax.transAxes, va='top')
#######        self.selected,  = ax.plot([xs[0]], [ys[0]], 'o', ms=12, alpha=0.4,
#######                                  color='yellow', visible=False)

#######    def onpress(self, event):
#######        if self.lastind is None: return
#######        if event.key not in ('n', 'p'): return
#######        if event.key=='n': inc = 1
#######        else:  inc = -1


#######        self.lastind += inc
#######        self.lastind = np.clip(self.lastind, 0, len(xs)-1)
#######        self.update()

#######    def onpick(self, event):

#######       if event.artist!=line: return True

#######       N = len(event.ind)
#######       if not N: return True

#######       # the click locations
#######       x = event.mouseevent.xdata
#######       y = event.mouseevent.ydata


#######       distances = np.hypot(x-xs[event.ind], y-ys[event.ind])
#######       indmin = distances.argmin()
#######       dataind = event.ind[indmin]

#######       self.lastind = dataind
#######       self.update()

#######    def update(self):
#######        if self.lastind is None: return

#######        dataind = self.lastind

#######        ax2.cla()
#######        ax2.plot(X[dataind])

#######        ax2.text(0.05, 0.9, 'mu=%1.3f\nsigma=%1.3f'%(xs[dataind], ys[dataind]),
#######                 transform=ax2.transAxes, va='top')
#######        ax2.set_ylim(-0.5, 1.5)
#######        self.selected.set_visible(True)
#######        self.selected.set_data(xs[dataind], ys[dataind])

#######        self.text.set_text('selected: %d'%dataind)
#######        fig.canvas.draw()


#######if __name__ == '__main__':
#######    import matplotlib.pyplot as plt

#######    X = np.random.rand(100, 200)
#######    xs = np.mean(X, axis=1)
#######    ys = np.std(X, axis=1)

#######    fig, (ax, ax2) = plt.subplots(2, 1)
#######    ax.set_title('click on point to plot time series')
#######    line, = ax.plot(xs, ys, 'o', picker=5)  # 5 points tolerance

#######    browser = PointBrowser()

#######    fig.canvas.mpl_connect('pick_event', browser.onpick)
#######    fig.canvas.mpl_connect('key_press_event', browser.onpress)

#######    plt.show()

