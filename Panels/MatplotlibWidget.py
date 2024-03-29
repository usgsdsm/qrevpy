
from PyQt5 import QtGui, QtCore

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as Canvas
from matplotlib.figure import Figure

class MatplotlibWidget(Canvas):        
    def __init__(self, parent=None, title='Title', xlabel='x label', ylabel='y label', dpi=70, hold=False):
        super(MatplotlibWidget, self).init__(Figure())

        self.setParent(parent)
        self.figure = Figure(figsize=(2, 1),dpi=dpi)
        self.canvas = Canvas(self.figure)
        self.theplot = self.figure.add_subplot(111)        

        self.theplot.set_title(title)
        self.theplot.set_xlabel(xlabel)
        self.theplot.set_ylabel(ylabel)

    def plotDataPoints(self, x, y):
        self.theplot.plot(x,y)
        self.draw()            