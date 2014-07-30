#!/usr/bin/python
# -*- coding: utf-8 -*-

from PyQt4 import QtCore, QtGui
from matplotlib import rcParams
rcParams['backend'] = 'QT4Agg'  #Fija el backend de las ventanas de matplotlib a qt4
rcParams['font.size'] = '8'
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg, NavigationToolbar2QT
from pylab import Figure, plot,  title,  figtext,  xlabel,  grid,  show
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import fsolve
from scipy import arange, linspace

from compuestos import Componente
from corriente import Corriente, Mezcla

class mpl(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=15, height=5, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        FigureCanvasQTAgg.__init__(self, self.fig)
        self.setParent(parent)
        self.axes2D = self.fig.add_subplot(111)
        FigureCanvasQTAgg.setSizePolicy(self, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        FigureCanvasQTAgg.updateGeometry(self)

    
    def config(self, xmin=None, xmax=None, ymin=None, ymax=None, scalex="linear", scaley="linear"):
        self.axes2D.clear()
        self.axes2D.set_autoscale_on(False)
        self.axes2D.set_xlabel("Pr", horizontalalignment='right', size='12')
        self.axes2D.set_ylabel("Z", size='14')
        self.axes2D.grid(True)
        if xmin!=None and xmax!=None:
            self.axes2D.set_xlim(xmin, xmax)
        else:
            self.axes2D.set_autoscalex_on(True)
        if ymin!=None and ymax!=None:
            self.axes2D.set_ylim(ymin, ymax)
        else:
            self.axes2D.set_autoscaley_on(True)
        
        self.axes2D.set_xscale(scalex)
        self.axes2D.set_yscale(scaley)
        
    def plot(self, x, y, color="#000000", grosor=1, linestyle="-"):
        self.axes2D.plot(x, y, color=color, lw=grosor, ls=linestyle)
        
    def data(self, *args, **kwargs):
        self.axes2D.plot(*args, **kwargs)
        self.draw()


class matplotlib(FigureCanvasQTAgg):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
    def __init__(self,  dim=2, parent=None):
        self.fig = Figure(figsize=(10, 10), dpi=100)
        FigureCanvasQTAgg.__init__(self, self.fig)
        self.setParent(parent)
        FigureCanvasQTAgg.setSizePolicy(self, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        FigureCanvasQTAgg.updateGeometry(self)
        
        if dim==2:
            self.ax = self.fig.add_subplot(111)
            self.ax.figure.subplots_adjust(left=0.08, right=0.98, bottom=0.08, top=0.92)

        else:
            self.ax = Axes3D(self.fig)

        FigureCanvasQTAgg.setSizePolicy(self, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        FigureCanvasQTAgg.updateGeometry(self)

    def plot_3D(self, labels, xdata, ydata, zdata, config=None):
        """Método que dibuja la matriz de datos"""
        self.ax.clear()
        self.data={"x": xdata[0], "y": ydata[:,0], "z": zdata}
        
        if config and config.getboolean("MEOS", "surface"):
            self.ax.plot_surface(xdata, ydata, zdata, rstride=1, cstride=1)
        else:
            self.ax.plot_wireframe(xdata, ydata, zdata, rstride=1, cstride=1)
            
        self.ax.set_xlabel(labels[0])
        self.ax.set_ylabel(labels[1])
        self.ax.set_zlabel(labels[2])
        self.ax.mouse_init(rotate_btn=1, zoom_btn=2)

#class PlotWidget(QtGui.QWidget):
#    def __init__(self, dim, parent=None):
#        super(PlotWidget, self).__init__(parent)
#        layout=QtGui.QVBoxLayout(self)
#        self.plot=matplotlib(dim)
#        layout.addWidget(self.plot)
#        
#        self.toolbar=NavigationToolbar2QT(self, self)
#        layout.addWidget(self.toolbar)
#        

class Plot(QtGui.QDialog):
    def __init__(self, parent=None):
        super(Plot, self).__init__(parent)
        gridLayout = QtGui.QGridLayout(self)

        self.plot=matplotlib()
        gridLayout.addWidget(self.plot,1,1,1,2)
        self.toolbar=NavigationToolbar2QT(self.plot, self.plot)
        gridLayout.addWidget(self.toolbar,2,1)
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Close)
        self.buttonBox.setSizePolicy(QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Maximum)
        self.buttonBox.rejected.connect(self.reject)
        gridLayout.addWidget(self.buttonBox,2,2)
        
    def addText(self, *args, **kwargs):
        self.plot.ax.text(*args, **kwargs)


    def addData(self, *args, **kwargs):
        self.plot.ax.plot(*args, **kwargs)
#        self.plot.draw()

    

if __name__ == '__main__':
    import sys
    t=[0.3, 0.45, 1., 1.5, 3.5, 7.5, 11.0, 24.0]
    k=[0.5, 0.6, 0.75, 0.8, 0.9, 0.95, 0.97, 0.99]
    
    app = QtGui.QApplication(sys.argv)
    grafico=Plot()
    grafico.data(t, k, 'ro')
    grafico.show()
    sys.exit(app.exec_())
