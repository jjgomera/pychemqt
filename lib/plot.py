#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2016, Juan José Gómez Romera <jjgomera@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''



from PyQt5 import QtCore, QtWidgets

from matplotlib import rcParams
rcParams['backend'] = 'QT5Agg'  #Fija el backend de las ventanas de matplotlib a qt5
rcParams['font.size'] = '9'
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from pylab import Figure
from mpl_toolkits.mplot3d import Axes3D


class mpl(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=15, height=5, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        FigureCanvasQTAgg.__init__(self, self.fig)
        self.setParent(parent)
        self.ax = self.fig.add_subplot(111)
        # self.ax.figure.subplots_adjust(left=0.08, right=0.92, bottom=0.1, top=0.95)
        FigureCanvasQTAgg.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        FigureCanvasQTAgg.updateGeometry(self)


    def config(self, xmin=None, xmax=None, ymin=None, ymax=None, scalex="linear", scaley="linear"):
        self.ax.clear()
        self.ax.set_autoscale_on(False)
        self.ax.set_xlabel("Pr", horizontalalignment='right', size='12')
        self.ax.set_ylabel("Z", size='14')
        self.ax.grid(True)
        if xmin!=None and xmax!=None:
            self.ax.set_xlim(xmin, xmax)
        else:
            self.ax.set_autoscalex_on(True)
        if ymin!=None and ymax!=None:
            self.ax.set_ylim(ymin, ymax)
        else:
            self.ax.set_autoscaley_on(True)

        self.ax.set_xscale(scalex)
        self.ax.set_yscale(scaley)

    def plot(self, *args, **kwargs):
        self.ax.plot(*args, **kwargs)

    def data(self, *args, **kwargs):
        self.ax.plot(*args, **kwargs)
        self.draw()


class matplotlib(FigureCanvasQTAgg):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
    def __init__(self,  dim=2, parent=None):
        self.fig = Figure(figsize=(10, 10), dpi=100)
        self.dim = dim
        FigureCanvasQTAgg.__init__(self, self.fig)
        self.setParent(parent)
        FigureCanvasQTAgg.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        FigureCanvasQTAgg.updateGeometry(self)

        if dim==2:
            self.ax = self.fig.add_subplot(111)
            self.ax.figure.subplots_adjust(left=0.08, right=0.98, bottom=0.08, top=0.92)

        else:
            self.ax = Axes3D(self.fig)
            self.ax.mouse_init(rotate_btn=1, zoom_btn=2)

        FigureCanvasQTAgg.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
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

class Plot(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super(Plot, self).__init__(parent)
        gridLayout = QtWidgets.QGridLayout(self)

        self.plot=matplotlib()
        gridLayout.addWidget(self.plot,1,1,1,2)
        self.toolbar=NavigationToolbar2QT(self.plot, self.plot)
        gridLayout.addWidget(self.toolbar,2,1)
        self.buttonBox = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Close)
        self.buttonBox.setSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Maximum)
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

    app = QtWidgets.QApplication(sys.argv)
    grafico=Plot()
    grafico.data(t, k, 'ro')
    grafico.show()
    sys.exit(app.exec_())
