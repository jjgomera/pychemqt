#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2022, Juan José Gómez Romera <jjgomera@gmail.com>

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


###############################################################################
# Plot library with all matplotlib functionality
###############################################################################


from matplotlib import rcParams
from matplotlib.backends import backend_qtagg
from matplotlib.figure import Figure

from tools.qt import QtWidgets


rcParams['font.size'] = '9'


class PlotWidget(backend_qtagg.FigureCanvasQTAgg):
    """QWidget with matplotlib integration"""
    def __init__(self, dim=2, width=15, height=5, parent=None):
        self.fig = Figure(figsize=(width, height), dpi=100)
        super().__init__(self.fig)

        self.dim = dim
        self.setParent(parent)
        self.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding)

        if dim == 2:
            self.ax = self.fig.add_subplot(111)
            self.ax.figure.subplots_adjust(
                left=0.08, right=0.98, bottom=0.08, top=0.92)

        else:
            self.ax = self.fig.add_subplot(projection="3d")
            self.ax.mouse_init(rotate_btn=1, zoom_btn=2)

    def plot(self, *args, **kwargs):
        """Direct accesst to ax plot procedure"""
        self.ax.plot(*args, **kwargs)

    def savePNG(self):
        """Save chart image to png file"""
        fmt = "Portable Network Graphics (*.png)"
        fname, ext = QtWidgets.QFileDialog.getSaveFileName(
            self,
            QtWidgets.QApplication.translate("pychemqt", "Save chart to file"),
            "./", fmt)
        if fname and ext == fmt:
            if fname.split(".")[-1] != "png":
                fname += ".png"
            self.fig.savefig(fname, facecolor='#fafafa')

    # def plot_3D(self, labels, xdata, ydata, zdata, config=None):
        # """Método que dibuja la matriz de datos"""
        # self.ax.clear()
        # self.data = {"x": xdata[0], "y": ydata[:, 0], "z": zdata}

        # if config and config.getboolean("MEOS", "surface"):
            # self.ax.plot_surface(xdata, ydata, zdata, rstride=1, cstride=1)
        # else:
            # self.ax.plot_wireframe(xdata, ydata, zdata, rstride=1, cstride=1)

        # self.ax.set_xlabel(labels[0])
        # self.ax.set_ylabel(labels[1])
        # self.ax.set_zlabel(labels[2])
        # self.ax.mouse_init(rotate_btn=1, zoom_btn=2)


class PlotDialog(QtWidgets.QDialog):
    """QDialog including Plotwidget, navigationtoolbar and a button to close"""
    def __init__(self, accept=False, cancel=True, parent=None):
        super().__init__(parent)
        gridLayout = QtWidgets.QGridLayout(self)

        self.plot = PlotWidget()
        gridLayout.addWidget(self.plot, 1, 1, 1, 2)
        self.toolbar = backend_qtagg.NavigationToolbar2QT(self.plot, self.plot)
        gridLayout.addWidget(self.toolbar, 2, 1)
        btonBox = QtWidgets.QDialogButtonBox()
        if accept:
            btonBox.addButton(QtWidgets.QDialogButtonBox.StandardButton.Ok)
            btonBox.accepted.connect(self.accept)
        if cancel:
            btonBox.addButton(QtWidgets.QDialogButtonBox.StandardButton.Cancel)
            btonBox.rejected.connect(self.reject)
        btonBox.setSizePolicy(QtWidgets.QSizePolicy.Policy.Maximum,
                              QtWidgets.QSizePolicy.Policy.Maximum)
        gridLayout.addWidget(btonBox, 2, 2)

    def addText(self, *args, **kwargs):
        """Direct access to ax text procedure"""
        self.plot.ax.text(*args, **kwargs)

    def addData(self, *args, **kwargs):
        """Direct access to ax plot procedure"""
        self.plot.ax.plot(*args, **kwargs)


if __name__ == '__main__':
    import sys
    t = [0.3, 0.45, 1., 1.5, 3.5, 7.5, 11.0, 24.0]
    k = [0.5, 0.6, 0.75, 0.8, 0.9, 0.95, 0.97, 0.99]

    app = QtWidgets.QApplication(sys.argv)
    grafico = PlotDialog()
    grafico.data(t, k, 'ro')
    grafico.show()
    sys.exit(app.exec())
