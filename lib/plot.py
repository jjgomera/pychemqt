#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2023, Juan José Gómez Romera <jjgomera@gmail.com>

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
from matplotlib import style
from numpy import arange

from tools.qt import QtWidgets
from lib.config import Preferences

# Load style defined in preferences
sty = Preferences.getint("Plot", 'style')
if sty == 0:
    style.use("default")
else:
    style.use(style.available[Preferences.getint("Plot", 'style')-1])

rcParams['font.size'] = '9'


class PlotWidget(backend_qtagg.FigureCanvasQTAgg):
    """QWidget with matplotlib integration"""
    def __init__(self, dim=2, width=15, height=5, dpi=100, parent=None):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super().__init__(self.fig)

        self.dim = dim
        self.setParent(parent)
        self.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding)

        if dim == 2:
            self.ax = self.fig.add_subplot(111)

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


class ConfPlot(QtWidgets.QDialog):
    """Matplotlib configuration"""
    def __init__(self, config=None, parent=None):
        super().__init__(parent)

        # TODO: Add support for custom Rcparams

        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Matplotlib Style:")), 1, 1)
        self.style = QtWidgets.QComboBox()
        layout.addWidget(self.style, 1, 2)
        self.style.addItem("default")
        for sty in style.available:
            self.style.addItem(sty)
        self.style.currentTextChanged.connect(self.updatePlot)

        layout.addItem(QtWidgets.QSpacerItem(
            10, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 10, 1, 1, 3)

        if config.has_section("Plot"):
            self.style.setCurrentIndex(config.getint("Plot", 'style'))

        self.updatePlot()

    def updatePlot(self, textStyle=None):
        """Update style of example plot"""
        if textStyle is None:
            textStyle = self.style.currentText()

        x = arange(-2, 8, .1)
        y = .1 * x ** 3 - x ** 2 + 3 * x + 2

        with style.context(textStyle):
            self.plot = PlotWidget(width=2, height=1)
            self.plot.plot(x, y)
            self.layout().addWidget(self.plot, 9, 1, 1, 3)

    def value(self, config):
        """Update ConfigParser instance with the config"""
        if not config.has_section("Plot"):
            config.add_section("Plot")
        config.set("Plot", "style", str(self.style.currentIndex()))

        return config


if __name__ == '__main__':
    import sys
    t = [0.3, 0.45, 1., 1.5, 3.5, 7.5, 11.0, 24.0]
    k = [0.5, 0.6, 0.75, 0.8, 0.9, 0.95, 0.97, 0.99]

    app = QtWidgets.QApplication(sys.argv)
    grafico = PlotDialog()
    grafico.data(t, k, 'ro')
    grafico.show()
    sys.exit(app.exec())
