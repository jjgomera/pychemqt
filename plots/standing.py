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
along with this program.  If not, see <http://www.gnu.org/licenses/>.


.. include:: standing.rst


The module include all related moody chart functionality
    * :class:`Standing_Katz`: Chart dialog
    * :func:`calculate`: Calculate procedure
    * :class:`CalculateDialog`: Dialog to calculate a specified point

and its configuration

    * :class:`Config`: Standing-Katz chart configuration
'''


import json
import os
import re

from numpy import arange
from scipy.optimize import fsolve

from lib.config import conf_dir, Preferences
from lib.crude import Z_list
from lib.plot import PlotWidget
from lib.utilities import formatLine
from tools.qt import QtCore, QtGui, QtWidgets
from UI.widgets import Entrada_con_unidades, GridConfig, LineConfig

from plots.ui import Chart


def calculate(config, dat=None):
    """Plot calculate procedure

    Parameters
    ----------
    config : Configparser
        pychemqt configparser configuration instance
    dat : dict
        dict with other method data

    Returns
    -------
    dat : dict
        dict with method data

    Notes
    -----
    This procedure is called when a new method of calculation of Z is
    necessary. Add the new data to the input dat with the other Method
    calculated yet
    """
    method = config.getint("Standing_Katz", "method")
    Tr = map(float, config.get("Standing_Katz", "Tr").split(","))
    Z = Z_list[method]

    Pr = arange(1e-3, 15, 0.01)

    if dat is None:
        dat = {}

    lines = {}
    for t in Tr:
        pr = []
        z = []
        for p in Pr:
            try:
                z.append(Z(t, p))
                pr.append(p)
            except NotImplementedError:
                pass
            except ValueError:
                pass

        # Only save the lines with data
        if pr:
            lines[t] = {"Pr": pr, "Z": z}

    dat[method] = lines

    # Save to file
    with open(conf_dir+"standing_katz.dat", "w", encoding="utf-8") as file:
        json.dump(dat, file, indent=4)


class Config(QtWidgets.QWidget):
    """Standing-Katz chart configuration"""
    TITLE = QtWidgets.QApplication.translate("Standing_Katz", "Standing-Katz chart")
    TITLECONFIG = QtWidgets.QApplication.translate("Standing_Katz", "Standing-Katz chart configuration")

    def __init__(self, config=None, parent=None):
        super().__init__(parent)
        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(QtWidgets.QLabel(self.tr("Method:")), 1, 1)
        self.method = QtWidgets.QComboBox()
        for Z in Z_list:
            name = Z.__name__[2:].replace("_", "-")
            year = re.search(r"((\d+))", Z.__doc__).group(0)
            doc = f"{name} ({year})"
            self.method.addItem(doc)
        layout.addWidget(self.method, 1, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed), 1, 3)

        layout.addWidget(QtWidgets.QLabel("Tr:"), 4, 1)
        self.Tr = QtWidgets.QLineEdit()
        layout.addWidget(self.Tr, 4, 2, 1, 2)
        self.lineconfig = LineConfig(
            "line", self.tr("Reduced temperature style line"))
        layout.addWidget(self.lineconfig, 5, 1, 1, 3)

        self.cruxconfig = LineConfig("crux", self.tr("Crux style line"))
        layout.addWidget(self.cruxconfig, 6, 1, 1, 3)

        self.gridconfig = GridConfig("grid", self.tr("Grid style line"))
        layout.addWidget(self.gridconfig, 7, 1, 1, 3)

        layout.addItem(QtWidgets.QSpacerItem(
            10, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 10, 1, 1, 3)

        if config and config.has_section("Standing_Katz"):
            self.method.setCurrentIndex(config.getint(
                "Standing_Katz", "method"))
            self.Tr.setText(config.get("Standing_Katz", "Tr"))
            self.lineconfig.setConfig(config, "Standing_Katz")
            self.cruxconfig.setConfig(config, "Standing_Katz")
            self.gridconfig.setConfig(config, "Standing_Katz")

    def value(self, config):
        """Update ConfigParser instance with the config"""
        if not config.has_section("Standing_Katz"):
            config.add_section("Standing_Katz")
        config.set("Standing_Katz", "method", str(self.method.currentIndex()))
        config.set("Standing_Katz", "Tr", self.Tr.text())
        config = self.lineconfig.value(config, "Standing_Katz")
        config = self.cruxconfig.value(config, "Standing_Katz")
        config = self.gridconfig.value(config, "Standing_Katz")
        return config


class Standing_Katz(Chart):
    """Standing-Katz chart dialog"""
    title = QtWidgets.QApplication.translate("Standing_Katz", "Standing and Katz compressivitity factors "
               "chart for natural gas")
    widgetConfig = Config
    locLogo = (0.8, 0.12, 0.1, 0.1)
    note = None

    def customUI(self):
        """Define custom UI element"""
        self.butonPNG.clicked.disconnect()
        self.butonPNG.clicked.connect(self.savePNG)

        self.plt2 = PlotWidget(parent=self, dpi=90)
        self.plt2.fig.canvas.mpl_connect('button_press_event', self.click)
        self.plotWidget.layout().addWidget(self.plt2, 1, 1)
        self.setMask()
        self.set_logo(self.plt2)

    def savePNG(self):
        """Save chart image to png file"""
        fmt = "Portable Network Graphics (*.png)"
        fname, ext = QtWidgets.QFileDialog.getSaveFileName(
            self, self.tr("Save chart to file"), "./", fmt)
        if fname and ext == fmt:
            if fname.split(".")[-1] != "png":
                fname += ".png"
            pix = self.plotWidget.grab()
            pix.save(fname, "png")

    def setMask(self):
        """Mask both plot to show only the region useful"""
        w = self.plt.width()
        h = self.plt.height()

        pol = QtGui.QPolygon()
        pol.append(QtCore.QPoint(int(0.05*w), int(0.05*h)))
        pol.append(QtCore.QPoint(int(w), int(0.05*h)))
        pol.append(QtCore.QPoint(int(w), int(0.23*h)))
        pol.append(QtCore.QPoint(int(0.9*w), int(0.228*h)))
        pol.append(QtCore.QPoint(int(0.3*w), int(0.698*h)))
        pol.append(QtCore.QPoint(int(0.125*w), int(0.698*h)))
        pol.append(QtCore.QPoint(int(0.125*w), int(0.7*h)))
        pol.append(QtCore.QPoint(int(0.05*w), int(0.7*h)))
        reg = QtGui.QRegion(pol)
        self.plt.setMask(reg)

        pol = QtGui.QPolygon()
        pol.append(QtCore.QPoint(int(0.05*w), int(0.7*h)))
        pol.append(QtCore.QPoint(int(0.3*w), int(0.7*h)))
        pol.append(QtCore.QPoint(int(0.90*w), int(0.23*h)))
        pol.append(QtCore.QPoint(int(0.90*w), int(0.228*h)))
        pol.append(QtCore.QPoint(int(w), int(0.228*h)))
        pol.append(QtCore.QPoint(int(w), int(h)))
        pol.append(QtCore.QPoint(int(0.05*w), int(h)))
        reg = QtGui.QRegion(pol)
        self.plt2.setMask(reg)

        x = (0, 1.8, 8)
        y = (0.276, 0.276, 0.95)
        self.plt.ax.plot(x, y, "black", lw=0.5)
        x = (7, 8.8, 15)
        y = (1.17, 1.17, 1.84)
        self.plt2.ax.plot(x, y, "black", lw=0.5)

        self.plt.draw()
        self.plt2.draw()

    def paintEvent(self, event):
        """Do redraw in each change of window size or position"""
        self.setMask()
        Chart.paintEvent(self, event)

    def click(self, event):
        """Update input and graph annotate when mouse click over chart"""
        Pr = event.xdata
        Z = event.ydata

        # Exit if click event if out of axis
        if Pr is None:
            self.clearCrux()
            return

        method = Preferences.getint("Standing_Katz", "method")
        f_Z = Z_list[method]

        def f(Tr):
            return f_Z(Tr, Pr) - Z

        Tr = 0
        for to in (1, 2, 3):
            try:
                rinput = fsolve(f, to, full_output=True)
                if rinput[2] == 1:
                    Tr = rinput[0][0]
                    break
            except ValueError:
                continue

        self.createCrux(Tr, Pr, Z)

    def plot(self):
        """Plot the Standing-Katz chart using the indicate method """
        method = Preferences.get("Standing_Katz", "method")

        self.plt.ax.clear()
        self.plt.ax.set_xlim(0, 8)
        self.plt.ax.set_ylim(0, 1.1)
        self.plt.ax.set_xlabel(r"$P_r=\frac{P}{P_c}$", ha='center', size='14')
        self.plt.ax.set_ylabel(r"$Z=\frac{PV}{nRT}$", va="bottom", size='14')

        self.plt2.ax.clear()
        self.plt2.ax.set_xlim(7, 15)
        self.plt2.ax.set_ylim(0.9, 2)
        self.plt2.ax.set_xlabel(r"$P_r=\frac{P}{P_c}$", ha='center', size='14')

        grid = Preferences.getboolean("Standing_Katz", "grid")
        kw = formatLine(Preferences, "Standing_Katz", "grid")
        del kw["marker"]
        if grid:
            self.plt.ax.grid(grid, **kw)
            self.plt2.ax.grid(grid, **kw)
        else:
            self.plt.ax.grid(grid)
            self.plt2.ax.grid(grid)

        if not os.path.isfile(conf_dir+"standing_katz.dat"):
            calculate(Preferences)

        load = False
        with open(conf_dir+"standing_katz.dat", "r", encoding="utf-8") as file:
            try:
                dat = json.load(file)
            except ValueError:
                calculate(Preferences)
                load = True

            if method not in dat:
                calculate(Preferences, dat)
                load = True

        # Reload file if it's created in last with statement
        if load:
            fname = conf_dir+"standing_katz.dat"
            with open(fname, "r", encoding="utf-8") as file:
                dat = json.load(file)

        # Define Crux
        kw = formatLine(Preferences, "Standing_Katz", "crux")
        self.plt.lx = self.plt.ax.axhline(**kw)  # the horiz line
        self.plt.ly = self.plt.ax.axvline(**kw)  # the vert line
        self.plt.lx.set_visible(False)
        self.plt.ly.set_visible(False)
        self.plt2.lx = self.plt2.ax.axhline(**kw)  # the horiz line
        self.plt2.ly = self.plt2.ax.axvline(**kw)  # the vert line
        self.plt2.lx.set_visible(False)
        self.plt2.ly.set_visible(False)
        self.plt.ax.xaxis.tick_top()
        self.plt2.ax.yaxis.tick_right()
        self.note = None

        # Plot data
        kw = formatLine(Preferences, "Standing_Katz", "line")
        pzmin = 4
        for Tr in sorted(dat[str(method)].keys()):
            line = dat[str(method)][Tr]
            self.plt.ax.plot(line["Pr"], line["Z"], **kw)
            self.plt2.ax.plot(line["Pr"], line["Z"], **kw)

            # Add Tr legend
            # Position as possible at minimum position of line
            zmin = min(line["Z"])
            if zmin < 1:
                pzmin = line["Pr"][line["Z"].index(zmin)]
            else:
                if 4 not in line["Pr"]:
                    line["Pr"].append(pzmin)
                    line["Pr"].sort()
                zmin = line["Z"][line["Pr"].index(pzmin)]
            self.plt.ax.text(pzmin, zmin, str(Tr),
                             size="x-small", ha='left', va='bottom')
            self.plt2.ax.text(line["Pr"][-1], line["Z"][-1], str(Tr),
                              size="x-small", ha='right', va='bottom')

        # Add explicative legend of isoline
        self.plt.ax.text(3, 1.07, r"$T_r$", size="12",
                         ha='left', va='center')

        self.plt.draw()

    def calculate(self):
        dlg = CalculateDialog()
        if dlg.exec():
            Tr = dlg.Tr.value
            Pr = dlg.Pr.value
            Z = dlg.Z.value
            self.createCrux(Tr, Pr, Z)

    def clearCrux(self):
        """Delete crux and note text"""
        if self.note:
            self.note.remove()
            self.note = None
        self.plt.lx.set_visible(False)
        self.plt.ly.set_visible(False)
        self.plt2.lx.set_visible(False)
        self.plt2.ly.set_visible(False)
        self.plt.draw()
        self.plt2.draw()

    def createCrux(self, Tr, Pr, Z):
        """Create a crux in selected point of plot and show data at bottom
        right corner"""

        self.clearCrux()

        if Pr < 8:
            self.plt.lx.set_visible(True)
            self.plt.ly.set_visible(True)
            self.plt.ly.set_xdata(Pr)
            self.plt.lx.set_ydata(Z)

        if Pr > 7:
            self.plt2.lx.set_visible(True)
            self.plt2.ly.set_visible(True)
            self.plt2.ly.set_xdata(Pr)
            self.plt2.lx.set_ydata(Z)

        if Tr:
            txt = f"$T_r$: {Tr:0.4g}\n$P_r$: {Pr:0.4g}\nZ: {Z:0.4g}"
        else:
            Tr = self.tr("Not converged")
            txt = f"$T_r$: {Tr}\n$P_r$: {Pr:0.4g}\nZ: {Z:0.4g}"
        self.note = self.plt2.fig.text(0.92, 0.05, txt, size="8", ha="right")

        self.plt.draw()
        self.plt2.draw()


class CalculateDialog(QtWidgets.QDialog):
    """Dialog to calculate a specified point"""
    def __init__(self, parent=None):
        super().__init__(parent)
        title = self.tr("Calculate compressibility factor of natural gas")
        self.setWindowTitle(title)
        layout = QtWidgets.QGridLayout(self)
        label = QtWidgets.QLabel(self.tr("Method:"))
        layout.addWidget(label, 1, 0)
        self.method = QtWidgets.QComboBox()
        for Z in Z_list:
            name = Z.__name__[2:].replace("_", "-")
            year = re.search(r"((\d+))", Z.__doc__).group(0)
            doc = f"{name} ({year})"
            self.method.addItem(doc)

        self.method.currentIndexChanged.connect(self.calculate)
        layout.addWidget(self.method, 1, 1, 1, 2)

        layout.addWidget(QtWidgets.QLabel("Tr"), 2, 1)
        self.Tr = Entrada_con_unidades(float, tolerancia=4)
        self.Tr.valueChanged.connect(self.calculate)
        layout.addWidget(self.Tr, 2, 2)
        layout.addWidget(QtWidgets.QLabel("Pr"), 3, 1)
        self.Pr = Entrada_con_unidades(float)
        self.Pr.valueChanged.connect(self.calculate)
        layout.addWidget(self.Pr, 3, 2)
        layout.addWidget(QtWidgets.QLabel("Z"), 4, 1)
        self.Z = Entrada_con_unidades(float, readOnly=True, decimales=8)
        layout.addWidget(self.Z, 4, 2)

        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 9, 1, 1, 3)

        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok
            | QtWidgets.QDialogButtonBox.StandardButton.Close)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox, 10, 1, 1, 3)

    def calculate(self):
        """Calculate point procedure"""
        index = self.method.currentIndex()
        Z = Z_list[index]
        Tr = self.Tr.value
        Pr = self.Pr.value
        if Pr and Tr is not None:
            z = Z(Tr, Pr)
            self.Z.setValue(z)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = Standing_Katz()
    Dialog.show()
    sys.exit(app.exec())
