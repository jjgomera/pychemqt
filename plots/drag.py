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
# Drag sphere chart
###############################################################################


from tools.qt import QtWidgets
from numpy import logspace

from lib import drag
from lib.utilities import formatLine
from UI.widgets import Entrada_con_unidades, GridConfig, LineConfig

from plots.ui import Chart


class Config(QtWidgets.QWidget):
    """Drag sphere chart configuration"""
    TITLE = QtWidgets.QApplication.translate("pychemqt", "Drag Sphere chart")

    def __init__(self, config=None, parent=None):
        super().__init__(parent)
        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Method:")), 1, 1)
        self.metodos = QtWidgets.QComboBox()
        for f in drag.f_list:
            self.metodos.addItem(f.__name__)
        layout.addWidget(self.metodos, 1, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed), 1, 3)

        self.lineconfig = LineConfig(
            "line", QtWidgets.QApplication.translate(
                "pychemqt", "Relative roughtness style line"))
        layout.addWidget(self.lineconfig, 4, 1, 1, 3)
        self.cruxconfig = LineConfig(
            "crux", QtWidgets.QApplication.translate(
                "pychemqt", "Crux style line"))
        layout.addWidget(self.cruxconfig, 5, 1, 1, 3)

        self.gridconfig = GridConfig(
            "grid", QtWidgets.QApplication.translate(
                "pychemqt", "Grid style line"))
        layout.addWidget(self.gridconfig, 6, 1, 1, 3)

        layout.addItem(QtWidgets.QSpacerItem(
            10, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 10, 1, 1, 3)

        if config and config.has_section("drag"):
            self.metodos.setCurrentIndex(config.getint("drag", 'method'))
            self.lineconfig.setConfig(config, "drag")
            self.cruxconfig.setConfig(config, "drag")
            self.gridconfig.setConfig(config, "drag")

    def value(self, config):
        """Update ConfigParser instance with the config"""
        if not config.has_section("drag"):
            config.add_section("drag")
        config.set("drag", "method", str(self.metodos.currentIndex()))
        config = self.lineconfig.value(config, "drag")
        config = self.cruxconfig.value(config, "drag")
        config = self.gridconfig.value(config, "drag")
        return config


class ConfigDialog(QtWidgets.QDialog):
    """Dialog to configure moody chart"""
    def __init__(self, config=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle(QtWidgets.QApplication.translate(
            "pychemqt", "Drag sphere diagram configuration"))
        layout = QtWidgets.QVBoxLayout(self)
        self.widget = Config(config)
        layout.addWidget(self.widget)
        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Cancel | QtWidgets.QDialogButtonBox.StandardButton.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox)

    def value(self, config):
        """Function result for wizard"""
        config = self.widget.value(config)
        return config


class Drag(Chart):
    """Drag sphere chart dialog"""
    title = QtWidgets.QApplication.translate("pychemqt", "Drag Sphere")
    configDialog = ConfigDialog
    locLogo = (0.8, 0.85, 0.1, 0.1)
    note = None

    def config(self):
        """Initialization action in plot don't neccesary to update in plot"""

        self.plt.fig.subplots_adjust(
            left=0.1, right=0.9, bottom=0.12, top=0.98)
        self.note = None

    def click(self, event):
        """Update input and graph annotate when mouse click over chart"""
        Re = event.xdata
        Cd = event.ydata

        # Exit if click event if out of axis
        if Re is None:
            self.clearCrux()
            return

        Cd = None
        method = self.Preferences.getint("drag", "method")
        F = drag.f_list[method]
        Cd = F(Re)
        self.createCrux(Re, Cd)

    @staticmethod
    def _txt(Re, Cd):
        txt = "Re: %0.4g\nCd: %0.4g" % (Re, Cd)
        return txt

    def clearCrux(self):
        """Remove crux and note if exist"""
        self.plt.lx.set_ydata(0)
        self.plt.ly.set_xdata(0)
        if self.note:
            self.note.remove()
            self.note = None
            self.plt.draw()

    def createCrux(self, Re, Cd):
        """Create a crux in selected point of plot and show data at bottom
        right corner"""
        self.clearCrux()
        txt = self._txt(Re, Cd)

        self.plt.lx.set_ydata(Cd)
        self.plt.ly.set_xdata(Re)

        self.note = self.plt.fig.text(0.92, 0.05, txt, size="6", va="top")
        self.plt.draw()

    def plot(self):
        """Plot the drag chart using the indicate method """
        method = self.Preferences.getint("drag", "method")
        f = drag.f_list[method]

        self.plt.ax.set_autoscale_on(False)
        self.plt.ax.clear()

        Re = logspace(-1, 6, 1000)
        Cd = []
        for re in Re:
            try:
                v = f(re)
            except NotImplementedError:
                v = None
            Cd.append(v)

        kw = formatLine(self.Preferences, "drag", "line")
        self.plt.ax.plot(Re, Cd, **kw)

        kw = formatLine(self.Preferences, "drag", "crux")
        self.plt.lx = self.plt.ax.axhline(**kw)  # the horiz line
        self.plt.ly = self.plt.ax.axvline(**kw)  # the vert line

        xlabel = QtWidgets.QApplication.translate(
            "pychemqt", "Reynolds number") + ", " + r"$Re=\frac{V\rho D}{\mu}$"
        self.plt.ax.set_xlabel(xlabel, ha='center', size='10')
        self.plt.ax.set_ylabel("Drag coefficient, $C_d$, [-]", size='10')

        grid = self.Preferences.getboolean("drag", "grid")
        kw = formatLine(self.Preferences, "drag", "grid")
        del kw["marker"]
        if grid:
            self.plt.ax.grid(grid, **kw)
        else:
            self.plt.ax.grid(grid)

        self.plt.ax.set_xlim(0.1, 1e6)
        self.plt.ax.set_ylim(0.06, 300)
        self.plt.ax.set_xscale("log")
        self.plt.ax.set_yscale("log")

        self.plt.draw()

    def calculate(self):
        dlg = CalculateDialog()
        if dlg.exec():
            Re = dlg.Re.value
            Cd = dlg.Cd.value
            self.createCrux(Re, Cd)


class CalculateDialog(QtWidgets.QDialog):
    """Dialog to calculate a specified point"""
    def __init__(self, parent=None):
        super().__init__(parent)
        title = QtWidgets.QApplication.translate(
            "pychemqt", "Calculate friction factor")
        self.setWindowTitle(title)
        layout = QtWidgets.QGridLayout(self)
        label = QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Method:"))
        layout.addWidget(label, 1, 0)
        self.metodos = QtWidgets.QComboBox()
        for f in drag.f_list:
            self.metodos.addItem(f.__name__)
        self.metodos.currentIndexChanged.connect(self.calculate)
        layout.addWidget(self.metodos, 1, 1, 1, 2)

        layout.addWidget(QtWidgets.QLabel("Re"), 3, 1)
        self.Re = Entrada_con_unidades(float, tolerancia=4)
        self.Re.valueChanged.connect(self.calculate)
        layout.addWidget(self.Re, 3, 2)
        layout.addWidget(QtWidgets.QLabel("Cd"), 5, 1)
        self.Cd = Entrada_con_unidades(float, readOnly=True, decimales=8)
        layout.addWidget(self.Cd, 5, 2)

        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok | QtWidgets.QDialogButtonBox.StandardButton.Close)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox, 10, 1, 1, 2)

    def calculate(self):
        """Calculate point procedure"""
        Re = self.Re.value
        index = self.metodos.currentIndex()
        F = drag.f_list[index]
        if Re:
            Cd = F(Re)
            self.Cd.setValue(Cd)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = Drag()
    Dialog.show()
    sys.exit(app.exec())
