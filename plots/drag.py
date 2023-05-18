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


.. include:: drag.rst


The module include all related moody chart functionality
    * :class:`Drag`: Chart dialog
    * :class:`CalculateDialog`: Dialog to calculate a specified point

and its configuration

    * :class:`Config`: Drag sphere chart configuration
'''


from numpy import logspace
from tools.qt import QtWidgets, tr

from lib import drag
from lib.config import Preferences
from lib.utilities import formatLine
from UI.widgets import Entrada_con_unidades, GridConfig, LineConfig

from plots.ui import Chart


class Config(QtWidgets.QWidget):
    """Drag sphere chart configuration"""
    TITLE = tr("pychemqt", "Drag Sphere chart")
    TITLECONFIG = tr("pychemqt", "Drag sphere diagram configuration")

    def __init__(self, config=None, parent=None):
        super().__init__(parent)
        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(QtWidgets.QLabel(tr("pychemqt", "Method:")), 1, 1)
        self.metodos = QtWidgets.QComboBox()
        for f in drag.f_list:
            self.metodos.addItem(f.__name__)
        layout.addWidget(self.metodos, 1, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed), 1, 3)

        self.lineconfig = LineConfig(
            "line", tr("pychemqt", "Drag coefficient style line"))
        layout.addWidget(self.lineconfig, 4, 1, 1, 3)
        self.cruxconfig = LineConfig(
            "crux", tr("pychemqt", "Crux style line"))
        layout.addWidget(self.cruxconfig, 5, 1, 1, 3)

        self.gridconfig = GridConfig(
            "grid", tr("pychemqt", "Grid style line"))
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


class Drag(Chart):
    """Drag sphere chart dialog"""
    title = tr("pychemqt", "Drag Sphere")
    widgetConfig = Config
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
        method = Preferences.getint("drag", "method")
        F = drag.f_list[method]
        Cd = F(Re)
        self.createCrux(Re, Cd)

    @staticmethod
    def _txt(Re, Cd):
        txt = f"Re: {Re:0.4g}\n$C_d$: {Cd:0.4g}"
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

        self.note = self.plt.fig.text(0.92, 0.05, txt, size="8", va="top")
        self.plt.draw()

    def plot(self):
        """Plot the drag chart using the indicate method """
        method = Preferences.getint("drag", "method")
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

        kw = formatLine(Preferences, "drag", "line")
        self.plt.ax.plot(Re, Cd, **kw)

        kw = formatLine(Preferences, "drag", "crux")
        self.plt.lx = self.plt.ax.axhline(**kw)  # the horiz line
        self.plt.ly = self.plt.ax.axvline(**kw)  # the vert line

        xlabel = tr("pychemqt", "Reynolds number") + ", " + \
            r"$Re=\frac{V\rho D}{\mu}$"
        self.plt.ax.set_xlabel(xlabel, ha='center', size='10')
        self.plt.ax.set_ylabel("Drag coefficient, $C_d$, [-]", size='10')

        grid = Preferences.getboolean("drag", "grid")
        kw = formatLine(Preferences, "drag", "grid")
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
        title = tr("pychemqt", "Calculate friction factor")
        self.setWindowTitle(title)
        layout = QtWidgets.QGridLayout(self)
        label = QtWidgets.QLabel(tr("pychemqt", "Method:"))
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
