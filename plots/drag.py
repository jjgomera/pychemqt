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


from PyQt5 import QtWidgets
from numpy import logspace

from lib import drag
from lib.utilities import formatLine
from UI.widgets import Entrada_con_unidades

from plots.ui import Chart
from plots.moody import ConfigDialog


class Drag(Chart):
    """Moody chart dialog"""
    title = QtWidgets.QApplication.translate("pychemqt", "Drag Sphere")
    configDialog = ConfigDialog
    locLogo = (0.8, 0.85, 0.1, 0.1)

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
        # method = self.Preferences.getint("Moody", "method")
        method = 0
        F = drag._all[method]
        Cd = F(Re)
        self.createCrux(Re, Cd)

    def _txt(self, Re, Cd):
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

        self.note = self.plt.fig.text(0.85, 0.08, txt, size="6", va="top")
        self.plt.draw()

    def plot(self):
        """Plot the drag chart using the indicate method """
        # method = self.Preferences.getint("Moody", "method")
        method = 0
        f = drag._all[method]

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

        kw = formatLine(self.Preferences, "Moody", "line")
        self.plt.ax.plot(Re, Cd, **kw)

        kw = formatLine(self.Preferences, "Moody", "crux")
        self.plt.lx = self.plt.ax.axhline(**kw)  # the horiz line
        self.plt.ly = self.plt.ax.axvline(**kw)  # the vert line

        xlabel = QtWidgets.QApplication.translate(
            "pychemqt", "Reynolds number") + ", " + r"$Re=\frac{V\rho D}{\mu}$"
        self.plt.ax.set_xlabel(xlabel, ha='center', size='10')
        self.plt.ax.set_ylabel("Drag coefficient, $C_d$, [-]", size='10')
        self.plt.ax.grid(True, "both", ls=":", lw=0.8)
        self.plt.ax.set_xlim(0.1, 1e6)
        self.plt.ax.set_ylim(0.06, 300)
        self.plt.ax.set_xscale("log")
        self.plt.ax.set_yscale("log")

        self.plt.draw()

    def calculate(self):
        dlg = CalculateDialog()
        if dlg.exec_():
            Re = dlg.Re.value
            Cd = dlg.Cd.value
            self.createCrux(Re, Cd)


class CalculateDialog(QtWidgets.QDialog):
    """Dialog to calculate a specified point"""
    def __init__(self, parent=None):
        super(CalculateDialog, self).__init__(parent)
        title = QtWidgets.QApplication.translate(
            "pychemqt", "Calculate friction factor")
        self.setWindowTitle(title)
        layout = QtWidgets.QGridLayout(self)
        label = QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Method:"))
        layout.addWidget(label, 1, 0)
        self.metodos = QtWidgets.QComboBox()
        for f in drag._all:
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
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Close)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox, 10, 1, 1, 2)

    def calculate(self, value):
        index = self.metodos.currentIndex()
        F = drag._all[index]
        Re = self.Re.value
        if Re:
            Cd = F(Re)
            self.Cd.setValue(Cd)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = Drag()
    Dialog.show()
    sys.exit(app.exec_())
