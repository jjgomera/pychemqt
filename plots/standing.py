#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2017, Juan José Gómez Romera <jjgomera@gmail.com>

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
# Standing Katz plot
###############################################################################


from configparser import ConfigParser
import json
import os
import re

from PyQt5 import QtGui, QtWidgets
from scipy import arange

from lib.config import conf_dir
from lib.petro import Z_list
from lib.plot import mpl
from lib.utilities import formatLine
from UI.widgets import Entrada_con_unidades


def calculate(config):
    """Plot calculate procedure"""
    method = config.getint("Standing_Katz", "method")
    Tr = eval(config.get("Standing_Katz", "Tr"))
    Prmin = config.getfloat("Standing_Katz", "Prmin")
    Prmax = config.getfloat("Standing_Katz", "Prmax")
    Z = Z_list[method]

    Pr = arange(Prmin, Prmax, 0.01)

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
        lines[t] = {"Pr": pr, "Z": z}

    dat[method] = lines

    # Save to file
    with open(conf_dir+"standing_katz.dat", "w") as file:
        json.dump(dat, file, indent=4)


class Standing_Katz(QtWidgets.QDialog):
    """Standing-Katz chart dialog"""
    title = QtWidgets.QApplication.translate(
        "pychemqt",
        "Standing and Katz compressivitity factors chart for natural gas")

    def __init__(self, parent=None):
        super(Standing_Katz, self).__init__(parent)
        self.showMaximized()
        self.setWindowTitle(self.title)
        layout = QtWidgets.QGridLayout(self)
        layout.setColumnStretch(3, 1)
        self.plt = mpl(self)
        layout.addWidget(self.plt, 2, 1, 1, 4)

        btBox = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Close)
        butonPNG = QtWidgets.QPushButton(QtGui.QIcon(
            os.environ["pychemqt"] +
            os.path.join("images", "button", "image.png")),
            QtWidgets.QApplication.translate("pychemqt", "Save as PNG"))
        butonPNG.clicked.connect(self.plt.savePNG)
        butonConfig = QtWidgets.QPushButton(QtGui.QIcon(
            os.environ["pychemqt"] +
            os.path.join("images", "button", "configure.png")),
            QtWidgets.QApplication.translate("pychemqt", "Configure"))
        butonConfig.clicked.connect(self.configure)
        butonCalculate = QtWidgets.QPushButton(QtGui.QIcon(
            os.environ["pychemqt"] +
            os.path.join("images", "button", "calculator.png")),
            QtWidgets.QApplication.translate("pychemqt", "Calculate point"))
        butonCalculate.clicked.connect(self.calculate)
        btBox.rejected.connect(self.reject)
        btBox.layout().insertWidget(0, butonPNG)
        btBox.layout().insertWidget(0, butonCalculate)
        btBox.layout().insertWidget(0, butonConfig)
        layout.addWidget(btBox, 3, 1, 1, 4)

        self.Preferences = ConfigParser()
        self.Preferences.read(conf_dir+"pychemqtrc")
        self.plot()

        self.note = None

    def configure(self):
        from UI.prefStandingKatz import Dialog
        dlg = Dialog(self.Preferences)
        if dlg.exec_():
            self.Preferences = dlg.value(self.Preferences)
            self.Preferences.write(open(conf_dir+"pychemqtrc", "w"))
            self.plot()

    def plot(self):
        """Plot the Standing-Katz chart using the indicate method """
        method = self.Preferences.getint("Standing_Katz", "method")

        Prmin = self.Preferences.getfloat("Standing_Katz", "Prmin")
        Prmax = self.Preferences.getfloat("Standing_Katz", "Prmax")

        self.plt.ax.clear()
        self.plt.ax.set_xlim(Prmin, Prmax)
        self.plt.ax.set_xlabel(r"$P_r=\frac{P}{P_c}$", ha='center', size='14')
        self.plt.ax.set_ylabel(r"$Z=\frac{PV}{nRT}$", va="bottom", size='14')
        self.plt.ax.grid(b=True, which='both', color='0.6', ls=':')

        if not os.path.isfile(conf_dir+"standing_katz.dat"):
            calculate(self.Preferences)

        load = False
        with open(conf_dir+"standing_katz.dat", "r") as file:
            try:
                dat = json.load(file)
            except ValueError:
                calculate(self.Preferences)
                load = True

            if method not in dat:
                calculate(self.Preferences)
                load = True

        # Reload file if it's created in last with statement
        if load:
            with open(conf_dir+"standing_katz.dat", "r") as file:
                dat = json.load(file)

        # Define Crux
        kw = formatLine(self.Preferences, "Standing_Katz", "crux")
        self.plt.lx = self.plt.ax.axhline(**kw)  # the horiz line
        self.plt.ly = self.plt.ax.axvline(**kw)  # the vert line

        # Plot data
        kw = formatLine(self.Preferences, "Standing_Katz", "line")
        for Tr in sorted(dat[str(method)].keys()):
            line = dat[str(method)][Tr]
            self.plt.ax.plot(line["Pr"], line["Z"], **kw)

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

        # Add explicative legend of isoline
        self.plt.ax.text(pzmin, zmin+0.1, r"$T_r$", size="12",
                         ha='left', va='center')

    def calculate(self):
        dlg = CalculateDialog()
        if dlg.exec_():
            Tr = dlg.Tr.value
            Pr = dlg.Pr.value
            Z = dlg.Z.value
            self.createCrux(Tr, Pr, Z)

    def createCrux(self, Tr, Pr, Z):
        """Create a crux in selected point of plot and show data at bottom
        right corner"""
        txt = "Tr: %0.4g\nPr: %0.4g\nZ: %0.4g" % (Tr, Pr, Z)
        self.plt.ly.set_xdata(Pr)
        self.plt.lx.set_ydata(Z)
        if self.note:
            self.note.remove()
            self.note = None
        self.note = self.plt.fig.text(0.85, 0.08, txt, size="8", va="top")
        self.plt.draw()


class CalculateDialog(QtWidgets.QDialog):
    """Dialog to calculate a specified point"""
    def __init__(self, parent=None):
        super(CalculateDialog, self).__init__(parent)
        title = QtWidgets.QApplication.translate(
            "pychemqt", "Calculate compressibility factor of natural gas")
        self.setWindowTitle(title)
        layout = QtWidgets.QGridLayout(self)
        label = QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Method:"))
        layout.addWidget(label, 1, 0)
        self.method = QtWidgets.QComboBox()
        for Z in Z_list:
            name = Z.__name__[2:].replace("_", "-")
            year = re.search("((\d+))", Z.__doc__).group(0)
            doc = "%s (%s)" % (name, year)
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

        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Close)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox, 10, 1, 1, 2)

    def calculate(self, value):
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
    sys.exit(app.exec_())
