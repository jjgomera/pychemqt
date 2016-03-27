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


###############################################################################
# Moody plot
###############################################################################


import os
import json

from PyQt5 import QtWidgets
from scipy import logspace, log10
from matplotlib.patches import ConnectionPatch

from lib.config import conf_dir, Preferences
from lib.friction import f_list
from lib.plot import mpl
from lib.utilities import representacion
from UI.widgets import Entrada_con_unidades


Re_laminar = [600, 2400]
Re_turbulent = logspace(log10(2400), 8, 50)
Re_fully = logspace(log10(4000), 8, 50)


def calculate(eD, F):
    fanning = Preferences.getboolean("Moody", "fanning")
    method = Preferences.getint("Moody", "method")

    dat = {}
    dat["fanning"] = fanning
    dat["method"] = method
    # laminar
    if fanning:
        dat["laminar"] = [16./R for R in Re_laminar]
        x = 4
    else:
        dat["laminar"] = [64./R for R in Re_laminar]
        x = 1

    # turbulent
    turb = {}
    for e in eD:
        turb[e] = [F(Rei, e)/x for Rei in Re_turbulent]
        dat["turbulent"] = turb

    # Line to define the fully desarrolled turbulent flux
    # dat["fully"] = [(1/(1.14-2*log10(3500/R)))**2 for R in Re_fully]

    # Save to file
    with open(conf_dir+"moody.dat", "w") as file:
        json.dump(dat, file, indent=4)


class Moody(QtWidgets.QDialog):
    """Moody chart dialog"""
    title = QtWidgets.QApplication.translate("pychemqt", "Moody Diagram")

    def __init__(self, parent=None):
        super(Moody, self).__init__(parent)
        self.showMaximized()
        self.setWindowTitle(self.title)
        layout = QtWidgets.QGridLayout(self)
        layout.setColumnStretch(3, 1)
        self.diagrama = mpl(self)
        layout.addWidget(self.diagrama, 2, 1, 1, 4)

        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Close)
        self.buttonBox.rejected.connect(self.reject)
        txt = QtWidgets.QApplication.translate("pychemqt", "Calculate point")
        self.buttonCalculation = QtWidgets.QPushButton(txt)
        self.buttonBox.addButton(self.buttonCalculation,
                                 QtWidgets.QDialogButtonBox.ActionRole)
        self.buttonCalculation.clicked.connect(self.calculate)
        layout.addWidget(self.buttonBox, 3, 1, 1, 4)

        self.__plot()

    def __plot(self):
        """Plot the Moody chart using the indicate method """
        fanning = Preferences.getboolean("Moody", "fanning")
        method = Preferences.getint("Moody", "method")

        if fanning:
            x = 4
        else:
            x = 1

        self.diagrama.ax.clear()
        self.diagrama.ax.set_autoscale_on(False)
        xlabel = QtWidgets.QApplication.translate(
            "pychemqt", "Reynolds number") + ", " + r"$Re=\frac{V\rho D}{\mu}$"
        self.diagrama.ax.set_xlabel(xlabel, ha='center', size='10')
        if fanning:
            ylabel = QtWidgets.QApplication.translate(
                "pychemqt", "Fanning Friction factor")
            formula = r"$f_f=\frac{2hDg}{LV^2}$"
        else:
            ylabel = QtWidgets.QApplication.translate(
                "pychemqt", "Darcy Friction factor")
            formula = r"$f_d=\frac{2hDg}{LV^2}$"
        self.diagrama.ax.set_ylabel(ylabel+",  " + formula, size='10')
        txt = QtWidgets.QApplication.translate(
            "pychemqt", "Relative roughness") + ", "+r"$r=\frac{\epsilon}{D}$"
        self.diagrama.fig.text(0.97, 0.5, txt, rotation=90, size='10',
                               va="center", ha="center")
        self.diagrama.ax.grid(True)
        self.diagrama.ax.set_xlim(600, 1e8)
        self.diagrama.ax.set_ylim(0.008/x, 0.11/x)
        self.diagrama.ax.set_xscale("log")
        self.diagrama.ax.set_yscale("log")
        xticks = [7e2, 8e2, 9e2, 1e3, 2e3, 3e3, 4e3, 5e3, 6e3, 7e3, 8e3, 9e3,
                  1e4, 2e4, 3e4, 4e4, 5e4, 6e4, 7e4, 8e4, 9e4, 1e5, 2e5, 3e5,
                  4e5, 5e5, 6e5, 7e5, 8e5, 9e5, 1e6, 2e6, 3e6, 4e6, 5e6, 6e6,
                  7e6, 8e6, 9e6, 1e7, 2e7, 3e7, 4e7, 5e7, 6e7, 7e7, 8e7, 9e7]
        self.diagrama.ax.set_xticks(xticks, minor=False)
        self.diagrama.ax.set_yticks([], minor=True)
        if fanning:
            yticks = [2e-3, 2.5e-3, 3e-3, 3.5e-3, 4e-3, 5e-3, 6e-3, 7e-3, 8e-3,
                      9e-3, 1e-2, 1.1e-2, 1.2e-2, 1.3e-2, 1.4e-2, 1.5e-2,
                      1.6e-2, 1.7e-2, 1.8e-2, 1.9e-2, 2e-2, 2.1e-2, 2.2e-2,
                      2.3e-2, 2.4e-2, 2.5e-2]
            ytickslabel = [2, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 9, r"$10^{-2}$", "",
                           1.2, "", 1.4, "", 1.6, "", 1.8, "", 2, "", "", "",
                           "", 2.5]
        else:
            yticks = [9e-3, 1e-2, 1.1e-2, 1.2e-2, 1.3e-2, 1.4e-2, 1.5e-2,
                      1.6e-2, 1.7e-2, 1.8e-2, 1.9e-2, 2e-2, 2.1e-2, 2.2e-2,
                      2.3e-2, 2.4e-2, 2.5e-2, 2.6e-2, 2.7e-2, 2.8e-2, 2.9e-2,
                      3e-2, 3.2e-2, 3.4e-2, 3.6e-2, 3.8e-2, 4e-2, 4.2e-2,
                      4.4e-2, 4.6e-2, 4.8e-2, 5e-2, 5.2e-2, 5.4e-2, 5.6e-2,
                      5.8e-2, 6e-2, 6.2e-2, 6.4e-2, 6.6e-2, 6.8e-2, 7e-2,
                      7.5e-2, 8e-2, 8.5e-2, 9e-2, 9.5e-2, 1e-1]
            ytickslabel = [9, r"$10^{-2}$", "", 1.2, "", 1.4, "", 1.6, "", 1.8,
                           "", 2, "", "", "", "", 2.5, "", "", "", "", 3, "",
                           "", "", "", 4, "", "", "", "", 5, "", "", "", "", 6,
                           "", "", "", "", 7, "", 8, "", 9, "", r"$10^{-1}$"]
        self.diagrama.ax.set_yticks(yticks)
        self.diagrama.ax.set_yticklabels(ytickslabel)

        eD = [0, 1e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 4e-4, 6e-4, 8e-4,
              .001, .0015, .002, .003, .004, .006, .008, .01, .0125, .015,
              .0175, .02, .025, .03, .035, .04, .045, .05, .06, .07]
        F = f_list[method]

        if not os.path.isfile(conf_dir+"moody.dat"):
            calculate(eD, F)

        load = False
        with open(conf_dir+"moody.dat", "r") as file:
            try:
                dat = json.load(file)
            except ValueError:
                calculate(eD, F)
                load = True

            if dat["fanning"] != fanning or dat["method"] != method:
                calculate(eD, F)
                load = True

        # Reload file if is create in las with statement
        if load:
            with open(conf_dir+"moody.dat", "r") as file:
                dat = json.load(file)

        # Plot data
        self.diagrama.ax.plot(Re_laminar, dat["laminar"], "k")
        for eD, f in dat["turbulent"].items():
            self.diagrama.ax.plot(Re_turbulent, f, "k")
            title = " " + representacion(eD, tol=4.5)
            if f[-1] > 0.008/x:
                self.diagrama.ax.text(1e8, f[-1], title, size="x-small",
                                      ha='left', va='center')
        # self.diagrama.ax.plot(Re_fully, dat["fully"], "k", lw=0.5, ls=":")

        # Add explicative legend
        # Laminar zone
        self.diagrama.ax.add_artist(
            ConnectionPatch((600, 0.009/x), (2400, 0.009/x), "data", "data",
                            arrowstyle="<|-|>", mutation_scale=20, fc="w"))
        txt = QtWidgets.QApplication.translate("pychemqt", "Laminar flux")
        self.diagrama.ax.text(1200, 0.0091/x, txt, size="small", va="bottom",
                              ha="center")

        # Critic zone
        self.diagrama.ax.add_artist(
            ConnectionPatch((2400, 0.009/x), (4000, 0.009/x), "data", "data",
                            arrowstyle="<|-|>", mutation_scale=20, fc="w"))
        txt = QtWidgets.QApplication.translate("pychemqt", "Critic\nzone")
        self.diagrama.ax.text(3200, 0.0091/x, txt, size="small", va="bottom",
                              ha="center")

        # Transition zone
        self.diagrama.ax.add_artist(
            ConnectionPatch((4000, 0.095/x), (10000, 0.095/x), "data", "data",
                            arrowstyle="<|-|>", mutation_scale=20, fc="w"))
        txt = QtWidgets.QApplication.translate("pychemqt", "Transition Zone")
        self.diagrama.ax.text(7000, 0.104/x, txt, size="small", va="top",
                              ha="center")

        # Turbulent zone
        self.diagrama.ax.add_artist(
            ConnectionPatch((10000, 0.095/x), (9.9e7, 0.095/x), "data", "data",
                            arrowstyle="<|-|>", mutation_scale=20, fc="w"))
        txt = QtWidgets.QApplication.translate(
            "pychemqt", "Turbulent flux fully developed")
        self.diagrama.ax.text(1e6, 0.104/x, txt, size="small", va="top",
                              ha="center")

        # Smooth tubes
        self.diagrama.ax.add_artist(
            ConnectionPatch((1e6, 0.009/x), (1.5e6, 0.011/x), "data", "data",
                            arrowstyle="<|-|>", mutation_scale=20, fc="w"))
        txt = QtWidgets.QApplication.translate("pychemqt", "Smooth tubes")
        self.diagrama.ax.text(1e6, 0.009/x, txt, size="small", va="top",
                              ha="right")

        # Laminar equation
        if fanning:
            txt = r"$f=\frac{16}{Re}$"
        else:
            txt = r"$f=\frac{64}{Re}$"
        self.diagrama.ax.text(1.4e3, 0.042/x, txt, size="12",
                              va="top", ha="center", rotation=-66)
        # TODO: Calculate the angle dynamically

    def calculate(self):
        dialog = CalculateDialog()
        dialog.exec_()


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
        for f in f_list:
            line = f.__doc__.split("\n")[1]
            year = line.split(" ")[-1]
            name = line.split(" ")[-3]
            doc = name + " " + year
            self.metodos.addItem(doc)
        self.metodos.currentIndexChanged.connect(self.calculate)
        layout.addWidget(self.metodos, 1, 1, 1, 2)
        self.fanning = QtWidgets.QCheckBox(QtWidgets.QApplication.translate(
            "pychemqt", "Calculate fanning friction factor"))
        self.fanning.toggled.connect(self.calculate)
        layout.addWidget(self.fanning, 2, 0, 1, 3)

        layout.addWidget(QtWidgets.QLabel("Re"), 3, 1)
        self.Re = Entrada_con_unidades(float, tolerancia=4)
        self.Re.valueChanged.connect(self.calculate)
        layout.addWidget(self.Re, 3, 2)
        layout.addWidget(QtWidgets.QLabel("e/D"), 4, 1)
        self.eD = Entrada_con_unidades(float)
        self.eD.valueChanged.connect(self.calculate)
        layout.addWidget(self.eD, 4, 2)
        layout.addWidget(QtWidgets.QLabel("f"), 5, 1)
        self.f = Entrada_con_unidades(float, readOnly=True, decimales=8)
        layout.addWidget(self.f, 5, 2)

        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Close)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox, 10, 1, 1, 2)

    def calculate(self, value):
        index = self.metodos.currentIndex()
        F = f_list[index]
        Re = self.Re.value
        eD = self.eD.value
        if Re and eD is not None:
            f = F(Re, eD)
            if self.fanning.isChecked():
                f /= 4
            self.f.setValue(f)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = Moody()
    Dialog.show()
    sys.exit(app.exec_())
