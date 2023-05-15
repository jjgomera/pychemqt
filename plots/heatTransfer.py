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
# Heat Transfer chart
###############################################################################

import os
from math import pi

from matplotlib import image
from numpy import arange, logspace, arctan
from tools.qt import QtWidgets, tr

from lib.config import IMAGE_PATH
from lib.heatTransfer import (efectividad, TemperatureEffectiveness,
                              CorrectionFactor, Fi)

from plots.ui import Chart


class ChartHeat(Chart):
    """Dialog to implement general heat exchanger plot"""
    flujo = []
    MIX = ("Cmin", "Cmax")
    locImage = None
    locLogo = None

    def customUI(self):
        """Define custom UI element for heat transfer plot"""
        self.butonConf.setEnabled(False)
        self.butonCalc.setEnabled(False)

        lyt = QtWidgets.QHBoxLayout()
        lyt.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Flow Arrangement")))
        self.flow = QtWidgets.QComboBox()
        for text in self.flujo:
            self.flow.addItem(text[0])
        self.flow.currentIndexChanged.connect(self.plot)
        lyt.addWidget(self.flow)

        self.mixed = QtWidgets.QComboBox()
        for text in self.MIX:
            self.mixed.addItem(text)
        self.mixed.currentIndexChanged.connect(self.plot)
        lyt.addWidget(self.mixed)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed))
        self.layout().addLayout(lyt, 0, 1, 1, 3)

        self.image = self.plt.fig.add_axes(self.locImage, anchor='SW')

    def set_image(self, flux):
        """Draw flux image in plot"""
        self.image.clear()
        logo = image.imread(os.path.join(
            IMAGE_PATH, "equipment", "heat transfer", flux)+".png")
        self.image.imshow(logo, alpha=0.9)
        self.image.axis('off')


class Efectividad(ChartHeat):
    """Heat Exchanger effectiveness plot"""
    title = tr(
        "pychemqt", "Heat Exchanger effectiveness")
    flujo = [
        (tr("pychemqt", "Counterflow"), "CF"),
        (tr("pychemqt", "Parallelflow"), "PF"),
        (tr(
            "pychemqt", "Crossflow, both fluids unmixed"), "CrFunMix"),
        (tr(
            "pychemqt", "Crossflow, one fluid mixed"), "CrFSMix"),
        (tr(
            "pychemqt", "Crossflow, both fluids mixed"), "CrFMix"),
        (tr(
            "pychemqt", "1-2 pass shell and tube exchanger"), "1-2TEMAE")]

    locImage = (0.7, 0.15, 0.2, 0.2)
    locLogo = (0.13, 0.77, 0.1, 0.1)

    def plot(self):
        self.plt.ax.clear()
        self.plt.ax.set_xlim(0, 6)
        self.plt.ax.set_ylim(0, 1)
        title = tr(
            "pychemqt", "Heat Transfer effectiveness")
        self.plt.ax.set_title(title, size='12')
        self.plt.ax.set_xlabel("NTU", size='12')
        self.plt.ax.set_ylabel("ε", size='14')

        index = self.flow.currentIndex()
        flujo = self.flujo[index][1]
        self.set_image(flujo)
        self.mixed.setVisible(flujo == "CrFSMix")
        kw = {}
        if flujo == "CrFSMix":
            kw["mixed"] = str(self.mixed.currentText())

        C = [0, 0.2, 0.4, 0.6, 0.8, 1.]

        NTU = arange(0, 6.1, 0.1)
        for ci in C:
            e = [0]
            for N in NTU[1:]:
                e.append(efectividad(N, ci, flujo, **kw))
            self.plt.plot(NTU, e, "k")

            fraccionx = (NTU[40]-NTU[30])/6
            fracciony = (e[40]-e[30])
            try:
                angle = arctan(fracciony/fraccionx)*360/2/pi
            except ZeroDivisionError:
                angle = 90

            self.plt.ax.annotate(
                "C*=%0.1f" % ci, (NTU[29], e[30]), rotation=angle,
                size="medium", ha="left", va="bottom")

        self.plt.draw()


class TemperatureEfectividad(ChartHeat):
    """Heat exchanger temperature-effectiveness plot"""

    title = tr(
        "pychemqt", "Heat Exchanger temperature effectiveness")
    flujo = [
        (tr("pychemqt", "Counterflow"), "CF"),
        (tr("pychemqt", "Parallelflow"), "PF"),
        (tr(
            "pychemqt", "Crossflow, both fluids unmixed"), "CrFunMix"),
        (tr(
            "pychemqt", "Crossflow, one fluid mixed"), "CrFSMix"),
        (tr(
            "pychemqt", "Crossflow, both fluids mixed"), "CrFMix"),
        (tr(
            "pychemqt", "1-2 TEMA E"), "1-2TEMAE"),
        (tr(
            "pychemqt", "1-2 TEMA E, shell fluid divided"), "1-2TEMAE2"),
        (tr(
            "pychemqt", "1-3 TEMA E"), "1-3TEMAE"),
        (tr(
            "pychemqt", "1-4 TEMA E"), "1-4TEMAE"),
        (tr(
            "pychemqt", "1-1 TEMA G"), "1-1TEMAG"),
        (tr(
            "pychemqt", "1-2 TEMA G"), "1-2TEMAG"),
        (tr(
            "pychemqt", "1-1 TEMA H"), "1-1TEMAH"),
        (tr(
            "pychemqt", "1-2 TEMA H"), "1-2TEMAH"),
        (tr(
            "pychemqt", "1-1 TEMA J"), "1-1TEMAJ"),
        (tr(
            "pychemqt", "1-2 TEMA J"), "1-2TEMAJ"),
        (tr(
            "pychemqt", "1-4 TEMA J"), "1-4TEMAJ")]

    locImage = (0.13, 0.6, 0.15, 0.15)
    locLogo = (0.13, 0.77, 0.1, 0.1)

    def plot(self):
        self.plt.ax.clear()
        self.plt.ax.set_xlim(0.1, 10)
        self.plt.ax.set_ylim(0, 1)
        self.plt.ax.set_xscale("log")

        self.plt.ax.set_title(tr(
            "pychemqt", "Heat Transfer Temperature Effectiveness"), size='12')
        self.plt.ax.set_xlabel("NTU", size='12')
        self.plt.ax.set_ylabel("P", size='14')
        # self.plt.ax.set_xticklabels([0.1, 1.0, 10])
        # xticklabels = [0.2, 0.3, None, 0.5, None, 0.7, None, None, 2.0, 3.0,
        #                None, 5.0, None, 7.0, None, None]
        # self.plt.ax.set_xticklabels(xticklabels, minor=True)

        index = self.flow.currentIndex()
        flujo = self.flujo[index][1]
        self.set_image(flujo)
        self.mixed.setVisible(flujo == "CrFSMix")
        kwargs = {}
        if flujo == "CrFSMix":
            kwargs["mixed"] = str(self.mixed.currentText())

        R = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6, 1.8,
             2., 2.5, 3., 4., 6., 8., 10., 15.]

        NTU = logspace(-1.5, 1, 100)
        for ri in R:
            e = [0]
            for N in NTU[1:]:
                e.append(TemperatureEffectiveness(N, ri, flujo, **kwargs))
            self.plt.plot(NTU, e, "k")
            self.plt.ax.annotate(" R=%0.1f" % ri, (NTU[-1], e[-1]),
                                 size="medium", ha="left", va="center")

        self.plt.draw()


class F(ChartHeat):
    """Heat Exchanger correction factor plot"""
    title = tr(
        "pychemqt", "ΔT Correction Factor")

    flujo = [
        (tr(
            "pychemqt", "Crossflow, both fluids unmixed"), "CrFunMix"),
        (tr(
            "pychemqt", "Crossflow, one fluid mixed"), "CrFSMix"),
        (tr(
            "pychemqt", "Crossflow, both fluids mixed"), "CrFMix"),
        (tr(
            "pychemqt", "1-2 pass shell and tube exchanger"), "1-2TEMAE")]

    locImage = (0.13, 0.6, 0.15, 0.15)
    locLogo = (0.13, 0.77, 0.1, 0.1)

    def plot(self):
        self.plt.ax.clear()
        self.plt.ax.set_xlim(0, 1)
        self.plt.ax.set_ylim(0, 1)
        title = r"$\Delta T_{ml}$ " + \
            tr("pychemqt", "Correction Factor")
        self.plt.ax.set_title(title, size='12')
        xlabel = "$P=\\frac{T_{1o}-T_{1i}}{T_{2i}-T_{1i}}$"
        self.plt.ax.set_xlabel(xlabel, size='12')
        self.plt.ax.set_ylabel("F", size='14')

        index = self.flow.currentIndex()
        flujo = self.flujo[index][1]
        self.set_image(flujo)
        self.mixed.setVisible(flujo == "CrFSMix")
        kwargs = {}
        if flujo == "CrFSMix":
            kwargs["mixed"] = str(self.mixed.currentText())

        R = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6,
             1.8, 2, 2.5, 3, 4, 6, 8, 10, 15, 20]

        P = arange(0, 1.01, 0.01)
        for ri in R:
            f = [CorrectionFactor(p, ri, flujo, **kwargs) for p in P]
            self.plt.plot(P, f, "k")

            # fraccionx=P[90]-P[80]
            # fracciony=f[90]-f[80]
            # try:
                # angle=arctan(fracciony/fraccionx)*360/2/pi
            # except ZeroDivisionError:
                # angle=90
            # self.plt.ax.annotate(
                # "R=%0.1f" %ri, (P[90], f[90]), rotation=angle, size="medium",
                # horizontalalignment="left", verticalalignment="bottom")

        self.plt.draw()


class Phi(ChartHeat):
    """Heat Exchanger correction factor plot"""
    title = tr("pychemqt", "ψ", None)

    flujo = [
        # (tr("pychemqt", "Crossflow, both fluids unmixed"), "CrFunMix"),
        # (tr("pychemqt", "Crossflow, one fluid mixed"), "CrFSMix"),
        # (tr("pychemqt", "Crossflow, both fluids mixed"), "CrFMix"),
        (tr("pychemqt", "1-2 pass shell and tube exchanger"), "1-2TEMAE")]

    locImage = (0.75, 0.6, 0.15, 0.15)
    locLogo = (0.8, 0.77, 0.1, 0.1)

    def plot(self):
        self.plt.ax.clear()
        self.plt.ax.set_xlim(0, 1)
        self.plt.ax.set_ylim(0, 1)
        title = r"$\Delta T_{ml}$ " + \
            tr("pychemqt", " Correction Factor")
        self.plt.ax.set_title(title, size='12')
        self.plt.ax.set_xlabel(
            "$P=\\frac{T_{1o}-T_{1i}}{T_{2i}-T_{1i}}$", size='12')
        self.plt.ax.set_ylabel("F", size='14')

        index = self.flow.currentIndex()
        flujo = self.flujo[index][1]
        self.set_image(flujo)
        self.mixed.setVisible(flujo == "CrFSMix")
        kwargs = {}
        if flujo == "CrFSMix":
            kwargs["mixed"] = str(self.mixed.currentText())

        R = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.4,
             1.6, 1.8, 2, 2.5, 3, 4, 6, 8, 10]

        P = arange(0, 1.01, 0.01)
        for ri in R:
            f = [Fi(p, ri, flujo, **kwargs) for p in P]
            self.plt.plot(P, f, "k")

        NTU = [0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.]
        for ntu in NTU:
            self.plt.plot([0, 1], [0, 1./ntu], "k", linestyle=":")

            # fraccionx=P[90]-P[80]
            # fracciony=f[90]-f[80]
            # try:
                # angle=arctan(fracciony/fraccionx)*360/2/pi
            # except ZeroDivisionError:
                # angle=90
            # self.plt.ax.annotate(
                # "R=%0.1f" %ri, (P[90], f[90]), rotation=angle, size="medium",
                # horizontalalignment="left", verticalalignment="bottom")

        self.plt.draw()


chartHE = (Efectividad, TemperatureEfectividad, F, Phi)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    dialogo = Efectividad()
    dialogo.show()
    sys.exit(app.exec())
