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


from PyQt5 import QtWidgets

from scipy import logspace, log10, pi, arctan
from matplotlib.patches import ConnectionPatch

from lib.physics import f_list
from lib.plot import mpl
from lib.utilities import representacion


class Moody(QtWidgets.QDialog):
    title = QtWidgets.QApplication.translate("pychemqt", "Moody Diagram")

    def __init__(self, parent=None):
        super(Moody, self).__init__(parent)
        self.showMaximized()
        self.setWindowTitle(self.title)
        layout = QtWidgets.QGridLayout(self)
        # layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Method:")), 1, 1)
        # self.metodos = QtWidgets.QComboBox()
        # self.metodos.addItem("Colebrook")
        # self.metodos.addItem("Chen (1979")
        # self.metodos.addItem("Romeo (2002)")
        # self.metodos.addItem("Goudar-Sonnad")
        # self.metodos.addItem("Manadilli (1997)")
        # self.metodos.addItem("Serghides")
        # self.metodos.addItem("Churchill (1977)")
        # self.metodos.addItem("Zigrang-Sylvester (1982)")
        # self.metodos.addItem("Swamee-Jain (1976)")

        # self.metodos.currentIndexChanged.connect(self.cambiar)
        # layout.addWidget(self.metodos, 1, 2)
        layout.setColumnStretch(3, 1)
        self.diagrama = mpl(self)
        layout.addWidget(self.diagrama, 2, 1, 1, 4)

        self.buttonBox = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Close)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox, 3, 1, 1, 4)

        self.cambiar(0)

    def cambiar(self, int):
        self.diagrama.ax.clear()
        self.diagrama.ax.set_autoscale_on(False)
        # title = QtWidgets.QApplication.translate("pychemqt", "Moody Diagram")
        # self.diagrama.ax.set_title(title, size='12')
        xlabel = QtWidgets.QApplication.translate("pychemqt", "Reynolds number") + \
            ",  " + r"$Re=\frac{V\rho D}{\mu}$"
        self.diagrama.ax.set_xlabel(xlabel, ha='center', size='10')
        ylabel = QtWidgets.QApplication.translate("pychemqt", "Friction factor") + \
            ",  " + r"$f=\frac{2hDg}{LV^2}$"
        self.diagrama.ax.set_ylabel(ylabel, size='10')
        txt = QtWidgets.QApplication.translate("pychemqt", "Relative roughness") + \
            ", "+r"$r=\frac{\epsilon}{D}$"
        self.diagrama.fig.text(0.95, 0.5, txt, rotation=90, size='10', va = "center", ha = "center")
        self.diagrama.ax.grid(True)
        self.diagrama.ax.set_xlim(600, 1e8)
        self.diagrama.ax.set_ylim(0.008, 0.1)
        self.diagrama.ax.set_xscale("log")
        self.diagrama.ax.set_yscale("log")
        xticks = [7e2, 8e2, 9e2, 1e3, 2e3, 3e3, 4e3, 5e3, 6e3, 7e3, 8e3, 9e3,
                  1e4, 2e4, 3e4, 4e4, 5e4, 6e4, 7e4, 8e4, 9e4, 1e5, 2e5, 3e5,
                  4e5, 5e5, 6e5, 7e5, 8e5, 9e5, 1e6, 2e6, 3e6, 4e6, 5e6, 6e6,
                  7e6, 8e6, 9e6, 1e7, 2e7, 3e7, 4e7, 5e7, 6e7, 7e7, 8e7, 9e7]
        yticks = [9e-3, 1e-2, 1.1e-2, 1.2e-2, 1.3e-2, 1.4e-2, 1.5e-2, 1.6e-2,
                  1.7e-2, 1.8e-2, 1.9e-2, 2e-2, 2.1e-2, 2.2e-2, 2.3e-2, 2.4e-2,
                  2.5e-2, 2.6e-2, 2.7e-2, 2.8e-2, 2.9e-2, 3e-2, 3.2e-2, 3.4e-2,
                  3.6e-2, 3.8e-2, 4e-2, 4.2e-2, 4.4e-2, 4.6e-2, 4.8e-2, 5e-2,
                  5.2e-2, 5.4e-2, 5.6e-2, 5.8e-2, 6e-2, 6.2e-2, 6.4e-2, 6.6e-2,
                  6.8e-2, 7e-2, 7.5e-2, 8e-2, 8.5e-2, 9e-2, 9.5e-2, 1e-1]
        self.diagrama.ax.set_xticks(xticks, minor=False)
        self.diagrama.ax.set_yticks(yticks)
        self.diagrama.ax.set_yticks([], minor=True)
        ytickslabel = [9, r"$10^{-2}$", "", 1.2, "", 1.4, "", 1.6, "", 1.8, "",
                       2, "", "", "", "", 2.5, "", "", "", "", 3, "", "", "",
                       "", 4, "", "", "", "", 5, "", "", "", "", 6, "", "", "",
                       "", 7, "", 8, "", 9, "", r"$10^{-1}$"]
        self.diagrama.ax.set_yticklabels(ytickslabel)
        self.__plot(int)
        self.diagrama.draw()

    def __plot(self, metodo=0, eD=[]):
        """Plot the Moody chart using the indicate method
        método de cálculo:
            0   -   Colebrook
            1   -   Chen (1979)
            2   -   Romeo (2002)
            3   -   Goudar-Sonnad
            4   -   Manadilli (1997)
            5   -   Serghides
            6   -   Churchill (1977)
            7   -   Zigrang-Sylvester (1982)
            8   -   Swamee-Jain (1976)")

        eD: lista con las líneas de rugosidades relativas a dibujar
        Prmin: escala del eje x, minimo valor de Pr a representar
        Prmax: escala del eje y, maximo valor de Pr a representar
        """
        if not eD:
            eD = [0, 1e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 4e-4, 6e-4,
                  8e-4, 0.001, 0.0015, 0.002, 0.003, 0.004, 0.006, 0.008,
                  0.01, 0.0125, 0.015, 0.0175, 0.02, 0.025, 0.03, 0.035,
                  0.04, 0.045, 0.05, 0.06, 0.07]
        F = f_list[metodo]

        # laminar
        Re = [600, 2400]
        f = [64./R for R in Re]
        self.diagrama.ax.plot(Re, f, "k")

        # turbulent
        Re = logspace(log10(2400), 8, 50)
        for e in eD:
            self.diagrama.ax.plot(Re, [F(Rei, e) for Rei in Re], "k")
            title = representacion(e, tol=4.5)
            angle = arctan((log10(F(Re[47], e))-log10(F(Re[35], e))) /
                           (log10(Re[47])-log10(Re[35])))*360/2/pi
            self.diagrama.ax.annotate(
                title, (Re[45], F(Re[45], e)), size="small", ha="center",
                va="bottom", rotation=angle)

        # Line to define the fully desarrolled turbulent flux
        f = [(1/(1.14-2*log10(3500/R)))**2 for R in Re]
        self.diagrama.ax.plot(Re, f, "k", lw=0.5, linestyle=":")

        self.diagrama.ax.add_artist(
            ConnectionPatch((600, 0.009), (2400, 0.009), "data", "data",
                            arrowstyle="<|-|>", mutation_scale=20, fc="w"))
        self.diagrama.ax.add_artist(
            ConnectionPatch((2400, 0.009), (6000, 0.009), "data", "data",
                            arrowstyle="<|-|>", mutation_scale=20, fc="w"))
        self.diagrama.ax.add_artist(
            ConnectionPatch((6000, 0.095), (40000, 0.095), "data", "data",
                            arrowstyle="<|-|>", mutation_scale=20, fc="w"))
        self.diagrama.ax.add_artist(
            ConnectionPatch((40000, 0.095), (9.9e7, 0.095), "data", "data",
                            arrowstyle="<|-|>", mutation_scale=20, fc="w"))
        txt = QtWidgets.QApplication.translate("pychemqt", "Transition Zone")
        self.diagrama.ax.text(15000, 0.094, txt, size="small", va="top", ha="center")
        txt = QtWidgets.QApplication.translate("pychemqt", "Turbulent flux fully desarrolled")
        self.diagrama.ax.text(2e6, 0.094, txt, size="small", va="top", ha="center")
        txt = QtWidgets.QApplication.translate("pychemqt", "Critic\nzone")
        self.diagrama.ax.text(4000, 0.0091, txt, size="small", va="bottom", ha="center")
        txt = QtWidgets.QApplication.translate("pychemqt", "Laminar flux")
        self.diagrama.ax.text(1200, 0.0091, txt, size="small", va="bottom", ha="center")

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = Moody()
    Dialog.show()
    sys.exit(app.exec_())

