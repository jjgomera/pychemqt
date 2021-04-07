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
# Heat Transfer chart
###############################################################################

import os

from matplotlib import image
from PyQt5 import QtGui, QtWidgets

from scipy import arctan, pi
from numpy import arange, logspace

from lib.plot import mpl
from lib.heatTransfer import (efectividad, TemperatureEffectiveness,
                              CorrectionFactor, Fi, NTU_fPR)


class Chart(QtWidgets.QDialog):
    """Dialog to implement general heat exchanger plot"""
    PosImage = None

    def __init__(self, parent=None):
        super(Chart, self).__init__(parent)
        self.setWindowTitle(self.title)
        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Flow Arrangement")), 1, 1)
        self.flow = QtWidgets.QComboBox()
        for text, flow in self.flujo:
            self.flow.addItem(text)
        self.flow.currentIndexChanged.connect(self.plot)
        layout.addWidget(self.flow, 1, 2)
        self.mixed = QtWidgets.QComboBox()
        for text in self.mezclado:
            self.mixed.addItem(text)
        self.mixed.currentIndexChanged.connect(self.changeMixed)
        layout.addWidget(self.mixed, 1, 3)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Fixed), 1, 4)

        self.diagrama = mpl(self, width=10, height=10)
        self.image = self.diagrama.fig.figimage([[0]], 0, 0, zorder=1)
        logo = image.imread('images/pychemqt_98.png')
        self.logo = self.diagrama.fig.figimage(logo, 0, 0, zorder=1)

        self.botonGuardar = QtWidgets.QToolButton()
        icon = os.environ["pychemqt"]+"/images/button/fileSave.png"
        self.botonGuardar.setIcon(QtGui.QIcon(QtGui.QPixmap(icon)))
        self.botonGuardar.setToolTip(
            QtWidgets.QApplication.translate("pychemqt", "Save chart to file"))
        self.botonGuardar.clicked.connect(self.save)
        layout.addWidget(self.botonGuardar, 5, 1)
        layout.addWidget(self.diagrama, 2, 1, 1, 6)
        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Close)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox, 5, 2, 1, 5)

        self.plot(0)
        self.showMaximized()
        self.refixImage()

    def resizeEvent(self, event):
        """Implement resizeEvent to precalculate the new position of image"""
        self.refixImage()
        QtWidgets.QDialog.resizeEvent(self, event)

    def changeMixed(self, indice):
        modelo = self.flow.currentIndex()
        self.plot(modelo)

    def refixImage(self, event=None):
        """Recalculate image position as resize/move window"""
        xmin, xmax = self.diagrama.ax.get_xlim()
        ymin, ymax = self.diagrama.ax.get_ylim()
        a = self.PosLogo[0]*xmax+(1-self.PosLogo[0])*xmin
        b = self.PosLogo[1]*ymax+(1-self.PosLogo[1])*ymin
        x, y = self.diagrama.ax.transData.transform_point((a, b))
        hlogo, wlogo = self.logo.get_size()
        himage, wimage = self.image.get_size()
        if self.PosLogo[0] == 1:
            x -= wlogo
            xi = x-wimage-20
        else:
            xi = x+wlogo+20
        if self.PosLogo[1] == 1:
            yi = y-himage-10
            y -= hlogo
        else:
            yi = y+10

        self.logo.ox = x
        self.logo.oy = y

        if self.PosImage:
            a = self.PosImage[0]*xmax+(1-self.PosImage[0])*xmin
            b = self.PosImage[1]*ymax+(1-self.PosImage[1])*ymin
            xi, yi = self.diagrama.ax.transData.transform_point((a, b))
            if self.PosImage[0] == 1:
                xi -= wimage+10
            else:
                xi += 10
            if self.PosImage[1] == 1:
                yi -= himage+10
            else:
                yi += 10

        self.image.ox = xi
        self.image.oy = yi
        self.diagrama.draw()

    def save(self):
        """Show the dialog to select name to save the image as file"""
        fname = str(QtWidgets[0].QFileDialog.getSaveFileName(
            self, QtWidgets.QApplication.translate("pychemqt", "Save chart to file"),
            "./", "Portable Network Graphics (*.png)"))
        self.diagrama.fig.savefig(fname, facecolor='#eeeeee')


class Efectividad(Chart):
    """Heat Exchanger effectiveness plot"""
    title = QtWidgets.QApplication.translate(
        "pychemqt", "Heat Exchanger effectiveness")
    flujo = [
        (QtWidgets.QApplication.translate("pychemqt", "Counterflow"), "CF"),
        (QtWidgets.QApplication.translate("pychemqt", "Parallelflow"), "PF"),
        (QtWidgets.QApplication.translate(
            "pychemqt", "Crossflow, both fluids unmixed"), "CrFunMix"),
        (QtWidgets.QApplication.translate(
            "pychemqt", "Crossflow, one fluid mixed, other unmixed"), "CrFSMix"),
        (QtWidgets.QApplication.translate(
            "pychemqt", "Crossflow, both fluids mixed"), "CrFMix"),
        (QtWidgets.QApplication.translate(
            "pychemqt", "1-2 pass shell and tube exchanger"), "1-2TEMAE")]

    mezclado = ("Cmin", "Cmax")
    PosLogo = (0, 1)
    PosImage = (1, 0)

    def plot(self, indice):
        self.diagrama.ax.clear()
        self.diagrama.ax.set_xlim(0, 6)
        self.diagrama.ax.set_ylim(0, 1)
        title = QtWidgets.QApplication.translate(
            "pychemqt", "Heat Transfer effectiveness")
        self.diagrama.ax.set_title(title, size='12')
        self.diagrama.ax.set_xlabel("NTU", size='12')
        self.diagrama.ax.set_ylabel("ε", size='14')

        flujo = self.flujo[indice][1]
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
            self.diagrama.plot(NTU, e, "k")

            fraccionx = (NTU[40]-NTU[30])/6
            fracciony = (e[40]-e[30])
            try:
                angle = arctan(fracciony/fraccionx)*360/2/pi
            except ZeroDivisionError:
                angle = 90

            self.diagrama.ax.annotate(
                "C*=%0.1f" % ci, (NTU[29], e[30]), rotation=angle,
                size="medium", ha="left", va="bottom")

        self.diagrama.draw()

        img = image.imread('images/equation/%s.png' % flujo)
        self.image.set_data(img)
        self.refixImage()


class TemperatureEfectividad(Chart):
    title = QtWidgets.QApplication.translate(
        "pychemqt", "Heat Exchanger temperature effectiveness")
    flujo = [
        (QtWidgets.QApplication.translate("pychemqt", "Counterflow"), "CF"),
        (QtWidgets.QApplication.translate("pychemqt", "Parallelflow"), "PF"),
        (QtWidgets.QApplication.translate(
            "pychemqt", "Crossflow, both fluids unmixed"), "CrFunMix"),
        (QtWidgets.QApplication.translate(
            "pychemqt", "Crossflow, one fluid mixed, other unmixed"), "CrFSMix"),
        (QtWidgets.QApplication.translate(
            "pychemqt", "Crossflow, both fluids mixed"), "CrFMix"),
        (QtWidgets.QApplication.translate("pychemqt", "1-2 TEMA E"), "1-2TEMAE"),
        (QtWidgets.QApplication.translate(
            "pychemqt", "1-2 TEMA E, shell fluid divided"), "1-2TEMAE2"),
        (QtWidgets.QApplication.translate("pychemqt", "1-3 TEMA E"), "1-3TEMAE"),
        (QtWidgets.QApplication.translate("pychemqt", "1-4 TEMA E"), "1-4TEMAE"),
        (QtWidgets.QApplication.translate("pychemqt", "1-1 TEMA G"), "1-1TEMAG"),
        (QtWidgets.QApplication.translate("pychemqt", "1-2 TEMA G"), "1-2TEMAG"),
        (QtWidgets.QApplication.translate("pychemqt", "1-1 TEMA H"), "1-1TEMAH"),
        (QtWidgets.QApplication.translate("pychemqt", "1-2 TEMA H"), "1-2TEMAH"),
        (QtWidgets.QApplication.translate("pychemqt", "1-1 TEMA J"), "1-1TEMAJ"),
        (QtWidgets.QApplication.translate("pychemqt", "1-2 TEMA J"), "1-2TEMAJ"),
        (QtWidgets.QApplication.translate("pychemqt", "1-4 TEMA J"), "1-4TEMAJ")]

    mezclado = ("1", "2")
    PosLogo = (0, 1)

    def plot(self, indice):
        self.diagrama.ax.clear()
        self.diagrama.ax.set_xlim(0.1, 10)
        self.diagrama.ax.set_ylim(0, 1)
        self.diagrama.ax.set_xscale("log")

        self.diagrama.ax.set_title(QtWidgets.QApplication.translate(
            "pychemqt", "Heat Transfer Temperature Effectiveness"), size='12')
        self.diagrama.ax.set_xlabel("NTU", size='12')
        self.diagrama.ax.set_ylabel("P", size='14')
        self.diagrama.ax.set_xticklabels(["0.1", "1.0", "10"])
        xticklabels = ["0.2", "0.3", "", "0.5", "", "0.7", "", "", "2.0",
                       "3.0", "", "5.0", "", "7.0", "", ""]
        self.diagrama.ax.set_xticklabels(xticklabels, minor=True)

        flujo = self.flujo[indice][1]
        self.mixed.setVisible(flujo == "CrFSMix")
        kwargs = {}
        if flujo == "CrFSMix":
            kwargs["mixed"] = str(self.mixed.currentText())

        R = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6, 1.8,
             2., 2.5, 3., 4., 6., 8., 10., 15.]

        NTU = logspace(-1.5, 1, 100)
        for ri in R:
            e = [0]+[TemperatureEffectiveness(N, ri, flujo, **kwargs) for N in NTU[1:]]
            self.diagrama.plot(NTU, e, "k")
            self.diagrama.ax.annotate(" R=%0.1f" % ri, (NTU[-1], e[-1]),
                                      size="medium", ha="left", va="center")

#        F=[0.3]
#        for f in F:
#            p=[]
#            NTU=[]
#            for r in R:
#                func=lambda P: CorrectionFactor(P, r, flujo)-f
#                print func(0.)
#                pi=fsolve(func, 0.2)
#                p.append(pi)
#                NTU.append(NTU_fPR(p, r, flujo))

#            p=[fsolve(lambda P: CorrectionFactor(P, r, flujo)-f, 0.5)[0] for r in R]
#            NTU=[NTU_fPR(pi, ri) for pi, ni in zip(p, R)]
#            self.diagrama.plot(NTU, p, "--")

        self.diagrama.draw()

        if flujo == "CrFSMix" and self.mixed.currentIndex():
            img = image.imread('images/equation/%s2.png' % flujo)
        else:
            img = image.imread('images/equation/%s.png' % flujo)
        self.image.set_data(img)
        self.refixImage()


class F(Chart):
    title = QtWidgets.QApplication.translate("pychemqt", "ΔT Correction Factor")

    flujo = [
        (QtWidgets.QApplication.translate(
            "pychemqt", "Crossflow, both fluids unmixed"), "CrFunMix"),
        (QtWidgets.QApplication.translate(
            "pychemqt", "Crossflow, one fluid mixed, other unmixed"), "CrFSMix"),
        (QtWidgets.QApplication.translate(
            "pychemqt", "Crossflow, both fluids mixed"), "CrFMix"),
        (QtWidgets.QApplication.translate(
            "pychemqt", "1-2 pass shell and tube exchanger"), "1-2TEMAE")]

    mezclado = ("Cmin", "Cmax")
    PosLogo = (1, 1)
    PosImage = (0, 0)

    def plot(self, indice):
        self.diagrama.ax.clear()
        self.diagrama.ax.set_xlim(0, 1)
        self.diagrama.ax.set_ylim(0, 1)
        self.diagrama.ax.set_title(QtWidgets.QApplication.translate("pychemqt", "$\Delta T_{ml}$ Correction Factor", None), size='12')
        self.diagrama.ax.set_xlabel("$P=\\frac{T_{1o}-T_{1i}}{T_{2i}-T_{1i}}$", size='12')
        self.diagrama.ax.set_ylabel("F", size='14')

        flujo = self.flujo[indice][1]
#        self.mixed.setVisible(flujo=="CrFSMix")
        kwargs = {}
        if flujo == "CrFSMix":
            kwargs["mixed"] = str(self.mixed.currentText())

        R = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6,
             1.8, 2, 2.5, 3, 4, 6, 8, 10, 15, 20]
#        R=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]

        P = arange(0, 1.01, 0.01)
        for ri in R:
            f = [CorrectionFactor(p, ri, flujo, **kwargs) for p in P]
            self.diagrama.plot(P, f, "k")

#            fraccionx=P[90]-P[80]
#            fracciony=f[90]-f[80]
#            try:
#                angle=arctan(fracciony/fraccionx)*360/2/pi
#            except ZeroDivisionError:
#                angle=90
#            self.diagrama.ax.annotate("R=%0.1f" %ri, (P[90], f[90]), rotation=angle, size="medium", horizontalalignment="left", verticalalignment="bottom")

        self.diagrama.draw()

        img = image.imread('images/equation/%s.png' % flujo)
        self.image.set_data(img)
        self.refixImage()


class Phi(Chart):
    title = QtWidgets.QApplication.translate("pychemqt", "ψ", None)

    flujo=[#(QtWidgets.QApplication.translate("pychemqt", "Crossflow, both fluids unmixed"), "CrFunMix"),
#                (QtGui.QApplication.translate("pychemqt", "Crossflow, one fluid mixed, other unmixed"), "CrFSMix"),
#                (QtGui.QApplication.translate("pychemqt", "Crossflow, both fluids mixed"), "CrFMix"),
                (QtWidgets.QApplication.translate("pychemqt", "1-2 pass shell and tube exchanger"), "1-2TEMAE")]

    mezclado = ("Cmin", "Cmax")
    PosLogo = (1, 1)

    def plot(self, indice):
        self.diagrama.ax.clear()
        self.diagrama.ax.set_xlim(0, 1)
        self.diagrama.ax.set_ylim(0, 1)
        self.diagrama.ax.set_title(QtWidgets.QApplication.translate("pychemqt", "$\Delta T_{ml}$ Correction Factor", None), size='12')
        self.diagrama.ax.set_xlabel("$P=\\frac{T_{1o}-T_{1i}}{T_{2i}-T_{1i}}$", size='12')
        self.diagrama.ax.set_ylabel("F", size='14')

        flujo = self.flujo[indice][1]
#        self.mixed.setVisible(flujo=="CrFSMix")
        kwargs = {}
        if flujo == "CrFSMix":
            kwargs["mixed"] = str(self.mixed.currentText())

        R = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.4,
             1.6, 1.8, 2, 2.5, 3, 4, 6, 8, 10]

        P = arange(0, 1.01, 0.01)
        for ri in R:
            f = [Fi(p, ri, flujo, **kwargs) for p in P]
            self.diagrama.plot(P, f, "k")

        NTU = [0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.]
        for ntu in NTU:
            self.diagrama.plot([0, 1], [0, 1./ntu], "k", linestyle=":")
#            fraccionx=P[90]-P[80]
#            fracciony=f[90]-f[80]
#            try:
#                angle=arctan(fracciony/fraccionx)*360/2/pi
#            except ZeroDivisionError:
#                angle=90
#            self.diagrama.ax.annotate("R=%0.1f" %ri, (P[90], f[90]), rotation=angle, size="medium", horizontalalignment="left", verticalalignment="bottom")

        self.diagrama.draw()

        img = image.imread('images/equation/%s.png' % flujo)
        self.image.set_data(img)
        self.refixImage()


chartHE = (Efectividad, TemperatureEfectividad, F, Phi)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    dialogo = F()
    dialogo.show()
    sys.exit(app.exec_())
