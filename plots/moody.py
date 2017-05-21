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
# Moody plot
###############################################################################


from configparser import ConfigParser
import json
import os

from PyQt5 import QtGui, QtWidgets
from scipy import logspace, log10
from matplotlib.patches import ConnectionPatch
from matplotlib import image

from lib.config import conf_dir, IMAGE_PATH
from lib.friction import f_list, eD
from lib.plot import mpl
from lib.utilities import formatLine, representacion
from UI.prefMoody import ConfigDialog
from UI.widgets import Entrada_con_unidades


Re_laminar = [600, 2400]
Re_turbulent = logspace(log10(2400), 8, 50)
Re_fully = logspace(log10(4000), 8, 50)


def calculate(config):
    fanning = config.getboolean("Moody", "fanning")
    method = config.getint("Moody", "method")
    eD = eval(config.get("Moody", "eD"))
    F = f_list[method]

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
    dat["fully"] = [(1/(1.14-2*log10(3500/R)))**2/x for R in Re_fully]

    # Save to file
    with open(conf_dir+"moody.dat", "w") as file:
        json.dump(dat, file, indent=4)


class Chart(QtWidgets.QDialog):
    """Generic chart dialog"""

    def __init__(self, parent=None):
        super(Chart, self).__init__(parent)
        self.setWindowTitle(self.title)
        layout = QtWidgets.QGridLayout(self)
        layout.setColumnStretch(3, 1)
        self.plt = mpl(self)
        self.plt.fig.canvas.mpl_connect('button_press_event', self.click)
        layout.addWidget(self.plt, 2, 1, 1, 4)

        btBox = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Close)
        butonPNG = QtWidgets.QPushButton(QtGui.QIcon(
            os.path.join(IMAGE_PATH, "button", "image.png")),
            QtWidgets.QApplication.translate("pychemqt", "Save as PNG"))
        butonPNG.clicked.connect(self.plt.savePNG)
        butonConfig = QtWidgets.QPushButton(QtGui.QIcon(
            os.path.join(IMAGE_PATH, "button", "configure.png")),
            QtWidgets.QApplication.translate("pychemqt", "Configure"))
        butonConfig.clicked.connect(self.configure)
        butonCalculate = QtWidgets.QPushButton(QtGui.QIcon(
            os.path.join(IMAGE_PATH, "button", "calculator.png")),
            QtWidgets.QApplication.translate("pychemqt", "Calculate point"))
        butonCalculate.clicked.connect(self.calculate)
        btBox.rejected.connect(self.reject)
        btBox.layout().insertWidget(0, butonPNG)
        btBox.layout().insertWidget(0, butonCalculate)
        btBox.layout().insertWidget(0, butonConfig)
        layout.addWidget(btBox, 3, 1, 1, 4)

        self.Preferences = ConfigParser()
        self.Preferences.read(conf_dir+"pychemqtrc")
        self.config()
        self.plot()

        logo = image.imread(os.path.join(IMAGE_PATH, "pychemqt_98.png"))
        self.logo = self.plt.fig.figimage(logo, 0, 0, zorder=1, alpha=0.5)
        self.showMaximized()

    def configure(self):
        dlg = self.configDialog(self.Preferences)
        if dlg.exec_():
            self.Preferences = dlg.value(self.Preferences)
            self.Preferences.write(open(conf_dir+"pychemqtrc", "w"))
            self.plot()

    def config(self):
        """Initialization action in plot don't neccesary to update in plot"""
        pass

    def click(self, event):
        """Define functionality when the mouse click in plot"""
        pass

    def plot(self):
        """Plot procedure"""
        pass

    def calculate(self):
        """Define the functionality when click the calculate point button"""
        pass

    def resizeEvent(self, event):
        """Implement resizeEvent to recalculate the new position of image"""
        self.refixImage()
        QtWidgets.QDialog.resizeEvent(self, event)

    def refixImage(self, event=None):
        """Recalculate image position as resize/move window"""
        xmin, xmax = self.plt.ax.get_xlim()
        ymin, ymax = self.plt.ax.get_ylim()
        x, y = self.plt.ax.transData.transform_point(self.PosLogo)
        hlogo, wlogo = self.logo.get_size()
        if self.PosLogo[0] == 1:
            x -= wlogo
        if self.PosLogo[1] == 1:
            y -= hlogo

        self.logo.ox = x
        self.logo.oy = y

        self.plt.draw()


class Moody(Chart):
    """Moody chart dialog"""
    title = QtWidgets.QApplication.translate("pychemqt", "Moody Diagram")
    configDialog = ConfigDialog
    PosLogo = 5e3, 1.2e-2

    def config(self):
        """Initialization action in plot don't neccesary to update in plot"""

        txt = QtWidgets.QApplication.translate(
            "pychemqt", "Relative roughness") + ", "+r"$r=\frac{\epsilon}{D}$"
        self.plt.fig.text(0.97, 0.5, txt, rotation=90, size='10',
                          va="center", ha="center")
        self.plt.fig.subplots_adjust(
            left=0.1, right=0.9, bottom=0.12, top=0.98)
        self.note = None

    def click(self, event):
        """Update input and graph annotate when mouse click over chart"""
        Re = event.xdata
        f = event.ydata
        ed = None
        method = self.Preferences.getint("Moody", "method")
        fanning = self.Preferences.getboolean("Moody", "fanning")
        F = f_list[method]

        if Re < 2400:
            if fanning:
                f = 16/Re
            else:
                f = 64/Re
        else:
            f_min = F(Re, 0)
            if f < f_min:
                Re = f = 0
            else:
                ed = eD(Re, f)

        self.createCrux(Re, f, ed)

    def _txt(self, Re, f, ed):
        if ed is None:
            txt = "Re: %0.4g\nf: %0.4g" % (Re, f)
        else:
            txt = "Re: %0.4g\ne/D: %0.4g\nf: %0.4g" % (Re, ed, f)
        return txt

    def createCrux(self, Re, f, ed):
        """Create a crux in selected point of plot and show data at bottom
        right corner"""
        if f and Re:
            txt = self._txt(Re, f, ed)

        self.plt.lx.set_ydata(f)
        self.plt.ly.set_xdata(Re)
        if self.note:
            self.note.remove()
            self.note = None
        if f and Re:
            self.note = self.plt.fig.text(0.85, 0.08, txt, size="6", va="top")
        self.plt.draw()

    def plot(self):
        """Plot the Moody chart using the indicate method """
        fanning = self.Preferences.getboolean("Moody", "fanning")
        method = self.Preferences.getint("Moody", "method")

        if fanning:
            x = 4
        else:
            x = 1

        self.plt.ax.set_autoscale_on(False)
        self.plt.ax.clear()

        kw = formatLine(self.Preferences, "Moody", "crux")
        self.plt.lx = self.plt.ax.axhline(**kw)  # the horiz line
        self.plt.ly = self.plt.ax.axvline(**kw)  # the vert line

        xlabel = QtWidgets.QApplication.translate(
            "pychemqt", "Reynolds number") + ", " + r"$Re=\frac{V\rho D}{\mu}$"
        self.plt.ax.set_xlabel(xlabel, ha='center', size='10')
        if fanning:
            ylabel = QtWidgets.QApplication.translate(
                "pychemqt", "Fanning Friction factor")
            formula = r"$f_f=\frac{2hDg}{LV^2}$"
        else:
            ylabel = QtWidgets.QApplication.translate(
                "pychemqt", "Darcy Friction factor")
            formula = r"$f_d=\frac{2hDg}{LV^2}$"
        self.plt.ax.set_ylabel(ylabel+",  " + formula, size='10')
        self.plt.ax.grid(True)
        self.plt.ax.set_xlim(600, 1e8)
        self.plt.ax.set_ylim(0.008/x, 0.11/x)
        self.plt.ax.set_xscale("log")
        self.plt.ax.set_yscale("log")

        xticks = [7e2, 8e2, 9e2, 1e3, 2e3, 3e3, 4e3, 5e3, 6e3, 7e3, 8e3, 9e3,
                  1e4, 2e4, 3e4, 4e4, 5e4, 6e4, 7e4, 8e4, 9e4, 1e5, 2e5, 3e5,
                  4e5, 5e5, 6e5, 7e5, 8e5, 9e5, 1e6, 2e6, 3e6, 4e6, 5e6, 6e6,
                  7e6, 8e6, 9e6, 1e7, 2e7, 3e7, 4e7, 5e7, 6e7, 7e7, 8e7, 9e7]
        self.plt.ax.set_xticks(xticks, minor=False)
        self.plt.ax.set_yticks([], minor=True)

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
        self.plt.ax.set_yticks(yticks)
        self.plt.ax.set_yticklabels(ytickslabel)

        if not os.path.isfile(conf_dir+"moody.dat"):
            calculate(self.Preferences)

        load = False
        with open(conf_dir+"moody.dat", "r") as file:
            try:
                dat = json.load(file)
            except ValueError:
                calculate(self.Preferences)
                load = True

            if dat["fanning"] != fanning or dat["method"] != method:
                calculate(self.Preferences)
                load = True

        # Reload file if it's created in last with statement
        if load:
            with open(conf_dir+"moody.dat", "r") as file:
                dat = json.load(file)

        # Plot data
        kw = formatLine(self.Preferences, "Moody", "line")
        self.plt.ax.plot(Re_laminar, dat["laminar"], **kw)
        for ed, f in dat["turbulent"].items():
            self.plt.ax.plot(Re_turbulent, f, **kw)
            title = " " + representacion(ed, tol=4.5)
            if f[-1] > 0.008/x:
                self.plt.ax.text(1e8, f[-1], title, size="x-small",
                                 ha='left', va='center')

        self.plt.ax.plot(Re_fully, dat["fully"], "k", lw=0.5, ls=":")

        # Add explicative legend
        # Laminar zone
        self.plt.ax.add_artist(
            ConnectionPatch((600, 0.01/x), (2400, 0.01/x), "data", "data",
                            arrowstyle="<|-|>", mutation_scale=15, fc="w"))
        txt = QtWidgets.QApplication.translate("pychemqt", "Laminar flux")
        self.plt.ax.text(1200, 0.0095/x, txt, size="small", va="top",
                         ha="center", backgroundcolor="#ffffff")

        # Critic zone
        self.plt.ax.add_artist(
            ConnectionPatch((2400, 0.01/x), (4000, 0.01/x), "data", "data",
                            arrowstyle="<|-|>", mutation_scale=15, fc="w"))
        txt = QtWidgets.QApplication.translate("pychemqt", "Critic\nzone")
        self.plt.ax.text(3200, 0.0095/x, txt, size="small", va="top",
                         ha="center", backgroundcolor="#ffffff")

        # Transition zone
        self.plt.ax.add_artist(
            ConnectionPatch((4000, 0.095/x), (40000, 0.095/x), "data", "data",
                            arrowstyle="<|-|>", mutation_scale=15, fc="w"))
        txt = QtWidgets.QApplication.translate("pychemqt", "Transition Zone")
        self.plt.ax.text(11000, 0.098/x, txt, size="small", va="bottom",
                         ha="center", backgroundcolor="#ffffff")

        # Turbulent zone
        self.plt.ax.add_artist(
            ConnectionPatch((40000, 0.095/x), (9.9e7, 0.095/x), "data", "data",
                            arrowstyle="<|-|>", mutation_scale=15, fc="w"))
        txt = QtWidgets.QApplication.translate(
            "pychemqt", "Turbulent flux fully developed")
        self.plt.ax.text(1e6, 0.098/x, txt, size="small", va="bottom",
                         ha="center", backgroundcolor="#ffffff")

        # Smooth tubes
        self.plt.ax.add_artist(
            ConnectionPatch((1e6, 0.0095/x), (1.5e6, 0.011/x), "data", "data",
                            arrowstyle="<|-|>", mutation_scale=15, fc="w"))
        txt = QtWidgets.QApplication.translate("pychemqt", "Smooth tubes")
        self.plt.ax.text(1e6, 0.009/x, txt, size="small", va="top",
                         ha="center", backgroundcolor="#ffffff")

        # Laminar equation
        if fanning:
            txt = r"$f=\frac{16}{Re}$"
        else:
            txt = r"$f=\frac{64}{Re}$"
        self.plt.ax.text(1.4e3, 0.042/x, txt, size="12",
                         va="top", ha="center", rotation=-66)
        # TODO: Calculate the angle dynamically

        self.plt.draw()

    def calculate(self):
        dlg = CalculateDialog()
        if dlg.exec_():
            Re = dlg.Re.value
            f = dlg.f.value
            eD = dlg.eD.value
            self.createCrux(Re, f, eD)


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
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Close)
        self.buttonBox.accepted.connect(self.accept)
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
