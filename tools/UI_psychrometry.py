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
# Psychrometric graphic tools
#   - PsychroPlot: Plot widget for psychrometric chart
#   - PsychroInput: Widget with input for psychrometric state
#   - UI_Psychrometry: Psychrometric chart
###############################################################################


from configparser import ConfigParser
from functools import partial
import logging
import os
import json

from PyQt5 import QtCore, QtGui, QtWidgets
from scipy import pi, arctan

from lib.psycrometry import PsyState, PsychroState, _Pbar, _height
from lib.config import conf_dir
from lib.plot import mpl
from lib.unidades import (Temperature, Pressure, Length, Mass,
                          SpecificVolume, Enthalpy)
from lib.utilities import formatLine
from UI.widgets import Entrada_con_unidades


class PsychroPlot(mpl):
    """
    Plot widget for psychrometric chart
        Add custom margins
        Define a point for text state properties, to easy remove and redraw
    """
    def __init__(self, *args, **kwargs):
        mpl.__init__(self, *args, **kwargs)
        self.state = None

    def config(self, config):
        self.ax.set_autoscale_on(False)
        chart = config.getboolean("Psychr", "chart")
        xlabel = "Tdb, " + Temperature.text()
        ylabel = "%s, %s/%s" % (
            QtWidgets.QApplication.translate("pychemqt", "Absolute humidity"),
            Mass.text(), Mass.text())

        tmin = Temperature(config.getfloat("Psychr", "isotdbStart")).config()
        tmax = Temperature(config.getfloat("Psychr", "isotdbEnd")).config()
        wmin = config.getfloat("Psychr", "isowStart")
        wmax = config.getfloat("Psychr", "isowEnd")

        if chart:
            self.ax.set_xlabel(xlabel, size="large")
            self.ax.set_ylabel(ylabel, size="large")
            self.ax.set_xlim(tmin, tmax)
            self.ax.set_ylim(wmin, wmax)
            self.ax.yaxis.set_ticks_position("right")
            self.ax.yaxis.set_label_position("right")
            self.ax.figure.subplots_adjust(left=0.05, top=0.95)
        else:
            self.ax.set_xlabel(ylabel, size="large")
            self.ax.set_ylabel(xlabel, size="large")
            self.ax.set_xlim(wmin, wmax)
            self.ax.set_ylim(tmin, tmax)
            self.ax.xaxis.set_ticks_position("top")
            self.ax.xaxis.set_label_position("top")
            self.ax.figure.subplots_adjust(right=0.95, bottom=0.05)

        kw = formatLine(config, "Psychr", "crux")
        self.lx = self.ax.axhline(**kw)  # the horiz line
        self.ly = self.ax.axvline(**kw)  # the vert line

    def createCrux(self, state, chart):
        """Update horizontal and vertical lines to show click point"""
        self.state = state
        if chart:
            self.lx.set_ydata(state.w)
            self.ly.set_xdata(state.tdb.config())
        else:
            self.lx.set_ydata(state.tdb.config())
            self.ly.set_xdata(state.w)
        self.showPointData(state, chart)

    def clearCrux(self):
        self.lx.set_ydata(0)
        self.ly.set_xdata(0)

    def showPointData(self, state, chart=True):
        """Update data of current cursor point in plot annotates"""
        self.clearPointData()

        txt = []
        for key in ("tdb", "tdp", "twb", "HR", "w", "h", "v", "rho"):
            txt.append(("%s: %s" % (key, state.__getattribute__(key).str),))

        if chart:
            loc = "upper left"
        else:
            loc = "lower right"
        self.ax.table(cellText=txt, loc=loc, cellLoc="left", colLoc="left")
        self.ax.tables[0].auto_set_column_width(0)
        self.draw()

    def clearPointData(self):
        while self.ax.tables:
            self.ax.tables.pop()
        self.draw()


class PsychroInput(QtWidgets.QWidget):
    """Widget with parameter for psychrometric state"""
    parameters = ["tdb", "twb", "tdp", "w", "HR", "v", "h"]
    stateChanged = QtCore.pyqtSignal(PsyState)
    pressureChanged = QtCore.pyqtSignal()

    def __init__(self, state=None, readOnly=False, parent=None):
        """
        constructor
        optional state parameter to assign initial psychrometric state
        """
        super(PsychroInput, self).__init__(parent)

        self.state = PsychroState(P=101325)

        layout = QtWidgets.QGridLayout(self)
        self.checkPresion = QtWidgets.QRadioButton(
                QtWidgets.QApplication.translate("pychemqt", "Pressure"))
        layout.addWidget(self.checkPresion, 1, 1, 1, 1)
        self.P = Entrada_con_unidades(Pressure, value=101325)
        self.P.valueChanged.connect(self.changePressure)
        layout.addWidget(self.P, 1, 2, 1, 1)
        self.checkAltitud = QtWidgets.QRadioButton(
                QtWidgets.QApplication.translate("pychemqt", "Altitude"))
        layout.addWidget(self.checkAltitud, 2, 1, 1, 1)
        self.z = Entrada_con_unidades(Length, value=0)
        self.checkPresion.toggled.connect(self.P.setEnabled)
        self.checkAltitud.toggled.connect(self.z.setEnabled)
        self.z.valueChanged.connect(self.changeAltitude)
        self.checkPresion.setChecked(True)
        self.z.setEnabled(False)
        layout.addWidget(self.z, 2, 2, 1, 1)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            3, 1, 1, 2)

        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Select point")), 4, 1, 1, 2)
        self.variables = QtWidgets.QComboBox()
        for txt in PsyState.TEXT_MODE:
            self.variables.addItem(txt)
        self.variables.currentIndexChanged.connect(self.updateInputs)
        layout.addWidget(self.variables, 5, 1, 1, 2)

        layout.addWidget(QtWidgets.QLabel("Tdb:"), 6, 1, 1, 1)
        self.tdb = Entrada_con_unidades(Temperature)
        self.tdb.valueChanged.connect(partial(self.updateKwargs, "tdb"))
        layout.addWidget(self.tdb, 6, 2, 1, 1)
        layout.addWidget(QtWidgets.QLabel("Twb:"), 7, 1, 1, 1)
        self.twb = Entrada_con_unidades(Temperature)
        self.twb.valueChanged.connect(partial(self.updateKwargs, "twb"))
        layout.addWidget(self.twb, 7, 2, 1, 1)
        layout.addWidget(QtWidgets.QLabel("Tdp:"), 8, 1, 1, 1)
        self.tdp = Entrada_con_unidades(Temperature)
        self.tdp.valueChanged.connect(partial(self.updateKwargs, "tdp"))
        layout.addWidget(self.tdp, 8, 2, 1, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Humidity Ratio:")), 9, 1, 1, 1)
        self.w = Entrada_con_unidades(float, textounidad="kgw/kgda")
        self.w.valueChanged.connect(partial(self.updateKwargs, "w"))
        layout.addWidget(self.w, 9, 2, 1, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Relative humidity:")), 10, 1, 1, 1)
        self.HR = Entrada_con_unidades(float, textounidad="%")
        self.HR.valueChanged.connect(partial(self.updateKwargs, "HR"))
        layout.addWidget(self.HR, 10, 2, 1, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Volume")), 11, 1, 1, 1)
        self.v = Entrada_con_unidades(SpecificVolume)
        self.v.valueChanged.connect(partial(self.updateKwargs, "v"))
        layout.addWidget(self.v, 11, 2, 1, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Enthalpy")), 12, 1, 1, 1)
        self.h = Entrada_con_unidades(Enthalpy)
        self.h.valueChanged.connect(partial(self.updateKwargs, "h"))
        layout.addWidget(self.h, 12, 2, 1, 1)
        layout.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            13, 1, 1, 2)

        self.setReadOnly(readOnly)
        self.updateInputs(0)
        if state:
            self.setState(state)

    def updateInputs(self, index):
        """Update inputs appearance to highlight active"""
        for par in self.parameters:
            self.__getattribute__(par).setReadOnly(True)
            self.__getattribute__(par).setResaltado(False)
        for par in PsyState.VAR_NAME[index]:
            self.__getattribute__(par).setReadOnly(False)
            self.__getattribute__(par).setResaltado(True)

        index = self.variables.currentIndex()
        kwargs = {"P": self.P.value}
        for par in PsyState.VAR_NAME[index]:
            if self.__getattribute__(par).value:
                kwargs[par] = self.state.__getattribute__(par)
        self.state = PsychroState(**kwargs)

    def setReadOnly(self, readOnly):
        self.checkPresion.setEnabled(not readOnly)
        self.checkAltitud.setEnabled(not readOnly)
        self.P.setReadOnly(readOnly)
        self.z.setReadOnly(readOnly)
        self.variables.setEnabled(not readOnly)
        for par in self.parameters:
            self.__getattribute__(par).setReadOnly(True)
            self.__getattribute__(par).setResaltado(False)

    def updateKwargs(self, key, value):
        """Update kwargs of state instance, if its correctly defined show it"""
        kwargs = {key: value}
        self.state(**kwargs)
        if self.state.status:
            self.setState(self.state)
            self.stateChanged.emit(self.state)

    def setState(self, state):
        """Fill data input with state properties"""
        self.state = state
        if state.w < state.ws:
            for p in self.parameters:
                self.__getattribute__(p).setValue(state.__getattribute__(p))

    def changePressure(self, value):
        """Change pressure to global plot and for states"""
        self.z.setValue(_height(value))
        self.state = PsychroState(P=value)
        self.pressureChanged.emit()

    def changeAltitude(self, value):
        """Change pressure through altitude and ICAO equation"""
        presion = _Pbar(value)
        self.P.setValue(presion)
        self.state = PsychroState(P=value)
        self.pressureChanged.emit()


class UI_Psychrometry(QtWidgets.QDialog):
    """Psychrometric charts tool"""
    def __init__(self, parent=None):
        super(UI_Psychrometry, self).__init__(parent)
        self.setWindowTitle(QtWidgets.QApplication.translate(
            "pychemqt", "Psychrometric chart"))
        self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"] + "/images/button/psychrometric.png")))
        self.showMaximized()

        layout = QtWidgets.QGridLayout(self)
        self.plt = PsychroPlot(self, width=100, height=1, dpi=90)
        self.plt.fig.canvas.mpl_connect('button_press_event', self.click)
        layout.addWidget(self.plt, 1, 3, 2, 2)

        self.inputs = PsychroInput()
        self.inputs.stateChanged.connect(self.plt.createCrux)
        self.inputs.pressureChanged.connect(self.plot)
        layout.addWidget(self.inputs, 1, 1, 1, 1)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 2, 1)

        self.buttonShowToolbox = QtWidgets.QToolButton()
        self.buttonShowToolbox.setCheckable(True)
        self.buttonShowToolbox.toggled.connect(self.showToolBar)
        layout.addWidget(self.buttonShowToolbox, 1, 2, 2, 1)
        self.line = QtWidgets.QFrame()
        self.line.setFrameShape(QtWidgets.QFrame.VLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        layout.addWidget(self.line, 1, 3, 3, 1)

        self.progressBar = QtWidgets.QProgressBar()
        self.progressBar.setVisible(False)
        layout.addWidget(self.progressBar, 3, 3)
        self.status = QtWidgets.QLabel()
        layout.addWidget(self.status, 3, 3)

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
        btBox.rejected.connect(self.reject)
        layout.addWidget(btBox, 3, 4)
        btBox.layout().insertWidget(0, butonPNG)
        btBox.layout().insertWidget(0, butonConfig)

        self.showToolBar(False)
        self.Preferences = ConfigParser()
        self.Preferences.read(conf_dir+"pychemqtrc")
        self.plot()
        logging.info(QtWidgets.QApplication.translate(
            "pychemqt", "Started psychrometric chart tool"))

    def configure(self):
        from UI.prefPsychrometric import Dialog
        dlg = Dialog(self.Preferences)
        if dlg.exec_():
            self.Preferences = dlg.value(self.Preferences)
            self.Preferences.write(open(conf_dir+"pychemqtrc", "w"))
            self.plot()

    def showToolBar(self, checked):
        """Show/Hide left toolbar with additional funcionality"""
        self.inputs.setVisible(checked)
        if checked:
            image = "arrow-left-double.png"
        else:
            image = "arrow-right-double.png"
        self.buttonShowToolbox.setIcon(QtGui.QIcon(
            os.environ["pychemqt"] + os.path.join("images", "button", image)))

    def drawlabel(self, name, t, W, label, unit):
        """
        Draw annotation for isolines
            name: name of isoline
            t: x array of line
            W: y array of line
            label: text value to draw
            unit: text units to draw
        """
        if self.Preferences.getboolean("Psychr", name+"label"):
            TMIN = self.Preferences.getfloat("Psychr", "isotdbStart")
            TMAX = self.Preferences.getfloat("Psychr", "isotdbEnd")
            tmin = Temperature(TMIN).config()
            tmax = Temperature(TMAX).config()
            wmin = self.Preferences.getfloat("Psychr", "isowStart")
            wmax = self.Preferences.getfloat("Psychr", "isowEnd")
            if self.Preferences.getboolean("Psychr", "chart"):
                x = tmax-tmin
                y = wmax-wmin
                i = 0
                for ti, wi in zip(t, W):
                    if tmin <= ti <= tmax and wmin <= wi <= wmax:
                        i += 1
            else:
                x = wmax-wmin
                y = tmax-tmin
                i = 0
                for ti, wi in zip(t, W):
                    if tmin <= wi <= tmax and wmin <= ti <= wmax:
                        i += 1

            if isinstance(label, float):
                label = "%4g" % label

            if self.Preferences.getboolean("Psychr", name+"units"):
                label += unit
            pos = self.Preferences.getfloat("Psychr", name+"position")
            p = int(i*pos/100-1)
            rot = arctan((W[p]-W[p-1])/y/(t[p]-t[p-1])*x)*360/2/pi
            self.plt.ax.annotate(label, (t[p], W[p]), rotation=rot,
                                 size="small", ha="center", va="center")

    def plot(self):
        """Plot chart"""
        self.plt.clearPointData()
        self.plt.ax.clear()
        chart = self.Preferences.getboolean("Psychr", "chart")
        self.plt.config(self.Preferences)
        filename = conf_dir+"%s_%i.json" % (
            PsychroState().__class__.__name__, self.inputs.P.value)
        if os.path.isfile(filename):
            with open(filename, "r") as archivo:
                data = json.load(archivo)
                self.status.setText(QtWidgets.QApplication.translate(
                    "pychemqt", "Loading cached data..."))
                QtWidgets.QApplication.processEvents()
        else:
            self.progressBar.setVisible(True)
            self.status.setText(QtWidgets.QApplication.translate(
                "pychemqt", "Calculating data..."))
            QtWidgets.QApplication.processEvents()
            data = PsychroState.calculatePlot(self)
            with open(filename, "w") as file:
                json.dump(data, file, indent=4)

            self.progressBar.setVisible(False)
        self.status.setText(
            QtWidgets.QApplication.translate("pychemqt", "Plotting..."))
        QtWidgets.QApplication.processEvents()

        # Saturation line
        t = [Temperature(ti).config() for ti in data["t"]]
        Hs = data["Hs"]
        format = formatLine(self.Preferences, "Psychr", "saturation")
        if chart:
            self.plt.plot(t, Hs, **format)
        else:
            self.plt.plot(Hs, t, **format)

        # Iso dew bulb temperaure lines, vertial lines in normal h-Tdb plot,
        # vertical in mollier diagram
        format = formatLine(self.Preferences, "Psychr", "isotdb")
        for i, T in enumerate(t):
            if chart:
                self.plt.plot([T, T], [0, Hs[i]], **format)
            else:
                self.plt.plot([0, Hs[i]], [T, T], **format)

        # Iso humidity lines, horizontal lines in normal h-Tdb plot, vertical
        # in mollier diagram
        H = data["H"]
        th = data["th"]
        tm = Temperature(self.Preferences.getfloat("Psychr", "isotdbEnd"))
        format = formatLine(self.Preferences, "Psychr", "isow")
        for i, H in enumerate(H):
            ts = Temperature(th[i]).config()
            if chart:
                self.plt.plot([ts, tm.config()], [H, H], **format)
            else:
                self.plt.plot([H, H], [ts, tm.config()], **format)

        # Iso relative humidity lines
        format = formatLine(self.Preferences, "Psychr", "isohr")
        for Hr, H0 in list(data["Hr"].items()):
            if chart:
                self.plt.plot(t, H0, **format)
                self.drawlabel("isohr", t, H0, Hr, "%")
            else:
                self.plt.plot(H0, t, **format)
                self.drawlabel("isohr", H0, t, Hr, "%")

        # Iso wet bulb temperature lines
        format = formatLine(self.Preferences, "Psychr", "isotwb")
        for T, (H, Tw) in list(data["Twb"].items()):
            value = Temperature(T).config()
            Tw_conf = [Temperature(Twi).config() for Twi in Tw]
            txt = Temperature.text()
            if chart:
                self.plt.plot(Tw_conf, H, **format)
                self.drawlabel("isotwb", Tw_conf, H, value, txt)
            else:
                self.plt.plot(H, Tw_conf, **format)
                self.drawlabel("isotwb", H, Tw_conf, value, txt)

        # Isochor lines
        format = formatLine(self.Preferences, "Psychr", "isochor")
        for v, (Td, H) in list(data["v"].items()):
            value = SpecificVolume(v).config()
            Td_conf = [Temperature(Tdi).config() for Tdi in Td]
            txt = SpecificVolume.text()
            if chart:
                self.plt.plot(Td_conf, H, **format)
                self.drawlabel("isochor", Td_conf, H, value, txt)
            else:
                self.plt.plot(H, Td_conf, **format)
                self.drawlabel("isochor", H, Td_conf, value, txt)

        if self.plt.state:
            self.plt.createCrux(self.plt.state, chart)
        self.plt.draw()
        self.status.setText("%s %s" % (
            QtWidgets.QApplication.translate("pychemqt", "Using"),
            PsychroState().__class__.__name__[3:]))

    def click(self, event):
        """Update input and graph annotate when mouse click over chart"""
        if event.xdata and event.ydata:
            chart = self.Preferences.getboolean("Psychr", "chart")
            if chart:
                state = self.createState(event.xdata, event.ydata)
            else:
                state = self.createState(event.ydata, event.xdata)
            if state.w <= state.ws:
                self.inputs.setState(state)
                self.plt.createCrux(state, chart)
            else:
                self.plt.clearCrux()
                self.plt.clearPointData()

    def createState(self, x, y):
        """Create psychrometric state from click or mouse position"""
        tdb = Temperature(x, "conf")
        punto = PsychroState(P=self.inputs.P.value, tdb=tdb, w=y)
        return punto

    def setProgressValue(self, value):
        self.progressBar.setValue(value)
        QtWidgets.QApplication.processEvents()


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    aireHumedo = UI_Psychrometry()
    sys.exit(app.exec_())
