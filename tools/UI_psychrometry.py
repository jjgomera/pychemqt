#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2025, Juan José Gómez Romera <jjgomera@gmail.com>

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


.. include:: UI_psychrometry.rst


The module include all related psychrometry chart functionality
    * :class:`UI_Psychrometry`: Psychrometric chart
    * :class:`PsychroPlot`: Plot widget for psychrometric chart
    * :class:`PsychroInput`: Widget with input for psychrometric state

and its configuration

    * :class:`Config`: Psychrometric chart configuration
    * :class:`ConfigDialog`: Dialog tool for standalone use
'''


from configparser import ConfigParser
from functools import partial
import json
import logging
from math import pi
import os

from numpy import arctan

from lib.config import conf_dir, IMAGE_PATH
from lib.plot import PlotWidget
from lib.psycrometry import PsyState, PsychroState, _Pbar, _height
from lib.unidades import (Temperature, Pressure, Length, Mass,
                          SpecificVolume, Enthalpy)
from lib.utilities import formatLine
from tools.qt import QtCore, QtGui, QtWidgets, translate
from tools.UI_Tables.prefMEOS import Isolinea
from UI.widgets import Entrada_con_unidades, LineConfig


class PsychroPlot(PlotWidget):
    """
    Plot widget for psychrometric chart
        Add custom margins
        Define a point for text state properties, to easy remove and redraw
    """
    def __init__(self, *args, **kwargs):
        PlotWidget.__init__(self, *args, **kwargs)
        self.state = None

    def config(self, config):
        """Apply configuration to plot"""
        self.ax.set_autoscale_on(False)
        chart = config.getboolean("Psychr", "chart")
        xlabel = "Tdb, " + Temperature.text()
        ylabel = f"{self.tr('Absolute humidity')}, {Mass.text()}/{Mass.text()}"

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
            self.lx.set_ydata([state.w])
            self.ly.set_xdata([state.tdb.config()])
        else:
            self.lx.set_ydata([state.tdb.config()])
            self.ly.set_xdata([state.w])
        self.showPointData(state, chart)

    def clearCrux(self):
        """Clear crux lines for click interaction in plot"""
        self.lx.set_ydata([0])
        self.ly.set_xdata([0])

    def showPointData(self, state, chart=True):
        """Update data of current cursor point in plot annotates"""
        self.clearPointData()

        txt = []
        for key in ("tdb", "tdp", "twb", "HR", "w", "h", "v", "rho"):
            txt.append((f"{key}: {getattr(state,key).str}",))

        if chart:
            loc = "upper left"
        else:
            loc = "lower right"
        self.ax.table(cellText=txt, loc=loc, cellLoc="left", colLoc="left",
                      edges="open", fontsize=8)
        self.ax.tables[0].auto_set_column_width(0)
        self.ax.tables[0].auto_set_font_size(True)
        self.draw()

    def clearPointData(self):
        """Delete point data from plot"""
        while self.ax.tables:
            self.ax.tables[0].remove()
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
        super().__init__(parent)

        self.state = PsychroState(P=101325)

        layout = QtWidgets.QGridLayout(self)
        self.checkPresion = QtWidgets.QRadioButton(self.tr("Pressure"))
        layout.addWidget(self.checkPresion, 1, 1, 1, 1)
        self.P = Entrada_con_unidades(Pressure, value=101325)
        self.P.valueChanged.connect(self.changePressure)
        layout.addWidget(self.P, 1, 2, 1, 1)
        self.checkAltitud = QtWidgets.QRadioButton(self.tr("Altitude"))
        layout.addWidget(self.checkAltitud, 2, 1, 1, 1)
        self.z = Entrada_con_unidades(Length, value=0)
        self.checkPresion.toggled.connect(self.P.setEnabled)
        self.checkAltitud.toggled.connect(self.z.setEnabled)
        self.z.valueChanged.connect(self.changeAltitude)
        self.checkPresion.setChecked(True)
        self.z.setEnabled(False)
        layout.addWidget(self.z, 2, 2, 1, 1)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 3, 1, 1, 2)

        layout.addWidget(QtWidgets.QLabel(self.tr("Select point")), 4, 1, 1, 2)
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
        layout.addWidget(QtWidgets.QLabel(
            self.tr("Humidity Ratio:")), 9, 1, 1, 1)
        self.w = Entrada_con_unidades(float, textounidad="kgw/kgda")
        self.w.valueChanged.connect(partial(self.updateKwargs, "w"))
        layout.addWidget(self.w, 9, 2, 1, 1)
        layout.addWidget(QtWidgets.QLabel(
            self.tr("Relative humidity:")), 10, 1, 1, 1)
        self.HR = Entrada_con_unidades(float, textounidad="%")
        self.HR.valueChanged.connect(partial(self.updateKwargs, "HR"))
        layout.addWidget(self.HR, 10, 2, 1, 1)
        layout.addWidget(QtWidgets.QLabel(
            self.tr("Volume")), 11, 1, 1, 1)
        self.v = Entrada_con_unidades(SpecificVolume)
        self.v.valueChanged.connect(partial(self.updateKwargs, "v"))
        layout.addWidget(self.v, 11, 2, 1, 1)
        layout.addWidget(QtWidgets.QLabel(
            self.tr("Enthalpy")), 12, 1, 1, 1)
        self.h = Entrada_con_unidades(Enthalpy)
        self.h.valueChanged.connect(partial(self.updateKwargs, "h"))
        layout.addWidget(self.h, 12, 2, 1, 1)
        layout.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 13, 1, 1, 2)

        self.setReadOnly(readOnly)
        self.updateInputs(0)
        if state:
            self.setState(state)

    def updateInputs(self, index):
        """Update inputs appearance to highlight active"""
        for par in self.parameters:
            getattr(self, par).setReadOnly(True)
            getattr(self, par).setResaltado(False)
        for par in PsyState.VAR_NAME[index]:
            getattr(self, par).setReadOnly(False)
            getattr(self, par).setResaltado(True)

        index = self.variables.currentIndex()
        kwargs = {"P": self.P.value}
        for par in PsyState.VAR_NAME[index]:
            if getattr(self, par).value:
                kwargs[par] = getattr(self.state, par)
        self.state = PsychroState(**kwargs)

    def setReadOnly(self, readOnly):
        """Set readOnly all widget"""
        self.checkPresion.setEnabled(not readOnly)
        self.checkAltitud.setEnabled(not readOnly)
        self.P.setReadOnly(readOnly)
        self.z.setReadOnly(readOnly)
        self.variables.setEnabled(not readOnly)
        for par in self.parameters:
            getattr(self, par).setReadOnly(True)
            getattr(self, par).setResaltado(False)

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
                getattr(self, p).setValue(getattr(state, p))

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
        super().__init__(parent)
        self.__TEXTSTATUS__ = self.tr("Launched humid air properties aplication")
        self.setWindowTitle(self.tr("Psychrometric chart"))
        self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(
            os.path.join(IMAGE_PATH, "button", "psychrometric.png"))))

        layout = QtWidgets.QGridLayout(self)
        self.plt = PsychroPlot(parent=self, width=100, height=1, dpi=90)
        self.plt.fig.canvas.mpl_connect('button_press_event', self.click)
        layout.addWidget(self.plt, 1, 3, 2, 2)

        self.inputs = PsychroInput()
        self.inputs.stateChanged.connect(self.plt.createCrux)
        self.inputs.pressureChanged.connect(self.plot)
        layout.addWidget(self.inputs, 1, 1, 1, 1)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 2, 1)

        self.buttonShowToolbox = QtWidgets.QToolButton()
        self.buttonShowToolbox.setCheckable(True)
        self.buttonShowToolbox.toggled.connect(self.showToolBar)
        layout.addWidget(self.buttonShowToolbox, 1, 2, 2, 1)
        self.line = QtWidgets.QFrame()
        self.line.setFrameShape(QtWidgets.QFrame.Shape.VLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Shadow.Sunken)
        layout.addWidget(self.line, 1, 3, 3, 1)

        self.progressBar = QtWidgets.QProgressBar()
        self.progressBar.setVisible(False)
        layout.addWidget(self.progressBar, 3, 3)
        self.status = QtWidgets.QLabel()
        layout.addWidget(self.status, 3, 3)

        btBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Close)
        butonPNG = QtWidgets.QPushButton(QtGui.QIcon(os.path.join(
            os.environ["pychemqt"], "images", "button", "image.png")),
            self.tr("Save as PNG"))
        butonPNG.clicked.connect(self.plt.savePNG)
        butonConfig = QtWidgets.QPushButton(QtGui.QIcon(os.path.join(
            os.environ["pychemqt"], "images", "button", "configure.png")),
            self.tr("Configure"))
        butonConfig.clicked.connect(self.configure)
        btBox.rejected.connect(self.reject)
        layout.addWidget(btBox, 3, 4)
        btBox.layout().insertWidget(0, butonPNG)
        btBox.layout().insertWidget(0, butonConfig)

        self.showToolBar(False)
        self.Preferences = ConfigParser()
        self.Preferences.read(conf_dir+"pychemqtrc")
        self.plot()
        logging.info(self.tr("Started psychrometric chart tool"))
        self.showMaximized()

    def configure(self):
        """Show configuration dialog"""
        dlg = ConfigDialog(self.Preferences)
        if dlg.exec():
            self.Preferences = dlg.value(self.Preferences)
            with open(conf_dir+"pychemqtrc", "w") as conf_file:
                self.Preferences.write(conf_file)
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
                label = f"{label:%4g}"

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
        filename = conf_dir + f"{PsychroState().__class__.__name__}_" \
            + f"{self.inputs.P.value:0.0f}.json"
        if os.path.isfile(filename):
            with open(filename, "r") as archivo:
                data = json.load(archivo)
                self.status.setText(self.tr("Loading cached data..."))
                QtWidgets.QApplication.processEvents()
        else:
            self.progressBar.setVisible(True)
            self.status.setText(self.tr("Calculating data..."))
            QtWidgets.QApplication.processEvents()
            data = PsychroState.calculatePlot(self)
            with open(filename, "w") as file:
                json.dump(data, file, indent=4)

            self.progressBar.setVisible(False)
        self.status.setText(self.tr("Plotting..."))
        QtWidgets.QApplication.processEvents()

        # Saturation line
        t = [Temperature(ti).config() for ti in data["t"]]
        Hs = data["Hs"]
        fmt = formatLine(self.Preferences, "Psychr", "saturation")
        if chart:
            self.plt.plot(t, Hs, **fmt)
        else:
            self.plt.plot(Hs, t, **fmt)

        # Iso dew bulb temperaure lines, vertial lines in normal h-Tdb plot,
        # vertical in mollier diagram
        fmt = formatLine(self.Preferences, "Psychr", "isotdb")
        for i, T in enumerate(t):
            if chart:
                self.plt.plot([T, T], [0, Hs[i]], **fmt)
            else:
                self.plt.plot([0, Hs[i]], [T, T], **fmt)

        # Iso humidity lines, horizontal lines in normal h-Tdb plot, vertical
        # in mollier diagram
        H = data["H"]
        th = data["th"]
        tm = Temperature(self.Preferences.getfloat("Psychr", "isotdbEnd"))
        fmt = formatLine(self.Preferences, "Psychr", "isow")
        for i, H in enumerate(H):
            ts = Temperature(th[i]).config()
            if chart:
                self.plt.plot([ts, tm.config()], [H, H], **fmt)
            else:
                self.plt.plot([H, H], [ts, tm.config()], **fmt)

        # Iso relative humidity lines
        fmt = formatLine(self.Preferences, "Psychr", "isohr")
        for Hr, H0 in list(data["Hr"].items()):
            if chart:
                self.plt.plot(t, H0, **fmt)
                self.drawlabel("isohr", t, H0, Hr, "%")
            else:
                self.plt.plot(H0, t, **fmt)
                self.drawlabel("isohr", H0, t, Hr, "%")

        # Iso wet bulb temperature lines
        fmt = formatLine(self.Preferences, "Psychr", "isotwb")
        for T, (H, Tw) in list(data["Twb"].items()):
            value = Temperature(T).config()
            Tw_conf = [Temperature(Twi).config() for Twi in Tw]
            txt = Temperature.text()
            if chart:
                self.plt.plot(Tw_conf, H, **fmt)
                self.drawlabel("isotwb", Tw_conf, H, value, txt)
            else:
                self.plt.plot(H, Tw_conf, **fmt)
                self.drawlabel("isotwb", H, Tw_conf, value, txt)

        # Isochor lines
        fmt = formatLine(self.Preferences, "Psychr", "isochor")
        for v, (Td, H) in list(data["v"].items()):
            value = SpecificVolume(v).config()
            Td_conf = [Temperature(Tdi).config() for Tdi in Td]
            txt = SpecificVolume.text()
            if chart:
                self.plt.plot(Td_conf, H, **fmt)
                self.drawlabel("isochor", Td_conf, H, value, txt)
            else:
                self.plt.plot(H, Td_conf, **fmt)
                self.drawlabel("isochor", H, Td_conf, value, txt)

        if self.plt.state:
            self.plt.createCrux(self.plt.state, chart)
        self.plt.draw()
        self.status.setText(
            f"{self.tr('Using')} {PsychroState().__class__.__name__[3:]}")

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
        """Update progress bar with new value"""
        self.progressBar.setValue(int(value))
        QtWidgets.QApplication.processEvents()


class Config(QtWidgets.QWidget):
    """Phychrometric chart configuration"""
    lineas = [
        ("IsoTdb", Temperature, translate("Psychrometry", "Iso dry bulb temperature")),
        ("IsoW", float, translate("Psychrometry", "Iso absolute humidity")),
        ("IsoHR", float, translate("Psychrometry", "Iso relative humidity")),
        ("IsoTwb", Temperature, translate("Psychrometry", "Iso wet bulb temperature")),
        ("Isochor", SpecificVolume, translate("Psychrometry", "Isochor"))]

    def __init__(self, config, parent=None):
        """constructor, config optional parameter to input project config"""
        super().__init__(parent)
        lyt = QtWidgets.QGridLayout(self)
        lyt.setContentsMargins(0, 0, 0, 0)
        scroll = QtWidgets.QScrollArea()
        scroll.setFrameStyle(QtWidgets.QFrame.Shape.NoFrame)
        lyt.addWidget(scroll)

        dlg = QtWidgets.QWidget()
        layout = QtWidgets.QGridLayout(dlg)

        groupType = QtWidgets.QGroupBox(self.tr("Chart type"))
        groupLayout = QtWidgets.QVBoxLayout(groupType)
        self.checkASHRAE = QtWidgets.QRadioButton(
            self.tr("ASHRAE Chart, W vs Tdb"))
        groupLayout.addWidget(self.checkASHRAE)
        self.checkMollier = QtWidgets.QRadioButton(self.tr("Mollier Chart ix"))
        groupLayout.addWidget(self.checkMollier)
        layout.addWidget(groupType, 0, 1, 1, 2)

        self.virial = QtWidgets.QCheckBox(
            self.tr("Use virial equation of state"))
        layout.addWidget(self.virial, 1, 1, 1, 2)
        self.coolProp = QtWidgets.QCheckBox(
            self.tr("Use external library coolProp (faster)"))
        self.coolProp.setEnabled(False)
        layout.addWidget(self.coolProp, 2, 2)
        self.refprop = QtWidgets.QCheckBox(
            self.tr("Use external library refprop (fastest)"))
        self.refprop.setEnabled(False)
        layout.addWidget(self.refprop, 3, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 4, 1)

        self.satlineconfig = LineConfig(
            "saturation", self.tr("Saturation Line Style"))
        layout.addWidget(self.satlineconfig, 5, 1, 1, 2)
        self.cruxlineconfig = LineConfig("crux", self.tr("Crux Line Style"))
        layout.addWidget(self.cruxlineconfig, 6, 1, 1, 2)
        group = QtWidgets.QGroupBox(self.tr("Isolines"))
        layout.addWidget(group, 7, 1, 1, 2)
        layoutgroup = QtWidgets.QGridLayout(group)
        self.comboIsolineas = QtWidgets.QComboBox()
        layoutgroup.addWidget(self.comboIsolineas, 1, 1)
        self.Isolineas = QtWidgets.QStackedWidget()
        self.comboIsolineas.currentIndexChanged.connect(
            self.Isolineas.setCurrentIndex)
        layoutgroup.addWidget(self.Isolineas, 2, 1)
        for name, unit, text in self.lineas:
            self.comboIsolineas.addItem(text)
            self.Isolineas.addWidget(Isolinea(unit, name, config, "Psychr"))
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 10, 2)

        scroll.setWidget(dlg)

        if os.environ["CoolProp"]:
            self.virial.toggled.connect(self.coolProp.setEnabled)
        if os.environ["refprop"]:
            self.virial.toggled.connect(self.refprop.setEnabled)

        if config.has_section("Psychr"):
            if config.getboolean("Psychr", 'chart'):
                self.checkASHRAE.setChecked(True)
            else:
                self.checkMollier.setChecked(True)
            self.virial.setChecked(config.getboolean("Psychr", 'virial'))
            self.coolProp.setChecked(config.getboolean("Psychr", 'coolprop'))
            self.refprop.setChecked(config.getboolean("Psychr", 'refprop'))
            self.satlineconfig.setConfig(config, "Psychr")
            self.cruxlineconfig.setConfig(config, "Psychr")

    def value(self, config):
        """Return value for main dialog"""
        if not config.has_section("Psychr"):
            config.add_section("Psychr")

        config.set("Psychr", "chart", str(self.checkASHRAE.isChecked()))
        config.set("Psychr", "virial", str(self.virial.isChecked()))
        config.set("Psychr", "coolprop", str(self.coolProp.isChecked()))
        config.set("Psychr", "refprop", str(self.refprop.isChecked()))
        config = self.satlineconfig.value(config, "Psychr")
        config = self.cruxlineconfig.value(config, "Psychr")

        for indice in range(self.Isolineas.count()):
            config = self.Isolineas.widget(indice).value(config)
        return config


class ConfigDialog(QtWidgets.QDialog):
    """Dialog to config thermal method calculations"""
    def __init__(self, config=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle(self.tr("Psychrometric chart configuration"))
        layout = QtWidgets.QVBoxLayout(self)
        self.widget = Config(config)
        layout.addWidget(self.widget)
        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Cancel
            | QtWidgets.QDialogButtonBox.StandardButton.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox)

    def value(self, config):
        """Function result for wizard"""
        config = self.widget.value(config)
        return config


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    aireHumedo = UI_Psychrometry()
    sys.exit(app.exec())
