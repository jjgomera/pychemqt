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
# Psychrometric graphic tools
#   - PsychroPlot: Plot widget for psychrometric chart
#   - PsychroInput: Widget with input for psychrometric state
#   - UI_Psychrometry: Psychrometric chart
###############################################################################


from configparser import ConfigParser
from functools import partial
import os
import pickle

from PyQt5 import QtCore, QtGui, QtWidgets
from scipy import pi, arctan

from lib.psycrometry import PsyState, PsychroState, _Pbar, _height
from lib.config import conf_dir
from lib.plot import mpl
from lib.unidades import (Temperature, Pressure, Length, Mass,
                          SpecificVolume, Enthalpy)
from UI.widgets import Entrada_con_unidades


class PsychroPlot(mpl):
    """
    Plot widget for psychrometric chart
        Add custom margins
        Define a pointer to text state properties, to remove and redraw
    """
    def __init__(self, *args, **kwargs):
        mpl.__init__(self, *args, **kwargs)
        self.notes = []

    def config(self):
        self.ax.set_autoscale_on(False)
        self.ax.set_xlabel("Tdb, " + Temperature.text())
        ylabel = "%s, %s/%s" % (
            QtWidgets.QApplication.translate("pychemqt", "Absolute humidity"),
            Mass.text(), Mass.text())
        self.ax.set_ylabel(ylabel)
        self.ax.yaxis.set_ticks_position("right")
        self.ax.yaxis.set_label_position("right")

        self.lx = self.ax.axhline(color='b')  # the horiz line
        self.ly = self.ax.axvline(color='b')  # the vert line

        tmin = Temperature(274).config()
        tmax = Temperature(329).config()

        self.ax.set_xlim(tmin, tmax)
        self.ax.set_ylim(0, 0.04)

    def showPointData(self, state):
        """Update data of current cursor point in plot annotates"""
        self.clearPointData()

        yi = 0.95
        for key in ("tdb", "tdp", "twb", "HR", "w", "h", "v", "rho"):
            self.notes.append(self.ax.annotate(
                "%s: %s" % (key, state.__getattribute__(key).str), (0.01, yi),
                xycoords='axes fraction', size="small", va="center"))
            yi -= 0.025
        self.draw()

    def clearPointData(self):
        while self.notes:
            anotation = self.notes.pop()
            anotation.remove()
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
        self.plt = PsychroPlot(self, dpi=90)
        self.plt.fig.canvas.mpl_connect('button_press_event', self.click)
        layout.addWidget(self.plt, 1, 3, 2, 2)

        self.inputs = PsychroInput()
        self.inputs.stateChanged.connect(self.createCrux)
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
        btBox.addButton(butonPNG, QtWidgets.QDialogButtonBox.AcceptRole)
        btBox.rejected.connect(self.reject)
        btBox.accepted.connect(self.savePNG)
        layout.addWidget(btBox, 3, 4)

        self.showToolBar(False)
        self.plot()

    def savePNG(self):
        """Save chart image to png file"""
        fmt = "Portable Network Graphics (*.png)"
        fname, ext = QtWidgets.QFileDialog.getSaveFileName(
            self,
            QtWidgets.QApplication.translate("pychemqt", "Save chart to file"),
            "./", fmt)
        if fname and ext == fmt:
            if fname.split(".")[-1] != "png":
                fname += ".png"
            self.plt.fig.savefig(fname, facecolor='#eeeeee')

    def showToolBar(self, checked):
        """Show/Hide left toolbar with additional funcionality"""
        self.inputs.setVisible(checked)
        if checked:
            image = "arrow-left-double.png"
        else:
            image = "arrow-right-double.png"
        self.buttonShowToolbox.setIcon(QtGui.QIcon(
            os.environ["pychemqt"] + os.path.join("images", "button", image)))

    def drawlabel(self, name, Preferences, t, W, label, unit):
        """
        Draw annotation for isolines
            name: name of isoline
            Preferences: Configparse instance of pychemqt preferences
            t: x array of line
            W: y array of line
            label: text value to draw
            unit: text units to draw
        """
        if Preferences.getboolean("Psychr", name+"label"):
            TMIN = Preferences.getfloat("Psychr", "isotdbStart")
            TMAX = Preferences.getfloat("Psychr", "isotdbEnd")
            tmin = Temperature(TMIN).config()
            tmax = Temperature(TMAX).config()
            x = tmax-tmin
            wmin = Preferences.getfloat("Psychr", "isowStart")
            wmax = Preferences.getfloat("Psychr", "isowEnd")
            y = wmax-wmin

            i = 0
            for ti, wi in zip(t, W):
                if tmin <= ti <= tmax and wmin <= wi <= wmax:
                    i += 1
            label = str(label)
            if Preferences.getboolean("Psychr", name+"units"):
                label += unit
            pos = Preferences.getfloat("Psychr", name+"position")
            p = int(i*pos/100-1)
            rot = arctan((W[p]-W[p-1])/y/(t[p]-t[p-1])*x)*360/2/pi
            self.plt.ax.annotate(label, (t[p], W[p]), rotation=rot,
                                 size="small", ha="center", va="center")

    def plot(self):
        """Plot chart"""
        Preferences = ConfigParser()
        Preferences.read(conf_dir+"pychemqtrc")

        self.plt.ax.clear()
        self.plt.config()
        filename = conf_dir+"%s_%i.pkl" % (
            PsychroState().__class__.__name__, self.inputs.P.value)
        if os.path.isfile(filename):
            with open(filename, "rb") as archivo:
                data = pickle.load(archivo)
                self.status.setText(QtWidgets.QApplication.translate(
                    "pychemqt", "Loading cached data..."))
                QtWidgets.QApplication.processEvents()
        else:
            self.progressBar.setVisible(True)
            self.status.setText(QtWidgets.QApplication.translate(
                "pychemqt", "Calculating data..."))
            QtWidgets.QApplication.processEvents()
            data = PsychroState.calculatePlot(self)
            pickle.dump(data, open(filename, "wb"))
            self.progressBar.setVisible(False)
        self.status.setText(
            QtWidgets.QApplication.translate("pychemqt", "Plotting..."))
        QtWidgets.QApplication.processEvents()

        tm = Temperature(Preferences.getfloat("Psychr", "isotdbEnd"))

        t = [Temperature(ti).config() for ti in data["t"]]
        Hs = data["Hs"]
        format = {}
        format["ls"] = Preferences.get("Psychr", "saturationlineStyle")
        format["lw"] = Preferences.getfloat("Psychr", "saturationlineWidth")
        format["color"] = Preferences.get("Psychr", "saturationColor")
        format["marker"] = Preferences.get("Psychr", "saturationmarker")
        format["markersize"] = 3
        self.plt.plot(t, Hs, **format)

        format = {}
        format["ls"] = Preferences.get("Psychr", "isotdblineStyle")
        format["lw"] = Preferences.getfloat("Psychr", "isotdblineWidth")
        format["color"] = Preferences.get("Psychr", "isotdbColor")
        format["marker"] = Preferences.get("Psychr", "isotdbmarker")
        format["markersize"] = 3
        for i, T in enumerate(t):
            self.plt.plot([T, T], [0, Hs[i]], **format)

        H = data["H"]
        th = data["th"]
        format = {}
        format["ls"] = Preferences.get("Psychr", "isowlineStyle")
        format["lw"] = Preferences.getfloat("Psychr", "isowlineWidth")
        format["color"] = Preferences.get("Psychr", "isowColor")
        format["marker"] = Preferences.get("Psychr", "isowmarker")
        format["markersize"] = 3
        for i, H in enumerate(H):
            self.plt.plot([th[i], tm.config()], [H, H], **format)

        format = {}
        format["ls"] = Preferences.get("Psychr", "isohrlineStyle")
        format["lw"] = Preferences.getfloat("Psychr", "isohrlineWidth")
        format["color"] = Preferences.get("Psychr", "isohrColor")
        format["marker"] = Preferences.get("Psychr", "isohrmarker")
        format["markersize"] = 3
        for Hr, H0 in list(data["Hr"].items()):
            self.plt.plot(t, H0, **format)
            self.drawlabel("isohr", Preferences, t, H0, Hr, "%")

        format = {}
        format["ls"] = Preferences.get("Psychr", "isotwblineStyle")
        format["lw"] = Preferences.getfloat("Psychr", "isotwblineWidth")
        format["color"] = Preferences.get("Psychr", "isotwbColor")
        format["marker"] = Preferences.get("Psychr", "isotwbmarker")
        format["markersize"] = 3
        for T, (H, Tw) in list(data["Twb"].items()):
            self.plt.plot(Tw, H, **format)
            value = Temperature(T).config()
            txt = Temperature.text()
            self.drawlabel("isotwb", Preferences, Tw, H, value, txt)

        format = {}
        format["ls"] = Preferences.get("Psychr", "isochorlineStyle")
        format["lw"] = Preferences.getfloat("Psychr", "isochorlineWidth")
        format["color"] = Preferences.get("Psychr", "isochorColor")
        format["marker"] = Preferences.get("Psychr", "isochormarker")
        format["markersize"] = 3
        for v, (Td, H) in list(data["v"].items()):
            self.plt.plot(Td, H, **format)
            value = SpecificVolume(v).config()
            txt = SpecificVolume.text()
            self.drawlabel("isochor", Preferences, Td, H, value, txt)

        self.plt.draw()
        self.status.setText("%s %s" % (
            QtWidgets.QApplication.translate("pychemqt", "Using"),
            PsychroState().__class__.__name__[3:]))

    def click(self, event):
        """Update input and graph annotate when mouse click over chart"""
        if event.xdata and event.ydata:
            state = self.createState(event.xdata, event.ydata)
            if state.w <= state.ws:
                self.inputs.setState(state)
                self.createCrux(state)
            else:
                self.plt.clearPointData()

    def createState(self, x, y):
        """Create psychrometric state from click or mouse position"""
        tdb = Temperature(x, "conf")
        punto = PsychroState(P=self.inputs.P.value, tdb=tdb, w=y)
        return punto

    def createCrux(self, state):
        """Update horizontal and vertical lines to show click point"""
        self.plt.lx.set_ydata(state.w)
        self.plt.ly.set_xdata(state.tdb.config())
        self.plt.showPointData(state)

    def setProgressValue(self, value):
        self.progressBar.setValue(value)
        QtWidgets.QApplication.processEvents()


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    aireHumedo = UI_Psychrometry()
    aireHumedo.show()
    sys.exit(app.exec_())
