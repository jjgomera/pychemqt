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
# Library with meos plugin functionality
#
#   plugin: Implement meos functionality to common use in menu and dialog
#   Menu: QMenu to add to mainwindow mainmenu with all meos addon functionality
#   Dialog: QDialog with all meos functionality
###############################################################################


from configparser import ConfigParser
from functools import partial
from math import log10
import os

from PyQt5 import QtGui, QtWidgets
from numpy import arange, append, concatenate, linspace, logspace
from scipy.optimize import fsolve

from lib import meos, mEoS, unidades, config
from lib.thermo import ThermoAdvanced
from UI.prefMEOS import Dialog as ConfDialog
from UI.widgets import createAction

from .chooseFluid import Ui_ChooseFluid
from .reference import Ui_ReferenceState, Ui_Properties
from .plot import PlotMEoS, Plot2D, Plot3D, calcIsoline, get_points, plot2D3D
from .table import Ui_Saturation, Ui_Isoproperty, createTabla
from .lib import getClassFluid, getMethod, calcPoint, N_PROP, KEYS, UNITS


class plugin(object):
    """Common functionality to add to menu and dialog in main window"""

    def _txt(self):
        """Common widget names
        fTxt: Fluid name, dynamic by configuration
        refTxt: Reference state name, dynamic by configuration
        propTxt: Properties option name, fixed
        confTxt: Configure option name, fixed
        """
        if self.config.has_option("MEoS", "fluid"):
            fTxt = mEoS.__all__[self.config.getint("MEoS", "fluid")].name
        else:
            fTxt = QtWidgets.QApplication.translate("pychemqt", "Fluid")
        if self.config.has_option("MEoS", "reference"):
            refTxt = self.config.get("MEoS", "reference")
        else:
            refTxt = QtWidgets.QApplication.translate(
                "pychemqt", "Reference State")
        propTxt = QtWidgets.QApplication.translate("pychemqt", "Properties")
        confTxt = QtWidgets.QApplication.translate("pychemqt", "Configure")

        return fTxt, refTxt, propTxt, confTxt

    def _menuCalculate(self):
        """QMenu for table actions"""
        menu = QtWidgets.QMenu(QtWidgets.QApplication.translate(
            "pychemqt", "Calculate"), parent=self)
        saturationAction = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Saturation"),
            slot=self.showSaturation, parent=self)
        menu.addAction(saturationAction)
        IsopropertyAction = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Isoproperty"),
            slot=self.showIsoproperty, parent=self)
        menu.addAction(IsopropertyAction)
        menu.addSeparator()
        SpecifyAction = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Specified point"),
            slot=self.addTableSpecified, parent=self)
        menu.addAction(SpecifyAction)
        return menu

    def _menuPlot(self):
        """QMenu for plot actions"""
        menu = QtWidgets.QMenu(
            QtWidgets.QApplication.translate("pychemqt", "Plot"), parent=self)
        Plot_T_s_Action = createAction(
            QtWidgets.QApplication.translate("pychemqt", "T-s diagram"),
            slot=partial(self.plot, "s", "T"), parent=self)
        menu.addAction(Plot_T_s_Action)
        Plot_T_rho_Action = createAction(
            QtWidgets.QApplication.translate("pychemqt", "T-rho diagram"),
            slot=partial(self.plot, "rho", "T"), parent=self)
        menu.addAction(Plot_T_rho_Action)
        Plot_P_h_Action = createAction(
            QtWidgets.QApplication.translate("pychemqt", "P-h diagram"),
            slot=partial(self.plot, "h", "P"), parent=self)
        menu.addAction(Plot_P_h_Action)
        Plot_P_v_Action = createAction(
            QtWidgets.QApplication.translate("pychemqt", "P-v diagram"),
            slot=partial(self.plot, "v", "P"), parent=self)
        menu.addAction(Plot_P_v_Action)
        Plot_P_T_Action = createAction(
            QtWidgets.QApplication.translate("pychemqt", "P-T diagram"),
            slot=partial(self.plot, "T", "P"), parent=self)
        menu.addAction(Plot_P_T_Action)
        Plot_h_s_Action = createAction(
            QtWidgets.QApplication.translate("pychemqt", "h-s diagram"),
            slot=partial(self.plot, "s", "h"), parent=self)
        menu.addAction(Plot_h_s_Action)
        Plot_v_u_Action = createAction(
            QtWidgets.QApplication.translate("pychemqt", "v-u diagram"),
            slot=partial(self.plot, "u", "v"), parent=self)
        menu.addAction(Plot_v_u_Action)
        Plot2DAction = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Other Plots"),
            slot=self.plot2D, parent=self)
        menu.addAction(Plot2DAction)
        menu.addSeparator()
        Plot3DAction = createAction(
            QtWidgets.QApplication.translate("pychemqt", "3D Plot"),
            slot=self.plot3D, parent=self)
        menu.addAction(Plot3DAction)
        return menu

    def showChooseFluid(self):
        """Show dialog to choose/view fluid"""
        dlg = Ui_ChooseFluid(self.config)
        if dlg.exec_():
            # Update configuration
            if not self.config.has_section("MEoS"):
                self.config.add_section("MEoS")
            self.config.set("MEoS", "fluid", str(dlg.id()))
            self.config.set("MEoS", "eq", str(dlg.eq.currentIndex()))
            self.config.set("MEoS", "PR", str(dlg.radioPR.isChecked()))
            self.config.set("MEoS", "Generalized",
                            str(dlg.generalized.isChecked()))
            self.config.set("MEoS", "visco", str(dlg.visco.currentIndex()))
            self.config.set("MEoS", "thermal", str(dlg.thermal.currentIndex()))
            self.checkProperties()
            self.parent().dirty[self.parent().idTab] = True
            self.parent().saveControl()

            # Update button text in dialog case
            if self.__class__.__name__ == "Dialog":
                fTxt = mEoS.__all__[dlg.lista.currentRow()].name
                self.fluido.setText(fTxt)

    def showReference(self):
        """Show dialog to choose reference state,
        use for enthalpy and entropy zero state
        Don't implemented yet"""
        dlg = Ui_ReferenceState(self.config)
        if dlg.exec_():
            # Get values
            if not self.config.has_section("MEoS"):
                self.config.add_section("MEoS")
            if dlg.OTO.isChecked():
                refName, refT, refP, refH, refS = "OTO", 298.15, 101325, 0, 0
            elif dlg.NBP.isChecked():
                Tb = mEoS.__all__[self.config.getint("MEoS", "fluid")].Tb
                refName, refT, refP, refH, refS = "NBP", Tb, 101325, 0, 0
            elif dlg.IIR.isChecked():
                refName, refT, refP, refH, refS = "IIR", 273.15, 101325, 200, 1
            elif dlg.ASHRAE.isChecked():
                refName, refT, refP, refH, refS = "ASHRAE", 233.15, 101325,
                refH, refS = 0, 0
            else:
                refName = "Custom"
                refT = dlg.T.value
                refP = dlg.P.value
                refH = dlg.h.value
                refS = dlg.s.value

            # Update configuration
            self.config.set("MEoS", "reference", refName)
            self.config.set("MEoS", "Tref", str(refT))
            self.config.set("MEoS", "Pref", str(refP))
            self.config.set("MEoS", "ho", str(refH))
            self.config.set("MEoS", "so", str(refS))
            self.checkProperties()
            self.parent().dirty[self.parent().idTab] = True
            self.parent().saveControl()

            # Update button text in dialog case
            if self.__class__.__name__ == "Dialog":
                self.reference.setText(refName)

    def checkProperties(self):
        """Add default properties to configuration automatic when choose
        fluid or reference state and properties are not defined"""
        if not self.config.has_option("MEoS", "properties"):
            self.config.set("MEoS", "properties", str(Ui_Properties._default))
            self.config.set("MEoS", "phase", "0")
            self.config.set("MEoS", "propertiesOrder",
                            str(list(range(N_PROP))))

    def showProperties(self):
        """Show dialog to choose/sort properties to show in tables"""
        dlg = Ui_Properties(self.config)
        if dlg.exec_():
            # Update configuration
            if not self.config.has_section("MEoS"):
                self.config.add_section("MEoS")
            self.config.set("MEoS", "properties", str(dlg.properties()))
            self.config.set("MEoS", "phase", str(dlg.checkFase.isChecked()))
            self.config.set("MEoS", "propertiesOrder", str(dlg.order))
            self.parent().dirty[self.parent().idTab] = True
            self.parent().saveControl()

    def configure(self):
        """Direct access to configuration"""
        Config = ConfigParser()
        Config.read(config.conf_dir + "pychemqtrc")
        dlg = ConfDialog(Config)
        if dlg.exec_():
            Config = dlg.value(Config)
            Config.write(open(config.conf_dir+"pychemqtrc", "w"))

    def showSaturation(self):
        """Show dialog to define input for a two-phase saturation table"""
        dlg = Ui_Saturation(self.config)
        if dlg.exec_():
            # Get values
            start = dlg.Inicial.value
            end = dlg.Final.value
            incr = dlg.Incremento.value
            fix = dlg.variableFix.value
            value = arange(start, end, incr)
            if (end-start) % incr == 0:
                value = append(value, end)
            fluid = getClassFluid(self.config)
            method = getMethod()

            fluidos = []
            if dlg.VL.isChecked():
                # Liquid-Gas line
                txt = QtWidgets.QApplication.translate(
                    "pychemqt", "Liquid-Gas Line")
                if dlg.VariarTemperatura.isChecked():
                    # Changing temperature
                    for val in value:
                        vconfig = unidades.Temperature(val).str
                        self.parent().statusbar.showMessage(
                            "%s: %s =%s, %s" % (fluid.name, "T", vconfig, txt))
                        fluidos.append(fluid._new(T=val, x=0.5))
                elif dlg.VariarPresion.isChecked():
                    # Changing pressure
                    for val in value:
                        vconfig = unidades.Temperature(val).str
                        self.parent().statusbar.showMessage(
                            "%s: %s =%s, %s" % (fluid.name, "P", vconfig, txt))
                        fluidos.append(fluid._new(P=val, x=0.5))
                elif dlg.VariarXconT.isChecked():
                    # Changing quality with fixed Temperature
                    fconfig = unidades.Temperature(fix).str
                    for val in value:
                        self.parent().statusbar.showMessage(
                            "%s: T =%s  x = %s, %s" % (
                                fluid.name, fconfig, val, txt))
                        fluidos.append(fluid._new(T=fix, x=val))
                elif dlg.VariarXconP.isChecked():
                    # Changing quality with fixed pressure
                    fconfig = unidades.Temperature(fix).str
                    for val in value:
                        self.parent().statusbar.showMessage(
                            "%s: P =%s  x = %s, %s" % (
                                fluid.name, fconfig, val, txt))
                        fluidos.append(fluid._new(P=fix, x=val))

            else:
                # Melting and sublimation line, only supported for meos
                # internal method
                if dlg.SL.isChecked():
                    func = fluid._Melting_Pressure
                    txt = QtWidgets.QApplication.translate(
                        "pychemqt", "Melting Line")
                elif dlg.SV.isChecked():
                    func = fluid._Sublimation_Pressure
                    txt = QtWidgets.QApplication.translate(
                        "pychemqt", "Sublimation Line")

                if dlg.VariarTemperatura.isChecked():
                    for val in value:
                        p = func(val)
                        fluidos.append(fluid._new(T=val, P=p))
                        self.parent().statusbar.showMessage(
                            "%s: %s=%0.2f, %s" % (fluid.name, "T", val, txt))
                else:
                    for p in value:
                        T = fsolve(lambda T: p-func(T), fluid.Tt)
                        fluidos.append(fluid._new(T=T, P=p))
                        self.parent().statusbar.showMessage(
                            "%s: %s=%0.2f, %s" % (fluid.name, "P", p, txt))

            title = QtWidgets.QApplication.translate(
                "pychemqt", "Table %s: %s changing %s (%s)" % (
                    fluid.name, txt, "T", method))
            self.addTable(fluidos, title)
            self.parent().statusbar.clearMessage()

    def showIsoproperty(self):
        """Show dialog to define input for isoproperty table calculations"""
        dlg = Ui_Isoproperty(self.parent())
        if dlg.exec_():
            self.parent().updateStatus(QtWidgets.QApplication.translate(
                "pychemqt", "Launch MEoS Isoproperty calculation..."))

            # Get data from dialog
            i = dlg.fix.currentIndex()
            j = dlg.vary.currentIndex()
            if j >= i:
                j += 1
            X = dlg.keys[i]
            keys = dlg.keys[:]
            Y = keys[j]
            value1 = dlg.variableFix.value
            start = dlg.Inicial.value
            end = dlg.Final.value
            incr = dlg.Incremento.value
            value2 = arange(start, end, incr)
            if (end-start) % incr == 0:
                value2 = append(value2, end)
            v1conf = dlg.unidades[i](value1).str

            fluid = getClassFluid(self.config)
            method = getMethod()

            kwarg = {}
            # Define option parameter for transport method, only available
            # for internal meos method
            if method == "MEOS":
                for key in ("eq", "visco", "thermal"):
                    kwarg[key] = self.config.getint("MEoS", key)

            fluidos = []
            for v2 in value2:
                kwarg[X] = value1
                kwarg[Y] = v2
                if dlg.unidades[j] == float:
                    v2conf = v2
                else:
                    v2conf = dlg.unidades[j](v2).str
                self.parent().statusbar.showMessage(
                    "%s: %s =%s, %s =%s" % (fluid.name, X, v1conf, Y, v2conf))
                fluidos.append(fluid._new(**kwarg))
            unitX = dlg.unidades[i].text()
            title = QtWidgets.QApplication.translate(
                "pychemqt", "%s: %s =%s %s changing %s (%s)" % (
                    fluid.name, X, v1conf, unitX, meos.propiedades[j],
                    method))
            self.addTable(fluidos, title)

    def addTable(self, fluidos, title):
        """Add table with properties to mainwindow
        fluidos: List with fluid instances
        title: Text title for window table"""
        tabla = createTabla(self.config, title, fluidos, self.parent())
        self.parent().centralwidget.currentWidget().addSubWindow(tabla)
        wdg = self.parent().centralwidget.currentWidget().subWindowList()[-1]
        wdg.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(tabla.icon)))
        tabla.show()

    def addTableSpecified(self):
        """Add blank table to mainwindow to calculata point data"""
        fluid = getClassFluid(self.config)
        name = fluid.name
        method = getMethod()
        title = "%s: %s (%s)" % (name, QtWidgets.QApplication.translate(
            "pychemqt", "Specified state points"), method)
        tabla = createTabla(self.config, title, None, self.parent())
        tabla.Point = fluid
        self.parent().centralwidget.currentWidget().addSubWindow(tabla)
        wdg = self.parent().centralwidget.currentWidget().subWindowList()[-1]
        wdg.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(tabla.icon)))
        tabla.show()

    def plot2D(self):
        """Add a generic 2D plot to project"""
        dlg = Plot2D(self.parent())
        if dlg.exec_():
            i = dlg.ejeX.currentIndex()
            j = dlg.ejeY.currentIndex()
            if j >= i:
                j += 1
            prop = ThermoAdvanced.propertiesKey()
            x = prop[i]
            y = prop[j]

            if dlg.Xscale.isChecked():
                xscale = "log"
            else:
                xscale = "linear"
            if dlg.Yscale.isChecked():
                yscale = "log"
            else:
                yscale = "linear"
            self.plot(x, y, xscale, yscale)

    def plot3D(self):
        """Add a generic 3D plot to project"""
        dlg = Plot3D(self.parent())
        if dlg.exec_():
            i = dlg.ejeX.currentIndex()
            j = dlg.ejeY.currentIndex()
            k = dlg.ejeZ.currentIndex()
            if k >= i:
                k += 1
            if k >= j:
                k += 1
            if j >= i:
                j += 1
            prop = ThermoAdvanced.propertiesKey()
            x = prop[i]
            y = prop[j]
            z = prop[k]
            self.plot(x, y, z=z)

    def plot(self, x, y, xscale=None, yscale=None, z=""):
        """Create a plot
        x: property for axes x
        y: property for axes y
        xscale: scale for axis x
        yscale: scale for axis y
        z: property for axis z, optional to 3D plot"""
        fluid = getClassFluid(self.config)
        method = getMethod()
        filename = "%s-%s.pkl" % (method, fluid.name)

        if z:
            title = QtWidgets.QApplication.translate(
                "pychemqt", "Plot %s: %s=f(%s,%s)" % (fluid.name, z, y, x))
            dim = 3
        else:
            title = QtWidgets.QApplication.translate(
                "pychemqt", "Plot %s: %s=f(%s)" % (fluid.name, y, x))
            dim = 2
        grafico = PlotMEoS(dim=dim, parent=self.parent(), filename=filename)
        grafico.setWindowTitle(title)
        grafico.x = x
        grafico.y = y
        grafico.z = z

        unitx = UNITS[KEYS.index(x)].magnitudes()[0][0]
        unity = UNITS[KEYS.index(y)].magnitudes()[0][0]
        i = self.config.getint("Units", unitx)
        j = self.config.getint("Units", unity)
        xtxt = "%s, %s" % (x, UNITS[KEYS.index(x)].__text__[i])
        ytxt = "%s, %s" % (y, UNITS[KEYS.index(y)].__text__[j])
        grafico.plot.ax.set_xlabel(xtxt)
        grafico.plot.ax.set_ylabel(ytxt)
        if z:
            grafico.z = z
            unitz = UNITS[KEYS.index(z)].magnitudes()[0][0]
            k = self.config.getint("Units", unitz)
            ztxt = "%s, %s" % (z, UNITS[KEYS.index(z)].__text__[k])
            grafico.plot.ax.set_zlabel(ztxt)

        self.parent().statusbar.showMessage(QtWidgets.QApplication.translate(
            "pychemqt", "Loading cached data..."))
        QtWidgets.QApplication.processEvents()
        data = grafico._getData()
        if not data:
            self.parent().progressBar.setValue(0)
            self.parent().progressBar.setVisible(True)
            self.parent().statusbar.showMessage(
                QtWidgets.QApplication.translate(
                    "pychemqt", "Calculating data, be patient..."))
            QtWidgets.QApplication.processEvents()
            data = self.calculatePlot(fluid)
            conf = {}
            conf["method"] = method
            conf["fluid"] = self.config.getint("MEoS", "fluid")
            conf["eq"] = self.config.getint("MEoS", "eq")
            conf["visco"] = self.config.getint("MEoS", "visco")
            conf["thermal"] = self.config.getint("MEoS", "thermal")
            data["config"] = conf
            grafico._saveData(data)
            self.parent().progressBar.setVisible(False)
        self.parent().statusbar.showMessage(
            QtWidgets.QApplication.translate("pychemqt", "Plotting..."))
        QtWidgets.QApplication.processEvents()
        grafico.config = data["config"]

        if z:
            plot2D3D(grafico, data, config.Preferences, x, y, z)
        else:
            plot2D3D(grafico, data, config.Preferences, x, y)

            if not xscale:
                if x in ["P", "rho", "v"]:
                    xscale = "log"
                else:
                    xscale = "linear"
            grafico.plot.ax.set_xscale(xscale)
            if not yscale:
                if y in ["P", "rho", "v"]:
                    yscale = "log"
                else:
                    yscale = "linear"
            grafico.plot.ax.set_yscale(yscale)

        grid = config.Preferences.getboolean("MEOS", "grid")
        grafico.plot.ax._gridOn = grid
        grafico.plot.ax.grid(grid)

        self.parent().centralwidget.currentWidget().addSubWindow(grafico)
        grafico.show()
        self.parent().statusbar.clearMessage()

    def calculatePlot(self, fluid):
        """Calculate data for plot
            fluid: class of meos fluid to calculate"""
        data = {}
        points = get_points(config.Preferences)
        method = getMethod()

        # Melting and sublimation line only supported in internal meos method
        if method == "MEOS":
            # Calculate melting line
            if fluid._melting:
                self.parent().statusbar.showMessage(
                    QtWidgets.QApplication.translate(
                        "pychemqt", "Calculating melting line..."))
                T = linspace(fluid._melting["Tmin"], fluid._melting["Tmax"],
                             points)
                fluidos = []
                for Ti in T:
                    P = fluid._Melting_Pressure(Ti)
                    fluido = calcPoint(fluid, self.config, T=Ti, P=P)
                    if fluido:
                        fluidos.append(fluido)
                    self.parent().progressBar.setValue(5*len(fluidos)/len(T))
                    QtWidgets.QApplication.processEvents()
                if fluidos:
                    data["melting"] = {}
                    for x in ThermoAdvanced.propertiesKey():
                        dat_propiedad = []
                        for fluido in fluidos:
                            num = fluido.__getattribute__(x)
                            if num is not None:
                                if x in ["fi", "f"]:
                                    num = num[0]
                                dat_propiedad.append(num._data)
                            else:
                                dat_propiedad.append(None)
                        data["melting"][x] = dat_propiedad

            # Calculate sublimation line
            if fluid._sublimation:
                self.parent().statusbar.showMessage(
                    QtWidgets.QApplication.translate(
                        "pychemqt", "Calculating sublimation line..."))
                T = linspace(fluid._sublimation["Tmin"],
                             fluid._sublimation["Tmax"], points)
                fluidos = []
                for Ti in T:
                    P = fluid._Sublimation_Pressure(Ti)
                    fluido = calcPoint(fluid, self.config, T=Ti, P=P)
                    if fluido:
                        fluidos.append(fluido)
                    self.parent().progressBar.setValue(5+5*len(fluidos)/len(T))
                    QtWidgets.QApplication.processEvents()
                if fluidos:
                    data["sublimation"] = {}
                    for x in ThermoAdvanced.propertiesKey():
                        dat_propiedad = []
                        for fluido in fluidos:
                            num = fluido.__getattribute__(x)
                            if num is not None:
                                if x in ["fi", "f"]:
                                    num = num[0]
                                dat_propiedad.append(num._data)
                            else:
                                dat_propiedad.append(None)
                        data["sublimation"][x] = dat_propiedad

        # Define the saturation temperature
        T = list(concatenate([linspace(fluid.Tt, 0.9*fluid.Tc, points),
                              linspace(0.9*fluid.Tc, 0.99*fluid.Tc, points),
                              linspace(0.99*fluid.Tc, fluid.Tc, points)]))
        for i in range(2, 0, -1):
            del T[points*i]

        # Calculate saturation
        self.parent().statusbar.showMessage(QtWidgets.QApplication.translate(
            "pychemqt", "Calculating Liquid-Vapour saturation line..."))
        for fase in [0, 1]:
            fluidos = []
            for Ti in T:
                fluidos.append(fluid._new(T=Ti, x=fase))
                self.parent().progressBar.setValue(
                    10+5*fase+5*len(fluidos)/len(T))
                QtWidgets.QApplication.processEvents()

            data["saturation_%i" % fase] = {}
            for key in ThermoAdvanced.propertiesKey():
                dat_propiedad = []
                for fluido in fluidos:
                    if fluido.status:
                        p = fluido.__getattribute__(key)
                        if key in ["fi", "f"]:
                            p = p[0]
                        dat_propiedad.append(p)
                data["saturation_%i" % fase][key] = dat_propiedad

        # Calculate isoquality lines
        data["x"] = {}
        self.parent().statusbar.showMessage(QtWidgets.QApplication.translate(
            "pychemqt", "Calculating isoquality lines..."))
        values = self.LineList("Isoquality", config.Preferences)
        for i, value in enumerate(values):
            fluidos = calcIsoline(fluid, self.config,
                                  "T", "x", T, value, 20, i, 20,
                                  len(values), self.parent().progressBar)

            data["x"][value] = {}
            for x in ThermoAdvanced.propertiesKey():
                dat_propiedad = []
                for fluido in fluidos:
                    if fluido is not None and fluido.status:
                        p = fluido.__getattribute__(key)
                        if key in ["fi", "f"]:
                            p = p[0]
                        dat_propiedad.append(p)
                    else:
                        dat_propiedad.append(None)
                data["x"][value][x] = dat_propiedad

        # Get limit equation
        if method == "MEOS":
            eq = fluid.eq[self.parent().currentConfig.getint("MEoS", "eq")]
            Tmin = eq["Tmin"]
            Tmax = eq["Tmax"]

            Tt = eq.get("Tt", fluid.Tt)
            if Tmin > Tt:
                Lt = fluid(T=Tmin, x=0)
            else:
                Lt = fluid(T=Tt, x=0)
            Pmin = Lt.P

            Pmax = eq["Pmax"]*1000
        elif method == "COOLPROP":
            Tmin = fluid.eq["Tmin"]
            Tmax = fluid.eq["Tmax"]
            Pmin = fluid.eq["Pmin"]
            Pmax = fluid.eq["Pmax"]
        elif method == "REFPROP":
            pass

        T = list(concatenate(
            [linspace(Tmin, 0.9*fluid.Tc, points),
             linspace(0.9*fluid.Tc, 0.99*fluid.Tc, points),
             linspace(0.99*fluid.Tc, fluid.Tc, points),
             linspace(fluid.Tc, 1.01*fluid.Tc, points),
             linspace(1.01*fluid.Tc, 1.1*fluid.Tc, points),
             linspace(1.1*fluid.Tc, Tmax, points)]))
        P = list(concatenate(
            [logspace(log10(Pmin), log10(0.9*fluid.Pc), points),
             linspace(0.9*fluid.Pc, 0.99*fluid.Pc, points),
             linspace(0.99*fluid.Pc, fluid.Pc, points),
             linspace(fluid.Pc, 1.01*fluid.Pc, points),
             linspace(1.01*fluid.Pc, 1.1*fluid.Pc, points),
             logspace(log10(1.1*fluid.Pc), log10(Pmax), points)]))
        for i in range(5, 0, -1):
            del T[points*i]
            del P[points*i]

        # Calculate isotherm lines
        data["T"] = {}
        self.parent().statusbar.showMessage(QtWidgets.QApplication.translate(
            "pychemqt", "Calculating isotherm lines..."))
        values = self.LineList("Isotherm", config.Preferences, fluid)
        for i, value in enumerate(values):
            fluidos = calcIsoline(fluid, self.config,
                                  "P", "T", P, value, 40, i, 10,
                                  len(values), self.parent().progressBar)
            data["T"][value] = {}
            for key in ThermoAdvanced.propertiesKey():
                dat_propiedad = []
                for fluido in fluidos:
                    if fluido is not None and fluido.status:
                        p = fluido.__getattribute__(key)
                        if key in ["fi", "f"]:
                            p = p[0]
                        dat_propiedad.append(p)
                    else:
                        dat_propiedad.append(None)
                data["T"][value][key] = dat_propiedad

        # Calculate isobar lines
        data["P"] = {}
        self.parent().statusbar.showMessage(QtWidgets.QApplication.translate(
            "pychemqt", "Calculating isobar lines..."))
        values = self.LineList("Isobar", config.Preferences, fluid)
        for i, value in enumerate(values):
            fluidos = calcIsoline(fluid, self.config,
                                  "T", "P", T, value, 50, i, 10,
                                  len(values), self.parent().progressBar)
            data["P"][value] = {}
            for key in ThermoAdvanced.propertiesKey():
                dat_propiedad = []
                for fluido in fluidos:
                    if fluido is not None and fluido.status:
                        p = fluido.__getattribute__(key)
                        if key in ["fi", "f"]:
                            p = p[0]
                        dat_propiedad.append(p)
                    else:
                        dat_propiedad.append(None)
                data["P"][value][key] = dat_propiedad

        # Calculate isochor lines
        data["v"] = {}
        self.parent().statusbar.showMessage(QtWidgets.QApplication.translate(
            "pychemqt", "Calculating isochor lines..."))
        values = self.LineList("Isochor", config.Preferences, fluid)
        for i, value in enumerate(values):
            fluidos = calcIsoline(fluid, self.config,
                                  "T", "v", T, value, 60, i, 10,
                                  len(values), self.parent().progressBar)
            data["v"][value] = {}
            for key in ThermoAdvanced.propertiesKey():
                dat_propiedad = []
                for fluido in fluidos:
                    if fluido is not None and fluido.status:
                        p = fluido.__getattribute__(key)
                        if key in ["fi", "f"]:
                            p = p[0]
                        dat_propiedad.append(p)
                    else:
                        dat_propiedad.append(None)
                data["v"][value][key] = dat_propiedad

        # Calculate isoenthalpic lines
        data["h"] = {}
        self.parent().statusbar.showMessage(QtWidgets.QApplication.translate(
            "pychemqt", "Calculating isoenthalpic lines..."))
        vals = self.LineList("Isoenthalpic", config.Preferences, fluid)
        for i, value in enumerate(vals):
            fluidos = calcIsoline(fluid, self.config,
                                  "P", "h", P, value, 70, i, 10,
                                  len(values), self.parent().progressBar)
            data["h"][value] = {}
            for key in ThermoAdvanced.propertiesKey():
                dat_propiedad = []
                for fluido in fluidos:
                    if fluido is not None and fluido.status:
                        p = fluido.__getattribute__(key)
                        if key in ["fi", "f"]:
                            p = p[0]
                        dat_propiedad.append(p)
                    else:
                        dat_propiedad.append(None)
                data["h"][value][key] = dat_propiedad

        # Calculate isoentropic lines
        data["s"] = {}
        self.parent().statusbar.showMessage(QtWidgets.QApplication.translate(
            "pychemqt", "Calculating isoentropic lines..."))
        values = self.LineList("Isoentropic", config.Preferences, fluid)
        for i, value in enumerate(values):
            fluidos = calcIsoline(fluid, self.config,
                                  "P", "s", P, value, 80, i, 20,
                                  len(values), self.parent().progressBar)
            data["s"][value] = {}
            for key in ThermoAdvanced.propertiesKey():
                dat_propiedad = []
                for fluido in fluidos:
                    if fluido is not None and fluido.status:
                        p = fluido.__getattribute__(key)
                        if key in ["fi", "f"]:
                            p = p[0]
                        dat_propiedad.append(p)
                    else:
                        dat_propiedad.append(None)
                data["s"][value][key] = dat_propiedad

        return data

    @staticmethod
    def LineList(name, Preferences, fluid=None):
        """Return a list with the values of isoline name to plot"""
        if Preferences.getboolean("MEOS", name+"Custom"):
            t = []
            for i in Preferences.get("MEOS", name+'List').split(','):
                if i:
                    t.append(float(i))
        else:
            start = Preferences.getfloat("MEOS", name+"Start")
            end = Preferences.getfloat("MEOS", name+"End")
            step = Preferences.getfloat("MEOS", name+"Step")
            t = list(arange(start, end+step, step))

        if fluid is not None and Preferences.getboolean("MEOS", name+"Critic"):
            if name == "Isotherm":
                t.append(fluid.Tc)
            elif name == "Isobar":
                t.append(fluid.Pc)
            elif name == "Isochor":
                t.append(1./fluid.rhoc)
            else:
                prop = {"Isoenthalpic": "h",
                        "Isoentropic": "s"}
                fc = fluid._new(T=fluid.Tc, rho=fluid.rhoc)
                t.append(fc.__getattribute__(prop[name]))
        return t


# Plugin to import in mainwindow, it implement all meos functionality as QMenu
class Menu(QtWidgets.QMenu, plugin):
    """QMenu to import in mainwindow with all meos addon functionality"""
    def __init__(self, parent=None):
        title = QtWidgets.QApplication.translate("pychemqt", "MEoS properties")
        super(Menu, self).__init__(title, parent)
        self.setIcon(QtGui.QIcon(
            os.path.join(config.IMAGE_PATH, "button", "tables.png")))
        self.aboutToShow.connect(self.aboutToShow_menu)

    def aboutToShow_menu(self):
        """Populate menu, check if fluid and reference state are defined to
        enable/disable calculation/plot option"""
        self.clear()
        self.config = self.parent().currentConfig

        fTxt, refTxt, propTxt, confTxt = self._txt()
        flAction = createAction(fTxt, slot=self.showChooseFluid, parent=self)
        refAction = createAction(refTxt, slot=self.showReference, parent=self)
        pAction = createAction(propTxt, slot=self.showProperties, parent=self)
        confAction = createAction(confTxt, slot=self.configure, parent=self)
        menuCalculate = self._menuCalculate()
        menuPlot = self._menuPlot()
        self.addAction(flAction)
        self.addAction(refAction)
        self.addAction(pAction)
        self.addAction(confAction)
        self.addSeparator()
        self.addAction(menuCalculate.menuAction())
        self.addAction(menuPlot.menuAction())
        self.addSeparator()

        # Disable calculation action if fluid and reference are not defined
        if not (self.config.has_option("MEoS", "fluid") and
                self.config.has_option("MEoS", "reference")):
            menuCalculate.setEnabled(False)
            menuPlot.setEnabled(False)


# Dialog with all meos functionality, to associate to a button in tools toolbar
class Dialog(QtWidgets.QDialog, plugin):
    """Dialog to choose fluid for meos plugins calculations"""
    def __init__(self, config=None, parent=None):
        super(Dialog, self).__init__(parent)
        if config is None:
            config = parent.currentConfig
        self.config = config
        self.setWindowTitle(
            QtWidgets.QApplication.translate("pychemqt", "Choose fluid"))
        layout = QtWidgets.QGridLayout(self)

        fTxt, refTxt, propTxt, confTxt = self._txt()
        fluid = QtWidgets.QPushButton(fTxt)
        fluid.clicked.connect(self.showChooseFluid)
        layout.addWidget(fluid, 1, 1)
        ref = QtWidgets.QPushButton(refTxt)
        ref.clicked.connect(self.showReference)
        layout.addWidget(ref, 2, 1)
        prop = QtWidgets.QPushButton(propTxt)
        prop.clicked.connect(self.showProperties)
        layout.addWidget(prop, 3, 1)
        conf = QtWidgets.QPushButton(confTxt)
        conf.clicked.connect(self.configure)
        layout.addWidget(conf, 4, 1)

        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Fixed,
            QtWidgets.QSizePolicy.Fixed), 5, 1)
        menuCalculate = self._menuCalculate()
        calculate = QtWidgets.QPushButton(menuCalculate.title())
        calculate.setMenu(menuCalculate)
        layout.addWidget(calculate, 6, 1)
        menuPlot = self._menuPlot()
        plot = QtWidgets.QPushButton(menuPlot.title())
        plot.setMenu(menuPlot)
        layout.addWidget(plot, 6, 2)

        # Disable calculation action if fluid and reference are not defined
        if not (self.config.has_option("MEoS", "fluid") and
                self.config.has_option("MEoS", "reference")):
            calculate.setEnabled(False)
            plot.setEnabled(False)

        layout.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 7, 1, 1, 3)
        btBox = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Close)
        btBox.clicked.connect(self.reject)
        layout.addWidget(btBox, 8, 1, 1, 3)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)

    conf = config.getMainWindowConfig()

    # SteamTables = AddPoint(conf)
    # SteamTables=AddLine(None)
    # SteamTables = Dialog(conf)
    SteamTables = Plot3D()

    SteamTables.show()
    sys.exit(app.exec_())
