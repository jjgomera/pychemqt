#!/usr/bin/python3
# -*- coding: utf-8 -*-
# pylint: disable=protected-access, not-callable

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


.. include:: UI_Tables.rst

The module include all UI related functionality of module
    * :class:`plugin`: Implement common functionality used in menu and dialog
    * :class:`Menu`: QMenu to add to mainwindow mainmenu
    * :class:`Dialog`: QDialog with all meos functionality

Dialogs for configuration:
    * :class:`reference.Ui_ReferenceState`: Dialog to select reference state
    * :class:`reference.Ui_Properties`: Dialog for select and sort shown \
    properties in tables
    * :class:`prefMEOS.ColorMapCombo`: Custom QComboBox to choose a \
    matplotlib colormap
    * :class:`prefMEOS.Isolinea`: Widget to configure isolines for mEoS
    * :class:`prefMEOS.Widget`: mEoS parameter configuration dialog
    * :class:`prefMEOS.Dialog`: Dialog tool for standalone use

Dialogs for fluid selection:
    * :class:`chooseFluid.Ui_ChooseFluid`: Dialog to choose fluid for \
    calculations
    * :class:`chooseFluid.DialogFilterFluid`: Dialog for filter compounds \
    family to show
    * :class:`chooseFluid.Dialog_InfoFluid`: Dialog to show parameter of \
    element with meos
    * :class:`chooseFluid.Widget_MEoS_Data`: Widget to show meos data
    * :class:`chooseFluid.transportDialog`: Dialog for transport and \
    ancillary equations
    * :class:`chooseFluid.Widget_Viscosity_Data`: Widget to show viscosity data
    * :class:`chooseFluid.Widget_Conductivity_Data`: Widget to show thermal \
    conductivity data

Library function for plugin
    * :func:`library.getMethod`: Return the thermo method name to use
    * :func:`library.getClassFluid`: Return the thermo class to calculate
    * :func:`library.calcPoint`: Calculate point state and check state in P-T \
    range of eq
    * :func:`library.get_propiedades`: Get the properties to show in tables
    * :func:`library._getData`: Get values of properties in fluid

Plot functionality:
    * :class:`plot.PlotMEoS`: Plot widget to show meos data as plot
    * :class:`plot.Plot2D`: Dialog for select a special 2D plot
    * :class:`plot.Plot3D`: Dialog for define a 3D plot
    * :class:`plot.EditPlot`: Dialog to edit plot
    * :class:`plot.AddLine`: Dialog to add new isoline to plot
    * :class:`plot.EditAxis`: Dialog to configure axes plot properties
    * :class:`plot.AxisWidget`: Dialog to configure axes plot properties
    * :func:`plot.calcIsoline`: Isoline calculation procedure
    * :func:`plot.get_points`: Get point number to plot lines from Preferences
    * :func:`plot.getLineFormat`: get matplotlib line format from preferences
    * :func:`plot.plotIsoline`: plot isoline procedure
    * :func:`plot.plot2D3D`: general procedure for plotting 2D and 3D
    * :func:`plot._getunitTransform`: Return the axis unit transform function \
    to map data to configurated unit

Table functionality:
    * :class:`table.TablaMEoS`: Tabla subclass to show meos data
    * :class:`table.Ui_Saturation`: Dialog to define a two-phase table
    * :class:`table.Ui_Isoproperty`: Dialog to define a isoproperty table
    * :class:`table.AddPoint`: Dialog to add new point to line2D
    * :func:`table.createTabla`: create TablaMEoS


The calculation library itself is in the `lib.meos <lib.meos.html>`__ module
with the compounds implemented in `lib.mEoS <lib.mEoS.html>`__
'''


from configparser import ConfigParser
from functools import partial
from math import log10
import os

from numpy import arange, append, concatenate, linspace, logspace
from scipy.optimize import fsolve

from lib import config, meos, mEoS, unidades
from lib.thermo import ThermoAdvanced
from UI.widgets import createAction
from tools.qt import QtGui, QtWidgets, QtCore, translate

from tools.UI_Tables.chooseFluid import Ui_ChooseFluid
from tools.UI_Tables.library import N_PROP, KEYS, UNITS
from tools.UI_Tables.library import (getClassFluid, getMethod, calcPoint,
                                     saveProperties)
from tools.UI_Tables.plot import (PlotMEoS, Plot2D, Plot3D, calcIsoline,
                                  get_points, plot2D3D, calcMesh)
from tools.UI_Tables.prefMEOS import Dialog as ConfDialog
from tools.UI_Tables.reference import Ui_ReferenceState, Ui_Properties
from tools.UI_Tables.table import Ui_Saturation, Ui_Isoproperty, createTabla


class plugin():
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
            fTxt = translate("meos", "Fluid")
        if self.config.has_option("MEoS", "reference"):
            refTxt = self.config.get("MEoS", "reference")
        else:
            refTxt = translate("meos", "Reference State")
        propTxt = translate("meos", "Properties")
        confTxt = translate("meos", "Configure")

        return fTxt, refTxt, propTxt, confTxt

    def _menuCalculate(self):
        """QMenu for table actions"""
        menu = QtWidgets.QMenu(translate("meos", "Calculate"), parent=self)
        saturationAction = createAction(
            translate("meos", "Saturation"),
            slot=self.showSaturation, parent=self)
        menu.addAction(saturationAction)
        IsopropertyAction = createAction(
            translate("meos", "Isoproperty"),
            slot=self.showIsoproperty, parent=self)
        menu.addAction(IsopropertyAction)
        menu.addSeparator()
        SpecifyAction = createAction(
            translate("meos", "Specified point"),
            slot=self.addTableSpecified, parent=self)
        menu.addAction(SpecifyAction)
        return menu

    def _menuPlot(self):
        """QMenu for plot actions"""
        menu = QtWidgets.QMenu(translate("meos", "Plot"), parent=self)
        Plot_T_s_Action = createAction(
            translate("meos", "T-s diagram"),
            slot=partial(self.plot, "s", "T"), parent=self)
        menu.addAction(Plot_T_s_Action)
        Plot_T_rho_Action = createAction(
            translate("meos", "T-rho diagram"),
            slot=partial(self.plot, "rho", "T"), parent=self)
        menu.addAction(Plot_T_rho_Action)
        Plot_P_h_Action = createAction(
            translate("meos", "P-h diagram"),
            slot=partial(self.plot, "h", "P"), parent=self)
        menu.addAction(Plot_P_h_Action)
        Plot_P_v_Action = createAction(
            translate("meos", "P-v diagram"),
            slot=partial(self.plot, "v", "P"), parent=self)
        menu.addAction(Plot_P_v_Action)
        Plot_P_T_Action = createAction(
            translate("meos", "P-T diagram"),
            slot=partial(self.plot, "T", "P"), parent=self)
        menu.addAction(Plot_P_T_Action)
        Plot_h_s_Action = createAction(
            translate("meos", "h-s diagram"),
            slot=partial(self.plot, "s", "h"), parent=self)
        menu.addAction(Plot_h_s_Action)
        Plot_v_u_Action = createAction(
            translate("meos", "v-u diagram"),
            slot=partial(self.plot, "u", "v"), parent=self)
        menu.addAction(Plot_v_u_Action)
        Plot2DAction = createAction(
            translate("meos", "Other Plots"), slot=self.plot2D, parent=self)
        menu.addAction(Plot2DAction)
        menu.addSeparator()
        Plot3DAction = createAction(
            translate("meos", "3D Plot"), slot=self.plot3D, parent=self)
        menu.addAction(Plot3DAction)
        return menu

    def showChooseFluid(self):
        """Show dialog to choose/view fluid"""
        dlg = Ui_ChooseFluid(self.config)
        if dlg.exec():
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
        if dlg.exec():
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
                refName, refT, refP = "ASHRAE", 233.15, 101325,
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
        if dlg.exec():
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
        if dlg.exec():
            Config = dlg.value(Config)
            with open(config.conf_dir+"pychemqtrc", "w") as file:
                Config.write(file)

    def showSaturation(self):
        """Show dialog to define input for a two-phase saturation table"""
        method = getMethod()
        index = self.config.getint("MEoS", "fluid")
        dlg = Ui_Saturation(method, index)
        if dlg.exec():
            # Get values
            start = dlg.Inicial.value
            end = dlg.Final.value
            incr = dlg.Incremento.value
            fix = dlg.variableFix.value
            value = arange(start, end, incr)
            if (end-start) % incr == 0:
                value = append(value, end)
            fluid = getClassFluid(method, index)

            fluidos = []
            if dlg.VL.isChecked():
                # Liquid-Gas line
                txt = translate("meos", "Liquid-Gas Line")
                if dlg.VariarTemperatura.isChecked():
                    # Changing temperature
                    for val in value:
                        vconfig = unidades.Temperature(val).str
                        self.parent().statusBar().showMessage(
                            "%s: %s =%s, %s" % (fluid.name, "T", vconfig, txt))
                        QtWidgets.QApplication.processEvents()
                        fluidos.append(fluid._new(T=val, x=0.5))
                elif dlg.VariarPresion.isChecked():
                    # Changing pressure
                    for val in value:
                        vconfig = unidades.Temperature(val).str
                        self.parent().statusBar().showMessage(
                            "%s: %s =%s, %s" % (fluid.name, "P", vconfig, txt))
                        QtWidgets.QApplication.processEvents()
                        fluidos.append(fluid._new(P=val, x=0.5))
                elif dlg.VariarXconT.isChecked():
                    # Changing quality with fixed Temperature
                    fconfig = unidades.Temperature(fix).str
                    for val in value:
                        self.parent().statusBar().showMessage(
                            "%s: T =%s  x = %s, %s" % (
                                fluid.name, fconfig, val, txt))
                        QtWidgets.QApplication.processEvents()
                        fluidos.append(fluid._new(T=fix, x=val))
                elif dlg.VariarXconP.isChecked():
                    # Changing quality with fixed pressure
                    fconfig = unidades.Temperature(fix).str
                    for val in value:
                        self.parent().statusBar().showMessage(
                            "%s: P =%s  x = %s, %s" % (
                                fluid.name, fconfig, val, txt))
                        QtWidgets.QApplication.processEvents()
                        fluidos.append(fluid._new(P=fix, x=val))

            else:
                # Melting and sublimation line, only supported for meos
                # internal method
                if dlg.SL.isChecked():
                    func = fluid._Melting_Pressure
                    txt = translate("meos", "Melting Line")
                elif dlg.SV.isChecked():
                    func = fluid._Sublimation_Pressure
                    txt = translate("meos", "Sublimation Line")

                if dlg.VariarTemperatura.isChecked():
                    for val in value:
                        p = func(val)
                        fluidos.append(fluid._new(T=val, P=p))
                        self.parent().statusBar().showMessage(
                            "%s: %s=%0.2f, %s" % (fluid.name, "T", val, txt))
                        QtWidgets.QApplication.processEvents()
                else:
                    for p in value:
                        T = fsolve(lambda T: p-func(T), fluid.Tt)
                        fluidos.append(fluid._new(T=T, P=p))
                        self.parent().statusBar().showMessage(
                            "%s: %s=%0.2f, %s" % (fluid.name, "P", p, txt))
                        QtWidgets.QApplication.processEvents()

            title = translate("meos", "Table %s: %s changing %s (%s)" % (
                    fluid.name, txt, "T", method))
            self.addTable(fluidos, title)
            self.parent().statusBar().clearMessage()

    def showIsoproperty(self):
        """Show dialog to define input for isoproperty table calculations"""
        dlg = Ui_Isoproperty(self.parent())
        if dlg.exec():
            self.parent().updateStatus(translate("meos",
        "Launch MEoS Isoproperty calculation..."))

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

            method = getMethod()
            index = self.config.getint("MEoS", "fluid")
            fluid = getClassFluid(method, index)

            kwarg = {}
            # Define option parameter for transport method, only available
            # for internal meos method
            if method == "meos":
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
                self.parent().statusBar().showMessage(
                    "%s: %s =%s, %s =%s" % (fluid.name, X, v1conf, Y, v2conf))
                QtWidgets.QApplication.processEvents()
                fluidos.append(fluid._new(**kwarg))
            unitX = dlg.unidades[i].text()
            title = translate("meos", "%s: %s =%s %s changing %s (%s)" % (
                    fluid.name, X, v1conf, unitX, meos.propiedades[j],
                    method.upper()))
            self.addTable(fluidos, title)

    def addTable(self, fluidos, title):
        """Add table with properties to mainwindow
        fluidos: List with fluid instances
        title: Text title for window table"""
        tabla = createTabla(self.config, title, fluidos, self.parent())
        method = getMethod()
        fluid = self.config.getint("MEoS", "fluid")
        tabla.Point = getClassFluid(method, fluid)
        wdg = self.parent().centralWidget().currentWidget().addSubWindow(tabla)
        wdg.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(tabla.icon)))
        self.parent().dirty[self.parent().idTab] = True
        self.parent().saveControl()
        tabla.show()

    def addTableSpecified(self):
        """Add blank table to mainwindow to calculata point data"""
        method = getMethod()
        index = self.config.getint("MEoS", "fluid")
        fluid = getClassFluid(method, index)
        name = fluid.name
        method = getMethod()
        title = "%s: %s (%s)" % (name, translate("meos",
    "Specified state points"), method.upper())
        tabla = createTabla(self.config, title, None, self.parent())
        tabla.Point = fluid
        wdg = self.parent().centralWidget().currentWidget().addSubWindow(tabla)
        wdg.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(tabla.icon)))
        self.parent().dirty[self.parent().idTab] = True
        self.parent().saveControl()
        tabla.show()

    def plot2D(self):
        """Add a generic 2D plot to project"""
        dlg = Plot2D(self.parent())
        if dlg.exec():
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
        if dlg.exec():
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
            mesh = dlg.checkMesh.isChecked()
            typeMesh = dlg.typeMesh.currentIndex()

            self.plot(x, y, z=z, mesh=mesh, typemesh=typeMesh)

    def plot(self, x, y, xscale=None, yscale=None, z="", mesh=False,
             typemesh=0):
        """Create a plot
        x: property for axes x
        y: property for axes y
        xscale: scale for axis x
        yscale: scale for axis y
        z: property for axis z, optional to 3D plot"""
        method = getMethod()
        index = self.config.getint("MEoS", "fluid")
        fluid = getClassFluid(method, index)
        method = getMethod()
        filename = "%s-%s.json" % (method, fluid.name.lower())

        if z:
            title = translate("meos",
        "Plot %s: %s=f(%s,%s)" % (fluid.name, z, y, x))
            dim = 3
        else:
            title = translate("meos", "Plot %s: %s=f(%s)" % (fluid.name, y, x))
            dim = 2
        grafico = PlotMEoS(dim=dim, parent=self.parent(), filename=filename)
        grafico.setWindowTitle(title)
        grafico.x = x
        grafico.y = y
        grafico.z = z

        if UNITS[KEYS.index(x)] != unidades.Dimensionless:
            unitx = UNITS[KEYS.index(x)].magnitudes()[0][0]
            i = self.config.getint("Units", unitx)
            xtxt = "%s, %s" % (x, UNITS[KEYS.index(x)].__text__[i])
        else:
            xtxt = "%s" % x
        grafico.plot.ax.set_xlabel(xtxt)

        if UNITS[KEYS.index(y)] != unidades.Dimensionless:
            unity = UNITS[KEYS.index(y)].magnitudes()[0][0]
            j = self.config.getint("Units", unity)
            ytxt = "%s, %s" % (y, UNITS[KEYS.index(y)].__text__[j])
        else:
            ytxt = "%s" % y
        grafico.plot.ax.set_ylabel(ytxt)

        if z:
            grafico.z = z
            if UNITS[KEYS.index(z)] != unidades.Dimensionless:
                unitz = UNITS[KEYS.index(z)].magnitudes()[0][0]
                k = self.config.getint("Units", unitz)
                ztxt = "%s, %s" % (z, UNITS[KEYS.index(z)].__text__[k])
            else:
                ztxt = "%s" % z
            grafico.plot.ax.set_zlabel(ztxt)

        self.parent().statusBar().showMessage(translate("meos",
    "Loading cached data..."))
        QtWidgets.QApplication.processEvents()
        data = grafico._getData()
        if not data:
            self.parent().progressBar.setValue(0)
            self.parent().progressBar.setVisible(True)
            self.parent().statusBar().showMessage(
                translate("meos", "Calculating data, be patient..."))
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
        self.parent().statusBar().showMessage(translate("meos", "Plotting..."))
        QtWidgets.QApplication.processEvents()
        grafico.config = data["config"]
        grafico.changeStatusThermo(data["config"])

        if z:
            kw = {"z": z,
                  "mesh": mesh,
                  "typemesh": typemesh}
            plot2D3D(grafico, data, config.Preferences, x, y, **kw)
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

        self.parent().centralWidget().currentWidget().addSubWindow(grafico)
        self.parent().dirty[self.parent().idTab] = True
        self.parent().saveControl()
        grafico.show()
        self.parent().statusBar().clearMessage()
        grafico.mouseMove.connect(grafico.updatePosition)

    def calculatePlot(self, fluid):
        """Calculate data for plot
            fluid: class of meos fluid to calculate"""
        data = {}
        points = get_points(config.Preferences)
        method = getMethod()

        # Melting and sublimation line only supported in internal meos method
        if method == "meos":
            # Calculate melting line
            if fluid._melting:
                self.parent().statusBar().showMessage(
                    translate("meos", "Calculating melting line..."))
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
                    data["melting"] = saveProperties(fluidos)

            # Calculate sublimation line
            if fluid._sublimation:
                self.parent().statusBar().showMessage(
                    translate("meos", "Calculating sublimation line..."))
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
                    data["sublimation"] = saveProperties(fluidos)

        # Define the saturation temperature
        T = list(concatenate([linspace(fluid.Tt, 0.9*fluid.Tc, points),
                              linspace(0.9*fluid.Tc, 0.99*fluid.Tc, points),
                              linspace(0.99*fluid.Tc, fluid.Tc, points)]))
        for i in range(2, 0, -1):
            del T[points*i]

        # Calculate saturation
        self.parent().statusBar().showMessage(translate("meos",
    "Calculating Liquid-Vapour saturation line..."))
        for fase in [0, 1]:
            fluidos = []
            for Ti in T:
                print("x = %i" % fase, Ti)
                try:
                    fluidos.append(fluid._new(T=Ti, x=fase))
                except:
                    pass
                self.parent().progressBar.setValue(
                    10+5*fase+5*len(fluidos)/len(T))
                QtWidgets.QApplication.processEvents()

            data["saturation_%i" % fase] = saveProperties(fluidos)

        # Calculate isoquality lines
        data["x"] = {}
        self.parent().statusBar().showMessage(translate("meos",
    "Calculating isoquality lines..."))
        values = self.LineList("Isoquality", config.Preferences)
        for i, value in enumerate(values):
            fluidos = calcIsoline(fluid, self.config,
                                  "T", "x", T, value, 20, i, 20,
                                  len(values), self.parent().progressBar)

            data["x"][value] = saveProperties(fluidos)

        # Get limit equation
        if method == "meos":
            eq = fluid.eq[self.parent().currentConfig.getint("MEoS", "eq")]
            Tmin = eq["Tmin"]
            Tmax = eq["Tmax"]

            Tt = eq.get("Tt", fluid.Tt)
            if Tmin > Tt:
                Lt = fluid._new(T=Tmin, x=0)
            else:
                Lt = fluid._new(T=Tt, x=0)
            Pmin = Lt.P

            Pmax = eq["Pmax"]*1000
        elif method == "coolprop":
            Tmin = fluid.eq["Tmin"]
            Tmax = fluid.eq["Tmax"]
            Pmin = fluid.eq["Pmin"]
            Pmax = fluid.eq["Pmax"]
        elif method == "refprop":
            import refprop
            refprop.setup("def", fluid.name)
            limit = refprop.limitx([1], t=-1)
            Tmin = limit["tmin"]
            try:
                Pmin = fluid(T=fluid.Tt, x=1).P
            except:
                Pmin = 100
            Tmax = limit["tmax"]
            Pmax = limit["pmax"]*1000

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
        self.parent().statusBar().showMessage(translate("meos",
    "Calculating isotherm lines..."))
        QtWidgets.QApplication.processEvents()
        values = self.LineList("Isotherm", config.Preferences, fluid)
        for i, value in enumerate(values):
            fluidos = calcIsoline(fluid, self.config,
                                  "P", "T", P, value, 40, i, 10,
                                  len(values), self.parent().progressBar)

            data["T"][value] = saveProperties(fluidos)

        # Calculate isobar lines
        data["P"] = {}
        self.parent().statusBar().showMessage(translate("meos",
    "Calculating isobar lines..."))
        QtWidgets.QApplication.processEvents()
        values = self.LineList("Isobar", config.Preferences, fluid)
        for i, value in enumerate(values):
            fluidos = calcIsoline(fluid, self.config,
                                  "T", "P", T, value, 50, i, 10,
                                  len(values), self.parent().progressBar)
            data["P"][value] = saveProperties(fluidos)

        # Calculate isochor lines
        data["v"] = {}
        self.parent().statusBar().showMessage(translate("meos",
    "Calculating isochor lines..."))
        QtWidgets.QApplication.processEvents()
        values = self.LineList("Isochor", config.Preferences, fluid)
        for i, value in enumerate(values):
            fluidos = calcIsoline(fluid, self.config,
                                  "T", "v", T, value, 60, i, 10,
                                  len(values), self.parent().progressBar)
            data["v"][value] = saveProperties(fluidos)

        # Calculate isoenthalpic lines
        data["h"] = {}
        self.parent().statusBar().showMessage(translate("meos",
    "Calculating isoenthalpic lines..."))
        QtWidgets.QApplication.processEvents()
        vals = self.LineList("Isoenthalpic", config.Preferences, fluid)
        for i, value in enumerate(vals):
            fluidos = calcIsoline(fluid, self.config,
                                  "P", "h", P, value, 70, i, 10,
                                  len(values), self.parent().progressBar)
            data["h"][value] = saveProperties(fluidos)

        # Calculate isoentropic lines
        data["s"] = {}
        self.parent().statusBar().showMessage(translate("meos",
    "Calculating isoentropic lines..."))
        QtWidgets.QApplication.processEvents()
        values = self.LineList("Isoentropic", config.Preferences, fluid)
        for i, value in enumerate(values):
            fluidos = calcIsoline(fluid, self.config,
                                  "P", "s", P, value, 80, i, 20,
                                  len(values), self.parent().progressBar)
            data["s"][value] = saveProperties(fluidos)

        # Calculate 3D mesh
        self.parent().statusBar().showMessage(translate("meos",
    "Calculating 3D mesh data..."))
        QtWidgets.QApplication.processEvents()
        fluidos = calcMesh(fluid, self.config, T, P)
        print(len(fluidos), len(fluidos[0]))
        data["mesh"] = saveProperties(fluidos)

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
            if start <= end and step:
                t = list(arange(start, end+step, step))
            else:
                t = []

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
        title = self.tr("MEoS properties")
        super().__init__(title, parent)
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
        if not (self.config.has_option("MEoS", "fluid")
                and self.config.has_option("MEoS", "reference")):
            menuCalculate.setEnabled(False)
            menuPlot.setEnabled(False)


# Dialog with all meos functionality, to associate to a button in tools toolbar
class Dialog(QtWidgets.QDialog, plugin):
    """Dialog to choose fluid for meos plugins calculations"""
    def __init__(self, config=None, parent=None):
        super().__init__(parent)
        if config is None:
            config = parent.currentConfig
        self.config = config
        self.setWindowTitle(self.tr("Choose fluid"))
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
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 5, 1)
        menuCalculate = self._menuCalculate()
        calculate = QtWidgets.QPushButton(menuCalculate.title())
        calculate.setMenu(menuCalculate)
        layout.addWidget(calculate, 6, 1)
        menuPlot = self._menuPlot()
        plot = QtWidgets.QPushButton(menuPlot.title())
        plot.setMenu(menuPlot)
        layout.addWidget(plot, 6, 2)

        # Disable calculation action if fluid and reference are not defined
        if not (self.config.has_option("MEoS", "fluid")
                and self.config.has_option("MEoS", "reference")):
            calculate.setEnabled(False)
            plot.setEnabled(False)

        layout.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 7, 1, 1, 3)
        btBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Close)
        btBox.clicked.connect(self.reject)
        layout.addWidget(btBox, 8, 1, 1, 3)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)

#     confi = config.getMainWindowConfig()

    SteamTables = Plot3D()

    SteamTables.show()
    sys.exit(app.exec())
