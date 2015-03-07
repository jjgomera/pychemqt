#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Library with meos plugin functionality
# 
#   - plugin: QMenu to import in mainwindow with all meos addon functionality
#
#   Dialogs for configuration:
#   - Ui_ChooseFluid: Dialog to choose fluid for meos plugins calculations
#   - Ui_ReferenceState: Dialog to select reference state
#   - Dialog_InfoFluid: Dialog to show parameter of element with meos
#       - Widget_MEoS_Data: Widget to show meos data
#   - transportDialog: Dialog to show parameters for transport and ancillary equations
#       - Widget_Viscosity_Data: Widget to show viscosity data
#       - Widget_Conductivity_Data: Widget to show thermal conductivity data
#   - Ui_Properties: Dialog for select and sort shown properties in tables
#
#   Table data:
#   - TablaMEoS: Tabla subclass to show meos data, add context menu options
#   - Ui_Saturation: Dialog to define input for a two-phase table calculation
#   - Ui_Isoproperty: Dialog to define input for isoproperty table calculations
#   - AddPoint: Dialog to add new point to line2D
#
#   Plot data:
#   - PlotMEoS: Plot widget to show meos plot data, add context menu options
#   - Plot2D: Dialog for select a special 2D plot
#   - EditPlot: Dialog to edit plot
#   - AddLine: Dialog to add new isoline to plot
#   - EditAxis: Dialog to configure axes plot properties
###############################################################################

import inspect, csv, os, cPickle, gzip
from functools import partial
from string import maketrans
from math import ceil, floor, log10, atan, pi
from ConfigParser import ConfigParser

from PyQt4 import QtCore, QtGui
from numpy import arange, append, concatenate, meshgrid, zeros, linspace, logspace, max, transpose, delete, insert, log
from matplotlib.lines import Line2D
from matplotlib.font_manager import FontProperties

from lib import meos, mEoS, unidades, plot, iapws, config
from lib.utilities import format2txt, representacion
from UI.widgets import Entrada_con_unidades, ClickableLabel, Tabla, createAction, LineStyleCombo, MarkerCombo, ColorSelector, InputFond, Status
from UI.delegate import CheckEditor
from tools.codeEditor import SimplePythonEditor
from tools.UI_Preferences import NumericFactor


prop_pickle = ["T", "P", "rho", "v", "h", "s", "u", "g", "a", "cv", "cp", "w", 
               "Z", "gamma", "fi", "alfav", "kappa", "alfap", "betap", "betas",
               "joule", "Gruneisen", "virialB", "virialC", "dpdT_rho", 
               "dpdrho_T", "drhodT_P", "drhodP_T", "dhdT_rho", "dhdP_T", 
               "dhdT_P", "dhdrho_T", "dhdP_rho", "kt", "ks", "Kt", "Ks", 
               "mu", "k", "nu", "alfa", "Prandt"]


# Plugin to inport in mainwindow, it implement of meos functionality as QMenu
class plugin(QtGui.QMenu):
    """QMenu to import in mainwindow with all meos addon functionality"""
    def __init__(self, parent=None):
        title = QtGui.QApplication.translate("pychemqt", "MEoS properties")
        super(plugin, self).__init__(title, parent)
        self.aboutToShow.connect(self.aboutToShow_menu)

    def aboutToShow_menu(self):
        """Populate menu, check if fluid and reference state are defined to 
        enable/disable calculation/plot option"""
        self.clear()

        if self.parent().currentConfig.has_option("MEoS", "fluid"):
            fluidTxt=mEoS.__all__[self.parent().currentConfig.getint("MEoS","fluid")].name
        else:
            fluidTxt=QtGui.QApplication.translate("pychemqt", "Fluid")
        if self.parent().currentConfig.has_option("MEoS", "reference"):
            refTxt=self.parent().currentConfig.get("MEoS","reference")
        else:
            refTxt=QtGui.QApplication.translate("pychemqt", "Reference State")
        fluidoAction=createAction(fluidTxt, slot=self.showChooseFluid, parent=self)
        referenciaAction=createAction(refTxt, slot=self.showReference, parent=self)
        propAction=createAction(QtGui.QApplication.translate("pychemqt", "Properties"), slot=self.showProperties, parent=self)

        menuCalculate=QtGui.QMenu(QtGui.QApplication.translate("pychemqt", "Calculate"), parent=self)
        saturationAction = createAction(QtGui.QApplication.translate("pychemqt", "Saturation"), slot=self.showSaturation, parent=self)
        menuCalculate.addAction(saturationAction)
        IsopropertyAction = createAction(QtGui.QApplication.translate("pychemqt", "Isoproperty"), slot=self.showIsoproperty, parent=self)
        menuCalculate.addAction(IsopropertyAction)
        menuCalculate.addSeparator()
        SpecifyAction = createAction(QtGui.QApplication.translate("pychemqt", "Specified point"), slot=self.addTableSpecified, parent=self)
        menuCalculate.addAction(SpecifyAction)

        menuPlot=QtGui.QMenu(QtGui.QApplication.translate("pychemqt", "Plot"), parent=self)
        Plot_T_s_Action = createAction(QtGui.QApplication.translate("pychemqt", "T-s diagram"), slot=partial(self.plot, "s", "T"), parent=self)
        menuPlot.addAction(Plot_T_s_Action)
        Plot_T_rho_Action = createAction(QtGui.QApplication.translate("pychemqt", "T-rho diagram"), slot=partial(self.plot, "rho", "T"), parent=self)
        menuPlot.addAction(Plot_T_rho_Action)
        Plot_P_h_Action = createAction(QtGui.QApplication.translate("pychemqt", "P-h diagram"), slot=partial(self.plot, "h", "P"), parent=self)
        menuPlot.addAction(Plot_P_h_Action)
        Plot_P_v_Action = createAction(QtGui.QApplication.translate("pychemqt", "P-v diagram"), slot=partial(self.plot, "v", "P"), parent=self)
        menuPlot.addAction(Plot_P_v_Action)
        Plot_P_T_Action = createAction(QtGui.QApplication.translate("pychemqt", "P-T diagram"), slot=partial(self.plot, "T", "P"), parent=self)
        menuPlot.addAction(Plot_P_T_Action)
        Plot_h_s_Action = createAction(QtGui.QApplication.translate("pychemqt", "h-s diagram"), slot=partial(self.plot, "s", "h"), parent=self)
        menuPlot.addAction(Plot_h_s_Action)
        Plot_v_u_Action = createAction(QtGui.QApplication.translate("pychemqt", "v-u diagram"), slot=partial(self.plot, "u", "v"), parent=self)
        menuPlot.addAction(Plot_v_u_Action)
        Plot2DAction = createAction(QtGui.QApplication.translate("pychemqt", "Other Plots"), slot=self.plot2D, parent=self)
        menuPlot.addAction(Plot2DAction)
        menuPlot.addSeparator()
        Plot3DAction = createAction(QtGui.QApplication.translate("pychemqt", "3D Plot"), slot=self.plot3D, parent=self)
        menuPlot.addAction(Plot3DAction)

        self.addAction(fluidoAction)
        self.addAction(referenciaAction)
        self.addSeparator()
        self.addAction(propAction)
        self.addSeparator()
        self.addAction(menuCalculate.menuAction())
        self.addAction(menuPlot.menuAction())
        self.addSeparator()

        if not (self.config.has_option("MEoS", "fluid") and
                self.config.has_option("MEoS", "reference")):
            menuCalculate.setEnabled(False)
            menuPlot.setEnabled(False)

    def showChooseFluid(self):
        """Show dialog to choose/view fluid"""
        dlg = Ui_ChooseFluid(self.config)
        if dlg.exec_():
            if not self.config.has_section("MEoS"):
                self.config.add_section("MEoS")
            self.config.set("MEoS", "fluid", str(dlg.lista.currentRow()))
            self.config.set("MEoS", "eq", str(dlg.eq.currentIndex()))
            self.config.set("MEoS", "PR", str(dlg.radioPR.isChecked()))
            self.config.set("MEoS", "Generalized", str(dlg.generalized.isChecked()))
            self.config.set("MEoS", "visco", str(dlg.visco.currentIndex()))
            self.config.set("MEoS", "thermal", str(dlg.thermal.currentIndex()))
            if not self.config.has_option("MEoS", "properties"):
                self.config.set("MEoS", "properties", str(Ui_Properties._default))
                self.config.set("MEoS", "phase", "0")
                self.config.set("MEoS", "propertiesOrder", range(64))
            self.parent().dirty[self.parent().idTab]=True
            self.parent().saveControl()

    def showReference(self):
        """Show dialog to choose reference state, for enthalpy and entropy zero state"""
        dlg=Ui_ReferenceState(self.config)
        if dlg.exec_():
            if not self.config.has_section("MEoS"):
                self.config.add_section("MEoS")
            if dlg.OTO.isChecked():
                refName, refT, refP, refH, refS = "OTO", 298.15, 101325, 0, 0
            elif dlg.NBP.isChecked():
                Tb=mEoS.__all__[self.config.getint("MEoS", "fluid")].Tb
                refName, refT, refP, refH, refS = "NBP", Tb, 101325, 0, 0
            elif dlg.IIR.isChecked():
                refName, refT, refP, refH, refS = "IIR", 273.15, 101325, 200, 1
            elif dlg.ASHRAE.isChecked():
                refName, refT, refP, refH, refS = "ASHRAE", 233.15, 101325, 0, 0
            else:
                refName = "Custom"
                refT = dlg.T.value
                refP = dlg.P.value
                refH = dlg.h.value
                refS = dlg.s.value
            self.config.set("MEoS", "reference", refName)
            self.config.set("MEoS", "Tref", str(refT))
            self.config.set("MEoS", "Pref", str(refP))
            self.config.set("MEoS", "ho", str(refH))
            self.config.set("MEoS", "so", str(refS))
            if not self.config.has_option("MEoS", "properties"):
                self.config.set("MEoS", "properties", str(Ui_Properties._default))
                self.config.set("MEoS", "phase", "0")
                self.config.set("MEoS", "propertiesOrder", range(64))
            self.parent().dirty[self.parent().idTab]=True
            self.parent().saveControl()

    def checkProperties(self):
        """Add default properties to show to configuration automatic when
        choose fluid or reference state and properties are not defined"""
        if not self.config.has_option("MEoS", "properties"):
            self.config.set("MEoS", "properties", str(Ui_Properties._default))
            self.config.set("MEoS", "phase", "0")
            self.config.set("MEoS", "propertiesOrder", range(64))

    def showProperties(self):
        """Show dialog to choose/sort properties to show in tables"""
        dlg=Ui_Properties(self.config)
        if dlg.exec_():
            if not self.config.has_section("MEoS"):
                self.config.add_section("MEoS")
            self.config.set("MEoS", "properties", str(dlg.properties))
            self.config.set("MEoS", "phase", str(dlg.checkFase.isChecked()))
            self.config.set("MEoS", "propertiesOrder", str(dlg.order))
            self.parent().dirty[self.parent().idTab]=True
            self.parent().saveControl()

    def showSaturation(self):
        """Show dialog to define input for a two-phase saturation table calculation"""
        dlg=Ui_Saturation(self.config)
        if dlg.exec_():
            start=dlg.Inicial.value
            end=dlg.Final.value
            fix=dlg.variableFix.value
            incr=dlg.Incremento.value
            value=arange(start, end+incr, incr)
            fluid=mEoS.__all__[self.config.getint("MEoS", "fluid")]

            fluidos=[]
            if dlg.VL.isChecked():
                txt=QtGui.QApplication.translate("pychemqt", "Liquid-Gas Line")
                if dlg.VariarTemperatura.isChecked():
                    for val in value:
                        vconfig = unidades.Temperature(val).str
                        self.parent().statusbar.showMessage(
                            "%s: %s =%s, %s" % (fluid.name, "T", vconfig, txt))
                        fluidos.append(fluid(T=val, x=0.5))
                elif dlg.VariarPresion.isChecked():
                    for val in value:
                        vconfig = unidades.Temperature(val).str
                        self.parent().statusbar.showMessage(
                            "%s: %s =%s, %s" % (fluid.name, "P", vconfig, txt))
                        fluidos.append(fluid(P=val, x=0.5))
                elif dlg.VariarXconT.isChecked():
                    fconfig = unidades.Temperature(fix).str
                    for val in value:
                        self.parent().statusbar.showMessage(
                            "%s: T =%s  x = %s, %s" % (fluid.name, fconfig, val, txt))
                        fluidos.append(fluid(T=fix, x=val))
                elif dlg.VariarXconP.isChecked():
                    fconfig = unidades.Temperature(fix).str
                    for val in value:
                        self.parent().statusbar.showMessage(
                            "%s: P =%s  x = %s, %s" % (fluid.name, fconfig, val, txt))
                        fluidos.append(fluid(P=fix, x=val))

            else:
                if dlg.SL.isChecked():
                    func=fluid._Melting_Pressure
                    txt=QtGui.QApplication.translate("pychemqt", "Melting Line")
                elif dlg.SV.isChecked():
                    func=fluid._Sublimation_Pressure
                    txt=QtGui.QApplication.translate("pychemqt", "Sublimation Line")

                if dlg.VariarTemperatura.isChecked():
                    for val in value:
                        p=func(val)
                        fluidos.append(fluid(T=val, P=p))
                        self.parent().statusbar.showMessage(
                            "%s: %s=%0.2f, %s" % (fluid.name, "T", val, txt))
                else:
                    for p in value:
                        T=fsolve(lambda T: p-func(T), fluid.Tt)
                        fluidos.append(fluid(T=T, P=p))
                        self.parent().statusbar.showMessage(
                            "%s: %s=%0.2f, %s" % (fluid.name, "P", p, txt))

            title=QtGui.QApplication.translate(
                "pychemqt", "Table %s: %s changing %s" %(fluid.name, txt, "T"))
            self.addTable(fluidos, title)
            self.parent().statusbar.clearMessage()
            
    def showIsoproperty(self):
        """Show dialog to define input for isoproperty table calculations"""
        dlg=Ui_Isoproperty(self.parent())
        if dlg.exec_():
            self.parent().updateStatus(QtGui.QApplication.translate(
                "pychemqt", "Launch MEoS Isoproperty calculation..."))
            indice1=dlg.fix.currentIndex()
            indice2=dlg.vary.currentIndex()
            if indice2 >= indice1:
                indice2 += 1
            var1=dlg.keys[indice1]
            keys=dlg.keys[:]
            var2=keys[indice2]
            value1=dlg.variableFix.value
            start=dlg.Inicial.value
            end=dlg.Final.value
            incr=dlg.Incremento.value
            value2=arange(start, end, incr)
            if (end-start)%incr == 0:
                value2=append(value2, end)
            fluid=mEoS.__all__[self.config.getint("MEoS", "fluid")]
            kwarg={}
            for key in ("eq", "visco", "thermal"):
                kwarg[key]=self.config.getint("MEoS", key)
            v1config = dlg.unidades[indice1](value1).str
            fluidos=[]
            for v2 in value2:
                kwarg[var1]=value1
                kwarg[var2]=v2
                if dlg.unidades[indice2] == float:
                    v2config = v2
                else:
                    v2config = dlg.unidades[indice2](v2).str
                self.parent().statusbar.showMessage(
                    "%s: %s =%s, %s =%s" % (fluid.name, var1, v1config, var2, v2config))
                fluidos.append(fluid(**kwarg))
            title=QtGui.QApplication.translate("pychemqt", "%s: %s =%s %s changing %s"
                %(fluid.formula, var1, v1config, dlg.unidades[indice1].text(), meos.propiedades[indice2]))
            self.addTable(fluidos, title)

    def addTable(self, fluidos, title):
        """Add table with properties to mainwindow
        fluidos: List with fluid instances
        title: Text title for window table"""
        tabla=createTabla(self.config, title, fluidos, self.parent())
        self.parent().centralwidget.currentWidget().addSubWindow(tabla)
        tabla.show()

    def addTableSpecified(self):
        """Add blank table to mainwindow to calculata point data"""
        fluid = mEoS.__all__[self.config.getint("MEoS", "fluid")]
        name = fluid.formula
        title="%s: %s" % (name, QtGui.QApplication.translate(
            "pychemqt", "Specified state points"))
        tabla=createTabla(self.config, title, None, self.parent())
        tabla.fluid = fluid
        tabla.Point = fluid()
        self.parent().centralwidget.currentWidget().addSubWindow(tabla)
        tabla.show()

    def plot2D(self):
        """Add a generic 2D plot to project"""
        dlg=Plot2D(self.parent())
        if dlg.exec_():
            i=dlg.ejeX.currentIndex()
            j=dlg.ejeY.currentIndex()
            if j>=i:
                j+=1
            x = prop_pickle[i]
            y = prop_pickle[j]

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
        """Add a generic 2D plot to project"""
        dlg=Plot3D(self.parent())
        if dlg.exec_():
            i=dlg.ejeX.currentIndex()
            j=dlg.ejeY.currentIndex()
            k=dlg.ejeZ.currentIndex()
            if k >= i:
                k+=1
            if k >= j:
                k+=1
            if j>=i:
                j+=1
            x = prop_pickle[i]
            y = prop_pickle[j]
            z = prop_pickle[k]
            self.plot(x, y, z=z)

    def plot(self, x, y, xscale=None, yscale=None, z=""):
        """Create a plot
        x: property for axes x
        y: property for axes y
        xscale: scale for axis x
        yscale: scale for axis y
        z: property for axis z, optional to 3D plot"""
        index = self.config.getint("MEoS", "fluid")
        fluid=mEoS.__all__[index]
        filename = "%s.pkl" % (fluid.formula)

        if z:
            title=QtGui.QApplication.translate(
                "pychemqt", "Plot %s: %s=f(%s,%s)" %(fluid.formula, z, y, x))
            dim = 3
        else:
            title=QtGui.QApplication.translate(
                "pychemqt", "Plot %s: %s=f(%s)" %(fluid.formula, y, x))
            dim = 2
        grafico=PlotMEoS(dim=dim, parent=self.parent(), filename=filename)
        grafico.setWindowTitle(title)
        grafico.x = x
        grafico.y = y
        grafico.z = z

        unitx = meos.units[meos.keys.index(x)].magnitudes()[0][0]
        unity = meos.units[meos.keys.index(y)].magnitudes()[0][0]
        i = self.config.getint("Units", unitx)
        j = self.config.getint("Units", unity)
        xtxt = "%s, %s" %(x, meos.units[meos.keys.index(x)].__text__[i])
        ytxt = "%s, %s" %(y, meos.units[meos.keys.index(y)].__text__[j])
        grafico.plot.ax.set_xlabel(xtxt)
        grafico.plot.ax.set_ylabel(ytxt)
        if z:
            grafico.z = z
            unitz = meos.units[meos.keys.index(z)].magnitudes()[0][0]
            k = self.config.getint("Units", unitz)
            ztxt = "%s, %s" %(z, meos.units[meos.keys.index(z)].__text__[k])
            grafico.plot.ax.set_zlabel(ztxt)

        self.parent().statusbar.showMessage(QtGui.QApplication.translate(
            "pychemqt", "Loading cached data..."))
        QtGui.QApplication.processEvents()
        data = grafico._getData()
        if not data:
            self.parent().progressBar.setValue(0)
            self.parent().progressBar.setVisible(True)
            self.parent().statusbar.showMessage(QtGui.QApplication.translate(
                "pychemqt", "Calculating data, be patient..."))
            QtGui.QApplication.processEvents()
            data = self.calculatePlot(fluid)
            conf = {}
            conf["fluid"] = index
            conf["eq"] = self.config.getint("MEoS", "eq")
            conf["visco"] = self.config.getint("MEoS", "visco")
            conf["thermal"] = self.config.getint("MEoS", "thermal")
            data["config"] = conf
            grafico._saveData(data)
            self.parent().progressBar.setVisible(False)
        self.parent().statusbar.showMessage(
            QtGui.QApplication.translate("pychemqt", "Plotting..."))
        QtGui.QApplication.processEvents()
        grafico.config=data["config"]

        if z:
            plot2D3D(grafico, data, self.parent().Preferences, x, y, z)
        else:
            plot2D3D(grafico, data, self.parent().Preferences, x, y)

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

        grid = self.parent().Preferences.getboolean("MEOS", "grid")
        grafico.plot.ax._gridOn = grid
        grafico.plot.ax.grid(grid)

        self.parent().centralwidget.currentWidget().addSubWindow(grafico)
        grafico.show()
        self.parent().statusbar.clearMessage()

    def calculatePlot(self, fluid):
        """Calculate data for plot
            fluid: class of meos fluid to calculate"""
        data = {}
        points = get_points(self.parent().Preferences)

        # Calculate melting line
        if fluid._melting:
            self.parent().statusbar.showMessage(QtGui.QApplication.translate(
                "pychemqt", "Calculating melting line..."))
            T = linspace(fluid._melting["Tmin"], fluid._melting["Tmax"], points)
            fluidos = []
            for Ti in T:
                P = fluid._Melting_Pressure(Ti)
                fluido = calcPoint(fluid, self.config, T=Ti, P=P)
                if fluido:
                    fluidos.append(fluido)
                self.parent().progressBar.setValue(5*len(fluidos)/len(T))
                QtGui.QApplication.processEvents()
            if fluidos:
                data["melting"]={}
                for x in prop_pickle:
                    dat_propiedad=[]
                    for fluido in fluidos:
                        num = fluido.__getattribute__(x)
                        if num is not None:
                            dat_propiedad.append(num._data)
                        else:
                            dat_propiedad.append(None)
                    data["melting"][x]=dat_propiedad

        # Calculate sublimation line
        if fluid._sublimation:
            self.parent().statusbar.showMessage(QtGui.QApplication.translate(
                "pychemqt", "Calculating sublimation line..."))
            T = linspace(fluid._sublimation["Tmin"], fluid._sublimation["Tmax"], points)
            fluidos = []
            for Ti in T:
                P = fluid._Sublimation_Pressure(Ti)
                fluido = calcPoint(fluid, self.config, T=Ti, P=P)
                if fluido:
                    fluidos.append(fluido)
                self.parent().progressBar.setValue(5+5*len(fluidos)/len(T))
                QtGui.QApplication.processEvents()
            if fluidos:
                data["sublimation"]={}
                for x in prop_pickle:
                    dat_propiedad=[]
                    for fluido in fluidos:
                        num = fluido.__getattribute__(x)
                        if num is not None:
                            dat_propiedad.append(num._data)
                        else:
                            dat_propiedad.append(None)
                    data["sublimation"][x]=dat_propiedad


        T = list(concatenate([linspace(fluid.Tt, 0.9*fluid.Tc, points),
                              linspace(0.9*fluid.Tc, 0.99*fluid.Tc, points),
                              linspace(0.99*fluid.Tc, fluid.Tc, points)]))
        for i in range(2, 0, -1):
            del T[points*i]

        # Calculate saturation
        self.parent().statusbar.showMessage(QtGui.QApplication.translate(
            "pychemqt", "Calculating Liquid-Vapour saturation line..."))
        if x == "P" and y == "T":
            fases = [0]
        else:
            fases = [0, 1]
        for fase in fases:
            fluidos = []
            for Ti in T:
                fluidos.append(fluid(T=Ti, x=fase))
                self.parent().progressBar.setValue(10+5*fase+5*len(fluidos)/len(T))
                QtGui.QApplication.processEvents()

            data["saturation_%i" %fase]={}
            for x in prop_pickle:
                dat_propiedad=[]
                for fluido in fluidos:
                    num = fluido.__getattribute__(x)
                    if num is not None:
                        dat_propiedad.append(num._data)
                    else:
                        dat_propiedad.append(None)
                data["saturation_%i" %fase][x]=dat_propiedad

        # Calculate isoquality lines
        data["x"]={}
        self.parent().statusbar.showMessage(QtGui.QApplication.translate(
            "pychemqt", "Calculating isoquality lines..."))
        values = self.LineList("Isoquality", self.parent().Preferences)
        for i, value in enumerate(values):
            fluidos = calcIsoline(fluid, self.config,
                                  "T", "x", T, value, 20, i, 20,
                                  len(values), self.parent().progressBar)

            data["x"][value]={}
            for x in prop_pickle:
                dat_propiedad=[]
                for fluido in fluidos:
                    num = fluido.__getattribute__(x)
                    if num is not None:
                        dat_propiedad.append(num._data)
                    else:
                        dat_propiedad.append(None)
                data["x"][value][x]=dat_propiedad

        eq = fluid.eq[self.parent().currentConfig.getint("MEoS", "eq")]
        T = list(concatenate([linspace(eq["Tmin"], 0.9*fluid.Tc, points),
                              linspace(0.9*fluid.Tc, 0.99*fluid.Tc, points),
                              linspace(0.99*fluid.Tc, fluid.Tc, points),
                              linspace(fluid.Tc, 1.01*fluid.Tc, points),
                              linspace(1.01*fluid.Tc, 1.1*fluid.Tc, points),
                              linspace(1.1*fluid.Tc, eq["Tmax"], points)]))
        Pmin = eq["Pmin"]*1000
        Pmax = eq["Pmax"]*1000
        P = list(concatenate([logspace(log10(Pmin), log10(0.9*fluid.Pc), points),
                              linspace(0.9*fluid.Pc, 0.99*fluid.Pc, points),
                              linspace(0.99*fluid.Pc, fluid.Pc, points),
                              linspace(fluid.Pc, 1.01*fluid.Pc, points),
                              linspace(1.01*fluid.Pc, 1.1*fluid.Pc, points),
                              logspace(log10(1.1*fluid.Pc), log10(Pmax), points)]))
        for i in range(5, 0, -1):
            del T[points*i]
            del P[points*i]

        # Calculate isotherm lines
        data["T"]={}
        self.parent().statusbar.showMessage(QtGui.QApplication.translate(
            "pychemqt", "Calculating isotherm lines..."))
        values = self.LineList("Isotherm", self.parent().Preferences, fluid)
        for i, value in enumerate(values):
            fluidos = calcIsoline(fluid, self.config,
                                  "P", "T", P, value, 40, i, 10,
                                  len(values), self.parent().progressBar)
            data["T"][value]={}
            for x in prop_pickle:
                dat_propiedad=[]
                for fluido in fluidos:
                    num = fluido.__getattribute__(x)
                    if num is not None:
                        dat_propiedad.append(num._data)
                    else:
                        dat_propiedad.append(None)
                data["T"][value][x]=dat_propiedad

        # Calculate isobar lines
        data["P"]={}
        self.parent().statusbar.showMessage(QtGui.QApplication.translate(
            "pychemqt", "Calculating isobar lines..."))
        values = self.LineList("Isobar", self.parent().Preferences, fluid)
        for i, value in enumerate(values):
            fluidos = calcIsoline(fluid, self.config,
                                  "T", "P", T, value, 50, i, 10,
                                  len(values), self.parent().progressBar)
            data["P"][value]={}
            for x in prop_pickle:
                dat_propiedad=[]
                for fluido in fluidos:
                    num = fluido.__getattribute__(x)
                    if num is not None:
                        dat_propiedad.append(num._data)
                    else:
                        dat_propiedad.append(None)
                data["P"][value][x]=dat_propiedad

        # Calculate isochor lines
        data["v"]={}
        self.parent().statusbar.showMessage(QtGui.QApplication.translate(
            "pychemqt", "Calculating isochor lines..."))
        values = self.LineList("Isochor", self.parent().Preferences, fluid)
        for i, value in enumerate(values):
            fluidos = calcIsoline(fluid, self.config,
                                  "T", "v", T, value, 60, i, 10,
                                  len(values), self.parent().progressBar)
            data["v"][value]={}
            for x in prop_pickle:
                dat_propiedad=[]
                for fluido in fluidos:
                    num = fluido.__getattribute__(x)
                    if num is not None:
                        dat_propiedad.append(num._data)
                    else:
                        dat_propiedad.append(None)
                data["v"][value][x]=dat_propiedad

        # Calculate isoenthalpic lines
        data["h"]={}
        self.parent().statusbar.showMessage(QtGui.QApplication.translate(
            "pychemqt", "Calculating isoenthalpic lines..."))
        values = self.LineList("Isoenthalpic", self.parent().Preferences, fluid)
        for i, value in enumerate(values):
            fluidos = calcIsoline(fluid, self.config,
                                  "T", "h", T, value, 70, i, 10,
                                  len(values), self.parent().progressBar)
            data["h"][value]={}
            for x in prop_pickle:
                dat_propiedad=[]
                for fluido in fluidos:
                    num = fluido.__getattribute__(x)
                    if num is not None:
                        dat_propiedad.append(num._data)
                    else:
                        dat_propiedad.append(None)
                data["h"][value][x]=dat_propiedad

        # Calculate isoentropic lines
        data["s"]={}
        self.parent().statusbar.showMessage(QtGui.QApplication.translate(
            "pychemqt", "Calculating isoentropic lines..."))
        values = self.LineList("Isoentropic", self.parent().Preferences, fluid)
        for i, value in enumerate(values):
            fluidos = calcIsoline(fluid, self.config,
                                  "T", "s", T, value, 80, i, 20,
                                  len(values), self.parent().progressBar)
            data["s"][value]={}
            for x in prop_pickle:
                dat_propiedad=[]
                for fluido in fluidos:
                    num = fluido.__getattribute__(x)
                    if num is not None:
                        dat_propiedad.append(num._data)
                    else:
                        dat_propiedad.append(None)
                data["s"][num][x]=dat_propiedad

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
                fc = fluid(T=fluid.Tc, rho=fluid.rhoc)
                t.append(fc.__getattribute__(prop[name]))
        return t


# Dialogs for configuration:
class Ui_ChooseFluid(QtGui.QDialog):
    """Dialog to choose fluid for meos plugins calculations"""
    def __init__(self, config=None, parent=None):
        """config: instance with project config to set initial values"""
        super(Ui_ChooseFluid, self).__init__(parent)
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Choose fluid"))
        layout = QtGui.QGridLayout(self)

        self.lista = QtGui.QListWidget()
        for fluido in mEoS.__all__:
            txt=fluido.name
            if fluido.synonym:
                txt+=" ("+fluido.synonym+")"
            self.lista.addItem(txt)
        self.lista.itemDoubleClicked.connect(self.accept)
        layout.addWidget(self.lista,1,1,3,1)

        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok|QtGui.QDialogButtonBox.Cancel, QtCore.Qt.Vertical)
        botonInfo=QtGui.QPushButton(QtGui.QApplication.translate("pychemqt", "Info"))
        self.buttonBox.addButton(botonInfo, QtGui.QDialogButtonBox.HelpRole)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        self.buttonBox.helpRequested.connect(self.info)
        layout.addWidget(self.buttonBox,1,2,1,1)

        self.widget=QtGui.QWidget(self)
        self.widget.setVisible(False)
        layout.addWidget(self.widget,4,1,1,2)
        gridLayout = QtGui.QGridLayout(self.widget)
        self.radioMEoS=QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "Use MEoS equation"))
        self.radioMEoS.setChecked(True)
        gridLayout.addWidget(self.radioMEoS,1,1,1,2)
        gridLayout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Equation")+": "),2,1)
        self.eq=QtGui.QComboBox()
        gridLayout.addWidget(self.eq,2,2)
        self.generalized=QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "Use generalizated expression"))
        gridLayout.addWidget(self.generalized,3,1,1,2)
        self.radioPR=QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "Use Peng-Robinson cubic equation"))
        gridLayout.addWidget(self.radioPR,4,1,1,2)

        gridLayout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),5,1)
        gridLayout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Viscosity")),6,1)
        self.visco=QtGui.QComboBox()
        gridLayout.addWidget(self.visco,6,2)
        gridLayout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Thermal")),7,1)
        self.thermal=QtGui.QComboBox()
        gridLayout.addWidget(self.thermal,7,2)
        gridLayout.addItem(QtGui.QSpacerItem(0,0,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Maximum),8,2)

        self.botonMore=QtGui.QPushButton(QtGui.QApplication.translate("pychemqt", "More..."))
        self.botonMore.setCheckable(True)
        self.botonMore.clicked.connect(self.widget.setVisible)
        layout.addWidget(self.botonMore,3,2,1,1)

        self.lista.currentRowChanged.connect(self.update)
        self.radioMEoS.toggled.connect(self.eq.setEnabled)

        if config and config.has_option("MEoS", "fluid"):
            self.lista.setCurrentRow(config.getint("MEoS", "fluid"))
            self.eq.setCurrentIndex(config.getint("MEoS", "eq"))
            self.radioPR.setChecked(config.getboolean("MEoS", "PR"))
            self.generalized.setChecked(config.getboolean("MEoS", "Generalized"))
            self.visco.setCurrentIndex(config.getint("MEoS", "visco"))
            self.thermal.setCurrentIndex(config.getint("MEoS", "thermal"))

    def info(self):
        """Show info dialog for fluid"""
        dialog=Dialog_InfoFluid(mEoS.__all__[self.lista.currentRow()])
        dialog.exec_()

    def update(self, indice):
        """Update data when selected fluid change"""
        fluido=mEoS.__all__[indice]
        self.eq.clear()
        for eq in fluido.eq:
            self.eq.addItem(eq["__name__"])

        self.visco.clear()
        if fluido._viscosity is not None:
            self.visco.setEnabled(True)
            for eq in fluido._viscosity:
                self.visco.addItem(eq["__name__"])
        else:
                self.visco.addItem(QtGui.QApplication.translate("pychemqt", "Undefined"))
                self.visco.setEnabled(False)

        self.thermal.clear()
        if fluido._thermal is not None:
            self.thermal.setEnabled(True)
            for eq in fluido._thermal:
                self.thermal.addItem(eq["__name__"])
        else:
            self.thermal.addItem(QtGui.QApplication.translate("pychemqt", "Undefined"))
            self.thermal.setEnabled(False)


class Ui_ReferenceState(QtGui.QDialog):
    """Dialog to select reference state"""
    def __init__(self, config=None, parent=None):
        """config: instance with project config to set initial values"""
        super(Ui_ReferenceState, self).__init__(parent)
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Select reference state"))
        layout = QtGui.QGridLayout(self)
        self.OTO=QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "OTO,  h,s=0 at 25ºC and 1 atm", None, QtGui.QApplication.UnicodeUTF8))
        layout.addWidget(self.OTO,0,1,1,7)
        self.NBP=QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "NBP,  h,s=0 saturated liquid at Tb"))
        layout.addWidget(self.NBP,1,1,1,7)
        self.IIR=QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "IIR,  h=200,s=1 saturated liquid at 0ºC", None, QtGui.QApplication.UnicodeUTF8))
        layout.addWidget(self.IIR,2,1,1,7)
        self.ASHRAE=QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "ASHRAE,  h,s=0 saturated liquid at -40ºC", None, QtGui.QApplication.UnicodeUTF8))
        layout.addWidget(self.ASHRAE,3,1,1,7)
        self.personalizado=QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "Custom"))
        self.personalizado.toggled.connect(self.setEnabled)
        layout.addWidget(self.personalizado,4,1,1,7)

        layout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed), 5,1,1,1)
        layout.addWidget(QtGui.QLabel("T:"),5,2,1,1)
        self.T = Entrada_con_unidades(unidades.Temperature, value=298.15)
        layout.addWidget(self.T,5,3,1,1)
        layout.addWidget(QtGui.QLabel("P:"),6,2,1,1)
        self.P = Entrada_con_unidades(unidades.Pressure, value=101325)
        layout.addWidget(self.P,6,3,1,1)
        layout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Fixed), 5,4,2,1)
        layout.addWidget(QtGui.QLabel("h:"),5,5,1,1)
        self.h = Entrada_con_unidades(unidades.Enthalpy, value=0)
        layout.addWidget(self.h,5,6,1,1)
        layout.addWidget(QtGui.QLabel("s:"),6,5,1,1)
        self.s = Entrada_con_unidades(unidades.SpecificHeat, value=0)
        layout.addWidget(self.s,6,6,1,1)
        layout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding), 7,7,1,1)

        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok|QtGui.QDialogButtonBox.Cancel)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox,8,1,1,7)

        if config and config.has_option("MEoS", "reference"):
            self.setEnabled(False)
            if config.get("MEoS", "reference")=="OTO":
                self.OTO.setChecked(True)
            elif config.get("MEoS", "reference")=="NBP":
                self.NBP.setChecked(True)
            elif config.get("MEoS", "reference")=="IIR":
                self.IIR.setChecked(True)
            elif config.get("MEoS", "reference")=="ASHRAE":
                self.ASHRAE.setChecked(True)
            else:
                self.personalizado.setChecked(True)
                self.setEnabled(True)
                self.T.setValue(config.getfloat("MEoS", "T"))
                self.P.setValue(config.getfloat("MEoS", "P"))
                self.h.setValue(config.getfloat("MEoS", "h"))
                self.s.setValue(config.getfloat("MEoS", "s"))
        else:
            self.OTO.setChecked(True)
            self.setEnabled(False)

    def setEnabled(self, bool):
        """Enable custom entriees"""
        self.T.setEnabled(bool)
        self.P.setEnabled(bool)
        self.h.setEnabled(bool)
        self.s.setEnabled(bool)


class Dialog_InfoFluid(QtGui.QDialog):
    """Dialog to show parameter of element with meos"""
    def __init__(self, element, parent=None):
        """element: class of element to show info"""
        super(Dialog_InfoFluid, self).__init__(parent)
        layout = QtGui.QGridLayout(self)
        self.element=element

        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Name")+":"),1,1)
        self.name = QtGui.QLabel()
        layout.addWidget(self.name,1,2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "R name")+":"),2,1)
        self.r_name = QtGui.QLabel()
        layout.addWidget(self.r_name,2,2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Formula")+":"),3,1)
        self.formula = QtGui.QLabel()
        layout.addWidget(self.formula,3,2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "CAS number")+":"),4,1)
        self.CAS = QtGui.QLabel()
        layout.addWidget(self.CAS,4,2)
        layout.addItem(QtGui.QSpacerItem(30,30,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),1,3,3,1)

        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "M")+":"),1,4)
        self.M = Entrada_con_unidades(float, textounidad="g/mol", readOnly=True)
        layout.addWidget(self.M,1,5)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Tc")+":"),2,4)
        self.Tc = Entrada_con_unidades(unidades.Temperature, readOnly=True)
        layout.addWidget(self.Tc,2,5)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Pc")+":"),3,4)
        self.Pc = Entrada_con_unidades(unidades.Pressure, readOnly=True)
        layout.addWidget(self.Pc,3,5)
        layout.addWidget(QtGui.QLabel(u"ρc"+":"),4,4)
        self.rhoc = Entrada_con_unidades(unidades.Density, "DenGas", readOnly=True)
        layout.addWidget(self.rhoc,4,5)
        layout.addItem(QtGui.QSpacerItem(30,30,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),1,6,3,1)

        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "T triple")+":"),1,7)
        self.Tt = Entrada_con_unidades(unidades.Temperature, readOnly=True)
        layout.addWidget(self.Tt,1,8)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "T boiling")+":"),2,7)
        self.Tb = Entrada_con_unidades(unidades.Temperature, readOnly=True)
        layout.addWidget(self.Tb,2,8)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Dipole moment")+":"),3,7)
        self.momento = Entrada_con_unidades(unidades.DipoleMoment, readOnly=True)
        layout.addWidget(self.momento,3,8)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "F acentric")+":"),4,7)
        self.f_acent = Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.f_acent,4,8)

        layout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),5,1)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Equation")+": "),6,1)
        self.eq = QtGui.QComboBox()
        layout.addWidget(self.eq,6,2,1,7)
        self.stacked = QtGui.QStackedWidget()
        layout.addWidget(self.stacked,7,1,1,8)
        self.eq.currentIndexChanged.connect(self.stacked.setCurrentIndex)

        self.moreButton=QtGui.QPushButton(QtGui.QApplication.translate("pychemqt", "More..."))
        self.moreButton.clicked.connect(self.more)
        layout.addWidget(self.moreButton,9,1,1,1)
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Close)
        self.buttonBox.clicked.connect(self.reject)
        layout.addWidget(self.buttonBox,9,2,1,7)

        self.fill(element)

    def fill(self, element):
        """Fill values"""
        self.name.setText(element.name)
        self.r_name.setText(element.synonym)
        self.formula.setText(element.formula)
        self.CAS.setText(element.CASNumber)
        self.M.setValue(element.M)
        self.Tc.setValue(element.Tc)
        self.Pc.setValue(element.Pc)
        self.rhoc.setValue(element.rhoc)
        self.Tb.setValue(element.Tb)
        self.Tt.setValue(element.Tt)
        self.momento.setValue(element.momentoDipolar)
        self.f_acent.setValue(element.f_acent)

        for eq in element.eq:
            widget=Widget_MEoS_Data(eq)
            self.stacked.addWidget(widget)
            self.eq.addItem(eq["__name__"])

    def more(self):
        """Show parameter for transpor and ancillary equations"""
        dialog = transportDialog(self.element, parent=self)
        dialog.show()


class Widget_MEoS_Data(QtGui.QWidget):
    """Widget to show meos data"""
    def __init__(self, eq, parent=None):
        """eq: dict with equation parameter"""
        super(Widget_MEoS_Data, self).__init__(parent)
        gridLayout = QtGui.QGridLayout(self)
        ref=QtGui.QLabel(eq["__doc__"])
        ref.setWordWrap(True)
        gridLayout.addWidget(ref,1,1)

        tabWidget = QtGui.QTabWidget()
        gridLayout.addWidget(tabWidget,3,1)

        # Cp tab
        if "ao_log" in eq["cp"]:
            pass
        else:
            tab1 = QtGui.QWidget()
            tabWidget.addTab(tab1,QtGui.QApplication.translate("pychemqt", "Cp"))
            gridLayout_Ideal=QtGui.QGridLayout(tab1)
            label=QtGui.QLabel()
            label.setAlignment(QtCore.Qt.AlignCenter)
            label.setPixmap(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equation/MEoS ideal.png"))
            gridLayout_Ideal.addWidget(label,1,1,1,3)
            self.Tabla_Cp_poly=Tabla(2, horizontalHeader=["n", "d"], stretch=False, readOnly=True)
            gridLayout_Ideal.addWidget(self.Tabla_Cp_poly,2,1)
            self.Tabla_Cp_exp=Tabla(2, horizontalHeader=["m", u"θ"], stretch=False, readOnly=True)
            gridLayout_Ideal.addWidget(self.Tabla_Cp_exp,2,2)
            self.Tabla_Cp_hyp=Tabla(2, horizontalHeader=["l", u"ψ"], stretch=False, readOnly=True)
            gridLayout_Ideal.addWidget(self.Tabla_Cp_hyp,2,3)

        if eq["__type__"]=="Helmholtz":
            label=QtGui.QLabel()
            label.setAlignment(QtCore.Qt.AlignCenter)
            label.setPixmap(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equation/MEoS.png"))
            gridLayout.addWidget(label,2,1)

            # Polinomial tab
            tab2 = QtGui.QWidget()
            tabWidget.addTab(tab2,QtGui.QApplication.translate("pychemqt", "Polinomial"))
            gridLayout_pol=QtGui.QGridLayout(tab2)
            label=QtGui.QLabel()
            label.setAlignment(QtCore.Qt.AlignCenter)
            label.setPixmap(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equation/MEoS lineal.png"))
            gridLayout_pol.addWidget(label,1,1)
            self.Tabla_lineal=Tabla(3, horizontalHeader=["n", "t", "d"], stretch=False, readOnly=True)
            gridLayout_pol.addWidget(self.Tabla_lineal,2,1)

            # Exponencial tab
            tab3 = QtGui.QWidget()
            tabWidget.addTab(tab3,QtGui.QApplication.translate("pychemqt", "Exponential"))
            gridLayout_Exp=QtGui.QGridLayout(tab3)
            label=QtGui.QLabel()
            label.setAlignment(QtCore.Qt.AlignCenter)
            label.setPixmap(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equation/MEoS exponential.png"))
            gridLayout_Exp.addWidget(label,1,1)
            self.Tabla_exponential=Tabla(5, horizontalHeader=["n", "t", "d", u"γ", "c"], stretch=False, readOnly=True)
            gridLayout_Exp.addWidget(self.Tabla_exponential,2,1)

            # Gaussian tab
            tab4 = QtGui.QWidget()
            tabWidget.addTab(tab4,QtGui.QApplication.translate("pychemqt", "Gaussian"))
            gridLayout_gauss=QtGui.QGridLayout(tab4)
            label=QtGui.QLabel()
            label.setAlignment(QtCore.Qt.AlignCenter)
            label.setPixmap(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equation/MEoS gaussian.png"))
            gridLayout_gauss.addWidget(label,1,1)
            self.Tabla_gauss=Tabla(7, horizontalHeader=["n", "t", "d", u"η", u"ε", u"β", u"γ"], stretch=False, readOnly=True)
            gridLayout_gauss.addWidget(self.Tabla_gauss,2,1)

            # Non analytic tab
            tab5 = QtGui.QWidget()
            tabWidget.addTab(tab5,QtGui.QApplication.translate("pychemqt", "Non analytic"))
            gridLayout_NA=QtGui.QGridLayout(tab5)
            label=QtGui.QLabel()
            label.setAlignment(QtCore.Qt.AlignCenter)
            label.setPixmap(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equation/MEoS non analitic.png"))
            gridLayout_NA.addWidget(label,1,1)
            label2=QtGui.QLabel()
            label2.setAlignment(QtCore.Qt.AlignCenter)
            label2.setPixmap(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equation/MEoS delta.png"))
            gridLayout_NA.addWidget(label2,2,1)
            self.Tabla_noanalytic=Tabla(8, horizontalHeader=["n", "a", "b", "A", "B", "C", "D", u"β"], stretch=False, readOnly=True)
            gridLayout_NA.addWidget(self.Tabla_noanalytic,3,1)

            # Hand Sphere tab
            tab6 = QtGui.QWidget()
            tabWidget.addTab(tab6,QtGui.QApplication.translate("pychemqt", "Hard Sphere"))
            gridLayout_HE=QtGui.QGridLayout(tab6)
            label=QtGui.QLabel()
            label.setAlignment(QtCore.Qt.AlignCenter)
            label.setPixmap(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equation/MEoS Hard Sphere.png"))
            gridLayout_HE.addWidget(label,1,1,1,2)
            gridLayout_HE.addWidget(QtGui.QLabel(u"φ:"),2,1)
            self.fi = Entrada_con_unidades(float, readOnly=True)
            gridLayout_HE.addWidget(self.fi,2,2)
            gridLayout_HE.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),3,1,1,2)

        elif eq["__type__"]=="MBWR":
            #Pestaña MBWR
            tab2 = QtGui.QWidget()
            tabWidget.addTab(tab2,QtGui.QApplication.translate("pychemqt", "MBWR"))
            gridLayout_MBWR=QtGui.QGridLayout(tab2)
            label=QtGui.QLabel()
            label.setAlignment(QtCore.Qt.AlignCenter)
            label.setPixmap(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equation/MEoS MBWR.png"))
            gridLayout_MBWR.addWidget(label,1,1)
            self.Tabla_MBWR=Tabla(1, horizontalHeader=["b"], stretch=False, readOnly=True)
            gridLayout_MBWR.addWidget(self.Tabla_MBWR,2,1)

        self.fill(eq)


    def fill(self, eq):
        format = {"format": 1, "total": 5}

        if "ao_log" in eq["cp"]:
            pass
        else:
            self.Tabla_Cp_poly.setColumn(0, [eq["cp"]["ao"]]+eq["cp"]["an"], **format)
            self.Tabla_Cp_poly.setColumn(1, [0]+eq["cp"]["pow"], **format)
            self.Tabla_Cp_poly.resizeColumnsToContents()
            self.Tabla_Cp_exp.setColumn(0, eq["cp"]["ao_exp"], **format)
            self.Tabla_Cp_exp.setColumn(1, eq["cp"]["exp"], **format)
            self.Tabla_Cp_exp.resizeColumnsToContents()
            self.Tabla_Cp_hyp.setColumn(0, eq["cp"]["ao_hyp"], **format)
            self.Tabla_Cp_hyp.setColumn(1, eq["cp"]["hyp"], **format)
            self.Tabla_Cp_hyp.resizeColumnsToContents()

        if eq["__type__"]=="Helmholtz":
            if eq.get("nr1", []):
                self.Tabla_lineal.setColumn(0, eq["nr1"], **format)
                self.Tabla_lineal.setColumn(1, eq["t1"], **format)
                self.Tabla_lineal.setColumn(2, eq["d1"], **format)
            if eq.get("nr2", []):
                self.Tabla_exponential.setColumn(0, eq["nr2"], **format)
                self.Tabla_exponential.setColumn(1, eq["t2"], **format)
                self.Tabla_exponential.setColumn(2, eq["d2"], **format)
                self.Tabla_exponential.setColumn(3, eq["gamma2"], **format)
                self.Tabla_exponential.setColumn(4, eq["c2"], **format)
            if eq.get("nr3", []):
                self.Tabla_gauss.setColumn(0, eq["nr3"], **format)
                self.Tabla_gauss.setColumn(1, eq["t3"], **format)
                self.Tabla_gauss.setColumn(2, eq["d3"], **format)
                self.Tabla_gauss.setColumn(3, eq["alfa3"], **format)
                self.Tabla_gauss.setColumn(4, eq["beta3"], **format)
                self.Tabla_gauss.setColumn(5, eq["gamma3"], **format)
                self.Tabla_gauss.setColumn(6, eq["epsilon3"], **format)
            if eq.get("nr4", []):
                self.Tabla_noanalytic.setColumn(0, eq["nr4"], **format)
                self.Tabla_noanalytic.setColumn(1, eq["a4"], **format)
                self.Tabla_noanalytic.setColumn(2, eq["b4"], **format)
                self.Tabla_noanalytic.setColumn(3, eq["A"], **format)
                self.Tabla_noanalytic.setColumn(4, eq["B"], **format)
                self.Tabla_noanalytic.setColumn(5, eq["C"], **format)
                self.Tabla_noanalytic.setColumn(6, eq["D"], **format)
                self.Tabla_noanalytic.setColumn(7, eq["beta4"], **format)
            self.Tabla_lineal.resizeColumnsToContents()
            self.Tabla_exponential.resizeColumnsToContents()
            self.Tabla_gauss.resizeColumnsToContents()
            self.Tabla_noanalytic.resizeColumnsToContents()

        elif eq["__type__"]=="MBWR":
            self.Tabla_MBWR.setColumn(0, eq["b"][1:], **format)
            self.Tabla_MBWR.resizeColumnsToContents()


class transportDialog(QtGui.QDialog):
    """Dialog to show parameters for transport and ancillary equations"""
    def __init__(self, element, parent=None):
        super(transportDialog, self).__init__(parent)
        gridLayout = QtGui.QGridLayout(self)
        self.element=element

        tabWidget = QtGui.QTabWidget()
        gridLayout.addWidget(tabWidget,1,1)

        #Tab viscosity
        tab3 = QtGui.QWidget()
        tabWidget.addTab(tab3,QtGui.QApplication.translate("pychemqt", "Viscosity"))
        gridLayout_viscosity=QtGui.QGridLayout(tab3)

        gridLayout_viscosity.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Equation")+": "),1,1)
        self.eqVisco = QtGui.QComboBox()
        gridLayout_viscosity.addWidget(self.eqVisco,1,2)
        gridLayout_viscosity.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Fixed),1,3)
        self.stackedVisco = QtGui.QStackedWidget()
        gridLayout_viscosity.addWidget(self.stackedVisco,2,1,1,3)
        self.eqVisco.currentIndexChanged.connect(self.stackedVisco.setCurrentIndex)

        if element._viscosity is not None:
            for eq in element._viscosity:
                widget=Widget_Viscosity_Data(element, eq)
                self.stackedVisco.addWidget(widget)
                self.eqVisco.addItem(eq["__name__"])
        else:
            self.eqVisco.addItem(QtGui.QApplication.translate("pychemqt", "Not Implemented"))


        #Tab thermal conductivity
        tab4 = QtGui.QWidget()
        tabWidget.addTab(tab4,QtGui.QApplication.translate("pychemqt", "Thermal Conductivity"))
        gridLayout_conductivity=QtGui.QGridLayout(tab4)

        gridLayout_conductivity.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Equation")+": "),1,1)
        self.eqThermo = QtGui.QComboBox()
        gridLayout_conductivity.addWidget(self.eqThermo,1,2)
        gridLayout_conductivity.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Fixed),1,3)
        self.stackedThermo = QtGui.QStackedWidget()
        gridLayout_conductivity.addWidget(self.stackedThermo,2,1,1,3)
        self.eqThermo.currentIndexChanged.connect(self.stackedThermo.setCurrentIndex)

        if element._thermal is not None:
            for eq in element._thermal:
                widget=Widget_Conductivity_Data(element, eq)
                self.stackedThermo.addWidget(widget)
                self.eqThermo.addItem(eq["__name__"])
        else:
            self.eqThermo.addItem(QtGui.QApplication.translate("pychemqt", "Not Implemented"))

        #Tab dielectric constant
        tab1 = QtGui.QWidget()
        tabWidget.addTab(tab1,QtGui.QApplication.translate("pychemqt", "Dielectric"))
        gridLayout_dielectric=QtGui.QGridLayout(tab1)

        if element._dielectric:
            label=QtGui.QLabel(element._Dielectric.__doc__)
            label.setWordWrap(True)
            gridLayout_dielectric.addWidget(label,1,1)

            self.Table_Dielectric=Tabla(1, verticalHeader=True, filas=5, stretch=False, readOnly=True)
            gridLayout_dielectric.addWidget(self.Table_Dielectric,2,1)
            i=0
            for key, valor in element._dielectric.iteritems():
                self.Table_Dielectric.setVerticalHeaderItem(i,QtGui.QTableWidgetItem(key))
                self.Table_Dielectric.setItem(0, i, QtGui.QTableWidgetItem(str(valor)))
                i+=1
            self.Table_Dielectric.resizeColumnsToContents()

        elif element._Dielectric != meos.MEoS._Dielectric:
            label=QtGui.QLabel(element._Dielectric.__doc__)
            label.setWordWrap(True)
            gridLayout_dielectric.addWidget(label,1,1)
            self.codigo_Dielectric=SimplePythonEditor()
            self.codigo_Dielectric.setText(inspect.getsource(element._Dielectric))
            gridLayout_dielectric.addWidget(self.codigo_Dielectric,2,1)
        else:
            gridLayout_dielectric.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Not Implemented")),1,1)
            gridLayout_dielectric.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),2,1)

        #Tab surface tension
        tab2 = QtGui.QWidget()
        tabWidget.addTab(tab2,QtGui.QApplication.translate("pychemqt", "Surface Tension"))
        gridLayout_surface=QtGui.QGridLayout(tab2)

        if element._surface:
            label=QtGui.QLabel(element._Surface.__doc__)
            label.setWordWrap(True)
            gridLayout_surface.addWidget(label,1,1)

            self.Table_Surface=Tabla(2, horizontalHeader=[u"σ", "n"], verticalHeader=True, stretch=False, readOnly=True)
            self.Table_Surface.setColumn(0, element._surface["sigma"])
            self.Table_Surface.setColumn(1, element._surface["exp"])
            gridLayout_surface.addWidget(self.Table_Surface,2,1)
            self.Table_Surface.resizeColumnsToContents()

        elif element._Surface != meos.MEoS._Surface:
            label=QtGui.QLabel(element._Surface.__doc__)
            label.setWordWrap(True)
            gridLayout_surface.addWidget(label,1,1)
            self.codigo_Surface=SimplePythonEditor()
            self.codigo_Surface.setText(inspect.getsource(element._Surface))
            gridLayout_surface.addWidget(self.codigo_Surface,2,1)
        else:
            gridLayout_surface.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Not Implemented")),1,1)
            gridLayout_surface.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),2,1)

        #Tab liquid density
        tab5 = QtGui.QWidget()
        tabWidget.addTab(tab5,QtGui.QApplication.translate("pychemqt", "Liquid Density"))
        gridLayout_liquid_density=QtGui.QGridLayout(tab5)

        if element._liquid_Density:
            label=QtGui.QLabel(element._Liquid_Density.__doc__)
            label.setWordWrap(True)
            gridLayout_liquid_density.addWidget(label,1,1)

            self.Table_Liquid_Density=Tabla(2, horizontalHeader=[u"ao", "n"], verticalHeader=True, stretch=False, readOnly=True)
            self.Table_Liquid_Density.setColumn(0, element._liquid_Density["ao"])
            self.Table_Liquid_Density.setColumn(1, element._liquid_Density["exp"])
            gridLayout_liquid_density.addWidget(self.Table_Liquid_Density,2,1)
            self.Table_Liquid_Density.resizeColumnsToContents()

        elif element._Liquid_Density != meos.MEoS._Liquid_Density:
            label=QtGui.QLabel(element._Liquid_Density.__doc__)
            label.setWordWrap(True)
            gridLayout_liquid_density.addWidget(label,1,1)
            self.codigo_Liquid_Density=SimplePythonEditor()
            self.codigo_Liquid_Density.setText(inspect.getsource(element._Liquid_Density))
            gridLayout_liquid_density.addWidget(self.codigo_Liquid_Density,2,1)
        else:
            gridLayout_liquid_density.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Not Implemented")),1,1)
            gridLayout_liquid_density.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),2,1)

        #Tab vapor density
        tab6 = QtGui.QWidget()
        tabWidget.addTab(tab6,QtGui.QApplication.translate("pychemqt", "Vapor Density"))
        gridLayout_vapor_density=QtGui.QGridLayout(tab6)

        if element._vapor_Density:
            label=QtGui.QLabel(element._Vapor_Density.__doc__)
            label.setWordWrap(True)
            gridLayout_vapor_density.addWidget(label,1,1)

            self.Table_Vapor_Density=Tabla(2, horizontalHeader=[u"ao", "n"], verticalHeader=True, stretch=False, readOnly=True)
            self.Table_Vapor_Density.setColumn(0, element._vapor_Density["ao"])
            self.Table_Vapor_Density.setColumn(1, element._vapor_Density["exp"])
            gridLayout_vapor_density.addWidget(self.Table_Vapor_Density,2,1)
            self.Table_Vapor_Density.resizeColumnsToContents()

        elif element._Vapor_Density != meos.MEoS._Vapor_Density:
            label=QtGui.QLabel(element._Vapor_Density.__doc__)
            label.setWordWrap(True)
            gridLayout_vapor_density.addWidget(label,1,1)
            self.codigo_Vapor_Density=SimplePythonEditor()
            self.codigo_Vapor_Density.setText(inspect.getsource(element._Vapor_Density))
            gridLayout_vapor_density.addWidget(self.codigo_Vapor_Density,2,1)
        else:
            gridLayout_vapor_density.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Not Implemented")),1,1)
            gridLayout_vapor_density.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),2,1)

        #Tab vapor presure
        tab7 = QtGui.QWidget()
        tabWidget.addTab(tab7,QtGui.QApplication.translate("pychemqt", "Vapor Pressure"))
        gridLayout_vapor_pressure=QtGui.QGridLayout(tab7)

        if element._vapor_Pressure:
            label=QtGui.QLabel(element._Vapor_Pressure.__doc__)
            label.setWordWrap(True)
            gridLayout_vapor_pressure.addWidget(label,1,1)

            self.Table_Vapor_Pressure=Tabla(2, horizontalHeader=[u"ao", "n"], verticalHeader=True, stretch=False, readOnly=True)
            self.Table_Vapor_Pressure.setColumn(0, element._vapor_Pressure["ao"])
            self.Table_Vapor_Pressure.setColumn(1, element._vapor_Pressure["exp"])
            gridLayout_vapor_pressure.addWidget(self.Table_Vapor_Pressure,2,1)
            self.Table_Vapor_Pressure.resizeColumnsToContents()

        elif element._Vapor_Pressure != meos.MEoS._Vapor_Pressure:
            label=QtGui.QLabel(element._Vapor_Pressure.__doc__)
            label.setWordWrap(True)
            gridLayout_vapor_pressure.addWidget(label,1,1)
            self.codigo_Vapor_Pressure=SimplePythonEditor()
            self.codigo_Vapor_Pressure.setText(inspect.getsource(element._Vapor_Pressure))
            gridLayout_vapor_pressure.addWidget(self.codigo_Vapor_Pressure,2,1)
        else:
            gridLayout_vapor_pressure.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Not Implemented")),1,1)
            gridLayout_vapor_pressure.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),2,1)

        #Tab melting presure
        tab8 = QtGui.QWidget()
        tabWidget.addTab(tab8,QtGui.QApplication.translate("pychemqt", "Melting Pressure"))
        gridLayout_melting_pressure=QtGui.QGridLayout(tab8)

        if element._melting:
            label=QtGui.QLabel(element._Melting_Pressure.__doc__)
            label.setWordWrap(True)
            gridLayout_melting_pressure.addWidget(label,1,1)

            self.Table_Melting_Pressure=Tabla(6, horizontalHeader=["a1", "n1", "a2", "n2", "a3", "n3"], verticalHeader=True, stretch=False, readOnly=True)
            self.Table_Melting_Pressure.setColumn(0, element._melting["a1"])
            self.Table_Melting_Pressure.setColumn(1, element._melting["exp1"])
            self.Table_Melting_Pressure.setColumn(2, element._melting["a2"])
            self.Table_Melting_Pressure.setColumn(3, element._melting["exp2"])
            self.Table_Melting_Pressure.setColumn(4, element._melting["a3"])
            self.Table_Melting_Pressure.setColumn(5, element._melting["exp3"])
            gridLayout_melting_pressure.addWidget(self.Table_Melting_Pressure,2,1)
            self.Table_Melting_Pressure.resizeColumnsToContents()

        elif element._Melting_Pressure != meos.MEoS._Melting_Pressure:
            label=QtGui.QLabel(element._Melting_Pressure.__doc__)
            label.setWordWrap(True)
            gridLayout_melting_pressure.addWidget(label,1,1)
            self.codigo_Melting_Pressure=SimplePythonEditor()
            self.codigo_Melting_Pressure.setText(inspect.getsource(element._Melting_Pressure))
            gridLayout_melting_pressure.addWidget(self.codigo_Melting_Pressure,2,1)
        else:
            gridLayout_melting_pressure.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Not Implemented")),1,1)
            gridLayout_melting_pressure.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),2,1)

        #Tab sublimation presure
        tab9 = QtGui.QWidget()
        tabWidget.addTab(tab9,QtGui.QApplication.translate("pychemqt", "Sublimation Pressure"))
        gridLayout__sublimation_pressure=QtGui.QGridLayout(tab9)

        if element._sublimation:
            label=QtGui.QLabel(element._Melting_Pressure.__doc__)
            label.setWordWrap(True)
            gridLayout__sublimation_pressure.addWidget(label,1,1)

            self.Table_Sublimation_Pressure=Tabla(6, horizontalHeader=["a1", "n1", "a2", "n2", "a3", "n3"], verticalHeader=True, stretch=False, readOnly=True)
            self.Table_Sublimation_Pressure.setColumn(0, element._sublimation["a1"])
            self.Table_Sublimation_Pressure.setColumn(1, element._sublimation["exp1"])
            self.Table_Sublimation_Pressure.setColumn(2, element._sublimation["a2"])
            self.Table_Sublimation_Pressure.setColumn(3, element._sublimation["exp2"])
            self.Table_Sublimation_Pressure.setColumn(4, element._sublimation["a3"])
            self.Table_Sublimation_Pressure.setColumn(5, element._sublimation["exp3"])
            gridLayout__sublimation_pressure.addWidget(self.Table_Sublimation_Pressure,2,1)
            self.Table_Sublimation_Pressure.resizeColumnsToContents()

        elif element._Sublimation_Pressure != meos.MEoS._Sublimation_Pressure:
            label=QtGui.QLabel(element._Sublimation_Pressure.__doc__)
            label.setWordWrap(True)
            gridLayout__sublimation_pressure.addWidget(label,1,1)
            self.codigo_Sublimation_Pressure=SimplePythonEditor()
            self.codigo_Sublimation_Pressure.setText(inspect.getsource(element._Sublimation_Pressure))
            gridLayout__sublimation_pressure.addWidget(self.codigo_Sublimation_Pressure,2,1)
        else:
            gridLayout__sublimation_pressure.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Not Implemented")),1,1)
            gridLayout__sublimation_pressure.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),2,1)

        #Tab Peng-Robinson
        tab10 = QtGui.QWidget()
        tabWidget.addTab(tab10,QtGui.QApplication.translate("pychemqt", "Peng-Robinson"))
        gridLayout_PengRobinson=QtGui.QGridLayout(tab10)

        if element._PR:
            label=QtGui.QLabel(element._PengRobinson.__doc__)
            label.setWordWrap(True)
            gridLayout_PengRobinson.addWidget(label,1,1,1,3)
            gridLayout_PengRobinson.addWidget(QtGui.QLabel("C"),2,1)
            self.PR=Entrada_con_unidades(float, decimales=6, value=element._PR, readOnly=True)
            gridLayout_PengRobinson.addWidget(self.PR,2,2)
            gridLayout_PengRobinson.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),3,1,1,3)
        else:
            gridLayout_PengRobinson.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "No Peneloux correction")),1,1)
            gridLayout_PengRobinson.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),2,1)

        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Close)
        self.buttonBox.clicked.connect(self.reject)
        gridLayout.addWidget(self.buttonBox,2,1)


class Widget_Viscosity_Data(QtGui.QWidget):
    """Widget to show viscosity data"""
    def __init__(self, element, eq, parent=None):
        """
        element: element class for code extract
        eq: dict with viscosity parameter"""
        super(Widget_Viscosity_Data, self).__init__(parent)
        gridLayout = QtGui.QGridLayout(self)
        if eq["eq"] == 0:
            doc = element.__getattribute__(element, eq["method"]).__doc__
            ref = QtGui.QLabel(doc)
        else:
            ref=QtGui.QLabel(eq["__doc__"])
        ref.setWordWrap(True)
        gridLayout.addWidget(ref,1,1,1,3)


        if eq["eq"] == 0:
            # Hardcoded method, show code
            self.codigo_Viscosity=SimplePythonEditor()
            code = ""
            for method in eq.get("__code__", ()):
                code += inspect.getsource(method)
                code += os.linesep
            code += inspect.getsource(element.__getattribute__(element, eq["method"]))
            self.codigo_Viscosity.setText(code)
            gridLayout.addWidget(self.codigo_Viscosity,2,1, 1, 3)
        elif eq["eq"] == 1:
            gridLayout.addWidget(QtGui.QLabel(u"ε/k"),4,1)
            self.ek=Entrada_con_unidades(float, value=eq["ek"], readOnly=True)
            gridLayout.addWidget(self.ek,4,2)
            gridLayout.addWidget(QtGui.QLabel(u"σ"),5,1)
            self.sigma=Entrada_con_unidades(float, value=eq["sigma"], readOnly=True)
            gridLayout.addWidget(self.sigma,5,2)
            tab = QtGui.QTabWidget()
            gridLayout.addWidget(tab, 6, 1, 1, 3)

            # Integral collision
            self.Tabla_Collision = Tabla(1, horizontalHeader=["b"], stretch=False, readOnly=True)
            if "collision" in eq:
                self.Tabla_Collision.setColumn(0, eq["collision"])
                self.Tabla_Collision.resizeColumnsToContents()
            else:
                self.Tabla_Collision.setDisabled(True)
            tab.addTab(self.Tabla_Collision, QtGui.QApplication.translate("pychemqt", "Collision"))

            # Virial
            self.Tabla_Virial = Tabla(2, horizontalHeader=["n", "t"], stretch=False, readOnly=True)
            if "n_virial" in eq:
                self.Tabla_Virial.setColumn(0, eq["n_virial"])
                self.Tabla_Virial.setColumn(1, eq["t_virial"])
                self.Tabla_Virial.resizeColumnsToContents()
            else:
                self.Tabla_Virial.setDisabled(True)
            tab.addTab(self.Tabla_Virial, QtGui.QApplication.translate("pychemqt", "Virial"))

            # Close-packed
            self.Tabla_Packed = Tabla(2, horizontalHeader=["n", "t"], stretch=False, readOnly=True)
            if "n_packed" in eq:
                self.Tabla_Packed.setColumn(0, eq["n_packed"])
                self.Tabla_Packed.setColumn(1, eq["t_packed"])
                self.Tabla_Packed.resizeColumnsToContents()
            else:
                self.Tabla_Packed.setDisabled(True)
            tab.addTab(self.Tabla_Packed, QtGui.QApplication.translate("pychemqt", "Close-packed density"))

            # polynomial term
            self.Tabla_Visco1 = Tabla(5, horizontalHeader=["n", "t", "d", "g", "c"], stretch=False, readOnly=True)
            if "n_poly" in eq:
                self.Tabla_Visco1.setColumn(0, eq["n_poly"])
                self.Tabla_Visco1.setColumn(1, eq["t_poly"])
                self.Tabla_Visco1.setColumn(2, eq["d_poly"])
                self.Tabla_Visco1.setColumn(3, eq["g_poly"])
                self.Tabla_Visco1.setColumn(4, eq["c_poly"])
                self.Tabla_Visco1.resizeColumnsToContents()
            else:
                self.Tabla_Visco1.setDisabled(True)
            tab.addTab(self.Tabla_Visco1, QtGui.QApplication.translate("pychemqt", "Polinomial"))

            # numerator of rational poly
            self.Tabla_numerator = Tabla(5, horizontalHeader=["n", "t", "d", "g", "c"], stretch=False, readOnly=True)
            if "n_num" in eq:
                self.Tabla_numerator.setColumn(0, eq["n_num"])
                self.Tabla_numerator.setColumn(1, eq["t_num"])
                self.Tabla_numerator.setColumn(2, eq["d_num"])
                self.Tabla_numerator.setColumn(3, eq["c_num"])
                self.Tabla_numerator.setColumn(4, eq["g_num"])
                self.Tabla_numerator.resizeColumnsToContents()
            else:
                self.Tabla_numerator.setDisabled(True)
            tab.addTab(self.Tabla_numerator, QtGui.QApplication.translate("pychemqt", "Numerator"))

            # denominator of rational poly
            self.Tabla_denominator = Tabla(5, horizontalHeader=["n", "t", "d", "g", "c"], stretch=False, readOnly=True)
            if "n_den" in eq:
                self.Tabla_denominator.setColumn(0, eq["n_den"])
                self.Tabla_denominator.setColumn(1, eq["t_den"])
                self.Tabla_denominator.setColumn(2, eq["d_den"])
                self.Tabla_denominator.setColumn(3, eq["c_den"])
                self.Tabla_denominator.setColumn(4, eq["g_den"])
                self.Tabla_denominator.resizeColumnsToContents()
            else:
                self.Tabla_denominator.setDisabled(True)
            tab.addTab(self.Tabla_denominator, QtGui.QApplication.translate("pychemqt", "Denominator"))
            gridLayout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),10,3)

        elif eq["eq"] == 2:
            gridLayout.addWidget(QtGui.QLabel(u"ε/k"),4,1)
            self.ek=Entrada_con_unidades(float, value=eq["ek"], readOnly=True)
            gridLayout.addWidget(self.ek,4,2)
            gridLayout.addWidget(QtGui.QLabel(u"σ"),5,1)
            self.sigma=Entrada_con_unidades(float, value=eq["sigma"], readOnly=True)
            gridLayout.addWidget(self.sigma,5,2)
            self.Tabla_Visco2 = Tabla(3, horizontalHeader=["b", "F", "E"], stretch=False, readOnly=True)
            if "collision" in eq:
                self.Tabla_Visco2.setColumn(0, eq["collision"])
            self.Tabla_Visco2.setColumn(1, eq["F"])
            self.Tabla_Visco2.setColumn(2, eq["E"])
            self.Tabla_Visco2.resizeColumnsToContents()
            gridLayout.addWidget(self.Tabla_Visco2, 6, 1, 1, 3)
            gridLayout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),10,3)

        elif eq["eq"] == 3:
            self.Tabla_Visco3 = Tabla(8, horizontalHeader=["n-poly", "t-poly", "n-num", "t-num", "d-num", "n-den", "t-den", "d-den"], stretch=False, readOnly=True)
            if "n_poly" in eq:
                self.Tabla_Visco3.setColumn(0, eq["n_poly"])
                self.Tabla_Visco3.setColumn(1, eq["t_poly"])
            if "n_num" in eq:
                self.Tabla_Visco3.setColumn(2, eq["n_num"])
                self.Tabla_Visco3.setColumn(3, eq["t_num"])
                self.Tabla_Visco3.setColumn(4, eq["d_num"])
            if "n_den" in eq:
                self.Tabla_Visco3.setColumn(5, eq["n_den"])
                self.Tabla_Visco3.setColumn(6, eq["t_den"])
                self.Tabla_Visco3.setColumn(7, eq["d_den"])
            self.Tabla_Visco3.resizeColumnsToContents()
            gridLayout.addWidget(self.Tabla_Visco3, 4, 1, 1, 3)
            gridLayout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),10,3)

        elif eq["eq"] == 4:
            gridLayout.addWidget(QtGui.QLabel(u"ε/k"),4,1)
            self.ek=Entrada_con_unidades(float, value=eq["ek"], readOnly=True)
            gridLayout.addWidget(self.ek,4,2)
            gridLayout.addWidget(QtGui.QLabel(u"σ"),5,1)
            self.sigma=Entrada_con_unidades(float, value=eq["sigma"], readOnly=True)
            gridLayout.addWidget(self.sigma,5,2)
            self.Tabla_Visco4 = Tabla(7, horizontalHeader=["a", "b", "c", "A", "B", "C", "D"], stretch=False, readOnly=True)
            format = {"format": 1, "decimales": 10}
            self.Tabla_Visco4.setColumn(0, eq["a"], **format)
            self.Tabla_Visco4.setColumn(1, eq["b"], **format)
            self.Tabla_Visco4.setColumn(2, eq["c"], **format)
            self.Tabla_Visco4.setColumn(3, eq["A"], **format)
            self.Tabla_Visco4.setColumn(4, eq["B"], **format)
            self.Tabla_Visco4.setColumn(5, eq["C"], **format)
            self.Tabla_Visco4.setColumn(6, eq["D"], **format)
            self.Tabla_Visco4.resizeColumnsToContents()
            gridLayout.addWidget(self.Tabla_Visco4, 6, 1, 1, 3)
            gridLayout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),10,3)

        elif eq["eq"] == 5:
            gridLayout.addWidget(QtGui.QLabel("w"),4,1)
            self.w=Entrada_con_unidades(float, value=eq["w"], readOnly=True)
            gridLayout.addWidget(self.w,4,2)
            gridLayout.addWidget(QtGui.QLabel("mur"),5,1)
            self.mur=Entrada_con_unidades(float, value=eq["mur"], readOnly=True)
            gridLayout.addWidget(self.mur,5,2)
            gridLayout.addWidget(QtGui.QLabel(u"ε/k"),6,1)
            self.k=Entrada_con_unidades(float, value=eq["k"], readOnly=True)
            gridLayout.addWidget(self.k,6,2)
            gridLayout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),10,3)


class Widget_Conductivity_Data(QtGui.QWidget):
    """Widget to show thermal conductivity data"""
    def __init__(self, element, eq, parent=None):
        """
        element: element class for code extract
        eq: dict with thermal conductivity parameter"""
        super(Widget_Conductivity_Data, self).__init__(parent)
        gridLayout = QtGui.QGridLayout(self)
        if eq["eq"] == 0:
            doc = element.__getattribute__(element, eq["method"]).__doc__
            ref = QtGui.QLabel(doc)
        else:
            ref=QtGui.QLabel(eq["__doc__"])
        ref.setWordWrap(True)
        gridLayout.addWidget(ref,1,1,1,3)

        if eq["eq"] == 0:
            # Hardcoded method, show code
            self.code=SimplePythonEditor()
            code = ""
            for method in eq.get("__code__", ()):
                code += inspect.getsource(method)
                code += os.linesep
            code += inspect.getsource(element.__getattribute__(element, eq["method"]))
            self.code.setText(code)
            gridLayout.addWidget(self.code,2,1, 1, 3)

        elif eq["eq"] == 1:
            self.Tabla_Therm1 = Tabla(11, horizontalHeader=["no", "co", "noden", "toden", "nb", "tb", "db", "cb", "nbden", "tbden", "dbden"], stretch=False, readOnly=True)
            if "no" in eq:
                self.Tabla_Therm1.setColumn(0, eq["no"])
                self.Tabla_Therm1.setColumn(1, eq["co"])
            if "noden" in eq:
                self.Tabla_Therm1.setColumn(2, eq["noden"])
                self.Tabla_Therm1.setColumn(3, eq["toden"])
            if "nb" in eq:
                self.Tabla_Therm1.setColumn(4, eq["nb"])
                self.Tabla_Therm1.setColumn(5, eq["tb"])
                self.Tabla_Therm1.setColumn(6, eq["db"])
                self.Tabla_Therm1.setColumn(7, eq["cb"])
            if "nbden" in eq:
                self.Tabla_Therm1.setColumn(8, eq["nbden"])
                self.Tabla_Therm1.setColumn(9, eq["tbden"])
                self.Tabla_Therm1.setColumn(10, eq["dbden"])
            self.Tabla_Therm1.resizeColumnsToContents()
            gridLayout.addWidget(self.Tabla_Therm1, 3, 1, 1, 3)

        elif eq["eq"] == 2:
            self.Tabla_Therm2 = Tabla(2, horizontalHeader=["E", "G"], stretch=False, readOnly=True)
            self.Tabla_Therm2.setColumn(0, eq["E"])
            self.Tabla_Therm2.setColumn(1, eq["G"])
            self.Tabla_Therm2.resizeColumnsToContents()
            gridLayout.addWidget(self.Tabla_Therm2, 3, 1, 1, 3)

        elif eq["eq"] == 3:
            self.Tabla_Therm3 = Tabla(3, horizontalHeader=["b", "F", "E"], stretch=False, readOnly=True)
            self.Tabla_Therm3.setColumn(0, eq["b"])
            self.Tabla_Therm3.setColumn(1, eq["F"])
            self.Tabla_Therm3.setColumn(2, eq["E"])
            self.Tabla_Therm3.resizeColumnsToContents()
            gridLayout.addWidget(self.Tabla_Therm3, 3, 1, 1, 3)

            parameter = QtGui.QWidget()
            gridLayout.addWidget(parameter, 4, 1, 1, 3)
            lyt = QtGui.QGridLayout(parameter)
            lyt.addWidget(QtGui.QLabel(u"ε/k"),1,1)
            self.ek=Entrada_con_unidades(float, value=eq["ek"], readOnly=True)
            lyt.addWidget(self.ek,1,2)
            lyt.addWidget(QtGui.QLabel(u"σ"),2,1)
            self.sigma=Entrada_con_unidades(float, value=eq["sigma"], readOnly=True)
            lyt.addWidget(self.sigma,2,2)
            lyt.addWidget(QtGui.QLabel(u"Nchapman"),3,1)
            self.Nchapman=Entrada_con_unidades(float, value=eq["Nchapman"], readOnly=True)
            lyt.addWidget(self.Nchapman,3,2)
            lyt.addWidget(QtGui.QLabel(u"tchapman"),4,1)
            self.tchapman=Entrada_con_unidades(float, value=eq["tchapman"], readOnly=True)
            lyt.addWidget(self.tchapman,4,2)
            lyt.addWidget(QtGui.QLabel(u"rhoc"),1,4)
            self.rhoc=Entrada_con_unidades(float, value=eq["rhoc"], readOnly=True)
            lyt.addWidget(self.rhoc,1,5)
            lyt.addWidget(QtGui.QLabel(u"ff"),2,4)
            self.ff=Entrada_con_unidades(float, value=eq["ff"], readOnly=True)
            lyt.addWidget(self.ff,2,5)
            lyt.addWidget(QtGui.QLabel(u"rm"),3,4)
            self.rm=Entrada_con_unidades(float, value=eq["rm"], readOnly=True)
            lyt.addWidget(self.rm,3,5)

        if "critical" in eq and eq["critical"]:
            gridLayout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),5,3)
            gridLayout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Critical enhancement")),6,1,1,2)
            if eq["critical"] == 3:
                gridLayout.addWidget(QtGui.QLabel(u"gnu"),7,1)
                self.gnu=Entrada_con_unidades(float, value=eq["gnu"], readOnly=True)
                gridLayout.addWidget(self.gnu,7,2)
                gridLayout.addWidget(QtGui.QLabel(u"γ"),8,1)
                self.gamma=Entrada_con_unidades(float, value=eq["gamma"], readOnly=True)
                gridLayout.addWidget(self.gamma,8,2)
                gridLayout.addWidget(QtGui.QLabel(u"Ro"),9,1)
                self.R0=Entrada_con_unidades(float, value=eq["R0"], readOnly=True)
                gridLayout.addWidget(self.R0,9,2)
                gridLayout.addWidget(QtGui.QLabel(u"ξo"),10,1)
                self.Xio=Entrada_con_unidades(float, value=eq["Xio"], readOnly=True)
                gridLayout.addWidget(self.Xio,10,2)
                gridLayout.addWidget(QtGui.QLabel(u"Γo"),11,1)
                self.gam0=Entrada_con_unidades(float, value=eq["gam0"], readOnly=True)
                gridLayout.addWidget(self.gam0,11,2)
                gridLayout.addWidget(QtGui.QLabel(u"qd"),12,1)
                self.qd=Entrada_con_unidades(float, value=eq["qd"], readOnly=True)
                gridLayout.addWidget(self.qd,12,2)
            elif eq["critical"] == 4:
                gridLayout.addWidget(QtGui.QLabel(u"γ"),7,1)
                self.gamma=Entrada_con_unidades(float, value=eq["gamma"], readOnly=True)
                gridLayout.addWidget(self.gamma,7,2)
                gridLayout.addWidget(QtGui.QLabel("v"),8,1)
                self.v=Entrada_con_unidades(float, value=eq["expo"], readOnly=True)
                gridLayout.addWidget(self.v,8,2)
                gridLayout.addWidget(QtGui.QLabel(u"α"),9,1)
                self.alfa=Entrada_con_unidades(float, value=eq["alfa"], readOnly=True)
                gridLayout.addWidget(self.alfa,9,2)
                gridLayout.addWidget(QtGui.QLabel(u"β"),10,1)
                self.beta=Entrada_con_unidades(float, value=eq["beta"], readOnly=True)
                gridLayout.addWidget(self.beta,10,2)
                gridLayout.addWidget(QtGui.QLabel(u"Γo"),11,1)
                self.Xio=Entrada_con_unidades(float, value=eq["Xio"], readOnly=True)
                gridLayout.addWidget(self.Xio,11,2)
            gridLayout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),15,3)


class Ui_Properties(QtGui.QDialog):
    """Dialog for select and sort shown properties in tables"""
    _default=[1, 0, 1, 0, 0, 1, 0, 1, 1]+[0]*55
    def __init__(self, config=None, parent=None):
        super(Ui_Properties, self).__init__(parent)
        if config and config.has_option("MEoS", "properties"):
            values = config.get("MEoS", "properties")
            if isinstance(values, str):
                values=eval(values)
            fase=config.getboolean("MEoS", "phase")
            self.order = config.get("MEoS", "propertiesOrder")
            if isinstance(self.order, str):
                self.order = eval(self.order)
        else:
            values=self._default
            fase=False
            self.order = range(64)

        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Select Properties"))
        layout = QtGui.QGridLayout(self)
        self.listaDisponibles=QtGui.QTableWidget(len(meos.propiedades), 2)
        self.listaDisponibles.verticalHeader().hide()
        self.listaDisponibles.horizontalHeader().hide()
        self.listaDisponibles.horizontalHeader().setStretchLastSection(True)
        self.listaDisponibles.setGridStyle(QtCore.Qt.NoPen)
        self.listaDisponibles.setColumnWidth(0, 18)
        self.listaDisponibles.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
        self.listaDisponibles.setEditTriggers(QtGui.QAbstractItemView.NoEditTriggers)
        self.listaDisponibles.setItemDelegateForColumn(0, CheckEditor(self))
        for i in range(64):
            self.listaDisponibles.setItem(i, 0, QtGui.QTableWidgetItem(str(values[i])))
            self.listaDisponibles.setItem(i, 1, QtGui.QTableWidgetItem(meos.propiedades[self.order[i]]))
            self.listaDisponibles.setRowHeight(i, 20)
            self.listaDisponibles.openPersistentEditor(self.listaDisponibles.item(i, 0))
        self.listaDisponibles.currentCellChanged.connect(self.comprobarBotones)
        self.listaDisponibles.cellDoubleClicked.connect(self.toggleCheck)
        layout.addWidget(self.listaDisponibles,1,1,6,1)

        self.ButtonArriba=QtGui.QToolButton()
        self.ButtonArriba.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/arrow-up.png")))
        self.ButtonArriba.clicked.connect(self.Up)
        layout.addWidget(self.ButtonArriba, 3, 2, 1, 1)
        self.ButtonAbajo=QtGui.QToolButton()
        self.ButtonAbajo.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/arrow-down.png")))
        self.ButtonAbajo.clicked.connect(self.Down)
        layout.addWidget(self.ButtonAbajo, 4, 2, 1, 1)

        self.checkFase=QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Show bulk, liquid and vapor properties"))
        self.checkFase.setChecked(fase)
        layout.addWidget(self.checkFase,7,1,1,2)
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Reset|QtGui.QDialogButtonBox.Ok|QtGui.QDialogButtonBox.Cancel)
        self.buttonBox.clicked.connect(self.buttonClicked)
        layout.addWidget(self.buttonBox,8,1,1,2)

    def toggleCheck(self, fila, columna):
        """Toggle check status with a doubleclick in row"""
        txt=self.listaDisponibles.item(fila, 0).text()
        if txt == "0":
            newtxt = "1"
        else:
            newtxt = "0"
        self.listaDisponibles.item(fila, 0).setText(newtxt)

    def Down(self):
        """Change current selected row with next row"""
        i=self.listaDisponibles.currentRow()
        txt=self.listaDisponibles.item(i, 0).text()
        self.listaDisponibles.item(i, 0).setText(self.listaDisponibles.item(i+1, 0).text())
        self.listaDisponibles.item(i+1, 0).setText(txt)
        item=self.listaDisponibles.takeItem(i, 1)
        self.listaDisponibles.setItem(i, 1, self.listaDisponibles.takeItem(i+1, 1))
        self.listaDisponibles.setItem(i+1, 1, item)
        self.listaDisponibles.setCurrentCell(i+1, 0)
        self.order[i], self.order[i+1] = self.order[i+1], self.order[i]

    def Up(self):
        """Change current selected row with previous row"""
        i=self.listaDisponibles.currentRow()
        txt=self.listaDisponibles.item(i, 0).text()
        self.listaDisponibles.item(i, 0).setText(self.listaDisponibles.item(i-1, 0).text())
        self.listaDisponibles.item(i-1, 0).setText(txt)
        item=self.listaDisponibles.takeItem(i, 1)
        self.listaDisponibles.setItem(i, 1, self.listaDisponibles.takeItem(i-1, 1))
        self.listaDisponibles.setItem(i-1, 1, item)
        self.listaDisponibles.setCurrentCell(i-1, 0)
        self.order[i], self.order[i-1] = self.order[i-1], self.order[i]

    def buttonClicked(self, boton):
        """Actions for dialogbuttonbox functionality"""
        if self.buttonBox.buttonRole(boton)==QtGui.QDialogButtonBox.AcceptRole:
            self.accept()
        elif self.buttonBox.buttonRole(boton)==QtGui.QDialogButtonBox.RejectRole:
            self.reject()
        elif self.buttonBox.buttonRole(boton)==QtGui.QDialogButtonBox.ResetRole:
            for i, propiedad in enumerate(self._default):
                self.listaDisponibles.item(i, 0).setText(str(propiedad))
            self.checkFase.setChecked(False)

    @property
    def properties(self):
        """Properties list"""
        value=[]
        for i in range(self.listaDisponibles.rowCount()):
            value.append(self.listaDisponibles.cellWidget(i, 0).isChecked())
        return value

    def comprobarBotones(self, fila):
        """Check if button are enabled or disabled"""
        self.ButtonArriba.setEnabled(fila>=1)
        self.ButtonAbajo.setEnabled(fila<self.listaDisponibles.rowCount()-1)

# Table data
def createTabla(config, title, fluidos=None, parent=None):
    """Create TablaMEoS to add to mainwindow
        config: configparser instance with project configuration
        title: title for the table
        fluidos: optional array with meos instances to fill de table
        parent: mainwindow pointer
        """
    propiedades, keys, units =get_propiedades(config)

    for i, unit in enumerate(units):
        sufx = unit.text()
        if not sufx:
            sufx = "[-]"
        propiedades[i]=propiedades[i]+os.linesep+sufx

    # Add two phases properties if requested
    if config.getboolean("MEoS", "phase"):
        for i in range(len(propiedades)-1, -1, -1):
            if keys[i] in meos._fase.__dict__:
                txt = [propiedades[i]]
                prefix = QtGui.QApplication.translate("pychemqt", "Liquid")
                txt.append(prefix+os.linesep+propiedades[i])
                prefix = QtGui.QApplication.translate("pychemqt", "Vapour")
                txt.append(prefix+os.linesep+propiedades[i])
                propiedades[i:i+1] = txt
                units[i:i+1] = [units[i]]*3

    if fluidos:
        data=[]
        for fluido in fluidos:
            fila = _getData(fluido, keys, config)
            data.append(fila)

        tabla = TablaMEoS(len(propiedades), horizontalHeader=propiedades, stretch=False, readOnly=True, units=units, parent=parent)
        tabla.setMatrix(data)
    else:
        columnInput=[]
        for key in keys:
            if key in ["P", "T", "x", "rho", "v", "h", "s"]:
                columnInput.append(False)
            else:
                columnInput.append(True)
            if config.getboolean("MEoS", "phase") and key in meos._fase.__dict__:
                columnInput.append(True)
                columnInput.append(True)

        if config.getboolean("MEoS", "phase"):
            for i in range(len(keys)-1, -1, -1):
                if keys[i] in meos._fase.__dict__:
                    keys[i:i+1] = [keys[i], "", ""]

        tabla = TablaMEoS(len(propiedades), horizontalHeader=propiedades, filas=1, units=units, keys=keys, stretch=False, columnReadOnly=columnInput, parent=parent)

    prefix=QtGui.QApplication.translate("pychemqt", "Table")
    tabla.setWindowTitle(prefix+": "+title)
    tabla.resizeColumnsToContents()
    return tabla


def get_propiedades(config):
    booleanos = config.get("MEoS", "properties")
    order = config.get("MEoS", "propertiesOrder")
    if isinstance(booleanos, str):
        booleanos = eval(booleanos)
    if isinstance(order, str):
        order = eval(order)
        
    propiedades=[]
    keys=[]
    units = []
    for indice, bool in zip(order, booleanos):
        if bool:
            propiedades.append(meos.propiedades[indice])
            keys.append(meos.keys[indice])
            units.append(meos.units[indice])
    return propiedades, keys, units


def _getData(fluido, keys, config, unit=None):
    fila=[]
    for i, key in enumerate(keys):
        if key:
            p = fluido.__getattribute__(key)
            if p is None:
                txt = QtGui.QApplication.translate("pychemqt", "undefined")
            else:
                if unit and unit[i]:
                    txt = p.__getattribute__(unit[i])
                else:
                    txt = p.config()
            fila.append(txt)
            if config.getboolean("MEoS", "phase") and key in meos._fase.__dict__:
                p = fluido.Liquido.__getattribute__(key)
                if p is None:
                    txt = QtGui.QApplication.translate("pychemqt", "undefined")
                elif isinstance(p, unidades.unidad):
                    if unit and unit[i]:
                        txt = p.__getattribute__(unit[i])
                    else:
                        txt = p.config()
                else:
                    txt = p
                fila.append(txt)
                p = fluido.Gas.__getattribute__(key)
                if p is None:
                    txt = QtGui.QApplication.translate("pychemqt", "undefined")
                elif isinstance(p, unidades.unidad):
                    if unit and unit[i]:
                        txt = p.__getattribute__(unit[i])
                    else:
                        txt = p.config()
                else:
                    txt = p
                fila.append(txt)
    return fila


class TablaMEoS(Tabla):
    """Tabla subclass to show meos data, add context menu options"""
    Plot = None
    def __init__(self, *args, **kwargs):
        if "keys" in kwargs:
            self.keys = kwargs["keys"]
            del kwargs["keys"]
        self.units = kwargs["units"]
        del kwargs["units"]
        if "orderUnit" in kwargs:
            self.orderUnit = kwargs["orderUnit"]
            del kwargs["orderUnit"]
        else:
            self.orderUnit = []
            for unit in self.units:
                if unit == unidades.Dimensionless:
                    self.orderUnit.append(0)
                else:
                    self.orderUnit.append(unit.Config.getint('Units', unit.__name__))

        if "format" in kwargs:
            self.format = kwargs["format"]
            del kwargs["format"]
        else:
            self.format=[{"format": 1, "decimales": 6, "signo": False}]*args[0]

        super(TablaMEoS, self).__init__(*args, **kwargs)
        self.horizontalHeader().setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.horizontalHeader().customContextMenuRequested.connect(self.horizontalHeaderClicked)
        self.verticalHeader().setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.verticalHeader().customContextMenuRequested.connect(self.verticalHeaderClicked)
        self.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
        self.itemSelectionChanged.connect(self.selectPoint)
        self.data=[]
        self.parent=kwargs.get("parent", None)
        if not self.readOnly:
            self.cellChanged.connect(self.calculatePoint)

    def _getPlot(self):
        """Return plot if it's loaded"""
        if not self.Plot:
            for window in self.parent.centralwidget.currentWidget().subWindowList():
                widget = window.widget()
                if isinstance(widget, PlotMEoS):
                    self.Plot = widget
                    break
        return self.Plot

    def horizontalHeaderClicked(self, event):
        """Show dialog to config format and unit"""
        column = self.horizontalHeader().logicalIndexAt(event)
        unit = self.units[column]
        dialog=NumericFactor(self.format[column], unit, self.orderUnit[column])
        if dialog.exec_():
            # Check unit change
            if unit != unidades.Dimensionless and \
                dialog.unit.currentIndex() != self.orderUnit[column]:
                for i, fila in enumerate(self.data):
                    value = unit(fila[column], unit.__units__[self.orderUnit[column]]).__getattribute__(unit.__units__[dialog.unit.currentIndex()])
                    self.data[i][column] = value
                self.orderUnit[column] = dialog.unit.currentIndex()
                txt = self.horizontalHeaderItem(column).text().split(os.linesep)[0]
                txt += os.linesep+unit.__text__[dialog.unit.currentIndex()]
                self.setHorizontalHeaderItem(column,QtGui.QTableWidgetItem(txt))

            # Check format change
            self.format[column]=dialog.args()
            self.setStr()
            self.resizeColumnToContents(column)
        self.setRangeSelected(QtGui.QTableWidgetSelectionRange(0, column, self.rowCount()-1, column), True)

    def verticalHeaderClicked(self, position):
        row = self.verticalHeader().logicalIndexAt(position)
        rows = []
        for item in self.selectedItems():
            if item.row() not in rows:
                rows.append(item.row())
        rows.sort(reverse=True)

        deleteAction=createAction(QtGui.QApplication.translate("pychemqt", "Delete Point"), icon=os.environ["pychemqt"]+"/images/button/editDelete", slot=partial(self.delete, rows), parent=self)
        inserPoint=createAction(QtGui.QApplication.translate("pychemqt", "Insert Point"), icon=os.environ["pychemqt"]+"/images/button/add", slot=partial(self.add, row), parent=self)

        menu = QtGui.QMenu()
        menu.addAction(deleteAction)
        menu.addSeparator()
        menu.addAction(inserPoint)
        menu.exec_(self.mapToGlobal(position))

    def delete(self, rows):
        self.parent.statusbar.showMessage(QtGui.QApplication.translate(
            "pychemqt", "Deleting point..."))
        QtGui.QApplication.processEvents()

        # Delete point from table
        for row in rows:
            self.removeRow(row)
            delete(self.data, row)

        # Delete point from data plot
        plot =  self._getPlot()
        data = plot._getData()
        title = self.windowTitle().split(QtGui.QApplication.translate("pychemqt", "Table from "))[1]
        for row in rows:
            if title == QtGui.QApplication.translate("pychemqt", "Melting Line"):
                for x in prop_pickle:
                    del data["melting"][x][row]
            elif title == QtGui.QApplication.translate("pychemqt", "Sublimation Line"):
                for x in prop_pickle:
                    del data["sublimation"][x][row]
            elif title == QtGui.QApplication.translate("pychemqt", "Saturation Line") or \
                    title == QtGui.QApplication.translate("pychemqt", "Liquid Saturation Line"):
                for x in prop_pickle:
                    del data["saturation_0"][x][row]
            elif title == QtGui.QApplication.translate("pychemqt", "Vapor Saturation Line"):
                for x in prop_pickle:
                    del data["saturation_1"][x][row]
            else:
                units = {"P": unidades.Pressure,
                         "T": unidades.Temperature,
                         "h": unidades.Enthalpy,
                         "s": unidades.Enthalpy,
                         "v": unidades.SpecificVolume,
                         "rho": unidades.Density}
                var = str(title.split(" = ")[0])
                txt = title.split(" = ")[1]
                unit = units[var]
                value = float(txt.split(" ")[0])
                stdValue = unit(value, "conf")
                for x in prop_pickle:
                    del data[var][stdValue][x][row]
        plot._saveData(data)

        # Delete point from data
        for line in plot.plot.ax.lines:
            if unicode(line.get_label()) == unicode(title):
                xdata = line._x
                ydata = line._y
                for row in rows:
                    xdata = delete(xdata, row)
                    ydata = delete(ydata, row)
                line.set_xdata(xdata)
                line.set_ydata(ydata)
                plot.plot.draw()
                break
        self.parent.statusbar.clearMessage()

    def add(self, row):
        dialog = AddPoint(self.parent.currentConfig, self.parent)
        if dialog.exec_():
            if dialog.checkBelow.isChecked():
                row += 1
            plot = self.Plot
            if plot is None:
                plot =  self._getPlot()

            datatoTable = []
            datatoTable.append(dialog.fluid.__getattribute__(plot.x).config())
            datatoTable.append(dialog.fluid.__getattribute__(plot.y).config())

            # Add point to table
            self.addRow(datatoTable, row)
            self.data=insert(self.data, row, datatoTable)

            # Add point to data plot
            data = plot._getData()
            title = self.windowTitle().split(QtGui.QApplication.translate("pychemqt", "Table from "))[1]
            if title == QtGui.QApplication.translate("pychemqt", "Melting Line"):
                for x in prop_pickle:
                    data["melting"][x].insert(row, dialog.fluid.__getattribute__(x))
            elif title == QtGui.QApplication.translate("pychemqt", "Sublimation Line"):
                for x in prop_pickle:
                    data["sublimation"].insert(row, dialog.fluid.__getattribute__(x))
            elif title == QtGui.QApplication.translate("pychemqt", "Saturation Line") or \
                    title == QtGui.QApplication.translate("pychemqt", "Liquid Saturation Line"):
                for x in prop_pickle:
                    data["saturation_0"].insert(row, dialog.fluid.__getattribute__(x))
            elif title == QtGui.QApplication.translate("pychemqt", "Vapor Saturation Line"):
                for x in prop_pickle:
                    data["saturation_1"].insert(row, dialog.fluid.__getattribute__(x))
            else:
                units = {"P": unidades.Pressure,
                         "T": unidades.Temperature,
                         "h": unidades.Enthalpy,
                         "s": unidades.Enthalpy,
                         "v": unidades.SpecificVolume,
                         "rho": unidades.Density}
                var = str(title.split(" = ")[0])
                txt = title.split(" = ")[1]
                unit = units[var]
                value = float(txt.split(" ")[0])
                stdValue = unit(value, "conf")

                for x in prop_pickle:
                    data[var][stdValue][x].insert(row, dialog.fluid.__getattribute__(x))
            plot._saveData(data)

            # Add point to data
            for line in plot.plot.ax.lines:
                if unicode(line.get_label()) == unicode(title):
                    xdata = insert(line._x, row, datatoTable[0])
                    ydata = insert(line._y, row, datatoTable[1])
                    line.set_xdata(xdata)
                    line.set_ydata(ydata)
                    plot.plot.draw()
                    break


    def selectPoint(self):
        plot = self._getPlot()
        if plot:
        #Remove old selected point if exist
            for i, line in enumerate(plot.plot.ax.lines):
                if line.get_label() == QtGui.QApplication.translate("pychemqt", "Selected Point"):
                    del line
                    del plot.plot.ax.lines[i]

        #Add new selected points
            x = []
            y = []
            for item in self.selectedItems():
                if item.column():
                    y.append(float(item.text()))
                else:
                    x.append(float(item.text()))
            label = QtGui.QApplication.translate("pychemqt", "Selected Point")
            plot.plot.ax.plot(x, y, 'ro', label=label)
            plot.plot.draw()



    def calculatePoint(self, row, column):
        """Add new value to kwargs for point, and show properties if is calculable
        row, column: index for modified cell in table"""
        key = self.keys[column]
        unit = self.units[column]
        if unit is unidades.Dimensionless:
            value = float(self.item(row, column).text())
        else:
            value = unit(float(self.item(row, column).text()), unit.__units__[self.orderUnit[column]])
        self.Point(**{key: value})
        if self.Point.status:
            units = []
            for ui, order in zip(self.units, self.orderUnit):
                if ui is unidades.Dimensionless:
                    units.append("")
                else:
                    units.append(ui.__units__[order])
            data = _getData(self.Point, self.keys, self.parent.currentConfig,
                            units)
            self.setRow(row, data)
            self.Point = self.fluid()


            self.addRow()
            self.setCurrentCell(row+1, column)

    def setMatrix(self, data):
        if self.readOnly:
            self.data=data
            self.setStr()
        else:
            for i, row in enumerate(data):
                self.setRow(i, row)
        self.resizeColumnsToContents()

    def setStr(self):
        for fila, array in enumerate(self.data):
            if fila>=self.rowCount():
                self.addRow()
            for columna, data in enumerate(array):
                if isinstance(data, QtCore.QString):
                    txt = data
                else:
                    txt=representacion(data, **self.format[columna])
                self.setValue(fila, columna, txt)

    def setRow(self, row, data):
        self.blockSignals(True)
        self.data.append(data)
        for column, data in enumerate(data):
            if isinstance(data, QtCore.QString):
                txt = data
            else:
                txt=representacion(data, **self.format[column])
            self.setValue(row, column, txt)
        self.resizeColumnsToContents()

        # Set calculate point readOnly
        if not self.readOnly:
            flags=QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsSelectable
            color=self.parent.Preferences.get("General", 'Color_ReadOnly')
            for i, bool in enumerate(self.columnReadOnly):
                if not bool:
                    self.item(row, i).setFlags(flags)
                    self.item(row, i).setBackground(QtGui.QColor(color))
        self.blockSignals(False)

    def contextMenuEvent(self, event):
        """Show context menu over cell"""
        menu = QtGui.QMenu()
        actionCopy = createAction(
            QtGui.QApplication.translate("pychemqt", "&Copy"),
            slot=partial(self.copy, event), shortcut=QtGui.QKeySequence.Copy,
            icon=os.environ["pychemqt"]+"/images/button/editCopy", parent=self)
        export = createAction(
            QtGui.QApplication.translate("pychemqt", "E&xport to csv"),
            self.exportCSV, icon=os.environ["pychemqt"]+"/images/button/export",
            tip=QtGui.QApplication.translate("pychemqt", "Export table to file"),
            parent=self)
        menu.addAction(actionCopy)
        menu.addSeparator()
        menu.addAction(export)
        menu.exec_(event.globalPos())

    def copy(self, event):
        """Copy selected value to clipboard"""
        widget=self.itemAt(self.viewport().mapFromGlobal(event.globalPos()))
        QtGui.QApplication.clipboard().setText(widget.text())

    def exportCSV(self):
        """Export data table as a csv file"""
        if self.parent.currentFilename:
            dir = os.path.dirname(str(self.parent.currentFilename))
        else:
            dir = "."

        pat=QtCore.QStringList()
        pat.append(QtGui.QApplication.translate("pychemqt", "CSV files") + " (*.csv)")
        if os.environ["ezodf"]:
            pat.append(QtGui.QApplication.translate(
            "pychemqt", "Libreoffice spreadsheet files")+ " (*.ods)")
        if os.environ["xlwt"]:
            pat.append(QtGui.QApplication.translate(
            "pychemqt", "Microsoft Excel 97/2000/XP/2003 XML")+ " (*.xls)")
        if os.environ["openpyxl"]:
            pat.append(QtGui.QApplication.translate(
            "pychemqt", "Microsoft Excel 2007/2010 XML")+ " (*.xlsx)")
        patron=pat.join(";;")

        dialog = QtGui.QFileDialog(self,
            QtGui.QApplication.translate("pychemqt", "Export table to file"),
            dir, patron)
        dialog.setAcceptMode(QtGui.QFileDialog.AcceptOpen)
        dialog.setFileMode(QtGui.QFileDialog.AnyFile)
        if dialog.exec_():
            fname = unicode(dialog.selectedFiles().join(""))
            ext = unicode(dialog.selectedNameFilter().split(".")[-1][:-1])
            
            if fname:
                exportTable(self.data, fname, ext, self.encabezadoHorizontal)

    def writeToStream(self, stream):
        stream.writeInt32(self.columnCount())

        # Save titles
        stream.writeString(self.windowTitle())
        for col in range(self.columnCount()):
            stream.writeString(self.horizontalHeaderItem(col).text())

        # Save units as index
        all = unidades._all
        all.append(unidades.Dimensionless)
        for unit in self.units:
            stream.writeInt32(all.index(unit))

        # Save keys if necesary
        stream.writeBool(self.readOnly)
        if not self.readOnly:
            stream.writeInt32(mEoS.__all__.index(self.fluid))
            for key in self.keys:
                stream.writeString(key)
            for boolean in self.columnReadOnly:
                stream.writeBool(boolean)

        # Save order unit
        for index in self.orderUnit:
            stream.writeInt32(index)

        # Save format
        for format in self.format:
            stream.writeInt32(format.get("format", 1))
            stream.writeInt32(format.get("decimales", 6))
            stream.writeBool(format.get("signo", False))
            stream.writeInt32(format.get("total", 10))
            stream.writeBool(format.get("exp", False))
            stream.writeInt32(format.get("tol", 6))
            stream.writeBool(format.get("thousand", False))

        # Save data if exist
        stream.writeInt32(len(self.data))
        for row in self.data:
            stream.writeInt32(len(row))
            for data in row:
                boolean = isinstance(data, str) or isinstance(data, QtCore.QString)
                stream.writeBool(boolean)
                if boolean:
                    stream.writeString(data)
                else:
                    stream.writeFloat(data)

    @classmethod
    def readFromStream(cls, stream, parent):
        columnCount = stream.readInt32()

        # Get titles
        title = stream.readString()
        propiedades = []
        for col in range(columnCount):
            propiedades.append(stream.readString())

        # Get units
        all = unidades._all
        all.append(unidades.Dimensionless)
        units = []
        for i in range(columnCount):
            index = stream.readInt32()
            units.append(all[index])

        # Get keys if neccesary
        readOnly = stream.readBool()
        if not readOnly:
            index = stream.readInt32()
            fluid = mEoS.__all__[index]
            keys = []
            for i in range(columnCount):
                keys.append(stream.readString())
            columnReadOnly = []
            for i in range(columnCount):
                columnReadOnly.append(stream.readBool())

        # Get OrderUnit
        orderUnit = []
        for i in range(columnCount):
            orderUnit.append(stream.readInt32())

        # Get format
        format = []
        for col in range(columnCount):
            fr = {}
            fr["format"] = stream.readInt32()
            fr["decimales"] = stream.readInt32()
            fr["signo"] = stream.readBool()
            fr["total"] = stream.readInt32()
            fr["exp"] = stream.readBool()
            fr["tol"] = stream.readInt32()
            fr["thousand"] = stream.readBool()
            format.append(fr)

        # Get data
        data = []
        datalen = stream.readInt32()
        for row in range(datalen):
            rowlen = stream.readInt32()
            lista = []
            for element in range(rowlen):
                boolean = stream.readBool()
                if boolean:
                    lista.append(stream.readString())
                else:
                    lista.append(stream.readFloat())
            data.append(lista)

        #Create Tabla
        args = (columnCount, )
        kwargs = {}
        kwargs["horizontalHeader"] = propiedades
        kwargs["format"] = format
        kwargs["stretch"] = False
        kwargs["parent"] = parent
        kwargs["units"] = units
        kwargs["orderUnit"] = orderUnit

        if readOnly:
            kwargs["readOnly"] = True
        else:
            kwargs["filas"] = datalen+1
            kwargs["keys"] = keys
            kwargs["columnReadOnly"] = columnReadOnly

        tabla = TablaMEoS(*args, **kwargs)
        tabla.setWindowTitle(title)
        tabla.setMatrix(data)
        if not readOnly:
            tabla.fluid = fluid
            tabla.Point = fluid()
        return tabla


class Ui_Saturation(QtGui.QDialog):
    """Dialog to define input for a two-phase saturation table calculation"""
    def __init__(self, config=None, parent=None):
        """config: instance with project config to set initial values"""
        super(Ui_Saturation, self).__init__(parent)
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Saturation Table"))
        layout = QtGui.QGridLayout(self)

        groupboxTypo=QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Interphase"))
        layout.addWidget(groupboxTypo,1,1,1,2)
        layoutg1=QtGui.QGridLayout(groupboxTypo)
        self.VL=QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "Vapor-Liquid (boiling line)" ))
        layoutg1.addWidget(self.VL,1,1)
        self.SL=QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "Solid-Liquid (melting line"))
        layoutg1.addWidget(self.SL,2,1)
        self.SV=QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "Solid-Vapor (Sublimation line)" ))
        layoutg1.addWidget(self.SV,3,1)

        groupboxVariar=QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Change"))
        layout.addWidget(groupboxVariar,1,3,1,2)
        layoutg2=QtGui.QGridLayout(groupboxVariar)
        self.VariarTemperatura=QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "Temperature"))
        self.VariarTemperatura.toggled.connect(self.updateVar)
        layoutg2.addWidget(self.VariarTemperatura,1,1)
        self.VariarPresion=QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "Pressure"))
        self.VariarPresion.toggled.connect(self.updateVar)
        layoutg2.addWidget(self.VariarPresion,2,1)
        self.VariarXconT=QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "Quality at fixed temperature"))
        self.VariarXconT.toggled.connect(self.updateVar)
        layoutg2.addWidget(self.VariarXconT,3,1)
        self.VariarXconP=QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "Quality at fixed pressure"))
        self.VariarXconP.toggled.connect(self.updateVar)
        layoutg2.addWidget(self.VariarXconP,4,1)

        line = QtGui.QFrame()
        line.setFrameShape(QtGui.QFrame.HLine)
        line.setFrameShadow(QtGui.QFrame.Sunken)
        layout.addWidget(line,2,1,1,4)

        self.labelFix=QtGui.QLabel()
        layout.addWidget(self.labelFix,4,3)
        self.variableFix=Entrada_con_unidades(float)
        layout.addWidget(self.variableFix,4,4)
        self.labelinicial=QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Initial"))
        layout.addWidget(self.labelinicial,4,1)
        self.Inicial=Entrada_con_unidades(float)
        layout.addWidget(self.Inicial,4,2)
        self.labelfinal=QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Final"))
        layout.addWidget(self.labelfinal,5,1)
        self.Final=Entrada_con_unidades(float)
        layout.addWidget(self.Final,5,2)
        self.labelincremento=QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Increment"))
        layout.addWidget(self.labelincremento,6,1)
        self.Incremento=Entrada_con_unidades(float)
        layout.addWidget(self.Incremento,6,2)

        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok|QtGui.QDialogButtonBox.Cancel)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox,10,1,1,4)

        if config:
            self.fluido=mEoS.__all__[config.getint("MEoS", "fluid")]
            if self.fluido._melting or self.fluido._Melting_Pressure != meos.MEoS._Melting_Pressure:
                self.SL.setEnabled(True)
            else:
                self.SL.setEnabled(False)
            if self.fluido._sublimation or self.fluido._Sublimation_Pressure != meos.MEoS._Sublimation_Pressure:
                self.SV.setEnabled(True)
            else:
                self.SV.setEnabled(False)

        self.VL.setChecked(True)
        self.VariarTemperatura.setChecked(True)
        self.updateVary()
        self.VL.toggled.connect(self.updateVary)


    def updateVary(self):
        """Update state for option to choose for properties to change"""
        self.VariarXconP.setEnabled(self.VL.isChecked())
        self.VariarXconT.setEnabled(self.VL.isChecked())
        self.VariarTemperatura.setChecked(not self.VL.isChecked())

    def updateVar(self, bool):
        """Update input values units and text"""
        if bool:
            # Select initial values
            fix, inicial, final, step = 0, 0, 0, 0
            if self.VL.isChecked():
                if self.sender() == self.VariarXconT:
                    fix = ceil((self.fluido.Tc-self.fluido.Tt)/2)
                    inicial = 0
                    final = 1
                    step = 0.1
                elif self.sender() == self.VariarXconP:
                    fix = ceil(self.fluido.Pc/2)
                    inicial = 0
                    final = 1
                    step = 0.1
                elif self.sender() == self.VariarTemperatura:
                    inicial = ceil(self.fluido.Tt)
                    final = floor(self.fluido.Tc)
                    step = 1.

            self.Inicial.deleteLater()
            self.Final.deleteLater()
            self.Incremento.deleteLater()
            if self.sender() == self.VariarXconT:
                self.labelFix.setVisible(True)
                self.labelFix.setText(unidades.Temperature.__title__)
                self.variableFix.deleteLater()
                self.variableFix=Entrada_con_unidades(unidades.Temperature, value=fix)
                self.layout().addWidget(self.variableFix,4,4)
                unidadVariable=float
                self.labelinicial.setText(QtGui.QApplication.translate("pychemqt", "Initial quality"))
                self.labelfinal.setText(QtGui.QApplication.translate("pychemqt", "Final quality"))

            elif self.sender() == self.VariarXconP:
                self.labelFix.setVisible(True)
                self.labelFix.setText(unidades.Pressure.__title__)
                self.variableFix.deleteLater()
                self.variableFix=Entrada_con_unidades(unidades.Pressure, value=fix)
                self.layout().addWidget(self.variableFix,4,4)
                unidadVariable=float
                self.labelinicial.setText(QtGui.QApplication.translate("pychemqt", "Initial quality"))
                self.labelfinal.setText(QtGui.QApplication.translate("pychemqt", "Final quality"))

            elif self.sender() == self.VariarTemperatura:
                self.labelFix.setVisible(False)
                self.variableFix.setVisible(False)
                unidadVariable=unidades.Temperature
                self.labelinicial.setText(QtGui.QApplication.translate("pychemqt", "Initial temperature"))
                self.labelfinal.setText(QtGui.QApplication.translate("pychemqt", "Final temperature"))

            else:
                self.labelFix.setVisible(False)
                self.variableFix.setVisible(False)
                unidadVariable=unidades.Pressure
                self.labelinicial.setText(QtGui.QApplication.translate("pychemqt", "Initial pressure"))
                self.labelfinal.setText(QtGui.QApplication.translate("pychemqt", "Final pressure"))

            self.Inicial=Entrada_con_unidades(unidadVariable, value=inicial)
            self.Final=Entrada_con_unidades(unidadVariable, value=final)
            if unidadVariable==unidades.Temperature:
                unidadDelta=unidades.DeltaT
            elif unidadVariable==unidades.Pressure:
                unidadDelta=unidades.DeltaP
            else:
                unidadDelta=unidadVariable

            self.Incremento=Entrada_con_unidades(unidadDelta, value=step)
            self.layout().addWidget(self.Inicial,4,2)
            self.layout().addWidget(self.Final,5,2)
            self.layout().addWidget(self.Incremento,6,2)


class Ui_Isoproperty(QtGui.QDialog):
    """Dialog to define input for isoproperty table calculations"""
    propiedades = [QtGui.QApplication.translate("pychemqt", "Temperature"),
                   QtGui.QApplication.translate("pychemqt", "Pressure"),
                   QtGui.QApplication.translate("pychemqt", "Density"),
                   QtGui.QApplication.translate("pychemqt", "Volume"),
                   QtGui.QApplication.translate("pychemqt", "Enthalpy"),
                   QtGui.QApplication.translate("pychemqt", "Entropy"),
                   QtGui.QApplication.translate("pychemqt", "Internal Energy")]
    unidades = [unidades.Temperature, unidades.Pressure, unidades.Density,
                unidades.SpecificVolume, unidades.Enthalpy,
                unidades.SpecificHeat, unidades.Enthalpy, float]
    keys = ["T", "P", "rho", "v", "h", "s", "u", "x"]

    def __init__(self, parent=None):
        super(Ui_Isoproperty, self).__init__(parent)
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Specify Isoproperty Table"))
        layout = QtGui.QGridLayout(self)

        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Hold constant")),1,1)
        self.fix=QtGui.QComboBox()
        for propiedad in self.propiedades:
            self.fix.addItem(propiedad)
        self.fix.currentIndexChanged.connect(self.actualizarUI)
        layout.addWidget(self.fix,1,2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Vary")),2,1)
        self.vary=QtGui.QComboBox()
        self.vary.currentIndexChanged.connect(self.actualizarVariable)
        layout.addWidget(self.vary,2,2)

        line = QtGui.QFrame()
        line.setFrameShape(QtGui.QFrame.HLine)
        line.setFrameShadow(QtGui.QFrame.Sunken)
        layout.addWidget(line,3,1,1,2)

        self.labelFix=QtGui.QLabel()
        layout.addWidget(self.labelFix,4,1)
        self.variableFix=Entrada_con_unidades(float)
        layout.addWidget(self.variableFix,4,2)
        self.labelinicial=QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Initial"))
        layout.addWidget(self.labelinicial,5,1)
        self.Inicial=Entrada_con_unidades(float)
        layout.addWidget(self.Inicial,5,2)
        self.labelfinal=QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Final"))
        layout.addWidget(self.labelfinal,6,1)
        self.Final=Entrada_con_unidades(float)
        layout.addWidget(self.Final,6,2)
        self.labelincremento=QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Increment"))
        layout.addWidget(self.labelincremento,7,1)
        self.Incremento=Entrada_con_unidades(float)
        layout.addWidget(self.Incremento,7,2)

        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok|QtGui.QDialogButtonBox.Cancel)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox,10,1,1,2)

        self.actualizarUI(0)

    def actualizarUI(self, indice):
        self.vary.clear()
        propiedades=self.propiedades[:]
        if indice <= 1:
            propiedades.append(QtGui.QApplication.translate("pychemqt", "Quality"))
        del propiedades[indice]
        for propiedad in propiedades:
            self.vary.addItem(propiedad)
        self.labelFix.setText(self.propiedades[indice])
        self.variableFix.deleteLater()
        self.variableFix=Entrada_con_unidades(self.unidades[indice])
        self.layout().addWidget(self.variableFix,4,2)

    def actualizarVariable(self, indice):
        self.Inicial.deleteLater()
        self.Final.deleteLater()
        self.Incremento.deleteLater()
        if indice>=self.fix.currentIndex():
            indice+=1
        self.Inicial=Entrada_con_unidades(self.unidades[indice])
        self.Final=Entrada_con_unidades(self.unidades[indice])
        self.Incremento=Entrada_con_unidades(self.unidades[indice])
        self.layout().addWidget(self.Inicial,5,2)
        self.layout().addWidget(self.Final,6,2)
        self.layout().addWidget(self.Incremento,7,2)

class AddPoint(QtGui.QDialog):
    """Dialog to add new point to line2D"""
    keys = ["T", "P", "x", "rho", "v", "h", "s", "u"]
    def __init__(self, config, parent=None):
        """config: configParser Instance with currentproject configuration"""
        super(AddPoint, self).__init__(parent)
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Add Point to line"))
        layout = QtGui.QGridLayout(self)
        fluid = mEoS.__all__[config.getint("MEoS", "fluid")]
        self.fluid = fluid()

        self.Inputs = []
        for i, (title, key, unit) in enumerate(meos.inputData):
            layout.addWidget(QtGui.QLabel(title),i,1)
            if unit is unidades.Dimensionless:
                entrada = Entrada_con_unidades(float)
            else:
                entrada = Entrada_con_unidades(unit)
            entrada.valueChanged.connect(partial(self.update, key))
            self.Inputs.append(entrada)
            layout.addWidget(entrada,i,2)

        self.status = Status(self.fluid.status, self.fluid.msg)
        layout.addWidget(self.status,i+1,1,1,2)

        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "To")),i+2,1)
        self.To = Entrada_con_unidades(unidades.Temperature)
        self.To.valueChanged.connect(partial(self.update, "To"))
        layout.addWidget(self.To,i+2,2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "rhoo")),i+3,1)
        self.rhoo = Entrada_con_unidades(unidades.Density)
        self.rhoo.valueChanged.connect(partial(self.update, "rhoo"))
        layout.addWidget(self.rhoo,i+3,2)

        self.checkBelow = QtGui.QCheckBox(
            QtGui.QApplication.translate("pychemqt", "Add below selected point"))
        layout.addWidget(self.checkBelow,i+4,1,1,2)

        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Reset|QtGui.QDialogButtonBox.Ok|QtGui.QDialogButtonBox.Cancel)
        self.buttonBox.clicked.connect(self.click)
        layout.addWidget(self.buttonBox,i+5,1,1,2)

    def click(self, button):
        """Manage mouse click event over buttonbox"""
        if QtGui.QDialogButtonBox.Reset == self.buttonBox.standardButton(button):
            self.reset()
        elif QtGui.QDialogButtonBox.Ok == self.buttonBox.standardButton(button):
            self.accept()
        elif QtGui.QDialogButtonBox.Cancel == self.buttonBox.standardButton(button):
            self.reject()

    def update(self, key, value):
        """Update fluid instance with new parameter key with value"""
        self.status.setState(4)
        QtGui.QApplication.processEvents()
        self.fluid(**{key: value})
        if self.fluid.status in (1, 3):
            self.fill(self.fluid)
        self.status.setState(self.fluid.status, self.fluid.msg)

    def fill(self, fluid):
        """Fill dialog widget with fluid properties values"""
        self.blockSignals(True)
        Config=ConfigParser()
        Config.read(config.conf_dir+"pychemqtrc")
        for key, input in zip(self.keys, self.Inputs):
            input.setValue(fluid.__getattribute__(key))
            if fluid.kwargs[key]:
                input.setResaltado(True)
            else:
                input.setResaltado(False)
        self.blockSignals(False)

    def reset(self):
        """Reset dialog widgets to initial clear status"""
        self.fluid = self.fluid.__class__()
        self.status.setState(self.fluid.status, self.fluid.msg)
        self.rhoo.clear()
        self.To.clear()
        for input in self.Inputs:
            input.clear()
            input.setResaltado(False)


# Plot data
class PlotMEoS(QtGui.QWidget):
    """Plot widget to show meos plot data, add context menu options"""
    def __init__(self, dim, toolbar=False, filename="", parent=None):
        super(PlotMEoS, self).__init__(parent)
        self.parent=parent
        self.dim=dim
        self.filename = filename
        self.notes = []

        layout=QtGui.QVBoxLayout(self)
        self.plot=plot.matplotlib(dim)

        format = {}
        self.plot.lx = self.plot.ax.axhline(c="#888888", ls=":")  # the horiz line
        self.plot.ly = self.plot.ax.axvline(c="#888888", ls=":")  # the vert line

        self.plot.lx.set_visible(False)
        self.plot.ly.set_visible(False)

        layout.addWidget(self.plot)
        self.toolbar=plot.NavigationToolbar2QT(self.plot, self.plot)
        self.toolbar.setVisible(toolbar)
        layout.addWidget(self.toolbar)

        self.editAxesAction=createAction(QtGui.QApplication.translate("pychemqt", "Edit &Axis"), icon=os.environ["pychemqt"]+"/images/button/editor", slot=self.editAxis, parent=self)
        self.editAction=createAction(QtGui.QApplication.translate("pychemqt", "Edit &Plot"), slot=self.edit, icon=os.environ["pychemqt"]+"/images/button/Regression", parent=self)
        self.editMarginAction=createAction(QtGui.QApplication.translate("pychemqt", "Edit &Margins"), slot=self.toolbar.configure_subplots, parent=self)
        self.saveAction=createAction(QtGui.QApplication.translate("pychemqt", "&Save Plot"), slot=self.toolbar.save_figure, icon=os.environ["pychemqt"]+"/images/button/fileSave", parent=self)
        self.toolbarVisibleAction = createAction(QtGui.QApplication.translate("pychemqt", "Toggle &Toolbar"), self.toolbar.setVisible, checkable=True, parent=self)
        self.gridToggleAction = createAction(QtGui.QApplication.translate("pychemqt", "Toggle &Grid"), self.grid, checkable=True, parent=self)
        self.gridToggleAction.setChecked(self.parent.Preferences.getboolean("MEOS", "grid"))

        if dim==2:
            self.plot.fig.canvas.mpl_connect('button_press_event', self.click)
        else:
            self.editMarginAction.setEnabled(False)
            
    def contextMenuEvent(self, event):
        menuTable=QtGui.QMenu(QtGui.QApplication.translate("pychemqt", "Tabulated data"))
        menuTable.setIcon(QtGui.QIcon(os.environ["pychemqt"]+"/images/button/table"))
        for linea in self.plot.ax.lines:
            action=createAction(linea.get_label(), slot=partial(self.table, linea), parent=self)
            menuTable.addAction(action)

        menu = QtGui.QMenu()
        menu.addAction(self.editAxesAction)
        menu.addAction(self.editAction)
        menu.addAction(self.editMarginAction)
        menu.addSeparator()
        menu.addAction(self.saveAction)
        menu.addAction(menuTable.menuAction())
        menu.addSeparator()
        menu.addAction(self.toolbarVisibleAction)
        menu.addAction(self.gridToggleAction)
        menu.exec_(event.globalPos())

        if self.plot.ax._gridOn:
            self.gridToggleAction.setChecked(True)

    def grid(self, bool):
        self.plot.ax.grid(bool)
        self.plot.ax._gridOn = bool
        self.plot.draw()

    def edit(self):
        dialog=EditPlot(self, self.parent)
        dialog.show()

    def editAxis(self):
        dialog=EditAxis(self.plot)
        dialog.exec_()

    def table(self, obj):
        """
        Export plot data to table
        obj: object (Line2D instance) to show data"""
        xtxt = meos.propiedades[meos.keys.index(self.x)]
        ytxt = meos.propiedades[meos.keys.index(self.y)]
        xunit = meos.units[meos.keys.index(self.x)]
        yunit = meos.units[meos.keys.index(self.y)]
        HHeader=[xtxt+os.linesep+xunit.text(), ytxt+os.linesep+yunit.text()]
        units=[xunit, yunit]
        if self.dim == 3:
            ztxt = meos.propiedades[meos.keys.index(self.z)]
            zunit = meos.units[meos.keys.index(self.z)]
            HHeader.append(ztxt+os.linesep+zunit.text())
            units.append(zunit)
            data = obj._verts3d
        else:
            data = obj.get_data(orig=True)
            
        tabla = TablaMEoS(self.dim, horizontalHeader=HHeader, units=units,
                          stretch=False, readOnly=True, parent=self.parent)
        tabla.setMatrix(transpose(data))
        tabla.verticalHeader().setContextMenuPolicy(QtCore.Qt.CustomContextMenu)

        title = QtGui.QApplication.translate("pychemqt", "Table from") + " " + \
            obj.get_label()
        tabla.setWindowTitle(title)
        self.parent.centralwidget.currentWidget().addSubWindow(tabla)
        tabla.show()

    def _getData(self):
        """Get data from file"""
        filenameHard = os.environ["pychemqt"]+"dat"+os.sep+"mEoS"+ \
                       os.sep+self.filename+".gz"
        filenameSoft = config.conf_dir+self.filename
        file = ""
        if os.path.isfile(filenameSoft):
            with open(filenameSoft) as archivo:
                data = cPickle.load(archivo)
            return data
        elif os.path.isfile(filenameHard):
            with gzip.GzipFile(filenameHard, 'rb') as archivo:
                data = cPickle.load(archivo)
            self._saveData(data)
            return data

    def _saveData(self, data):
        """Save changes in data to file"""
        with open(config.conf_dir+self.filename, 'wb') as file:
            cPickle.dump(data, file)

    def click(self, event):
        """Update input and graph annotate when mouse click over chart"""
        # Accept only left click
        if event.button != 1:
            return
        units = {"x": unidades.Dimensionless,
                 "T": unidades.Temperature,
                 "P": unidades.Pressure,
                 "h": unidades.Enthalpy,
                 "u": unidades.Enthalpy,
                 "s": unidades.SpecificHeat,
                 "v": unidades.SpecificVolume,
                 "rho": unidades.Density}
        if self.x in units and self.y in units:
            x = units[self.x](event.xdata, "conf")
            y = units[self.y](event.ydata, "conf")

            fluid = mEoS.__all__[self.config["fluid"]]
            kwargs = {self.x: x, self.y: y}
            fluido = calcPoint(fluid, self.config, **kwargs)
            if fluido and fluido.status and \
                    fluido._constants["Tmin"] <= fluido.T <= fluido._constants["Tmax"] and \
                    fluido._constants["Pmin"] <= fluido.P.kPa <= fluido._constants["Pmax"]:
                self.plot.lx.set_ydata(event.ydata)
                self.plot.ly.set_xdata(event.xdata)
                self.plot.lx.set_visible(True)
                self.plot.ly.set_visible(True)
                self.showPointData(fluido)
            else:
                self.plot.lx.set_visible(False)
                self.plot.ly.set_visible(False)
                self.clearPointData()

    def showPointData(self, state):
        self.clearPointData()
        yi = 0.98
        for key in ("T", "P", "x", "v", "rho", "h", "s", "u"):
            self.notes.append(self.plot.ax.annotate(
                "%s: %s" % (key, state.__getattribute__(key).str), (0.01, yi),
                xycoords='axes fraction', size="small", va="center"))
            yi -= 0.025
        self.plot.draw()

    def clearPointData(self):
        while self.notes:
            anotation = self.notes.pop()
            anotation.remove()
        self.plot.draw()

    def writeToStream(self, stream):
        stream.writeString(self.filename)
        stream.writeString(self.windowTitle())
        stream.writeString(self.x)
        stream.writeString(self.y)
        stream.writeString(self.z)

        # TODO: Add support for save font properties
        stream.writeQString(self.plot.ax.get_title())
        stream.writeString(self.plot.ax.title.get_color())
        stream.writeQString(QtCore.QString(self.plot.ax.get_xlabel()))
        stream.writeString(self.plot.ax.xaxis.get_label().get_color())
        stream.writeQString(QtCore.QString(self.plot.ax.get_ylabel()))
        stream.writeString(self.plot.ax.yaxis.get_label().get_color())
        if self.z:
            stream.writeQString(QtCore.QString(self.plot.ax.get_zlabel()))
            stream.writeString(self.plot.ax.zaxis.get_label().get_color())
        stream.writeBool(self.plot.ax._gridOn)
        stream.writeString(self.plot.ax.get_xscale())
        stream.writeString(self.plot.ax.get_yscale())
        xmin, xmax=self.plot.ax.get_xlim()
        stream.writeFloat(xmin)
        stream.writeFloat(xmax)
        ymin, ymax=self.plot.ax.get_ylim()
        stream.writeFloat(ymin)
        stream.writeFloat(ymax)
        if self.z:
            zmin, zmax=self.plot.ax.get_zlim()
            stream.writeFloat(zmin)
            stream.writeFloat(zmax)

        # Margins
        stream.writeFloat(self.plot.fig.subplotpars.left)
        stream.writeFloat(self.plot.fig.subplotpars.bottom)
        stream.writeFloat(self.plot.fig.subplotpars.right)
        stream.writeFloat(self.plot.fig.subplotpars.top)

        # Config
        stream.writeInt32(self.config["fluid"])
        stream.writeInt32(self.config["eq"])
        stream.writeInt32(self.config["visco"])
        stream.writeInt32(self.config["thermal"])

    @classmethod
    def readFromStream(cls, stream, parent):
        filename = stream.readString()
        title = stream.readString()
        x = stream.readString()
        y = stream.readString()
        z = stream.readString()
        if z:
            dim = 3
        else:
            dim = 2
        grafico = PlotMEoS(dim=dim, parent=parent, filename=filename)
        grafico.x = x
        grafico.y = y
        grafico.z = z
        grafico.setWindowTitle(title)

        xtxt = "%s, %s" %(x, meos.units[meos.keys.index(x)].text())
        ytxt = "%s, %s" %(y, meos.units[meos.keys.index(y)].text())
        grafico.plot.ax.set_xlabel(xtxt)
        grafico.plot.ax.set_ylabel(ytxt)
        if z:
            ztxt = "%s, %s" %(z, meos.units[meos.keys.index(z)].text())
            grafico.plot.ax.set_zlabel(ztxt)

        plotTitle = stream.readQString()
        if plotTitle:
            grafico.plot.ax.set_title(unicode(plotTitle))
        titleColor = stream.readString()
        grafico.plot.ax.title.set_color(titleColor)
        xlabel = stream.readQString()
        if xlabel:
            grafico.plot.ax.set_xlabel(unicode(xlabel))
        xlabelColor = stream.readString()
        grafico.plot.ax.xaxis.get_label().set_color(xlabelColor)
        ylabel = stream.readQString()
        if ylabel:
            grafico.plot.ax.set_ylabel(unicode(ylabel))
        ylabelColor = stream.readString()
        grafico.plot.ax.yaxis.get_label().set_color(ylabelColor)
        if z:
            zlabel = stream.readQString()
            if zlabel:
                grafico.plot.ax.set_zlabel(unicode(zlabel))
            zlabelColor = stream.readString()
            grafico.plot.ax.zaxis.get_label().set_color(zlabelColor)
            
        grid = stream.readBool()
        grafico.plot.ax._gridOn = grid
        grafico.plot.ax.grid(grid)
        xscale = stream.readString()
        yscale = stream.readString()

        xmin = stream.readFloat()
        xmax = stream.readFloat()
        grafico.plot.ax.set_xlim(xmin, xmax)
        ymin = stream.readFloat()
        ymax = stream.readFloat()
        grafico.plot.ax.set_ylim(ymin, ymax)
        if z:
            zmin = stream.readFloat()
            zmax = stream.readFloat()
            grafico.plot.ax.set_zlim(zmin, zmax)

        data = grafico._getData()
        if z:
            plot2D3D(grafico, data, parent.Preferences, x, y, z)
        else:
            plot2D3D(grafico, data, parent.Preferences, x, y)
            
        if xscale:
            grafico.plot.ax.set_xscale(xscale)
        if yscale:
            grafico.plot.ax.set_yscale(yscale)

        # Load margins
        left = stream.readFloat()
        bottom = stream.readFloat()
        right = stream.readFloat()
        top = stream.readFloat()
        grafico.plot.fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

        # Load config
        conf = {}
        conf["fluid"] = stream.readInt32()
        conf["eq"] = stream.readInt32()
        conf["visco"] = stream.readInt32()
        conf["thermal"] = stream.readInt32()
        grafico.config = conf

        return grafico

class Plot2D(QtGui.QDialog):
    """Dialog for select a special 2D plot"""
    def __init__(self, config=None, parent=None):
        super(Plot2D, self).__init__(parent)
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Setup 2D Plot"))
        layout = QtGui.QGridLayout(self)
        group_Ejex=QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Axis X"))
        layout.addWidget(group_Ejex,1,1)
        layout_GroupX=QtGui.QGridLayout(group_Ejex)
        self.ejeX = QtGui.QComboBox()
        layout_GroupX.addWidget(self.ejeX,1,1)
        self.Xscale=QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Logarithmic scale"))
        layout_GroupX.addWidget(self.Xscale,2,1)
        for prop in prop_pickle:
            self.ejeX.addItem(meos.properties[prop])

        group_Ejey=QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Axis Y"))
        layout.addWidget(group_Ejey,2,1)
        layout_GroupY=QtGui.QGridLayout(group_Ejey)
        self.ejeY = QtGui.QComboBox()
        layout_GroupY.addWidget(self.ejeY,1,1)
        self.Yscale=QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Logarithmic scale"))
        layout_GroupY.addWidget(self.Yscale,2,1)

        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok|QtGui.QDialogButtonBox.Cancel)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox,10,1)

        self.ejeXChanged(0)
        self.ejeX.currentIndexChanged.connect(self.ejeXChanged)

    def ejeXChanged(self, int):
        """Rellena las variables disponibles para el ejeY en el gráfico 2D, todos menos el que este activo en el ejeX"""
        self.ejeY.clear()
        prop2 = prop_pickle[:]
        del prop2[int]
        for prop in prop2:
            self.ejeY.addItem(meos.properties[prop])


class Plot3D(QtGui.QDialog):
    """Widget de configuracion inicial de graficos 3D"""

    def __init__(self, parent=None):
        super(Plot3D, self).__init__(parent)
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Setup 3D Plot"))
        layout = QtGui.QGridLayout(self)

        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Axis X")),1,1)
        self.ejeX = QtGui.QComboBox()
        for prop in prop_pickle:
            self.ejeX.addItem(meos.properties[prop])

        layout.addWidget(self.ejeX,1,2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Axis Y")),2,1)
        self.ejeY = QtGui.QComboBox()
        layout.addWidget(self.ejeY,2,2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Axis Z")),3,1)
        self.ejeZ = QtGui.QComboBox()
        layout.addWidget(self.ejeZ,3,2)

        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok|QtGui.QDialogButtonBox.Cancel)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox,4,1,1,2)

        self.ejeX.currentIndexChanged.connect(self.ejeXChanged)
        self.ejeY.currentIndexChanged.connect(self.ejeYChanged)
        self.ejeXChanged(0)

    def ejeXChanged(self, int):
        """Rellena las variables disponibles para el ejeY en el gráfico 3D, todos menos el que este activo en el ejeX"""
        self.ejeY.clear()
        prop2 = prop_pickle[:]
        del prop2[int]
        for prop in prop2:
            self.ejeY.addItem(meos.properties[prop])

    def ejeYChanged(self, int):
        """Rellena las variables disponibles para el ejeZ en el gráfico 3D, todos menos los activos en el ejeX y ejeY"""
        self.ejeZ.clear()
        prop2 = prop_pickle[:]
        intX = self.ejeX.currentIndex()
        del prop2[intX]
        del prop2[int]
        for prop in prop2:
            self.ejeZ.addItem(meos.properties[prop])


class EditPlot(QtGui.QWidget):
    """Dialog to edit plot"""
    def __init__(self, plotMEoS, mainwindow, parent=None):
        super(EditPlot, self).__init__(parent)
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Edit Plot"))
        layout = QtGui.QGridLayout(self)
        self.plotMEoS=plotMEoS
        self.fig=plotMEoS.plot
        self.mainwindow=mainwindow

        self.lista=QtGui.QListWidget()
        layout.addWidget(self.lista,0,1,1,3)

        lytTitle = QtGui.QHBoxLayout()
        label=QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Label"))
        lytTitle.addWidget(label)
        self.label=QtGui.QLineEdit()
        lytTitle.addWidget(self.label)
        layout.addLayout(lytTitle,1,1,1,3)

        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Line Width")),2,1)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Line Style")),2,2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Color")),2,3)
        self.Grosor = QtGui.QDoubleSpinBox()
        self.Grosor.setAlignment(QtCore.Qt.AlignRight)
        self.Grosor.setRange(0.1, 5)
        self.Grosor.setDecimals(1)
        self.Grosor.setSingleStep(0.1)
        layout.addWidget(self.Grosor,3,1)
        self.Linea = LineStyleCombo()
        layout.addWidget(self.Linea,3,2)
        self.ColorButton = ColorSelector()
        layout.addWidget(self.ColorButton,3,3)

        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Marker")),4,1)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Marker Size")),4,2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Marker Color")),4,3)
        self.Marca = MarkerCombo()
        layout.addWidget(self.Marca,5,1)
        self.markerSize = QtGui.QDoubleSpinBox()
        self.markerSize.setAlignment(QtCore.Qt.AlignRight)
        self.markerSize.setDecimals(1)
        self.markerSize.setSingleStep(0.1)
        layout.addWidget(self.markerSize,5,2)
        self.markerfacecolor = ColorSelector()
        layout.addWidget(self.markerfacecolor,5,3)

        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Marker edge")),7,1)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Width")),6,2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Color")),6,3)
        self.markerEdgeSize = QtGui.QDoubleSpinBox()
        self.markerEdgeSize.setAlignment(QtCore.Qt.AlignRight)
        self.markerEdgeSize.setDecimals(1)
        self.markerEdgeSize.setSingleStep(0.1)
        layout.addWidget(self.markerEdgeSize,7,2)
        self.markeredgecolor = ColorSelector()
        layout.addWidget(self.markeredgecolor,7,3)

        self.visible=QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Visible"))
        layout.addWidget(self.visible,8,1,1,3)
        self.antialiases=QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Antialiases"))
        layout.addWidget(self.antialiases,9,1,1,3)

        layoutButton=QtGui.QHBoxLayout()
        layout.addLayout(layoutButton,10,1,1,3)
        self.botonAdd=QtGui.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/add.png")), "")
        self.botonAdd.clicked.connect(self.add)
        layoutButton.addWidget(self.botonAdd)
        self.botonRemove=QtGui.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/remove.png")), "")
        self.botonRemove.clicked.connect(self.remove)
        layoutButton.addWidget(self.botonRemove)
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Close)
        self.buttonBox.rejected.connect(self.close)
        layoutButton.addWidget(self.buttonBox)

        for linea in self.fig.ax.lines[2:]:
            self.lista.addItem(linea._label)

        self.lista.currentRowChanged.connect(self.update)
        self.label.textChanged.connect(partial(self.changeValue, "label"))
        self.Grosor.valueChanged.connect(partial(self.changeValue, "lw"))
        self.Linea.valueChanged.connect(partial(self.changeValue, "ls"))
        self.Linea.currentIndexChanged.connect(self.ColorButton.setEnabled)
        self.ColorButton.valueChanged.connect(partial(self.changeValue, "color"))
        self.Marca.valueChanged.connect(partial(self.changeValue, "marker"))
        self.Marca.currentIndexChanged.connect(self.markerSize.setEnabled)
        self.Marca.currentIndexChanged.connect(self.markerfacecolor.setEnabled)
        self.Marca.currentIndexChanged.connect(self.markerEdgeSize.setEnabled)
        self.Marca.currentIndexChanged.connect(self.markeredgecolor.setEnabled)
        self.markerSize.valueChanged.connect(partial(self.changeValue, "ms"))
        self.markerfacecolor.valueChanged.connect(partial(self.changeValue, "mfc"))
        self.markerEdgeSize.valueChanged.connect(partial(self.changeValue, "mew"))
        self.markeredgecolor.valueChanged.connect(partial(self.changeValue, "mec"))
        self.visible.toggled.connect(partial(self.changeValue, "visible"))
        self.antialiases.toggled.connect(partial(self.changeValue, "antialiases"))
        self.lista.setCurrentRow(0)

    def update(self, i):
        """Fill format widget with value of selected line"""
        linea=self.fig.ax.lines[i+2]
        self.label.setText(linea.get_label())
        self.Grosor.setValue(linea.get_lw())
        self.Linea.setCurrentValue(linea.get_ls())
        self.ColorButton.setColor(linea.get_color())
        self.Marca.setCurrentValue(linea.get_marker())
        self.Marca.currentIndexChanged.emit(self.Marca.currentIndex())
        self.markerSize.setValue(linea.get_ms())
        self.markerfacecolor.setColor(linea.get_mfc())
        self.markerEdgeSize.setValue(linea.get_mew())
        self.markeredgecolor.setColor(linea.get_mec())
        self.visible.setChecked(linea.get_visible())
        self.antialiases.setChecked(linea.get_antialiased())

    def changeValue(self, key, value):
        """Actualiza datos del grafico, cambios hechos al vuelo en el grafico"""
        linea=self.fig.ax.lines[self.lista.currentRow()+2]
        func={"label": linea.set_label,
                    "lw": linea.set_lw,
                    "ls": linea.set_ls,
                    "marker": linea.set_marker,
                    "color": linea.set_color,
                    "ms": linea.set_ms,
                    "mfc": linea.set_mfc,
                    "mew": linea.set_mew,
                    "mec": linea.set_mec,
                    "visible": linea.set_visible,
                    "antialiases": linea.set_antialiased}
        if key in ("ls", "marker", "color", "mfc", "mec"):
            value=str(value)
        func[key](value)
        if key=="label":
            self.lista.currentItem().setText(value)
        else:
            self.fig.draw()

    def add(self):
        """Add a isoline to plot"""
        dialog=AddLine()
        if dialog.exec_():
            points = get_points(self.mainwindow.Preferences)
            self.mainwindow.progressBar.setVisible(True)
            fluid=mEoS.__all__[self.mainwindow.currentConfig.getint("MEoS", "fluid")]
            prop=dialog.tipo.currentIndex()
            value=dialog.input[prop].value

            eq = fluid.eq[self.mainwindow.currentConfig.getint("MEoS", "eq")]
            T = list(concatenate([linspace(eq["Tmin"], 0.9*fluid.Tc, points),
                                  linspace(0.9*fluid.Tc, 0.99*fluid.Tc, points),
                                  linspace(0.99*fluid.Tc, fluid.Tc, points),
                                  linspace(fluid.Tc, 1.01*fluid.Tc, points),
                                  linspace(1.01*fluid.Tc, 1.1*fluid.Tc, points),
                                  linspace(1.1*fluid.Tc, eq["Tmax"], points)]))
            Pmin = eq["Pmin"]*1000
            Pmax = eq["Pmax"]*1000
            P = list(concatenate([logspace(log10(Pmin), log10(0.9*fluid.Pc), points),
                                  linspace(0.9*fluid.Pc, 0.99*fluid.Pc, points),
                                  linspace(0.99*fluid.Pc, fluid.Pc, points),
                                  linspace(fluid.Pc, 1.01*fluid.Pc, points),
                                  linspace(1.01*fluid.Pc, 1.1*fluid.Pc, points),
                                  logspace(log10(1.1*fluid.Pc), log10(Pmax), points)]))
            for i in range(5, 0, -1):
                del T[points*i]
                del P[points*i]

            if prop == 0:
                # Calcualte isotherm line
                self.mainwindow.statusbar.showMessage(QtGui.QApplication.translate(
                    "pychemqt", "Adding isotherm line..."))
                fluidos = calcIsoline(fluid, self.mainwindow.currentConfig,
                                      "P", "T", P, value, 0, 0, 100,
                                      1, self.mainwindow.progressBar)
                var = "T"
                name = "Isotherm"
                unidad = unidades.Temperature
            elif prop == 1:
                # Calculate isobar line
                self.mainwindow.statusbar.showMessage(QtGui.QApplication.translate(
                    "pychemqt", "Adding isobar line..."))
                fluidos = calcIsoline(fluid, self.mainwindow.currentConfig,
                                      "T", "P", T, value, 0, 0, 100,
                                      1, self.mainwindow.progressBar)
                var = "P"
                name = "Isobar"
                unidad = unidades.Pressure
            elif prop == 2:
                # Calculate isoenthalpic line
                self.mainwindow.statusbar.showMessage(QtGui.QApplication.translate(
                    "pychemqt", "Adding isoenthalpic line..."))
                fluidos = calcIsoline(fluid, self.mainwindow.currentConfig,
                                      "P", "h", P, value, 0, 0, 100,
                                      1, self.mainwindow.progressBar)
                var = "h"
                name = "Isoenthalpic"
                unidad = unidades.Enthalpy
            elif prop == 3:
                # Calculate isoentropic line
                self.mainwindow.statusbar.showMessage(QtGui.QApplication.translate(
                    "pychemqt", "Adding isoentropic line..."))
                fluidos = calcIsoline(fluid, self.mainwindow.currentConfig,
                                      "T", "s", T, value, 0, 0, 100,
                                      1, self.mainwindow.progressBar)
                var = "s"
                name = "Isoentropic"
                unidad = unidades.SpecificHeat
            elif prop == 4:
                # Calculate isochor line
                self.mainwindow.statusbar.showMessage(QtGui.QApplication.translate(
                    "pychemqt", "Adding isochor line..."))
                fluidos = calcIsoline(fluid, self.mainwindow.currentConfig,
                                      "T", "v", T, value, 0, 0, 100,
                                      1, self.mainwindow.progressBar)
                var = "v"
                name = "Isochor"
                unidad = unidades.SpecificVolume
            elif prop == 5:
                # Calculate isodensity line
                self.mainwindow.statusbar.showMessage(QtGui.QApplication.translate(
                    "pychemqt", "Adding isodensity line..."))
                fluidos = calcIsoline(fluid, self.mainwindow.currentConfig,
                                      "T", "rho", T, value, 0, 0, 100,
                                      1, self.mainwindow.progressBar)
                var = "rho"
                name = "Isochor"
                unidad = unidades.Density
            elif prop == 6:
                # Calculate isoquality line
                self.mainwindow.statusbar.showMessage(QtGui.QApplication.translate(
                    "pychemqt", "Adding isoquality line..."))
                T = T[:3*points-2]
                fluidos = calcIsoline(fluid, self.mainwindow.currentConfig,
                                      "T", "x", T, value, 0, 0, 100,
                                      1, self.mainwindow.progressBar)
                var = "x"
                name = "Isoquality"
                unidad = unidades.Dimensionless

            line = {value: {}}
            for x in prop_pickle:
                dat_propiedad=[]
                for fluido in fluidos:
                    num = fluido.__getattribute__(x)
                    if num is not None:
                        dat_propiedad.append(num._data)
                    else:
                        dat_propiedad.append(None)
                line[value][x]=dat_propiedad

            format = getLineFormat(self.mainwindow.Preferences, name)
            transform = _getunitTransform((self.plotMEoS.x, self.plotMEoS.y))
            if self.plotMEoS.dim==3:
                plotIsoline(line, (self.plotMEoS.x, self.plotMEoS.y, self.plotMEoS.z), var, unidad, self.plotMEoS, transform, **format)
            else:
                plotIsoline(line, (self.plotMEoS.x, self.plotMEoS.y), var, unidad, self.plotMEoS, transform, **format)

            self.plotMEoS.plot.draw()
            self.mainwindow.progressBar.setVisible(False)
            self.lista.addItem(self.fig.ax.lines[-1].get_label())
            self.lista.setCurrentRow(self.lista.count()-1)

            # Save new line to file
            data = self.plotMEoS._getData()
            if var not in data:
                data[var] = {}
            data[var][value] = line[value]
            self.plotMEoS._saveData(data)

    def remove(self):
        """Remove a line from plot"""
        self.mainwindow.statusbar.showMessage(QtGui.QApplication.translate(
            "pychemqt", "Deleting line..."))
        QtGui.QApplication.processEvents()

        # Remove data from file
        data = self.plotMEoS._getData()
        txt = unicode(self.lista.currentItem().text()).split()
        var = txt[0]
        units = {"T": unidades.Temperature,
                 "P": unidades.Pressure,
                 "v": unidades.SpecificVolume,
                 "rho": unidades.Density,
                 "h": unidades.Enthalpy,
                 "s": unidades.SpecificHeat,
                 "x": unidades.Dimensionless}
        if var in units:
            unit = units[var]
            for key in data[var]:
                str = unit(key).str
                if str[1:] == " ".join(txt[2:]):
                    del data[var][key]
                    self.plotMEoS._saveData(data)
                    break

        # Remove line to plot and update list element
        index = self.lista.currentRow()
        del self.fig.ax.lines[index+2]
        if index == 0:
            self.lista.setCurrentRow(1)
        else:
            self.lista.setCurrentRow(index-1)
        self.lista.takeItem(index)
        self.fig.draw()
        self.mainwindow.statusbar.clearMessage()


class AddLine(QtGui.QDialog):
    """Dialog to add new isoline to plot"""
    lineas = [(QtGui.QApplication.translate("pychemqt", "Isotherm"), unidades.Temperature, None),
              (QtGui.QApplication.translate("pychemqt", "Isobar"), unidades.Pressure, None),
              (QtGui.QApplication.translate("pychemqt", "Isoenthalpic"), unidades.Enthalpy, None),
              (QtGui.QApplication.translate("pychemqt", "Isoentropic"), unidades.SpecificHeat, "SpecificEntropy"),
              (QtGui.QApplication.translate("pychemqt", "Isochor"), unidades.SpecificVolume, None),
              (QtGui.QApplication.translate("pychemqt", "Isodensity"), unidades.Density, None),
              (QtGui.QApplication.translate("pychemqt", "Isoquality"), float, None)]

    def __init__(self, parent=None):
        super(AddLine, self).__init__(parent)
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Add Line to Plot"))
        layout = QtGui.QGridLayout(self)

        self.tipo=QtGui.QComboBox()
        layout.addWidget(self.tipo,1,1,1,2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Value")),2,1)

        self.input=[]
        for title, unidad, magnitud in self.lineas:
            self.input.append(Entrada_con_unidades(unidad, magnitud))
            layout.addWidget(self.input[-1],2,2)
            self.tipo.addItem(title)

        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok|QtGui.QDialogButtonBox.Cancel)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox,10,1,1,2)

        self.isolineaChanged(0)
        self.tipo.currentIndexChanged.connect(self.isolineaChanged)

    def isolineaChanged(self, int):
        """Deja visible solo entrada seleccionada en la lista"""
        for i in self.input:
            i.setVisible(False)
        self.input[int].setVisible(True)


class EditAxis(QtGui.QDialog):
    """Dialog to configure axes plot properties"""
    def __init__(self, fig=None, parent=None):
        super(EditAxis, self).__init__(parent)
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Edit Axis"))
        layout = QtGui.QGridLayout(self)
        self.fig=fig
        
        lytTitle = QtGui.QHBoxLayout()
        label=QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Title"))
        label.setSizePolicy(QtGui.QSizePolicy.Maximum,QtGui.QSizePolicy.Maximum)
        lytTitle.addWidget(label)
        self.title=InputFond()
        lytTitle.addWidget(self.title)
        layout.addLayout(lytTitle,1,1,1,self.fig.dim)

        self.axisX = AxisWidget("x", self)
        layout.addWidget(self.axisX,2,1)
        self.axisY = AxisWidget("y", self)
        layout.addWidget(self.axisY,2,2)

        if self.fig.dim == 3:
            self.axisZ = AxisWidget("z", self)
            layout.addWidget(self.axisZ,2,3)
            self.axisX.scale.setEnabled(False)
            self.axisY.scale.setEnabled(False)
            self.axisZ.scale.setEnabled(False)
            
        self.gridCheckbox=QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Show Grid"))
        layout.addWidget(self.gridCheckbox,3,1,1,self.fig.dim)
        layout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),5,1,1,self.fig.dim)
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Close)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox,10,1,1,self.fig.dim)

        if fig:
            self.populate()

        self.title.textChanged.connect(partial(self.update, "title"))
        self.title.colorChanged.connect(partial(self.update, "titlecolor"))
        self.title.fontChanged.connect(partial(self.update, "titlefont"))
        self.axisX.label.textChanged.connect(partial(self.update, "xlabel"))
        self.axisX.label.colorChanged.connect(partial(self.update, "xlabelcolor"))
        self.axisX.label.fontChanged.connect(partial(self.update, "xlabelfont"))
        self.axisY.label.textChanged.connect(partial(self.update, "ylabel"))
        self.axisY.label.colorChanged.connect(partial(self.update, "ylabelcolor"))
        self.axisY.label.fontChanged.connect(partial(self.update, "ylabelfont"))
        self.gridCheckbox.toggled.connect(partial(self.update, "grid"))
        self.axisX.scale.toggled.connect(partial(self.update, "xscale"))
        self.axisY.scale.toggled.connect(partial(self.update, "yscale"))
        self.axisX.min.valueChanged.connect(partial(self.update, "xmin"))
        self.axisY.min.valueChanged.connect(partial(self.update, "ymin"))
        self.axisX.max.valueChanged.connect(partial(self.update, "xmax"))
        self.axisY.max.valueChanged.connect(partial(self.update, "ymax"))
        if self.fig.dim == 3:
            self.axisZ.label.textChanged.connect(partial(self.update, "zlabel"))
            self.axisZ.label.colorChanged.connect(partial(self.update, "zlabelcolor"))
            self.axisZ.label.fontChanged.connect(partial(self.update, "zlabelfont"))
            self.axisZ.min.valueChanged.connect(partial(self.update, "zmin"))
            self.axisZ.max.valueChanged.connect(partial(self.update, "zmax"))


    def populate(self):
        """Rellena widgets con los valores del gráfico"""
        self.title.setText(self.fig.ax.get_title())
        self.title.setColor(QtGui.QColor(self.fig.ax.title.get_color()))
        self.axisX.label.setText(self.fig.ax.get_xlabel())
        self.axisX.label.setColor(QtGui.QColor(self.fig.ax.xaxis.get_label().get_color()))
        self.axisY.label.setText(self.fig.ax.get_ylabel())
        self.axisY.label.setColor(QtGui.QColor(self.fig.ax.yaxis.get_label().get_color()))
        self.gridCheckbox.setChecked(self.fig.ax._gridOn)
        self.axisX.scale.setChecked(self.fig.ax.get_xscale()=="log")
        self.axisY.scale.setChecked(self.fig.ax.get_yscale()=="log")
        xmin, xmax=self.fig.ax.get_xlim()
        self.axisX.min.setValue(xmin)
        self.axisX.max.setValue(xmax)
        ymin, ymax=self.fig.ax.get_ylim()
        self.axisY.min.setValue(ymin)
        self.axisY.max.setValue(ymax)
        if self.fig.dim == 3:
            self.axisZ.label.setText(self.fig.ax.get_zlabel())
            self.axisZ.label.setColor(QtGui.QColor(self.fig.ax.zaxis.get_label().get_color()))
            zmin, zmax=self.fig.ax.get_zlim()
            self.axisZ.min.setValue(zmin)
            self.axisZ.max.setValue(zmax)


    def update(self, key, value):
        """Actualiza datos del grafico, cambios hechos al vuelo en el grafico"""
        func={"xlabel": self.fig.ax.set_xlabel,
                "xlabelcolor": self.fig.ax.xaxis.get_label().set_color,
                "xlabelfont": self.fig.ax.xaxis.get_label().set_fontproperties,
                "ylabel": self.fig.ax.set_ylabel,
                "ylabelcolor": self.fig.ax.yaxis.get_label().set_color,
                "ylabelfont": self.fig.ax.yaxis.get_label().set_fontproperties,
                "title": self.fig.ax.set_title,
                "titlecolor": self.fig.ax.title.set_color,
                "titlefont": self.fig.ax.title.set_fontproperties,
                "xscale": self.fig.ax.set_xscale,
                "yscale": self.fig.ax.set_yscale,
                "grid": self.fig.ax.grid}
        
        if self.fig.dim == 3:
            func["zlabel"] = self.fig.ax.set_zlabel
            func["zlabelcolor"] = self.fig.ax.zaxis.get_label().set_color
            func["zlabelfont"] = self.fig.ax.zaxis.get_label().set_fontproperties
            
        
        if key in ("xscale", "yscale"):
            if value:
                value="log"
            else:
                value="linear"
        if key == "grid":
            self.fig.ax._gridOn = value
        if key in ("titlecolor", "xlabelcolor", "ylabelcolor"):
            value=str(value)
        if key in ("titlefont", "xlabelfont", "ylabelfont"):
            value=self.convertFont(value)

        if key in ("xmin", "xmax"):
            xmin=self.axisX.min.value
            xmax=self.axisX.max.value
            self.fig.ax.set_xlim(xmin, xmax)
        elif key in ("ymin", "ymax"):
            ymin=self.axisY.min.value
            ymax=self.axisY.max.value
            self.fig.ax.set_ylim(ymin, ymax)
        elif key in ("zmin", "zmax"):
            ymin=self.axisZ.min.value
            ymax=self.axisZ.max.value
            self.fig.ax.set_zlim(ymin, ymax)
        else:
            func[key](value)
        self.fig.draw()

    def convertFont(self, font):
        """Convierte la QFont devuelta por QFontDialog en una FontProperties usada por matplotlib"""
        family=str(font.family())
        if str(font.style()) in ("normal", "italic", "oblique"):
            style=str(font.style())
        else:
            style=None
        font=FontProperties(family, style, None, font.stretch(), font.weight(), font.pointSize())
        return font


class AxisWidget(QtGui.QGroupBox):
    """Dialog to configure axes plot properties"""
    def __init__(self, name, parent=None):
        title = name+" "+QtGui.QApplication.translate("pychemqt", "Axis")
        super(AxisWidget, self).__init__(title, parent)
        lyt=QtGui.QGridLayout(self)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Label")),1,1)
        self.label=InputFond()
        lyt.addWidget(self.label,1,2)
        self.scale=QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Logarithmic scale"))
        lyt.addWidget(self.scale,2,1,1,2)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "from")),3,1)
        self.min=Entrada_con_unidades(float, min=float("-inf"))
        lyt.addWidget(self.min,3,2)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "to")),4,1)
        self.max=Entrada_con_unidades(float, min=float("-inf"))
        lyt.addWidget(self.max,4,2)


def calcIsoline(fluid, config, var, fix, valuevar, valuefix, ini, step, end, total, bar):
    fluidos = []
    fase = None
    rhoo = 0
    To = 0
    for Ti in valuevar:
        kwargs = {var: Ti, fix: valuefix, "rho0": rhoo, "T0": To}
        fluido = calcPoint(fluid, config, **kwargs)
        if fluido and fluido.status and (fluido.rho != rhoo or fluido.T != To):
            if var not in ("T", "P") or fix not in ("T", "P"):
                rhoo = fluido.rho
                To = fluido.T

            fluidos.append(fluido)
#            if var in ("T", "P") and fix in ("T", "P"):
#                if fase is None:
#                    fase = fluido.x
#                if fase != fluido.x and fase <= 0:
#                    if fluido.P < fluid.Pc and fluido.T < fluid.Tc:
#                        fluido_x0 = calcPoint(fluid, config, **{fix: valuefix, "x": 0.})
#                        fluidos.insert(-1, fluido_x0)
#                elif fase != fluido.x and fase >= 1:
#                    if fluido.P < fluid.Pc and fluido.T < fluid.Tc:
#                        fluido_x1 = calcPoint(fluid, config, **{fix: valuefix, "x": 1.})
#                        fluidos.insert(-1, fluido_x1)
#                if fase != fluido.x and fluido.x >= 1:
#                    if fluido.P < fluid.Pc and fluido.T < fluid.Tc:
#                        fluido_x1 = calcPoint(fluid, config, **{fix: valuefix, "x": 1.})
#                        fluidos.insert(-1, fluido_x1)
##                        rhoo = fluido_x1.rho
##                        To = fluido_x1.T
#                elif fase != fluido.x and fluido.x <= 0:
#                    if fluido.P < fluid.Pc and fluido.T < fluid.Tc:
#                        fluido_x0 = calcPoint(fluid, config, **{fix: valuefix, "x": 0.})
#                        fluidos.insert(-1, fluido_x0)
##                        rhoo = fluido_x0.rho
##                        To = fluido_x0.T
#                fase = fluido.x

        bar.setValue(ini+end*step/total+end/total*len(fluidos)/len(valuevar))
        QtGui.QApplication.processEvents()
    return fluidos


def get_points(Preferences):
    """Get point number to plot lines from Preferences"""
    definition = Preferences.getint("MEOS", "definition")
    if definition == 1:
        points = 10
    elif definition == 2:
        points = 25
    elif definition == 3:
        points = 50
    elif definition == 4:
        points = 100
    else:
        points = 5
    return points


def getLineFormat(Preferences, name):
    """get matplotlib line format from preferences
        Preferences: configparser instance with pycheqmt preferences
        name: name of isoline"""
    format={}
    format["ls"]=Preferences.get("MEOS", name+"lineStyle")
    format["lw"]=Preferences.getfloat("MEOS", name+"lineWidth")
    format["color"]=Preferences.get("MEOS", name+"Color")
    format["marker"]=Preferences.get("MEOS", name+"marker")
    format["ms"]=3

    # Anotation
    if name != "saturation":
        format["annotate"]=Preferences.getboolean("MEOS", name+"label")
        format["pos"]=Preferences.getint("MEOS", name+"position")
        format["unit"]=Preferences.getboolean("MEOS", name+"units")
        format["variable"]=Preferences.getboolean("MEOS", name+"variable")

    return format


def plotIsoline(data, axis, title, unidad, grafico, transform, **format):
    x, y, z = axis
    fx, fy, fz = transform
    xscale = grafico.plot.ax.get_xscale()
    yscale = grafico.plot.ax.get_yscale()
    annotate = format.pop("annotate")
    pos = format.pop("pos")
    unit = format.pop("unit")
    variable = format.pop("variable")
    for key in sorted(data.keys()):
        xi = map(fx, data[key][x])
        yi = map(fy, data[key][y])
        label = "%s =%s" % (title, unidad(key).str)
        if z:
            zi = map(fz, data[key][z])
            grafico.plot.ax.plot(xi, yi, zi, label=label, **format)
        else:
            grafico.plot.ax.plot(xi, yi, label=label, **format)

        # Add annotate for isolines
        if annotate and not z:
            if variable and unit:
                txt = label
            elif variable:
                txt =  "%s =%s" % (title, unidad(key).config())
            elif unit:
                txt = unidad(key).str
            else:
                txt = unidad(key).config()

            xmin, xmax = grafico.plot.ax.get_xlim()
            ymin, ymax = grafico.plot.ax.get_ylim()

            i = int(len(xi)*pos/100)
            if pos > 50:
                j = i-10
            else:
                j = i+10
            if xscale == "log":
                fx=(log(xi[i])-log(xi[j]))/(log(xmax)-log(xmin))
            else:
                fx=(xi[i]-xi[j])/(xmax-xmin)
            if yscale == "log":
                fy=(log(yi[i])-log(yi[j]))/(log(ymax)-log(ymin))
            else:
                fy=(yi[i]-yi[j])/(ymax-ymin)

            rot = atan(fy/fx)*360/2/pi
            grafico.plot.ax.annotate(txt, (xi[i], yi[i]),
                    rotation=rot, size="small", ha="center", va="center")


def plot2D3D(grafico, data, Preferences, x, y, z=None):
    """Plot procedure
    Parameters:
        grafico: plot
        data: data to plot
        Preferences: ConfigParser instance from mainwindow preferencesChanged
        x: Key for x axis
        y: Key for y axis
        z: Key for z axis Optional for 3D plot"""

    functionx = _getunitTransform(x)
    functiony = _getunitTransform(y)
    functionz = _getunitTransform(z)
    transform = (functionx, functiony, functionz)

    # Plot saturation lines
    format = getLineFormat(Preferences, "saturation")
    if x == "P" and y == "T":
        satLines = [(QtGui.QApplication.translate("pychemqt", "Saturation Line"), 0)]
    else:
        satLines = [(QtGui.QApplication.translate("pychemqt", "Liquid Saturation Line"), 0),
                    (QtGui.QApplication.translate("pychemqt", "Vapor Saturation Line"), 1)]
    for label, fase in satLines:
        xsat = map(functionx, data["saturation_%i" %fase][x])
        ysat = map(functiony, data["saturation_%i" %fase][y])
        if z:
            zsat = map(functionz, data["saturation_%i" %fase][z])
            grafico.plot.ax.plot(xsat, ysat, zsat, label=label, **format)
        else:
            grafico.plot.ax.plot(xsat, ysat, label=label, **format)

    # Plot melting and sublimation lines
    if "melting" in data:
        label = QtGui.QApplication.translate("pychemqt", "Melting Line")
        xmel = map(functionx, data["melting"][x])
        ymel = map(functiony, data["melting"][y])
        if z:
            zmel = map(functionz, data["melting"][z])
            grafico.plot.ax.plot(xmel, ymel, zmel, label=label, **format)
        else:
            grafico.plot.ax.plot(xmel, ymel, label=label, **format)
    if "sublimation" in data:
        xsub = map(functionx, data["sublimation"][x])
        ysub = map(functiony, data["sublimation"][y])
        label = QtGui.QApplication.translate("pychemqt", "Sublimation Line")
        if z:
            zmel = map(functionz, data["melting"][z])
            grafico.plot.ax.plot(xmel, ymel, zmel, label=label, **format)
        else:
            grafico.plot.ax.plot(xsub, ysub, label=label, **format)

    # Plot quality isolines
    if x not in ["P", "T"] or y not in ["P", "T"] or z:
        format = getLineFormat(Preferences, "Isoquality")
        plotIsoline(data["x"], (x, y, z), "x", unidades.Dimensionless, grafico, transform, **format)

    # Plot isotherm lines
    if x != "T" and y != "T" or z:
        format = getLineFormat(Preferences, "Isotherm")
        plotIsoline(data["T"], (x, y, z), "T", unidades.Temperature, grafico, transform, **format)

    # Plot isobar lines
    if x != "P" and y != "P" or z:
        format = getLineFormat(Preferences, "Isobar")
        plotIsoline(data["P"], (x, y, z), "P", unidades.Pressure, grafico, transform, **format)

    # Plot isochor lines
    if x not in ["rho", "v"] and y not in ["rho", "v"] or z:
        format = getLineFormat(Preferences, "Isochor")
        plotIsoline(data["v"], (x, y, z), "v", unidades.SpecificVolume, grafico, transform, **format)
        # Plot isodeensity lines
        if "rho" in data:
            plotIsoline(data["rho"], (x, y, z), "rho", unidades.Density, grafico, transform, **format)

    # Plot isoenthalpic lines
    if x != "h" and y != "h" or z:
        format = getLineFormat(Preferences, "Isoenthalpic")
        plotIsoline(data["h"], (x, y, z), "h", unidades.Enthalpy, grafico, transform, **format)

    # Plot isoentropic lines
    if x != "s" and y != "s" or z:
        format = getLineFormat(Preferences, "Isoentropic")
        plotIsoline(data["s"], (x, y, z), "s", unidades.SpecificHeat, grafico, transform, **format)


def _getunitTransform(eje):
    """Return the axis transform function to map all plot data to configurated unit
        Parameters:
            seq: list with axis property keys
    """
    if not eje:
        return None
    elif eje == "T":
        index = config.getMainWindowConfig().getint("Units", "Temperature")
        return [None, unidades.K2C, unidades.K2R, unidades.K2F, unidades.K2Re][index]
    else:
        unit = meos.units[meos.keys.index(eje)]
        factor = unit(1.).config()
        return lambda val: val*factor if val is not None else nan
    

def calcPoint(fluid, config, **kwargs):
    if isinstance(config, dict):
        option = config
    else:
        option = {}
        option["eq"] = config.getint("MEoS", "eq")
        option["visco"] = config.getint("MEoS", "visco")
        option["thermal"] = config.getint("MEoS", "thermal")
    kwargs.update(option)
    Tmin = fluid.eq[option["eq"]]["Tmin"]
    Tmax = fluid.eq[option["eq"]]["Tmax"]
    Pmin = fluid.eq[option["eq"]]["Pmin"]*1000
    Pmax = fluid.eq[option["eq"]]["Pmax"]*1000
    if "T" in kwargs:
        if kwargs["T"] < Tmin or kwargs["T"] > Tmax:
            return None
    if "P" in kwargs:
        if kwargs["P"] < Pmin-1 or kwargs["P"] > Pmax+1:
            return None
    fluido = fluid(**kwargs)

    if fluido.status not in [1, 3]:
        return None
    if fluido._melting and fluido._melting["Tmin"] < fluido.T < fluido._melting["Tmax"]:
        Pmel = fluido._Melting_Pressure(fluido.T)
        Pmax = min(Pmax, Pmel)

    if fluido.P < Pmin-1 or fluido.P > Pmax+1 or fluido.T < Tmin or fluido.T > Tmax:
        return None
    return fluido




if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)

    SteamTables = EditAxis()
#    SteamTables=AddLine(None)
#    SteamTables=transportDialog(mEoS.__all__[2])

    SteamTables.show()
    sys.exit(app.exec_())

