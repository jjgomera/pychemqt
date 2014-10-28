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
#   Dialogs for input data:
#   - Ui_Saturation: Dialog to define input for a two-phase table calculation
#   - Ui_Isoproperty: Dialog to define input for isoproperty table calculations
###############################################################################

import inspect, csv, os
from functools import partial
from string import maketrans

from PyQt4 import QtCore, QtGui
from numpy import arange, append, concatenate, meshgrid, zeros, linspace, logspace, min, max
from matplotlib.lines import Line2D
from matplotlib.font_manager import FontProperties

from lib import meos, mEoS, unidades, plot, iapws
from lib.utilities import format2txt, representacion
from UI.widgets import Entrada_con_unidades, ClickableLabel, Tabla, createAction, LineStyleCombo, MarkerCombo, ColorSelector, InputFond
from UI.delegate import CheckEditor
from tools.codeEditor import SimplePythonEditor
from tools.UI_Preferences import NumericFactor

if os.environ["freesteam"]:
    import freesteam


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
        SaturationSpecifyAction = createAction(QtGui.QApplication.translate("pychemqt", "Saturation Specified point"), slot=self.addTableSpecified, parent=self)
        menuCalculate.addAction(SaturationSpecifyAction)

        menuPlot=QtGui.QMenu(QtGui.QApplication.translate("pychemqt", "Plot"), parent=self)
        Plot2DAction = createAction(QtGui.QApplication.translate("pychemqt", "2D Plot"), slot=self.plot2D, parent=self)
        menuPlot.addAction(Plot2DAction)
        Plot3DAction = createAction(QtGui.QApplication.translate("pychemqt", "3D Plot"), slot=self.plot3D, parent=self)
        menuPlot.addAction(Plot3DAction)
        menuPlot.addSeparator()
        Plot_T_s_Action = createAction(QtGui.QApplication.translate("pychemqt", "T-s diagram"), slot=self.showChooseFluid, parent=self)
        menuPlot.addAction(Plot_T_s_Action)
        Plot_T_h_Action = createAction(QtGui.QApplication.translate("pychemqt", "T-h diagram"), slot=self.showChooseFluid, parent=self)
        menuPlot.addAction(Plot_T_h_Action)

        self.addAction(fluidoAction)
        self.addAction(referenciaAction)
        self.addSeparator()
        self.addAction(propAction)
        self.addSeparator()
        self.addAction(menuCalculate.menuAction())
        self.addAction(menuPlot.menuAction())
        self.addSeparator()

        if not (self.parent().currentConfig.has_option("MEoS", "fluid") and self.parent().currentConfig.has_option("MEoS", "reference")):
            menuCalculate.setEnabled(False)
            menuPlot.setEnabled(False)

    def showChooseFluid(self):
        """Show dialog to choose/view fluid"""
        dialog=Ui_ChooseFluid(self.parent().currentConfig)
        if dialog.exec_():
            if not self.parent().currentConfig.has_section("MEoS"):
                self.parent().currentConfig.add_section("MEoS")
            self.parent().currentConfig.set("MEoS", "fluid", str(dialog.lista.currentRow()))
            self.parent().currentConfig.set("MEoS", "eq", str(dialog.eq.currentIndex()))
            self.parent().currentConfig.set("MEoS", "PR", str(dialog.radioPR.isChecked()))
            self.parent().currentConfig.set("MEoS", "Generalized", str(dialog.radioGeneralized.isChecked()))
            self.parent().currentConfig.set("MEoS", "visco", str(dialog.visco.currentIndex()))
            self.parent().currentConfig.set("MEoS", "thermal", str(dialog.thermal.currentIndex()))
            self.checkProperties()
            self.parent().dirty[self.parent().idTab]=True
            self.parent().saveControl()
    
    def showReference(self):
        """Show dialog to choose reference state, for enthalpy and entropy zero state"""
        dialog=Ui_ReferenceState(self.parent().currentConfig)
        if dialog.exec_():
            if not self.parent().currentConfig.has_section("MEoS"):
                self.parent().currentConfig.add_section("MEoS")
            if dialog.OTO.isChecked():
                referencia=["OTO", 298.15, 101325, 0, 0]
            elif dialog.NBP.isChecked():
                Tb=mEoS.__all__[self.parent().currentConfig.getint("MEoS", "fluid")].Tb
                referencia=["NBP", Tb, 101325, 0, 0]
            elif dialog.IIR.isChecked():
                referencia=["IIR", 273.15, 101325, 200, 1]
            elif dialog.ASHRAE.isChecked():
                referencia=["ASHRAE", 233.15, 101325, 0, 0]
            else:
                referencia=["Custom", dialog.T.value, dialog.P.value, dialog.h.value, dialog.s.value]
            self.parent().currentConfig.set("MEoS", "reference", referencia[0])
            self.parent().currentConfig.set("MEoS", "T", str(referencia[1]))
            self.parent().currentConfig.set("MEoS", "P", str(referencia[2]))
            self.parent().currentConfig.set("MEoS", "h", str(referencia[3]))
            self.parent().currentConfig.set("MEoS", "s", str(referencia[4]))
            self.checkProperties()
            self.parent().dirty[self.parent().idTab]=True
            self.parent().saveControl()
            
    def checkProperties(self):
        """Add default properties to show to configuration automatic when
        choose fluid or reference state and properties are not defined"""
        if not self.parent().currentConfig.has_option("MEoS", "properties"):
            self.parent().currentConfig.set("MEoS", "properties", str(Ui_Properties._default))
            self.parent().currentConfig.set("MEoS", "phase", "False")
            self.parent().currentConfig.set("MEoS", "propertiesOrder", range(64))

    def showProperties(self):
        """Show dialog to choose/sort properties to show in tables"""
        dialog=Ui_Properties(self.parent().currentConfig)
        if dialog.exec_():
            if not self.parent().currentConfig.has_section("MEoS"):
                self.parent().currentConfig.add_section("MEoS")
            self.parent().currentConfig.set("MEoS", "properties", str(dialog.properties))
            self.parent().currentConfig.set("MEoS", "phase", str(dialog.checkFase.isChecked()))
            self.parent().currentConfig.set("MEoS", "propertiesOrder", str(dialog.order))
            self.parent().dirty[self.parent().idTab]=True
            self.parent().saveControl()

    def showSaturation(self):
        """Show dialog to define input for a two-phase saturation table calculation"""
        dialog=Ui_Saturation(self.parent().currentConfig)
        if dialog.exec_():
            start=dialog.Inicial.value
            end=dialog.Final.value
            fix=dialog.variableFix.value
            incr=dialog.Incremento.value
            value=arange(start, end, incr)
            if (end-start)%incr == 0:
                value=append(value, end)
            fluid=mEoS.__all__[self.parent().currentConfig.getint("MEoS", "fluid")]

            fluidos=[]
            if dialog.VL.isChecked():
                txt=QtGui.QApplication.translate("pychemqt", "Liquid-Gas Line")
                if dialog.VariarTemperatura.isChecked():
                    for val in value:
                        vconfig = unidades.Temperature(val).str
                        self.parent().statusbar.showMessage("%s: %s=%s, %s" % (fluid.name, "T", vconfig, txt), 3000)
                        fluidos.append(fluid(T=val, x=1))
                elif dialog.VariarPresion.isChecked():
                    for val in value:
                        vconfig = unidades.Temperature(val).str
                        self.parent().statusbar.showMessage("%s: %s=%s, %s" % (fluid.name, "P", vconfig, txt), 3000)
                        fluidos.append(fluid(P=val, x=1))
                elif dialog.VariarXconT.isChecked():
                    fconfig = unidades.Temperature(fix).str
                    for val in value:
                        self.parent().statusbar.showMessage("%s: %s=%0.2f %s=%s, %s" % (fluid.name, "T", fconfig, "x", val, txt), 3000)
                        fluidos.append(fluid(T=fix, x=val))
                elif dialog.VariarXconP.isChecked():
                    fconfig = unidades.Temperature(fix).str
                    for val in value:
                        self.parent().statusbar.showMessage("%s: %s=%0.2f %s=%s, %s" % (fluid.name, "P", fconfig, "x", val, txt), 3000)
                        fluidos.append(fluid(P=fix, x=val))

            else:
                if dialog.SL.isChecked():
                    func=fluid._Melting_Pressure
                    txt=QtGui.QApplication.translate("pychemqt", "Melting Line")
                elif dialog.SV.isChecked():
                    func=fluid._Sublimation_Pressure
                    txt=QtGui.QApplication.translate("pychemqt", "Sublimation Line")

                if dialog.VariarTemperatura.isChecked():
                    for val in value:
                        p=func(val)
                        fluidos.append(fluid(T=val, P=p.MPa))
                        self.parent().statusbar.showMessage("%s: %s=%0.2f, %s" % (fluid.name, "T", val, txt), 3000)
                else:
                    for p in value:
                        T=fsolve(lambda T: p-func(T), fluid.Tt)
                        fluidos.append(fluid(T=T, P=p*1e-6))
                        self.parent().statusbar.showMessage("%s: %s=%0.2f, %s" % (fluid.name, "P", p, txt), 3000)

            title=QtGui.QApplication.translate("pychemqt", "Table %s: %s changing %s" %(fluid.name, txt, "T"))
            self.addTable(fluidos, title)

    def showIsoproperty(self):
        dialog=Ui_Isoproperty(self.parent())
        if dialog.exec_():
            self.parent().updateStatus(QtGui.QApplication.translate("pychemqt", "Launch MEoS Isoproperty calculation..."))
            indice1=dialog.fix.currentIndex()
            indice2=dialog.vary.currentIndex()
            var1=meos.keys[indice1]
            keys=meos.keys[:]
            del keys[indice1]
            var2=keys[indice2]
            if var1=="P":
                value1=dialog.variableFix.value.MPa
            else:
                value1=dialog.variableFix.value
            if var2=="P":
                start=dialog.Inicial.value.MPa
                end=dialog.Final.value.MPa
                incr=dialog.Incremento.value.MPa
            else:
                start=dialog.Inicial.value
                end=dialog.Final.value
                incr=dialog.Incremento.value
            value2=arange(start, end, incr)
            if (end-start)%incr == 0:
                value2=append(value2, end)
            fluid=mEoS.__all__[self.parent().currentConfig.getint("MEoS", "fluid")]
            kwarg={}
            for key in ("eq", "visco", "thermal"):
                kwarg[key]=self.parent().currentConfig.getint("MEoS", key)
            fluidos=[]
            for v2 in value2:
                kwarg[var1]=value1
                kwarg[var2]=v2
                fluidos.append(fluid(**kwarg))
                self.parent().statusbar.showMessage("%s: %s=%0.2f, %s=%0.2f" % (fluid.name, var1, value1, var2, v2), 3000)
            title=QtGui.QApplication.translate("pychemqt", "%s: %s=%0.2f %s changing %s" %(fluid.formula, var1, value1, dialog.unidades[indice1].text(), meos.propiedades[indice2]))
            self.addTable(fluidos, title)


    def plot2D(self):
        dialog=Plot2D(self.parent())
        if dialog.exec_():
            self.parent().progressBar.setVisible(True)
            i=dialog.ejeX.currentIndex()
            j=dialog.ejeY.currentIndex()
            if j>=i:
                j+=1
            xini=dialog.ejeX_min[i].value
            xfin=dialog.ejeX_max[i].value
            yini=dialog.ejeY_min[j].value
            yfin=dialog.ejeY_max[j].value
            c1=dialog.var[i][0]
            c2=dialog.var[j][0]

            fluid=mEoS.__all__[self.parent().currentConfig.getint("MEoS", "fluid")]
            sufx=configSufx(fluid, self.parent())
            title=QtGui.QApplication.translate("pychemqt", "Plot %s: %s=f(%s) %s" %(fluid.formula, dialog.ejeY.currentText(), dialog.ejeX.currentText(), sufx))
            grafico=PlotMEoS(dim=2, parent=self.parent())
            grafico.setWindowTitle(title)
            var=configVariables(self.parent())
            grafico.plot.ax.set_xlabel(var[str(dialog.ejeX.currentText())])
            grafico.plot.ax.set_ylabel(var[str(dialog.ejeY.currentText())])
            grafico.plot.ax.set_title("")
            grafico.plot.ax.c1=c1
            grafico.plot.ax.c2=c2
            grafico.plot.ax.property=None

            calcularSaturacion(self.parent().Preferences, grafico, fluid, dialog.metodo, xini, xfin, yini, yfin, c1, c2)
            calcularIsolineas(self.parent().Preferences, grafico, fluid, dialog.metodo, xini, xfin, yini, yfin, c1, c2)

            grafico.plot.ax.set_xlim(xini, xfin)
            grafico.plot.ax.set_ylim(yini, yfin)
            grafico.plot.ax.grid(self.parent().Preferences.getboolean("MEOS", "grid"))
            if dialog.ejeX_escala.isChecked():
                grafico.plot.ax.set_xscale("log")
            if dialog.ejeY_escala.isChecked():
                grafico.plot.ax.set_yscale("log")
            self.parent().centralwidget.currentWidget().addSubWindow(grafico)
            grafico.show()
            self.parent().progressBar.setVisible(False)

    def plot3D(self):
        dialog=Plot3D(self.parent())
        if dialog.exec_():
            self.parent().progressBar.setVisible(True)
            i, j=dialog.currentIndex()
            xini=dialog.abscisaInicio[i].value
            xfin=dialog.abscisaFin[i].value
            xsalto=dialog.abscisaIntervalo[i].value
            xn=int((xfin-xini)/xsalto+1)
            yini=dialog.ordenadaInicio[j].value
            yfin=dialog.ordenadaFin[j].value
            ysalto=dialog.ordenadaIntervalo[j].value
            yn=int((yfin-yini)/ysalto+1)
            xi=arange(xini, xfin, xsalto)
            if (xfin-xini)/xsalto==float(int((xfin-xini)/xsalto)):
                xi=concatenate((xi, [xfin]))
            yi=arange(yini, yfin, ysalto)
            if (yfin-yini)/ysalto==float(int((yfin-yini)/ysalto)):
                yi=concatenate((yi, [yfin]))

            #python 2.6 compatibility
            inv_dict={}
            for k, v in configVariables(self.parent()).items():
                inv_dict[v]=k
            #inv_dict = {v:k for k, v in configVariables(self.parent()).items()}
            property=inv_dict[dialog.variableTabla.currentText()]
            c1, c2=map(str, dialog.ejesTabla.currentText().split(","))

            fluid=mEoS.__all__[self.parent().currentConfig.getint("MEoS", "fluid")]
            sufx=configSufx(fluid, self.parent())
            title=QtGui.QApplication.translate("pychemqt", "Plot %s: %s=f(%s) %s" %(fluid.formula, property, dialog.ejesTabla.currentText(), sufx))
            labels=[dialog.label_ejeX.text(), dialog.label_ejeY.text(), "z"]
            grafico=PlotMEoS(dim=3, parent=self.parent())
            grafico.setWindowTitle(title)
            xdata, ydata, zdata=calculate(self.parent(), xi, yi, c1, c2, property, dialog)
            grafico.plot.plot_3D(labels, xdata, ydata, zdata, self.parent().Preferences)
            calcularSaturacion(self.parent().Preferences, grafico, fluid, dialog.metodo, xini, xfin, yini, yfin, c1, c2, property)
            calcularIsolineas(self.parent().Preferences, grafico, fluid, dialog.metodo, xini, xfin, yini, yfin, c1, c2, property)

            grafico.plot.ax.set_title("")
            grafico.plot.ax.set_xlim3d(xini, xfin)
            grafico.plot.ax.set_ylim3d(yini, yfin)
            grafico.plot.ax.set_zlim3d(min(zdata), max(zdata))
            grafico.plot.ax.grid(self.parent().Preferences.getboolean("MEOS", "grid"))
            grafico.plot.ax.c1=c1
            grafico.plot.ax.c2=c2
            grafico.plot.ax.property=property

            self.parent().centralwidget.currentWidget().addSubWindow(grafico)
            grafico.show()
            self.parent().progressBar.setVisible(False)






    def addPlot(self, labels, title, xdata, ydata, zdata):
        grafico=PlotMEoS(dim=3, parent=self.parent())
        grafico.setWindowTitle(title)
        grafico.plot.plot_3D(labels, xdata, ydata, zdata, self.parent().Preferences)
        self.parent().centralwidget.currentWidget().addSubWindow(grafico)
        grafico.show()

    def addTable(self, fluidos, title, tabla=None):
        tabla=createTabla(self.parent(), title, fluidos)
        self.parent().centralwidget.currentWidget().addSubWindow(tabla)
        tabla.show()

    def addTableSpecified(self):
        name=mEoS.__all__[self.parent().currentConfig.getint("MEoS", "fluid")].formula
        title="%s: %s" % (name, QtGui.QApplication.translate("pychemqt", "Specified state points"))
        tabla=createTabla(self.parent(), title)
        self.parent().centralwidget.currentWidget().addSubWindow(tabla)
        tabla.show()


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
        self.radioGeneralized=QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "Use generalizated expression"))
        gridLayout.addWidget(self.radioGeneralized,3,1,1,2)
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
            self.radioGeneralized.setChecked(config.getboolean("MEoS", "Generalized"))
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
            values=eval(config.get("MEoS", "properties"))
            fase=config.getboolean("MEoS", "phase")
            self.order = eval(config.get("MEoS", "propertiesOrder"))
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
        print self.order

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
        print self.order

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
            value.append(str(self.listaDisponibles.cellWidget(i, 0).isChecked()))
        return value

    def comprobarBotones(self, fila):
        """Check if button are enabled or disabled"""
        self.ButtonArriba.setEnabled(fila>=1)
        self.ButtonAbajo.setEnabled(fila<self.listaDisponibles.rowCount()-1)


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

        self.VL.setChecked(True)
        self.VariarTemperatura.setChecked(True)
        self.updateVary()
        self.VL.toggled.connect(self.updateVary)

        if config:
            fluido=mEoS.__all__[config.getint("MEoS", "fluid")]
            if fluido._melting or fluido._Melting_Pressure != meos.MEoS._Melting_Pressure:
                self.SL.setEnabled(True)
            else:
                self.SL.setEnabled(False)
            if fluido._sublimation or fluido._Sublimation_Pressure != meos.MEoS._Sublimation_Pressure:
                self.SV.setEnabled(True)
            else:
                self.SV.setEnabled(False)

    def updateVary(self):
        """Update state for option to choose for properties to change"""
        self.VariarXconP.setEnabled(self.VL.isChecked())
        self.VariarXconT.setEnabled(self.VL.isChecked())
        self.VariarTemperatura.setChecked(not self.VL.isChecked())

    def updateVar(self, bool):
        """Update input values units and text"""
        if bool:
            self.Inicial.deleteLater()
            self.Final.deleteLater()
            self.Incremento.deleteLater()
            if self.sender() == self.VariarXconT:
                self.labelFix.setVisible(True)
                self.labelFix.setText(unidades.Temperature.__title__)
                self.variableFix.deleteLater()
                self.variableFix=Entrada_con_unidades(unidades.Temperature)
                self.layout().addWidget(self.variableFix,4,4)
                unidadVariable=float
                self.labelinicial.setText(QtGui.QApplication.translate("pychemqt", "Initial quality"))
                self.labelfinal.setText(QtGui.QApplication.translate("pychemqt", "Final quality"))

            elif self.sender() == self.VariarXconP:
                self.labelFix.setVisible(True)
                self.labelFix.setText(unidades.Pressure.__title__)
                self.variableFix.deleteLater()
                self.variableFix=Entrada_con_unidades(unidades.Pressure)
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

            self.Inicial=Entrada_con_unidades(unidadVariable)
            self.Final=Entrada_con_unidades(unidadVariable)
            if unidadVariable==unidades.Temperature:
                unidadDelta=unidades.DeltaT
            elif unidadVariable==unidades.Pressure:
                unidadDelta=unidades.DeltaP
            else:
                unidadDelta=unidadVariable

            self.Incremento=Entrada_con_unidades(unidadDelta)
            self.layout().addWidget(self.Inicial,4,2)
            self.layout().addWidget(self.Final,5,2)
            self.layout().addWidget(self.Incremento,6,2)


class Ui_Isoproperty(QtGui.QDialog):
    """Dialog to define input for isoproperty table calculations"""
    propiedades=[QtGui.QApplication.translate("pychemqt", "Temperature"),
                                QtGui.QApplication.translate("pychemqt", "Pressure"),
                                QtGui.QApplication.translate("pychemqt", "Density"),
                                QtGui.QApplication.translate("pychemqt", "Volume"),
                                QtGui.QApplication.translate("pychemqt", "Enthalpy"),
                                QtGui.QApplication.translate("pychemqt", "Entropy")]
    unidades=[unidades.Temperature, unidades.Pressure, unidades.Density, unidades.SpecificVolume, unidades.Enthalpy, unidades.SpecificHeat]

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
        propiedades=self.propiedades[:3]
        if indice<3:
            del propiedades[indice]
        elif indice==3:
            del propiedades[2]
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


def method(parent):
    if parent.currentConfig.getfloat("MEoS", "fluid"):
        return "meos"
    elif parent.Preferences.getboolean("MEOS", "iapws") and parent.Preferences.getboolean("MEOS", "freesteam"):
        return "freesteam"
    else:
        return "iapws"

def configVariables(parent):
    metodo=method(parent)
    if metodo=="freesteam":
        variables={"p": QtGui.QApplication.translate("pychemqt", "Pressure"),
               "T": QtGui.QApplication.translate("pychemqt", "Temperature"),
               "v": QtGui.QApplication.translate("pychemqt", "Specific Volume"),
               "h": QtGui.QApplication.translate("pychemqt", "Enghalpy"),
               "s": QtGui.QApplication.translate("pychemqt", "Entropy"),
               "u": QtGui.QApplication.translate("pychemqt", "Internal Energy"),
               "cp": QtGui.QApplication.translate("pychemqt", "Cp"),
               "cv": QtGui.QApplication.translate("pychemqt", "Cv"),
               "rho": QtGui.QApplication.translate("pychemqt", "Density"),
               "k": QtGui.QApplication.translate("pychemqt", "Thermal Conductivity"),
               "mu": QtGui.QApplication.translate("pychemqt", "Viscosity"),
               "x": QtGui.QApplication.translate("pychemqt", "Quality"),
               "w": QtGui.QApplication.translate("pychemqt", "Speed of sound")}
    elif metodo=="iapws":
        variables=iapws.properties
    else:
        variables=meos.properties
    return variables

def configUnidades(parent):
    metodo=method(parent)
    if metodo=="freesteam":
        variables=[("p", unidades.Pressure, None),
               ("T", unidades.Temperature, None),
               ("v", unidades.SpecificVolume, None),
               ("h", unidades.Enthalpy, None),
               ("s", unidades.SpecificHeat, "Entropy"),
               ("u", unidades.Enthalpy, None),
               ("cp", unidades.SpecificHeat, None),
               ("cv", unidades.SpecificHeat, None),
               ("rho", unidades.Density, None),
               ("k", unidades.ThermalConductivity, None),
               ("mu", unidades.Viscosity, None),
               ("x", float, None),
               ("w", unidades.Speed, None)]

    elif metodo=="iapws":
        variables=iapws.units
    else:
        variables=meos.units
    return variables


def configSufx(fluid, parent):
    metodo=method(parent)
    if metodo=="freesteam":
        sufx=QtGui.QApplication.translate("pychemqt", "using")+" freesteam"
    elif metodo=="iapws":
        sufx=QtGui.QApplication.translate("pychemqt", "using")+" iapws97"
    elif fluid.formula == "H2O":
        sufx=QtGui.QApplication.translate("pychemqt", "using")+" iapws95"
    else:
        sufx=""

    return sufx


def get_propiedades(parent):
    metodo=method(parent)
    if metodo=="freesteam":
        variables=configVariables(parent)
        propiedades=[]
        keys=[]
        for key, propiedad in variables.iteritems():
            if key != "w":
                propiedades.append(propiedad)
                keys.append(key)
    else:
        booleanos = eval(parent.currentConfig.get("MEoS", "properties"))
        order = eval(parent.currentConfig.get("MEoS", "propertiesOrder"))
        propiedades=[]
        keys=[]
        for indice, bool in zip(order, booleanos):
            if bool=="True":
                propiedades.append(meos.propiedades[indice])
                keys.append(meos.keys[indice])
    return propiedades, keys

def createTabla(parent, title, fluidos=None):
    propiedades, keys=get_propiedades(parent)
    if fluidos:
        for i, key in enumerate(keys):
            sufx = fluidos[0].__getattribute__(key).text()
            if not sufx:
                sufx = "[-]"
            propiedades[i]=propiedades[i]+"\n"+sufx
        tabla = TablaMEoS(len(propiedades), horizontalHeader=propiedades, stretch=False, readOnly=True, parent=parent)
        data=[]
        for fluido in fluidos:
            fila=[]
            for key in keys:
                fila.append(fluido.__getattribute__(key).config())
            data.append(fila)
        tabla.setMatrix(data)
    else:
        columnInput=[]
        for key in keys:
            if key in ["P", "T", "x", "rho", "v", "h", "s"]:
                columnInput.append(False)
            else:
                columnInput.append(True)
        tabla = TablaMEoS(len(propiedades), horizontalHeader=propiedades, filas=1, dinamica=True, stretch=False, columnReadOnly=columnInput, parent=parent)

    if fluidos:
        sufx=configSufx(fluidos[0], parent)
    else:
        sufx=""
    prefix=QtGui.QApplication.translate("pychemqt", "Table")
    tabla.setWindowTitle(prefix+": "+title+" "+sufx)
    tabla.resizeColumnsToContents()
    return tabla


def calculate(mainwindow, xi, yi, c1, c2, property, dialog):
    """Calculo de mesh de datos en los graficos 3D
    soporte terminado para freesteam, iapws y mEoS"""
    xdata,ydata = meshgrid(xi, yi)
    zdata = zeros(xdata.shape)
    matriz=[]
    contador=0
    progressfactor=100./len(xi)/len(yi)
    kwarg={}
    if dialog.metodo=="freesteam":
        func=[freesteam.steam_pT, freesteam.steam_ph, freesteam.steam_ps, freesteam.steam_pv, freesteam.steam_Ts, freesteam.steam_Tx][dialog.ejesTabla.currentIndex()]

        for x in xi:
            fila=[]
            for y in yi:
                mainwindow.progressBar.setValue(contador*progressfactor)
                mainwindow.statusbar.showMessage("%s %s=%f,%s=%f" % (QtGui.QApplication.translate("pychemqt", "Calculating..."), c1, x, c2, y), 3000)
                QtGui.QApplication.processEvents()
                fila.append(func(x, y))
                contador+=1
            matriz.append(fila)

    else:
        if dialog.metodo=="iapws":
            func=iapws.IAPWS97
        else:
            func=mEoS.__all__[mainwindow.currentConfig.getint("MEoS", "fluid")]
            for key in ("eq", "visco", "thermal"):
                kwarg[key]=mainwindow.currentConfig.getint("MEoS", key)

        factor={"P": 1.0e6, "s": 1.0e3, "h": 1.0e3, "u": 1.0e3}
        factor1=factor.get(c1, 1.)
        factor2=factor.get(c2, 1.)
        for x in xi:
            fila=[]
            for y in yi:
                mainwindow.progressBar.setValue(contador*progressfactor)
                QtGui.QApplication.processEvents()
                mainwindow.statusbar.showMessage("%s %s=%f,%s=%f" % (QtGui.QApplication.translate("pychemqt", "Calculating..."), c1, x, c2, y), 3000)
                kwarg[c1]=x/factor1
                kwarg[c2]=y/factor2
                fluido=func(**kwarg)
                fila.append(fluido)
                contador+=1
            matriz.append(fila)

    for i, fila in enumerate(matriz):
        for j, fluid in enumerate(fila):
            zdata[j, i]=fluid.__getattribute__(property)

    return xdata, ydata, zdata


def calcularSaturacion(Preferences, grafico, fluid, metodo, xini, xfin, yini, yfin, c1, c2, property=None):
    """Método que calcula datos de la línea de saturación
    Soporte para freesteam y iapws"""
    T = linspace(fluid.Tt, fluid.Tc, 100)
    factor1, factor2, factorProperty=1., 1., 1.
    if metodo=="freesteam":
        func=freesteam.region4_Tx
    elif metodo=="iapws":
        func=iapws.IAPWS97_Tx
        factor={"P": 1.0e6, "s": 1.0e3, "h": 1.0e3, "u": 1.0e3}
        factor1=factor.get(c1, 1.)
        factor2=factor.get(c2, 1.)
        if property:
            factorProperty=factor.get(property, 1.)
    else:
        pass

    format={}
    format["ls"]=Preferences.get("MEOS", "saturationlineStyle")
    format["lw"]=Preferences.getfloat("MEOS", "saturationlineWidth")
    format["color"]=Preferences.get("MEOS", "saturationColor")
    format["marker"]=Preferences.get("MEOS", "saturationmarker")
    format["ms"]=3

    label=[QtGui.QApplication.translate("pychemqt", "Liquid Saturation Line"),
                QtGui.QApplication.translate("pychemqt", "Vapor Saturation Line")]
    for fase in (0, 1):
        fluidos=[func(Ti, fase) for Ti in T]
        xsat=[fluido.__getattribute__(c1)*factor1 for fluido in fluidos]
        ysat=[fluido.__getattribute__(c2)*factor2 for fluido in fluidos]
        if property:
            zsat=[fluido.__getattribute__(property)*factorProperty for fluido in fluidos]
        else:
            zsat=None
        plotLine(grafico, xsat, ysat, zsat, xini, xfin, yini, yfin, format, label[fase])
        grafico.plot.ax.lines[-1].__setattr__("fluids", fluidos)


def calcularIsolineas(Preferences, grafico, fluid, metodo, xini, xfin, yini, yfin, c1, c2, property=None, add=None):
    """Método que calcula datos para isolineas
    add: parametro opcional que permite añadir lineas especificadas en forma de array:
        indice de la isolinea
        valor fijo
        soporte completo para freesteam e iapws"""

    factor1, factor2, factorProperty=1., 1., 1.

    if metodo=="freesteam":

        isolineas=[("Isotherm", "T", freesteam.steam_pT, QtGui.QApplication.translate("pychemqt", "Isotherm")),
                ("Isobar", "p", freesteam.steam_pT, QtGui.QApplication.translate("pychemqt", "Isobar")),
                ("Isoenthalpic", "h", freesteam.steam_ph, QtGui.QApplication.translate("pychemqt", "Isoenthalpic")),
                ("Isoentropic", "s", freesteam.steam_ps, QtGui.QApplication.translate("pychemqt", "Isoentropic")),
                ("Isochor", "v", freesteam.steam_pv, QtGui.QApplication.translate("pychemqt", "Isochor")),
                ("Isoquality", "x", freesteam.steam_Tx, QtGui.QApplication.translate("pychemqt", "Isoquality"))]
    else:
        factor={"P": 1.0e6, "s": 1.0e3, "h": 1.0e3, "u": 1.0e3}
        factor1=factor.get(c1, 1.)
        factor2=factor.get(c2, 1.)
        if property:
            factorProperty=factor.get(property, 1.)

        if metodo=="iapws":
            isolineas=[("Isotherm", "T", iapws.IAPWS97_PT, QtGui.QApplication.translate("pychemqt", "Isotherm")),
                    ("Isobar", "P", iapws.IAPWS97_PT, QtGui.QApplication.translate("pychemqt", "Isobar")),
#                    ("Isoenthalpic", "h", iapws.IAPWS97_Ph, QtGui.QApplication.translate("pychemqt", "Isoenthalpic")),
#                    ("Isoentropic", "s", iapws.IAPWS97_Ps, QtGui.QApplication.translate("pychemqt", "Isoentropic")),
                    ("Isochor", "v", iapws.IAPWS97_Pv, QtGui.QApplication.translate("pychemqt", "Isochor")),
                    ("Isoquality", "x", iapws.IAPWS97_Tx, QtGui.QApplication.translate("pychemqt", "Isoquality"))]

        else:
            pass


    if add:
        isolineas=[isolineas[add[0]]]

    for name, prop, func, title in isolineas:
        format={}
        format["ls"]=Preferences.get("MEOS", name+"lineStyle")
        format["lw"]=Preferences.getfloat("MEOS", name+"lineWidth")
        format["color"]=Preferences.get("MEOS", name+"Color")
        format["marker"]=Preferences.get("MEOS", name+"marker")
        format["ms"]=3

        factorProperty=factor.get(property, 1.)

        if add:
            lineas=[add[1]]
        else:
            lineas=[]
            if Preferences.getboolean("MEOS", name+"Custom"):
                for i in Preferences.get("MEOS", name+'List').split(','):
                    if i:
                        lineas.append(float(i))
            else:
                start=Preferences.getboolean("MEOS", name+"Start")
                end=Preferences.getboolean("MEOS", name+"End")
                step=Preferences.getboolean("MEOS", name+"Step")
                lineas=arange(start, end, step)

            if prop != "x" and Preferences.getboolean("MEOS", name+"Critic"):
                if metodo=="freesteam":
                    fcri=freesteam.steam_pT
                    Pc=fluid.Pc
                elif metodo=="iapws":
                    fcri=iapws.IAPWS97_PT
                    Pc=fluid.Pc.MPa
                PropCRIT=fcri(Pc, fluid.Tc).__getattribute__(prop)/factorProperty
                lineas.append(PropCRIT)

        if prop=="x":
            X=linspace(fluid.Tt, fluid.Tc, 100)
        elif prop=="P":
            X=linspace(fluid.Tt, 1000, 100)
        else:
            X= logspace(-3, 3, 100)/10
        for linea in lineas:
            xi, yi, zi, fluidos=calcIsolinea(prop, func, linea, X, c1, c2, property, factor1, factor2, factorProperty)
            sucess=plotLine(grafico, xi, yi, zi, xini, xfin, yini, yfin, format, "%s %s=%0.4f" %(title, prop, linea))
            if sucess and prop in (c1, c2):
                grafico.plot.ax.lines[-1].set_visible(False)
            grafico.plot.ax.lines[-1].__setattr__("fluids", fluidos)



def calcIsolinea(prop, func, linea, X, c1, c2, property, factor1, factor2, factorProperty):
    if prop in ["p", "P"]:
        fluidos=[func(linea, Xi) for Xi in X]
    else:
        fluidos=[]
        for Xi in X:
#            print prop, func, linea, Xi
            fluidos.append(func(Xi, linea))
#        fluidos=[func(Xi, linea) for Xi in X]
    xi=[fluido.__getattribute__(c1)*factor1 for fluido in fluidos]
    yi=[fluido.__getattribute__(c2)*factor2 for fluido in fluidos]
    if property:
        zi=[fluido.__getattribute__(property)*factorProperty for fluido in fluidos]
    else:
        zi=None
    return xi, yi, zi, fluidos


def plotLine(grafico, xi=None, yi=None, zi=None, xini=None, xfin=None, yini=None, yfin=None, format=None, label=None):
    if grafico.dim==2:
        list=(xi, yi)
    else:
        list=(xi, yi, zi)
        for i in range(len(xi)-1, -1, -1):
            if xi[i]<xini or xi[i]>xfin or yi[i]<yini or yi[i]>yfin:
                for lista in list:
                    del lista[i]

    if xi:
        grafico.plot.ax.plot(*list, label=label, **format)
        return True




















class EditPlot(QtGui.QWidget):
    """Dialogo de configuracion de datos del grafico"""
    def __init__(self, plotMEoS, mainwindow, parent=None):
        super(EditPlot, self).__init__(parent)
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Edit Plot"))
        layout = QtGui.QGridLayout(self)
        self.plotMEoS=plotMEoS
        self.fig=plotMEoS.plot
        self.mainwindow=mainwindow

        self.lista=QtGui.QListWidget()
        layout.addWidget(self.lista,0,1,1,3)

        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Label")),1,1)
        self.label=QtGui.QLineEdit()
        layout.addWidget(self.label,1,2,1,2)

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

        self.populateList()

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

    def populateList(self):
        """Rellena la lista con los label de los elementos disponibles en el grafico"""
        self.lista.clear()
        for linea in self.fig.ax.lines:
            self.lista.addItem(linea._label)

    def update(self, i):
        """Rellena widgets con los valores del elemento del grafico seleccionado en la lista"""
        linea=self.fig.ax.lines[i]
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
        linea=self.fig.ax.lines[self.lista.currentRow()]
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
        """Añade una isolinea al grafico"""
        dialog=AddLine()
        if dialog.exec_():
            fluid=mEoS.__all__[self.mainwindow.currentConfig.getint("MEoS", "fluid")]
            metodo=method(self.mainwindow)
            Preferences=self.mainwindow.Preferences
            xini, xfin=self.fig.ax.get_xlim()
            yini, yfin=self.fig.ax.get_ylim()
            c1=self.fig.ax.c1
            c2=self.fig.ax.c2
            property=self.fig.ax.property
            prop=dialog.tipo.currentIndex()
            value=dialog.input[prop].value
            calcularIsolineas(Preferences, self.plotMEoS, fluid, metodo, xini, xfin, yini, yfin, c1, c2, property, (prop, value))
            self.fig.draw()
            self.lista.addItem(self.fig.ax.lines[-1].get_label())
            self.lista.setCurrentRow(self.lista.count()-1)

    def remove(self):
        """Elimina el elemento seleccionado en la lista del grafico"""
        del self.fig.ax.lines[self.lista.currentRow()]
        self.populateList()
        self.fig.draw()


class AddLine(QtGui.QDialog):
    """Dialogo de definicion de nuevas isolineas"""
    lineas=[(QtGui.QApplication.translate("pychemqt", "Isotherm"), unidades.Temperature, None),
                (QtGui.QApplication.translate("pychemqt", "Isobar"), unidades.Pressure, None),
                (QtGui.QApplication.translate("pychemqt", "Isoenthalpic"), unidades.Enthalpy, None),
                (QtGui.QApplication.translate("pychemqt", "Isoentropic"), unidades.SpecificHeat, "Entropy"),
                (QtGui.QApplication.translate("pychemqt", "Isochor"), unidades.SpecificVolume, None),
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
    """Dialogo de configuracion de parametros generales del grafico"""
    def __init__(self, fig, parent=None):
        super(EditAxis, self).__init__(parent)
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Edit Axis"))
        layout = QtGui.QGridLayout(self)
        self.fig=fig

        label=QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Title"))
        label.setSizePolicy(QtGui.QSizePolicy.Maximum,QtGui.QSizePolicy.Maximum)
        layout.addWidget(label,1,1)
        self.title=InputFond()
        layout.addWidget(self.title,1,2,1,3)

        groupX=QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "x Axis"))
        layout.addWidget(groupX,2,1,1,2)
        groupXlayout=QtGui.QGridLayout(groupX)
        groupXlayout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Label")),1,1)
        self.labelX=InputFond()
        groupXlayout.addWidget(self.labelX,1,2)
        self.scaleX=QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Logarithmic scale"))
        groupXlayout.addWidget(self.scaleX,2,1,1,2)
        groupXlayout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "from")),3,1)
        self.xmin=Entrada_con_unidades(float)
        groupXlayout.addWidget(self.xmin,3,2)
        groupXlayout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "to")),4,1)
        self.xmax=Entrada_con_unidades(float)
        groupXlayout.addWidget(self.xmax,4,2)

        groupY=QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "y Axis"))
        layout.addWidget(groupY,2,3,1,2)
        groupYlayout=QtGui.QGridLayout(groupY)
        groupYlayout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Label")),1,1)
        self.labelY=InputFond()
        groupYlayout.addWidget(self.labelY,1,2)
        self.scaleY=QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Logarithmic scale"))
        groupYlayout.addWidget(self.scaleY,2,1,1,2)
        groupYlayout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "from")),3,1)
        self.ymin=Entrada_con_unidades(float)
        groupYlayout.addWidget(self.ymin,3,2)
        groupYlayout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "to")),4,1)
        self.ymax=Entrada_con_unidades(float)
        groupYlayout.addWidget(self.ymax,4,2)

        self.gridCheckbox=QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Show Grid"))
        layout.addWidget(self.gridCheckbox,3,1,1,4)
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Close)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox,10,1,1,4)

        self.populate()

        self.title.textChanged.connect(partial(self.update, "title"))
        self.title.colorChanged.connect(partial(self.update, "titlecolor"))
        self.title.fontChanged.connect(partial(self.update, "titlefont"))
        self.labelX.textChanged.connect(partial(self.update, "xlabel"))
        self.labelX.colorChanged.connect(partial(self.update, "xlabelcolor"))
        self.labelX.fontChanged.connect(partial(self.update, "xlabelfont"))
        self.labelY.textChanged.connect(partial(self.update, "ylabel"))
        self.labelY.colorChanged.connect(partial(self.update, "ylabelcolor"))
        self.labelY.fontChanged.connect(partial(self.update, "ylabelfont"))
        self.gridCheckbox.toggled.connect(partial(self.update, "grid"))
        self.scaleX.toggled.connect(partial(self.update, "xscale"))
        self.scaleY.toggled.connect(partial(self.update, "yscale"))
        self.xmin.valueChanged.connect(partial(self.update, "xmin"))
        self.xmax.valueChanged.connect(partial(self.update, "xmax"))
        self.ymin.valueChanged.connect(partial(self.update, "ymin"))
        self.ymax.valueChanged.connect(partial(self.update, "ymax"))

    def populate(self):
        """Rellena widgets con los valores del gráfico"""
        self.title.setText(self.fig.ax.get_title())
        self.title.setColor(QtGui.QColor(self.fig.ax.title.get_color()))
        self.labelX.setText(self.fig.ax.get_xlabel())
        self.labelX.setColor(QtGui.QColor(self.fig.ax.xaxis.get_label().get_color()))
        self.labelY.setText(self.fig.ax.get_ylabel())
        self.labelY.setColor(QtGui.QColor(self.fig.ax.yaxis.get_label().get_color()))
        self.gridCheckbox.setChecked(self.fig.ax.get_xgridlines()[0].get_visible())
        self.scaleX.setChecked(self.fig.ax.get_xscale()=="log")
        self.scaleY.setChecked(self.fig.ax.get_yscale()=="log")
        xmin, xmax=self.fig.ax.get_xlim()
        self.xmin.setValue(xmin)
        self.xmax.setValue(xmax)
        ymin, ymax=self.fig.ax.get_ylim()
        self.ymin.setValue(ymin)
        self.ymax.setValue(ymax)

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
                    "grid": self.fig.ax.grid,
                    "xscale": self.fig.ax.set_xscale,
                    "yscale": self.fig.ax.set_yscale}

        if key in ("xscale", "yscale"):
            if value:
                value="log"
            else:
                value="linear"
        if key in ("titlecolor", "xlabelcolor", "ylabelcolor"):
            value=str(value)
        if key in ("titlefont", "xlabelfont", "ylabelfont"):
            value=self.convertFont(value)

        if key in ("xmin", "xmax"):
            xmin=self.xmin.value
            xmax=self.xmax.value
            self.fig.ax.set_xlim(xmin, xmax)
        elif key in ("ymin", "ymax"):
            ymin=self.ymin.value
            ymax=self.ymax.value
            self.fig.ax.set_ylim(ymin, ymax)
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


class Plot2D(QtGui.QDialog):
    """Widget de configuracion inicial de graficos 2D"""

    def __init__(self, parent=None):
        super(Plot2D, self).__init__(parent)
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Setup 2D Plot"))
        layout = QtGui.QGridLayout(self)
        self.metodo=method(self.parent())

        self.var=configUnidades(self.parent())
        self.unit=[var[0] for var in self.var]
        group_Ejex=QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Axis X"))
        layout.addWidget(group_Ejex,1,1)
        layout_GroupX=QtGui.QGridLayout(group_Ejex)
        self.ejeX = QtGui.QComboBox()
        layout_GroupX.addWidget(self.ejeX,1,1)
        self.ejeX_escala=QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Logarithmic scale"))
        layout_GroupX.addWidget(self.ejeX_escala,2,1)
        layout_GroupX.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Start")),1,2)
        layout_GroupX.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "End")),2,2)
        self.ejeX_min=[]
        self.ejeX_max=[]
        for nombre, unidad, magnitud in self.var:
            self.ejeX_min.append(Entrada_con_unidades(unidad, magnitud))
            layout_GroupX.addWidget(self.ejeX_min[-1],1,3,1,2)
            self.ejeX_max.append(Entrada_con_unidades(unidad, magnitud))
            layout_GroupX.addWidget(self.ejeX_max[-1],2,3,1,2)
            self.ejeX.addItem(nombre)

        group_Ejey=QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Axis Y"))
        layout.addWidget(group_Ejey,2,1)
        layout_GroupY=QtGui.QGridLayout(group_Ejey)
        self.ejeY = QtGui.QComboBox()
        layout_GroupY.addWidget(self.ejeY,1,1)
        self.ejeY_escala=QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Logarithmic scale"))
        layout_GroupY.addWidget(self.ejeY_escala,2,1)
        layout_GroupY.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Start")),1,2)
        layout_GroupY.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "End")),2,2)
        self.ejeY_min=[]
        self.ejeY_max=[]
        for nombre, unidad, magnitud in self.var:
            self.ejeY_min.append(Entrada_con_unidades(unidad, magnitud))
            layout_GroupY.addWidget(self.ejeY_min[-1],1,3,1,2)
            self.ejeY_max.append(Entrada_con_unidades(unidad, magnitud))
            layout_GroupY.addWidget(self.ejeY_max[-1],2,3,1,2)
            self.ejeY.addItem(nombre)

        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok|QtGui.QDialogButtonBox.Cancel)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox,10,1)

        self.ejeXChanged(0)
        self.ejeX.currentIndexChanged.connect(self.ejeXChanged)
        self.ejeY.currentIndexChanged.connect(self.ejeYChanged)
        self.ejeYChanged()

    def ejeXChanged(self, int):
        """Rellena las variables disponibles para el ejeY en el gráfico 2D, todos menos el que este activo en el ejeX"""
        self.ejeY.clear()
        Ejes2D=self.unit[:]
        del Ejes2D[int]
        for nombre in Ejes2D:
            self.ejeY.addItem(nombre)

        for i in self.ejeX_min:
            i.setVisible(False)
        self.ejeX_min[int].setVisible(True)
        for i in self.ejeX_max:
            i.setVisible(False)
        self.ejeX_max[int].setVisible(True)


    def ejeYChanged(self):
        try:
            int=self.unit.index(self.ejeY.currentText())
        except ValueError:
            int=0

        for i in self.ejeY_min:
            i.setVisible(False)
        self.ejeY_min[int].setVisible(True)
        for i in self.ejeY_max:
            i.setVisible(False)
        self.ejeY_max[int].setVisible(True)


class Plot3D(QtGui.QDialog):
    """Widget de configuracion inicial de graficos 3D"""

    def __init__(self, parent=None):
        super(Plot3D, self).__init__(parent)
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Setup 3D Plot"))
        layout = QtGui.QGridLayout(self)

        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Axis")),1,1)
        self.ejesTabla = QtGui.QComboBox()
        self.ejesTabla.setToolTip(QtGui.QApplication.translate("pychemqt", "p\tPressure\nT\tTemperature\nh\tEnthalpy\ns\tEntropy\nv\tSpecific Volume\nx\tQuality"))
        layout.addWidget(self.ejesTabla,1,2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Variable")),2,1)
        self.variableTabla = QtGui.QComboBox()
        self.variableTabla.setToolTip(QtGui.QApplication.translate("pychemqt", "Define variable to draw"))
        layout.addWidget(self.variableTabla,2,2,1,3)
        self.label_ejeX = QtGui.QLabel()
        layout.addWidget(self.label_ejeX,4,1,1,1)
        self.label_ejeY = QtGui.QLabel()
        layout.addWidget(self.label_ejeY,5,1,1,1)
        label=QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Start"))
        label.setAlignment(QtCore.Qt.AlignCenter|QtCore.Qt.AlignBottom)
        layout.addWidget(label,3,2)
        label=QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "End"))
        label.setAlignment(QtCore.Qt.AlignCenter|QtCore.Qt.AlignBottom)
        layout.addWidget(label,3,3)
        label=QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Step")+"         ")
        label.setAlignment(QtCore.Qt.AlignCenter|QtCore.Qt.AlignBottom)
        layout.addWidget(label,3,4)
        self.metodo=method(self.parent())
        if self.metodo=="freesteam":
            Ejes=["p,T", "p,h", "p,s", "p,v", "T,s", "T,x"]

            self.abscisaInicio = [Entrada_con_unidades(unidades.Pressure, texto=False), Entrada_con_unidades(unidades.Temperature, texto=False)]
            self.abscisaFin = [Entrada_con_unidades(unidades.Pressure, texto=False), Entrada_con_unidades(unidades.Temperature, texto=False)]
            self.abscisaIntervalo = [Entrada_con_unidades(unidades.DeltaP), Entrada_con_unidades(unidades.DeltaT)]
            self.ordenadaInicio = [Entrada_con_unidades(unidades.Temperature, texto=False), Entrada_con_unidades(unidades.Enthalpy, texto=False), Entrada_con_unidades(unidades.SpecificHeat, "Entropy", texto=False, title=QtGui.QApplication.translate("pychemqt", "Entropy")), Entrada_con_unidades(unidades.SpecificVolume, texto=False), Entrada_con_unidades(unidades.Dimensionless, texto=False, title=QtGui.QApplication.translate("pychemqt", "Quality"))]
            self.ordenadaFin = [Entrada_con_unidades(unidades.Temperature, texto=False), Entrada_con_unidades(unidades.Enthalpy, texto=False), Entrada_con_unidades(unidades.SpecificHeat, "Entropy", texto=False), Entrada_con_unidades(unidades.SpecificVolume, texto=False), Entrada_con_unidades(unidades.Dimensionless, texto=False)]
            self.ordenadaIntervalo = [Entrada_con_unidades(unidades.DeltaT), Entrada_con_unidades(unidades.Enthalpy), Entrada_con_unidades(unidades.SpecificHeat, "Entropy"), Entrada_con_unidades(unidades.SpecificVolume), Entrada_con_unidades(unidades.Dimensionless, texto=False)]

        elif self.metodo=="iapws":
            Ejes=["T,P", "P,h", "P,s", "h,s", "T,x", "P,x"]

            self.abscisaInicio = [Entrada_con_unidades(unidades.Temperature, texto=False), Entrada_con_unidades(unidades.Pressure, texto=False), Entrada_con_unidades(unidades.Enthalpy, texto=False)]
            self.abscisaFin = [Entrada_con_unidades(unidades.Temperature, texto=False), Entrada_con_unidades(unidades.Pressure, texto=False), Entrada_con_unidades(unidades.Enthalpy, texto=False)]
            self.abscisaIntervalo = [Entrada_con_unidades(unidades.DeltaT), Entrada_con_unidades(unidades.DeltaP), Entrada_con_unidades(unidades.Enthalpy)]
            self.ordenadaInicio = [Entrada_con_unidades(unidades.Pressure, texto=False), Entrada_con_unidades(unidades.Enthalpy, texto=False), Entrada_con_unidades(unidades.SpecificHeat, "Entropy", texto=False, title=QtGui.QApplication.translate("pychemqt", "Entropy")), Entrada_con_unidades(unidades.Dimensionless, texto=False, title=QtGui.QApplication.translate("pychemqt", "Quality"))]
            self.ordenadaFin = [Entrada_con_unidades(unidades.Pressure, texto=False), Entrada_con_unidades(unidades.Enthalpy, texto=False), Entrada_con_unidades(unidades.SpecificHeat, "Entropy", texto=False), Entrada_con_unidades(unidades.Dimensionless, texto=False)]
            self.ordenadaIntervalo = [Entrada_con_unidades(unidades.DeltaP), Entrada_con_unidades(unidades.Enthalpy), Entrada_con_unidades(unidades.SpecificHeat, "Entropy"), Entrada_con_unidades(unidades.Dimensionless, texto=False)]

        else:
            Ejes=["T,P", "T,rho", "T,v", "T,h", "T,s", "T,u", "T,x", "P,rho", "P,v", "P,h", "P,s", "P,u", "P,x", "rho,h", "rho,s", "rho,u", "v,h", "v,s", "v,u", "h,s", "h,u", "s,u"]
            self.abscisaInicio = [Entrada_con_unidades(unidades.Temperature, texto=False), Entrada_con_unidades(unidades.Pressure, texto=False), Entrada_con_unidades(unidades.Density, texto=False), Entrada_con_unidades(unidades.SpecificVolume, texto=False), Entrada_con_unidades(unidades.Enthalpy, texto=False), Entrada_con_unidades(unidades.SpecificHeat, "Entropy", texto=False, title=QtGui.QApplication.translate("pychemqt", "Entropy"))]
            self.abscisaFin = [Entrada_con_unidades(unidades.Temperature, texto=False), Entrada_con_unidades(unidades.Pressure, texto=False), Entrada_con_unidades(unidades.Density, texto=False), Entrada_con_unidades(unidades.SpecificVolume, texto=False), Entrada_con_unidades(unidades.Enthalpy, texto=False), Entrada_con_unidades(unidades.SpecificHeat, "Entropy", texto=False)]
            self.abscisaIntervalo = [Entrada_con_unidades(unidades.DeltaT), Entrada_con_unidades(unidades.DeltaP), Entrada_con_unidades(unidades.Density), Entrada_con_unidades(unidades.SpecificVolume), Entrada_con_unidades(unidades.Enthalpy), Entrada_con_unidades(unidades.SpecificHeat, "Entropy")]
            self.ordenadaInicio = [Entrada_con_unidades(unidades.Pressure, texto=False), Entrada_con_unidades(unidades.Density, texto=False), Entrada_con_unidades(unidades.SpecificVolume, texto=False), Entrada_con_unidades(unidades.Enthalpy, texto=False), Entrada_con_unidades(unidades.SpecificHeat, "Entropy", texto=False), Entrada_con_unidades(unidades.Enthalpy, texto=False, title=QtGui.QApplication.translate("pychemqt", "Internal Energy")), Entrada_con_unidades(unidades.Dimensionless, texto=False, title=QtGui.QApplication.translate("pychemqt", "Quality"))]
            self.ordenadaFin = [Entrada_con_unidades(unidades.Pressure, texto=False), Entrada_con_unidades(unidades.Density, texto=False), Entrada_con_unidades(unidades.SpecificVolume, texto=False), Entrada_con_unidades(unidades.Enthalpy, texto=False), Entrada_con_unidades(unidades.SpecificHeat, "Entropy", texto=False), Entrada_con_unidades(unidades.Enthalpy, texto=False), Entrada_con_unidades(unidades.Dimensionless, texto=False)]
            self.ordenadaIntervalo = [Entrada_con_unidades(unidades.DeltaP), Entrada_con_unidades(unidades.Density), Entrada_con_unidades(unidades.SpecificVolume), Entrada_con_unidades(unidades.Enthalpy), Entrada_con_unidades(unidades.SpecificHeat, "Entropy"), Entrada_con_unidades(unidades.Enthalpy), Entrada_con_unidades(unidades.Dimensionless, texto=False)]

        for widget in self.abscisaInicio:
            layout.addWidget(widget,4,2)
        for widget in self.abscisaFin:
            layout.addWidget(widget,4,3)
        for widget in self.abscisaIntervalo:
            layout.addWidget(widget,4,4)
        for widget in self.ordenadaInicio:
            layout.addWidget(widget,5,2)
        for widget in self.ordenadaFin:
            layout.addWidget(widget,5,3)
        for widget in self.ordenadaIntervalo:
            layout.addWidget(widget,5,4)

        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok|QtGui.QDialogButtonBox.Cancel)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox,10,1,1,4)

        for i in Ejes:
            self.ejesTabla.addItem(i)

        self.ejesTabla.currentIndexChanged.connect(self.ejesTabla_currentIndexChanged)
        self.ejesTabla_currentIndexChanged(0)


    def ejesTabla_currentIndexChanged(self, indice):
        """Hace los cambios pertinentes en la gui cuando se cambian los ejes de la tabla 3D:
            Actualiza unidades mostradas en la entrada de datos de la tabla
            Actualiza los checkbox de isolineas habilitados (todos menos los de los ejes x e y)"""

        coord=self.ejesTabla.currentText().split(",")
        dict=configVariables(self.parent()).copy()
        for entry in coord:
            del dict[str(entry)]

        self.variableTabla.clear()
        for key, value in dict.iteritems():
            self.variableTabla.addItem(value)

        i, j=self.currentIndex()
        self.label_ejeX.setText(self.abscisaInicio[i].unidad.__title__)
        self.label_ejeY.setText(self.ordenadaInicio[j].unidad.__title__)
        for lista in (self.abscisaInicio, self.abscisaFin, self.abscisaIntervalo):
            for widget in lista:
                widget.setVisible(False)
            lista[i].setVisible(True)
        for lista in (self.ordenadaInicio, self.ordenadaFin, self.ordenadaIntervalo):
            for widget in lista:
                widget.setVisible(False)
            lista[j].setVisible(True)

    def currentIndex(self):
        if self.metodo=="freesteam":
            i=[0, 0, 0, 0, 1, 1][self.ejesTabla.currentIndex()]
            j=[0, 1, 2, 3, 2, 4][self.ejesTabla.currentIndex()]

        elif self.metodo=="iapws":
            i=[0, 1, 1, 2, 0, 1][self.ejesTabla.currentIndex()]
            j=[0, 1, 2, 2, 3, 3][self.ejesTabla.currentIndex()]

        else:
            i=[0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 5][self.ejesTabla.currentIndex()]
            j=[0, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 3, 4, 5, 3, 4, 5, 4, 5, 5][self.ejesTabla.currentIndex()]

        return i, j


class PlotMEoS(QtGui.QWidget):
    """Gráfico de datos de las ecuaciones multiparametro, reescribe las acciones de menu contextual."""
    def __init__(self, dim, toolbar=False, parent=None):
        super(PlotMEoS, self).__init__(parent)
        self.parent=parent
        self.dim=dim

        layout=QtGui.QVBoxLayout(self)
        self.plot=plot.matplotlib(dim)
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

    def grid(self, bool):
        self.plot.ax.grid(bool)
        self.plot.draw()

    def edit(self):
        dialog=EditPlot(self, self.parent)
        dialog.show()

    def editAxis(self):
        dialog=EditAxis(self.plot)
        dialog.exec_()

    def table(self, obj):
        """Crea una tabla con los datos del grafico"""
        title=QtGui.QApplication.translate("pychemqt", "Table from")+" "+obj.get_label()
        tabla=createTabla(self.parent, title, obj.fluids)
        self.parent.centralwidget.currentWidget().addSubWindow(tabla)
        tabla.show()

#        HHeader=[str(x) for x in self.plot.data["x"]]
#        VHeader=[str(y) for y in self.plot.data["y"]]
#        tabla = TablaMEoS(len(self.plot.data["x"]), horizontalHeader=HHeader, verticalHeaderLabels=VHeader, stretch=False, readOnly=True, parent=self.parent)
#        tabla.setMatrix(self.plot.data["z"])
#        prefix=QtGui.QApplication.translate("pychemqt", "Table from")+" "
#        tabla.setWindowTitle(prefix+obj.get_label())
#        self.parent.centralwidget.currentWidget().addSubWindow(tabla)
#        tabla.show()



class TablaMEoS(Tabla):
    """Tabla de datos de las ecuaciones multiparametro, reescribe las acciones de menu contextual."""

    def __init__(self, *args, **kwargs):
        super(TablaMEoS, self).__init__(*args, **kwargs)
        self.horizontalHeader().setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.horizontalHeader().customContextMenuRequested.connect(self.show_header_context_menu)
        self.format=[{"format": 1, "decimales": 4, "signo": False}]*args[0]
        self.data=[]
        self.parent=kwargs.get("parent", None)

    def show_header_context_menu(self, position):
        column = self.horizontalHeader().logicalIndexAt(position)
        dialog=NumericFactor(self.format[column])
        if dialog.exec_():
            self.format[column]=dialog.args()
            self.setStr()
            self.resizeColumnToContents(column)
        self.setRangeSelected(QtGui.QTableWidgetSelectionRange(0, column, self.rowCount()-1, column), True)

    def setMatrix(self, data):
        self.data=data
        self.setStr()

    def setStr(self):
        for fila, array in enumerate(self.data):
            if fila>=self.rowCount():
                self.addRow()
            for columna, data in enumerate(array):
                str=representacion(data, **self.format[columna])
                self.setValue(fila, columna, str)


    def contextMenuEvent(self, event):
        menu = QtGui.QMenu()
        actionCopy=createAction(QtGui.QApplication.translate("pychemqt", "&Copy"), slot=partial(self.copy, event), shortcut=QtGui.QKeySequence.Copy, icon=os.environ["pychemqt"]+"/images/button/editCopy", parent=self)
        export = createAction(QtGui.QApplication.translate("pychemqt", "E&xport to csv"), self.export, icon=os.environ["pychemqt"]+"/images/button/export", tip=QtGui.QApplication.translate("pychemqt", "Export table to csv"), parent=self)
        menu.addAction(actionCopy)
        menu.addSeparator()
        menu.addAction(export)
        menu.exec_(event.globalPos())

    def copy(self, event):
        widget=self.itemAt(self.viewport().mapFromGlobal(event.globalPos()))
        QtGui.QApplication.clipboard().setText(widget.text())

    def export(self):
        dir = self.parent.currentFilename if self.parent.currentFilename else "."
        fname = unicode(QtGui.QFileDialog.getSaveFileName(self,
                            QtGui.QApplication.translate("pychemqt", "Export table to csv"), dir,
                            "csv files (*.csv)"))
        if fname:
            if fname.split(".")[-1]!="csv":
                fname+=".csv"

            cambio=maketrans(".", ",")
            with open(fname, "w") as archivo:
                writer = csv.writer(archivo, delimiter='\t', quotechar='"', quoting=csv.QUOTE_NONE)
                for row in self.data:
                    writer.writerow([str(data).translate(cambio) for data in row])







if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)

    SteamTables = Ui_Saturation()
#    SteamTables=AddLine(None)
#    SteamTables=transportDialog(mEoS.__all__[2])

    SteamTables.show()
    sys.exit(app.exec_())

