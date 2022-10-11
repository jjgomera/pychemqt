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
# Module with petroleum fraction pseudocomponent definition UI interfaces
#   -Definicion_Petro: Definition of crude and oil fraction
#   -View_Petro: Dialog to show the properties of a petroleum fraction
###############################################################################


import os
from functools import partial

from qt import QtGui, QtWidgets
from scipy import arange


from lib.config import IMAGE_PATH, Preferences
from lib.crude import Crudo
from lib import sql
from lib.plot import Plot
from lib.petro import Petroleo, curve_Predicted, _Tb_Predicted
from lib import unidades
from lib.unidades import Temperature, Pressure, Diffusivity
from UI import prefPetro
from UI.inputTable import InputTableWidget
from UI.newComponent import newComponent
from UI.widgets import Entrada_con_unidades


class View_Petro(QtWidgets.QDialog):
    """Dialog to show the properties of a petroleum fractions"""

    def __init__(self, petroleo=None, parent=None):
        super(View_Petro, self).__init__(parent)
        self.setWindowTitle(QtWidgets.QApplication.translate(
            "pychemqt", "Petrol assay characteristics"))
        layout = QtWidgets.QGridLayout(self)
        self.nombre = QtWidgets.QLabel()
        layout.addWidget(self.nombre, 1, 1, 1, 5)
        label = QtWidgets.QLabel("M")
        label.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Molecular Weight"))
        layout.addWidget(label, 2, 1)
        self.M = Entrada_con_unidades(float, textounidad="g/mol")
        layout.addWidget(self.M, 2, 2)
        label = QtWidgets.QLabel("Tb")
        label.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Boiling Temperature"))
        layout.addWidget(label, 3, 1)
        self.Tb = Entrada_con_unidades(unidades.Temperature)
        layout.addWidget(self.Tb, 3, 2)
        label = QtWidgets.QLabel("SG")
        label.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Specific gravity at 60ºF"))
        layout.addWidget(label, 4, 1)
        self.gravity = Entrada_con_unidades(float)
        layout.addWidget(self.gravity, 4, 2)
        label = QtWidgets.QLabel("API")
        label.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "API Specific gravity"))
        layout.addWidget(label, 5, 1)
        self.API = Entrada_con_unidades(float)
        layout.addWidget(self.API, 5, 2)
        label = QtWidgets.QLabel("Kw")
        label.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Watson characterization factor"))
        layout.addWidget(label, 6, 1)
        self.watson = Entrada_con_unidades(float)
        layout.addWidget(self.watson, 6, 2)
        label = QtWidgets.QLabel("n")
        label.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Refractive Index"))
        layout.addWidget(label, 7, 1)
        self.n = Entrada_con_unidades(float)
        layout.addWidget(self.n, 7, 2)
        label = QtWidgets.QLabel("I")
        label.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Huang parameter"))
        layout.addWidget(label, 8, 1)
        self.I = Entrada_con_unidades(float)
        layout.addWidget(self.I, 8, 2)
        label = QtWidgets.QLabel("ν<sub>100F</sub>")
        label.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Kinematic viscosity at 100ºF"))
        layout.addWidget(label, 9, 1)
        self.v100 = Entrada_con_unidades(unidades.Diffusivity)
        layout.addWidget(self.v100, 9, 2)
        label = QtWidgets.QLabel("ν<sub>210F</sub>")
        label.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Kinematic viscosity at 210ºF"))
        layout.addWidget(label, 10, 1)
        self.v210 = Entrada_con_unidades(unidades.Diffusivity)
        layout.addWidget(self.v210, 10, 2)

        layout.addWidget(QtWidgets.QLabel("Tc"), 2, 4)
        self.Tc = Entrada_con_unidades(unidades.Temperature)
        layout.addWidget(self.Tc, 2, 5)
        layout.addWidget(QtWidgets.QLabel("Pc"), 3, 4)
        self.Pc = Entrada_con_unidades(unidades.Pressure)
        layout.addWidget(self.Pc, 3, 5)
        layout.addWidget(QtWidgets.QLabel("Vc"), 4, 4)
        self.Vc = Entrada_con_unidades(unidades.SpecificVolume)
        layout.addWidget(self.Vc, 4, 5)
        layout.addWidget(QtWidgets.QLabel("Zc"), 5, 4)
        self.Zc = Entrada_con_unidades(float)
        layout.addWidget(self.Zc, 5, 5)
        label = QtWidgets.QLabel("ω")
        label.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Acentric factor"))
        layout.addWidget(label, 6, 4)
        self.f_acent = Entrada_con_unidades(float)
        layout.addWidget(self.f_acent, 6, 5)
        label = QtWidgets.QLabel("m")
        label.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Refractivity Intercept"))
        layout.addWidget(label, 7, 4)
        self.refractivity = Entrada_con_unidades(float)
        layout.addWidget(self.refractivity, 7, 5)
        layout.addWidget(QtWidgets.QLabel("CH"), 8, 4)
        self.CH = Entrada_con_unidades(float)
        layout.addWidget(self.CH, 8, 5)
        layout.addWidget(QtWidgets.QLabel("%S"), 9, 4)
        self.S = Entrada_con_unidades(float)
        layout.addWidget(self.S, 9, 5)
        layout.addWidget(QtWidgets.QLabel("%H"), 10, 4)
        self.H = Entrada_con_unidades(float)
        layout.addWidget(self.H, 10, 5)

        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "VGC")), 2, 7)
        self.VGC = Entrada_con_unidades(float)
        layout.addWidget(self.VGC, 2, 8)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Cetane index")), 3, 7)
        self.cetane = Entrada_con_unidades(float)
        layout.addWidget(self.cetane, 3, 8)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Pour point")), 4, 7)
        self.pour = Entrada_con_unidades(unidades.Temperature)
        layout.addWidget(self.pour, 4, 8)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Aniline point")), 5, 7)
        self.aniline = Entrada_con_unidades(unidades.Temperature)
        layout.addWidget(self.aniline, 5, 8)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Freezing point")), 6, 7)
        self.freezing = Entrada_con_unidades(unidades.Temperature)
        layout.addWidget(self.freezing, 6, 8)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Cloud point")), 7, 7)
        self.cloud = Entrada_con_unidades(unidades.Temperature)
        layout.addWidget(self.cloud, 7, 8)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Smoke point")), 8, 7)
        self.smoke = Entrada_con_unidades(unidades.Length)
        layout.addWidget(self.smoke, 8, 8)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Flash point (open)")), 9, 7)
        self.flashOpen = Entrada_con_unidades(unidades.Temperature)
        layout.addWidget(self.flashOpen, 9, 8)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Flash point (closed)")), 10, 7)
        self.flashClosed = Entrada_con_unidades(unidades.Temperature)
        layout.addWidget(self.flashClosed, 10, 8)

        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 15, 8)
        button = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Close)
        button.rejected.connect(self.reject)
        layout.addWidget(button, 16, 1, 1, 8)

        self.setReadOnly(True)
        if petroleo:
            self.rellenar(petroleo)

    def setReadOnly(self, bool):
        self.M.setReadOnly(bool)
        self.Tb.setReadOnly(bool)
        self.gravity.setReadOnly(bool)
        self.API.setReadOnly(bool)
        self.watson.setReadOnly(bool)

        self.Tc.setReadOnly(bool)
        self.Pc.setReadOnly(bool)
        self.Vc.setReadOnly(bool)
        self.Zc.setReadOnly(bool)
        self.f_acent.setReadOnly(bool)
        self.refractivity.setReadOnly(bool)
        self.CH.setReadOnly(bool)
        self.S.setReadOnly(bool)
        self.H.setReadOnly(bool)

        self.n.setReadOnly(bool)
        self.I.setReadOnly(bool)
        self.cetane.setReadOnly(bool)
        self.aniline.setReadOnly(bool)
        self.cloud.setReadOnly(bool)
        self.pour.setReadOnly(bool)
        self.freezing.setReadOnly(bool)
        self.smoke.setReadOnly(bool)
        self.v100.setReadOnly(bool)
        self.v210.setReadOnly(bool)
        self.VGC.setReadOnly(bool)
        self.flashOpen.setReadOnly(bool)
        self.flashClosed.setReadOnly(bool)

    def rellenar(self, petroleo):
        self.nombre.setText(petroleo.name)
        self.M.setValue(petroleo.M)
        self.Tb.setValue(petroleo.Tb)
        self.gravity.setValue(petroleo.SG)
        self.API.setValue(petroleo.API)
        self.watson.setValue(petroleo.Kw)

        self.Tc.setValue(petroleo.Tc)
        self.Pc.setValue(petroleo.Pc)
        self.Vc.setValue(petroleo.Vc)
        self.Zc.setValue(petroleo.Zc)
        self.f_acent.setValue(petroleo.f_acent)
        self.refractivity.setValue(petroleo.Ri)
        self.CH.setValue(petroleo.CH)
        self.S.setValue(petroleo.S)
        self.H.setValue(petroleo.H)

        self.n.setValue(petroleo.n)
        self.I.setValue(petroleo.I)
        self.cetane.setValue(petroleo.CetaneI)
        self.aniline.setValue(petroleo.AnilineP)
        self.cloud.setValue(petroleo.CloudP)
        self.pour.setValue(petroleo.PourP)
        self.freezing.setValue(petroleo.FreezingP)
        self.smoke.setValue(petroleo.SmokeP)
        self.v100.setValue(petroleo.v100)
        self.v210.setValue(petroleo.v210)
        # self.VGC.setValue(petroleo.VGC)
        if petroleo.hasCurve:
            self.flashOpen.setValue(petroleo.self.FlashPo)
            self.flashClosed.setValue(petroleo.self.FlashPc)


class Definicion_Petro(newComponent):
    """Dialog for define hypothetical crude and oil fraction"""
    ViewDetails = View_Petro

    def __init__(self, parent=None):
        super(Definicion_Petro, self).__init__(parent)
        self.setWindowTitle(QtWidgets.QApplication.translate(
            "pychemqt", "Petrol component definition"))

        layout = QtWidgets.QVBoxLayout(self)
        self.toolBox = QtWidgets.QTabWidget()
        self.toolBox.setTabPosition(QtWidgets.QTabWidget.TabPosition.South)
        layout.addWidget(self.toolBox)

        # Distillation data definition
        distilationPage = QtWidgets.QWidget()
        self.toolBox.addTab(
            distilationPage,
            QtWidgets.QApplication.translate("pychemqt", "Distillation data"))
        lyt = QtWidgets.QGridLayout(distilationPage)

        # Widget with curve functionality
        curveWidget = QtWidgets.QWidget()
        lytcurve = QtWidgets.QGridLayout(curveWidget)
        lytcurve.addWidget(QtWidgets.QLabel("Curve type"), 1, 1)
        self.tipoCurva = QtWidgets.QComboBox()
        for method in Petroleo.CURVE_TYPE:
            self.tipoCurva.addItem(method)
        self.tipoCurva.currentIndexChanged.connect(self.curveIndexChanged)
        lytcurve.addWidget(self.tipoCurva, 1, 2)
        self.curvaDestilacion = InputTableWidget(2)
        self.curvaDestilacion.tabla.horizontalHeader().show()
        self.curvaDestilacion.tabla.rowFinished.connect(self.checkStatusCurve)
        lytcurve.addWidget(self.curvaDestilacion, 2, 1, 3, 3)
        self.regresionButton = QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(
            os.path.join(IMAGE_PATH, "button", "fit.png"))),
            QtWidgets.QApplication.translate("pychemqt", "Regression"))
        self.regresionButton.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Calculate missing required values from a curve fit"))
        self.regresionButton.clicked.connect(self.regresionCurve)
        lytcurve.addWidget(self.regresionButton, 2, 3)
        self.finishButton = QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(
            os.path.join(IMAGE_PATH, "button", "arrow-right.png"))),
            QtWidgets.QApplication.translate("pychemqt", "Finish"))
        self.finishButton.clicked.connect(self.finishCurva)
        lytcurve.addWidget(self.finishButton, 5, 3)
        lytcurve.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Pressure")), 5, 1)
        self.presion = Entrada_con_unidades(Pressure, value=101325.)
        self.presion.valueChanged.connect(partial(
            self.changeParams, "P_curve"))
        lytcurve.addWidget(self.presion, 5, 2)
        lytcurve.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 6, 4)

        # Widget with crude functionality
        crudeWidget = QtWidgets.QWidget()
        lytcrude = QtWidgets.QGridLayout(crudeWidget)
        self.crudo = QtWidgets.QComboBox()
        self.crudo.addItem("")
        query = "SELECT name, location, API, sulfur FROM CrudeOil"
        sql.databank.execute(query)
        for name, location, API, sulfur in sql.databank:
            self.crudo.addItem("%s (%s)  API: %s %s: %s" % (
                name, location, API, "%S", sulfur))
        self.crudo.currentIndexChanged.connect(partial(
            self.changeParams, "index"))
        lytcrude.addWidget(self.crudo, 1, 1, 1, 2)
        lytcrude.addWidget(QtWidgets.QLabel("Pseudo C+"), 2, 1)
        self.Cplus = Entrada_con_unidades(int, width=50)
        self.Cplus.valueChanged.connect(partial(self.changeParams, "Cplus"))
        lytcrude.addWidget(self.Cplus, 2, 2)

        self.checkCurva = QtWidgets.QRadioButton(
            QtWidgets.QApplication.translate(
                "pychemqt", "Define destillation curve"))
        self.checkCurva.toggled.connect(curveWidget.setEnabled)
        curveWidget.setEnabled(False)
        lyt.addWidget(self.checkCurva, 1, 1, 1, 2)
        lyt.addWidget(curveWidget, 2, 1, 1, 2)
        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed, QtWidgets.QSizePolicy.Policy.Fixed),
            3, 1)
        self.checkCrudo = QtWidgets.QRadioButton(
            QtWidgets.QApplication.translate(
                "pychemqt", "Use petrol fraction from list"))
        self.checkCrudo.toggled.connect(self.changeUnknown)
        self.checkCrudo.toggled.connect(crudeWidget.setEnabled)
        crudeWidget.setEnabled(False)
        lyt.addWidget(self.checkCrudo, 4, 1, 1, 2)
        lyt.addWidget(crudeWidget, 5, 1, 1, 2)
        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed, QtWidgets.QSizePolicy.Policy.Fixed),
            6, 1, 1, 2)
        self.checkBlend = QtWidgets.QCheckBox(QtWidgets.QApplication.translate(
            "pychemqt", "Blend if its necessary"))
        lyt.addWidget(self.checkBlend, 7, 1, 1, 2)
        self.cutButton = QtWidgets.QPushButton(
            QtWidgets.QApplication.translate("pychemqt", "Define cut ranges"))
        self.cutButton.setEnabled(False)
        self.cutButton.clicked.connect(self.showCutRange)
        lyt.addWidget(self.cutButton, 7, 2)
        self.checkBlend.toggled.connect(self.cutButton.setEnabled)
        lyt.addItem(QtWidgets.QSpacerItem(
            5, 5, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 8, 1)

        # Definition with bulk properties
        definitionPage = QtWidgets.QWidget()
        self.toolBox.addTab(
            definitionPage,
            QtWidgets.QApplication.translate("pychemqt", "Bulk Definition"))

        lyt = QtWidgets.QGridLayout(definitionPage)
        txt = QtWidgets.QLabel("Tb")
        txt.setToolTip(
            QtWidgets.QApplication.translate("pychemqt", "Boiling point"))
        lyt.addWidget(txt, 1, 1)
        self.Tb = Entrada_con_unidades(Temperature)
        self.Tb.valueChanged.connect(partial(self.changeParams, "Tb"))
        lyt.addWidget(self.Tb, 1, 2)
        txt = QtWidgets.QLabel("M")
        txt.setToolTip(
            QtWidgets.QApplication.translate("pychemqt", "Molecular weight"))
        lyt.addWidget(txt, 2, 1)
        self.M = Entrada_con_unidades(float, textounidad="g/mol")
        self.M.valueChanged.connect(partial(self.changeParams, "M"))
        lyt.addWidget(self.M, 2, 2)
        txt = QtWidgets.QLabel("SG")
        txt.setToolTip(
            QtWidgets.QApplication.translate("pychemqt", "Specific Gravity"))
        lyt.addWidget(txt, 3, 1)
        self.SG = Entrada_con_unidades(float)
        self.SG.valueChanged.connect(partial(self.changeParams, "SG"))
        lyt.addWidget(self.SG, 3, 2)
        txt = QtWidgets.QLabel("API")
        txt.setToolTip(
            QtWidgets.QApplication.translate("pychemqt", "API Gravity"))
        lyt.addWidget(txt, 4, 1)
        self.API = Entrada_con_unidades(float)
        self.API.valueChanged.connect(partial(self.changeParams, "API"))
        lyt.addWidget(self.API, 4, 2)
        txt = QtWidgets.QLabel("Kw")
        txt.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Watson characterization factor"))
        lyt.addWidget(txt, 5, 1)
        self.Kw = Entrada_con_unidades(float)
        self.Kw.valueChanged.connect(partial(self.changeParams, "Kw"))
        lyt.addWidget(self.Kw, 5, 2)
        lyt.addWidget(QtWidgets.QLabel("C/H"), 6, 1)
        self.CH = Entrada_con_unidades(float)
        self.CH.valueChanged.connect(partial(self.changeParams, "CH"))
        lyt.addWidget(self.CH, 6, 2)
        txt = QtWidgets.QLabel("ν<sub>100F</sub>")
        txt.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Kinematic viscosity at 100ºF"))
        lyt.addWidget(txt, 7, 1)
        self.v100 = Entrada_con_unidades(Diffusivity)
        self.v100.valueChanged.connect(partial(self.changeParams, "v100"))
        lyt.addWidget(self.v100, 7, 2)
        txt = QtWidgets.QLabel("ν<sub>210F</sub>")
        txt.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Kinematic viscosity at 210ºF"))
        lyt.addWidget(txt, 8, 1)
        self.v210 = Entrada_con_unidades(Diffusivity)
        self.v210.valueChanged.connect(partial(self.changeParams, "v210"))
        lyt.addWidget(self.v210, 8, 2)
        txt = QtWidgets.QLabel("n")
        txt.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Refractive index"))
        lyt.addWidget(txt, 9, 1)
        self.n = Entrada_con_unidades(float)
        self.n.valueChanged.connect(partial(self.changeParams, "n"))
        lyt.addWidget(self.n, 9, 2)
        txt = QtWidgets.QLabel("I")
        txt.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Huang Parameter"))
        lyt.addWidget(txt, 10, 1)
        self.I = Entrada_con_unidades(float)
        self.I.valueChanged.connect(partial(self.changeParams, "I"))
        lyt.addWidget(self.I, 10, 2)
        lyt.addWidget(QtWidgets.QLabel("%S"), 11, 1)
        self.S = Entrada_con_unidades(float, spinbox=True, step=1.0, max=100)
        self.S.valueChanged.connect(partial(self.changeParams, "S"))
        lyt.addWidget(self.S, 11, 2)
        lyt.addWidget(QtWidgets.QLabel("%H"), 12, 1)
        self.H = Entrada_con_unidades(float, spinbox=True, step=1.0, max=100)
        self.H.valueChanged.connect(partial(self.changeParams, "H"))
        lyt.addWidget(self.H, 12, 2)
        lyt.addWidget(QtWidgets.QLabel("%N"), 13, 1)
        self.N = Entrada_con_unidades(float, spinbox=True, step=1.0, max=100)
        self.N.valueChanged.connect(partial(self.changeParams, "N"))
        lyt.addWidget(self.N, 13, 2)

        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed, QtWidgets.QSizePolicy.Policy.Fixed),
            14, 1, 1, 2)
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Alternate definition, poor accuracy")), 15, 1, 1, 2)
        txt = QtWidgets.QLabel("Nc")
        txt.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Carbon number"))
        lyt.addWidget(txt, 16, 1)
        self.Nc = Entrada_con_unidades(int, width=50)
        self.Nc.valueChanged.connect(partial(self.changeParams, "Nc"))
        lyt.addWidget(self.Nc, 16, 2)

        # Configuration
        configPage = prefPetro.Widget(Preferences)
        self.toolBox.addTab(
            configPage,
            QtGui.QIcon(IMAGE_PATH + "button/configure.png"),
            QtWidgets.QApplication.translate("pychemqt", "Configuration"))

        # Initialization section
        newComponent.loadUI(self)
        self.curveParameters = None  # Fitting parameter for distillation curve

        self.Petroleo = Petroleo()
        self.Crudo = Crudo()
        self.curveIndexChanged(0)
        self.checkStatusCurve()

    @property
    def unknown(self):
        if self.checkCrudo.isChecked():
            return self.Crudo
        else:
            return self.Petroleo

    def changeUnknown(self):
        self.status.setState(self.unknown.status, self.unknown.msg)
        self.buttonShowDetails.setEnabled(self.unknown.status)
        self.buttonBox.button(QtWidgets.QDialogButtonBox.StandardButton.Save).setEnabled(
            self.unknown.status)

    # Curve distillation definition
    def curveIndexChanged(self, index):
        """Show the composition unit appropiated to the new curve selected"""
        if index == 3:
            header = ["wt.%", "Tb, " + Temperature.text()]
        else:
            header = ["Vol.%", "Tb, " + Temperature.text()]
        self.curvaDestilacion.tabla.setHorizontalHeaderLabels(header)

    def finishCurva(self):
        """End the curve distillation definition and add the data to the
        Petroleo instance"""
        kwargs = {}
        curve = Petroleo.CURVE_TYPE[self.tipoCurva.currentIndex()]
        kwargs["curveType"] = curve
        kwargs["X_curve"] = self.curvaDestilacion.column(0)
        kwargs["T_curve"] = self.curvaDestilacion.column(1)
        kwargs["fit_curve"] = self.curveParameters
        self.calculo(**kwargs)

    def checkStatusCurve(self):
        """Check curren data of curve to check completeness of its definition
        and enable/disable accordly the buttons"""
        X = self.curvaDestilacion.column(0)
        self.regresionButton.setEnabled(len(X) > 3)

        defined = True
        for xi in [0.1, 0.5]:
            defined = defined and xi in X
        regresion = self.curveParameters is not None
        self.finishButton.setEnabled(defined or regresion)

    def regresionCurve(self):
        dlg = Plot(accept=True)
        x = self.curvaDestilacion.column(0)
        T = self.curvaDestilacion.column(1, Temperature)
        dlg.addData(x, T, color="black", ls="None", marker="s", mfc="red")
        parameters, r2 = curve_Predicted(x, T)
        xi = arange(0, 1, 0.01)
        Ti = [_Tb_Predicted(parameters, x_i) for x_i in xi]
        dlg.addData(xi, Ti, color="black", lw=0.5)

        # Add equation formula to plot
        txt = r"$\frac{T-T_{o}}{T_{o}}=\left[\frac{A}{B}\ln\left(\frac{1}{1-x}"
        txt += r"\right)\right]^{1/B}$"
        To = Temperature(parameters[0])
        txt2 = "\n\n\n$T_o=%s$" % To.str
        txt2 += "\n$A=%0.4f$" % parameters[1]
        txt2 += "\n$B=%0.4f$" % parameters[2]
        txt2 += "\n$r^2=%0.6f$" % r2
        dlg.plot.ax.text(0, T[-1], txt, size="14", va="top", ha="left")
        dlg.plot.ax.text(0, T[-1], txt2, size="10", va="top", ha="left")
        if dlg.exec():
            self.curveParameters = parameters
            self.checkStatusCurve()

    def showCutRange(self):
        pass


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    # petroleo = Petroleo(name="Petroleo", API=22.5, M=339.7)
    petroleo = Petroleo(name="Petroleo", Nc=20)
    Dialog = View_Petro(petroleo)
    Dialog.show()
    sys.exit(app.exec())
