#!/usr/bin/python
# -*- coding: utf-8 -*-

from functools import partial
import os

from PyQt5 import QtCore, QtGui, QtWidgets

from scipy import array, exp, optimize, linspace

from lib.plot import Plot
from lib.compuestos import Componente
from lib import unidades, sql
from .entrada_datos import Entrada_Datos, eqDIPPR
from UI.delegate import SpinEditor
from UI.widgets import Entrada_con_unidades, Tabla, okToContinue


class View_Petro(QtWidgets.QDialog):
    """Clase que define el ventana con las propiedades de fracciones petrolíferas"""
    def __init__(self, petroleo=None, parent=None):
        super(View_Petro, self).__init__(parent)
        self.setWindowTitle(QtCore.QCoreApplication.translate("pychemqt", "Petrol assay characteristics"))
        layout = QtWidgets.QGridLayout(self)
        self.nombre=QtWidgets.QLabel()
        layout.addWidget(self.nombre,1,1,1,5)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Molecular Weight")),2,1)
        self.M=Entrada_con_unidades(float, readOnly=True, textounidad="g/mol")
        layout.addWidget(self.M,2,2)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Boiling Point")),3,1)
        self.Tb=Entrada_con_unidades(unidades.Temperature, readOnly=True)
        layout.addWidget(self.Tb,3,2)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "SG 60ºF", None)),4,1)
        self.gravity=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.gravity,4,2)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "API gravity")),5,1)
        self.API=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.API,5,2)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Watson Factor")),6,1)
        self.watson=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.watson,6,2)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Refractive Index")),7,1)
        self.n=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.n,7,2)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Huang Parameter")),8,1)
        self.I=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.I,8,2)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "v100")),9,1)
        self.v100=Entrada_con_unidades(unidades.Diffusivity, readOnly=True)
        layout.addWidget(self.v100,9,2)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "v210")),10,1)
        self.v210=Entrada_con_unidades(unidades.Diffusivity, readOnly=True)
        layout.addWidget(self.v210,10,2)

        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Tc")),2,4)
        self.Tc=Entrada_con_unidades(unidades.Temperature, readOnly=True)
        layout.addWidget(self.Tc,2,5)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Pc")),3,4)
        self.Pc=Entrada_con_unidades(unidades.Pressure, readOnly=True)
        layout.addWidget(self.Pc,3,5)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Vc")),4,4)
        self.Vc=Entrada_con_unidades(unidades.SpecificVolume, readOnly=True)
        layout.addWidget(self.Vc,4,5)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Zc")),5,4)
        self.Zc=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.Zc,5,5)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Acentric factor")),6,4)
        self.f_acent=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.f_acent,6,5)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Refractivity Intercept")),7,4)
        self.refractivity=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.refractivity,7,5)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "C-H ratio")),8,4)
        self.CH=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.CH,8,5)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "% Sulfur")),9,4)
        self.S=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.S,9,5)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "%H")),10,4)
        self.H=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.H,10,5)

        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "VGC")),2,7)
        self.VGC=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.VGC,2,8)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Cetane index")),3,7)
        self.cetane=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.cetane,3,8)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Pour point")),4,7)
        self.pour=Entrada_con_unidades(unidades.Temperature, readOnly=True)
        layout.addWidget(self.pour,4,8)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Aniline point")),5,7)
        self.aniline=Entrada_con_unidades(unidades.Temperature, readOnly=True)
        layout.addWidget(self.aniline,5,8)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Freezing point")),6,7)
        self.freezing=Entrada_con_unidades(unidades.Temperature, readOnly=True)
        layout.addWidget(self.freezing,6,8)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Cloud point")),7,7)
        self.cloud=Entrada_con_unidades(unidades.Temperature, readOnly=True)
        layout.addWidget(self.cloud,7,8)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Smoke point")),8,7)
        self.smoke=Entrada_con_unidades(unidades.Length, readOnly=True)
        layout.addWidget(self.smoke,8,8)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Flash point (open)")),9,7)
        self.flashOpen=Entrada_con_unidades(unidades.Temperature, readOnly=True)
        layout.addWidget(self.flashOpen,9,8)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Flash point (closed)")),10,7)
        self.flashClosed=Entrada_con_unidades(unidades.Temperature, readOnly=True)
        layout.addWidget(self.flashClosed,10,8)
        
        layout.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Expanding),15,3,1,1)
        self.boton = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Close)
        self.boton.rejected.connect(self.accept)
        layout.addWidget(self.boton,16,1,1,8)
        
        if petroleo:
            self.rellenar(petroleo)
        
    def rellenar(self, petroleo):
        self.nombre.setText(petroleo.name)
        self.M.setValue(petroleo.M)
        self.Tb.setValue(petroleo.Tb)
        self.gravity.setValue(petroleo.SG)
        self.API.setValue(petroleo.API)
        self.watson.setValue(petroleo.watson)
        
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
        self.cetane.setValue(petroleo.CI)
        self.aniline.setValue(petroleo.AnilineP)
        self.cloud.setValue(petroleo.CloudP)
        self.pour.setValue(petroleo.PourP)
        self.freezing.setValue(petroleo.FreezingP)
        self.smoke.setValue(petroleo.SmokeP)
        self.v100.setValue(petroleo.v100)
        self.v210.setValue(petroleo.v210)
        self.VGC.setValue(petroleo.VGC)
        if petroleo.hasCurve:
            self.flashOpen.setValue(petroleo.self.FlashPo)
            self.flashClosed.setValue(petroleo.self.FlashPc)

class View_Contribution(QtWidgets.QDialog):
    """Ventana que muestra las propiedades del nuevo componente definido mediante los metodos de contribucion de grupo"""
    def __init__(self, petroleo=None, parent=None):
        super(View_Contribution, self).__init__(parent)
        self.setWindowTitle(QtCore.QCoreApplication.translate("pychemqt", "Group Contribution new component"))
        layout = QtWidgets.QGridLayout(self)

        self.nombre=QtWidgets.QLabel()
        layout.addWidget(self.nombre,1,1,1,5)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Molecular Weight")),2,1)
        self.M=Entrada_con_unidades(float, readOnly=True, textounidad="g/mol")
        layout.addWidget(self.M,2,2)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Boiling Point")),3,1)
        self.Tb=Entrada_con_unidades(unidades.Temperature, readOnly=True)
        layout.addWidget(self.Tb,3,2)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Melting Point")),4,1)
        self.Tf=Entrada_con_unidades(unidades.Temperature, readOnly=True)
        layout.addWidget(self.Tf,4,2)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Tc")),5,1)
        self.Tc=Entrada_con_unidades(unidades.Temperature, readOnly=True)
        layout.addWidget(self.Tc,5,2)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Pc")),6,1)
        self.Pc=Entrada_con_unidades(unidades.Pressure, readOnly=True)
        layout.addWidget(self.Pc,6,2)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Vc")),7,1)
        self.Vc=Entrada_con_unidades(unidades.SpecificVolume, readOnly=True)
        layout.addWidget(self.Vc,7,2)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "ΔHf", None)),8,1)
        self.Hf=Entrada_con_unidades(unidades.Enthalpy, readOnly=True)
        layout.addWidget(self.Hf,8,2)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "ΔGf", None)),9,1)
        self.Gf=Entrada_con_unidades(unidades.Enthalpy, readOnly=True)
        layout.addWidget(self.Gf,9,2)

        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "ΔHm", None)),2,4)
        self.Hm=Entrada_con_unidades(unidades.Enthalpy, readOnly=True)
        layout.addWidget(self.Hm,2,5)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "ΔHv", None)),3,4)
        self.Hv=Entrada_con_unidades(unidades.Enthalpy, readOnly=True)
        layout.addWidget(self.Hv,3,5)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Cpa")),4,4)
        self.cpa=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.cpa,4,5)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Cpb")),5,4)
        self.cpb=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.cpb,5,5)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Cpc")),6,4)
        self.cpc=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.cpc,6,5)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Cpd")),7,4)
        self.cpd=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.cpd,7,5)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "μa", None)),8,4)
        self.mua=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.mua,8,5)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "μb", None)),9,4)
        self.mub=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.mub,9,5)

        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "SG 60ºF", None)),2,7)
        self.gravity=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.gravity,2,8)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "API gravity")),3,7)
        self.API=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.API,3,8)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Watson Factor")),4,7)
        self.watson=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.watson,4,8)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Acentric factor")),5,7)
        self.f_acent=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.f_acent,5,8)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Rackett")),6,7)
        self.rackett=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.rackett,6,8)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Volume Liquid Constant")),7,7)
        self.Vliq=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.Vliq,7,8)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Solubility Parameter")),8,7)
        self.Parametro_solubilidad=Entrada_con_unidades(unidades.SolubilityParameter, readOnly=True)
        layout.addWidget(self.Parametro_solubilidad,8,8)

        layout.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Expanding),15,3,1,1)
        self.boton = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Close)
        self.boton.rejected.connect(self.accept)
        layout.addWidget(self.boton,16,1,1,8)
        
        if petroleo:
            self.rellenar(petroleo)
        
    def rellenar(self, petroleo):
        self.nombre.setText(petroleo.name)
        self.M.setValue(petroleo.M)
        self.Tb.setValue(petroleo.Tb)
        self.Tf.setValue(petroleo.Tf)
        self.Tc.setValue(petroleo.Tc)
        self.Pc.setValue(petroleo.Pc)
        self.Vc.setValue(petroleo.Vc)
        self.Hf.setValue(petroleo.Hf)
        self.Gf.setValue(petroleo.Gf)
        
        self.Hm.setValue(petroleo.Hm)
        self.Hv.setValue(petroleo.Hv)
        self.cpa.setValue(petroleo.cp[0])
        self.cpb.setValue(petroleo.cp[1])
        self.cpc.setValue(petroleo.cp[2])
        self.cpd.setValue(petroleo.cp[3])
        
        self.gravity.setValue(petroleo.SG)
        self.API.setValue(petroleo.API)
        self.watson.setValue(petroleo.watson)
        self.f_acent.setValue(petroleo.f_acent)
        self.rackett.setValue(petroleo.rackett)
        self.Vliq.setValue(petroleo.Vliq)
        self.Parametro_solubilidad.setValue(petroleo.Parametro_solubilidad)


class DIPPR_widget(QtWidgets.QGroupBox):
    """Clase que define los widgets comunes"""
    valueChanged = QtCore.pyqtSignal()
    def __init__(self, title, indice=0, propiedad="", parent=None):
        super(DIPPR_widget, self).__init__(title, parent)
        self.propiedad=propiedad
        self.title, self.unit=title.split(", ")
        self.parent=parent
        self.t=[]
        self.data=[]
        
        layout = QtWidgets.QGridLayout(self)
        self.buttonRegression = QtWidgets.QToolButton()
        self.buttonRegression.setToolTip(QtCore.QCoreApplication.translate("pychemqt", "Fit parameters from experimental data"))
        self.buttonRegression.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/Regression.gif")))
        self.buttonRegression.setIconSize(QtCore.QSize(32,32))
        self.buttonRegression.setFixedSize(QtCore.QSize(32,32))
        self.buttonRegression.clicked.connect(self.regresion)
        layout.addWidget(self.buttonRegression,1,1,2,1)
        self.buttonPlot = QtWidgets.QToolButton()
        self.buttonPlot.setToolTip(QtCore.QCoreApplication.translate("pychemqt", "Plot equation vs temperature"))
        self.buttonPlot.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/plot.png")))
        self.buttonPlot.setIconSize(QtCore.QSize(32,32))
        self.buttonPlot.setFixedSize(QtCore.QSize(32,32))
        self.buttonPlot.setEnabled(False)
        self.buttonPlot.clicked.connect(self.plot)
        layout.addWidget(self.buttonPlot,1,2,2,1)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "T min")),4,1)
        self.tmin =Entrada_con_unidades(unidades.Temperature, width=70)
        self.tmin.valueChanged.connect(self.valueChanged.emit)
        layout.addWidget(self.tmin,4,2)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "T max")),5,1)
        self.tmax =Entrada_con_unidades(unidades.Temperature, width=70)
        self.tmax.valueChanged.connect(self.valueChanged.emit)
        layout.addWidget(self.tmax,5,2)
        
        txt=["A", "B", "C", "D", "E"]
        self.entradas=[]
        for i in range(5):
            layout.addWidget(QtWidgets.QLabel(txt[i]),1+i,4)
            self.entradas.append(Entrada_con_unidades(float))
            self.entradas[-1].valueChanged.connect(self.valueChanged.emit)
            layout.addWidget(self.entradas[-1],1+i,5)
        layout.addItem(QtWidgets.QSpacerItem(0,0,QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Expanding),6,3)

        self.eq=eqDIPPR(1)
        layout.addWidget(self.eq,0,0,1,6)
        self.eqformula = QtWidgets.QLabel()
        self.eqformula.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.eqformula.setFrameShadow(QtWidgets.QFrame.Plain)
        self.eqformula.setScaledContents(False)
        self.eqformula.setAlignment(QtCore.Qt.AlignCenter)
        layout.addWidget(self.eqformula, 0, 0, 1, 6)

        self.changeIndice(indice)
        
        
    def changeIndice(self, indice):
        self.indice=indice
        if indice == 0 or indice > 1000:
            self.eqformula.setVisible(False)
            self.setReadOnly(False)
            self.buttonRegression.setEnabled(True)
        else:
            self.setReadOnly(True)
            self.eq.setVisible(True)
            self.buttonRegression.setEnabled(False)

    def setReadOnly(self, bool):
        for entrada in self.entradas:
            entrada.setReadOnly(bool)
        self.tmin.setReadOnly(bool)
        self.tmax.setReadOnly(bool)
        self.buttonRegression.setDisabled(bool)
        
    @property
    def value(self):
        valor=[self.eq.value()]
        for entrada in self.entradas:
            elemento=entrada.value
            if elemento:
                valor.append(elemento)
            else:
                valor.append(0)
        if valor[1:]==[0]*5:
            return []
        valor.append(self.tmin.value)
        valor.append(self.tmax.value)
        return valor

    def clear(self):
        self.tmin.clear()
        self.tmax.clear()
        self.eq.clear()
        self.eqformula.clear()
        for entrada in self.entradas:
            entrada.clear()
        self.t=[]
        self.data=[]
        self.buttonPlot.setEnabled(False)

    def rellenar(self, array, indice=None):
        if indice:
            self.changeIndice(indice)
        if array[0]!=0:
            for valor, entrada in zip(array[1:6], self.entradas):
                entrada.setValue(valor)
            self.tmin.setValue(array[6])
            self.tmax.setValue(array[7])
            self.eq.setValue(array[0])
            self.eqformula.setPixmap(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equation/DIPPR%i.gif" %array[0]))
            self.eq.setVisible(False)
            self.eqformula.setVisible(True)
            self.buttonPlot.setEnabled(True)
            self.equation=array[0]


    def formula_DIPPR(self, ecuacion, parametros):
        if ecuacion == 1:
            string="$%s = " % self.propiedad
            if parametros[0]:
                string+="%0.3f" % parametros[0]
            if parametros[1]:
                string+="%+0.3fT" % parametros[1]
            if parametros[2]:
                string+="%+0.3fT^2" % parametros[2]
            if parametros[3]:
                string+="%+0.3fT^3" % parametros[3]
            if parametros[4]:
                string+="%+0.3fT^4" % parametros[4]
            string+="$"
            return string
        elif ecuacion == 2:
            string="$%s = e^{" % self.propiedad
            if parametros[0]:
                string+="%0.3f" % parametros[0]
            if parametros[1]:
                string+="%+0.3f/T" % parametros[1]
            if parametros[2]:
                string+="%+0.3f \ln(T)" % parametros[2]
            if parametros[3]:
                string+="%+0.3fT^{%0.3f}" % (parametros[3], parametros[4])
            string+="}$"
            return string
        elif ecuacion == 3:
            string="$%s = %0.3f T^{\\frac{%0.3f}{1%+0.3f T" % (self.propiedad, parametros[0], parametros[1], parametros[2])
            if parametros[3]:
                string+="%+0.3f^2" % parametros[3]
            string+="}}$"
            return string
        elif ecuacion == 4:
            string="$%s = %0.3f%+0.3f*exp(" % (self.propiedad, parametros[0], parametros[1])
            if parametros[2] <0:
                string+="-"
            string+="%0.3f/T^{%0.3f})$" % (parametros[2], parametros[3])
            return string
        elif ecuacion == 5:
            string="$%s = " % self.propiedad
            if parametros[0]:
                string+="%0.3f" % parametros[0]
            if parametros[1]:
                string+="+\\frac{%0.3f}{T}" % parametros[1]
            if parametros[2]:
                string+="+\\frac{%0.3f}{T^3}" % parametros[2]
            if parametros[3]:
                string+="+\\frac{%0.3f}{T^8}" % parametros[3]
            if parametros[4]:
                string+="+\\frac{%0.3f}{T^9}" % parametros[4]
            string+="$"
            return string
        elif ecuacion == 6:
            string="$%s = \\frac{%0.3f}{%0.3f^{1+\\left(1-T/%0.3f\\right)^{%0.3f}}}$" % (self.propiedad, parametros[0], parametros[1], parametros[2], parametros[3])
            return string
        elif ecuacion == 7:
            string="$%s = %0.3f (1-T_r)^{" % (self.propiedad, parametros[0])
            if parametros[1]:
                string+="%0.3f" % parametros[1]
            if parametros[2]:
                string+="%+0.3fT_r" % parametros[2]
            if parametros[3]:
                string+="%+0.3fT_r^2" % parametros[3]
            if parametros[4]:
                string+="%+0.3fT_r^3" % parametros[4]
            string+="}$"
            return string
        elif ecuacion == 8:
            string="$%s = " % self.propiedad
            if parametros[0]:
                string+="%0.3f" % parametros[0]
            if parametros[1]:
                string+="%+0.3f\\left(\\frac{%0.3f/T}{sinh(%0.3f/T)}\\right)^2" % (parametros[1], parametros[2], parametros[2])
            if parametros[3]:
                string+="%+0.3f\\left(\\frac{%0.3f/T}{cosh(%0.3f/T)}\\right)^2" % (parametros[3], parametros[4], parametros[4])
            string+="$"    
            return string
        elif ecuacion == 9:
            string="$%s = " % self.propiedad
            if parametros[0]:
                string+="\\frac{%0.3f}{T_r}" % parametros[0]
            if parametros[1]:
                string+="%+0.3f" % parametros[1]
            if parametros[0] and parametros[2]:
                string+="%+0.3fT_r" % -2*parametros[0]*parametros[2]
            if parametros[0] and parametros[3]:
                string+="%+0.3fT_r^2" % -parametros[0]*parametros[3]
            if parametros[2]:
                string+="%+0.3fT_r^3" % -parametros[2]**2/3
            if parametros[2] and parametros[3]:
                string+="%+0.3fT_r^4" % -parametros[2]*parametros[3]/2
            if parametros[3]:
                string+="%+0.3fT_r^5" % -parametros[3]**2/5                
            string+="$"
            return string


    def plot(self):
        array=self.value
        t=linspace(array[-2], array[-1], 100)
        var=self.parent.compuesto.DIPPR(t, array[:-2])
        dialog=Plot()
        dialog.addData(t, var)
        if self.t and self.data:
            dialog.addData(self.t, self.data, "ro")
        dialog.plot.ax.grid(True)
        dialog.plot.ax.set_title(self.title, size="14")
        formula=self.formula_DIPPR(array[0], array[1:-2])
        dialog.addText(min(t), max(var), formula, size="14")
        dialog.plot.ax.set_xlabel("T, K", horizontalalignment='right', size="12")
        dialog.plot.ax.set_ylabel("$"+self.propiedad+",\;"+self.unit+"$", horizontalalignment='right', size="12")
        dialog.exec_()


    def regresion(self):
        dialogo=Entrada_Datos(title=self.title, DIPPR=True, horizontalHeader=["T, K", QtCore.QCoreApplication.translate("pychemqt", "Property")+", "+self.unit], tc=True, tcValue=self.parent.tempCritica.value, t=self.t, property=self.data, eq=self.eq.value())
        if dialogo.exec_():
            t=dialogo.tabla.getColumn(0, fill=False)
            p=dialogo.tabla.getColumn(1, fill=False)
            ecuacion=dialogo.eqDIPPR.value()
            self.parent.tempCritica.setValue(dialogo.tc.value)
#            self.parent.compuesto.Tc=dialogo.tc.value

            inicio=[1, 1, 1, 1, 1]
        
            def resto(parametros, ecuacion, t, f):
                var=array([ecuacion]+list(parametros))
                return f-array([self.parent.compuesto.DIPPR(ti, var) for ti in t])

            # Realizamos los ajustes
            ajuste=optimize.leastsq(resto,inicio,args=(ecuacion, t, p))
                
            dialog=Plot()
            dialog.addData(t, p, "ro")
            if ajuste[1] in range(1, 5):
                self.rellenar([ecuacion]+list(ajuste[0])+[min(t), max(t)])
                self.t=t
                self.data=p
                t2=linspace(min(t), max(t), 100)
                var=[self.parent.compuesto.DIPPR(ti, [ecuacion]+list(ajuste[0])) for ti in t2]
                dialog.addData(t2, var)
                dialog.plot.ax.grid(True)
                dialog.plot.ax.set_title(self.title, size="14")
                valores=self.value
                formula=self.formula_DIPPR(valores[0], valores[1:-2])
                dialog.addText(t[0]*0.95+t[-1]*0.05, max(var)*0.9+min(var)*0.1, formula, size="14")
                dialog.addText(t[0]*0.95+t[-1]*0.05, max(var)*0.85+min(var)*0.15, "residuo=%0.2e" % sum(resto(ajuste[0], ecuacion, t, p)), size="10")
                dialog.plot.ax.set_xlabel("T, K", horizontalalignment='right', size="12")
                dialog.plot.ax.set_ylabel("$"+self.propiedad+",\;"+self.unit+"$", horizontalalignment='right', size="12")
                self.valueChanged.emit()
            else:
                aviso=QtWidgets.QMessageBox.warning(self, QtCore.QCoreApplication.translate("pychemqt", "Warning"), QtCore.QCoreApplication.translate("pychemqt", "Fit unsuccessfully"))
            dialog.exec_()


class Parametric_widget(QtWidgets.QGroupBox):
    """Clase que define los widgets comunes"""
    valueChanged = QtCore.pyqtSignal()
    def __init__(self, title, indice=0, propiedad="", numero=2, imagen="", parent=None):
        super(Parametric_widget, self).__init__(title, parent)
        self.propiedad=propiedad
        self.title, self.unit=title.split(", ")
        self.indice=indice
        self.array=[]
        self.numero=numero
        self.imagen=imagen
        self.parent=parent
        self.t=[]
        self.data=[]

        layout = QtWidgets.QGridLayout(self)
        self.formula = QtWidgets.QLabel()
        self.formula.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.formula.setFrameShadow(QtWidgets.QFrame.Plain)
        self.formula.setAlignment(QtCore.Qt.AlignCenter)
        self.formula.setPixmap(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equation/%s.gif" % self.imagen))
        self.formula.setScaledContents(False)
        layout.addWidget(self.formula,0,1,1,5)
        self.buttonRegression = QtWidgets.QToolButton()
        self.buttonRegression.setToolTip(QtCore.QCoreApplication.translate("pychemqt", "Fit parameters from experimental data"))
        self.buttonRegression.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/Regression.gif")))
        self.buttonRegression.setIconSize(QtCore.QSize(32,32))
        self.buttonRegression.setFixedSize(QtCore.QSize(32,32))
        self.buttonRegression.clicked.connect(self.regresion)
        layout.addWidget(self.buttonRegression,1,1,2,1)
        self.buttonPlot = QtWidgets.QToolButton()
        self.buttonPlot.setToolTip(QtCore.QCoreApplication.translate("pychemqt", "Plot equation vs temperature"))
        self.buttonPlot.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/plot.png")))
        self.buttonPlot.setIconSize(QtCore.QSize(32,32))
        self.buttonPlot.setFixedSize(QtCore.QSize(32,32))
        self.buttonPlot.clicked.connect(self.plot)
        self.buttonPlot.setEnabled(False)
        layout.addWidget(self.buttonPlot,1,2,2,1)
        
        txt=["A", "B", "C", "D"]
        self.entradas=[]
        for i in range(numero):
            layout.addWidget(QtWidgets.QLabel(txt[i]),1+i,4)
            self.entradas.append(Entrada_con_unidades(float))
            self.entradas[-1].valueChanged.connect(self.valueChanged.emit)
            layout.addWidget(self.entradas[-1],1+i,5)
        layout.addItem(QtWidgets.QSpacerItem(0,0,QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Expanding),2+numero,3)

        self.changeIndice(indice)
        
    def changeIndice(self, indice):
        self.indice=indice
        if indice == 0 or indice > 1000:
            self.setReadOnly(False)
            self.buttonRegression.setEnabled(True)
        else:
            self.setReadOnly(True)
            self.buttonRegression.setEnabled(False)

    def setReadOnly(self, bool):
        for entrada in self.entradas:
            entrada.setReadOnly(bool)
        self.buttonRegression.setDisabled(bool)

    @property
    def value(self):
        valor=[]
        for entrada in self.entradas:
            elemento=entrada.value
            if elemento:
                valor.append(elemento)
            else:
                valor.append(0)
        if valor!=[0]*self.numero:
            return valor
        else:
            return []

    def clear(self):
        for entrada in self.entradas:
            entrada.clear()
        self.t=[]
        self.data=[]
        self.buttonPlot.setEnabled(False)

    
    def rellenar(self, array, indice=None):
        if indice:
            self.changeIndice(indice)

        if array!=(0, )*self.numero:
            self.array=array
            for i, valor in enumerate(array):
                self.entradas[i].setValue(valor)
            self.buttonPlot.setEnabled(True)


    def plot(self):
        array=self.value
        if self.parent.compuesto.Tf:
            tmin=self.parent.compuesto.Tf
        elif self.parent.compuesto.presion_vapor!=[0]*8:
            tmin=self.parent.compuesto.presion_vapor[-2]
        else:
            tmin=300
        if self.parent.compuesto.Tb:
            tmax=self.parent.compuesto.Tb
        elif self.parent.compuesto.presion_vapor!=[0]*8:
            tmax=self.parent.compuesto.presion_vapor[-1]
        else:
            tmax=500
        t=linspace(tmin, tmax, 100)
        
        funcion={"viscosidad": self.parent.compuesto.Mu_Liquido_Parametrica, 
                        "antoine": self.parent.compuesto.Pv_Antoine, 
                        "tensionsuperficial": self.parent.compuesto.Tension_Parametrica, 
                        "henry": self.parent.compuesto.constante_Henry}

        formula={"viscosidad": "$\\log\\mu=%0.2f\\left(\\frac{1}{T}-\\frac{1}{%0.2f}\\right)$" % self.parent.compuesto.viscosidad_parametrica, 
                        "antoine": "$P_{v}=e^{%0.2f-\\frac{%0.2f}{T%+0.2f}}$" % self.parent.compuesto.antoine, 
                        "tensionsuperficial": "$\\sigma=%0.2f\\left(1-T_{r}\\right)^{%0.2f}$" % self.parent.compuesto.tension_superficial_parametrica, 
                        "henry": "$\\lnH=\\frac{%0.2f}{T}%+0.2f\\cdot\\lnT%+0.2f\\cdot T%+0.2f$" % self.parent.compuesto.henry}

        var=[funcion[self.imagen](ti) for ti in t]
        dialog=Plot()
        dialog.addData(t, var)
        if self.t and self.data:
            dialog.addData(self.t, self.data, "ro")
        dialog.plot.ax.grid(True)
        dialog.plot.ax.set_title(self.title, size="14")
        dialog.addText(min(t), max(var), formula[self.imagen], size="14")
        dialog.plot.ax.set_xlabel("T, K", horizontalalignment='right', size="12")
        dialog.plot.ax.set_ylabel("$"+self.propiedad+",\;"+self.unit+"$", horizontalalignment='right', size="12")
        dialog.exec_()


    def regresion(self):
        dialogo=Entrada_Datos(title=self.title, horizontalHeader=["T, K", QtCore.QCoreApplication.translate("pychemqt", "Property")+", "+self.unit], tc=self.imagen=="tensionsuperficial", tcValue=self.parent.tempCritica.value)
        if dialogo.exec_():
            t=dialogo.tabla.getColumn(0, fill=False)
            p=dialogo.tabla.getColumn(1, fill=False)
            self.parent.tempCritica.setValue(dialogo.tc.value)
            self.parent.compuesto.Tc=dialogo.tc.value
            
            #TODO: Añadir soporte para configurar valores iniciales
            inicio=[1.]*self.numero
            funcion={"viscosidad": self.parent.compuesto.Mu_Liquido_Parametrica, 
                            "antoine": self.parent.compuesto.Pv_Antoine, 
                            "tensionsuperficial": self.parent.compuesto.Tension_Parametrica, 
                            "henry": self.parent.compuesto.constante_Henry}

            def resto(parametros, t, f):
                return f-array([funcion[self.imagen](ti, parametros) for ti in t])

            # Realizamos los ajustes
            ajuste=optimize.leastsq(resto,inicio,args=(t, p))
                
            if ajuste[1] in range(1, 5):
                self.rellenar(list(ajuste[0]))
                self.t=t
                self.data=p
                dialog=Plot()
                dialog.addData(t, p, "ro")
                t2=linspace(t[0], t[-1], 100)
                var=[funcion[self.imagen](ti, ajuste[0]) for ti in t2]
                dialog.addData(t2, var)
                dialog.plot.ax.grid(True)
                dialog.plot.ax.set_title(self.title, size="14")

                if self.imagen=="viscosidad":
                    formula="$\\log\\mu=%0.2f\\left(\\frac{1}{T}-\\frac{1}{%0.2f}\\right)$" % tuple(ajuste[0])
                elif self.imagen=="antoine":
                    formula="$P_{v}=e^{%0.2f-\\frac{%0.2f}{T%+0.2f}}$" % tuple(ajuste[0])
                elif self.imagen=="tensionsuperficial":
                    formula="$\\sigma=%0.2f\\left(1-T_{r}\\right)^{%0.2f}$" % tuple(ajuste[0])
                elif self.imagen=="henry":
                    formula="$\\lnH=\\frac{%0.2f}{T}%+0.2f\\cdot\\lnT%+0.2f\\cdot T%+0.2f$" % tuple(ajuste[0])
                dialog.addText(t[0]*0.95+t[-1]*0.05, max(var)*0.9+min(var)*0.1, formula, size="14")
                dialog.addText(t[0]*0.95+t[-1]*0.05, max(var)*0.85+min(var)*0.15, "residuo=%0.2e" % sum(resto(ajuste[0], t, p)), size="10")
                dialog.plot.ax.set_xlabel("T, K", horizontalalignment='right', size="12")
                dialog.plot.ax.set_ylabel("$"+self.propiedad+",\;"+self.unit+"$", horizontalalignment='right', size="12")
                dialog.exec_()
                self.valueChanged.emit()
            else:
                aviso=QtWidgets.QMessageBox.warning(self, QtCore.QCoreApplication.translate("pychemqt", "Warning"), QtCore.QCoreApplication.translate("pychemqt", "Fit unsuccessfully"))


class View_Component(QtWidgets.QDialog):
    compuesto=Componente()
    
    def __init__(self, indice=0, parent=None):
        super(View_Component, self).__init__(parent)
        gridLayout = QtWidgets.QGridLayout(self)
        
        layoutTitle=QtWidgets.QHBoxLayout()
        layoutTitle.setSpacing(5)
        if indice == 0:
            self.setWindowTitle(QtCore.QCoreApplication.translate("pychemqt", "Define New Component"))
            layoutTitle.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Name")+": ", self))
            self.nombre = QtWidgets.QLineEdit(self)
            self.nombre.editingFinished.connect(self.setDirty)
            layoutTitle.addWidget(self.nombre)
            layoutTitle.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Fixed))
        else:
            if indice > 1000:
                self.setWindowTitle(QtCore.QCoreApplication.translate("pychemqt", "Custom Component Properties"))
            else:
                self.setWindowTitle(QtCore.QCoreApplication.translate("pychemqt", "Component Properties"))
            self.buttonFirst=QtWidgets.QToolButton(self)
            self.buttonFirst.setToolTip(QtCore.QCoreApplication.translate("pychemqt", "Go to first element"))
            self.buttonFirst.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/arrow-left-double.png")))
            self.buttonFirst.clicked.connect(partial(self.change, 1))
            layoutTitle.addWidget(self.buttonFirst)
            self.buttonPrevious=QtWidgets.QToolButton(self)
            self.buttonPrevious.setToolTip(QtCore.QCoreApplication.translate("pychemqt", "Go to previous element"))
            self.buttonPrevious.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/arrow-left.png")))
            self.buttonPrevious.clicked.connect(partial(self.change, "-1"))
            layoutTitle.addWidget(self.buttonPrevious)
            self.nombre = QtWidgets.QLabel(self)
            layoutTitle.addWidget(self.nombre)
            self.buttonNext=QtWidgets.QToolButton(self)
            self.buttonNext.setToolTip(QtCore.QCoreApplication.translate("pychemqt", "Go to next element"))
            self.buttonNext.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/arrow-right.png")))
            self.buttonNext.clicked.connect(partial(self.change, "+1"))
            layoutTitle.addWidget(self.buttonNext)
            self.buttonLast=QtWidgets.QToolButton(self)
            self.buttonLast.setToolTip(QtCore.QCoreApplication.translate("pychemqt", "Go to last element"))
            self.buttonLast.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/arrow-right-double.png")))
            self.buttonLast.clicked.connect(partial(self.change, "last"))
            layoutTitle.addWidget(self.buttonLast)
        gridLayout.addItem(layoutTitle,1,1,1,3)
        
        tabWidget = QtWidgets.QTabWidget()
        gridLayout.addWidget(tabWidget,2,1,1,3)
        
        self.buttonBox = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Discard|QtWidgets.QDialogButtonBox.Save|QtWidgets.QDialogButtonBox.Apply|QtWidgets.QDialogButtonBox.Close)
        self.buttonBox.clicked.connect(self.buttonClicked)
        gridLayout.addWidget(self.buttonBox,3,1,1,3)
        
        #Pestaña general
        tab1 = QtWidgets.QWidget()
        tabWidget.addTab(tab1,QtCore.QCoreApplication.translate("pychemqt", "&General"))
        gridLayout_General=QtWidgets.QGridLayout(tab1)
        gridLayout_General.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Formula")),1,1)
        self.formula1 = QtWidgets.QLineEdit()
        self.formula1.editingFinished.connect(self.setDirty)
        gridLayout_General.addWidget(self.formula1,1,2)
        gridLayout_General.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "CAS Number")),2,1)
        self.CAS_number = QtWidgets.QLineEdit()
        self.CAS_number.editingFinished.connect(self.setDirty)
        gridLayout_General.addWidget(self.CAS_number,2,2)
        gridLayout_General.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Alternative Name")),3,1)
        self.nombre_alternativo = QtWidgets.QLineEdit()
        self.nombre_alternativo.editingFinished.connect(self.setDirty)
        gridLayout_General.addWidget(self.nombre_alternativo,3,2)
        
        self.labelSmile=QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Smile Code"))
        gridLayout_General.addWidget(self.labelSmile,4,1,1,1)
        self.smile=QtWidgets.QLineEdit()
        self.smile.editingFinished.connect(self.setDirty)
        gridLayout_General.addWidget(self.smile,4,2,1,1)
        self.labelFormula2=QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Expanded Formula"))
        gridLayout_General.addWidget(self.labelFormula2,5,1,1,1)
        self.formula2 = QtWidgets.QLabel()
        self.formula2.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.formula2.setFrameShadow(QtWidgets.QFrame.Plain)
        self.formula2.setAlignment(QtCore.Qt.AlignCenter)
        self.formula2.setScaledContents(False)
        gridLayout_General.addWidget(self.formula2,5,2,4,1)

        label_UNIFAC=QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "UNIFAC Groups" ))
        label_UNIFAC.setAlignment(QtCore.Qt.AlignTop)
        gridLayout_General.addWidget(label_UNIFAC,9,1,1,1)
        
        self.UNIFAC=Tabla(2, filas=1, horizontalHeader=[QtCore.QCoreApplication.translate("pychemqt", "Group" ), QtCore.QCoreApplication.translate("pychemqt", "Contribution")], verticalHeader=False, dinamica=True, readOnly=indice)
        self.UNIFAC.setItemDelegateForColumn(0, SpinEditor(self))
        self.UNIFAC.setItemDelegateForColumn(1, SpinEditor(self))    
        self.UNIFAC.setColumnWidth(0, 80)
        self.UNIFAC.setFixedWidth(160)
        gridLayout_General.addWidget(self.UNIFAC, 9,2,3,1)

        gridLayout_General.addItem(QtWidgets.QSpacerItem(30,30,QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Expanding),9,3,3,1)
        
        gridLayout_General.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Molecular Weight")),1,4)
        self.M = Entrada_con_unidades(float,  textounidad="g/mol")
        self.M.valueChanged.connect(self.setDirty)
        gridLayout_General.addWidget(self.M,1,5)
        gridLayout_General.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "T critic ")),2,4)
        self.tempCritica = Entrada_con_unidades(unidades.Temperature)
        self.tempCritica.valueChanged.connect(self.setDirty)
        gridLayout_General.addWidget(self.tempCritica,2,5)
        gridLayout_General.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "P critic")),3,4)
        self.Pc = Entrada_con_unidades(unidades.Pressure)
        self.Pc.valueChanged.connect(self.setDirty)
        gridLayout_General.addWidget(self.Pc,3,5)
        gridLayout_General.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "V critic")),4,4)
        self.volumenCritico = Entrada_con_unidades(unidades.SpecificVolume)
        self.volumenCritico.valueChanged.connect(self.setDirty)
        gridLayout_General.addWidget(self.volumenCritico,4,5)
        gridLayout_General.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Acentric factor")),5,4)
        self.factorAcentrico = Entrada_con_unidades(float)
        self.factorAcentrico.valueChanged.connect(self.setDirty)
        gridLayout_General.addWidget(self.factorAcentrico,5,5)
        gridLayout_General.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "SG 60ºF", None)),6,4)
        self.SG = Entrada_con_unidades(float)
        self.SG.valueChanged.connect(self.setDirty)
        gridLayout_General.addWidget(self.SG,6,5)
        gridLayout_General.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Melting Point")),7,4)
        self.temp_fusion = Entrada_con_unidades(unidades.Temperature)
        self.temp_fusion.valueChanged.connect(self.setDirty)
        gridLayout_General.addWidget(self.temp_fusion,7,5)
        gridLayout_General.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Boiling Point")),8,4)
        self.temp_ebullicion = Entrada_con_unidades(unidades.Temperature)
        self.temp_ebullicion.valueChanged.connect(self.setDirty)
        gridLayout_General.addWidget(self.temp_ebullicion,8,5)
        
        gridLayout_General.addWidget(QtWidgets.QLabel("ΔH<sub>f</sub>"),9,4)
        self.calorFormacionGas = Entrada_con_unidades(unidades.Enthalpy)
        self.calorFormacionGas.valueChanged.connect(self.setDirty)
        gridLayout_General.addWidget(self.calorFormacionGas,9,5)
        gridLayout_General.addWidget(QtWidgets.QLabel("ΔG<sub>f</sub>"),10,4)
        self.energiaGibbsGas = Entrada_con_unidades(unidades.Enthalpy)
        self.energiaGibbsGas.valueChanged.connect(self.setDirty)
        gridLayout_General.addWidget(self.energiaGibbsGas,10,5)
        gridLayout_General.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Expanding),11,4,1,2)
        gridLayout_General.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Expanding),13,1,1,6)
        
        
        #pestaña cp
        tab2 = QtWidgets.QWidget()
        tabWidget.addTab(tab2,QtCore.QCoreApplication.translate("pychemqt", "&Cp"))
        gridLayout_cp = QtWidgets.QGridLayout(tab2)
        group_CpGasIdeal = QtWidgets.QGroupBox(QtCore.QCoreApplication.translate("pychemqt", "Cp ideal gas")+", cal/mol·K")
        gridLayout_cp.addWidget(group_CpGasIdeal,1,1,1,1)
        lyt_cpGasIdeal = QtWidgets.QGridLayout(group_CpGasIdeal)
        self.eqCpGasIdeal = QtWidgets.QLabel()
        self.eqCpGasIdeal.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.eqCpGasIdeal.setFrameShadow(QtWidgets.QFrame.Plain)
        self.eqCpGasIdeal.setPixmap(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equation/cp.gif"))
        self.eqCpGasIdeal.setScaledContents(False)
        lyt_cpGasIdeal.addWidget(self.eqCpGasIdeal,0,1,1,5)
        lyt_cpGasIdeal.addWidget(QtWidgets.QLabel("A"),1,1)
        self.cpa = Entrada_con_unidades(float)
        self.cpa.valueChanged.connect(self.setDirty)
        lyt_cpGasIdeal.addWidget(self.cpa,1,2)
        lyt_cpGasIdeal.addWidget(QtWidgets.QLabel("B"),2,1)
        self.cpb = Entrada_con_unidades(float)
        self.cpb.valueChanged.connect(self.setDirty)
        lyt_cpGasIdeal.addWidget(self.cpb,2,2)
        lyt_cpGasIdeal.addWidget(QtWidgets.QLabel("C"),3,1)
        self.cpc = Entrada_con_unidades(float)
        self.cpc.valueChanged.connect(self.setDirty)
        lyt_cpGasIdeal.addWidget(self.cpc,3,2)
        lyt_cpGasIdeal.addWidget(QtWidgets.QLabel("D"),1,4)
        self.cpd = Entrada_con_unidades(float)
        self.cpd.valueChanged.connect(self.setDirty)
        lyt_cpGasIdeal.addWidget(self.cpd,1,5)
        lyt_cpGasIdeal.addWidget(QtWidgets.QLabel("E"),2,4)
        self.cpe = Entrada_con_unidades(float)
        self.cpe.valueChanged.connect(self.setDirty)
        lyt_cpGasIdeal.addWidget(self.cpe,2,5)
        lyt_cpGasIdeal.addWidget(QtWidgets.QLabel("F"),3,4)
        self.cpf = Entrada_con_unidades(float)
        self.cpf.valueChanged.connect(self.setDirty)
        lyt_cpGasIdeal.addWidget(self.cpf,3,5)
        lyt_cpGasIdeal.addItem(QtWidgets.QSpacerItem(0,0,QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Expanding),4,1,1,4)
        
        self.cpgasDIPPR = DIPPR_widget(QtCore.QCoreApplication.translate("pychemqt", "Cp ideal gas DIPPR")+", cal/mol K", indice, "Cp_g", parent=self)
        self.cpgasDIPPR.valueChanged.connect(self.setDirty)
        gridLayout_cp.addWidget(self.cpgasDIPPR,1,2)
        self.cpLiquidoDIPPR = DIPPR_widget(QtCore.QCoreApplication.translate("pychemqt", "Cp liquid DIPPR")+", cal/mol K", indice, "Cp_l", parent=self)
        self.cpLiquidoDIPPR.valueChanged.connect(self.setDirty)
        gridLayout_cp.addWidget(self.cpLiquidoDIPPR,2,1)
        self.cpSolidoDIPPR = DIPPR_widget(QtCore.QCoreApplication.translate("pychemqt", "Cp solid DIPPR")+", cal/mol K", indice, "Cp_s", parent=self)
        self.cpSolidoDIPPR.valueChanged.connect(self.setDirty)
        gridLayout_cp.addWidget(self.cpSolidoDIPPR,2,2)
        gridLayout_cp.addItem(QtWidgets.QSpacerItem(0,0,QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Expanding),3,1)


        #pestaña densidad
        tab3 = QtWidgets.QWidget()
        tabWidget.addTab(tab3,QtCore.QCoreApplication.translate("pychemqt", "&Density"))
        gridLayout_rho = QtWidgets.QGridLayout(tab3)
        
        self.RhoSolido = DIPPR_widget(QtCore.QCoreApplication.translate("pychemqt", "Solid Density DIPPR")+", kmol/m^3", indice, "\\rho_s", parent=self)
        self.RhoSolido.valueChanged.connect(self.setDirty)
        gridLayout_rho.addWidget(self.RhoSolido,1,1)
        self.RhoLiquido = DIPPR_widget(QtCore.QCoreApplication.translate("pychemqt", "Liquid Density DIPPR")+", kmol/m^3", indice, "\\rho_l", parent=self)
        self.RhoLiquido.valueChanged.connect(self.setDirty)
        gridLayout_rho.addWidget(self.RhoLiquido,1,2)
        gridLayout_rho.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Expanding),2,1)


        #Pestaña viscosidad
        tab4 = QtWidgets.QWidget()
        tabWidget.addTab(tab4,QtCore.QCoreApplication.translate("pychemqt", "&Viscosity"))
        gridLayout_mu = QtWidgets.QGridLayout(tab4)

        self.muLiquido = DIPPR_widget(QtCore.QCoreApplication.translate("pychemqt", "Liquid Viscosity DIPPR")+", Pa s", indice, "\\mu_l", parent=self)
        self.muLiquido.valueChanged.connect(self.setDirty)
        gridLayout_mu.addWidget(self.muLiquido,1,1)
        self.muGas = DIPPR_widget(QtCore.QCoreApplication.translate("pychemqt", "Gas Viscosity DIPPR")+", Pa s", indice, "\\mu_g", parent=self)
        self.muGas.valueChanged.connect(self.setDirty)
        gridLayout_mu.addWidget(self.muGas,1,2)
        self.muParametric = Parametric_widget(QtCore.QCoreApplication.translate("pychemqt", "Viscosity")+", cP", indice, "\\mu_g", 2, "viscosidad", parent=self)
        self.muParametric.valueChanged.connect(self.setDirty)
        gridLayout_mu.addWidget(self.muParametric,2,1)
        gridLayout_mu.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Expanding),3,1)

        
        #Pestaña Pv, Hv
        tab5 = QtWidgets.QWidget()
        tabWidget.addTab(tab5,QtCore.QCoreApplication.translate("pychemqt", "P&v && Hv"))
        gridLayout_vapor = QtWidgets.QGridLayout(tab5)

        self.boilingHeat = DIPPR_widget(QtCore.QCoreApplication.translate("pychemqt", "Heat of vaporization DIPPR")+", J/kmol", indice, "H_v", parent=self)
        self.boilingHeat.valueChanged.connect(self.setDirty)
        gridLayout_vapor.addWidget(self.boilingHeat,1,1)
        self.vaporPressure = DIPPR_widget(QtCore.QCoreApplication.translate("pychemqt", "Vapor Pressure DIPPR")+", Pa", indice, "P_v", parent=self)
        self.vaporPressure.valueChanged.connect(self.setDirty)
        gridLayout_vapor.addWidget(self.vaporPressure,1,2)
        self.PvAntoine = Parametric_widget(QtCore.QCoreApplication.translate("pychemqt", "Antoine Vapor Pressure")+", mmHg", indice, "P_v", 3, "antoine", parent=self)
        self.PvAntoine.valueChanged.connect(self.setDirty)
        gridLayout_vapor.addWidget(self.PvAntoine,2,1)
        gridLayout_vapor.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Expanding),3,1)


        #Pestaña conductividad térmica y tension
        tab6 = QtWidgets.QWidget()
        tabWidget.addTab(tab6,QtCore.QCoreApplication.translate("pychemqt", "&Tension && Conductivity"))
        gridLayout_conductividad = QtWidgets.QGridLayout(tab6)

        self.conductivityLiquido = DIPPR_widget(QtCore.QCoreApplication.translate("pychemqt", "Liquid Thermal Conductivity DIPPR")+", W/m K", indice, "\\lambda_l", parent=self)
        self.conductivityLiquido.valueChanged.connect(self.setDirty)
        gridLayout_conductividad.addWidget(self.conductivityLiquido,1,1)
        self.conductivityGas = DIPPR_widget(QtCore.QCoreApplication.translate("pychemqt", "Gas Thermal Conductivity DIPPR")+", W/m K", indice, "\lambda_v", parent=self)
        self.conductivityGas.valueChanged.connect(self.setDirty)
        gridLayout_conductividad.addWidget(self.conductivityGas,1,2)
        self.tension = DIPPR_widget(QtCore.QCoreApplication.translate("pychemqt", "Surface Tension DIPPR")+", N/m", indice, "\sigma", parent=self)
        self.tension.valueChanged.connect(self.setDirty)
        gridLayout_conductividad.addWidget(self.tension,2,1)
        self.tensionParametrica = Parametric_widget(QtCore.QCoreApplication.translate("pychemqt", "Surface Tension")+", N/m", indice, "\sigma", 2, "tensionsuperficial", parent=self)
        self.tensionParametrica.valueChanged.connect(self.setDirty)
        gridLayout_conductividad.addWidget(self.tensionParametrica,2,2)
        gridLayout_conductividad.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Expanding),3,1)

        #Pestaña parametricos
        tab7 = QtWidgets.QWidget()
        tabWidget.addTab(tab7,QtCore.QCoreApplication.translate("pychemqt", "&Parametrics"))
        gridLayout_Parametricos=QtWidgets.QGridLayout(tab7)

        self.Henry = Parametric_widget(QtCore.QCoreApplication.translate("pychemqt", "Henry Costant")+", psia/Xg", indice, "K_H", 4, "henry", parent=self)
        self.Henry.valueChanged.connect(self.setDirty)
        gridLayout_Parametricos.addWidget(self.Henry,1,1,2,1)

        group_MSRK = QtWidgets.QGroupBox(QtCore.QCoreApplication.translate("pychemqt", "MSRK Coefficients"))
        gridLayout_Parametricos.addWidget(group_MSRK,1,2)
        lyt_MSRK = QtWidgets.QGridLayout(group_MSRK)
        lyt_MSRK.addWidget(QtWidgets.QLabel("A"),1,1)
        self.MSRKa = Entrada_con_unidades(float)
        self.MSRKa.valueChanged.connect(self.setDirty)
        lyt_MSRK.addWidget(self.MSRKa,1,2)
        lyt_MSRK.addWidget(QtWidgets.QLabel("B"),2,1)
        self.MSRKb = Entrada_con_unidades(float)
        self.MSRKb.valueChanged.connect(self.setDirty)
        lyt_MSRK.addWidget(self.MSRKb,2,2)
        gridLayout_Parametricos.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Fixed),1,3)
        gridLayout_Parametricos.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Expanding),3,1)


        #Pestaña otras
        tab8 = QtWidgets.QWidget()
        tabWidget.addTab(tab8,QtCore.QCoreApplication.translate("pychemqt", "&Others"))
        gridLayout_Otros=QtWidgets.QGridLayout(tab8)
                
        gridLayout_Otros.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Solubility Parameter")),1,1)
        self.parametroSolubilidad = Entrada_con_unidades(unidades.SolubilityParameter)
        self.parametroSolubilidad.valueChanged.connect(self.setDirty)
        gridLayout_Otros.addWidget(self.parametroSolubilidad,1,2)
        gridLayout_Otros.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Dipole Moment")),2,1)
        self.momentoDipolar = Entrada_con_unidades(unidades.DipoleMoment)
        self.momentoDipolar.valueChanged.connect(self.setDirty)
        gridLayout_Otros.addWidget(self.momentoDipolar,2,2)
        gridLayout_Otros.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Molecular Diameter")),3,1)
        self.diametroMolecular = Entrada_con_unidades(float, textounidad="Å")
        self.diametroMolecular.valueChanged.connect(self.setDirty)
        gridLayout_Otros.addWidget(self.diametroMolecular,3,2)
        gridLayout_Otros.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Net calorific value")),4,1)
        self.calorCombustionNeto = Entrada_con_unidades(unidades.MolarEnthalpy)
        self.calorCombustionNeto.valueChanged.connect(self.setDirty)
        gridLayout_Otros.addWidget(self.calorCombustionNeto,4,2)
        gridLayout_Otros.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Gross calorific value")),5,1)
        self.calorCombustionBruto = Entrada_con_unidades(unidades.MolarEnthalpy)
        self.calorCombustionBruto.valueChanged.connect(self.setDirty)
        gridLayout_Otros.addWidget(self.calorCombustionBruto,5,2)
        gridLayout_Otros.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Liquid Volume Costant")),6,1)
        self.constanteVolumenLiquido = Entrada_con_unidades(unidades.MolarVolume)
        self.constanteVolumenLiquido.valueChanged.connect(self.setDirty)
        gridLayout_Otros.addWidget(self.constanteVolumenLiquido,6,2)
        gridLayout_Otros.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "API Gravity")),7,1)
        self.gravedadAPI = Entrada_con_unidades(float)
        self.gravedadAPI.valueChanged.connect(self.setDirty)
        gridLayout_Otros.addWidget(self.gravedadAPI,7,2)
        gridLayout_Otros.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Modified acentric factor")),8,1)
        self.fAcentricoModificado = Entrada_con_unidades(float)
        self.fAcentricoModificado.valueChanged.connect(self.setDirty)
        gridLayout_Otros.addWidget(self.fAcentricoModificado,8,2)
        gridLayout_Otros.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Fixed,QtWidgets.QSizePolicy.Fixed),1,3)
        gridLayout_Otros.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "UNIQUAC area")),1,4)
        self.UNIQUACArea = Entrada_con_unidades(float)
        self.UNIQUACArea.valueChanged.connect(self.setDirty)
        gridLayout_Otros.addWidget(self.UNIQUACArea,1,5)
        gridLayout_Otros.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "UNIQUAC volume")),2,4)
        self.UNIQUACVolumen = Entrada_con_unidades(float)
        self.UNIQUACVolumen.valueChanged.connect(self.setDirty)
        gridLayout_Otros.addWidget(self.UNIQUACVolumen,2,5)
        gridLayout_Otros.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Watson volume")),3,4)
        self.wilson = Entrada_con_unidades(float)
        self.wilson.valueChanged.connect(self.setDirty)
        gridLayout_Otros.addWidget(self.wilson,3,5)
        gridLayout_Otros.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Stiehl polar factor")),4,4)
        self.stiehl = Entrada_con_unidades(float)
        self.stiehl.valueChanged.connect(self.setDirty)
        gridLayout_Otros.addWidget(self.stiehl,4,5)
        gridLayout_Otros.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Rackett constant")),5,4)
        self.rackett = Entrada_con_unidades(float)
        self.rackett.valueChanged.connect(self.setDirty)
        gridLayout_Otros.addWidget(self.rackett,5,5)
        gridLayout_Otros.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Polar parameter")),6,4)
        self.parametroPolar = Entrada_con_unidades(float)
        self.parametroPolar.valueChanged.connect(self.setDirty)
        gridLayout_Otros.addWidget(self.parametroPolar,6,5)
        gridLayout_Otros.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Eps/k")),7,4)
        self.EpsK = Entrada_con_unidades(float)
        self.EpsK.valueChanged.connect(self.setDirty)
        gridLayout_Otros.addWidget(self.EpsK,7,5)
        gridLayout_Otros.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Watson factor")),8,4)
        self.watson = Entrada_con_unidades(float)
        self.watson.valueChanged.connect(self.setDirty)
        gridLayout_Otros.addWidget(self.watson,8,5)
        gridLayout_Otros.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Fixed,QtWidgets.QSizePolicy.Fixed),9,1)
        
        group_electrolitos = QtWidgets.QGroupBox(QtCore.QCoreApplication.translate("pychemqt", "Electrolitics properties"))
        gridLayout_Otros.addWidget(group_electrolitos,10,1,1,5)
        lyt_electrolitos=QtWidgets.QGridLayout(group_electrolitos)
        lyt_electrolitos.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Heat of solution")),1,1)
        self.calorDisolucion = Entrada_con_unidades(unidades.MolarEnthalpy)
        lyt_electrolitos.addWidget(self.calorDisolucion,1,2)
        lyt_electrolitos.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Solubility, 20ºC", None)),2,1)
        self.solubilidad = Entrada_con_unidades(float, textounidad="%pp")
        lyt_electrolitos.addWidget(self.solubilidad,2,2)
        lyt_electrolitos.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Fixed,QtWidgets.QSizePolicy.Fixed),1,3)
        lyt_electrolitos.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Acid-Base Characteristic")),1,4)
        self.tipo_ph = QtWidgets.QComboBox()
        self.tipo_ph.addItem(QtCore.QCoreApplication.translate("pychemqt", "Neutral"))
        self.tipo_ph.addItem(QtCore.QCoreApplication.translate("pychemqt", "Strong acid"))
        self.tipo_ph.addItem(QtCore.QCoreApplication.translate("pychemqt", "Weak acid"))
        self.tipo_ph.addItem(QtCore.QCoreApplication.translate("pychemqt", "Strong base"))
        self.tipo_ph.addItem(QtCore.QCoreApplication.translate("pychemqt", "Weak base"))
        self.tipo_ph.addItem(QtCore.QCoreApplication.translate("pychemqt", "Amphoteric"))
        lyt_electrolitos.addWidget(self.tipo_ph,1,5)
        lyt_electrolitos.addWidget(QtWidgets.QLabel("pKa"),2,4)
        self.pka = Entrada_con_unidades(float)
        lyt_electrolitos.addWidget(self.pka,2,5)
        lyt_electrolitos.addWidget(QtWidgets.QLabel("pKb"),3,4)
        self.pkb = Entrada_con_unidades(float)
        lyt_electrolitos.addWidget(self.pkb,3,5)
        lyt_electrolitos.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate("pychemqt", "Electrolitic charge")),4,4)
        self.cargaElectrolitica = Entrada_con_unidades(float)
        lyt_electrolitos.addWidget(self.cargaElectrolitica,4,5)
        gridLayout_Otros.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Expanding),11,1)

        self.change(indice)
        self.dirty=False

    def change(self, indice):            
        if indice=="-1":
            if self.indice==1001:
                indice=sql.N_comp
            else:
                indice=self.indice-1
        elif indice=="+1":
            if self.indice==sql.N_comp:
                indice=1001
            else:
                indice=self.indice+1
        elif indice=="last":
            if sql.N_comp_Custom:
                indice=1000+sql.N_comp_Custom
            else:
                indice=sql.N_comp

        if indice == 0:
            self.buttonBox.button(QtWidgets.QDialogButtonBox.Discard).setVisible(True)
            self.buttonBox.button(QtWidgets.QDialogButtonBox.Save).setVisible(True)
            self.buttonBox.button(QtWidgets.QDialogButtonBox.Apply).setVisible(False)
            self.buttonBox.button(QtWidgets.QDialogButtonBox.Close).setVisible(False)
        elif indice < 1000:
            self.buttonBox.button(QtWidgets.QDialogButtonBox.Discard).setVisible(False)
            self.buttonBox.button(QtWidgets.QDialogButtonBox.Save).setVisible(False)
            self.buttonBox.button(QtWidgets.QDialogButtonBox.Apply).setVisible(False)
            self.buttonBox.button(QtWidgets.QDialogButtonBox.Close).setVisible(True)
        else:
            self.buttonBox.button(QtWidgets.QDialogButtonBox.Discard).setVisible(True)
            self.buttonBox.button(QtWidgets.QDialogButtonBox.Save).setVisible(False)
            self.buttonBox.button(QtWidgets.QDialogButtonBox.Apply).setVisible(True)
            self.buttonBox.button(QtWidgets.QDialogButtonBox.Close).setVisible(False)
            
        if indice == 0 or indice >1000:
            self.setReadOnly(False)
        else:
            self.setReadOnly(True)

        if indice > 0:
            self.clear()
            self.buttonFirst.setDisabled(indice==1)
            self.buttonPrevious.setDisabled(indice==1)
            if sql.N_comp_Custom:
                last=1000+sql.N_comp_Custom
            else:
                last=sql.N_comp
            self.buttonNext.setDisabled(indice==last)
            self.buttonLast.setDisabled(indice==last)
            self.rellenar(indice)
        self.indice=indice

    def setDirty(self):
        self.dirty=True
        
    def clear(self):
        self.nombre.clear()
        self.nombre_alternativo.clear()
        self.CAS_number.clear()
        self.formula1.clear()
        self.formula2.clear()
        self.UNIFAC.clear()
        self.M.clear()
        self.tempCritica.clear()
        self.Pc.clear()
        self.volumenCritico.clear()
        self.temp_fusion.clear()
        self.temp_ebullicion.clear()
        self.factorAcentrico.clear()
        self.SG.clear()
        self.calorFormacionGas.clear()
        self.energiaGibbsGas.clear()
        self.cpa.clear()
        self.cpb.clear()
        self.cpc.clear()
        self.cpd.clear()
        self.cpe.clear()
        self.cpf.clear()
        self.cpgasDIPPR.clear()
        self.cpLiquidoDIPPR.clear()
        self.cpSolidoDIPPR.clear()
        self.RhoSolido.clear()
        self.RhoLiquido.clear()
        self.muLiquido.clear()
        self.muGas.clear()
        self.boilingHeat.clear()
        self.vaporPressure.clear()
        self.conductivityLiquido.clear()
        self.conductivityGas.clear()
        self.tension.clear()
        self.muParametric.clear()
        self.PvAntoine.clear()
        self.tensionParametrica.clear()
        self.Henry.clear()
        self.MSRKa.clear()
        self.MSRKb.clear()
        self.parametroSolubilidad.clear()
        self.momentoDipolar.clear()
        self.diametroMolecular.clear()
        self.calorCombustionNeto.clear()
        self.calorCombustionBruto.clear()
        self.constanteVolumenLiquido.clear()
        self.gravedadAPI.clear()
        self.fAcentricoModificado.clear()
        self.UNIQUACArea.clear()
        self.UNIQUACVolumen.clear()
        self.wilson.clear()
        self.stiehl.clear()
        self.rackett.clear()
        self.EpsK.clear()
        self.watson.clear()

        
            
    def rellenar(self, indice):
        self.indice=indice
        self.compuesto=Componente(indice)
        self.nombre.setText("%i - %s" % (self.compuesto.indice, self.compuesto.nombre))
        self.nombre_alternativo.setText(self.compuesto.nombre_alternativo)
        self.CAS_number.setText(self.compuesto.CASNumber)
        self.formula1.setText(self.compuesto.formula)
        if self.compuesto.smile!="":
            self.formula2.setPixmap(QtGui.QPixmap(self.compuesto.archivo_imagen.name))
            self.smile.setText(self.compuesto.smile)
        if self.compuesto.UNIFAC!=[]:
            self.UNIFAC.setRowCount(len(self.compuesto.UNIFAC))
            for i in range(len(self.compuesto.UNIFAC)):
                self.UNIFAC.setRowHeight(i, 25)
                self.UNIFAC.setItem(i, 0, QtWidgets.QTableWidgetItem(str(self.compuesto.UNIFAC[i][0])))
                self.UNIFAC.setItem(i, 1, QtWidgets.QTableWidgetItem(str(self.compuesto.UNIFAC[i][1])))
                self.UNIFAC.item(i, 1).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
            
        if self.compuesto.M:
            self.M.setValue(self.compuesto.M)
        if self.compuesto.Tc:
            self.tempCritica.setValue(self.compuesto.Tc)
        if self.compuesto.Pc:
            self.Pc.setValue(self.compuesto.Pc)
        if self.compuesto.Vc:
            self.volumenCritico.setValue(self.compuesto.Vc)
        if self.compuesto.Tf:
            self.temp_fusion.setValue(self.compuesto.Tf)
        if self.compuesto.Tb:
            self.temp_ebullicion.setValue(self.compuesto.Tb)
        if self.compuesto.f_acent:
            self.factorAcentrico.setValue(self.compuesto.f_acent)
        if self.compuesto.SG:
            self.SG.setValue(self.compuesto.SG)
        if self.compuesto.calor_formacion:
            self.calorFormacionGas.setValue(self.compuesto.calor_formacion)
        if self.compuesto.energia_formacion:
            self.energiaGibbsGas.setValue(self.compuesto.energia_formacion)
        self.cpa.setValue(self.compuesto.cp[0])
        self.cpb.setValue(self.compuesto.cp[1])
        self.cpc.setValue(self.compuesto.cp[2])
        self.cpd.setValue(self.compuesto.cp[3])
        self.cpe.setValue(self.compuesto.cp[4])
        self.cpf.setValue(self.compuesto.cp[5])
        
        self.cpgasDIPPR.rellenar(self.compuesto.capacidad_calorifica_gas, indice)
        self.cpLiquidoDIPPR.rellenar(self.compuesto.capacidad_calorifica_liquido, indice)
        self.cpSolidoDIPPR.rellenar(self.compuesto.capacidad_calorifica_solido, indice)
        self.RhoSolido.rellenar(self.compuesto.densidad_solido, indice)
        self.RhoLiquido.rellenar(self.compuesto.densidad_liquido, indice)
        self.muLiquido.rellenar(self.compuesto.viscosidad_liquido, indice)
        self.muGas.rellenar(self.compuesto.viscosidad_gas, indice)
        self.boilingHeat.rellenar(self.compuesto.calor_vaporizacion, indice)
        self.vaporPressure.rellenar(self.compuesto.presion_vapor, indice)
        self.conductivityLiquido.rellenar(self.compuesto.conductividad_liquido, indice)
        self.conductivityGas.rellenar(self.compuesto.conductividad_gas, indice)
        self.tension.rellenar(self.compuesto.tension_superficial, indice)
        
        self.muParametric.rellenar(self.compuesto.viscosidad_parametrica, indice)
        self.PvAntoine.rellenar(self.compuesto.antoine, indice)
        self.tensionParametrica.rellenar(self.compuesto.tension_superficial_parametrica, indice)
        self.Henry.rellenar(self.compuesto.henry, indice)

        if self.compuesto.MSRK[0]!=0 and self.compuesto.MSRK[1]!=0:
            self.MSRKa.setValue(self.compuesto.MSRK[0])
            self.MSRKb.setValue(self.compuesto.MSRK[1])
        if self.compuesto.parametro_solubilidad:
            self.parametroSolubilidad.setValue(self.compuesto.parametro_solubilidad)
        if self.compuesto.momento_dipolar:
            self.momentoDipolar.setValue(self.compuesto.momento_dipolar)
        if self.compuesto.diametro_molecular:
            self.diametroMolecular.setValue(self.compuesto.diametro_molecular)
        if self.compuesto.calor_combustion_neto:
            self.calorCombustionNeto.setValue(self.compuesto.calor_combustion_neto)
        if self.compuesto.calor_combustion_bruto:
            self.calorCombustionBruto.setValue(self.compuesto.calor_combustion_bruto)
        if self.compuesto.Vliq:
            self.constanteVolumenLiquido.setValue(self.compuesto.Vliq)
        if self.compuesto.API:
            self.gravedadAPI.setValue(self.compuesto.API)
        if self.compuesto.f_acent_mod:
            self.fAcentricoModificado.setValue(self.compuesto.f_acent_mod)
        if self.compuesto.UNIQUAC_area:
            self.UNIQUACArea.setValue(self.compuesto.UNIQUAC_area)
        if self.compuesto.UNIQUAC_volumen:
            self.UNIQUACVolumen.setValue(self.compuesto.UNIQUAC_volumen)
        if self.compuesto.wilson:
            self.wilson.setValue(self.compuesto.wilson)
        if self.compuesto.stiehl:
            self.stiehl.setValue(self.compuesto.stiehl)
        if self.compuesto.rackett:
            self.rackett.setValue(self.compuesto.rackett)
    #        if self.compuesto.[53]:
    #            self.parametroPolar.setValueself.compuesto.[53])
        if self.compuesto.ek:
            self.EpsK.setValue(self.compuesto.ek)
        if self.compuesto.Kw:
            self.watson.setValue(self.compuesto.Kw)
        #        if self.compuesto.[52]:
        #            self.cargaElectrolitica.setText(str(self.compuesto.[52]))
        #        if self.compuesto.[54]:
        #            self.calorDisolucion.setText(str(self.compuesto.[54]))
        #        if self.compuesto.[55]:
        #            self.solubilidad.setText(str(self.compuesto.[55]))
        #        if self.compuesto.[57]:
        #            self.pka.setText(str(self.compuesto.[57]))
        #        if self.compuesto.[58]:
        #            self.pkb.setText(str(self.compuesto.[58]))
        #        self.tipo_ph.setCurrentIndex(self.compuesto.[56])

    def getComponent(self):
        nuevo_elemento=[]
        nuevo_elemento.append(str(self.formula1.text()))
        if self.indice==0:
            nuevo_elemento.append(str(self.nombre.text()))
        else:
            nuevo_elemento.append(str(self.nombre.text()).split(" - ")[1])
        
        nuevo_elemento.append(self.M.value)
        nuevo_elemento.append(self.tempCritica.value)
        nuevo_elemento.append(self.Pc.value)
        nuevo_elemento.append(self.volumenCritico.value)
        nuevo_elemento.append(self.gravedadAPI.value)
            
        cpideal=[]
        cpideal.append(self.cpa.value)
        cpideal.append(self.cpb.value)
        cpideal.append(self.cpc.value)
        cpideal.append(self.cpd.value)
        cpideal.append(self.cpe.value)
        cpideal.append(self.cpf.value)
        nuevo_elemento.append(cpideal)
        
        nuevo_elemento.append(self.PvAntoine.value)
        nuevo_elemento.append(self.Henry.value)
        nuevo_elemento.append(self.muParametric.value)
        nuevo_elemento.append(self.tensionParametrica.value)
        
        nuevo_elemento.append(self.RhoSolido.value)
        nuevo_elemento.append(self.RhoLiquido.value)
        nuevo_elemento.append(self.vaporPressure.value)
        nuevo_elemento.append(self.boilingHeat.value)
        nuevo_elemento.append(self.cpSolidoDIPPR.value)
        nuevo_elemento.append(self.cpLiquidoDIPPR.value)
        nuevo_elemento.append(self.cpgasDIPPR.value)
        nuevo_elemento.append(self.muLiquido.value)
        nuevo_elemento.append(self.muGas.value)
        nuevo_elemento.append(self.conductivityLiquido.value)
        nuevo_elemento.append(self.conductivityGas.value)
        nuevo_elemento.append(self.tension.value)

        nuevo_elemento.append(self.momentoDipolar.value)
        nuevo_elemento.append(self.constanteVolumenLiquido.value)
        nuevo_elemento.append(self.rackett.value)
        nuevo_elemento.append(self.SG.value)
        nuevo_elemento.append(self.factorAcentrico.value)
        nuevo_elemento.append(self.parametroSolubilidad.value)
        nuevo_elemento.append(self.watson.value)
        
        msrk=[]
        msrk.append(self.MSRKa.value)
        msrk.append(self.MSRKb.value)
        nuevo_elemento.append(msrk)
        
        nuevo_elemento.append(self.stiehl.value)
        nuevo_elemento.append(self.temp_ebullicion.value)
        nuevo_elemento.append(self.temp_fusion.value)
        nuevo_elemento.append(str(self.CAS_number.text()))
        nuevo_elemento.append(str(self.formula2.text()))

        UNIFAC=self.UNIFAC.getMatrix()
        nuevo_elemento.append(UNIFAC)

        nuevo_elemento.append(self.diametroMolecular.value)
        nuevo_elemento.append(self.EpsK.value)
        nuevo_elemento.append(self.UNIQUACArea.value)
        nuevo_elemento.append(self.UNIQUACVolumen.value)
        nuevo_elemento.append(self.fAcentricoModificado.value)
        nuevo_elemento.append(self.calorFormacionGas.value)
        nuevo_elemento.append(self.energiaGibbsGas.value)
        nuevo_elemento.append(self.wilson.value)
        nuevo_elemento.append(self.calorCombustionNeto.value)
        nuevo_elemento.append(self.calorCombustionBruto.value)

        nuevo_elemento.append(str(self.nombre_alternativo.text()))
        nuevo_elemento.append(0)
        nuevo_elemento.append(0)
        nuevo_elemento.append(0)

        nuevo_elemento.append(self.cargaElectrolitica.value)
        nuevo_elemento.append(self.parametroPolar.value)
        nuevo_elemento.append(self.calorDisolucion.value)
        nuevo_elemento.append(self.solubilidad.value)

        nuevo_elemento.append(self.tipo_ph.currentIndex())
        nuevo_elemento.append(self.pka.value)
        nuevo_elemento.append(self.pkb.value)
        nuevo_elemento.append(str(self.smile.text()))
        
        return nuevo_elemento

    def buttonClicked(self, boton):
        role=self.buttonBox.buttonRole(boton)
        if role==QtWidgets.QDialogButtonBox.AcceptRole:
            componente=self.getComponent()
            sql.inserElementsFromArray(sql.databank_Custom_name, [componente])
            self.accept()
        elif role==QtWidgets.QDialogButtonBox.ApplyRole:
            componente=self.getComponent()
            sql.updateElement(componente, self.indice)
            self.accept()
        elif role==QtWidgets.QDialogButtonBox.DestructiveRole:
            componente=self.getComponent()
            if self.indice==0:
                func=sql.inserElementsFromArray
                arg=(sql.databank_Custom_name, [componente])
            elif self.indice>1000:
                func=sql.updateElement
                arg=(componente, self.indice)
            if okToContinue(self, self.dirty, func, arg):
                self.accept()
            else:
                self.reject()
        else:
                self.reject()
                


    def setReadOnly(self, bool):
        self.formula1.setReadOnly(bool)
        self.CAS_number.setReadOnly(bool)
        self.nombre_alternativo.setReadOnly(bool)
        self.formula1.setReadOnly(bool)
        self.smile.setReadOnly(bool)
        self.M.setReadOnly(bool)
        self.tempCritica.setReadOnly(bool)
        self.Pc.setReadOnly(bool)
        self.volumenCritico.setReadOnly(bool)
        self.factorAcentrico.setReadOnly(bool)
        self.SG.setReadOnly(bool)
        self.temp_ebullicion.setReadOnly(bool)
        self.temp_fusion.setReadOnly(bool)
        self.calorFormacionGas.setReadOnly(bool)
        self.energiaGibbsGas.setReadOnly(bool)
        self.cpa.setReadOnly(bool)
        self.cpb.setReadOnly(bool)
        self.cpc.setReadOnly(bool)
        self.cpd.setReadOnly(bool)
        self.cpe.setReadOnly(bool)
        self.cpf.setReadOnly(bool)
        self.MSRKa.setReadOnly(bool)
        self.MSRKb.setReadOnly(bool)
        self.parametroSolubilidad.setReadOnly(bool)
        self.momentoDipolar.setReadOnly(bool)
        self.diametroMolecular.setReadOnly(bool)
        self.calorCombustionNeto.setReadOnly(bool)
        self.calorCombustionBruto.setReadOnly(bool)
        self.constanteVolumenLiquido.setReadOnly(bool)
        self.gravedadAPI.setReadOnly(bool)
        self.fAcentricoModificado.setReadOnly(bool)
        self.UNIQUACArea.setReadOnly(bool)
        self.UNIQUACVolumen.setReadOnly(bool)
        self.wilson.setReadOnly(bool)
        self.stiehl.setReadOnly(bool)
        self.rackett.setReadOnly(bool)
        self.parametroPolar.setReadOnly(bool)
        self.EpsK.setReadOnly(bool)
        self.watson.setReadOnly(bool)
        self.calorDisolucion.setReadOnly(bool)
        self.solubilidad.setReadOnly(bool)
        self.pka.setReadOnly(bool)
        self.pkb.setReadOnly(bool)
        self.cargaElectrolitica.setReadOnly(bool)
        self.tipo_ph.setDisabled(bool)

        
if __name__ == "__main__":
    import sys
    from lib.petro import Petroleo
    app = QtWidgets.QApplication(sys.argv)
#    Dialog = View_Component(10)
#    petroleo=Petroleo(name="Petroleo", API=22.5, M=339.7)
#    Dialog = View_Petro(petroleo)
    Dialog = View_Contribution()
    Dialog.show()
    sys.exit(app.exec_())
