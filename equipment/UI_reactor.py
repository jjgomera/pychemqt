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


###################################
###   Diálogo de definición de reactores UI_reactor  ###
###################################

import pickle
from functools import partial
from math import exp, log
import os

from PyQt5 import QtCore, QtGui, QtWidgets

from scipy.optimize import leastsq
from scipy import array

from lib import unidades, reaction, plot
from lib.thread import Evaluate
from lib.config import getComponents
from UI.widgets import Status
from equipment.parents import UI_equip
from equipment.reactor import Reactor
from UI import UI_corriente, inputTable
from UI.widgets import Entrada_con_unidades, Tabla, QLabelMath


class widgetReacciones(QtWidgets.QWidget):
    """Widget con la tabla de reacciones y los botones para modificar la lista de reacciones"""
    changed = QtCore.pyqtSignal()
    reacciones=[]
    reaccion=None
    activo=None
    ajuste=None

    def __init__(self, parent=None):
        super(widgetReacciones, self).__init__(parent)
        self.indices, self.nombres, M=getComponents()
        gridLayout = QtWidgets.QGridLayout(self)

        self.TablaReacciones=Tabla(5, horizontalHeader=[QtWidgets.QApplication.translate("pychemqt", "Reaction"), "ΔHr, %s" %unidades.MolarEnthalpy(None).text(), QtWidgets.QApplication.translate("pychemqt", "Type"), QtWidgets.QApplication.translate("pychemqt", "Phase"), QtWidgets.QApplication.translate("pychemqt", "Description")], dinamica=False, verticalHeader=True, orientacion=QtCore.Qt.AlignmentFlag.AlignLeft)
        self.TablaReacciones.setMinimumWidth(500)
        self.TablaReacciones.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.TablaReacciones.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        self.TablaReacciones.horizontalHeader().setStretchLastSection(True)
        self.TablaReacciones.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.TablaReacciones.itemSelectionChanged.connect(self.actualizarBotones)
        gridLayout.addWidget(self.TablaReacciones,1,1,6,4)

        self.botonAbrir=QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/fileOpen.png")), QtWidgets.QApplication.translate("pychemqt", "Open"))
        self.botonAbrir.clicked.connect(self.botonAbrirClicked)
        gridLayout.addWidget(self.botonAbrir,1,5)
        self.botonGuardar=QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/fileSave.png")), QtWidgets.QApplication.translate("pychemqt", "Save"))
        self.botonGuardar.clicked.connect(self.botonGuardarClicked)
        self.botonGuardar.setSizePolicy(QtWidgets.QSizePolicy.Policy.Fixed,QtWidgets.QSizePolicy.Policy.Fixed)
        self.botonGuardar.setEnabled(False)
        gridLayout.addWidget(self.botonGuardar,2,5)

        self.botonNew=QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/fileNew.png")), QtWidgets.QApplication.translate("pychemqt", "New"))
        self.botonNew.clicked.connect(self.botonNewClicked)
        gridLayout.addWidget(self.botonNew,3,5)
        self.botonEdit=QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/editor.png")), QtWidgets.QApplication.translate("pychemqt", "Edit"))
        self.botonEdit.setEnabled(False)
        self.botonEdit.setCheckable(True)
        self.botonEdit.clicked.connect(self.botonEditClicked)
        gridLayout.addWidget(self.botonEdit,4,5)
        self.botonDelete=QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/editDelete.png")), QtWidgets.QApplication.translate("pychemqt", "Delete"))
        self.botonDelete.setEnabled(False)
        self.botonDelete.clicked.connect(self.botonDeleteClicked)
        gridLayout.addWidget(self.botonDelete,5,5)
        self.botonClear=QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/clear.png")), QtWidgets.QApplication.translate("pychemqt", "Clear"))
        self.botonClear.clicked.connect(self.botonClearClicked)
        gridLayout.addWidget(self.botonClear,6,5)
        gridLayout.addItem(QtWidgets.QSpacerItem(10,10,QtWidgets.QSizePolicy.Policy.Expanding,QtWidgets.QSizePolicy.Policy.Expanding),10,1)


    def actualizarBotones(self, bool=True):
        self.botonEdit.setEnabled(bool)
        self.botonDelete.setEnabled(bool)

    def botonAbrirClicked(self):
        fname = str(QtWidgets.QFileDialog.getOpenFileName(self, QtWidgets.QApplication.translate("pychemqt", "Open reaction file"), "./", QtWidgets.QApplication.translate("pychemqt", "reaction file")+" (*.rec);;"+QtWidgets.QApplication.translate("pychemqt", "All files")+" (*.*)")[0])
        if fname:
            with open(fname, "r") as archivo:
                reacciones=pickle.load(archivo)
            print(reacciones)
            self.reacciones=reacciones
            self.botonGuardar.setEnabled(True)
            for fila, reaccion in enumerate(reacciones):
                self.TablaReacciones.addRow()
                self.TablaReacciones.setValue(fila, 0, reaccion.text)
                self.TablaReacciones.setValue(fila, 1, "%0.4e" %reaccion.Hr.config(), QtCore.Qt.AlignmentFlag.AlignRight)
                self.TablaReacciones.setValue(fila, 2, str(reaccion.tipo+1)+" - "+reaction.Reaction.TEXT_TYPE[reaccion.tipo])
                self.TablaReacciones.setValue(fila, 3, reaction.Reaction.TEXT_PHASE[reaccion.fase])
                self.TablaReacciones.item(fila, 4).setFlags(QtCore.Qt.ItemFlag.ItemIsEditable|QtCore.Qt.ItemFlag.ItemIsEnabled|QtCore.Qt.ItemFlag.ItemIsSelectable)
            for i in range(4):
                self.TablaReacciones.resizeColumnToContents(i)
        self.changed.emit()

    def botonGuardarClicked(self):
        fname = str(QtWidgets.QFileDialog.getSaveFileName(self, QtWidgets.QApplication.translate("pychemqt", "Save reaction to file"), "./", QtWidgets.QApplication.translate("pychemqt", "reaction file")+" (*.rec)")[0])
        if fname:
            if fname.split(".")[-1]!="rec":
                fname+=".rec"
            pickle.dump(self.reacciones, open(fname, "w"))

    def botonNewClicked(self):
        dialog=UI_reacciones(parent=self)
        if dialog.exec():
            pass


    def botonEditClicked(self, bool):
        if bool:
            indice=self.TablaReacciones.currentRow()
            reaccion=self.reacciones[indice]
            dialogo=UI_reacciones(reaccion, self)
            dialogo.exec()
#            self.rellenar(self.reaccion)
#            self.activo=indice
        else:
            self.botonAddClicked(self.activo, False)
            self.reacciones[self.activo]=self.reaccion
            self.TablaReacciones.setCurrentCell(self.activo, 0)
            self.activo=-1
            self.changed.emit()

        self.botonNew.setEnabled(not bool)
        self.botonDelete.setEnabled(not bool)
        self.botonClear.setEnabled(not bool)
        self.botonAdd.setEnabled(not bool)
        self.botonAbrir.setEnabled(not bool)
        self.botonGuardar.setEnabled(not bool)


    def botonDeleteClicked(self):
        indice=self.TablaReacciones.currentRow()
        self.TablaReacciones.removeRow(indice)
        del self.reacciones[indice]
        self.TablaReacciones.clearSelection()
        self.actualizarBotones(False)
        self.changed.emit()

    def botonClearClicked(self):
        if self.reacciones:
            self.reacciones=[]
            self.TablaReacciones.setRowCount(0)
            self.botonGuardar.setEnabled(False)

    def botonAddClicked(self, fila, add=True):
        if add:
            fila=self.TablaReacciones.rowCount()
            self.TablaReacciones.addRow()
        self.TablaReacciones.setValue(fila, 0, self.Formula.text())
        self.TablaReacciones.setValue(fila, 1, "%0.4e" %self.Hr.value.config(), QtCore.Qt.AlignmentFlag.AlignRight)
        self.TablaReacciones.setValue(fila, 2, str(self.tipo.currentIndex()+1)+" - "+self.tipo.currentText())
        self.TablaReacciones.setValue(fila, 3, self.Fase.currentText())
        self.TablaReacciones.item(fila, 4).setFlags(QtCore.Qt.ItemFlag.ItemIsEditable|QtCore.Qt.ItemFlag.ItemIsEnabled|QtCore.Qt.ItemFlag.ItemIsSelectable)
        for i in range(4):
            self.TablaReacciones.resizeColumnToContents(i)
        self.reacciones.insert(fila, self.reaccion)
        self.botonGuardar.setEnabled(True)
        self.changed.emit()

class UI_reacciones(QtWidgets.QDialog):
    reaction=reaction.Reaction()
    def __init__(self, reaccion=None, parent=None):
        super(UI_reacciones, self).__init__(parent)
        self.evaluate=Evaluate()
        self.evaluate.finished.connect(self.rellenar)
        self.indices, self.nombres, M=getComponents()
        gridLayout = QtWidgets.QGridLayout(self)

        lyt=QtWidgets.QHBoxLayout()
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Key component")))
        self.key=QtWidgets.QComboBox()
        for i, nombre in enumerate(self.nombres):
            self.key.addItem("%i - %s" %(i+1, nombre))
        self.key.currentIndexChanged.connect(partial(self.changeParams, "key"))
        lyt.addWidget(self.key)
        lyt.addItem(QtWidgets.QSpacerItem(10,10,QtWidgets.QSizePolicy.Policy.Expanding,QtWidgets.QSizePolicy.Policy.Fixed))
        gridLayout.addLayout(lyt,1,1,1,5)

        lyt=QtWidgets.QHBoxLayout()
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Phase")))
        self.fase=QtWidgets.QComboBox()
        for txt in reaction.Reaction.TEXT_PHASE:
            self.fase.addItem(txt)
        self.fase.currentIndexChanged.connect(partial(self.changeParams, "fase"))
        lyt.addWidget(self.fase)
        self.Formula=QtWidgets.QLabel()
        self.Formula.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.Formula.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding,QtWidgets.QSizePolicy.Policy.Fixed)
        lyt.addWidget(self.Formula)
        gridLayout.addLayout(lyt,2,1,1,5)

        lyt=QtWidgets.QVBoxLayout()
        title=self.nombres[:]
        title.append("")
        self.Estequiometria=Tabla(1, verticalHeaderLabels=title, horizontalHeader=[QtWidgets.QApplication.translate("pychemqt", "Coefficients")], filas=len(self.indices))
        self.Estequiometria.setFixedHeight(22*len(self.indices)+22+4+22)
        lyt.addWidget(self.Estequiometria)
        self.Estequiometria.addRow()
        brush=QtGui.QBrush(QtGui.QColor("#eaeaea"))
        self.Estequiometria.item(len(self.indices), 0).setBackground(brush)
        self.Estequiometria.item(len(self.indices), 0).setFlags(QtCore.Qt.ItemFlag.NoItemFlags)
        self.Estequiometria.cellChanged.connect(self.reaccionCambiada)
        lyt.addItem(QtWidgets.QSpacerItem(10,10,QtWidgets.QSizePolicy.Policy.Fixed,QtWidgets.QSizePolicy.Policy.Expanding))
        gridLayout.addLayout(lyt,3,1,1,2)

        lyt=QtWidgets.QGridLayout()
        lyt.addItem(QtWidgets.QSpacerItem(10,10,QtWidgets.QSizePolicy.Policy.Fixed,QtWidgets.QSizePolicy.Policy.Fixed),1,1)
        self.formula=QtWidgets.QCheckBox(QtWidgets.QApplication.translate("pychemqt", "Use name in formula"))
        self.formula.toggled.connect(partial(self.changeParams, "formula"))
        lyt.addWidget(self.formula,1,2,1,2)
        self.customHr=QtWidgets.QCheckBox("ΔHr "+QtWidgets.QApplication.translate("pychemqt", "user specified"))
        self.customHr.toggled.connect(self.changeHr)
        lyt.addWidget(self.customHr,2,2,1,2)
        lyt.addWidget(QtWidgets.QLabel("ΔHr<sup>o</sup>"),3,2)
        self.Hr=Entrada_con_unidades(unidades.MolarEnthalpy, readOnly=True)
        self.Hr.valueChanged.connect(partial(self.changeParams, "Hr"))
        lyt.addWidget(self.Hr,3,3)
        lyt.addItem(QtWidgets.QSpacerItem(10,10,QtWidgets.QSizePolicy.Policy.Fixed,QtWidgets.QSizePolicy.Policy.Expanding))
        gridLayout.addLayout(lyt,3,3,1,2)

        gridLayout.addItem(QtWidgets.QSpacerItem(10,10,QtWidgets.QSizePolicy.Policy.Fixed,QtWidgets.QSizePolicy.Policy.Fixed),4,2)

        lyt=QtWidgets.QHBoxLayout()
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Type")))
        self.tipo=QtWidgets.QComboBox()
        for txt in reaction.Reaction.TEXT_TYPE:
            self.tipo.addItem(txt)
        self.tipo.currentIndexChanged.connect(partial(self.changeParams, "tipo"))
        lyt.addWidget(self.tipo)
        lyt.addItem(QtWidgets.QSpacerItem(20,10,QtWidgets.QSizePolicy.Policy.Fixed,QtWidgets.QSizePolicy.Policy.Fixed))
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Concentration")))
        self.base=QtWidgets.QComboBox()
        for txt in reaction.Reaction.TEXT_BASE:
            self.base.addItem(txt)
        self.base.currentIndexChanged.connect(partial(self.changeParams, "base"))
        lyt.addWidget(self.base)
        lyt.addItem(QtWidgets.QSpacerItem(10,10,QtWidgets.QSizePolicy.Policy.Expanding,QtWidgets.QSizePolicy.Policy.Fixed))
        gridLayout.addLayout(lyt,5,1,1,5)

        self.stacked = QtWidgets.QStackedWidget()
        self.tipo.currentIndexChanged.connect(self.stacked.setCurrentIndex)
        gridLayout.addWidget(self.stacked,6,1,1,5)
        gridLayout.addItem(QtWidgets.QSpacerItem(10,10,QtWidgets.QSizePolicy.Policy.Expanding,QtWidgets.QSizePolicy.Policy.Expanding),7,1,1,5)

        widget=QtWidgets.QWidget()
        self.stacked.addWidget(widget)
        lyt=QtWidgets.QGridLayout(widget)
        lyt.addWidget(QtWidgets.QLabel("<h3>"+QtWidgets.QApplication.translate("pychemqt", "Estequiometric reaction")+"</h3>"),1,1,1,4)
        self.Conversion=Tabla(1, verticalHeaderModel="C", filas=3)
        self.Conversion.setConnected()
        self.Conversion.setFixedWidth(100)
        lyt.addWidget(self.Conversion,2,1,3,1)
        mathTex = r"$Conversion = C_o + C_1T + C_2T^2 + \cdots + C_nT^n$"
        label = QLabelMath(mathTex)
        lyt.addWidget(label,2,2,1,3)
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Temperature unit")),3,2)
        self.unidadesTemperatura=QtWidgets.QComboBox()
        for i in unidades.Temperature.__text__:
            self.unidadesTemperatura.addItem(i)
        lyt.addWidget(self.unidadesTemperatura,3,3)
        lyt.addItem(QtWidgets.QSpacerItem(10,10,QtWidgets.QSizePolicy.Policy.Expanding,QtWidgets.QSizePolicy.Policy.Expanding),4,4)
        lyt.addItem(QtWidgets.QSpacerItem(10,10,QtWidgets.QSizePolicy.Policy.Expanding,QtWidgets.QSizePolicy.Policy.Expanding),5,1,1,5)


        widget=QtWidgets.QWidget()
        self.stacked.addWidget(widget)
        lyt=QtWidgets.QGridLayout(widget)
        self.check_KFijo=QtWidgets.QRadioButton(QtWidgets.QApplication.translate("pychemqt", "Fixed"))
        self.check_KFijo.toggled.connect(self.KeqChanged)
        lyt.addWidget(self.check_KFijo,1,1,1,2)
        lyt.addWidget(QtWidgets.QLabel("K<sub>eq</sub>"),1,3)
        self.Keq=Entrada_con_unidades(float)
        lyt.addWidget(self.Keq,1,4)
        mathTex = r"$aA + bB \rightleftharpoons cC + dD \therefore "
        mathTex += r"K_{eq} = \frac{[C]^c [D]^d}{[A]^a [B]^b}$"
        label = QLabelMath(mathTex)
        lyt.addWidget(label,1,5,1,4)

        self.check_KEq=QtWidgets.QRadioButton(QtWidgets.QApplication.translate("pychemqt", "Equation"))
        self.check_KEq.toggled.connect(self.KeqChanged)
        lyt.addWidget(self.check_KEq,2,1,1,2)
        self.check_KTabla=QtWidgets.QRadioButton(QtWidgets.QApplication.translate("pychemqt", "Table"))
        self.check_KTabla.toggled.connect(self.KeqChanged)
        lyt.addWidget(self.check_KTabla,2,5,1,2)
        self.KEq_Dat=Tabla(1, verticalHeaderLabels=["A", "B", "C", "D", "E", "F", "G", "H"], filas=8)
        self.KEq_Dat.setFixedHeight(22*8+4)
        self.KEq_Dat.setFixedWidth(120)
        lyt.addWidget(self.KEq_Dat,3,3,1,2)
        self.KEq_Tab=Tabla(4, horizontalHeader=["T, K", "Keq", "Kcalc", "%Error"], verticalHeader=False, columnReadOnly=[False, False, True, True])
        self.KEq_Tab.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding,QtWidgets.QSizePolicy.Policy.Expanding)
        self.KEq_Tab.setFixedWidth(400)
        self.KEq_Tab.setConnected()
        self.KEq_Tab.rowFinished.connect(self.Regresion)
        self.KEq_Tab.setAlternatingRowColors(False)
        lyt.addWidget(self.KEq_Tab,3,5,1,4)
        lyt.addWidget(QtWidgets.QLabel("r²"),4,5)
        self.r2=Entrada_con_unidades(float, readOnly=True)
        lyt.addWidget(self.r2,4,6)
        self.botonTablaPlot=QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/plot.png")), QtWidgets.QApplication.translate("pychemqt", "Plot"))
        self.botonTablaPlot.clicked.connect(self.Plot)
        lyt.addWidget(self.botonTablaPlot,4,7)
        self.botonTablaClear=QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/clear.png")), QtWidgets.QApplication.translate("pychemqt", "Clear"))
        self.botonTablaClear.clicked.connect(self.KEq_Tab.clear)
        lyt.addWidget(self.botonTablaClear,4,8)
        mathTex = r"$\lnK_eq = A+B/T+C\lnT+DT+ET^2+FT^3+GT^4+HT^5$"
        label=QLabelMath(mathTex)
        lyt.addWidget(label,5,1,1,8)

        self.checkGibbs=QtWidgets.QRadioButton(QtWidgets.QApplication.translate("pychemqt", "From Gibbs free energy minimization"))
        lyt.addWidget(self.checkGibbs,6,1,1,4)

        self.check_KFijo.setChecked(True)


        widget=QtWidgets.QWidget()
        self.stacked.addWidget(widget)
        lyt=QtWidgets.QGridLayout(widget)

        widget=QtWidgets.QWidget()
        self.stacked.addWidget(widget)
        lyt=QtWidgets.QGridLayout(widget)

        self.status=Status()
        gridLayout.addWidget(self.status, 10,1)
        self.buttonBox = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Cancel|QtWidgets.QDialogButtonBox.StandardButton.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        gridLayout.addWidget(self.buttonBox,10,2,1,4)

        if reaccion:
            self.setReaction(reaccion)

    def changeParams(self, parametro, valor):
        self.calculo(**{parametro: valor})

    def calculo(self, **kwargs):
        self.status.setState(4)
        self.evaluate.start(self.reaction, kwargs)

    def changeHr(self, bool):
        self.Hr.setReadOnly(not bool)
        self.changeParams("customHr", bool)

    def reaccionCambiada(self):
        kwargs={"componentes": self.indices,
                    "coeficientes": self.Estequiometria.getColumn(0)[:-1]}
        self.calculo(**kwargs)

    def setReaction(self, reaction):
        self.reaction=reaction
        self.rellenar()

#        if self.Estequiometria.getValue(0, self.Base.currentIndex()):
#            reaccion=reaction.Reaction(self.indices, self.Estequiometria.getColumn(0), base=self.Base.currentIndex(), estequiometria=[0, 0, 0.5], formulas=self.checkFormula.isChecked(), calor=self.checkCalorEspecificado.isChecked(), Hr=self.Hr.value, tipo=self.tipo.currentIndex(), conversion=self.Conversion.getColumn(0)[-1::-1])
#            self.Balance.setValue(reaccion.error)
#            if reaccion.state:
#                self.Formula.setText(reaccion._txt(self.checkFormula.isChecked()))
#                self.Hr.setValue(reaccion.Hr)
#            else:
#                self.Formula.clear()
#                self.Hr.clear()
#            self.botonAdd.setEnabled(reaccion.state and not self.botonEdit.isChecked())
#            self.reaccion=reaccion

    def rellenar(self):
        self.blockSignals(True)
        for variable in self.reaction.kwargsValue:
            self.__getattribute__(variable).setValue(self.reaction.kwargs[variable])
        for combo in self.reaction.kwargsList:
            self.__getattribute__(combo).setCurrentIndex(self.reaction.kwargs[combo])
        for check in self.reaction.kwargsCheck:
            self.__getattribute__(check).setChecked(self.reaction.kwargs[check])

        self.Estequiometria.setColumn(0, self.reaction.kwargs["coeficientes"])
#        self.Conversion.setColumn(0, self.reaction.estequiometria[-1::-1])
        self.blockSignals(False)

        self.status.setState(self.reaction.status, self.reaction.msg)
        self.Estequiometria.item(len(self.indices), 0).setText(str(self.reaction.error))
        if self.reaction.status:
            self.Formula.setText(self.reaction._txt())
            self.Hr.setValue(self.reaction.Hr)

    def KeqChanged(self):
        self.Keq.setReadOnly(not self.check_KFijo.isChecked())
        self.KEq_Dat.setEnabled(self.check_KEq.isChecked())
        self.KEq_Tab.setEnabled(self.check_KTabla.isChecked())
        self.botonTablaClear.setEnabled(self.check_KTabla.isChecked())
        self.botonTablaPlot.setEnabled(self.check_KTabla.isChecked())

    def Regresion(self):
        t=array(self.KEq_Tab.getColumn(0)[:-1])
        k=array(self.KEq_Tab.getColumn(1)[:-1])
        if len(t)>=4:
            if 4<=len(t)<8:
                inicio=r_[0, 0, 0, 0]
                f=lambda par, T: exp(par[0]+par[1]/T+par[2]*log(T)+par[3]*T)
                resto=lambda par, T, k: k-f(par, T)
            else:
                inicio=r_[0, 0, 0, 0, 0, 0, 0, 0]
                f=lambda par, T: exp(par[0]+par[1]/T+par[2]*log(T)+par[3]*T+par[4]*T**2+par[5]*T**3+par[6]*T**4+par[7]*T**5)
                resto=lambda par, T, k: k-f(par, T)

            ajuste=leastsq(resto,inicio,args=(t, k))
            kcalc=f(ajuste[0], t)
            error=(k-kcalc)/k*100
            self.KEq_Dat.setColumn(0, ajuste[0])
            self.KEq_Tab.setColumn(2, kcalc)
            self.KEq_Tab.setColumn(3, error)

            if ajuste[1] in [1, 2, 3, 4]:
                self.ajuste=ajuste[0]

    def Plot(self):
        if self.ajuste!=None:
            t=array(self.KEq_Tab.getColumn(0)[:-1])
            k=array(self.KEq_Tab.getColumn(1)[:-1])
            if 4<=len(t)<8:
                f=lambda par, T: exp(par[0]+par[1]/T+par[2]*log(T)+par[3]*T)
            else:
                f=lambda par, T: exp(par[0]+par[1]/T+par[2]*log(T)+par[3]*T+par[4]*T**2+par[5]*T**3+par[6]*T**4+par[7]*T**5)

            grafico=plot.Plot()
            grafico.data(t, k, 'ro', t, f(self.ajuste, t))
            grafico.exec()



class UI_equipment(UI_equip):
    """Diálogo de definición de reactores"""
    profile_T=None
    Equipment=Reactor()
    def __init__(self, equipment=None, parent=None):
        """entrada: Parametro opcional de clase corriente que indica la corriente de entrada en el reactor"""
        super(UI_equipment, self).__init__(Reactor, entrada=False, salida=False, parent=parent)

        #Pestaña reacciones
        self.Reacciones= widgetReacciones()
        self.Reacciones.changed.connect(self.calculo)
        self.tabWidget.insertTab(1, self.Reacciones, QtWidgets.QApplication.translate("pychemqt", "Reactions"))

        #Pestaña calculo
        gridLayout_Calculo = QtWidgets.QGridLayout(self.tabCalculo)
        gridLayout_Calculo.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Pressure")),2,0,1,1)
        self.P=Entrada_con_unidades(unidades.Pressure)
        self.P.valueChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.P,2,1,1,1)
        gridLayout_Calculo.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Pressure drop")),3,0,1,1)
        self.DeltaP=Entrada_con_unidades(unidades.Pressure)
        self.DeltaP.valueChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.DeltaP,3,1,1,1)
        gridLayout_Calculo.addItem(QtWidgets.QSpacerItem(10,10,QtWidgets.QSizePolicy.Policy.Fixed,QtWidgets.QSizePolicy.Policy.Fixed),4,0,1,5)
        lyt=QtWidgets.QHBoxLayout()
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Type")))
        self.tipo=QtWidgets.QComboBox()
        self.tipo.addItem(QtWidgets.QApplication.translate("pychemqt", "CSTR, continuous stirred-tank"))
        self.tipo.addItem(QtWidgets.QApplication.translate("pychemqt", "PFR, plug flow"))
        self.tipo.currentIndexChanged.connect(self.tipoCambiado)
        lyt.addWidget(self.tipo)
        lyt.addItem(QtWidgets.QSpacerItem(10,10,QtWidgets.QSizePolicy.Policy.Expanding,QtWidgets.QSizePolicy.Policy.Fixed))
        gridLayout_Calculo.addLayout(lyt,5,0,1,5)
        gridLayout_Calculo.addItem(QtWidgets.QSpacerItem(10,10,QtWidgets.QSizePolicy.Policy.Fixed,QtWidgets.QSizePolicy.Policy.Fixed),6,0,1,5)

        groupbox=QtWidgets.QGroupBox(QtWidgets.QApplication.translate("pychemqt", "Thermal mode"))
        layout=QtWidgets.QGridLayout(groupbox)
        self.checkAdiabatico=QtWidgets.QRadioButton(QtWidgets.QApplication.translate("pychemqt", "Adiabatic"))
        self.checkAdiabatico.toggled.connect(self.heatChanged)
        layout.addWidget(self.checkAdiabatico, 1, 1, 1, 1)
        self.checkIsotermico=QtWidgets.QRadioButton(QtWidgets.QApplication.translate("pychemqt", "Isothermal"))
        self.checkIsotermico.toggled.connect(self.heatChanged)
        layout.addWidget(self.checkIsotermico, 2, 1, 1, 1)
        self.checkFlux=QtWidgets.QRadioButton(QtWidgets.QApplication.translate("pychemqt", "Heat duty"))
        self.checkFlux.toggled.connect(self.heatChanged)
        layout.addWidget(self.checkFlux, 3, 1, 1, 1)
        self.checkIntercambio=QtWidgets.QRadioButton(QtWidgets.QApplication.translate("pychemqt", "Heat transfer"))
        self.checkIntercambio.toggled.connect(self.heatChanged)
        layout.addWidget(self.checkIntercambio, 4, 1, 1, 1)
        self.checkPerfil=QtWidgets.QRadioButton(QtWidgets.QApplication.translate("pychemqt", "PFR temperature profile"))
        self.checkPerfil.toggled.connect(self.heatChanged)
        layout.addWidget(self.checkPerfil, 5, 1, 1, 1)
        self.T=Entrada_con_unidades(unidades.Temperature)
        self.T.valueChanged.connect(self.calculo)
        layout.addWidget(self.T, 2, 2, 1, 2)
        self.Q=Entrada_con_unidades(unidades.Power)
        self.Q.valueChanged.connect(self.calculo)
        layout.addWidget(self.Q, 3, 2, 1, 2)
        self.T_ext=Entrada_con_unidades(unidades.Temperature)
        self.T_ext.valueChanged.connect(self.calculo)
        layout.addWidget(self.T_ext,4,2,1,2)
        self.U=Entrada_con_unidades(unidades.HeatTransfCoef)
        self.U.valueChanged.connect(self.calculo)
        layout.addWidget(self.U,4,4)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Flow")),4,5)
        self.direccion=QtWidgets.QComboBox()
        self.direccion.addItem(QtWidgets.QApplication.translate("pychemqt", "Countercurrent"))
        self.direccion.addItem(QtWidgets.QApplication.translate("pychemqt", "Cocurrent"))
        layout.addWidget(self.direccion,4,6)
        self.botonPerfil=QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/table.png")), QtWidgets.QApplication.translate("pychemqt", "Add Profile"))
        self.botonPerfil.clicked.connect(self.editorPerfil)
        layout.addWidget(self.botonPerfil,5,2,1,1)
        gridLayout_Calculo.addWidget(groupbox, 7, 0, 1, 5)

        self.groupBox_Diseno= QtWidgets.QGroupBox(QtWidgets.QApplication.translate("pychemqt", "Design"))
        gridLayout_Calculo.addWidget(self.groupBox_Diseno,8,0,1,5)
        lyt = QtWidgets.QGridLayout(self.groupBox_Diseno)
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Mode")),1,1)
        self.modo=QtWidgets.QComboBox()
        self.modo.addItem(QtWidgets.QApplication.translate("pychemqt", "Rating: calculate conversión"))
        self.modo.addItem(QtWidgets.QApplication.translate("pychemqt", "Design, calculate volumen"))
        self.modo.currentIndexChanged.connect(self.calculo)
        lyt.addWidget(self.modo,1,2,1,3)
        lyt.addItem(QtWidgets.QSpacerItem(10,10,QtWidgets.QSizePolicy.Policy.Expanding,QtWidgets.QSizePolicy.Policy.Fixed),1,5)
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Reactor Volume")),2,1)
        self.V=Entrada_con_unidades(unidades.Volume, "VolLiq")
        lyt.addWidget(self.V,2,2)
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Key Component")),3,1)
        self.key=QtWidgets.QComboBox()
#        for i, nombre in enumerate(self.nombres):
#            self.key.addItem("%i - %s" %(i+1, nombre))
        lyt.addWidget(self.key,3,2)
        lyt.addItem(QtWidgets.QSpacerItem(20,10,QtWidgets.QSizePolicy.Policy.Fixed,QtWidgets.QSizePolicy.Policy.Fixed),3,3)
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Conversion")),3,4)
        self.conversion=Entrada_con_unidades(float, max=1)
        lyt.addWidget(self.conversion,3,5)

        gridLayout_Calculo.addItem(QtWidgets.QSpacerItem(10,10,QtWidgets.QSizePolicy.Policy.Expanding,QtWidgets.QSizePolicy.Policy.Expanding),10,0,1,5)

        groupBox_Calculo = QtWidgets.QGroupBox(QtWidgets.QApplication.translate("pychemqt", "Results"))
        gridLayout_Calculo.addWidget(groupBox_Calculo,11,0,1,5)
        gridLayout_1 = QtWidgets.QGridLayout(groupBox_Calculo)
        gridLayout_1.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "T output")),0,1,1,1)
        self.TCalc=Entrada_con_unidades(unidades.Temperature, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.TCalc,0,2,1,1)
        gridLayout_1.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Heat")),1,1,1,1)
        self.HeatCalc=Entrada_con_unidades(unidades.Power, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.HeatCalc,1,2,1,1)

        gridLayout_1.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Conversion")),0,4)
        self.conversionCalc=Entrada_con_unidades(float, readOnly=True)
        gridLayout_1.addWidget(self.conversionCalc,0,5)
        gridLayout_1.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Reactor Volume")),1,4)
        self.VCalc=Entrada_con_unidades(unidades.Volume, "VolLiq", readOnly=True)
        gridLayout_1.addWidget(self.VCalc,1,5)

        self.checkAdiabatico.setChecked(True)
        self.tipoCambiado(0)


    def heatChanged(self):
        self.T.setReadOnly(not self.checkIsotermico.isChecked())
        self.Q.setReadOnly(not self.checkFlux.isChecked())
        self.T_ext.setReadOnly(not self.checkIntercambio.isChecked())
        self.U.setReadOnly(not self.checkIntercambio.isChecked())
        self.direccion.setEnabled(self.checkIntercambio.isChecked() and self.tipo.currentIndex())
        self.botonPerfil.setEnabled(self.checkPerfil.isChecked())

    def tipoCambiado(self, ind):
        self.checkPerfil.setEnabled(ind)
        self.direccion.setEnabled(self.checkIntercambio.isChecked() and ind)
        self.calculo()


    def editorPerfil(self):
        dialog=inputTable.InputTableDialog(2, data=self.profile_T, title=QtWidgets.QApplication.translate("pychemqt", "Temperature profile"), horizontalHeader=["x", "T, "+unidades.Temperature(None).text()])
        if dialog.exec():
            self.profile_T=dialog.data

#    def cambiar_entrada(self, corriente):
#        self.entrada=corriente
#        self.calculo()
#
#    def todos_datos(self):
#        if self.checkAdiabatico.isChecked():
#            heat=True
#        elif self.checkIsotermico.isChecked():
#            heat=self.T.value
#        elif self.checkFlux.isChecked():
#            heat=self.Q.value
#        else: heat=self.T_ext.value and self.U.value
#        return  self.Reacciones.reacciones and heat

#    def calculo(self):
#        if self.todos_datos():
#            self.status.setState(4)
#            if self.checkAdiabatico.isChecked():
#                thermal=0
#            elif self.checkIsotermico.isChecked():
#                thermal=1
#            elif self.checkFlux.isChecked():
#                thermal=2
#            else: thermal=3
#
#            if self.P.value:
#                P=self.P.value.atm
#            else:
#                P=None
#
#            self.Equipment(entrada=self.entrada, thermal=thermal, reaccion=self.Reacciones.reacciones, P=P, T=self.T.value, Q=self.Q.value, Text=self.T_ext.value, U=self.U.value)
#
#            self.rellenoSalida()
#            self.status.setState(1)
#
#
#    def rellenoSalida(self):
#        self.TCalc.setValue(self.Equipment.Salida.T)
#        self.HeatCalc.setValue(self.Equipment.Heat)
#        self.Salida.rellenar(self.Equipment.Salida)


if __name__ == "__main__":
    import sys
    from lib.corriente import Corriente, Mezcla, Solid
    from numpy import r_
    app = QtWidgets.QApplication(sys.argv)
#    dialogo = UI_reacciones()
    dialogo = UI_equipment()
    dialogo.show()
    sys.exit(app.exec())
