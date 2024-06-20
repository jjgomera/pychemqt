#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2023, Juan José Gómez Romera <jjgomera@gmail.com>

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
# Diálogo de definición de reactores UI_reactor                               #
###############################################################################

from functools import partial
from math import exp, log
import os
import pickle

from numpy import array
from scipy.optimize import leastsq

from lib import unidades, reaction, plot
from lib.config import getComponents, IMAGE_PATH
from lib.thread import Evaluate
from equipment.parents import UI_equip
from equipment.reactor import Reactor
from tools.qt import QtCore, QtGui, QtWidgets
from UI import inputTable
from UI.widgets import Entrada_con_unidades, Tabla, QLabelMath, Status


class ReactionWidget(QtWidgets.QWidget):
    """Show defined reactions and let user configure"""
    changed = QtCore.pyqtSignal()
    reacciones = []
    reaccion = None
    activo = None
    ajuste = None

    def __init__(self, parent=None):
        super().__init__(parent)
        self.indices, self.nombres, M = getComponents()
        gridLayout = QtWidgets.QGridLayout(self)

        htitle = [self.tr("Reaction"),
                  "ΔHr, %s" % unidades.MolarEnthalpy(None).text(),
                  self.tr("Type"),
                  self.tr("Phase"),
                  self.tr("Description")]
        self.Table = Tabla(
            5, horizontalHeader=htitle, dinamica=False, verticalHeader=True,
            orientacion=QtCore.Qt.AlignmentFlag.AlignLeft)
        self.Table.setMinimumWidth(500)
        self.Table.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.Table.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        self.Table.horizontalHeader().setStretchLastSection(True)
        self.Table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.Table.itemSelectionChanged.connect(self.updateUI)
        gridLayout.addWidget(self.Table, 1, 1, 6, 4)

        self.btnOpen = QtWidgets.QPushButton(
            QtGui.QIcon(os.path.join(IMAGE_PATH, "button", "fileOpen.png")),
            self.tr("Open"))
        self.btnOpen.clicked.connect(self.btnOpenClicked)
        gridLayout.addWidget(self.btnOpen, 1, 5)
        self.btnSave = QtWidgets.QPushButton(
            QtGui.QIcon(os.path.join(IMAGE_PATH, "button", "fileSave.png")),
            self.tr("Save"))
        self.btnSave.clicked.connect(self.btnSaveClicked)
        self.btnSave.setSizePolicy(QtWidgets.QSizePolicy.Policy.Fixed,
                                   QtWidgets.QSizePolicy.Policy.Fixed)
        self.btnSave.setEnabled(False)
        gridLayout.addWidget(self.btnSave, 2, 5)

        self.btnNew = QtWidgets.QPushButton(
            QtGui.QIcon(os.path.join(IMAGE_PATH, "button", "fileNew.png")),
            self.tr("New"))
        self.btnNew.clicked.connect(self.btnNewClicked)
        gridLayout.addWidget(self.btnNew, 3, 5)
        self.btnEdit = QtWidgets.QPushButton(
            QtGui.QIcon(os.path.join(IMAGE_PATH, "button", "editor.png")),
            self.tr("Edit"))
        self.btnEdit.setEnabled(False)
        self.btnEdit.setCheckable(True)
        self.btnEdit.clicked.connect(self.btnEditClicked)
        gridLayout.addWidget(self.btnEdit, 4, 5)
        self.btnDel = QtWidgets.QPushButton(
            QtGui.QIcon(os.path.join(IMAGE_PATH, "button", "editDelete.png")),
            self.tr("Delete"))
        self.btnDel.setEnabled(False)
        self.btnDel.clicked.connect(self.btnDelClicked)
        gridLayout.addWidget(self.btnDel, 5, 5)
        self.btnClear = QtWidgets.QPushButton(
            QtGui.QIcon(os.path.join(IMAGE_PATH, "button", "clear.png")),
            self.tr("Clear"))
        self.btnClear.clicked.connect(self.btnClearClicked)
        gridLayout.addWidget(self.btnClear, 6, 5)
        gridLayout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 10, 1)

    def updateUI(self, bool=True):
        """Update button state"""
        self.btnEdit.setEnabled(bool)
        self.btnDel.setEnabled(bool)

    def btnOpenClicked(self):
        patron = self.tr("reaction file")
        patron += " (*.rec)"
        title = self.tr("Open reaction file")
        fname = QtWidgets.QFileDialog.getOpenFileName(
            self, title, "./", patron)[0]
        if fname:
            with open(fname, "rb") as archivo:
                reacciones = pickle.load(archivo)

            self.reacciones = reacciones
            self.btnSave.setEnabled(True)
            for fila, reaccion in enumerate(reacciones):
                reactiontxt = "%i - %s" % (
                    reaccion.tipo+1,
                    reaction.Reaction.TEXT_TYPE[reaccion.tipo])
                self.Table.addRow()
                self.Table.setValue(fila, 0, reaccion.text)
                self.Table.setValue(fila, 1, "%0.4e" % reaccion.Hr.config())
                self.Table.setValue(fila, 2, reactiontxt)
                self.Table.setValue(
                    fila, 3, reaction.Reaction.TEXT_PHASE[reaccion.fase])
                self.Table.item(fila, 4).setFlags(
                    QtCore.Qt.ItemFlag.ItemIsEditable
                    | QtCore.Qt.ItemFlag.ItemIsEnabled
                    | QtCore.Qt.ItemFlag.ItemIsSelectable)
            for i in range(4):
                self.Table.resizeColumnToContents(i)
        self.changed.emit()

    def btnSaveClicked(self):
        title = self.tr("Save reaction to file")
        patron = self.tr("reaction file")
        patron += " (*.rec)"
        fname = QtWidgets.QFileDialog.getSaveFileName(
            self, title, "./", patron)[0]
        if fname:
            if fname.split(".")[-1] != "rec":
                fname += ".rec"
            pickle.dump(self.reacciones, open(fname, "w"))

    def btnNewClicked(self):
        dialog = UI_reacciones(parent=self)
        if dialog.exec():
            # TODO: Add new reaction
            pass

    def btnEditClicked(self, bool):
        if bool:
            indice = self.Table.currentRow()
            reaccion = self.reacciones[indice]
            dialog = UI_reacciones(reaccion, self)
            dialog.exec()
            # self.rellenar(self.reaccion)
            # self.activo=indice
        else:
            self.botonAddClicked(self.activo, False)
            self.reacciones[self.activo] = self.reaccion
            self.Table.setCurrentCell(self.activo, 0)
            self.activo = -1
            self.changed.emit()

        self.btnNew.setEnabled(not bool)
        self.btnDel.setEnabled(not bool)
        self.btnClear.setEnabled(not bool)
        self.botonAdd.setEnabled(not bool)
        self.btnOpen.setEnabled(not bool)
        self.btnSave.setEnabled(not bool)

    def btnDelClicked(self):
        indice = self.Table.currentRow()
        self.Table.removeRow(indice)
        del self.reacciones[indice]
        self.Table.clearSelection()
        self.updateUI(False)
        self.changed.emit()

    def btnClearClicked(self):
        if self.reacciones:
            self.reacciones = []
            self.Table.setRowCount(0)
            self.btnSave.setEnabled(False)

    def botonAddClicked(self, fila, add=True):
        if add:
            fila = self.Table.rowCount()
            self.Table.addRow()
        self.Table.setValue(fila, 0, self.Formula.text())
        self.Table.setValue(fila, 1, "%0.4e" % self.Hr.value.config())
        self.Table.setValue(fila, 2, "%i - %s" % (
            self.tipo.currentIndex()+1, self.tipo.currentText()))
        self.Table.setValue(fila, 3, self.Fase.currentText())
        self.Table.item(fila, 4).setFlags(
            QtCore.Qt.ItemFlag.ItemIsEditable
            | QtCore.Qt.ItemFlag.ItemIsEnabled
            | QtCore.Qt.ItemFlag.ItemIsSelectable)
        for i in range(4):
            self.Table.resizeColumnToContents(i)
        self.reacciones.insert(fila, self.reaccion)
        self.btnSave.setEnabled(True)
        self.changed.emit()


class UI_reacciones(QtWidgets.QDialog):
    reaction = reaction.Reaction()

    def __init__(self, reaccion=None, parent=None):
        super().__init__(parent)
        self.evaluate = Evaluate()
        self.evaluate.finished.connect(self.rellenar)
        self.indices, self.nombres, M = getComponents()
        gridLayout = QtWidgets.QGridLayout(self)

        lyt = QtWidgets.QHBoxLayout()
        lyt.addWidget(QtWidgets.QLabel(
            self.tr("Key component")))
        self.key = QtWidgets.QComboBox()
        for i, nombre in enumerate(self.nombres):
            self.key.addItem("%i - %s" % (i+1, nombre))
        self.key.currentIndexChanged.connect(partial(self.changeParams, "key"))
        lyt.addWidget(self.key)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed))
        gridLayout.addLayout(lyt, 1, 1, 1, 5)

        lyt = QtWidgets.QHBoxLayout()
        lyt.addWidget(QtWidgets.QLabel(
            self.tr("Phase")))
        self.fase = QtWidgets.QComboBox()
        for txt in reaction.Reaction.TEXT_PHASE:
            self.fase.addItem(txt)
        self.fase.currentIndexChanged.connect(
            partial(self.changeParams, "fase"))
        lyt.addWidget(self.fase)
        self.Formula = QtWidgets.QLabel()
        self.Formula.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.Formula.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding,
                                   QtWidgets.QSizePolicy.Policy.Fixed)
        lyt.addWidget(self.Formula)
        gridLayout.addLayout(lyt, 2, 1, 1, 5)

        lyt = QtWidgets.QVBoxLayout()
        title = self.nombres[:]
        title.append("")
        self.Estequiometria = Tabla(
            1, verticalHeaderLabels=title,
            horizontalHeader=[self.tr("Coefficients")],
            filas=len(self.indices))
        self.Estequiometria.setFixedHeight(22*len(self.indices)+22+4+22)
        lyt.addWidget(self.Estequiometria)
        self.Estequiometria.addRow()
        brush = QtGui.QBrush(QtGui.QColor("#eaeaea"))
        self.Estequiometria.item(len(self.indices), 0).setBackground(brush)
        self.Estequiometria.item(len(self.indices), 0).setFlags(
            QtCore.Qt.ItemFlag.NoItemFlags)
        self.Estequiometria.cellChanged.connect(self.reaccionCambiada)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Expanding))
        gridLayout.addLayout(lyt, 3, 1, 1, 2)

        lyt = QtWidgets.QGridLayout()
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 1, 1)
        self.formula = QtWidgets.QCheckBox(
            self.tr("Use name in formula"))
        self.formula.toggled.connect(partial(self.changeParams, "formula"))
        lyt.addWidget(self.formula, 1, 2, 1, 2)
        self.customHr = QtWidgets.QCheckBox("ΔHr " + self.tr("user specified"))
        self.customHr.toggled.connect(self.changeHr)
        lyt.addWidget(self.customHr, 2, 2, 1, 2)
        lyt.addWidget(QtWidgets.QLabel("ΔHr<sup>o</sup>"), 3, 2)
        self.Hr = Entrada_con_unidades(unidades.MolarEnthalpy, readOnly=True)
        self.Hr.valueChanged.connect(partial(self.changeParams, "Hr"))
        lyt.addWidget(self.Hr, 3, 3)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Expanding))
        gridLayout.addLayout(lyt, 3, 3, 1, 2)

        gridLayout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 4, 2)

        lyt = QtWidgets.QHBoxLayout()
        lyt.addWidget(QtWidgets.QLabel(self.tr("Type")))
        self.tipo = QtWidgets.QComboBox()
        for txt in reaction.Reaction.TEXT_TYPE:
            self.tipo.addItem(txt)
        self.tipo.currentIndexChanged.connect(
            partial(self.changeParams, "tipo"))
        lyt.addWidget(self.tipo)
        lyt.addItem(QtWidgets.QSpacerItem(
            20, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed))
        lyt.addWidget(QtWidgets.QLabel(self.tr("Concentration")))
        self.base = QtWidgets.QComboBox()
        for txt in reaction.Reaction.TEXT_BASE:
            self.base.addItem(txt)
        self.base.currentIndexChanged.connect(
            partial(self.changeParams, "base"))
        lyt.addWidget(self.base)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed))
        gridLayout.addLayout(lyt, 5, 1, 1, 5)

        self.stacked = QtWidgets.QStackedWidget()
        self.tipo.currentIndexChanged.connect(self.stacked.setCurrentIndex)
        gridLayout.addWidget(self.stacked, 6, 1, 1, 5)
        gridLayout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 7, 1, 1, 5)

        widget = QtWidgets.QWidget()
        self.stacked.addWidget(widget)
        lyt = QtWidgets.QGridLayout(widget)
        lyt.addWidget(QtWidgets.QLabel(
            "<h3>"+self.tr("Estequiometric reaction")+"</h3>"),1,1,1,4)
        self.Conversion = Tabla(1, verticalHeaderModel="C", filas=3)
        self.Conversion.setConnected()
        self.Conversion.setFixedWidth(100)
        lyt.addWidget(self.Conversion, 2, 1, 3, 1)
        mathTex = r"$Conversion = C_o + C_1T + C_2T^2 + \cdots + C_nT^n$"
        label = QLabelMath(mathTex)
        lyt.addWidget(label, 2, 2, 1, 3)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Temperature unit")),3,2)
        self.unidadesTemperatura = QtWidgets.QComboBox()
        for i in unidades.Temperature.__text__:
            self.unidadesTemperatura.addItem(i)
        lyt.addWidget(self.unidadesTemperatura, 3, 3)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 4, 4)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 5, 1, 1, 5)

        widget = QtWidgets.QWidget()
        self.stacked.addWidget(widget)
        lyt = QtWidgets.QGridLayout(widget)
        self.check_KFijo = QtWidgets.QRadioButton(self.tr("Fixed"))
        self.check_KFijo.toggled.connect(self.KeqChanged)
        lyt.addWidget(self.check_KFijo, 1, 1, 1, 2)
        lyt.addWidget(QtWidgets.QLabel("K<sub>eq</sub>"), 1, 3)
        self.Keq = Entrada_con_unidades(float)
        lyt.addWidget(self.Keq, 1, 4)
        mathTex = r"$aA + bB \rightleftharpoons cC + dD \therefore "
        mathTex += r"K_{eq} = \frac{[C]^c [D]^d}{[A]^a [B]^b}$"
        label = QLabelMath(mathTex)
        lyt.addWidget(label, 1, 5, 1, 4)

        self.check_KEq = QtWidgets.QRadioButton(self.tr("Equation"))
        self.check_KEq.toggled.connect(self.KeqChanged)
        lyt.addWidget(self.check_KEq, 2, 1, 1, 2)
        self.check_KTabla = QtWidgets.QRadioButton(self.tr("Table"))
        self.check_KTabla.toggled.connect(self.KeqChanged)
        lyt.addWidget(self.check_KTabla, 2, 5, 1, 2)
        self.KEq_Dat = Tabla(1, verticalHeaderLabels=["A", "B", "C", "D", "E", "F", "G", "H"], filas=8)
        self.KEq_Dat.setFixedHeight(22*8+4)
        self.KEq_Dat.setFixedWidth(120)
        lyt.addWidget(self.KEq_Dat, 3, 3, 1, 2)
        self.KEq_Tab = Tabla(4, horizontalHeader=["T, K", "Keq", "Kcalc", "%Error"], verticalHeader=False, columnReadOnly=[False, False, True, True])
        self.KEq_Tab.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding,
                                   QtWidgets.QSizePolicy.Policy.Expanding)
        self.KEq_Tab.setFixedWidth(400)
        self.KEq_Tab.setConnected()
        self.KEq_Tab.rowFinished.connect(self.Regresion)
        self.KEq_Tab.setAlternatingRowColors(False)
        lyt.addWidget(self.KEq_Tab, 3, 5, 1, 4)
        lyt.addWidget(QtWidgets.QLabel("r²"), 4, 5)
        self.r2 = Entrada_con_unidades(float, readOnly=True)
        lyt.addWidget(self.r2, 4, 6)
        self.botonTablaPlot = QtWidgets.QPushButton(
            QtGui.QIcon(os.path.join(IMAGE_PATH, "button", "plot.png")),
            self.tr("Plot"))
        self.botonTablaPlot.clicked.connect(self.Plot)
        lyt.addWidget(self.botonTablaPlot, 4, 7)
        self.botonTablaClear = QtWidgets.QPushButton(
            QtGui.QIcon(os.path.join(IMAGE_PATH, "button", "clear.png")),
            self.tr("Clear"))
        self.botonTablaClear.clicked.connect(self.KEq_Tab.clear)
        lyt.addWidget(self.botonTablaClear, 4, 8)
        mathTex = r"$\ln K_eq = A+B/T+C\ln T+DT+ET^2+FT^3+GT^4+HT^5$"
        label = QLabelMath(mathTex)
        lyt.addWidget(label, 5, 1, 1, 8)

        self.checkGibbs = QtWidgets.QRadioButton(
            self.tr("From Gibbs free energy minimization"))
        lyt.addWidget(self.checkGibbs, 6, 1, 1, 4)
        self.check_KFijo.setChecked(True)

        widget = QtWidgets.QWidget()
        self.stacked.addWidget(widget)
        lyt = QtWidgets.QGridLayout(widget)

        widget = QtWidgets.QWidget()
        self.stacked.addWidget(widget)
        lyt = QtWidgets.QGridLayout(widget)

        self.status = Status()
        gridLayout.addWidget(self.status, 10, 1)
        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Cancel
            | QtWidgets.QDialogButtonBox.StandardButton.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        gridLayout.addWidget(self.buttonBox, 10, 2, 1, 4)

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
        kwargs = {"componentes": self.indices,
                  "coeficientes": self.Estequiometria.getColumn(0)[:-1]}
        self.calculo(**kwargs)

    def setReaction(self, reaction):
        self.reaction=reaction
        self.rellenar()

        # if self.Estequiometria.getValue(0, self.Base.currentIndex()):
            # reaccion = reaction.Reaction(
                # self.indices, self.Estequiometria.getColumn(0),
                # base=self.Base.currentIndex(), estequiometria=[0, 0, 0.5],
                # formulas=self.checkFormula.isChecked(),
                # calor=self.checkCalorEspecificado.isChecked(),
                # Hr=self.Hr.value, tipo=self.tipo.currentIndex(),
                # conversion=self.Conversion.getColumn(0)[-1::-1])
            # self.Balance.setValue(reaccion.error)
            # if reaccion.state:
                # self.Formula.setText(
                    # reaccion._txt(self.checkFormula.isChecked()))
                # self.Hr.setValue(reaccion.Hr)
            # else:
                # self.Formula.clear()
                # self.Hr.clear()
            # self.botonAdd.setEnabled(
                # reaccion.state and not self.btnEdit.isChecked())
            # self.reaccion=reaccion

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
        self.Estequiometria.item(len(self.indices), 0).setText(
            str(self.reaction.error))
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
        t = array(self.KEq_Tab.getColumn(0)[:-1])
        k = array(self.KEq_Tab.getColumn(1)[:-1])
        if len(t) >= 4:
            if 4 <= len(t) < 8:
                inicio = r_[0, 0, 0, 0]
                f = lambda par, T: exp(par[0]+par[1]/T+par[2]*log(T)+par[3]*T)
                resto = lambda par, T, k: k-f(par, T)
            else:
                inicio = r_[0, 0, 0, 0, 0, 0, 0, 0]
                f = lambda par, T: exp(par[0]+par[1]/T+par[2]*log(T)+par[3]*T+par[4]*T**2+par[5]*T**3+par[6]*T**4+par[7]*T**5)
                resto = lambda par, T, k: k-f(par, T)

            ajuste = leastsq(resto,inicio,args=(t, k))
            kcalc = f(ajuste[0], t)
            error = (k-kcalc)/k*100
            self.KEq_Dat.setColumn(0, ajuste[0])
            self.KEq_Tab.setColumn(2, kcalc)
            self.KEq_Tab.setColumn(3, error)

            if ajuste[1] in [1, 2, 3, 4]:
                self.ajuste = ajuste[0]

    def Plot(self):
        if self.ajuste != None:
            t = array(self.KEq_Tab.getColumn(0)[:-1])
            k = array(self.KEq_Tab.getColumn(1)[:-1])
            if 4 <= len(t) < 8:
                f = lambda par, T: exp(par[0]+par[1]/T+par[2]*log(T)+par[3]*T)
            else:
                f = lambda par, T: exp(par[0]+par[1]/T+par[2]*log(T)+par[3]*T+par[4]*T**2+par[5]*T**3+par[6]*T**4+par[7]*T**5)

            grafico = plot.PlotDialog()
            grafico.data(t, k, 'ro', t, f(self.ajuste, t))
            grafico.exec()


class UI_equipment(UI_equip):
    """Diálogo de definición de reactores"""
    profile_T = None
    Equipment = Reactor()
    def __init__(self, equipment=None, parent=None):
        """entrada: Parametro opcional de clase corriente que indica la
        corriente de entrada en el reactor"""
        super().__init__(Reactor, entrada=False, salida=False, parent=parent)

        # Reactions tab
        self.Reacciones = ReactionWidget()
        self.Reacciones.changed.connect(self.calculo)
        self.tabWidget.insertTab(1, self.Reacciones, self.tr("Reactions"))

        # Calculate tab
        gridLayout_Calculo = QtWidgets.QGridLayout(self.tabCalculo)
        gridLayout_Calculo.addWidget(QtWidgets.QLabel(
            self.tr("Pressure")), 2, 0, 1, 1)
        self.P = Entrada_con_unidades(unidades.Pressure)
        self.P.valueChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.P, 2, 1, 1, 1)
        gridLayout_Calculo.addWidget(QtWidgets.QLabel(
            self.tr("Pressure drop")), 3, 0, 1, 1)
        self.DeltaP = Entrada_con_unidades(unidades.Pressure)
        self.DeltaP.valueChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.DeltaP, 3, 1, 1, 1)
        gridLayout_Calculo.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 4, 0, 1, 5)
        lyt = QtWidgets.QHBoxLayout()
        lyt.addWidget(QtWidgets.QLabel(self.tr("Type")))
        self.tipo = QtWidgets.QComboBox()
        self.tipo.addItem(self.tr("CSTR, continuous stirred-tank"))
        self.tipo.addItem(self.tr("PFR, plug flow"))
        self.tipo.currentIndexChanged.connect(self.tipoCambiado)
        lyt.addWidget(self.tipo)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed))
        gridLayout_Calculo.addLayout(lyt, 5, 0, 1, 5)
        gridLayout_Calculo.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 6, 0, 1, 5)

        groupbox = QtWidgets.QGroupBox(self.tr("Thermal mode"))
        layout = QtWidgets.QGridLayout(groupbox)
        self.checkAdiabatico = QtWidgets.QRadioButton(self.tr("Adiabatic"))
        self.checkAdiabatico.toggled.connect(self.heatChanged)
        layout.addWidget(self.checkAdiabatico, 1, 1, 1, 1)
        self.checkIsotermico = QtWidgets.QRadioButton(self.tr("Isothermal"))
        self.checkIsotermico.toggled.connect(self.heatChanged)
        layout.addWidget(self.checkIsotermico, 2, 1, 1, 1)
        self.checkFlux = QtWidgets.QRadioButton(self.tr("Heat duty"))
        self.checkFlux.toggled.connect(self.heatChanged)
        layout.addWidget(self.checkFlux, 3, 1, 1, 1)
        self.checkIntercambio = QtWidgets.QRadioButton(
            self.tr("Heat transfer"))
        self.checkIntercambio.toggled.connect(self.heatChanged)
        layout.addWidget(self.checkIntercambio, 4, 1, 1, 1)
        self.checkPerfil = QtWidgets.QRadioButton(
            self.tr("PFR temperature profile"))
        self.checkPerfil.toggled.connect(self.heatChanged)
        layout.addWidget(self.checkPerfil, 5, 1, 1, 1)
        self.T = Entrada_con_unidades(unidades.Temperature)
        self.T.valueChanged.connect(self.calculo)
        layout.addWidget(self.T, 2, 2, 1, 2)
        self.Q = Entrada_con_unidades(unidades.Power)
        self.Q.valueChanged.connect(self.calculo)
        layout.addWidget(self.Q, 3, 2, 1, 2)
        self.T_ext = Entrada_con_unidades(unidades.Temperature)
        self.T_ext.valueChanged.connect(self.calculo)
        layout.addWidget(self.T_ext, 4, 2, 1, 2)
        self.U = Entrada_con_unidades(unidades.HeatTransfCoef)
        self.U.valueChanged.connect(self.calculo)
        layout.addWidget(self.U, 4, 4)
        layout.addWidget(QtWidgets.QLabel(self.tr("Flow")), 4, 5)
        self.direccion = QtWidgets.QComboBox()
        self.direccion.addItem(self.tr("Countercurrent"))
        self.direccion.addItem(self.tr("Cocurrent"))
        layout.addWidget(self.direccion, 4, 6)
        self.botonPerfil = QtWidgets.QPushButton(
            QtGui.QIcon(os.path.join(IMAGE_PATH, "button", "table.png")),
            self.tr("Add Profile"))
        self.botonPerfil.clicked.connect(self.editorPerfil)
        layout.addWidget(self.botonPerfil, 5, 2, 1, 1)
        gridLayout_Calculo.addWidget(groupbox, 7, 0, 1, 5)

        self.groupBox_Diseno = QtWidgets.QGroupBox(self.tr("Design"))
        gridLayout_Calculo.addWidget(self.groupBox_Diseno, 8, 0, 1, 5)
        lyt = QtWidgets.QGridLayout(self.groupBox_Diseno)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Mode")), 1, 1)
        self.modo = QtWidgets.QComboBox()
        self.modo.addItem(self.tr("Rating: calculate conversión"))
        self.modo.addItem(self.tr("Design, calculate volumen"))
        self.modo.currentIndexChanged.connect(self.calculo)
        lyt.addWidget(self.modo, 1, 2, 1, 3)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed), 1, 5)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Reactor Volume")), 2, 1)
        self.V = Entrada_con_unidades(unidades.Volume, "VolLiq")
        lyt.addWidget(self.V, 2, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Key Component")), 3, 1)
        self.key = QtWidgets.QComboBox()
        for i, nombre in enumerate(self.Reacciones.nombres):
            self.key.addItem("%i - %s" % (i+1, nombre))
        lyt.addWidget(self.key, 3, 2)
        lyt.addItem(QtWidgets.QSpacerItem(
            20, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 3, 3)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Conversion")), 3, 4)
        self.conversion = Entrada_con_unidades(float, max=1)
        lyt.addWidget(self.conversion, 3, 5)

        gridLayout_Calculo.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 10, 0, 1, 5)

        groupBox_Calculo = QtWidgets.QGroupBox(self.tr("Results"))
        gridLayout_Calculo.addWidget(groupBox_Calculo, 11, 0, 1, 5)
        gridLayout_1 = QtWidgets.QGridLayout(groupBox_Calculo)
        gridLayout_1.addWidget(QtWidgets.QLabel(self.tr("T output")),
            0, 1, 1, 1)
        self.TCalc = Entrada_con_unidades(
            unidades.Temperature, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.TCalc, 0, 2, 1, 1)
        gridLayout_1.addWidget(QtWidgets.QLabel(self.tr("Heat")), 1, 1, 1, 1)
        self.HeatCalc = Entrada_con_unidades(
            unidades.Power, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.HeatCalc, 1, 2, 1, 1)

        gridLayout_1.addWidget(QtWidgets.QLabel(self.tr("Conversion")), 0, 4)
        self.conversionCalc = Entrada_con_unidades(float, readOnly=True)
        gridLayout_1.addWidget(self.conversionCalc, 0, 5)
        gridLayout_1.addWidget(QtWidgets.QLabel(
            self.tr("Reactor Volume")), 1, 4)
        self.VCalc = Entrada_con_unidades(
            unidades.Volume, "VolLiq", readOnly=True)
        gridLayout_1.addWidget(self.VCalc, 1, 5)

        self.checkAdiabatico.setChecked(True)
        self.tipoCambiado(0)

    def heatChanged(self):
        self.T.setReadOnly(not self.checkIsotermico.isChecked())
        self.Q.setReadOnly(not self.checkFlux.isChecked())
        self.T_ext.setReadOnly(not self.checkIntercambio.isChecked())
        self.U.setReadOnly(not self.checkIntercambio.isChecked())
        self.direccion.setEnabled(self.checkIntercambio.isChecked()
                                  and self.tipo.currentIndex())
        self.botonPerfil.setEnabled(self.checkPerfil.isChecked())

    def tipoCambiado(self, ind):
        self.checkPerfil.setEnabled(ind)
        self.direccion.setEnabled(self.checkIntercambio.isChecked() and ind)
        self.calculo()

    def editorPerfil(self):
        dialog = inputTable.InputTableDialog(
            2, data=self.profile_T,
            title=self.tr("Temperature profile"),
            horizontalHeader=["x", "T, "+unidades.Temperature(None).text()])
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
    dialogo = UI_equipment()
    # dialogo = UI_reacciones()
    # dialogo = ReactionWidget()
    dialogo.show()
    sys.exit (app.exec())
