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
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


###############################################################################
# Pipe equipment dialog and functionality

# - Dialog: Generic dialog for definition of special fitting
#
###############################################################################

from functools import partial
from math import pi
import os

from equipment.parents import UI_equip
from equipment.pipe import Pipe
from lib.config import Preferences, IMAGE_PATH
from lib.utilities import representacion
from lib.pipeDatabase import CATALOG, CATALOG_TRANSLATE, FITTING, FITTING_DESC
from lib.fitting import (K_contraction_crane, K_enlargement_crane, K_flush,
                         K_MitreBend, Ft, K_longBend)
from lib.unidades import (Length, Temperature, HeatTransfCoef, Pressure, Power,
                          Speed, Currency)
from tools.costIndex import CostData
from tools.qt import QtCore, QtGui, QtWidgets
from UI.delegate import SpinEditor, CellEditor
from UI.widgets import Entrada_con_unidades


class Dialog(QtWidgets.QDialog):
    """Generic dialog for definition of special fitting

    Parameters
    ----------
    tipo : int
        Define the kind of dialog
        0 - Contraction, enlargement
        1 - Codo segmentado de angulo arbitrario
        2 - Codo largo
        3 - Entrada redondeada
    titulo : str
        Title of dialog
    icon : str
        Path of icon file
    textos : List
    """

    def __init__(self, tipo=0, titulo="", icon=None, textos=None, parent=None):
        super().__init__(parent)

        if textos is None:
            textos = ("", "")

        self.tipo = tipo
        self.setWindowTitle(titulo)
        self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(icon)))
        gridLayout = QtWidgets.QGridLayout(self)
        self.label = QtWidgets.QLabel(textos[0])
        self.label.setAlignment(
            QtCore.Qt.AlignmentFlag.AlignRight
            | QtCore.Qt.AlignmentFlag.AlignTrailing
            | QtCore.Qt.AlignmentFlag.AlignVCenter)
        gridLayout.addWidget(self.label, 0, 0)
        self.D1 = Entrada_con_unidades(Length)
        self.D1.valueChanged.connect(self.calculaK)
        gridLayout.addWidget(self.D1, 0, 1, 1, 2)
        lb_2 = QtWidgets.QLabel(textos[1])
        lb_2.setAlignment(
            QtCore.Qt.AlignmentFlag.AlignRight
            | QtCore.Qt.AlignmentFlag.AlignTrailing
            | QtCore.Qt.AlignmentFlag.AlignVCenter)
        gridLayout.addWidget(lb_2, 1, 0)
        self.D2 = Entrada_con_unidades(Length)
        self.D2.valueChanged.connect(self.calculaK)
        gridLayout.addWidget(self.D2, 1, 1, 1, 2)
        self.checkBox = QtWidgets.QCheckBox()
        self.checkBox.setText(self.tr("Gradual"))
        self.checkBox.toggled.connect(self.check_toggled)
        self.checkBox.toggled.connect(self.calculaK)
        gridLayout.addWidget(self.checkBox, 2, 0)
        lb_3 = QtWidgets.QLabel(self.tr("Angle"))
        lb_3.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight
                          | QtCore.Qt.AlignmentFlag.AlignTrailing
                          | QtCore.Qt.AlignmentFlag.AlignVCenter)
        gridLayout.addWidget(lb_3, 2, 1)
        self.angulo = QtWidgets.QSpinBox()
        self.angulo.setFixedWidth(50)
        self.angulo.setAlignment(
            QtCore.Qt.AlignmentFlag.AlignRight
            | QtCore.Qt.AlignmentFlag.AlignTrailing
            | QtCore.Qt.AlignmentFlag.AlignVCenter)
        self.angulo.setRange(1, 90)
        gridLayout.addWidget(self.angulo, 2, 2)
        self.regla = QtWidgets.QSlider()
        self.regla.setMaximum(90)
        self.regla.setOrientation(QtCore.Qt.Orientation.Horizontal)
        self.regla.setTickPosition(QtWidgets.QSlider.TickPosition.TicksAbove)
        self.regla.setTickInterval(5)
        self.regla.valueChanged.connect(self.angulo.setValue)
        self.regla.valueChanged.connect(self.calculaK)
        self.angulo.valueChanged.connect(self.regla.setValue)
        self.angulo.valueChanged.connect(self.calculaK)
        gridLayout.addWidget(self.regla, 3, 0, 1, 3)
        self.K_text = QtWidgets.QLabel("K=")
        self.K_text.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter
                                 | QtCore.Qt.AlignmentFlag.AlignVCenter)
        gridLayout.addWidget(self.K_text, 4, 0, 1, 3)
        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok
            | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        self.buttonBox.setCenterButtons(True)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        gridLayout.addWidget(self.buttonBox, 5, 0, 1, 3)

        self.check_toggled(False)

        if tipo == 1:
            lb_2.setVisible(False)
            self.D2.setVisible(False)
            self.checkBox.setVisible(False)
            self.check_toggled(True)

        if tipo >= 2:
            self.checkBox.setVisible(False)
            lb_3.setVisible(False)
            self.angulo.setVisible(False)
            self.regla.setVisible(False)

    def check_toggled(self, bool):
        self.angulo.setEnabled(bool)
        self.regla.setEnabled(bool)

    def todos_datos(self):
        if self.tipo == 1:
            return True
        else:
            return self.D1.value and self.D2.value

    def calculaK(self):
        if self.tipo == 1:
            tita = self.angulo.value()
            f = Ft(self.D1.value.mm)
            K = K_MitreBend(tita)*f
            self.K_text.setText(f"K = {K:0.3f}")
            self.K = K
        else:
            D1 = self.D1.value
            D2 = self.D2.value
            if D1 and D2:
                if self.checkBox.isChecked():
                    tita = self.angulo.value()
                else:
                    tita = 180
                if self.tipo == 0:
                    if D2 > D1:
                        K = K_enlargement_crane(D1, D2, tita)
                    else:
                        K = K_contraction_crane(D1, D2, tita)
                elif self.tipo == 2:
                    f = Ft(self.D2.value.mm)
                    rD = self.D1.value/self.D2.value
                    K = K_longBend(rD)*f
                else:
                    K = K_flush(D1/D2)
                self.K_text.setText(f"K = {K:0.3f}")
                self.K = K


class Catalogo_Materiales(QtWidgets.QWidget):
    """Diálogo del catalogo de materiales y tamaños de tuberías"""
    valueChanged = QtCore.pyqtSignal(list)

    def __init__(self, parent=None):
        super().__init__(parent)
        gridLayout = QtWidgets.QGridLayout(self)

        self.TablaMaterial = QtWidgets.QTableWidget()
        self.TablaMaterial.verticalHeader().hide()
        self.TablaMaterial.setRowCount(1)
        self.TablaMaterial.setColumnCount(3)
        titles = [self.tr("Material"),
                  self.tr("Class"),
                  self.tr("Roughness") + ", mm"]
        self.TablaMaterial.setHorizontalHeaderLabels(titles)
        self.TablaMaterial.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.NoSelection)
        self.TablaMaterial.setColumnWidth(0, 200)
        self.TablaMaterial.setColumnWidth(1, 80)
        self.TablaMaterial.setColumnWidth(2, 100)
        self.TablaMaterial.setRowHeight(0, 20)
        self.TablaMaterial.setFixedHeight(47)
        self.TablaMaterial.setFixedWidth(384)
        self.TablaMaterial.setItemDelegateForColumn(2, CellEditor(self))
        self.TablaMaterial.itemChanged.connect(self.emitirSignal)
        gridLayout.addWidget(self.TablaMaterial, 1, 0, 1, 2)

        self.TablaDiametro = QtWidgets.QTableWidget()
        self.TablaDiametro.verticalHeader().hide()
        self.TablaDiametro.setRowCount(1)
        self.TablaDiametro.setColumnCount(7)
        titles = [
            self.tr("Nº", None),
            self.tr("D int")+", mm",
            self.tr("Thickness") + ", mm",
            self.tr("D ext")+", mm",
            self.tr("Weight")+", kg/m",
            "V"+", m³/100m", "S"+", m²/100m"]
        self.TablaDiametro.setHorizontalHeaderLabels(titles)
        self.TablaDiametro.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.NoSelection)
        self.TablaDiametro.horizontalHeader().setStretchLastSection(True)
        self.TablaDiametro.setRowHeight(0, 20)
        self.TablaDiametro.setFixedHeight(47)
        self.TablaDiametro.setFixedWidth(704)
        self.TablaDiametro.setItemDelegateForColumn(1, CellEditor(self))
        self.TablaDiametro.setItemDelegateForColumn(2, CellEditor(self))
        self.TablaDiametro.setItemDelegateForColumn(3, CellEditor(self))
        self.TablaDiametro.setItemDelegateForColumn(4, CellEditor(self))
        self.TablaDiametro.itemChanged.connect(self.emitirSignal)
        gridLayout.addWidget(self.TablaDiametro, 2, 0, 1, 2)

        path = os.path.join(IMAGE_PATH, "button", "clear.png")
        self.botonClear = QtWidgets.QPushButton(
            QtGui.QIcon(QtGui.QPixmap(path)), self.tr("Reset"))
        self.botonClear.clicked.connect(self.clear)
        gridLayout.addWidget(self.botonClear, 1, 2)

        path = os.path.join(IMAGE_PATH, "button", "filenew.png")
        self.botonEdit = QtWidgets.QPushButton(
            QtGui.QIcon(QtGui.QPixmap(path)), self.tr("Edit"))
        self.botonEdit.clicked.connect(self.edit)
        self.botonEdit.setCheckable(True)
        gridLayout.addWidget(self.botonEdit, 2, 2)

        self.line = QtWidgets.QFrame()
        self.line.setFrameShape(QtWidgets.QFrame.Shape.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Shadow.Sunken)
        gridLayout.addWidget(self.line, 3, 0, 1, 3)

        label = QtWidgets.QLabel(self.tr("Material:"))
        label.setSizePolicy(QtWidgets.QSizePolicy.Policy.Fixed,
                            QtWidgets.QSizePolicy.Policy.Fixed)
        gridLayout.addWidget(label, 4, 0)
        self.ComboMaterial = QtWidgets.QComboBox()
        gridLayout.addWidget(self.ComboMaterial, 4, 1)

        self.TablaDiametro2 = QtWidgets.QTableWidget()
        self.TablaDiametro2.verticalHeader().hide()
        self.TablaDiametro2.setColumnCount(7)
        self.TablaDiametro2.setHorizontalHeaderLabels(titles)
        self.TablaDiametro2.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.TablaDiametro2.setVerticalScrollBarPolicy(
            QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.TablaDiametro2.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.TablaDiametro2.setFixedWidth(720)
        gridLayout.addWidget(self.TablaDiametro2, 5, 0, 1, 2)

        path = os.path.join(IMAGE_PATH, "button", "arrow-up-double.png")
        self.botonSelect = QtWidgets.QPushButton(
            QtGui.QIcon(QtGui.QPixmap(path)), self.tr("Select"))
        self.botonSelect.setDisabled(True)
        self.TablaDiametro2.currentCellChanged.connect(self.selectedChanged)
        self.TablaDiametro2.cellDoubleClicked.connect(self.add)
        self.botonSelect.clicked.connect(self.add)
        gridLayout.addWidget(self.botonSelect, 4, 2)
        gridLayout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Preferred,
            QtWidgets.QSizePolicy.Policy.Preferred), 6, 0, 1, 3)

        materiales = []
        diametros = []
        for material in CATALOG:
            txt = list(material[1:3])
            if txt not in materiales:
                materiales.append(txt)
                diametros.append([])
            indice = materiales.index(txt)
            diametros[indice].append(material[3:])
        for material in materiales:
            self.ComboMaterial.addItem(material[0]+" "+material[1])
            if material[0] in CATALOG_TRANSLATE:
                material[0] = CATALOG_TRANSLATE[material[0]]
        self.ComboMaterial.currentIndexChanged.connect(self.rellenarDiametros)
        self.diametros = diametros
        self.materiales = materiales
        self.rellenarDiametros(0)

        self.material = None
        self.diametro = None
        self.rellenarActivo()
        self.edit(False)

    def rellenarActivo(self):
        self.TablaDiametro.setItem(0, 0, QtWidgets.QTableWidgetItem(""))
        for i in range(1, 7):
            self.TablaDiametro.setItem(0, i, QtWidgets.QTableWidgetItem(""))
            self.TablaDiametro.item(0, i).setTextAlignment(
                QtCore.Qt.AlignmentFlag.AlignRight
                | QtCore.Qt.AlignmentFlag.AlignVCenter)
        for i in range(3):
            self.TablaMaterial.setItem(0, i, QtWidgets.QTableWidgetItem(""))
        self.TablaMaterial.item(0, 2).setTextAlignment(
            QtCore.Qt.AlignmentFlag.AlignRight
            | QtCore.Qt.AlignmentFlag.AlignVCenter)

    def selectedChanged(self):
        self.botonSelect.setEnabled(True)

    def clear(self):
        self.material = None
        self.diametro = None
        for i in range(7):
            self.TablaDiametro.item(0, i).setText("")
        for i in range(3):
            self.TablaMaterial.item(0, i).setText("")

    def add(self, i=None, j=None, load=True):
        if load:
            i = self.ComboMaterial.currentIndex()
            j = self.TablaDiametro2.currentRow()
        else:
            self.ComboMaterial.setCurrentIndex(i)
            self.TablaDiametro2.setCurrentRow(j)
        self.blockSignals(True)
        self.TablaMaterial.item(0, 0).setText(self.materiales[i][0])
        self.TablaMaterial.item(0, 1).setText(self.materiales[i][1])
        self.TablaMaterial.item(0, 2).setText(str(self.diametros[i][j][0]))
        for columna in range(7):
            self.TablaDiametro.item(0, columna).setText(
                self.TablaDiametro2.item(j, columna).text())
        self.material = self.materiales[i]
        self.diametro = self.diametros[j]
        self.blockSignals(False)
        self.emitirSignal()

    def addCustom(self, lista):
        self.blockSignals(True)
        self.TablaMaterial.item(0, 0).setText(lista[0])
        self.TablaMaterial.item(0, 1).setText(lista[1])
        self.TablaMaterial.item(0, 2).setText(str(lista[2]))
        for columna in range(7):
            self.TablaDiametro.item(0, columna).setText(str(lista[columna+3]))
        self.blockSignals(False)

    def edit(self, bool):
        self.blockSignals(True)
        inactivo = QtGui.QColor(Preferences.get("General", 'Color_ReadOnly'))
        activo = QtGui.QColor(Preferences.get("General", 'Color_Resaltado'))
        if bool:
            flags = QtCore.Qt.ItemFlag.ItemIsEditable \
               | QtCore.Qt.ItemFlag.ItemIsEnabled \
               | QtCore.Qt.ItemFlag.ItemIsSelectable
            color = activo
        else:
            flags = QtCore.Qt.ItemFlag.ItemIsEnabled \
               | QtCore.Qt.ItemFlag.ItemIsSelectable
            color = inactivo

        for i in [0, 1, 2]:
            self.TablaMaterial.item(0, i).setFlags(flags)
            self.TablaMaterial.item(0, i).setBackground(color)
        for i in range(0, 5):
            self.TablaDiametro.item(0, i).setFlags(flags)
            self.TablaDiametro.item(0, i).setBackground(color)

        self.TablaDiametro.item(0, 5).setBackground(inactivo)
        self.TablaDiametro.item(0, 6).setBackground(inactivo)
        self.TablaDiametro.item(0, 5).setFlags(
            QtCore.Qt.ItemFlag.ItemIsEnabled
            | QtCore.Qt.ItemFlag.ItemIsSelectable)
        self.TablaDiametro.item(0, 6).setFlags(
            QtCore.Qt.ItemFlag.ItemIsEnabled
            | QtCore.Qt.ItemFlag.ItemIsSelectable)
        self.blockSignals(False)

    def rellenarDiametros(self, indice):
        filas = len(self.diametros[indice])
        self.TablaDiametro2.setRowCount(filas)
        for i in range(filas):
            di = self.diametros[indice][i][4]-self.diametros[indice][i][3]*2
            V = pi/4*di**2*100/1e6
            S = pi*self.diametros[indice][i][4]*100/1e3
            self.TablaDiametro2.setRowHeight(i, 20)
            self.TablaDiametro2.setItem(i, 0, QtWidgets.QTableWidgetItem(""))
            self.TablaDiametro2.item(i, 0).setText(self.diametros[indice][i][2])
            self.TablaDiametro2.setItem(i, 1, QtWidgets.QTableWidgetItem(""))
            self.TablaDiametro2.item(i, 1).setTextAlignment(
                QtCore.Qt.AlignmentFlag.AlignRight
                | QtCore.Qt.AlignmentFlag.AlignVCenter)
            self.TablaDiametro2.item(i, 1).setText(str(round(di, 2)))
            self.TablaDiametro2.setItem(i, 2, QtWidgets.QTableWidgetItem(""))
            self.TablaDiametro2.item(i, 2).setTextAlignment(
                QtCore.Qt.AlignmentFlag.AlignRight
                | QtCore.Qt.AlignmentFlag.AlignVCenter)
            self.TablaDiametro2.item(i, 2).setText(
                str(self.diametros[indice][i][3]))
            self.TablaDiametro2.setItem(i, 3, QtWidgets.QTableWidgetItem(""))
            self.TablaDiametro2.item(i, 3).setTextAlignment(
                QtCore.Qt.AlignmentFlag.AlignRight
                | QtCore.Qt.AlignmentFlag.AlignVCenter)
            self.TablaDiametro2.item(i, 3).setText(
                str(self.diametros[indice][i][4]))
            self.TablaDiametro2.setItem(i, 4, QtWidgets.QTableWidgetItem(""))
            self.TablaDiametro2.item(i, 4).setTextAlignment(
                QtCore.Qt.AlignmentFlag.AlignRight
                | QtCore.Qt.AlignmentFlag.AlignVCenter)
            self.TablaDiametro2.item(i, 4).setText(
                str(round(self.diametros[indice][i][5], 3)))
            self.TablaDiametro2.setItem(i, 5, QtWidgets.QTableWidgetItem(""))
            self.TablaDiametro2.item(i, 5).setTextAlignment(
                QtCore.Qt.AlignmentFlag.AlignRight
                | QtCore.Qt.AlignmentFlag.AlignVCenter)
            self.TablaDiametro2.item(i, 5).setText(str(round(V, 3)))
            self.TablaDiametro2.setItem(i, 6, QtWidgets.QTableWidgetItem(""))
            self.TablaDiametro2.item(i, 6).setTextAlignment(
                QtCore.Qt.AlignmentFlag.AlignRight
                | QtCore.Qt.AlignmentFlag.AlignVCenter)
            self.TablaDiametro2.item(i, 6).setText(str(round(S, 3)))

    def getMaterial(self):
        try:
            lista = []
            lista.append(self.TablaMaterial.item(0, 0).text())
            lista.append(self.TablaMaterial.item(0, 1).text())
            lista.append(float(self.TablaMaterial.item(0, 2).text()))
            lista.append(self.TablaDiametro.item(0, 0).text())
            for i in range(1, 4):
                lista.append(float(self.TablaDiametro.item(0, i).text()))
            for i in range(4, 7):
                lista.append(float(self.TablaDiametro.item(0, i).text()))
            if self.botonEdit.isChecked():
                lista.append(self.ComboMaterial.currentIndex())
                lista.append(self.TablaDiametro2.currentRow())
            else:
                lista.append(-1)
                lista.append(-1)
        except:
            return False
        else:
            return lista

    def emitirSignal(self):
        material = self.getMaterial()
        if material:
            self.valueChanged.emit(material)


class Catalogo_Materiales_Dialog(QtWidgets.QDialog, Catalogo_Materiales):
    valueChanged = QtCore.pyqtSignal(list)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(self.tr("Select Pipe from Database"))
        buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Cancel
            | QtWidgets.QDialogButtonBox.StandardButton.Ok)
        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)
        self.layout().addWidget(buttonBox, 10, 0, 1, 3)


class Catalogo_Accesorios(QtWidgets.QWidget):
    """Widget to select/edit/add fitting to pipe"""
    valueChanged = QtCore.pyqtSignal(list)
    titles = [
        QtWidgets.QApplication.translate("equipment", "Type"),
        "D,mm", "D,in", "K", "Nº",
        QtWidgets.QApplication.translate("equipment", "Description")]

    def __init__(self, parent=None):
        super().__init__(parent)
        self.semaforo = QtCore.QSemaphore(1)
        self.accesorios = []
        self.K = 0
        gridLayout = QtWidgets.QGridLayout(self)
        self.Accesorios = QtWidgets.QTableWidget()
        self.Accesorios.verticalHeader().hide()
        self.Accesorios.setRowCount(0)
        self.Accesorios.setColumnCount(6)
        self.Accesorios.setHorizontalHeaderLabels(self.titles)
        self.Accesorios.setItemDelegateForColumn(4, SpinEditor(self))
        self.Accesorios.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.Accesorios.setSortingEnabled(True)
        self.Accesorios.horizontalHeader().setStretchLastSection(True)
        self.Accesorios.setColumnWidth(0, 80)
        self.Accesorios.setColumnWidth(1, 60)
        self.Accesorios.setColumnWidth(2, 60)
        self.Accesorios.setColumnWidth(3, 50)
        self.Accesorios.setColumnWidth(4, 50)
        self.Accesorios.cellChanged.connect(self.CalcularK)
        self.Accesorios.currentCellChanged.connect(
            self.selectedChangedAccesorios)
        gridLayout.addWidget(self.Accesorios, 0, 0, 3, 1)

        self.k_total = QtWidgets.QLabel()
        self.k_total.setText("K = 0")
        self.k_total.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter
                                  | QtCore.Qt.AlignmentFlag.AlignTop)
        self.k_total.setFixedSize(100, 25)
        gridLayout.addWidget(self.k_total, 0, 1)
        path = os.path.join(IMAGE_PATH, "button", "editDelete.png")
        self.botonBorrar = QtWidgets.QPushButton(
            QtGui.QIcon(QtGui.QPixmap(path)), self.tr("Delete"))
        path = os.path.join(IMAGE_PATH, "button", "clear.png")
        self.botonBorrar.clicked.connect(self.borrar)
        gridLayout.addWidget(self.botonBorrar, 1, 1)
        self.botonClear = QtWidgets.QPushButton(
            QtGui.QIcon(QtGui.QPixmap(path)), self.tr("Clear"))
        self.botonClear.clicked.connect(self.clear)
        gridLayout.addWidget(self.botonClear, 2, 1)

        self.line = QtWidgets.QFrame()
        self.line.setFrameShape(QtWidgets.QFrame.Shape.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Shadow.Sunken)
        gridLayout.addWidget(self.line, 3, 0, 1, 2)

        self.TablaAccesorios = QtWidgets.QTableWidget()
        self.TablaAccesorios.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.TablaAccesorios.horizontalHeader().setStretchLastSection(True)
        self.TablaAccesorios.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.TablaAccesorios.currentCellChanged.connect(self.selectedChanged)
        self.TablaAccesorios.cellDoubleClicked.connect(self.add)
        gridLayout.addWidget(self.TablaAccesorios, 4, 0, 7, 1)

        path = os.path.join(IMAGE_PATH, "button", "add.png")
        self.botonAdd = QtWidgets.QPushButton(
            QtGui.QIcon(QtGui.QPixmap(path)), self.tr("Add"))
        self.botonAdd.setDisabled(True)
        self.botonAdd.clicked.connect(self.add)
        gridLayout.addWidget(self.botonAdd, 4, 1)

        path = os.path.join(IMAGE_PATH, "equipment", "pipe", "ER.png")
        self.botonEntrada = QtWidgets.QPushButton(
            QtGui.QIcon(QtGui.QPixmap(path)), self.tr("Entrance"))
        self.botonEntrada.clicked.connect(self.EntradaRedondeada)
        gridLayout.addWidget(self.botonEntrada, 6, 1)

        path = os.path.join(IMAGE_PATH, "equipment", "pipe", "GE.png")
        self.botonExpansion = QtWidgets.QPushButton(
            QtGui.QIcon(QtGui.QPixmap(path)), self.tr("Enlargement"))
        self.botonExpansion.clicked.connect(self.expansion)
        gridLayout.addWidget(self.botonExpansion, 7, 1)

        path = os.path.join(IMAGE_PATH, "equipment", "pipe", "GC.png")
        self.botonContraccion = QtWidgets.QPushButton(
            QtGui.QIcon(QtGui.QPixmap(path)), self.tr("Contraction"))
        self.botonContraccion.clicked.connect(self.contraccion)
        gridLayout.addWidget(self.botonContraccion, 8, 1)

        path = os.path.join(IMAGE_PATH, "equipment", "pipe", "LB.png")
        self.botonCodoLargo = QtWidgets.QPushButton(
            QtGui.QIcon(QtGui.QPixmap(path)), self.tr("Long bend"))
        self.botonCodoLargo.clicked.connect(self.CodoLargo)
        gridLayout.addWidget(self.botonCodoLargo, 9, 1)

        path = os.path.join(IMAGE_PATH, "equipment", "pipe", "MB45.png")
        self.botonCodoSegmentado = QtWidgets.QPushButton(
            QtGui.QIcon(QtGui.QPixmap(path)), self.tr("Mitre bend"))
        self.botonCodoSegmentado.clicked.connect(self.CodoSegmentado)
        gridLayout.addWidget(self.botonCodoSegmentado, 10, 1)

        self.rellenarTabla()

    def rellenarTabla(self):
        self.TablaAccesorios.verticalHeader().hide()
        self.TablaAccesorios.setRowCount(0)
        self.TablaAccesorios.setColumnCount(5)
        titles = self.titles[:]
        del titles[4]
        self.TablaAccesorios.setHorizontalHeaderLabels(titles)
        self.ListaAccesorios = []
        self.diametros = []
        self.pulgadas = []
        self.k = []
        self.comboxPulgadas = []
        self.comboxmm = []

        for i, Di, key, Di_in, K in FITTING:
            if key in self.ListaAccesorios:
                indice = self.ListaAccesorios.index(key)
                self.diametros[indice].append(Di)
                self.pulgadas[indice].append(Di_in)
                self.k[indice].append(K)
                self.comboxPulgadas[indice].addItem(Di_in)
                self.comboxmm[indice].addItem(str(int(Di)))
            else:
                indice = self.TablaAccesorios.rowCount()
                self.ListaAccesorios.append(key)
                self.diametros.append([Di])
                self.pulgadas.append([Di_in])
                self.k.append([K])
                self.TablaAccesorios.setRowCount(indice+1)
                icon = os.path.join(IMAGE_PATH, "equipment", "pipe", key) \
                    + ".png"
                self.TablaAccesorios.setItem(
                    indice, 0, QtWidgets.QTableWidgetItem(QtGui.QIcon(
                        QtGui.QPixmap(icon)), key))
                self.comboxmm.append(QtWidgets.QComboBox())
                self.TablaAccesorios.setCellWidget(indice, 1, self.comboxmm[-1])
                self.comboxmm[-1].addItem(str(int(Di)))
                self.comboxmm[-1].currentIndexChanged.connect(
                    self.combommChanged)
                self.comboxPulgadas.append(QtWidgets.QComboBox())
                self.TablaAccesorios.setCellWidget(
                    indice, 2, self.comboxPulgadas[-1])
                self.comboxPulgadas[-1].addItem(Di_in)
                self.comboxPulgadas[-1].currentIndexChanged.connect(
                    self.comboinchChanged)
                self.TablaAccesorios.setItem(
                    indice, 3, QtWidgets.QTableWidgetItem("%0.3f" % K))
                self.TablaAccesorios.item(indice, 3).setTextAlignment(
                    QtCore.Qt.AlignmentFlag.AlignRight
                    | QtCore.Qt.AlignmentFlag.AlignVCenter)
                self.TablaAccesorios.setItem(
                    indice, 4, QtWidgets.QTableWidgetItem(FITTING_DESC[key]))
                self.TablaAccesorios.setRowHeight(
                    self.TablaAccesorios.rowCount()-1, 20)

        self.TablaAccesorios.resizeColumnToContents(0)
        self.TablaAccesorios.setColumnWidth(1, 60)
        self.TablaAccesorios.setColumnWidth(2, 60)
        self.TablaAccesorios.setColumnWidth(3, 50)

    def expansion(self):
        title = self.tr("Expansion")
        icon = os.path.join(IMAGE_PATH, "equipment", "pipe", "GE.png")
        parameter = [
            self.tr("Input diameter"),
            self.tr("Output diameter")]
        dialog = Dialog(0, title, icon, parameter)
        if dialog.exec():
            indice = self.Accesorios.rowCount()
            self.Accesorios.setRowCount(indice+1)
            if dialog.checkBox.isChecked():
                type = "GE"
                txt = self.tr("Gradual enlargement")
                dia = "%0.3f->%0.3f θ=%iº" % (
                    dialog.D1.value, dialog.D2.value, dialog.angulo.value())
            else:
                icon = QtGui.QIcon(QtGui.QPixmap(
                    os.path.join(IMAGE_PATH, "equipment", "pipe", "SE.png")))
                type = "SE"
                txt = self.tr("Sudden enlargement")
                dia = "%0.3f->%0.3f" % (dialog.D1.value, dialog.D2.value)
            self.Accesorios.setItem(
                indice, 0, QtWidgets.QTableWidgetItem(icon, type))
            self.Accesorios.setSpan(indice, 1, 1, 2)
            self.Accesorios.setItem(indice, 1, QtWidgets.QTableWidgetItem(dia))
            self.Accesorios.setItem(
                indice, 3, QtWidgets.QTableWidgetItem(representacion(dialog.K, 3)))
            self.Accesorios.item(indice, 3).setTextAlignment(
                QtCore.Qt.AlignmentFlag.AlignRight
                | QtCore.Qt.AlignmentFlag.AlignVCenter)
            self.Accesorios.setItem(indice, 4, QtWidgets.QTableWidgetItem(str(1)))
            self.Accesorios.item(indice, 4).setTextAlignment(
                QtCore.Qt.AlignmentFlag.AlignRight
                | QtCore.Qt.AlignmentFlag.AlignVCenter)
            self.Accesorios.setItem(indice, 5, QtWidgets.QTableWidgetItem(txt))
            self.Accesorios.setRowHeight(indice, 20)
            self.CalcularK()

    def contraccion(self):
        title = self.tr("Contraction")
        icon = os.path.join(IMAGE_PATH, "equipment", "pipe", "GC.png")
        parameter = [
            self.tr("Input diameter"),
            self.tr("Output diameter")]
        dialog = Dialog(0, title, icon, parameter)
        if dialog.exec():
            indice = self.Accesorios.rowCount()
            self.Accesorios.setRowCount(indice+1)
            if dialog.checkBox.isChecked():
                type = "GC"
                txt = self.tr("Gradual contraction")
                dia = "%0.3f->%0.3f θ=%iº" % (
                    dialog.D1.value, dialog.D2.value, dialog.angulo.value())
                icon = QtGui.QIcon(QtGui.QPixmap(
                    os.path.join(IMAGE_PATH, "equipment", "pipe", "GC.png")))
            else:
                icon = QtGui.QIcon(QtGui.QPixmap(
                    os.path.join(IMAGE_PATH, "equipment", "pipe", "SC.png")))
                type = "SC"
                txt = self.tr("Sudden contraction")
                dia = "%0.3f->%0.3f" % (dialog.D1.value, dialog.D2.value)
            self.Accesorios.setItem(indice, 0,
                                    QtWidgets.QTableWidgetItem(icon, type))
            self.Accesorios.setSpan(indice, 1, 1, 2)
            self.Accesorios.setItem(indice, 1, QtWidgets.QTableWidgetItem(dia))
            self.Accesorios.setItem(
                indice, 3, QtWidgets.QTableWidgetItem(representacion(dialog.K, 3)))
            self.Accesorios.item(indice, 3).setTextAlignment(
                QtCore.Qt.AlignmentFlag.AlignRight
                | QtCore.Qt.AlignmentFlag.AlignVCenter)
            self.Accesorios.setItem(indice, 4, QtWidgets.QTableWidgetItem(str(1)))
            self.Accesorios.item(indice, 4).setTextAlignment(
                QtCore.Qt.AlignmentFlag.AlignRight
                | QtCore.Qt.AlignmentFlag.AlignVCenter)
            self.Accesorios.setItem(indice, 5, QtWidgets.QTableWidgetItem(txt))
            self.Accesorios.setRowHeight(indice, 20)
            self.CalcularK()

    def EntradaRedondeada(self):
        title = self.tr("Rounded entrance")
        icon = os.path.join(IMAGE_PATH, "equipment", "pipe", "ER.png")
        parameter = [
            self.tr("Exit radio"),
            self.tr("Pipe Diameter")]
        dialog = Dialog(3, title, icon, parameter)
        if dialog.exec():
            indice = self.Accesorios.rowCount()
            self.Accesorios.setRowCount(indice+1)
            self.Accesorios.setItem(indice, 0, QtWidgets.QTableWidgetItem(
                QtGui.QIcon(QtGui.QPixmap(icon)), "ER"))
            self.Accesorios.setSpan(indice, 1, 1, 2)
            self.Accesorios.setItem(indice, 1, QtWidgets.QTableWidgetItem(
                f"r={dialog.D1.value.mm:0.3f}, D={dialog.D2.value.mm:0.3f}"))
            self.Accesorios.setItem(indice, 3, QtWidgets.QTableWidgetItem(
                representacion(dialog.K, 3)))
            self.Accesorios.item(indice, 3).setTextAlignment(
                QtCore.Qt.AlignmentFlag.AlignRight
                | QtCore.Qt.AlignmentFlag.AlignVCenter)
            self.Accesorios.setItem(indice, 4, QtWidgets.QTableWidgetItem(str(1)))
            self.Accesorios.item(indice, 4).setTextAlignment(
                QtCore.Qt.AlignmentFlag.AlignRight
                | QtCore.Qt.AlignmentFlag.AlignVCenter)
            self.Accesorios.setItem(indice, 5, QtWidgets.QTableWidgetItem(
                self.tr("Pipe entrance rounded intake")))
            self.Accesorios.setRowHeight(indice, 20)
            self.CalcularK()

    def CodoLargo(self):
        title = self.tr("Long Pipe Bend")
        icon = os.path.join(IMAGE_PATH, "equipment", "pipe", "LB.png")
        parameter = [
            self.tr("Bend radio"),
            self.tr("Pipe diameter")]
        dialog = Dialog(2, title, icon, parameter)
        if dialog.exec():
            indice = self.Accesorios.rowCount()
            self.Accesorios.setRowCount(indice+1)
            self.Accesorios.setItem(indice, 0, QtWidgets.QTableWidgetItem(
                QtGui.QIcon(QtGui.QPixmap(icon)), "LB"))
            self.Accesorios.setSpan(indice, 1, 1, 2)
            self.Accesorios.setItem(indice, 1, QtWidgets.QTableWidgetItem(
                "r=%0.3f, D=%0.3f" % (dialog.D1.value.mm, dialog.D2.value.mm)))
            self.Accesorios.setItem(indice, 3, QtWidgets.QTableWidgetItem(
                representacion(dialog.K, 3)))
            self.Accesorios.item(indice, 3).setTextAlignment(
                QtCore.Qt.AlignmentFlag.AlignRight
                | QtCore.Qt.AlignmentFlag.AlignVCenter)
            self.Accesorios.setItem(indice, 4, QtWidgets.QTableWidgetItem(str(1)))
            self.Accesorios.item(indice, 4).setTextAlignment(
                QtCore.Qt.AlignmentFlag.AlignRight
                | QtCore.Qt.AlignmentFlag.AlignVCenter)
            self.Accesorios.setItem(indice, 5, QtWidgets.QTableWidgetItem(title))
            self.Accesorios.setRowHeight(indice, 20)
            self.CalcularK()

    def CodoSegmentado(self):
        title = self.tr("Mitre bend with custom angle")
        icon = os.path.join(IMAGE_PATH, "equipment", "pipe", "MB45.png")
        parameter = [self.tr("Pipe diameter"), ""]
        dialog = Dialog(1, title, icon, parameter)
        if dialog.exec():
            indice = self.Accesorios.rowCount()
            self.Accesorios.setRowCount(indice+1)
            self.Accesorios.setItem(indice, 0, QtWidgets.QTableWidgetItem(
                QtGui.QIcon(QtGui.QPixmap(icon)), "MBx"))
            self.Accesorios.setItem(indice, 1, QtWidgets.QTableWidgetItem(
                "%0.3f" % dialog.D1.value.mm))
            self.Accesorios.setItem(indice, 2, QtWidgets.QTableWidgetItem(
                "θ=%iº" % dialog.angulo.value()))
            self.Accesorios.setItem(indice, 3, QtWidgets.QTableWidgetItem(
                representacion(dialog.K, 3)))
            self.Accesorios.item(indice, 3).setTextAlignment(
                QtCore.Qt.AlignmentFlag.AlignRight
                | QtCore.Qt.AlignmentFlag.AlignVCenter)
            self.Accesorios.setItem(indice, 4, QtWidgets.QTableWidgetItem(str(1)))
            self.Accesorios.item(indice, 4).setTextAlignment(
                QtCore.Qt.AlignmentFlag.AlignRight
                | QtCore.Qt.AlignmentFlag.AlignVCenter)
            self.Accesorios.setItem(indice, 5, QtWidgets.QTableWidgetItem(
                self.tr("Mitre Bend")))
            self.Accesorios.setRowHeight(indice, 20)
            self.CalcularK()

    def combommChanged(self, i):
        if self.semaforo.available() > 0:
            self.semaforo.acquire(1)
            indice = self.comboxmm.index(self.focusWidget())
            self.comboxPulgadas[indice].setCurrentIndex(i)
            self.TablaAccesorios.item(indice, 3).setText(
                "%0.3f" % self.k[indice][i])
            self.TablaAccesorios.selectRow(indice)
            self.semaforo.release(1)

    def comboinchChanged(self, i):
        if self.semaforo.available() > 0:
            self.semaforo.acquire(1)
            indice = self.comboxPulgadas.index(self.focusWidget())
            self.comboxmm[indice].setCurrentIndex(i)
            self.TablaAccesorios.item(indice, 3).setText(
                "%0.3f" % self.k[indice][i])
            self.TablaAccesorios.selectRow(indice)
            self.semaforo.release(1)

    def borrar(self):
        del self.accesorios[self.Accesorios.currentRow()]
        self.Accesorios.removeRow(self.Accesorios.currentRow())
        if self.accesorios:
            self.CalcularK()
        else:
            self.clear()

    def clear(self):
        self.Accesorios.clearContents()
        self.Accesorios.setRowCount(0)
        self.k_total.setText("K = 0")
        self.K = 0
        self.accesorios = []

    def CalcularK(self):
        k = 0
        for i in range(self.Accesorios.rowCount()):
            count = float(self.Accesorios.item(i, 3).text())
            ki = float(self.Accesorios.item(i, 4).text())
            k += count*ki
        self.k_total.setText("K = %0.3f" % k)
        self.K = k

    def selectedChanged(self, i):
        self.botonAdd.setEnabled(i > -1)

    def selectedChangedAccesorios(self, newi, newj, oldi, oldj):
        self.Accesorios.openPersistentEditor(self.Accesorios.item(newi, 4))
        self.Accesorios.closePersistentEditor(self.Accesorios.item(oldi, 4))

    def add(self, indicetabla=None, indicecombo=None, emit=True):
        self.Accesorios.blockSignals(True)
        if indicetabla is None:
            indicetabla = self.TablaAccesorios.currentRow()
            indicecombo = self.comboxmm[indicetabla].currentIndex()

        if [indicetabla, indicecombo] not in self.accesorios:
            indice = self.Accesorios.rowCount()
            self.accesorios.append([indicetabla, indicecombo])
            self.Accesorios.setRowCount(indice+1)
            self.Accesorios.setItem(indice, 0, QtWidgets.QTableWidgetItem(
                self.TablaAccesorios.item(indicetabla, 0)))
            self.Accesorios.item(indice, 0).setFlags(
                QtCore.Qt.ItemFlag.ItemIsEnabled | QtCore.Qt.ItemFlag.ItemIsSelectable)
            self.Accesorios.setItem(indice, 1, QtWidgets.QTableWidgetItem(
                self.comboxmm[indicetabla].currentText()))
            self.Accesorios.item(indice, 1).setTextAlignment(
                QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter)
            self.Accesorios.item(indice, 1).setFlags(
                QtCore.Qt.ItemFlag.ItemIsEnabled | QtCore.Qt.ItemFlag.ItemIsSelectable)
            self.Accesorios.setItem(indice, 2, QtWidgets.QTableWidgetItem(
                self.comboxPulgadas[indicetabla].currentText()))
            self.Accesorios.item(indice, 2).setTextAlignment(
                QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter)
            self.Accesorios.item(indice, 2).setFlags(
                QtCore.Qt.ItemFlag.ItemIsEnabled | QtCore.Qt.ItemFlag.ItemIsSelectable)
            self.Accesorios.setItem(indice, 3, QtWidgets.QTableWidgetItem(
                self.TablaAccesorios.item(indicetabla, 3)))
            self.Accesorios.item(indice, 3).setFlags(
                QtCore.Qt.ItemFlag.ItemIsEnabled | QtCore.Qt.ItemFlag.ItemIsSelectable)
            self.Accesorios.setItem(indice, 4, QtWidgets.QTableWidgetItem("1"))
            self.Accesorios.item(indice, 4).setTextAlignment(
                QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter)
            self.Accesorios.setItem(indice, 5, QtWidgets.QTableWidgetItem(
                self.TablaAccesorios.item(indicetabla, 4)))
            self.Accesorios.item(indice, 5).setFlags(
                QtCore.Qt.ItemFlag.ItemIsEnabled | QtCore.Qt.ItemFlag.ItemIsSelectable)
            self.Accesorios.setRowHeight(indice, 20)
        else:
            indice = self.accesorios.index([indicetabla, indicecombo])
            self.Accesorios.item(indice, 4).setText(
                str(int(self.Accesorios.item(indice, 4).text())+1))
        self.CalcularK()
        self.Accesorios.blockSignals(False)

        if emit:
            array = []
            for i, [indicetabla, indicecombo] in enumerate(self.accesorios):
                lista = [indicetabla, indicecombo]
                lista.append(float(self.Accesorios.item(i, 3).text()))
                lista.append(int(self.Accesorios.item(i, 4).text()))
                lista.append(str(self.Accesorios.item(i, 0).text()))
                lista.append(str(self.Accesorios.item(i, 1).text()))
                lista.append(str(self.Accesorios.item(i, 2).text()))
                lista.append(str(self.Accesorios.item(i, 5).text()))
                array.append(lista)
            self.valueChanged.emit(array)


class UI_equipment(UI_equip):
    """Pipe equipment edition dialog"""
    Equipment = Pipe()

    def __init__(self, equipment=None, parent=None):
        """
        equipment: Initial equipment instance to model
        """
        super().__init__(Pipe, entrada=False, salida=False, parent=parent)

        # Catalog tab
        self.tabCatalogo = Catalogo_Materiales()
        self.tabCatalogo.valueChanged.connect(
            partial(self.changeParams, "material"))
        self.tabWidget.insertTab(1, self.tabCatalogo, self.tr("Catalog"))

        # Fitting tab
        self.tabAccesorios = Catalogo_Accesorios()
        self.tabAccesorios.valueChanged.connect(
            partial(self.changeParams, "accesorios"))
        self.tabWidget.insertTab(2, self.tabAccesorios, self.tr("Fittings"))

        # Calculate tab
        lyt = QtWidgets.QGridLayout(self.tabCalculo)
#        lyt.addWidget(QtGui.QLabel(self.tr("Modo:")),0,0)
#        self.Modo=QtGui.QComboBox()
#        lyt.addWidget(self.Modo,0,1,1,1)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Method")), 1, 0)
        self.metodo = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_METODO:
            self.metodo.addItem(txt)
        self.metodo.currentIndexChanged.connect(self.metodoCambiado)
        lyt.addWidget(self.metodo, 1, 1, 1, 5)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 1, 0, 1, 6)

        lyt.addWidget(QtWidgets.QLabel(self.tr("Length")), 2, 0)
        self.l = Entrada_con_unidades(Length, resaltado=True)
        self.l.valueChanged.connect(partial(self.changeParams, "l"))
        lyt.addWidget(self.l, 2, 1)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Elevation")), 3, 0)
        self.h = Entrada_con_unidades(Length, min=float("-inf"))
        self.h.valueChanged.connect(partial(self.changeParams, "h"))
        lyt.addWidget(self.h, 3, 1)
        self.labelC = QtWidgets.QLabel(self.tr("C Factor"))
        lyt.addWidget(self.labelC, 2, 3)
        self.C = Entrada_con_unidades(float)
        self.C.valueChanged.connect(partial(self.changeParams, "C"))
        lyt.addWidget(self.C, 2, 4)
        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 4, 0, 1, 6)

        lyt.addWidget(QtWidgets.QLabel(self.tr("Thermal Mode")), 11, 0)
        self.thermal = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_THERMAL:
            self.thermal.addItem(txt)
        self.thermal.currentIndexChanged.connect(self.Heat_visible)
        lyt.addWidget(self.thermal, 11, 1, 1, 5)

        self.groupBox_Heat = QtWidgets.QGroupBox(
            self.tr("Heat transfer to surroundings"))
        lyt.addWidget(self.groupBox_Heat, 12, 0, 3, 2)
        layout = QtWidgets.QGridLayout(self.groupBox_Heat)
        layout.addWidget(QtWidgets.QLabel(
            self.tr("External temperature")), 1, 1)
        self.T_ext = Entrada_con_unidades(Temperature)
        self.T_ext.valueChanged.connect(partial(self.changeParams, "T_ext"))
        layout.addWidget(self.T_ext, 1, 2)
        layout.addWidget(QtWidgets.QLabel(
            self.tr("Heat transfer coefficient")), 2, 1)
        self.U = Entrada_con_unidades(HeatTransfCoef)
        self.U.valueChanged.connect(partial(self.changeParams, "U"))
        layout.addWidget(self.U, 2, 2)

        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 12, 0, 1, 6)
        self.labelQ = QtWidgets.QLabel(self.tr("Heat Flux"))
        lyt.addWidget(self.labelQ, 13, 0)
        self.Calor = Entrada_con_unidades(Power)
        self.Calor.valueChanged.connect(partial(self.changeParams, "Q"))
        lyt.addWidget(self.Calor, 13, 1)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 14, 0, 1, 6)

        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 15, 0, 1, 6)

        group = QtWidgets.QGroupBox(self.tr("Results"))
        lyt.addWidget(group, 16, 0, 1, 6)
        layout = QtWidgets.QGridLayout(group)
        layout.addWidget(QtWidgets.QLabel(self.tr("ΔP total")), 0, 0)
        self.DeltaP = Entrada_con_unidades(
            Pressure, retornar=False, readOnly=True)
        layout.addWidget(self.DeltaP, 0, 1)
        layout.addWidget(QtWidgets.QLabel(self.tr("ΔP friction")), 1, 0)
        self.DeltaP_f = Entrada_con_unidades(
            Pressure, retornar=False, readOnly=True)
        layout.addWidget(self.DeltaP_f, 1, 1)
        layout.addWidget(QtWidgets.QLabel(self.tr("ΔP fittings")), 2, 0)
        self.DeltaP_ac = Entrada_con_unidades(
            Pressure, retornar=False, readOnly=True)
        layout.addWidget(self.DeltaP_ac, 2, 1)
        layout.addWidget(QtWidgets.QLabel(self.tr("ΔP elevation")), 3, 0)
        self.DeltaP_h = Entrada_con_unidades(
            Pressure, retornar=False, readOnly=True)
        layout.addWidget(self.DeltaP_h, 3, 1)
        layout.addWidget(QtWidgets.QLabel(self.tr("ΔP acceleration")),
                         4, 0)
        self.DeltaP_v = Entrada_con_unidades(
            Pressure, retornar=False, readOnly=True)
        layout.addWidget(self.DeltaP_v, 4, 1)

        layout.addWidget(QtWidgets.QLabel("ΔP/100ft"), 0, 4)
        self.DeltaP_100ft = Entrada_con_unidades(Pressure, retornar=False,
                                                 readOnly=True)
        layout.addWidget(self.DeltaP_100ft, 0, 5)
        layout.addWidget(QtWidgets.QLabel(self.tr("Speed")), 1, 4)
        self.V = Entrada_con_unidades(Speed, retornar=False, readOnly=True)
        layout.addWidget(self.V, 1, 5)
        layout.addWidget(QtWidgets.QLabel(self.tr("Friction factor")),
                         2, 4)
        self.f = Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.f, 2, 5)
        layout.addWidget(QtWidgets.QLabel(self.tr("Reynolds")), 3, 4)
        self.Re = Entrada_con_unidades(float, tolerancia=5, decimales=1,
                                       readOnly=True)
        layout.addWidget(self.Re, 3, 5)
        layout.addWidget(QtWidgets.QLabel(self.tr("Output Temperature")), 4, 4)
        self.Tout = Entrada_con_unidades(
            Temperature, decimales=2, retornar=False, readOnly=True)
        layout.addWidget(self.Tout, 4, 5)

        layout.addWidget(QtWidgets.QLabel(self.tr("Flow region")), 0, 7)
        self.regimen = Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.regimen, 0, 8)
        layout.addWidget(QtWidgets.QLabel(self.tr("Reynolds liquid")), 1, 7)
        self.Re_liq = Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.Re_liq, 1, 8)
        layout.addWidget(QtWidgets.QLabel(self.tr("Reynolds gas")), 2, 7)
        self.Re_gas = Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.Re_gas, 2, 8)
        layout.addWidget(QtWidgets.QLabel(self.tr("Friction liquid")), 3, 7)
        self.f_liq = Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.f_liq, 3, 8)
        layout.addWidget(QtWidgets.QLabel(self.tr("Friction gas")), 4, 7)
        self.f_gas = Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.f_gas, 4, 8)

        layout.addItem(QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Minimum), 0, 3)
        layout.addItem(QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Minimum), 0, 6)

        # Cost tab
        lyt = QtWidgets.QGridLayout(self.tabCostos)
        self.labelAvisoCostos = QtWidgets.QLabel(self.tr(
    "Costs only available for steel pipes"))
        self.labelAvisoCostos.setVisible(False)
        self.labelAvisoCostos.setStyleSheet("color: #FF0000;")
        lyt.addWidget(self.labelAvisoCostos, 0, 0, 1, 2)

        self.Costos = CostData(self.Equipment)
        self.Costos.valueChanged.connect(self.changeParamsCoste)
        lyt.addWidget(self.Costos, 1, 0, 1, 2)

        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 2, 0, 1, 2)
        group = QtWidgets.QGroupBox(self.tr("Stimated Cost"))
        lyt.addWidget(group, 3, 0, 1, 2)
        layout = QtWidgets.QGridLayout(group)
        layout.addWidget(QtWidgets.QLabel(self.tr("Purchase Costs")), 0, 1)
        self.C_adq = Entrada_con_unidades(
            Currency, retornar=False, readOnly=True)
        layout.addWidget(self.C_adq, 0, 2)
        layout.addWidget(QtWidgets.QLabel(self.tr("Installed Costs")), 1, 1)
        self.C_inst = Entrada_con_unidades(
            Currency, retornar=False, readOnly=True)
        layout.addWidget(self.C_inst, 1, 2)
        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 4, 0, 1, 2)

        self.Heat_visible(0)
        self.metodoCambiado(0)
        if equipment:
            self.setEquipment(equipment)

    def Heat_visible(self, indice):
        if indice == 1:
            self.groupBox_Heat.setVisible(False)
            self.labelQ.setVisible(True)
            self.Calor.setVisible(True)
        elif indice == 2:
            self.groupBox_Heat.setVisible(True)
            self.labelQ.setVisible(False)
            self.Calor.setVisible(False)
        else:
            self.groupBox_Heat.setVisible(False)
            self.labelQ.setVisible(False)
            self.Calor.setVisible(False)
        self.changeParams("thermal", indice)

    def metodoCambiado(self, indice):
        self.C.setEnabled(indice == 1)
        self.labelC.setEnabled(indice == 1)
        self.changeParams("metodo", indice)

    def rellenar(self):
        UI_equip.rellenar(self)
        if self.Equipment.status == 1:
            if self.Equipment.statusCoste:
                self.labelAvisoCostos.setVisible(False)
            else:
                self.labelAvisoCostos.setVisible(True)

    def rellenarInput(self):
        self.blockSignals(True)
        UI_equip.rellenarInput(self)
        if self.Equipment.kwargs["material"]:
            if self.Equipment.kwargs["material"][10] >= 0:
                self.tabCatalogo.add(
                    self.Equipment.kwargs["material"][10],
                    self.Equipment.kwargs["material"][11], False)
            else:
                self.tabCatalogo.addCustom(self.Equipment.kwargs["material"])

        if self.Equipment.kwargs["accesorios"]:
            self.tabAccesorios.clear()
            for accesorio in self.Equipment.kwargs["accesorios"]:
                for numero in range(accesorio[3]):
                    self.tabAccesorios.add(accesorio[0], accesorio[1], False)
        self.blockSignals(False)


if __name__ == "__main__":
    import sys
    from lib.corriente import Corriente
    app = QtWidgets.QApplication(sys.argv)
    # agua = Corriente(T=300, P=101325, caudalMasico=1, fraccionMolar=[1.])
    # tuberia = Pipe(entrada=agua, metodo=0, l=5, material=["Steel (ANSI)", "Sch. 40", 0.12192, '6"', 152.4, 11.43, 175.26, 42.07, 1.824, 55.06, -1, 2])
    dialogo = UI_equipment()
    dialogo.show()
    sys.exit(app.exec())
