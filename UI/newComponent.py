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
# Module to implement new component by group contribution methods
#   -newComponent: Main dialog class with common functionality
#   -UI_Contribution: Definition for group contribution
#   -View_Contribution: Dialog to show the properties of a group contribution
###############################################################################


import os
from functools import partial

from PyQt5 import QtCore, QtGui, QtWidgets

from lib.compuestos import (Joback, Constantinou_Gani, Wilson_Jasperson,
                            Marrero_Pardillo, Elliott, Ambrose)
from lib import sql
from lib.unidades import Temperature
from UI.delegate import SpinEditor
from UI.viewComponents import View_Component
from UI.widgets import Entrada_con_unidades, Status


class newComponent(QtWidgets.QDialog):
    """Main dialog class with common functionality"""

    def loadUI(self):
        """Define common widget for chid class"""
        layoutBottom = QtWidgets.QHBoxLayout()
        self.status = Status()
        layoutBottom.addWidget(self.status)
        self.buttonShowDetails = QtWidgets.QPushButton(
            QtWidgets.QApplication.translate("pychemqt", "Show Details"))
        self.buttonShowDetails.clicked.connect(self.showDetails)
        self.buttonShowDetails.setEnabled(False)
        layoutBottom.addWidget(self.buttonShowDetails)
        self.btonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Cancel |
            QtWidgets.QDialogButtonBox.Save)
        self.btonBox.button(QtWidgets.QDialogButtonBox.Save).setEnabled(False)
        self.btonBox.accepted.connect(self.save)
        self.btonBox.rejected.connect(self.reject)
        layoutBottom.addWidget(self.btonBox)
        # self.layout().addLayout(layoutBottom, 30, 0, 1, 6)
        self.layout().addLayout(layoutBottom)

    def save(self):
        """Save new componente in user database"""
        elemento = self.unknown.export2Component()
        sql.inserElementsFromArray(sql.databank_Custom_name, [elemento])
        Dialog = View_Component(1001+sql.N_comp_Custom)
        Dialog.show()
        QtWidgets.QDialog.accept(self)

    def changeParams(self, parametro, valor):
        self.calculo(**{parametro: valor})

    def calculo(self, **kwargs):
        self.status.setState(4)
        self.unknown(**kwargs)
        self.status.setState(self.unknown.status, self.unknown.msg)
        self.buttonShowDetails.setEnabled(self.unknown.status)
        self.btonBox.button(QtWidgets.QDialogButtonBox.Save).setEnabled(
            self.unknown.status)

    def showDetails(self):
        """Show details of new component"""
        dialog = self.ViewDetails(self.unknown)
        dialog.exec_()


class View_Contribution(QtWidgets.QDialog):
    """Ventana que muestra las propiedades del nuevo componente definido mediante los metodos de contribucion de grupo"""
    def __init__(self, petroleo=None, parent=None):
        super(View_Contribution, self).__init__(parent)
        self.setWindowTitle(QtWidgets.QApplication.translate("pychemqt", "Group Contribution new component"))
        layout = QtWidgets.QGridLayout(self)

        self.nombre=QtWidgets.QLabel()
        layout.addWidget(self.nombre,1,1,1,5)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Molecular Weight")),2,1)
        self.M=Entrada_con_unidades(float, readOnly=True, textounidad="g/mol")
        layout.addWidget(self.M,2,2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Boiling Point")),3,1)
        self.Tb=Entrada_con_unidades(unidades.Temperature, readOnly=True)
        layout.addWidget(self.Tb,3,2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Melting Point")),4,1)
        self.Tf=Entrada_con_unidades(unidades.Temperature, readOnly=True)
        layout.addWidget(self.Tf,4,2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Tc")),5,1)
        self.Tc=Entrada_con_unidades(unidades.Temperature, readOnly=True)
        layout.addWidget(self.Tc,5,2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Pc")),6,1)
        self.Pc=Entrada_con_unidades(unidades.Pressure, readOnly=True)
        layout.addWidget(self.Pc,6,2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Vc")),7,1)
        self.Vc=Entrada_con_unidades(unidades.SpecificVolume, readOnly=True)
        layout.addWidget(self.Vc,7,2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "ΔHf", None)),8,1)
        self.Hf=Entrada_con_unidades(unidades.Enthalpy, readOnly=True)
        layout.addWidget(self.Hf,8,2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "ΔGf", None)),9,1)
        self.Gf=Entrada_con_unidades(unidades.Enthalpy, readOnly=True)
        layout.addWidget(self.Gf,9,2)

        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "ΔHm", None)),2,4)
        self.Hm=Entrada_con_unidades(unidades.Enthalpy, readOnly=True)
        layout.addWidget(self.Hm,2,5)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "ΔHv", None)),3,4)
        self.Hv=Entrada_con_unidades(unidades.Enthalpy, readOnly=True)
        layout.addWidget(self.Hv,3,5)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Cpa")),4,4)
        self.cpa=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.cpa,4,5)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Cpb")),5,4)
        self.cpb=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.cpb,5,5)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Cpc")),6,4)
        self.cpc=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.cpc,6,5)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Cpd")),7,4)
        self.cpd=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.cpd,7,5)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "μa", None)),8,4)
        self.mua=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.mua,8,5)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "μb", None)),9,4)
        self.mub=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.mub,9,5)

        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "SG 60ºF", None)),2,7)
        self.gravity=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.gravity,2,8)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "API gravity")),3,7)
        self.API=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.API,3,8)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Watson Factor")),4,7)
        self.watson=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.watson,4,8)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Acentric factor")),5,7)
        self.f_acent=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.f_acent,5,8)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Rackett")),6,7)
        self.rackett=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.rackett,6,8)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Volume Liquid Constant")),7,7)
        self.Vliq=Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.Vliq,7,8)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Solubility Parameter")),8,7)
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



class Ui_Contribution(newComponent):
    """Dialog to define hypotethical new component with several group
    contribucion methods"""
    ViewDetails = View_Contribution

    def __init__(self, metodo, parent=None):
        """Metodo: name of group contribution method:
            Joback
            Constantinou-Gani
            Wilson-Jasperson
            Marrero-Pardillo
            Elliott
            Ambrose
        """
        super(Ui_Contribution, self).__init__(parent)
        self.setWindowTitle(QtWidgets.QApplication.translate(
            "pychemqt", "Select the component group for method") +" "+ metodo)

        self.grupo = []
        self.indices = []
        self.contribucion = []
        self.metodo = metodo
        layout = QtWidgets.QGridLayout(self)
        self.Grupos = QtWidgets.QTableWidget()
        self.Grupos.verticalHeader().hide()
        self.Grupos.setRowCount(0)
        self.Grupos.setColumnCount(2)
        self.Grupos.setHorizontalHeaderItem(0, QtWidgets.QTableWidgetItem("Nk"))
        self.Grupos.setHorizontalHeaderItem(1, QtWidgets.QTableWidgetItem(
            QtWidgets.QApplication.translate("pychemqt", "Group")))
        self.Grupos.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.Grupos.setSortingEnabled(True)
        self.Grupos.horizontalHeader().setStretchLastSection(True)
        self.Grupos.setColumnWidth(0, 50)
        self.Grupos.setItemDelegateForColumn(0, SpinEditor(self))
        self.Grupos.cellChanged.connect(self.cellChanged)
        self.Grupos.setEditTriggers(QtWidgets.QAbstractItemView.AllEditTriggers)
        layout.addWidget(self.Grupos, 0, 0, 3, 3)

        self.Formula = QtWidgets.QLabel()
        font = QtGui.QFont()
        font.setPointSize(12)
        self.Formula.setFont(font)
        self.Formula.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
        self.Formula.setFixedHeight(50)
        layout.addWidget(self.Formula, 0, 3)
        self.botonBorrar = QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"]+"/images/button/editDelete.png")),
            QtWidgets.QApplication.translate("pychemqt", "Delete"))
        self.botonBorrar.clicked.connect(self.borrar)
        layout.addWidget(self.botonBorrar, 1, 3)
        self.botonClear = QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"]+"/images/button/clear.png")),
            QtWidgets.QApplication.translate("pychemqt", "Clear"))
        self.botonClear.clicked.connect(self.clear)
        layout.addWidget(self.botonClear, 2, 3)

        self.line = QtWidgets.QFrame()
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        layout.addWidget(self.line, 3, 0, 1, 4)

        self.TablaContribuciones = QtWidgets.QListWidget()
        self.TablaContribuciones.currentItemChanged.connect(self.selectedChanged)
        self.TablaContribuciones.itemDoubleClicked.connect(self.add)
        layout.addWidget(self.TablaContribuciones, 4, 0, 7, 3)
        self.botonAdd = QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"]+"/images/button/add.png")),
            QtWidgets.QApplication.translate("pychemqt", "Add"))
        self.botonAdd.setDisabled(True)
        self.botonAdd.clicked.connect(self.add)
        layout.addWidget(self.botonAdd, 4, 3)
        layout.addItem(QtWidgets.QSpacerItem(20, 20, QtWidgets.QSizePolicy.Expanding,
                                         QtWidgets.QSizePolicy.Expanding), 5, 1, 1, 1)

        # Show widget for specific method
        if metodo in ["Constantinou", "Wilson"]:
            self.Order1 = QtWidgets.QRadioButton(
                QtWidgets.QApplication.translate("pychemqt", "1st order"))
            self.Order1.setChecked(True)
            self.Order1.toggled.connect(self.Order)
            layout.addWidget(self.Order1, 6, 3)
            self.Order2 = QtWidgets.QRadioButton(
                QtWidgets.QApplication.translate("pychemqt", "2nd order"))
            layout.addWidget(self.Order2, 7, 3)

        if metodo == "Wilson":
            layout.addWidget(QtWidgets.QLabel(
                QtWidgets.QApplication.translate("pychemqt", "Rings")), 8, 3)
            self.anillos = QtWidgets.QSpinBox()
            self.anillos.valueChanged.connect(partial(self.changeParams, "ring"))
            layout.addWidget(self.anillos, 9, 3)

        if metodo == "Marrero":
            layout.addWidget(QtWidgets.QLabel(
                QtWidgets.QApplication.translate("pychemqt", "Atoms")), 8, 3)
            self.Atomos = QtWidgets.QSpinBox()
            self.Atomos.valueChanged.connect(partial(self.changeParams, "atomos"))
            layout.addWidget(self.Atomos, 9, 3)

        if metodo == "Ambrose":
            layout.addWidget(QtWidgets.QLabel(
                QtWidgets.QApplication.translate("pychemqt", "Platt number")), 8, 3)
            self.Platt = QtWidgets.QSpinBox()
            self.Platt.setToolTip(QtWidgets.QApplication.translate(
                "pychemqt", "The Platt number is the number of pairs of carbon \
                atoms which are separated \nby three carbon-carbon bonds and \
                is an indicator of the degree of branching in the molecule.\n\
                The Platt number of an n-alkane is equal to the number of \
                carbons minus three"))
            self.Platt.valueChanged.connect(partial(self.changeParams, "platt"))
            layout.addWidget(self.Platt, 9, 3)

        layout.addItem(QtWidgets.QSpacerItem(20, 20, QtWidgets.QSizePolicy.Expanding,
                                         QtWidgets.QSizePolicy.Expanding), 10, 1, 1, 1)
        layout.addItem(QtWidgets.QSpacerItem(20, 20, QtWidgets.QSizePolicy.Expanding,
                                         QtWidgets.QSizePolicy.Expanding), 11, 0, 1, 4)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Name")), 12, 0)
        self.nombre = QtWidgets.QLineEdit()
        self.nombre.textChanged.connect(partial(self.changeParams, "name"))
        layout.addWidget(self.nombre, 12, 1, 1, 3)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Molecular Weight")), 13, 0)
        self.M = Entrada_con_unidades(float, textounidad="g/mol")
        self.M.valueChanged.connect(partial(self.changeParams, "M"))
        layout.addWidget(self.M, 13, 1)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Boiling point")), 14, 0)
        self.Tb = Entrada_con_unidades(Temperature)
        self.Tb.valueChanged.connect(partial(self.changeParams, "Tb"))
        layout.addWidget(self.Tb, 14, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Specific Gravity")), 15, 0)
        self.SG = Entrada_con_unidades(float)
        self.SG.valueChanged.connect(partial(self.changeParams, "SG"))
        layout.addWidget(self.SG, 15, 1)

        newComponent.loadUI(self)

        func = {"Constantinou": Constantinou_Gani,
                "Wilson": Wilson_Jasperson,
                "Joback": Joback,
                "Ambrose": Ambrose,
                "Elliott": Elliott,
                "Marrero": Marrero_Pardillo}
        self.unknown = func[self.metodo]()

        for i, nombre in enumerate(self.unknown.coeff["txt"]):
            self.TablaContribuciones.addItem(nombre[0])

        if metodo in ["Constantinou", "Wilson"]:
            self.Order()

    def Order(self):
        """Show/Hide group of undesired order"""
        for i in range(self.unknown.FirstOrder):
            self.TablaContribuciones.item(i).setHidden(self.Order2.isChecked())
        for i in range(self.unknown.FirstOrder, self.unknown.SecondOrder):
            self.TablaContribuciones.item(i).setHidden(self.Order1.isChecked())

    def borrar(self, indice=None):
        """Remove some group contribution from list"""
        if not indice:
            indice = self.Grupos.currentRow()
        if indice != -1:
            self.Grupos.removeRow(indice)
            del self.grupo[indice]
            del self.indices[indice]
            del self.contribucion[indice]
            self.calculo(**{"group": self.indices,
                            "contribution": self.contribucion})

    def clear(self):
        """Clear widgets from dialog"""
        self.Grupos.clearContents()
        self.Grupos.setRowCount(0)
        self.grupo = []
        self.indices = []
        self.contribucion = []
        self.Formula.clear()
        self.M.clear()
        self.nombre.clear()
        self.Tb.clear()
        self.SG.clear()
        self.unknown.clear()
        self.status.setState(self.unknown.status, self.unknown.msg)

    def cellChanged(self, i, j):
        if j == 0:
            valor = int(self.Grupos.item(i, j).text())
            if valor <= 0:
                self.borrar(i)
            else:
                self.contribucion[i] = int(valor)
        self.calculo(**{"group": self.indices, "contribution": self.contribucion})

    def selectedChanged(self, i):
        self.botonAdd.setEnabled(i != -1)

    def add(self):
        indice = self.Grupos.rowCount()
        grupo = self.TablaContribuciones.currentItem().text()
        if grupo not in self.grupo:
            self.grupo.append(grupo)
            self.indices.append(self.TablaContribuciones.currentRow())
            self.contribucion.append(1)
            self.Grupos.setRowCount(indice+1)
            self.Grupos.setItem(indice, 0, QtWidgets.QTableWidgetItem("1"))
            self.Grupos.item(indice, 0).setTextAlignment(
                QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
            self.Grupos.setItem(indice, 1, QtWidgets.QTableWidgetItem(grupo))
            self.Grupos.item(indice, 1).setFlags(
                QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsEnabled)
            self.Grupos.setRowHeight(indice, 20)
        else:
            indice = self.grupo.index(grupo)
            self.contribucion[indice] += 1
            self.Grupos.item(indice, 0).setText(str(int(
                self.Grupos.item(indice, 0).text())+1))
        self.calculo(**{"group": self.indices, "contribution": self.contribucion})

    def calculo(self, **kwargs):
        """Calculate function"""
        newComponent.calculo(self, **kwargs)
        if self.unknown.status:
            self.Formula.setText(self.unknown.formula)

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = Definicion_Petro()
    # Dialog = Ui_Contribution("Ambrose")
    Dialog.show()
    sys.exit(app.exec_())
