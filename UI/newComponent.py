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
# Module to implement new component by group contribution methods UI
#   -newComponent: Main dialog class with common functionality
#   -UI_Contribution: Definition for group contribution
#   -View_Contribution: Dialog to show the properties of a group contribution
###############################################################################


import os
from functools import partial

from PyQt5 import QtCore, QtGui, QtWidgets

from lib import sql
from lib.config import IMAGE_PATH
from lib.newComponent import _methods, _group2Order
from lib.unidades import (Temperature, Pressure, SpecificVolume, Enthalpy,
                          SolubilityParameter)
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
    """Dialog to show the properties of a group contribution"""
    def __init__(self, cmp=None, parent=None):
        """Constructor
        cmp: optional new component to show the properties"""
        super(View_Contribution, self).__init__(parent)
        self.setWindowTitle(QtWidgets.QApplication.translate(
            "pychemqt", "Group Contribution new component"))
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
        self.Tb = Entrada_con_unidades(Temperature)
        layout.addWidget(self.Tb, 3, 2)
        label = QtWidgets.QLabel("Tm")
        label.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Melting Temperature"))
        layout.addWidget(label, 4, 1)
        self.Tf = Entrada_con_unidades(Temperature)
        layout.addWidget(self.Tf, 4, 2)
        layout.addWidget(QtWidgets.QLabel("Tc"), 5, 1)
        self.Tc = Entrada_con_unidades(Temperature)
        layout.addWidget(self.Tc, 5, 2)
        layout.addWidget(QtWidgets.QLabel("Pc"), 6, 1)
        self.Pc = Entrada_con_unidades(Pressure)
        layout.addWidget(self.Pc, 6, 2)
        layout.addWidget(QtWidgets.QLabel("Vc"), 7, 1)
        self.Vc = Entrada_con_unidades(SpecificVolume)
        layout.addWidget(self.Vc, 7, 2)

        label = QtWidgets.QLabel("ΔHf")
        label.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Enthalpy of formation of ideal gas"))
        layout.addWidget(label, 8, 1)
        self.Hf = Entrada_con_unidades(Enthalpy)
        layout.addWidget(self.Hf, 8, 2)
        label = QtWidgets.QLabel("ΔGf")
        label.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Gibbs free energy of formation of ideal gas"))
        layout.addWidget(label, 9, 1)
        self.Gf = Entrada_con_unidades(Enthalpy)
        layout.addWidget(self.Gf, 9, 2)

        label = QtWidgets.QLabel("ΔHf")
        label.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Enthalpy of fusion"))
        layout.addWidget(label, 2, 4)
        self.Hm = Entrada_con_unidades(Enthalpy)
        layout.addWidget(self.Hm, 2, 5)
        label = QtWidgets.QLabel("ΔHv")
        label.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Enthalpy of vaporization"))
        layout.addWidget(label, 3, 4)
        self.Hv = Entrada_con_unidades(Enthalpy)
        layout.addWidget(self.Hv, 3, 5)
        layout.addWidget(QtWidgets.QLabel("Cpa"), 4, 4)
        self.cpa = Entrada_con_unidades(float)
        layout.addWidget(self.cpa, 4, 5)
        layout.addWidget(QtWidgets.QLabel("Cpb"), 5, 4)
        self.cpb = Entrada_con_unidades(float)
        layout.addWidget(self.cpb, 5, 5)
        layout.addWidget(QtWidgets.QLabel("Cpc"), 6, 4)
        self.cpc = Entrada_con_unidades(float)
        layout.addWidget(self.cpc, 6, 5)
        layout.addWidget(QtWidgets.QLabel("Cpd"), 7, 4)
        self.cpd = Entrada_con_unidades(float)
        layout.addWidget(self.cpd, 7, 5)
        layout.addWidget(QtWidgets.QLabel("μa"), 8, 4)
        self.mua = Entrada_con_unidades(float)
        layout.addWidget(self.mua, 8, 5)
        layout.addWidget(QtWidgets.QLabel("μb"), 9, 4)
        self.mub = Entrada_con_unidades(float)
        layout.addWidget(self.mub, 9, 5)

        label = QtWidgets.QLabel("SG")
        label.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Specific gravity at 60ºF"))
        layout.addWidget(label, 2, 7)
        self.gravity = Entrada_con_unidades(float)
        layout.addWidget(self.gravity, 2, 8)
        label = QtWidgets.QLabel("API")
        label.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "API Specific gravity"))
        layout.addWidget(label, 3, 7)
        self.API = Entrada_con_unidades(float)
        layout.addWidget(self.API, 3, 8)
        label = QtWidgets.QLabel("Kw")
        label.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Watson characterization factor"))
        layout.addWidget(label, 4, 7)
        self.watson = Entrada_con_unidades(float)
        layout.addWidget(self.watson, 4, 8)

        label = QtWidgets.QLabel("w")
        label.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Acentric factor"))
        layout.addWidget(label, 5, 7)
        self.f_acent = Entrada_con_unidades(float)
        layout.addWidget(self.f_acent, 5, 8)
        label = QtWidgets.QLabel("Ra")
        label.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Rackett constant"))
        layout.addWidget(label, 6, 7)
        self.rackett = Entrada_con_unidades(float)
        layout.addWidget(self.rackett, 6, 8)
        label = QtWidgets.QLabel("Vliq")
        label.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Volume Liquid Constant"))
        layout.addWidget(label, 7, 7)
        self.Vliq = Entrada_con_unidades(float)
        layout.addWidget(self.Vliq, 7, 8)
        label = QtWidgets.QLabel("Sol")
        label.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Solubility Parameter"))
        layout.addWidget(label, 8, 7)
        self.Parametro_solubilidad = Entrada_con_unidades(SolubilityParameter)
        layout.addWidget(self.Parametro_solubilidad, 8, 8)

        layout.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 15, 8)
        btn = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Close)
        btn.rejected.connect(self.accept)
        layout.addWidget(btn, 16, 1, 1, 8)

        self.setReadOnly(True)
        if cmp:
            self.rellenar(cmp)

    def setReadOnly(self, bool):
        self.M.setReadOnly(bool)
        self.Tb.setReadOnly(bool)
        self.Tf.setReadOnly(bool)
        self.Tc.setReadOnly(bool)
        self.Pc.setReadOnly(bool)
        self.Vc.setReadOnly(bool)
        self.Hf.setReadOnly(bool)
        self.Gf.setReadOnly(bool)

        self.Hm.setReadOnly(bool)
        self.Hv.setReadOnly(bool)
        self.cpa.setReadOnly(bool)
        self.cpb.setReadOnly(bool)
        self.cpc.setReadOnly(bool)
        self.cpd.setReadOnly(bool)
        self.mua.setReadOnly(bool)
        self.mub.setReadOnly(bool)

        self.gravity.setReadOnly(bool)
        self.API.setReadOnly(bool)
        self.watson.setReadOnly(bool)
        self.f_acent.setReadOnly(bool)
        self.rackett.setReadOnly(bool)
        self.Vliq.setReadOnly(bool)
        self.Parametro_solubilidad.setReadOnly(bool)

    def rellenar(self, cmp):
        self.nombre.setText(cmp.name)
        self.M.setValue(cmp.M)
        self.Tb.setValue(cmp.Tb)
        self.Tf.setValue(cmp.Tf)
        self.Tc.setValue(cmp.Tc)
        self.Pc.setValue(cmp.Pc)
        self.Vc.setValue(cmp.Vc)
        self.Hf.setValue(cmp.Hf)
        self.Gf.setValue(cmp.Gf)

        self.Hm.setValue(cmp.Hm)
        self.Hv.setValue(cmp.Hv)
        self.cpa.setValue(cmp.cp[0])
        self.cpb.setValue(cmp.cp[1])
        self.cpc.setValue(cmp.cp[2])
        self.cpd.setValue(cmp.cp[3])

        self.gravity.setValue(cmp.SG)
        self.API.setValue(cmp.API)
        self.watson.setValue(cmp.Kw)
        self.f_acent.setValue(cmp.f_acent)
        self.rackett.setValue(cmp.rackett)
        self.Vliq.setValue(cmp.Vliq)
        self.Parametro_solubilidad.setValue(cmp.Parametro_solubilidad)


class Ui_Contribution(newComponent):
    """Dialog to define hypotethical new component with several group
    contribucion methods"""
    ViewDetails = View_Contribution

    def __init__(self, metodo, parent=None):
        """Metodo: name of group contribution method:
            Joback, Constantinou, Wilson, Marrero, Elliott, Ambrose, Klincewicz
        """
        super(Ui_Contribution, self).__init__(parent)

        self.grupo = []
        self.indices = []
        self.contribucion = []
        self.metodo = metodo
        lyt = QtWidgets.QVBoxLayout(self)
        widget = QtWidgets.QWidget()
        layout = QtWidgets.QGridLayout(widget)

        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Name")), 1, 0)
        self.nombre = QtWidgets.QLineEdit()
        self.nombre.textChanged.connect(partial(self.changeParams, "name"))
        layout.addWidget(self.nombre, 1, 1, 1, 3)

        self.Group = QtWidgets.QTableWidget()
        self.Group.verticalHeader().hide()
        self.Group.setRowCount(0)
        self.Group.setColumnCount(2)
        self.Group.setHorizontalHeaderItem(0, QtWidgets.QTableWidgetItem("Nk"))
        self.Group.setHorizontalHeaderItem(1, QtWidgets.QTableWidgetItem(
            QtWidgets.QApplication.translate("pychemqt", "Group")))
        self.Group.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.Group.setSortingEnabled(True)
        self.Group.horizontalHeader().setStretchLastSection(True)
        self.Group.setColumnWidth(0, 50)
        self.Group.setItemDelegateForColumn(0, SpinEditor(self))
        self.Group.cellChanged.connect(self.cellChanged)
        self.Group.setEditTriggers(QtWidgets.QAbstractItemView.AllEditTriggers)
        layout.addWidget(self.Group, 2, 0, 3, 3)

        self.Formula = QtWidgets.QLabel()
        font = QtGui.QFont()
        font.setPointSize(12)
        self.Formula.setFont(font)
        self.Formula.setAlignment(
                QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
        self.Formula.setFixedHeight(50)
        layout.addWidget(self.Formula, 2, 3)
        self.btnDelete = QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(
            os.path.join(IMAGE_PATH, "button", "editDelete.png"))),
            QtWidgets.QApplication.translate("pychemqt", "Delete"))
        self.btnDelete.clicked.connect(self.borrar)
        layout.addWidget(self.btnDelete, 3, 3)
        self.btnClear = QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(
            os.path.join(IMAGE_PATH, "button", "clear.png"))),
            QtWidgets.QApplication.translate("pychemqt", "Clear"))
        self.btnClear.clicked.connect(self.clear)
        layout.addWidget(self.btnClear, 4, 3)

        self.line = QtWidgets.QFrame()
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        layout.addWidget(self.line, 5, 0, 1, 4)

        self.groupContributions = QtWidgets.QListWidget()
        self.groupContributions.currentItemChanged.connect(self.selectChanged)
        self.groupContributions.itemDoubleClicked.connect(self.add)
        layout.addWidget(self.groupContributions, 6, 0, 7, 3)
        self.btnAdd = QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(
            os.path.join(IMAGE_PATH, "button", "add.png"))),
            QtWidgets.QApplication.translate("pychemqt", "Add"))
        self.btnAdd.setDisabled(True)
        self.btnAdd.clicked.connect(self.add)
        layout.addWidget(self.btnAdd, 6, 3)
        layout.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 7, 1)

        # Show selection for method with several order contributions
        if metodo in _group2Order:
            self.Order1 = QtWidgets.QRadioButton(
                QtWidgets.QApplication.translate("pychemqt", "1st order"))
            self.Order1.setChecked(True)
            self.Order1.toggled.connect(self.Order)
            layout.addWidget(self.Order1, 10, 3)
            self.Order2 = QtWidgets.QRadioButton(
                QtWidgets.QApplication.translate("pychemqt", "2nd order"))
            self.Order2.toggled.connect(self.Order)
            layout.addWidget(self.Order2, 11, 3)

            if metodo == "MarreroGani":
                self.Order3 = QtWidgets.QRadioButton(
                    QtWidgets.QApplication.translate("pychemqt", "3rd order"))
                layout.addWidget(self.Order3, 12, 3)
                self.Order3.toggled.connect(self.Order)

        # layout.addItem(QtWidgets.QSpacerItem(
            # 10, 10, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            # 12, 1)
        labelTb = QtWidgets.QLabel("Tb")
        labelTb.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Experimental Boiling Temperature"))
        layout.addWidget(labelTb, 13, 0)
        self.Tb = Entrada_con_unidades(Temperature)
        self.Tb.valueChanged.connect(partial(self.changeParams, "Tb"))
        layout.addWidget(self.Tb, 13, 1)
        labelM = QtWidgets.QLabel("M")
        labelM.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Molecular Weight"))
        layout.addWidget(labelM, 14, 0)
        self.M = Entrada_con_unidades(float, textounidad="g/mol")
        self.M.valueChanged.connect(partial(self.changeParams, "M"))
        layout.addWidget(self.M, 14, 1)
        label = QtWidgets.QLabel("SG")
        label.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Experimental Specific Gravity at 60ºF"))
        layout.addWidget(label, 15, 0)
        self.SG = Entrada_con_unidades(float)
        self.SG.valueChanged.connect(partial(self.changeParams, "SG"))
        layout.addWidget(self.SG, 15, 1)
        lyt.addWidget(widget)

        # Show widget for specific method
        if metodo == "Constantinou":
            # Disable Tb as input parameter
            labelTb.setEnabled(False)
            self.Tb.setEnabled(False)

        elif metodo == "Wilson":
            self.Tb.setResaltado(True)
            layout.addWidget(QtWidgets.QLabel(
                QtWidgets.QApplication.translate("pychemqt", "Rings")), 16, 0)
            self.ring = QtWidgets.QSpinBox()
            self.ring.setFixedWidth(80)
            self.ring.valueChanged.connect(partial(self.changeParams, "ring"))
            layout.addWidget(self.ring, 16, 1)

        elif metodo == "Elliott":
            self.M.setResaltado(True)

        elif metodo == "Ambrose":
            self.Tb.setResaltado(True)
            layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
                "pychemqt", "Platt number")), 16, 0)
            self.plat = QtWidgets.QSpinBox()
            self.plat.setFixedWidth(80)
            tr = "The Platt number is the number of pairs of carbon atoms "
            tr += "which are separated by three carbon-carbon bonds and is an "
            tr += "indicator of the degree of branching in the molecule. The "
            tr += "Platt number of an n-alkane is equal to the number of "
            tr += "carbons minus three"
            self.plat.setToolTip(QtWidgets.QApplication.translate(
                "pychemqt", tr))
            self.plat.valueChanged.connect(partial(self.changeParams, "platt"))
            layout.addWidget(self.plat, 16, 1)

        elif metodo == "Klincewicz":
            self.Tb.setResaltado(True)

            self.nogroupKlincewicz = QtWidgets.QCheckBox(
                    QtWidgets.QApplication.translate(
                        "pychemqt",
                        "Use the simple no group contribution method"))
            self.nogroupKlincewicz.toggled.connect(self.nogroupCheckToggled)
            layout.addWidget(self.nogroupKlincewicz, 16, 0, 1, 4)
            layout.addWidget(QtWidgets.QLabel(
                QtWidgets.QApplication.translate("pychemqt", "Atoms")), 17, 0)
            self.atoms = QtWidgets.QSpinBox()
            self.atoms.setFixedWidth(80)
            self.atoms.valueChanged.connect(
                    partial(self.changeParams, "atoms"))
            layout.addWidget(self.atoms, 17, 1)

        newComponent.loadUI(self)

        # Initialization variables
        func = {}
        for f in _methods:
            func[f.__name__] = f
        self.unknown = func[self.metodo]()

        title = self.unknown.__title__
        title += " " + QtWidgets.QApplication.translate(
                "pychemqt", "new component definition")
        self.setWindowTitle(title)

        for i, nombre in enumerate(self.unknown.coeff["txt"]):
            self.groupContributions.addItem(nombre[0])

        if metodo in _group2Order:
            self.Order()
        if metodo == "Klincewicz":
            self.nogroupCheckToggled(False)

    def Order(self):
        """Show/Hide group of undesired order"""
        for i in range(self.unknown.FirstOrder):
            item = self.groupContributions.item(i)
            item.setHidden(not self.Order1.isChecked())
        for i in range(self.unknown.FirstOrder, self.unknown.SecondOrder):
            item = self.groupContributions.item(i)
            item.setHidden(not self.Order2.isChecked())
        if self.unknown.ThirdOrder:
            for i in range(self.unknown.SecondOrder, self.unknown.ThirdOrder):
                item = self.groupContributions.item(i)
                item.setHidden(not self.Order3.isChecked())

    def nogroupCheckToggled(self, boolean):
        """Set advanced properties input status for Klincewitcz method"""
        self.changeParams("nogroup", boolean)
        self.M.setResaltado(boolean)
        self.atoms.setEnabled(boolean)
        self.Group.setDisabled(boolean)
        self.Formula.setDisabled(boolean)
        self.btnDelete.setDisabled(boolean)
        self.btnClear.setDisabled(boolean)
        self.btnAdd.setDisabled(boolean)
        self.groupContributions.setDisabled(boolean)

    def borrar(self, indice=None):
        """Remove some group contribution from added group list"""
        if not indice:
            indice = self.Group.currentRow()
        if indice != -1:
            self.Group.removeRow(indice)
            del self.grupo[indice]
            del self.indices[indice]
            del self.contribucion[indice]
            self.calculo(**{"group": self.indices,
                            "contribution": self.contribucion})

    def clear(self):
        """Clear widgets from dialog"""
        self.Group.clearContents()
        self.Group.setRowCount(0)
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
        """Process the user manual count of group contribution changed"""
        if j == 0:
            valor = int(self.Group.item(i, j).text())
            if valor <= 0:
                self.borrar(i)
            else:
                self.contribucion[i] = int(valor)
        kw = {"group": self.indices, "contribution": self.contribucion}
        self.calculo(**kw)

    def selectChanged(self, i):
        """The add button is only enabled when the group list have any selected
        row"""
        self.btnAdd.setEnabled(i != -1)

    def add(self):
        """Add the current selected item to the list of group"""
        indice = self.Group.rowCount()
        grupo = self.groupContributions.currentItem().text()
        if grupo not in self.grupo:
            self.grupo.append(grupo)
            self.indices.append(self.groupContributions.currentRow())
            self.contribucion.append(1)
            self.Group.setRowCount(indice+1)
            self.Group.setItem(indice, 0, QtWidgets.QTableWidgetItem("1"))
            self.Group.item(indice, 0).setTextAlignment(
                QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
            self.Group.setItem(indice, 1, QtWidgets.QTableWidgetItem(grupo))
            self.Group.item(indice, 1).setFlags(
                QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsEnabled)
            self.Group.setRowHeight(indice, 20)
        else:
            indice = self.grupo.index(grupo)
            self.contribucion[indice] += 1
            self.Group.item(indice, 0).setText(str(int(
                self.Group.item(indice, 0).text())+1))
        kw = {"group": self.indices, "contribution": self.contribucion}
        self.calculo(**kw)

    def calculo(self, **kwargs):
        """Calculate function"""
        newComponent.calculo(self, **kwargs)
        if self.unknown.status in (1, 3):
            self.Formula.setText(self.unknown.formula)

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = Ui_Contribution("MarreroGani")
    Dialog.show()
    sys.exit(app.exec_())
