#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Module to implement new component
#   -newComponent: Main dialog class with common functionality
#   -UI_Contribution: Definition for group contribution
#   -Definicion_Petro: Definition of crude and oil fraction
###############################################################################

import os
from functools import partial

from PyQt4 import QtCore, QtGui

from UI.widgets import Entrada_con_unidades, Tabla, Status
from UI.delegate import SpinEditor
from UI.viewComponents import View_Component, View_Petro, View_Contribution
from lib.unidades import Temperature, Pressure, Diffusivity
from lib.petro import Crudo, crudo, Petroleo
from lib.compuestos import (Joback, Constantinou_Gani, Wilson_Jasperson,
                            Marrero_Pardillo, Elliott, Ambrose)
from lib import sql


class newComponent(QtGui.QDialog):
    """Main dialog class with common functionality"""

    def loadUI(self):
        """Define common widget for chid class"""
        layoutBottom = QtGui.QHBoxLayout()
        self.status = Status()
        layoutBottom.addWidget(self.status)
        self.buttonShowDetails = QtGui.QPushButton(
            QtGui.QApplication.translate("pychemqt", "Show Details"))
        self.buttonShowDetails.clicked.connect(self.showDetails)
        self.buttonShowDetails.setEnabled(False)
        layoutBottom.addWidget(self.buttonShowDetails)
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Cancel |
                                                QtGui.QDialogButtonBox.Save)
        self.buttonBox.button(QtGui.QDialogButtonBox.Save).setEnabled(False)
        self.buttonBox.accepted.connect(self.save)
        self.buttonBox.rejected.connect(self.reject)
        layoutBottom.addWidget(self.buttonBox)
        self.layout().addLayout(layoutBottom, 30, 0, 1, 6)

    def save(self):
        """Save new componente in user database"""
        elemento = self.unknown.export2Component()
        sql.inserElementsFromArray(sql.databank_Custom_name, [elemento])
        Dialog = View_Component(1001+sql.N_comp_Custom)
        Dialog.show()
        QtGui.QDialog.accept(self)

    def changeParams(self, parametro, valor):
        self.calculo(**{parametro: valor})

    def calculo(self, **kwargs):
        self.status.setState(4)
        self.unknown(**kwargs)
        self.status.setState(self.unknown.status, self.unknown.msg)
        self.buttonShowDetails.setEnabled(self.unknown.status)
        self.buttonBox.button(QtGui.QDialogButtonBox.Save).setEnabled(
            self.unknown.status)

    def showDetails(self):
        """Show details of new component"""
        dialog = self.ViewDetails(self.unknown)
        dialog.exec_()


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
        self.setWindowTitle(QtGui.QApplication.translate(
            "pychemqt", "Select the component group for method") +" "+ metodo)

        self.grupo = []
        self.indices = []
        self.contribucion = []
        self.metodo = metodo
        layout = QtGui.QGridLayout(self)
        self.Grupos = QtGui.QTableWidget()
        self.Grupos.verticalHeader().hide()
        self.Grupos.setRowCount(0)
        self.Grupos.setColumnCount(2)
        self.Grupos.setHorizontalHeaderItem(0, QtGui.QTableWidgetItem("Nk"))
        self.Grupos.setHorizontalHeaderItem(1, QtGui.QTableWidgetItem(
            QtGui.QApplication.translate("pychemqt", "Group")))
        self.Grupos.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
        self.Grupos.setSortingEnabled(True)
        self.Grupos.horizontalHeader().setStretchLastSection(True)
        self.Grupos.setColumnWidth(0, 50)
        self.Grupos.setItemDelegateForColumn(0, SpinEditor(self))
        self.Grupos.cellChanged.connect(self.cellChanged)
        self.Grupos.setEditTriggers(QtGui.QAbstractItemView.AllEditTriggers)
        layout.addWidget(self.Grupos, 0, 0, 3, 3)

        self.Formula = QtGui.QLabel()
        font = QtGui.QFont()
        font.setPointSize(12)
        self.Formula.setFont(font)
        self.Formula.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
        self.Formula.setFixedHeight(50)
        layout.addWidget(self.Formula, 0, 3)
        self.botonBorrar = QtGui.QPushButton(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"]+"/images/button/editDelete.png")),
            QtGui.QApplication.translate("pychemqt", "Delete"))
        self.botonBorrar.clicked.connect(self.borrar)
        layout.addWidget(self.botonBorrar, 1, 3)
        self.botonClear = QtGui.QPushButton(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"]+"/images/button/clear.png")),
            QtGui.QApplication.translate("pychemqt", "Clear"))
        self.botonClear.clicked.connect(self.clear)
        layout.addWidget(self.botonClear, 2, 3)

        self.line = QtGui.QFrame()
        self.line.setFrameShape(QtGui.QFrame.HLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        layout.addWidget(self.line, 3, 0, 1, 4)

        self.TablaContribuciones = QtGui.QListWidget()
        self.TablaContribuciones.currentItemChanged.connect(self.selectedChanged)
        self.TablaContribuciones.itemDoubleClicked.connect(self.add)
        layout.addWidget(self.TablaContribuciones, 4, 0, 7, 3)
        self.botonAdd = QtGui.QPushButton(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"]+"/images/button/add.png")),
            QtGui.QApplication.translate("pychemqt", "Add"))
        self.botonAdd.setDisabled(True)
        self.botonAdd.clicked.connect(self.add)
        layout.addWidget(self.botonAdd, 4, 3)
        layout.addItem(QtGui.QSpacerItem(20, 20, QtGui.QSizePolicy.Expanding,
                                         QtGui.QSizePolicy.Expanding), 5, 1, 1, 1)

        # Show widget for specific method
        if metodo in ["Constantinou", "Wilson"]:
            self.Order1 = QtGui.QRadioButton(
                QtGui.QApplication.translate("pychemqt", "1st order"))
            self.Order1.setChecked(True)
            self.Order1.toggled.connect(self.Order)
            layout.addWidget(self.Order1, 6, 3)
            self.Order2 = QtGui.QRadioButton(
                QtGui.QApplication.translate("pychemqt", "2nd order"))
            layout.addWidget(self.Order2, 7, 3)

        if metodo == "Wilson":
            layout.addWidget(QtGui.QLabel(
                QtGui.QApplication.translate("pychemqt", "Rings")), 8, 3)
            self.anillos = QtGui.QSpinBox()
            self.anillos.valueChanged.connect(partial(self.changeParams, "ring"))
            layout.addWidget(self.anillos, 9, 3)

        if metodo == "Marrero":
            layout.addWidget(QtGui.QLabel(
                QtGui.QApplication.translate("pychemqt", "Atoms")), 8, 3)
            self.Atomos = QtGui.QSpinBox()
            self.Atomos.valueChanged.connect(partial(self.changeParams, "atomos"))
            layout.addWidget(self.Atomos, 9, 3)

        if metodo == "Ambrose":
            layout.addWidget(QtGui.QLabel(
                QtGui.QApplication.translate("pychemqt", "Platt number")), 8, 3)
            self.Platt = QtGui.QSpinBox()
            self.Platt.setToolTip(QtGui.QApplication.translate(
                "pychemqt", "The Platt number is the number of pairs of carbon \
                atoms which are separated \nby three carbon-carbon bonds and \
                is an indicator of the degree of branching in the molecule.\n\
                The Platt number of an n-alkane is equal to the number of \
                carbons minus three"))
            self.Platt.valueChanged.connect(partial(self.changeParams, "platt"))
            layout.addWidget(self.Platt, 9, 3)

        layout.addItem(QtGui.QSpacerItem(20, 20, QtGui.QSizePolicy.Expanding,
                                         QtGui.QSizePolicy.Expanding), 10, 1, 1, 1)
        layout.addItem(QtGui.QSpacerItem(20, 20, QtGui.QSizePolicy.Expanding,
                                         QtGui.QSizePolicy.Expanding), 11, 0, 1, 4)
        layout.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Name")), 12, 0)
        self.nombre = QtGui.QLineEdit()
        self.nombre.textChanged.connect(partial(self.changeParams, "name"))
        layout.addWidget(self.nombre, 12, 1, 1, 3)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Molecular Weight")), 13, 0)
        self.M = Entrada_con_unidades(float, textounidad="g/mol")
        self.M.valueChanged.connect(partial(self.changeParams, "M"))
        layout.addWidget(self.M, 13, 1)
        layout.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Boiling point")), 14, 0)
        self.Tb = Entrada_con_unidades(Temperature)
        self.Tb.valueChanged.connect(partial(self.changeParams, "Tb"))
        layout.addWidget(self.Tb, 14, 1)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
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
            self.Grupos.setItem(indice, 0, QtGui.QTableWidgetItem("1"))
            self.Grupos.item(indice, 0).setTextAlignment(
                QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
            self.Grupos.setItem(indice, 1, QtGui.QTableWidgetItem(grupo))
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


class Definicion_Petro(newComponent):
    """Dialog for define hypothetical crude and oil fraction"""
    ViewDetails = View_Petro

    def __init__(self, parent=None):
        super(Definicion_Petro, self).__init__(parent)
        self.setWindowTitle(QtGui.QApplication.translate(
            "pychemqt", "Petrol component definition"))
        layout = QtGui.QGridLayout(self)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Name")), 1, 1)
        self.nombre = QtGui.QLineEdit()
        self.nombre.textChanged.connect(partial(self.changeParams, "name"))
        layout.addWidget(self.nombre, 1, 2, 1, 4)
        layout.addItem(QtGui.QSpacerItem(10, 10, QtGui.QSizePolicy.Fixed,
                                         QtGui.QSizePolicy.Fixed), 2, 1, 1, 2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Boiling point")), 3, 1)
        self.Tb = Entrada_con_unidades(Temperature)
        self.Tb.valueChanged.connect(partial(self.changeParams, "Tb"))
        layout.addWidget(self.Tb, 3, 2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Molecular Weight")), 4, 1)
        self.M = Entrada_con_unidades(float, textounidad="g/mol")
        self.M.valueChanged.connect(partial(self.changeParams, "M"))
        layout.addWidget(self.M, 4, 2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Specific Gravity")), 5, 1)
        self.SG = Entrada_con_unidades(float)
        self.SG.valueChanged.connect(partial(self.changeParams, "SG"))
        layout.addWidget(self.SG, 5, 2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "API gravity")), 6, 1)
        self.API = Entrada_con_unidades(float)
        self.API.valueChanged.connect(partial(self.changeParams, "API"))
        layout.addWidget(self.API, 6, 2)
        layout.addWidget(QtGui.QLabel("K watson:"), 7, 1)
        self.Kw = Entrada_con_unidades(float)
        self.Kw.valueChanged.connect(partial(self.changeParams, "Kw"))
        layout.addWidget(self.Kw, 7, 2)
        layout.addWidget(QtGui.QLabel("C/H:"), 8, 1)
        self.CH = Entrada_con_unidades(float)
        self.CH.valueChanged.connect(partial(self.changeParams, "CH"))
        layout.addWidget(self.CH, 8, 2)
        layout.addWidget(QtGui.QLabel(u"ν<sub>100F</sub>:"), 9, 1)
        self.v100 = Entrada_con_unidades(Diffusivity)
        self.v100.valueChanged.connect(partial(self.changeParams, "v100"))
        layout.addWidget(self.v100, 9, 2)
        layout.addWidget(QtGui.QLabel(u"ν<sub>210F</sub>:"), 10, 1)
        self.v210 = Entrada_con_unidades(Diffusivity)
        self.v210.valueChanged.connect(partial(self.changeParams, "v210"))
        layout.addWidget(self.v210, 10, 2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Refraction index")), 11, 1)
        self.n = Entrada_con_unidades(float)
        self.n.valueChanged.connect(partial(self.changeParams, "n"))
        layout.addWidget(self.n, 11, 2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Huang parameter")), 12, 1)
        self.I = Entrada_con_unidades(float)
        self.I.valueChanged.connect(partial(self.changeParams, "I"))
        layout.addWidget(self.I, 12, 2)
        layout.addWidget(QtGui.QLabel("%S"), 13, 1)
        self.S = Entrada_con_unidades(float, spinbox=True, step=1.0, max=100)
        self.S.valueChanged.connect(partial(self.changeParams, "S"))
        layout.addWidget(self.S, 13, 2)
        layout.addWidget(QtGui.QLabel("%H"), 14, 1)
        self.H = Entrada_con_unidades(float, spinbox=True, step=1.0, max=100)
        self.H.valueChanged.connect(partial(self.changeParams, "H"))
        layout.addWidget(self.H, 14, 2)
        layout.addWidget(QtGui.QLabel("%N"), 15, 1)
        self.N = Entrada_con_unidades(float, spinbox=True, step=1.0, max=100)
        self.N.valueChanged.connect(partial(self.changeParams, "N"))
        layout.addWidget(self.N, 15, 2)

        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Carbons number")), 19, 1)
        self.carbonos = Entrada_con_unidades(int, width=50, spinbox=True,
                                             step=1, start=7, min=5, max=100)
        self.N.valueChanged.connect(partial(self.changeParams, "Nc"))
        layout.addWidget(self.carbonos, 19, 2)

        layout.addItem(QtGui.QSpacerItem(10, 10, QtGui.QSizePolicy.Expanding,
                                         QtGui.QSizePolicy.Expanding), 3, 3, 15, 1)
        self.checkCurva = QtGui.QCheckBox(QtGui.QApplication.translate(
            "pychemqt", "Define destillation curve"))
        layout.addWidget(self.checkCurva, 3, 4, 1, 2)
        self.tipoCurva = QtGui.QComboBox()
        self.tipoCurva.addItem("ASTM D86")
        self.tipoCurva.addItem("TBP")
        self.tipoCurva.addItem("EFV")
        self.tipoCurva.addItem("ASTM D1186")
        self.tipoCurva.addItem("ASTM D2887 (SD)")
        self.tipoCurva.setEnabled(False)
        layout.addWidget(self.tipoCurva, 4, 4, 1, 2)
        self.textoPresion = QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Pressure"))
        self.textoPresion.setEnabled(False)
        layout.addWidget(self.textoPresion, 5, 4)
        self.presion = Entrada_con_unidades(Pressure, value=101325.)
        self.presion.setEnabled(False)
        self.presion.valueChanged.connect(partial(self.changeParams, "P_dist"))
        layout.addWidget(self.presion, 5, 5)
        self.curvaDestilacion = Tabla(
            2, filas=1, horizontalHeader=["%dist", "Tb, "+Temperature.text()],
            verticalHeader=False, dinamica=True)
        self.curvaDestilacion.setEnabled(False)
        self.curvaDestilacion.editingFinished.connect(self.changeCurva)
        layout.addWidget(self.curvaDestilacion, 6, 4, 13, 2)
        self.checkBlend = QtGui.QCheckBox(QtGui.QApplication.translate(
            "pychemqt", "Blend if its necessary"))
        layout.addWidget(self.checkBlend, 19, 4, 1, 2)

        self.checkCurva.toggled.connect(self.tipoCurva.setEnabled)
        self.checkCurva.toggled.connect(self.presion.setEnabled)
        self.checkCurva.toggled.connect(self.textoPresion.setEnabled)
        self.checkCurva.toggled.connect(self.curvaDestilacion.setEnabled)

        layout.addItem(QtGui.QSpacerItem(20, 20, QtGui.QSizePolicy.Fixed,
                                         QtGui.QSizePolicy.Fixed), 20, 1, 1, 2)
        self.checkCrudo = QtGui.QCheckBox(QtGui.QApplication.translate(
            "pychemqt", "Use petrol fraction from list"))
        self.checkCrudo.toggled.connect(self.changeUnknown)
        layout.addWidget(self.checkCrudo, 21, 1, 1, 2)
        self.crudo = QtGui.QComboBox()
        self.crudo.setEnabled(False)
        self.crudo.addItem("")
        for i in crudo[1:]:
            self.crudo.addItem("%s (%s)  API: %s %S: %s" % (i[0], i[1], i[3], i[4]))
            # i[0]+" ("+i[1]+")"+"   API: "+str(i[3])+"  %S: "+str(i[4]))
        self.crudo.currentIndexChanged.connect(partial(
            self.changeParams, "indice"))
        layout.addWidget(self.crudo, 23, 1, 1, 5)
        layout.addWidget(QtGui.QLabel("Pseudo C+"), 24, 1)
        self.Cplus = Entrada_con_unidades(int, width=50, spinbox=True, step=1,
                                          min=6)
        self.Cplus.valueChanged.connect(partial(self.changeParams, "Cplus"))
        layout.addWidget(self.Cplus, 24, 2)
        self.checkCrudo.toggled.connect(self.crudo.setEnabled)
        self.checkCrudo.toggled.connect(self.Cplus.setEnabled)
        layout.addItem(QtGui.QSpacerItem(5, 5, QtGui.QSizePolicy.Expanding,
                                         QtGui.QSizePolicy.Expanding), 29, 1, 1, 2)

        newComponent.loadUI(self)

        self.Petroleo = Petroleo()
        self.Crudo = Crudo()

    @property
    def unknown(self):
        if self.checkCrudo.isChecked():
            return self.Crudo
        else:
            return self.Petroleo

    def changeUnknown(self):
        self.status.setState(self.unknown.status, self.unknown.msg)
        self.buttonShowDetails.setEnabled(self.unknown.status)
        self.buttonBox.button(QtGui.QDialogButtonBox.Save).setEnabled(
            self.unknown.status)

    def changeCurva(self):
        temp = self.curvaDestilacion.getColumn(0)
        T_dist = self.curvaDestilacion.getColumn(1)
        curvas = ["D86", "TBP", "EFV", "D1160", "SD"]
        x = curvas.pop(self.tipoCurva.currentIndex())
        kwargs = {}
        kwargs[x] = temp
        kwargs["Tdist"] = T_dist
        for curva in curvas:
            kwargs[curva] = []
        self.calculo(**kwargs)


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    Dialog = Definicion_Petro()
    # Dialog = Ui_Contribution("Ambrose")
    Dialog.show()
    sys.exit(app.exec_())
