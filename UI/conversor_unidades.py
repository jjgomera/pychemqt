#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
from PyQt4 import QtCore, QtGui

from lib.unidades import Currency, getdata
from UI.delegate import CellEditor


class UI_conversorUnidades(QtGui.QDialog):
    def __init__(self, unidad, valor=None, parent=None):
        super(UI_conversorUnidades, self).__init__(parent)
        self.unidad=unidad
        self.magnitud=unidad.__name__
        self.texto=unidad.__text__
        self.unit=unidad.__units__
        if unidad.__tooltip__:
            self.tooltip=unidad.__tooltip__
        else:
            self.tooltip=unidad.__text__
        self.value=self.unidad(valor)
        self.setWindowTitle(unidad.__title__)
        self.gridLayout = QtGui.QGridLayout(self)
        self.tabla=QtGui.QTableWidget()
        self.tabla.setRowCount(len(self.texto))
        self.tabla.setColumnCount(1)
        self.tabla.setItemDelegateForColumn(0, CellEditor(self))
        self.tabla.horizontalHeader().setVisible(False)
        self.tabla.horizontalHeader().setStretchLastSection(True)
        if self.magnitud in ["SpecificVolume", "Density", "MassFlow", "VolFlow", "ThermalConductivity", "HeatTransfCoef"]:
            self.resize(215, self.minimumHeight())
        elif self.magnitud=="Currency":
            self.resize(250, 500)
        else:
            self.resize(self.minimumSize())

        if self.magnitud in ["Temperature", "Area", "Volume",  "Length", "Angle", "Time"]:
            x=15
        elif self.magnitud in["ThermalConductivity"]:
            x=10
        elif self.magnitud in ["Speed", "Mass", "Acceleration", "Energy", "Enthalpy", "MassFlow", "Diffusivity", "Tension", "Solubility_parameter", "HeatTransfCoef"]:
            x=5
        else: x=0
        self.gridLayout.addItem(QtGui.QSpacerItem(x,15,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),2,0)
        self.gridLayout.addItem(QtGui.QSpacerItem(x,15,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),2,2)

        for i in range(len(self.texto)):
            self.tabla.setVerticalHeaderItem(i, QtGui.QTableWidgetItem(self.texto[i]))
            self.tabla.setRowHeight(i,24)
            self.tabla.setItem(i, 0, QtGui.QTableWidgetItem(""))
            self.tabla.item(i, 0).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)

        for i in range(len(self.tooltip)):
            self.tabla.item(i, 0).setToolTip(QtGui.QApplication.translate("pychemqt", self.tooltip[i]))

        if valor:
            self.rellenarTabla(self.value)
            self.tabla.resizeColumnsToContents()
        if self.magnitud!="Currency":
            self.tabla.setFixedHeight(len(self.texto)*24+4)
        self.gridLayout.addWidget(self.tabla, 2, 1, 1, 1)
        self.tabla.cellChanged.connect(self.actualizar)

        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setCenterButtons(True)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        self.gridLayout.addWidget(self.buttonBox,3,0,1,3)

    def rellenarTabla(self, valor):
        for i, key in enumerate(self.unit):
            self.tabla.item(i, 0).setText(valor.format(key))

    def actualizar(self, fila, columna):
        self.tabla.blockSignals(True)
        self.value=self.unidad(float(self.tabla.item(fila, columna).text()), self.unit[fila])
        self.rellenarTabla(self.value)
        self.tabla.blockSignals(False)


class moneda(UI_conversorUnidades):
    def __init__(self, valor=None, parent=None):
        super(moneda, self).__init__(Currency, valor=valor, parent=parent)

        self.fecha=QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Date::")+self.value.fecha)
        self.gridLayout.addWidget(self.fecha, 0, 1)
        self.botonActualizar=QtGui.QPushButton(QtGui.QApplication.translate("pychemqt", "Update"))
        self.botonActualizar.clicked.connect(self.getdata)
        self.gridLayout.addWidget(self.botonActualizar, 1, 1)

        for i in range(len(Currency.__units__)):
            self.tabla.verticalHeaderItem(i).setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/flag/%s.gif" % Currency.__units__[i])))

    def getdata(self):
        getdata()
        self.value=self.unidad(self.value)
        self.fecha.setText(QtGui.QApplication.translate("pychemqt", "Date:")+self.value.fecha)
        if self.value!=0:
            self.actualizar(0, 0)


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    dialogo = moneda(300)
    dialogo.show()
    sys.exit(app.exec_())

