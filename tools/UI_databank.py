#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2016, Juan José Gómez Romera <jjgomera@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>."""


###############################################################################
# Database tools
###############################################################################

import os

from PyQt5 import QtGui, QtWidgets

from lib import sql
from UI import viewComponents


class UI_databank_widget(QtWidgets.QWidget):
    """Database widget, to use in whatever need: database dialog, proyect
    component list definnition"""
    def __init__(self, parent=None):
        super(UI_databank_widget, self).__init__(parent)
        gridLayout = QtWidgets.QGridLayout(self)

        layoutTitle = QtWidgets.QHBoxLayout()
        layoutTitle.setSpacing(5)
        self.buttonNew = QtWidgets.QToolButton(self)
        self.buttonNew.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Create new element"))
        self.buttonNew.setIcon(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"]+"/images/button/fileNew.png")))
        self.buttonNew.clicked.connect(self.newComponent)
        layoutTitle.addWidget(self.buttonNew)
        self.buttonCopy = QtWidgets.QToolButton(self)
        self.buttonCopy.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Clone this element"))
        self.buttonCopy.setIcon(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"]+"/images/button/editCopy.png")))
        self.buttonCopy.clicked.connect(self.copyComponent)
        layoutTitle.addWidget(self.buttonCopy)
        self.buttonDelete = QtWidgets.QToolButton(self)
        self.buttonDelete.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Delete element"))
        self.buttonDelete.setIcon(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"]+"/images/button/editDelete.png")))
        self.buttonDelete.clicked.connect(self.deleteComponent)
        self.buttonDelete.setEnabled(False)
        layoutTitle.addWidget(self.buttonDelete)
        gridLayout.addItem(layoutTitle, 1, 1)

        gridLayout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Find")+": "), 1, 2)
        self.Busqueda = QtWidgets.QLineEdit()
        self.Busqueda.textChanged.connect(self.buscar)
        gridLayout.addWidget(self.Busqueda, 1, 3)
        self.siguiente = QtWidgets.QPushButton(
            QtWidgets.QApplication.translate("pychemqt", "Next"))
        self.siguiente.clicked.connect(self.Next)
        gridLayout.addWidget(self.siguiente, 1, 4)

        self.BaseDatos = QtWidgets.QTableWidget()
        self.BaseDatos.setMinimumWidth(375)
        self.BaseDatos.verticalHeader().hide()
        self.BaseDatos.setEditTriggers(
            QtWidgets.QAbstractItemView.NoEditTriggers)
        self.BaseDatos.setShowGrid(False)
        self.BaseDatos.setRowCount(0)
        self.BaseDatos.setColumnCount(3)
        self.BaseDatos.setHorizontalHeaderItem(
            0, QtWidgets.QTableWidgetItem("Id"))
        self.BaseDatos.setHorizontalHeaderItem(1, QtWidgets.QTableWidgetItem(
            QtWidgets.QApplication.translate("pychemqt", "Name")))
        self.BaseDatos.setHorizontalHeaderItem(2, QtWidgets.QTableWidgetItem(
            QtWidgets.QApplication.translate("pychemqt", "Formula")))
        self.BaseDatos.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectRows)
        self.BaseDatos.horizontalHeader().setStretchLastSection(True)
        self.BaseDatos.currentCellChanged.connect(self.checkButton)
        self.BaseDatos.doubleClicked.connect(self.mostrarPropiedades)
        gridLayout.addWidget(self.BaseDatos, 2, 1, 1, 4)
        self.correctos = []
        self.indice = 0
        self.rellenar()

    def rellenar(self):
        """Fill in list with component from database"""
        self.BaseDatos.setRowCount(0)
        sql.databank.execute("select * from compuestos")
        for i in sql.databank:
            self.BaseDatos.setRowCount(self.BaseDatos.rowCount()+1)
            self.BaseDatos.setItem(
                i[0]-1, 0, QtWidgets.QTableWidgetItem(str(i[0])))
            self.BaseDatos.setItem(i[0]-1, 1, QtWidgets.QTableWidgetItem(i[2]))
            self.BaseDatos.setItem(i[0]-1, 2, QtWidgets.QTableWidgetItem(i[1]))
            self.BaseDatos.setRowHeight(self.BaseDatos.rowCount()-1, 20)

        sql.databank_Custom.execute("select * from compuestos")
        for i in sql.databank_Custom:
            filas = self.BaseDatos.rowCount()
            self.BaseDatos.setRowCount(filas+1)
            self.BaseDatos.setItem(
                filas, 0, QtWidgets.QTableWidgetItem(str(i[0])))
            self.BaseDatos.setItem(filas, 1, QtWidgets.QTableWidgetItem(i[2]))
            self.BaseDatos.setItem(filas, 2, QtWidgets.QTableWidgetItem(i[1]))
            self.BaseDatos.setRowHeight(self.BaseDatos.rowCount()-1, 20)

        self.BaseDatos.resizeColumnsToContents()

    def buscar(self):
        """Search str at database"""
        self.indice = 0
        texto = "%"+self.Busqueda.text()+"%"
        query = "select * from compuestos where "
        query += "nombre LIKE '%s' or formula LIKE '%s'" % (texto, texto)
        sql.databank.execute(query)
        self.correctos = []
        for i in sql.databank:
            self.correctos.append(i[0])
        self.BaseDatos.setCurrentCell(self.correctos[self.indice]-1, 0)

    def Next(self):
        """Show next coincidence with search string"""
        if self.indice < len(self.correctos)-1:
            self.indice += 1
        else:
            self.indice = 0
        self.BaseDatos.setCurrentCell(self.correctos[self.indice]-1, 0)

    def checkButton(self, indice):
        """Edit action are only available in custom database elements"""
        if indice >= sql.N_comp:
            self.buttonDelete.setEnabled(True)
        else:
            self.buttonDelete.setEnabled(False)

    def mostrarPropiedades(self):
        """Show properties of selected component"""
        indice = self.currentIndice
        if indice > 0:
            Dialog = viewComponents.View_Component(indice)
            Dialog.exec_()

    @property
    def currentIndice(self):
        value = self.BaseDatos.item(self.BaseDatos.currentRow(), 0).text()
        return float(value)

    def currentRow(self):
        return self.BaseDatos.currentRow()

    def newComponent(self):
        Dialog = viewComponents.View_Component(0)
        if Dialog.exec_():
            self.rellenar()

    def copyComponent(self):
        sql.copyElement(self.currentIndice)
        self.rellenar()

    def deleteComponent(self):
        sql.deleteElement(self.currentIndice)
        self.rellenar()


class UI_databank(QtWidgets.QDialog):
    """Database dialog"""
    def __init__(self, parent=None):
        super(UI_databank, self).__init__(parent)
        self.setWindowTitle(QtWidgets.QApplication.translate(
            "pychemqt", "Components database"))
        layout = QtWidgets.QVBoxLayout(self)
        self.databank = UI_databank_widget()
        layout.addWidget(self.databank)
        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Cancel | QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = UI_databank()
    Dialog.show()
    sys.exit(app.exec_())
