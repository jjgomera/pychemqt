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
# Tools with units converter UI
#    - UI_conversorUnidades: Dialog to show all values of a magnitud in every
#        supported unit
#    - moneda: Customized unit converter for currency values
#    - UI_unitConverter: Dialog to choose between the availables unit
###############################################################################


import logging
import os

from PyQt5 import QtCore, QtGui, QtWidgets

from lib.unidades import Currency, getrates, _all
from lib.config import conf_dir
from UI.delegate import CellEditor


class UI_conversorUnidades(QtWidgets.QDialog):
    """Dialog to show all values of a magnitud in every supported unit"""
    def __init__(self, unidad, valor=None, parent=None):
        """
        unidad: unidades class
        valor: optional value to show
        """
        super(UI_conversorUnidades, self).__init__(parent)
        self.setWindowTitle(unidad.__title__)
        self.mutex = QtCore.QMutex()

        self.unidad = unidad
        if unidad.__tooltip__:
            self.tooltip = unidad.__tooltip__
        else:
            self.tooltip = unidad.__text__
        self.value = self.unidad(valor)

        lyt = QtWidgets.QGridLayout(self)
        self.tabla = QtWidgets.QTableWidget()
        self.tabla.setRowCount(len(unidad.__text__))
        self.tabla.setColumnCount(1)
        self.tabla.setItemDelegateForColumn(0, CellEditor(self))
        self.tabla.horizontalHeader().setVisible(False)
        self.tabla.horizontalHeader().setStretchLastSection(True)
        self.tabla.setMaximumHeight(len(unidad.__text__)*24+4)

        for i, txt in enumerate(unidad.__text__):
            item = QtWidgets.QTableWidgetItem(txt)
            self.tabla.setVerticalHeaderItem(i, item)
            self.tabla.setRowHeight(i, 24)
            self.tabla.setItem(i, 0, QtWidgets.QTableWidgetItem(""))
            self.tabla.item(i, 0).setTextAlignment(
                QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        for i, tip in enumerate(self.tooltip):
            self.tabla.item(i, 0).setToolTip(tip)
        self.tabla.cellChanged.connect(self.update)
        lyt.addWidget(self.tabla, 2, 1)

        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Cancel | QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setCenterButtons(True)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        lyt.addWidget(self.buttonBox, 3, 0, 1, 3)

        if valor:
            self.fill(self.value)

    def fill(self, valor):
        for i, key in enumerate(self.unidad.__units__):
            self.tabla.item(i, 0).setText(valor.format(key))

    def update(self, fila, columna):
        if self.mutex.tryLock():
            new = self.tabla.item(fila, columna).text()
            self.value = self.unidad(float(new), self.unidad.__units__[fila])
            self.fill(self.value)
            self.mutex.unlock()


class moneda(UI_conversorUnidades):
    """Customized unit converter for currency values
        - Add date and exchange update button at top
        - Add country flags to easy recognize
    """
    def __init__(self, valor=None, parent=None):
        super(moneda, self).__init__(Currency, valor=valor, parent=parent)

        self.fecha = QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Date") + ": " + self.value.date)
        self.layout().addWidget(self.fecha, 0, 1)
        self.botonActualizar = QtWidgets.QPushButton(
            QtWidgets.QApplication.translate("pychemqt", "Update"))
        self.botonActualizar.clicked.connect(self.getrates)
        self.layout().addWidget(self.botonActualizar, 1, 1)

        for i in range(len(Currency.__units__)):
            header = self.tabla.verticalHeaderItem(i)
            header.setIcon(QtGui.QIcon(QtGui.QPixmap(
                os.environ["pychemqt"] +
                "/images/flag/%s.png" % Currency.__units__[i])))
            # Set backgroundcolor to better look or rare currencies
            # Use olimpic continent color code
            main = len(Currency._uMain)
            Europe = main + len(Currency._uEurope)
            America = Europe + len(Currency._uAmerica)
            Africa = America + len(Currency._uAfrica)
            Asia = Africa + len(Currency._uAsia)
            Oceania = Asia + len(Currency._uOceania)

            if i < main:
                color = "#FFFFFF"
            elif i < Europe:
                color = "#FF5555"
            elif i < America:
                color = "#55FF55"
            elif i < Africa:
                color = "#888888"
            elif i < Asia:
                color = "#FFFF00"
            elif i < Oceania:
                color = "#5555FF"
            else:
                color = "#FF8800"
            header.setBackground(QtGui.QBrush(QtGui.QColor(color)))

    def getrates(self):
        filename = conf_dir + "moneda.dat"
        getrates(filename)
        self.value = self.unidad(self.value)
        self.fecha.setText(QtWidgets.QApplication.translate(
            "pychemqt", "Date") + ": " + self.value.date)
        if self.value != 0:
            self.update(0, 0)


class UI_unitConverter(QtWidgets.QDialog):
    """Dialog to choose between the availables unit,
    used as main window tool
    """

    def __init__(self, parent=None):
        super(UI_unitConverter, self).__init__(parent)
        self.setWindowTitle(
            QtWidgets.QApplication.translate("pychemqt", "Units converter"))

        self.verticalLayout = QtWidgets.QVBoxLayout(self)
        self.lista = QtWidgets.QListWidget()
        self.lista.itemDoubleClicked.connect(self.showChildWindow)
        self.verticalLayout.addWidget(self.lista)
        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Close)
        self.buttonBox.rejected.connect(self.reject)
        self.verticalLayout.addWidget(self.buttonBox)
        for unidad in _all:
            self.lista.addItem(unidad.__title__)

        self.lista.setCurrentRow(-1)
        logging.info(QtWidgets.QApplication.translate(
            "pychemqt", "Starting unit converte tool"))

    def showChildWindow(self):
        """Show child window with selected unit converter"""
        indice = self.lista.currentRow()
        if _all[indice].__name__ == "Currency":
            dialog = moneda()
        else:
            dialog = UI_conversorUnidades(_all[indice])
        dialog.exec_()


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    # conversion_unidades = UI_unitConverter()
    conversion_unidades = moneda()
    conversion_unidades.show()
    sys.exit(app.exec_())
