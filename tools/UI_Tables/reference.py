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
#   Dialogs for configuration:
#   - Ui_ReferenceState: Dialog to select reference state
#   - Ui_Properties: Dialog for select and sort shown properties in tables
###############################################################################


import os

from PyQt5 import QtCore, QtGui, QtWidgets

from lib import unidades
from lib.thermo import ThermoAdvanced
from UI.delegate import CheckEditor
from UI.widgets import Entrada_con_unidades

from tools.UI_Tables.library import N_PROP


class Ui_ReferenceState(QtWidgets.QDialog):
    """Dialog for select reference state"""
    def __init__(self, config=None, parent=None):
        """config: instance with project config to set initial values"""
        super(Ui_ReferenceState, self).__init__(parent)
        self.setWindowTitle(QtWidgets.QApplication.translate(
            "pychemqt", "Select reference state"))
        layout = QtWidgets.QGridLayout(self)
        self.OTO = QtWidgets.QRadioButton(QtWidgets.QApplication.translate(
            "pychemqt", "OTO,  h,s=0 at 298K and 1 atm"))
        layout.addWidget(self.OTO, 0, 1, 1, 7)
        self.NBP = QtWidgets.QRadioButton(QtWidgets.QApplication.translate(
            "pychemqt", "NBP,  h,s=0 saturated liquid at Tb"))
        layout.addWidget(self.NBP, 1, 1, 1, 7)
        self.IIR = QtWidgets.QRadioButton(QtWidgets.QApplication.translate(
            "pychemqt", "IIR,  h=200,s=1 saturated liquid at 273K"))
        layout.addWidget(self.IIR, 2, 1, 1, 7)
        self.ASHRAE = QtWidgets.QRadioButton(QtWidgets.QApplication.translate(
            "pychemqt", "ASHRAE,  h,s=0 saturated liquid at 243K"))
        layout.addWidget(self.ASHRAE, 3, 1, 1, 7)
        self.personalizado = QtWidgets.QRadioButton(
            QtWidgets.QApplication.translate("pychemqt", "Custom"))
        self.personalizado.toggled.connect(self.setEnabled)
        layout.addWidget(self.personalizado, 4, 1, 1, 7)

        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed, QtWidgets.QSizePolicy.Policy.Fixed),
            5, 1)
        layout.addWidget(QtWidgets.QLabel("T:"), 5, 2)
        self.T = Entrada_con_unidades(unidades.Temperature, value=298.15)
        layout.addWidget(self.T, 5, 3)
        layout.addWidget(QtWidgets.QLabel("P:"), 6, 2)
        self.P = Entrada_con_unidades(unidades.Pressure, value=101325)
        layout.addWidget(self.P, 6, 3)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed), 5, 4, 2, 1)
        layout.addWidget(QtWidgets.QLabel("h:"), 5, 5)
        self.h = Entrada_con_unidades(unidades.Enthalpy, value=0)
        layout.addWidget(self.h, 5, 6)
        layout.addWidget(QtWidgets.QLabel("s:"), 6, 5)
        self.s = Entrada_con_unidades(unidades.SpecificHeat, value=0)
        layout.addWidget(self.s, 6, 6)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 7, 7)

        buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)
        layout.addWidget(buttonBox, 8, 1, 1, 7)

        if config and config.has_option("MEoS", "reference"):
            self.setEnabled(False)
            if config.get("MEoS", "reference") == "OTO":
                self.OTO.setChecked(True)
            elif config.get("MEoS", "reference") == "NBP":
                self.NBP.setChecked(True)
            elif config.get("MEoS", "reference") == "IIR":
                self.IIR.setChecked(True)
            elif config.get("MEoS", "reference") == "ASHRAE":
                self.ASHRAE.setChecked(True)
            else:
                self.personalizado.setChecked(True)
                self.setEnabled(True)
                self.T.setValue(config.getfloat("MEoS", "T"))
                self.P.setValue(config.getfloat("MEoS", "P"))
                self.h.setValue(config.getfloat("MEoS", "h"))
                self.s.setValue(config.getfloat("MEoS", "s"))
        else:
            self.OTO.setChecked(True)
            self.setEnabled(False)

    def setEnabled(self, bool):
        """Enable custom entriees"""
        self.T.setEnabled(bool)
        self.P.setEnabled(bool)
        self.h.setEnabled(bool)
        self.s.setEnabled(bool)


class Ui_Properties(QtWidgets.QDialog):
    """Dialog for select and sort shown properties in tables"""
    _default = [1, 0, 1, 0, 0, 1, 0, 1, 1]+[0]*(N_PROP-9)

    def __init__(self, config=None, parent=None):
        super(Ui_Properties, self).__init__(parent)
        if config and config.has_option("MEoS", "properties"):
            values = config.get("MEoS", "properties")
            if isinstance(values, str):
                values = eval(values)
            fase = config.getboolean("MEoS", "phase")
            self.order = config.get("MEoS", "propertiesOrder")
            if isinstance(self.order, str):
                self.order = eval(self.order)
        else:
            values = self._default
            fase = False
            self.order = list(range(N_PROP))

        self.setWindowTitle(
            QtWidgets.QApplication.translate("pychemqt", "Select Properties"))
        layout = QtWidgets.QGridLayout(self)
        self.prop = QtWidgets.QTableWidget(len(ThermoAdvanced.properties()), 2)
        self.prop.verticalHeader().hide()
        self.prop.horizontalHeader().hide()
        self.prop.horizontalHeader().setStretchLastSection(True)
        self.prop.setGridStyle(QtCore.Qt.PenStyle.NoPen)
        self.prop.setColumnWidth(0, 18)
        self.prop.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.prop.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.prop.setItemDelegateForColumn(0, CheckEditor(self))
        for i, value in enumerate(values):
            if value == 1:
                val = "1"
            else:
                val = ""
            self.prop.setItem(i, 0, QtWidgets.QTableWidgetItem(val))
            name = ThermoAdvanced.propertiesName()[self.order[i]]
            self.prop.setItem(i, 1, QtWidgets.QTableWidgetItem(name))
            self.prop.setRowHeight(i, 20)
            self.prop.openPersistentEditor(self.prop.item(i, 0))
        self.prop.currentCellChanged.connect(self.comprobarBotones)
        self.prop.cellDoubleClicked.connect(self.toggleCheck)
        layout.addWidget(self.prop, 1, 1, 6, 1)

        self.ButtonArriba = QtWidgets.QToolButton()
        self.ButtonArriba.setIcon(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"] +
            os.path.join("images", "button", "arrow-up.png"))))
        self.ButtonArriba.clicked.connect(self.Up)
        layout.addWidget(self.ButtonArriba, 3, 2, 1, 1)
        self.ButtonAbajo = QtWidgets.QToolButton()
        self.ButtonAbajo.setIcon(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"] +
            os.path.join("images", "button", "arrow-down.png"))))
        self.ButtonAbajo.clicked.connect(self.Down)
        layout.addWidget(self.ButtonAbajo, 4, 2, 1, 1)

        self.checkFase = QtWidgets.QCheckBox(QtWidgets.QApplication.translate(
            "pychemqt", "Show bulk, liquid and vapor properties"))
        self.checkFase.setChecked(fase)
        layout.addWidget(self.checkFase, 7, 1, 1, 2)
        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Reset | QtWidgets.QDialogButtonBox.StandardButton.Ok |
            QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        self.buttonBox.addButton(
            QtWidgets.QApplication.translate("pychemqt", "Mark all"),
            QtWidgets.QDialogButtonBox.ButtonRole.ResetRole)
        self.buttonBox.addButton(
            QtWidgets.QApplication.translate("pychemqt", "No Mark"),
            QtWidgets.QDialogButtonBox.ButtonRole.ResetRole)
        self.btYes = QtWidgets.QPushButton
        self.buttonBox.clicked.connect(self.buttonClicked)
        layout.addWidget(self.buttonBox, 8, 1, 1, 2)

    def toggleCheck(self, fila, columna):
        """Toggle check status with a doubleclick in row"""
        txt = self.prop.item(fila, 0).text()
        if txt == "0":
            newtxt = "1"
        else:
            newtxt = ""
        self.prop.item(fila, 0).setText(newtxt)

    def Down(self):
        """Change current selected row with next row"""
        i = self.prop.currentRow()
        txt = self.prop.item(i, 0).text()
        self.prop.item(i, 0).setText(self.prop.item(i+1, 0).text())
        self.prop.item(i+1, 0).setText(txt)
        item = self.prop.takeItem(i, 1)
        self.prop.setItem(i, 1, self.prop.takeItem(i+1, 1))
        self.prop.setItem(i+1, 1, item)
        self.prop.setCurrentCell(i+1, 0)
        self.order[i], self.order[i+1] = self.order[i+1], self.order[i]

    def Up(self):
        """Change current selected row with previous row"""
        i = self.prop.currentRow()
        txt = self.prop.item(i, 0).text()
        self.prop.item(i, 0).setText(self.prop.item(i-1, 0).text())
        self.prop.item(i-1, 0).setText(txt)
        item = self.prop.takeItem(i, 1)
        self.prop.setItem(i, 1, self.prop.takeItem(i-1, 1))
        self.prop.setItem(i-1, 1, item)
        self.prop.setCurrentCell(i-1, 0)
        self.order[i], self.order[i-1] = self.order[i-1], self.order[i]

    def buttonClicked(self, boton):
        """Actions for dialogbuttonbox functionality"""
        if self.buttonBox.buttonRole(boton) == \
                QtWidgets.QDialogButtonBox.ButtonRole.AcceptRole:
            self.accept()
        elif self.buttonBox.buttonRole(boton) == \
                QtWidgets.QDialogButtonBox.ButtonRole.RejectRole:
            self.reject()
        else:
            if boton == \
                    self.buttonBox.button(QtWidgets.QDialogButtonBox.StandardButton.Reset):
                values = self._default
            elif boton.text() == \
                    QtWidgets.QApplication.translate("pychemqt", "No Mark"):
                values = [0]*N_PROP
            else:
                values = [1]*N_PROP

            for i, propiedad in enumerate(values):
                if propiedad == 1:
                    val = "1"
                else:
                    val = ""
                self.prop.item(i, 0).setText(val)

    def properties(self):
        """Properties list"""
        value = []
        for i in range(self.prop.rowCount()):
            value.append(self.prop.cellWidget(i, 0).isChecked())
        return value

    def comprobarBotones(self, fila):
        """Check if button are enabled or disabled"""
        self.ButtonArriba.setEnabled(fila >= 1)
        self.ButtonAbajo.setEnabled(fila < self.prop.rowCount()-1)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)

    SteamTables = Ui_Properties()
    # SteamTables = Ui_ReferenceState()

    SteamTables.show()
    sys.exit(app.exec())
