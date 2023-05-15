#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>."""


###############################################################################
# Spreadsheet interaction equipment dialog
###############################################################################


import os

from tools.qt import QtCore, QtGui, QtWidgets, tr

try:
    import ezodf
except:
    pass

from UI.widgets import PathConfig, Tabla
from equipment.parents import UI_equip
from equipment.spreadsheet import Spreadsheet


class TableDelegate(QtWidgets.QItemDelegate):
    """Delegate table with combobox options"""
    def __init__(self, owner, items=None):
        super(TableDelegate, self).__init__(owner)
        if not items:
            items = {}
            for ind in range(4):
                items[ind] = []
        self.setItems(items)

    def setItems(self, items):
        self.items = items

    def setItemsByIndex(self, index, items):
        self.items[index] = items

    def createEditor(self, parent, option, index):
        if index.column() < 4:
            self.editor = QtWidgets.QComboBox(parent)
            self.editor.addItems(self.items[index.column()])
        else:
            self.editor = QtWidgets.QLineEdit(parent)
            regExp = QtCore.QRegExp("[A-Z]{1,3}\\d{1,5}")
            validator = QtGui.QRegExpValidator(regExp)
            self.editor.setValidator(validator)
        return self.editor

    def setEditorData(self, editor, index):
        value = index.data(QtCore.Qt.ItemDataRole.DisplayRole)
        if index.column() < 4:
            try:
                num = self.items[index.column()].index(value)
            except ValueError:
                num = -1
            editor.setCurrentIndex(num)
        else:
            editor.setText(value)

    def setModelData(self, editor, model, index):
        if index.column() < 4:
            value = editor.currentText()
        else:
            value = editor.text().upper()

        model.setData(index, QtCore.QVariant(value), QtCore.Qt.ItemDataRole.DisplayRole)


class UI_equipment(UI_equip):
    """Spreadsheet interaction equipment edition dialog"""
    Equipment = Spreadsheet()

    def __init__(self, equipment=None, project=None, parent=None):
        """
        equipment: Initial equipment instance to model
        """
        super().__init__(Spreadsheet, entrada=True, salida=True,
                         calculo=False, parent=parent)
        self.project = project

        # Calculate tab
        layout = QtWidgets.QGridLayout(self.Entrada)
        label = tr(
            "pychemqt", "Spreadsheet path")+":"
        msg = tr(
            "pychemqt", "Select Spreadsheet")
        patrones = []
        if os.environ["ezodf"]:
            patrones.append(tr(
                "pychemqt", "Libreoffice spreadsheet files")+" (*.ods)")
        if os.environ["xlwt"]:
            patrones.append(tr(
                "pychemqt", "Microsoft Excel 97/2000/XP/2003 XML")+" (*.xls)")
        if os.environ["openpyxl"]:
            patrones.append(tr(
                "pychemqt", "Microsoft Excel 2007/2010 XML")+" (*.xlsx)")
        patron = ";;".join(patrones)
        self.filename = PathConfig(label, msg=msg, patron=patron)
        self.filename.valueChanged.connect(self.changeSpreadsheet)
        layout.addWidget(self.filename, 1, 1)
        header = [tr("pychemqt", "Entity"),
                  tr("pychemqt", "Variable"),
                  tr("pychemqt", "Unit value"),
                  tr("pychemqt", "Sheet"),
                  tr("pychemqt", "Cell")]
        self.datamap = Tabla(
            5, filas=1, dinamica=True, horizontalHeader=header,
            verticalHeader=False, orientacion=QtCore.Qt.AlignmentFlag.AlignLeft,
            delegate=None, delegateforRow=TableDelegate, parent=self)
        self.datamap.setEnabled(False)
        self.datamap.cellChanged.connect(self.cellChanged)
        self.datamap.rowFinished.connect(self.addRow)
        layout.addWidget(self.datamap, 2, 1)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 10, 1)

        entitys = []
        for stream in list(self.project.streams.keys()):
            entitys.append("s%i" % stream)
        for equip in list(self.project.items.keys()):
            if equip[0] == "e":
                entitys.append(equip)
        self.datamap.itemDelegateForRow(0).setItemsByIndex(0, entitys)
        self.entitys = entitys
        if equipment:
            self.setEquipment(equipment)

    def changeSpreadsheet(self, path):
        self.datamap.setEnabled(bool(path))
        self.changeParams("filename", str(path))
        self.datamap.blockSignals(True)
        self.datamap.clear()
        self.datamap.blockSignals(False)
        spreadsheet = ezodf.opendoc(path)
        sheets = [name for name in spreadsheet.sheets.names()]
        self.datamap.itemDelegateForRow(0).setItemsByIndex(3, sheets)

    def rellenarInput(self):
        self.blockSignals(True)
        self.datamap.itemDelegateForRow(
            self.datamap.rowCount()-1).setItemsByIndex(0, self.entitys)
        if self.Equipment.status:
            self.datamap.setEnabled(True)
            self.filename.setText(self.Equipment.kwargs["filename"])
            self.datamap.itemDelegateForRow(0).setItemsByIndex(
                3, self.Equipment.sheets)

        self.datamap.blockSignals(True)
        self.datamap.clear()
        if self.Equipment.kwargs["datamap"]:
            for i, data in enumerate(self.Equipment.kwargs["datamap"]):
                self.datamap.addRow()
                self.datamap.itemDelegateForRow(i).setItemsByIndex(
                    0, self.entitys)
                self.datamap.itemDelegateForRow(i).setItemsByIndex(
                    3, self.Equipment.sheets)
                self.datamap.setItem(
                    i, 0, QtWidgets.QTableWidgetItem(data["entity"]))
                self.datamap.setItem(
                    i, 1, QtWidgets.QTableWidgetItem(data["property"]))
                self.datamap.setItem(
                    i, 2, QtWidgets.QTableWidgetItem(data["unit"]))
                self.datamap.setItem(
                    i, 3, QtWidgets.QTableWidgetItem(data["sheet"]))
                self.datamap.setItem(
                    i, 4, QtWidgets.QTableWidgetItem(data["cell"]))
            self.datamap.itemDelegateForRow(
                self.datamap.rowCount()-1).setItemsByIndex(0, self.entitys)
            self.datamap.itemDelegateForRow(
                self.datamap.rowCount()-1).setItemsByIndex(
                    3, self.Equipment.sheets)
        self.datamap.blockSignals(False)
        self.blockSignals(False)

    def rellenar(self):
        self.rellenarInput()
        self.status.setState(self.Equipment.status, self.Equipment.msg)

    def cellChanged(self, i, j):
        obj = self.project.getObject(str(self.datamap.item(i, 0).text()))
        properties = [prop[0] for prop in obj.propertiesNames()]
        if j == 0:  # Entity cambiado, cambiar variables disponibles
            self.datamap.itemDelegateForRow(i).setItemsByIndex(1, properties)
            editor = QtWidgets.QComboBox()
            editor.addItems(self.datamap.itemDelegateForRow(i).items[1])
            self.datamap.setColumnWidth(1, editor.sizeHint().width())
        elif j == 1:   # Variable cambiada, cambiar unidades disponibles
            value = self.datamap.item(i, 1).text()
            ind = properties.index(value)
            if obj.propertiesUnit()[ind] == str:
                self.datamap.itemDelegateForRow(i).setItemsByIndex(2, [" "])
                self.datamap.item(i, 2).setText(" ")
            else:
                self.datamap.itemDelegateForRow(i).setItemsByIndex(
                    2, obj.propertiesNames()[ind][2].__text__)
        elif j == 3:
            self.datamap.item(i, 4).setText("")

    def addRow(self, fila):
        datamap = self.Equipment.kwargs["datamap"][:]
        data = {}
        data["entity"] = str(fila[0])
        data["property"] = str(fila[1])
        data["unit"] = str(fila[2])
        data["sheet"] = str(fila[3])
        data["cell"] = str(fila[4])
        datamap.append(data)
        self.changeParams("datamap", datamap)


if __name__ == "__main__":
    import sys
    from lib.corriente import Corriente
    from lib.project import Project
    from equipment.heatExchanger import Hairpin
    project = Project()
    project.addItem("i1", Corriente())
    project.addItem("i2", Corriente())
    Cambiador = Hairpin()
    project.addItem("e1", Cambiador)
    project.addStream(1, "i1", "e1", ind_down=0)
    project.addStream(2, "i2", "e1", ind_down=1)
    project.addItem("o1", Corriente())
    project.addStream(3, "e1", "o1", ind_up=0)
    project.addItem("o2", Corriente())
    project.addStream(4, "e1", "o2", ind_up=1)
    caliente = Corriente(T=140+273.15, P=361540., caudalMasico=1.36,
                         ids=[62], fraccionMolar=[1.])
    project.setInput(1, caliente)

    spreadsheet = Spreadsheet(filename="/home/jjgomera/Programacion/pychemqt/Samples/ejemplo.ods",
                              project=project)
    app = QtWidgets.QApplication(sys.argv)
    dialogo = UI_equipment(spreadsheet, project=project)
    dialogo.show()
    sys.exit(app.exec())
