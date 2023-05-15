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
# Divider equipment dialog
###############################################################################


from functools import partial

from tools.qt import QtCore, QtWidgets, tr

from lib.unidades import Pressure, MassFlow
from equipment.parents import UI_equip
from equipment.flux import Divider
from UI import UI_corriente
from UI.widgets import Entrada_con_unidades, Tabla
from UI.delegate import CellEditor


class UI_equipment(UI_equip):
    """Divider equipment edition dialog"""
    Equipment = Divider()

    def __init__(self, equipment=None, salidas=0, parent=None):
        """
        equipment: Initial equipment instance to model
        salidas: Stream Output number to equipment
        """
        super().__init__(Divider, entrada=False, parent=parent)

        # Calculate tab
        lyt_Calc = QtWidgets.QGridLayout(self.tabCalculo)
        lyt_Calc.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Separation")), 1, 1, 1, 1)
        self.criterio = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_CRITERIO:
            self.criterio.addItem(txt)
        self.criterio.currentIndexChanged.connect(self.criterio_Changed)
        lyt_Calc.addWidget(self.criterio, 1, 2, 1, 1)

        self.fracciones = Tabla(1, horizontalHeader=[True], stretch=False)
        self.fracciones.setItemDelegateForColumn(0, CellEditor(self))
        lyt_Calc.addWidget(self.fracciones, 2, 1, 1, 2)

        lyt_Calc.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Pressure lost")), 3, 1, 1, 1)
        self.deltaP = Entrada_con_unidades(Pressure, value=0)
        self.deltaP.valueChanged.connect(partial(self.changeParams, "deltaP"))
        lyt_Calc.addWidget(self.deltaP, 3, 2, 1, 1)

        if equipment and salidas:
            equipment(salidas=salidas)
        elif equipment:
            salidas = equipment.kwargs["salidas"]
        else:
            self.Equipment = Divider(salidas=salidas)

        self.fracciones.setRowCount(salidas)
        for i in range(salidas):
            itm = QtWidgets.QTableWidgetItem("%i" % i)
            itm.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter)
            self.fracciones.setItem(0, i, itm)
            self.fracciones.setRowHeight(i, 20)
            widget = UI_corriente.Ui_corriente(readOnly=True)
            self.Salida.addTab(widget, str(i+1))

        self.criterio_Changed(0)
        self.fracciones.editingFinished.connect(
            partial(self.changeParams, "split"))
        self.setEquipment(equipment)

    def criterio_Changed(self, int):
        if int:
            item = QtWidgets.QTableWidgetItem(tr(
                "pychemqt", "Flow")+", "+MassFlow.text())
            self.fracciones.setHorizontalHeaderItem(0, item)
            self.fracciones.item(self.fracciones.rowCount()-1, 0).setFlags(
                QtCore.Qt.ItemFlag.ItemIsEditable | QtCore.Qt.ItemFlag.ItemIsEnabled |
                QtCore.Qt.ItemFlag.ItemIsSelectable)
        else:
            item = QtWidgets.QTableWidgetItem(tr(
                "pychemqt", "Flow")+", "+MassFlow.text())
            self.fracciones.setHorizontalHeaderItem(0, item)
            self.fracciones.item(self.fracciones.rowCount()-1, 0).setFlags(
                QtCore.Qt.ItemFlag.NoItemFlags)
        self.changeParams("criterio", int)

    def changeParams(self, parametro, valor=None):
        if parametro == "split":
            valor = self.fracciones.getColumn(0, False)
            if self.criterio.currentIndex() == 0:
                if len(valor)+1 < self.fracciones.rowCount():
                    return
                elif len(valor)+1 == self.fracciones.rowCount():
                    valor.append(1-sum(valor))
                elif len(valor) == self.fracciones.rowCount():
                    valor[-1] = 1-sum(valor[:-1])
        self.calculo(**{parametro: valor})

    def rellenar(self):
        UI_equip.rellenar(self)
        if self.Equipment.status == 1 and self.criterio.currentIndex() == 1:
                self.entrada.setCorriente(self.Equipment.entrada)

    def rellenarInput(self):
        UI_equip.rellenarInput(self)
        self.fracciones.setColumn(0, self.Equipment.kwargs["split"])


if __name__ == "__main__":
    import sys
    from lib.corriente import Corriente
    app = QtWidgets.QApplication(sys.argv)
    agua = Corriente(T=300, P=101325, caudalMasico=1, ids=[62],
                     fraccionMolar=[0.1])
    divisor = Divider(entrada=agua, salidas=2, split=[0.3, 0.7])
    dialogo = UI_equipment(divisor)
    dialogo.show()
    sys.exit(app.exec())
