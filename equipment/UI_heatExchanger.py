#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>."""


###############################################################################
# Generic heat exchanger equipment dialog
###############################################################################


from functools import partial

from PyQt5 import QtWidgets

from lib.unidades import (Temperature, DeltaT, DeltaP, Power, Area,
                          HeatTransfCoef)
from UI.widgets import Entrada_con_unidades
from equipment.heatExchanger import Heat_Exchanger
from equipment.parents import UI_equip


class UI_equipment(UI_equip):
    """Generic heat exchanger equipment edition dialog"""
    Equipment = Heat_Exchanger()

    def __init__(self, equipment=None, parent=None):
        """
        equipment: Initial equipment instance to model
        """
        super().__init__(Heat_Exchanger, entrada=False, salida=False,
                         parent=parent)

        # Calculate tab
        lyt = QtWidgets.QGridLayout(self.tabCalculo)
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Output temperature")), 1, 1)
        self.Tout = Entrada_con_unidades(Temperature)
        self.Tout.valueChanged.connect(partial(self.changeParams, "Tout"))
        lyt.addWidget(self.Tout, 1, 2)
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Temperature increase")), 2, 1)
        self.DeltaT = Entrada_con_unidades(DeltaT)
        self.DeltaT.valueChanged.connect(partial(self.changeParams, "DeltaT"))
        lyt.addWidget(self.DeltaT, 2, 2)
        lyt.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Heat Duty")), 3, 1)
        self.Heat = Entrada_con_unidades(Power)
        self.Heat.valueChanged.connect(partial(self.changeParams, "Heat"))
        lyt.addWidget(self.Heat, 3, 2)
        group = QtWidgets.QGroupBox(
            QtWidgets.QApplication.translate("pychemqt", "Heat Transfer"))
        lyt.addWidget(group, 4, 1, 1, 2)
        lyt1 = QtWidgets.QGridLayout(group)
        lyt1.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Area")), 1, 1)
        self.A = Entrada_con_unidades(Area)
        self.A.valueChanged.connect(partial(self.changeParams, "A"))
        lyt1.addWidget(self.A, 1, 2)
        lyt1.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Heat Transfer Coefficient")), 2, 1)
        self.U = Entrada_con_unidades(HeatTransfCoef)
        self.U.valueChanged.connect(partial(self.changeParams, "U"))
        lyt1.addWidget(self.U, 2, 2)
        lyt1.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "External Temperature")), 3, 1)
        self.Text = Entrada_con_unidades(Temperature)
        self.Text.valueChanged.connect(partial(self.changeParams, "Text"))
        lyt1.addWidget(self.Text, 3, 2)

        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            5, 0, 1, 3)
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Pressure loss")), 6, 1)
        self.DeltaP = Entrada_con_unidades(DeltaP, value=0)
        self.DeltaP.valueChanged.connect(partial(self.changeParams, "DeltaP"))
        lyt.addWidget(self.DeltaP, 6, 2)
        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 7, 0, 1, 3)

        group = QtWidgets.QGroupBox(
            QtWidgets.QApplication.translate("pychemqt", "Results"))
        lyt.addWidget(group, 8, 1, 1, 5)
        layout = QtWidgets.QGridLayout(group)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Heat Duty")), 0, 1)
        self.HeatCalc = Entrada_con_unidades(Power, retornar=False)
        self.HeatCalc.setReadOnly(True)
        layout.addWidget(self.HeatCalc, 0, 2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Output Temperature")), 1, 1)
        self.ToutCalc = Entrada_con_unidades(Temperature, retornar=False)
        self.ToutCalc.setReadOnly(True)
        layout.addWidget(self.ToutCalc, 1, 2)

        lyt.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            9, 0, 1, 3)

        if equipment:
            self.setEquipment(equipment)

    def changeParams(self, parametro, valor):
        if parametro == "Tout":
            self.Heat.clear()
            self.DeltaT.clear()
        elif parametro == "DeltaT":
            self.Heat.clear()
            self.Tout.clear()
        elif parametro == "Heat":
            self.DeltaT.clear()
            self.Tout.clear()
        self.calculo(**{parametro: valor})


if __name__ == "__main__":
    import sys
    from lib.corriente import Corriente
    app = QtWidgets.QApplication(sys.argv)
    agua = Corriente(T=300, P=101325, caudalMasico=1, ids=[62],
                     fraccionMolar=[1.])
    cambiador = Heat_Exchanger(entrada=agua, Tout=90+273.15)
    dialogo = UI_equipment(cambiador)
    dialogo.show()
    sys.exit(app.exec_())
