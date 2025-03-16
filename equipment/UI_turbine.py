#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2025, Juan José Gómez Romera <jjgomera@gmail.com>

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
# Turbine equipment dialog
###############################################################################


from functools import partial

from equipment.compressor import Turbine
from equipment.parents import UI_equip
from lib.unidades import Pressure, Power, Currency
from tools.qt import QtWidgets
from tools.costIndex import CostData
from UI.widgets import Entrada_con_unidades


class UI_equipment(UI_equip):
    """Turbine equipment edition dialog"""
    Equipment = Turbine()

    def __init__(self, equipment=None, parent=None):
        """
        equipment: Initial equipment instance to model
        """
        super().__init__(Turbine, entrada=False, salida=False, parent=parent)

        # Calculate tab
        lyt_Calc = QtWidgets.QGridLayout(self.tabCalculo)
        lyt_Calc.addWidget(QtWidgets.QLabel(self.tr("Method:")), 1, 1)
        self.metodo = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_METODO:
            self.metodo.addItem(txt)
        # self.metodo.addItem(self.tr(
            # "pychemqt", "Especificar curva de funcionamiento"))
        self.metodo.currentIndexChanged.connect(
            self.on_tipoCalculo_currentIndexChanged)
        lyt_Calc.addWidget(self.metodo, 1, 2, 1, 2)
        lyt_Calc.addItem(QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 2, 0, 1, 4)
        lyt_Calc.addWidget(QtWidgets.QLabel(self.tr("Thermodynamic:")), 3, 1)
        self.termodinamica = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_TERMODINAMICA:
            self.termodinamica.addItem(txt)
        self.termodinamica.currentIndexChanged.connect(
            partial(self.changeParams, "termodinamica"))
        lyt_Calc.addWidget(self.termodinamica, 3, 2, 1, 2)
        lyt_Calc.addItem(QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 4, 0, 1, 4)
        lyt_Calc.addWidget(QtWidgets.QLabel(self.tr("Out Pressure")), 5, 1)
        self.Pout = Entrada_con_unidades(Pressure)
        self.Pout.valueChanged.connect(partial(self.changeParams, "Pout"))
        lyt_Calc.addWidget(self.Pout, 5, 2)
        lyt_Calc.addWidget(QtWidgets.QLabel(self.tr("Pressure ratio")), 6, 1)
        self.razon = Entrada_con_unidades(float)
        self.razon.valueChanged.connect(partial(self.changeParams, "razon"))
        lyt_Calc.addWidget(self.razon, 6, 2)
        lyt_Calc.addWidget(QtWidgets.QLabel(self.tr("Efficiency")), 7, 1)
        self.rendimiento = Entrada_con_unidades(float)
        self.rendimiento.valueChanged.connect(
            partial(self.changeParams, "rendimiento"))
        lyt_Calc.addWidget(self.rendimiento, 7, 2)
        lyt_Calc.addWidget(QtWidgets.QLabel(self.tr("Actual Power")), 8, 1)
        self.trabajo = Entrada_con_unidades(Power)
        self.trabajo.valueChanged.connect(
            partial(self.changeParams, "trabajo"))
        lyt_Calc.addWidget(self.trabajo, 8, 2)
        lyt_Calc.setRowStretch(10, 1)

        group = QtWidgets.QGroupBox()
        group.setTitle(self.tr("Results"))
        lyt_Calc.addWidget(group, 12, 1, 1, 2)
        lyt = QtWidgets.QGridLayout(group)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Power")), 1, 1)
        self.power = Entrada_con_unidades(Power, retornar=False, readOnly=True)
        lyt.addWidget(self.power, 1, 2)
        lyt.addWidget(QtWidgets.QLabel("Cp/Cv"), 2, 1)
        self.cp_cv = Entrada_con_unidades(float, retornar=False, readOnly=True)
        lyt.addWidget(self.cp_cv, 2, 2)
        lyt.setColumnStretch(3, 1)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Pressure ratio")), 1, 4)
        self.razonCalculada = Entrada_con_unidades(float, readOnly=True)
        lyt.addWidget(self.razonCalculada, 1, 5)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Efficiency")), 2, 4)
        self.rendimientoCalculado = Entrada_con_unidades(float, readOnly=True)
        lyt.addWidget(self.rendimientoCalculado, 2, 5)

        # Cost tab
        lyt_Cost = QtWidgets.QGridLayout(self.tabCostos)
        self.Costos = CostData(self.Equipment)
        self.Costos.valueChanged.connect(self.changeParamsCoste)
        lyt_Cost.addWidget(self.Costos, 1, 0, 1, 2)

        lyt_Cost.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 2, 0, 1, 2)
        group = QtWidgets.QGroupBox(self.tr("Stimated Costs"))
        lyt_Cost.addWidget(group, 3, 0, 1, 2)
        lyt = QtWidgets.QGridLayout(group)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Purchase cost")), 0, 1)
        self.C_adq = Entrada_con_unidades(Currency, retornar=False,
                                          readOnly=True)
        lyt.addWidget(self.C_adq, 0, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Installed cost")), 1, 1)
        self.C_inst = Entrada_con_unidades(Currency, retornar=False,
                                           readOnly=True)
        lyt.addWidget(self.C_inst, 1, 2)
        lyt_Cost.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 4, 0, 1, 2)

        self.on_tipoCalculo_currentIndexChanged(0)
        if equipment:
            self.setEquipment(equipment)

    def on_tipoCalculo_currentIndexChanged(self, idx):
        """Enabled or disabled widget for data entry to calculate"""
        if idx == 0:
            self.trabajo.setReadOnly(True)
            self.razon.setReadOnly(True)
            self.Pout.setReadOnly(False)
            self.rendimiento.setReadOnly(False)
            self.Pout.setResaltado(True)
            self.rendimiento.setResaltado(True)
            self.razon.setResaltado(False)
            self.trabajo.setResaltado(False)
        elif idx == 1:
            self.Pout.setReadOnly(True)
            self.razon.setReadOnly(False)
            self.trabajo.setReadOnly(True)
            self.rendimiento.setReadOnly(False)
            self.Pout.setResaltado(False)
            self.rendimiento.setResaltado(True)
            self.razon.setResaltado(True)
            self.trabajo.setResaltado(False)
        elif idx == 2:
            self.trabajo.setReadOnly(False)
            self.Pout.setReadOnly(True)
            self.razon.setReadOnly(True)
            self.rendimiento.setReadOnly(False)
            self.Pout.setResaltado(False)
            self.rendimiento.setResaltado(True)
            self.razon.setResaltado(False)
            self.trabajo.setResaltado(True)
        elif idx == 3:
            self.rendimiento.setReadOnly(True)
            self.razon.setReadOnly(True)
            self.trabajo.setReadOnly(False)
            self.Pout.setReadOnly(False)
            self.Pout.setResaltado(True)
            self.rendimiento.setResaltado(False)
            self.razon.setResaltado(False)
            self.trabajo.setResaltado(True)
        elif idx == 4:
            self.rendimiento.setReadOnly(True)
            self.razon.setReadOnly(False)
            self.trabajo.setReadOnly(False)
            self.Pout.setReadOnly(True)
            self.Pout.setResaltado(False)
            self.rendimiento.setResaltado(False)
            self.razon.setResaltado(True)
            self.trabajo.setResaltado(True)
        else:
            self.rendimiento.setReadOnly(False)
            self.razon.setReadOnly(False)
            self.trabajo.setReadOnly(False)
            self.Pout.setReadOnly(False)
            self.Pout.setResaltado(True)
            self.rendimiento.setResaltado(True)
            self.razon.setResaltado(True)
            self.trabajo.setResaltado(True)

        self.changeParams("metodo", idx)

    def rellenar(self):
        UI_equip.rellenar(self)
        if self.Equipment.status == 1 and self.metodo.currentIndex() == 5:
            self.entrada.setCorriente(self.Equipment.entrada)


if __name__ == "__main__":
    import sys
    from lib.corriente import Corriente
    app = QtWidgets.QApplication(sys.argv)
    agua = Corriente(T=500, P=3*101325, caudalMasico=1, fraccionMasica=[1.])
    turbina = Turbine(entrada=agua, metodo=1, razon=0.3, rendimiento=1.)
    dialogo = UI_equipment(turbina)
    dialogo.show()
    sys.exit(app.exec())
