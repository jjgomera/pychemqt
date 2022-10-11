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
# Fireheater equipment dialog
###############################################################################


from functools import partial

from qt import QtWidgets

from lib.unidades import Temperature, Pressure, Power, VolFlow, Currency
from tools.costIndex import CostData
from equipment.parents import UI_equip
from equipment.heatExchanger import Fired_Heater
from UI.widgets import Entrada_con_unidades


class UI_equipment(UI_equip):
    """Fireheater equipment edition dialog"""
    Equipment = Fired_Heater()

    def __init__(self, equipment=None, parent=None):
        """
        equipment: Initial equipment instance to model
        """
        super().__init__(Fired_Heater, entrada=False, salida=False,
                         parent=parent)

        # Calculate tab
        layout = QtWidgets.QGridLayout(self.tabCalculo)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Output Temperature")), 1, 1)
        self.Tout = Entrada_con_unidades(Temperature, resaltado=True)
        self.Tout.valueChanged.connect(partial(self.changeParams, "Tout"))
        layout.addWidget(self.Tout, 1, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy. Fixed, QtWidgets.QSizePolicy.Policy.Fixed),
            2, 0, 1, 6)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Pressure drop")), 3, 1)
        self.deltaP = Entrada_con_unidades(Pressure)
        self.deltaP.valueChanged.connect(partial(self.changeParams, "deltaP"))
        layout.addWidget(self.deltaP, 3, 2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Maximum heat flux")), 4, 1)
        self.Hmax = Entrada_con_unidades(Power)
        self.Hmax.valueChanged.connect(partial(self.changeParams, "Hmax"))
        layout.addWidget(self.Hmax, 4, 2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Fuel calorific value")), 5, 1)
        self.poderCalorifico = Entrada_con_unidades(float)
        self.poderCalorifico.valueChanged.connect(
            partial(self.changeParams, "poderCalorifico"))
        layout.addWidget(self.poderCalorifico, 5, 2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Efficiency")), 6, 1)
        self.eficiencia = Entrada_con_unidades(float, spinbox=True)
        self.eficiencia.valueChanged.connect(
            partial(self.changeParams, "eficiencia"))
        layout.addWidget(self.eficiencia, 6, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 7, 0, 1, 6)

        group = QtWidgets.QGroupBox(
            QtWidgets.QApplication.translate("pychemqt", "Results"))
        layout.addWidget(group, 8, 1, 1, 5)
        lyt = QtWidgets.QGridLayout(group)
        lyt.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Heat")), 0, 1)
        self.Heat = Entrada_con_unidades(Power, retornar=False, readOnly=True)
        lyt.addWidget(self.Heat, 0, 2)
        lyt.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Fuel")), 1, 1)
        self.CombustibleRequerido = Entrada_con_unidades(
            VolFlow, "QLiq", retornar=False, readOnly=True)
        lyt.addWidget(self.CombustibleRequerido, 1, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 9, 0, 1, 6)

        # Cost tab
        lyt_Cost = QtWidgets.QGridLayout(self.tabCostos)
        lyt_Cost.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Type")), 1, 1)
        self.tipo = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_TIPO:
            self.tipo.addItem(txt)
        self.tipo.currentIndexChanged.connect(self.mostrarSubclasificacion)
        lyt_Cost.addWidget(self.tipo, 1, 2)
        self.label = QtWidgets.QLabel()
        lyt_Cost.addWidget(self.label, 2, 1)
        self.subtipoBox = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_SUBTIPOBOX:
            self.subtipoBox.addItem(txt)
        self.subtipoBox.currentIndexChanged.connect(
            partial(self.changeParamsCoste, "subtipoBox"))
        lyt_Cost.addWidget(self.subtipoBox, 2, 2)
        self.subtipoCylindrical = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_SUBTIPOCYLINDRICAL:
            self.subtipoCylindrical.addItem(txt)
        self.subtipoCylindrical.currentIndexChanged.connect(
            partial(self.changeParamsCoste, "subtipoCylindrical"))
        lyt_Cost.addWidget(self.subtipoCylindrical, 2, 2)
        lyt_Cost.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Material")), 3, 1)
        self.material = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_MATERIAL:
            self.material.addItem(txt)
        self.material.currentIndexChanged.connect(
            partial(self.changeParamsCoste, "material"))
        lyt_Cost.addWidget(self.material, 3, 2)
        lyt_Cost.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Design pressure")), 4, 1)
        self.P_dis = Entrada_con_unidades(Pressure)
        self.P_dis.valueChanged.connect(
            partial(self.changeParamsCoste, "P_dis"))
        lyt_Cost.addWidget(self.P_dis, 4, 2)
        lyt_Cost.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 5, 1, 1, 6)

        self.Costos = CostData(self.Equipment)
        self.Costos.valueChanged.connect(self.changeParamsCoste)
        lyt_Cost.addWidget(self.Costos, 6, 1, 2, 5)

        lyt_Cost.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 8, 1, 1, 6)
        group = QtWidgets.QGroupBox(
            QtWidgets.QApplication.translate("pychemqt", "Stimated Costs"))
        lyt_Cost.addWidget(group, 9, 1, 1, 6)
        layout = QtWidgets.QGridLayout(group)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Purchase costs")), 0, 1)
        self.C_adq = Entrada_con_unidades(Currency, retornar=False)
        self.C_adq.setReadOnly(True)
        layout.addWidget(self.C_adq, 0, 2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Installed costs")), 1, 1)
        self.C_inst = Entrada_con_unidades(Currency, retornar=False)
        self.C_inst.setReadOnly(True)
        layout.addWidget(self.C_inst, 1, 2)
        lyt_Cost.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 10, 1, 1, 6)

        self.mostrarSubclasificacion(0)
        if equipment:
            self.setEquipment(equipment)

    def mostrarSubclasificacion(self, ind):
        if ind:
            txt = QtWidgets.QApplication.translate(
                "pychemqt", "Cylindrical heater type")
        else:
            txt = QtWidgets.QApplication.translate(
                "pychemqt", "Box heater type")
        self.label.setText(txt)
        self.subtipoBox.setVisible(not ind)
        self.subtipoCylindrical.setVisible(ind)
        self.changeParamsCoste("tipo", ind)


if __name__ == "__main__":
    import sys
    from lib.corriente import Corriente
    from equipment.heatExchanger import Fired_Heater
    app = QtWidgets.QApplication(sys.argv)
    agua = Corriente(T=300, P=101325., caudalMasico=0.01, fraccionMasica=[1.])
    fireheater = Fired_Heater(entrada=agua, Tout=450)
    dialogo = UI_equipment(fireheater)
    dialogo.show()
    sys.exit(app.exec())
