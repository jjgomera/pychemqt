#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


###############################################################################
# Grinder equipment dialog
###############################################################################


from functools import partial


from equipment.parents import UI_equip
from equipment.solids import Grinder
from lib import unidades
from tools import costIndex
from tools.qt import QtWidgets
from UI.widgets import Entrada_con_unidades


class UI_equipment(UI_equip):
    """Diálogo de definición de molinos trituradores de sólidos"""
    Equipment = Grinder()

    def __init__(self, equipment=None, parent=None):
        """equipment: Initial equipment instance to model"""
        super().__init__(Grinder, entrada=False, salida=False, parent=parent)

        # Calculate tab
        lyt = QtWidgets.QGridLayout(self.tabCalculo)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Bond work index")), 1, 0, 1, 1)
        self.Material = QtWidgets.QComboBox()
        self.Material.addItem(self.tr("User defined"))
        for key, idx in self.Equipment.BOND_INDEX:
            self.Material.addItem(key)
        self.Material.currentIndexChanged.connect(self.cambiarBondWordIndex)
        lyt.addWidget(self.Material, 1, 1, 1, 1)
        self.BondIndex = Entrada_con_unidades(float)
        self.BondIndex.valueChanged.connect(
            partial(self.changeParams, "BondIndex"))
        lyt.addWidget(self.BondIndex, 1, 2, 1, 1)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Exponent")), 2, 0, 1, 1)
        self.exponent = Entrada_con_unidades(float)
        self.exponent.valueChanged.connect(
            partial(self.changeParams, "exponent"))
        lyt.addWidget(self.exponent, 2, 1, 1, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("D80")), 3, 0, 1, 1)
        self.D80 = Entrada_con_unidades(unidades.Length, "ParticleDiameter")
        self.D80.valueChanged.connect(
            partial(self.changeParams, "D80"))
        lyt.addWidget(self.D80, 3, 1, 1, 2)

        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 4, 0, 1, 5)

        group = QtWidgets.QGroupBox(self.tr("Results"))
        layout = QtWidgets.QGridLayout(group)
        layout.addWidget(QtWidgets.QLabel(self.tr("Power")), 1, 0)
        self.power = Entrada_con_unidades(
            unidades.Power, retornar=False, readOnly=True)
        layout.addWidget(self.power, 1, 1)
        layout.addWidget(QtWidgets.QLabel(self.tr("Solid Mass Flow")), 2, 0)
        self.solidflow = Entrada_con_unidades(
            unidades.MassFlow, retornar=False)
        self.solidflow.setReadOnly(True)
        layout.addWidget(self.solidflow, 2, 1)
        lyt.addWidget(group, 5, 0, 1, 5)

        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 10, 0, 1, 5)

        # Cost tab
        lyt = QtWidgets.QGridLayout(self.tabCostos)
        lyt.addWidget(
            QtWidgets.QLabel(self.tr("Type:")), 1, 1, 1, 1)
        self.tipo = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_TIPO_COSTOS:
            self.tipo.addItem(txt)
        self.tipo.currentIndexChanged.connect(
            partial(self.changeParamsCoste, "tipoCoste"))
        lyt.addWidget(self.tipo, 1, 2, 1, 1)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 2, 1, 1, 2)

        self.Costos = costIndex.CostData(self.Equipment)
        self.Costos.valueChanged.connect(self.changeParamsCoste)
        lyt.addWidget(self.Costos, 4, 1, 2, 5)

        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 6, 1, 1, 6)
        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 10, 1, 1, 6)
        groupBox_Costos = QtWidgets.QGroupBox(self.tr("Stimated costs"))
        lyt.addWidget(groupBox_Costos, 7, 1, 1, 6)
        lytgroup = QtWidgets.QGridLayout(groupBox_Costos)
        lytgroup.addWidget(QtWidgets.QLabel(self.tr("Purchase cost")), 0, 1)
        self.C_adq = Entrada_con_unidades(
            unidades.Currency, retornar=False, readOnly=True)
        lytgroup.addWidget(self.C_adq, 0, 2)
        lytgroup.addWidget(QtWidgets.QLabel(self.tr("Installed cost")), 1, 1)
        self.C_inst = Entrada_con_unidades(
            unidades.Currency, retornar=False, readOnly=True)
        lytgroup.addWidget(self.C_inst, 1, 2)

        if equipment:
            self.setEquipment(equipment)

    def cambiarBondWordIndex(self, idx):
        try:
            value = self.Equipment.BOND_INDEX[idx+1][1]
        except KeyError:
            self.BondIndex.setReadOnly(False)
        else:
            self.BondIndex.setValue(value)
            self.BondIndex.setReadOnly(True)
            self.changeParams("BondIndex", value)


if __name__ == "__main__":
    import sys
    from lib.corriente import Corriente
    from lib.solids import Solid
    app = QtWidgets.QApplication(sys.argv)

    dm = [17.5, 22.4, 26.2, 31.8, 37, 42.4, 48, 54,
          60, 69, 81.3, 96.5, 109, 127]
    fracciones = [0.02, 0.03, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                  0.05, 0.03, 0.02]
    solid = Solid(caudalSolido=[1], distribucion_diametro=dm,
                  distribucion_fraccion=fracciones, solids=[638])
    grinder = Grinder(entrada=Corriente(solido=solid), D80=1e-5, BondIndex=10)
    dialogo = UI_equipment(grinder)
    dialogo.show()
    sys.exit(app.exec())
