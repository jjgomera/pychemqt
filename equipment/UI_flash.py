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
# Flash phase separator equipment dialog
###############################################################################

from functools import partial


from equipment.distillation import Flash
from equipment.parents import UI_equip
from lib.unidades import Length, Mass, Volume, Density, Currency
from tools.qt import QtWidgets
from tools.costIndex import CostData
from UI.widgets import Entrada_con_unidades


class UI_equipment (UI_equip):
    """Flash phase separator equipment edition dialog"""
    Equipment = Flash()

    def __init__(self, equipment=None, parent=None):
        """
        equipment: Initial equipment instance to model
        """
        super(UI_equipment, self).__init__(Flash, entrada=False, parent=parent)

        # Calculate tab
        lyt_Calc = QtWidgets.QGridLayout(self.tabCalculo)
        lyt_Calc.addWidget(QtWidgets.QLabel(self.tr("Method")), 0, 1)
        self.flash = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_FLASH:
            self.flash.addItem(txt)
        self.flash.currentIndexChanged.connect(
            partial(self.changeParams, "metodo"))
        lyt_Calc.addWidget(self.flash, 0, 2)
        lyt_Calc.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 1, 1, 1, 6)

        # Cost tab
        lyt_Cost = QtWidgets.QGridLayout(self.tabCostos)
        lyt_Cost.addWidget(QtWidgets.QLabel(self.tr("Orientation")), 0, 1)
        self.orientacion = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_ORIENTATION:
            self.orientacion.addItem(txt)
        self.orientacion.currentIndexChanged.connect(
            partial(self.changeParamsCoste, "orientacion"))
        lyt_Cost.addWidget(self.orientacion, 0, 2)
        lyt_Cost.addWidget(QtWidgets.QLabel(self.tr("Material")), 1, 1)
        self.material = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_MATERIAL:
            self.material.addItem(txt)
        self.material.currentIndexChanged.connect(
            partial(self.changeParamsCoste, "material"))
        lyt_Cost.addWidget(self.material, 1, 2, 1, 4)
        lyt_Cost.addWidget(QtWidgets.QLabel(self.tr("Density")), 2, 4)
        self.Densidad = Entrada_con_unidades(Density, "DenLiq")
        self.Densidad.valueChanged.connect(
            partial(self.changeParamsCoste, "densidad"))
        lyt_Cost.addWidget(self.Densidad, 2, 5)
        lyt_Cost.addWidget(QtWidgets.QLabel(self.tr("Diameter")), 2, 1)
        self.diametro = Entrada_con_unidades(Length)
        self.diametro.valueChanged.connect(
            partial(self.changeParamsCoste, "diametro"))
        lyt_Cost.addWidget(self.diametro, 2, 2)
        lyt_Cost.addWidget(QtWidgets.QLabel(self.tr("Length")), 3, 1)
        self.longitud = Entrada_con_unidades(Length)
        self.longitud.valueChanged.connect(
            partial(self.changeParamsCoste, "longitud"))
        lyt_Cost.addWidget(self.longitud, 3, 2)
        lyt_Cost.addWidget(QtWidgets.QLabel(self.tr("Thickness")), 4, 1)
        self.espesor = Entrada_con_unidades(Length, "Thickness")
        self.espesor.valueChanged.connect(
            partial(self.changeParamsCoste, "espesor"))
        lyt_Cost.addWidget(self.espesor, 4, 2)
        lyt_Cost.addWidget(QtWidgets.QLabel(self.tr("Head type")), 5, 1)
        self.cabeza = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_HEAD:
            self.cabeza.addItem(txt)
        self.cabeza.currentIndexChanged.connect(
            partial(self.changeParamsCoste, "cabeza"))
        lyt_Cost.addWidget(self.cabeza, 5, 2)
        lyt_Cost.addWidget(QtWidgets.QLabel(self.tr("Head Thickness")), 6, 1)
        self.espesor_cabeza = Entrada_con_unidades(Length, "Thickness")
        self.espesor_cabeza.valueChanged.connect(
            partial(self.changeParamsCoste, "espesor_cabeza"))
        lyt_Cost.addWidget(self.espesor_cabeza, 6, 2)
        lyt_Cost.addWidget(QtWidgets.QLabel(
            self.tr("Straight flange length")), 7, 1)
        self.reborde = Entrada_con_unidades(Length)
        self.reborde.valueChanged.connect(
            partial(self.changeParamsCoste, "reborde"))
        lyt_Cost.addWidget(self.reborde, 7, 2)
        lyt_Cost.addWidget(QtWidgets.QLabel(self.tr("Volume")), 6, 4)
        self.Volumen = Entrada_con_unidades(Volume, "VolLiq", retornar=False)
        self.Volumen.setReadOnly(True)
        lyt_Cost.addWidget(self.Volumen, 6, 5)
        lyt_Cost.addWidget(QtWidgets.QLabel(self.tr("Weight")), 7, 4)
        self.Peso = Entrada_con_unidades(Mass, readOnly=True)
        lyt_Cost.addWidget(self.Peso, 7, 5)
        lyt_Cost.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 2, 3, 6, 1)
        lyt_Cost.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 8, 0, 1, 6)

        self.Costos = CostData(self.Equipment)
        self.Costos.valueChanged.connect(self.changeParamsCoste)
        lyt_Cost.addWidget(self.Costos, 9, 1, 2, 5)

        lyt_Cost.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 11, 0, 1, 6)
        group = QtWidgets.QGroupBox(self.tr("Stimated Costs"))
        lyt_Cost.addWidget(group, 12, 1, 1, 5)
        layout = QtWidgets.QGridLayout(group)
        layout.addWidget(QtWidgets.QLabel(self.tr("Purchase costs")), 0, 1)
        self.C_adq = Entrada_con_unidades(
            Currency, retornar=False, tolerancia=8, decimales=2)
        self.C_adq.setReadOnly(True)
        layout.addWidget(self.C_adq, 0, 2)
        layout.addWidget(QtWidgets.QLabel(self.tr("Installed costs")), 1, 1)
        self.C_inst = Entrada_con_unidades(
            Currency, retornar=False, tolerancia=8, decimales=2)
        self.C_inst.setReadOnly(True)
        layout.addWidget(self.C_inst, 1, 2)
        lyt_Cost.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 13, 0, 1, 6)

        # Output tab
        self.addSalida(self.tr("Destilate"))
        self.addSalida(self.tr("Residue"))

        if equipment:
            self.setEquipment(equipment)


if __name__ == "__main__":
    import sys
    from lib.corriente import Corriente
    app = QtWidgets.QApplication(sys.argv)
    entrada = Corriente(T=340, P=101325, caudalMasico=0.01,
                        ids=[10, 38, 22, 61],
                        fraccionMolar=[.3, 0.5, 0.05, 0.15])
    flash = Flash(entrada=entrada)
    dialogo = UI_equipment(flash)
    dialogo.show()
    sys.exit(app.exec())
