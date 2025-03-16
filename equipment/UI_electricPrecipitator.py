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
# Electric precipitator equipment dialog
###############################################################################


from functools import partial

from tools.qt import QtWidgets

from equipment.gas_solid import ElectricPrecipitator
from equipment.parents import UI_equip
from lib.unidades import DeltaP, PotencialElectric, Area
from UI.widgets import Entrada_con_unidades


class UI_equipment(UI_equip):
    """Electric precipitator equipment edition dialog"""
    Equipment = ElectricPrecipitator()

    def __init__(self, equipment=None, parent=None):
        """
        equipment: Initial equipment instance to model
        """
        super().__init__(ElectricPrecipitator, entrada=False, parent=parent)

        # Calculate tab
        lyt_Calc = QtWidgets.QGridLayout(self.tabCalculo)
        lyt_Calc.addWidget(QtWidgets.QLabel(self.tr("Mode")), 1, 1)
        self.metodo = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_TIPO:
            self.metodo.addItem(txt)
        self.metodo.currentIndexChanged.connect(self.tipoCalculoCambiado)
        lyt_Calc.addWidget(self.metodo, 1, 2, 1, 4)
        lyt_Calc.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 2, 1, 1, 6)

        lyt_Calc.addWidget(QtWidgets.QLabel(self.tr("Area")), 3, 1)
        self.area = Entrada_con_unidades(Area, resaltado=True)
        self.area.valueChanged.connect(partial(self.changeParams, "area"))
        lyt_Calc.addWidget(self.area, 3, 2)
        lyt_Calc.addWidget(QtWidgets.QLabel(
            self.tr("Allowable efficiency")), 4, 1)
        self.rendimientoAdmisible = Entrada_con_unidades(float,  readOnly=True)
        self.rendimientoAdmisible.valueChanged.connect(
            partial(self.changeParams, "rendimientoAdmisible"))
        lyt_Calc.addWidget(self.rendimientoAdmisible, 4, 2)
        lyt_Calc.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 5, 1, 1, 6)

        lyt_Calc.addWidget(QtWidgets.QLabel(
            self.tr("Dielectric constant")), 6, 1)
        self.epsilon = Entrada_con_unidades(float)
        self.epsilon.valueChanged.connect(
            partial(self.changeParams, "epsilon"))
        lyt_Calc.addWidget(self.epsilon, 6, 2)
        lyt_Calc.addWidget(QtWidgets.QLabel(self.tr("Charging field")), 7, 1)
        self.potencialCarga = Entrada_con_unidades(PotencialElectric)
        self.potencialCarga.valueChanged.connect(
            partial(self.changeParams, "potencialCarga"))
        lyt_Calc.addWidget(self.potencialCarga, 7, 2)
        lyt_Calc.addWidget(QtWidgets.QLabel(self.tr("Collecting field")), 8, 1)
        self.potencialDescarga = Entrada_con_unidades(PotencialElectric)
        self.potencialDescarga.valueChanged.connect(
            partial(self.changeParams, "potencialDescarga"))
        lyt_Calc.addWidget(self.potencialDescarga, 8, 2)
        lyt_Calc.addWidget(QtWidgets.QLabel(self.tr("Pressure drop")), 9, 1)
        self.deltaP = Entrada_con_unidades(DeltaP)
        self.deltaP.valueChanged.connect(partial(self.changeParams, "deltaP"))
        lyt_Calc.addWidget(self.deltaP, 9, 2)
        lyt_Calc.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 10, 1, 1, 6)

        groupbox = QtWidgets.QGroupBox(self.tr("Result"))
        lyt_Calc.addWidget(groupbox, 11, 1, 1, 5)
        lyt = QtWidgets.QGridLayout(groupbox)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Area")), 0, 1)
        self.areaCalculada = Entrada_con_unidades(Area, retornar=False)
        self.areaCalculada.setReadOnly(True)
        lyt.addWidget(self.areaCalculada, 0, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Efficiency")), 1, 1)
        self.rendimiento = Entrada_con_unidades(float, readOnly=True)
        lyt.addWidget(self.rendimiento, 1, 2)
        lyt_Calc.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 12, 1, 1, 6)

        # Output tab
        self.addSalida(
            self.tr("Filtered gas"))
        self.addSalida(
            self.tr("Collected solids"))

        if equipment:
            self.setEquipment(equipment)

    def tipoCalculoCambiado(self, tipo_calculo):
        self.area.setReadOnly(tipo_calculo)
        self.area.setResaltado(not tipo_calculo)
        self.rendimientoAdmisible.setReadOnly(not tipo_calculo)
        self.rendimientoAdmisible.setResaltado(tipo_calculo)
        self.changeParams("metodo", tipo_calculo)


if __name__ == "__main__":
    import sys
    from lib.corriente import Corriente, Solid
    app = QtWidgets.QApplication(sys.argv)
    diametros = [17.5e-6, 22.4e-6, 26.2e-6, 31.8e-6, 37e-6, 42.4e-6, 48e-6,
                 54e-6, 60e-6, 69e-6, 81.3e-6, 96.5e-6, 109e-6, 127e-6]
    fracciones = [0.02, 0.03, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                  0.05, 0.03, 0.02]
    solido = Solid(caudalSolido=[0.1], distribucion_diametro=diametros,
                   distribucion_fraccion=fracciones)
    corriente = Corriente(T=300, P=101325, caudalMasico=1.,
                          fraccionMolar=[1.], solido=solido)
    precipitador = ElectricPrecipitator(entrada=corriente, metodo=1,
                                        rendimientoAdmisible=0.9)
    dialogo = UI_equipment(precipitador)
    dialogo.show()
    sys.exit(app.exec())
