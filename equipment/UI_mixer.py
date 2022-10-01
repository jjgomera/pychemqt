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
# Mixer equipment dialog
###############################################################################


from functools import partial

from PyQt5 import QtWidgets

from lib.unidades import Pressure
from equipment.parents import UI_equip
from equipment.flux import Mixer
from UI import UI_corriente
from UI.widgets import Entrada_con_unidades


class UI_equipment(UI_equip):
    """Mixer equipment edition dialog"""
    Equipment = Mixer()

    def __init__(self, equipment=None, entradas=1, parent=None):
        """
        equipment: Initial equipment instance to model
        entradas: Stream Input number to equipment
        """
        super().__init__(Mixer, salida=False, parent=parent)

        # Input tab
        for i in range(entradas):
            entrada = UI_corriente.Ui_corriente()
            entrada.Changed.connect(partial(self.cambiarEntrada, i))
            self.Entrada.addTab(entrada, str(i+1))

        # Calculate tab
        lyt_Calc = QtWidgets.QGridLayout(self.tabCalculo)
        lyt_Calc.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Output Pressure Method")), 1, 1)
        self.criterio = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_METODO:
            self.criterio.addItem(txt)
        self.criterio.currentIndexChanged.connect(self.criterio_Changed)
        lyt_Calc.addWidget(self.criterio, 1, 2)

        lyt_Calc.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            2, 1, 1, 3)
        lyt_Calc.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Output Pressure")), 3, 1)
        self.Pout = Entrada_con_unidades(Pressure)
        self.Pout.valueChanged.connect(partial(self.changeParams, "Pout"))
        lyt_Calc.addWidget(self.Pout, 3, 2)
        lyt_Calc.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 4, 1, 1, 3)

        self.criterio_Changed(0)

        if equipment:
            self.setEquipment(equipment)
        else:
            self.Equipment = Mixer(entradas=entradas)

    def criterio_Changed(self, int):
        self.Pout.setEnabled(int == 2)
        self.changeParams("criterio", int)

    def cambiarEntrada(self, ind, corriente):
        self.Equipment(id_entrada=ind, entrada=corriente)

    def rellenarInput(self):
        UI_equip.rellenarInput(self)
        for i, entrada in enumerate(self.Equipment.kwargs["entrada"]):
            if entrada:
                self.Entrada.widget(i).setCorriente(entrada)


if __name__ == "__main__":
    import sys
    from lib.corriente import Corriente
    app = QtWidgets.QApplication(sys.argv)
    agua = Corriente(T=300, P=101325., caudalMasico=1, ids=[62],
                     fraccionMasica=[1.])
    agua2 = Corriente(T=300, P=101325.*2, caudalMasico=2, ids=[62],
                      fraccionMasica=[1.])
#    mezclador = Mixer(entrada=[agua, agua2], criterio=0)
    mezclador = Mixer(criterio=0)
    dialogo = UI_equipment(mezclador, entradas=2)
    dialogo.show()
    sys.exit(app.exec())
