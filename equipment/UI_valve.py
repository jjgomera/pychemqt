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
# Valve equipment definition
###############################################################################


from functools import partial

from tools.qt import QtWidgets

from lib.unidades import Temperature, Pressure
from equipment.parents import UI_equip
from equipment.flux import Valve
from UI.widgets import Entrada_con_unidades


class UI_equipment (UI_equip):
    """Valve equipment edition dialog"""
    Equipment = Valve()

    def __init__(self, equipment=None, parent=None):
        """
        equipment: Initial equipment instance to model
        """
        super().__init__(Valve, entrada=False, salida=False, parent=parent)

        # Calculate tab
        lyt_Calc = QtWidgets.QGridLayout(self.tabCalculo)
        lyt_Calc.addWidget(QtWidgets.QLabel(self.tr("Valve operation")), 1, 1)
        self.off = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_WORKING:
            self.off.addItem(txt)
        self.off.currentIndexChanged.connect(self.criterio_Changed)
        lyt_Calc.addWidget(self.off, 1, 2)
        lyt_Calc.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 2, 1, 1, 6)
        lyt_Calc.addWidget(QtWidgets.QLabel(self.tr("Pout")), 3, 1)
        self.Pout = Entrada_con_unidades(Pressure)
        self.Pout.valueChanged.connect(partial(self.changeParams, "Pout"))
        lyt_Calc.addWidget(self.Pout, 3, 2)
        lyt_Calc.addWidget(QtWidgets.QLabel(self.tr("DeltaP")), 4, 1)
        self.DeltaP = Entrada_con_unidades(Pressure)
        self.DeltaP.valueChanged.connect(partial(self.changeParams, "DeltaP"))
        lyt_Calc.addWidget(self.DeltaP, 4, 2)
        lyt_Calc.addWidget(QtWidgets.QLabel(self.tr("T dew point")), 5, 1)
        self.Dew = Entrada_con_unidades(Temperature)
        self.Dew.valueChanged.connect(partial(self.changeParams, "Dew"))
        lyt_Calc.addWidget(self.Dew, 5, 2)
        lyt_Calc.addWidget(QtWidgets.QLabel(self.tr("T bubble point")), 6, 1)
        self.Bubble = Entrada_con_unidades(Temperature)
        self.Bubble.valueChanged.connect(partial(self.changeParams, "Bubble"))
        lyt_Calc.addWidget(self.Bubble, 6, 2)
        lyt_Calc.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 10, 1, 1, 6)

        self.criterio_Changed(0)
        if equipment:
            self.setEquipment(equipment)

    def criterio_Changed(self, int):
        self.Pout.setEnabled(int == 1)
        self.DeltaP.setEnabled(int == 1)
        self.Dew.setEnabled(int == 1)
        self.Bubble.setEnabled(int == 1)
        self.calculo(off=int)


if __name__ == "__main__":
    import sys
    from lib.corriente import Corriente
    app = QtWidgets.QApplication(sys.argv)
    agua = Corriente(T=300, P=101325, caudalMasico=1, fraccionMasica=[1.])
    valvula = Valve(entrada=agua, off=1, DeltaP=1000)
    dialogo = UI_equipment(valvula)
    dialogo.show()
    sys.exit(app.exec())
