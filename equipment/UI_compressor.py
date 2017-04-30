#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2017, Juan José Gómez Romera <jjgomera@gmail.com>

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
# Compressor equipment dialog
###############################################################################


from functools import partial

from PyQt5 import QtWidgets

from lib.unidades import Pressure, Power, Currency
from tools.costIndex import CostData
from equipment.parents import UI_equip
from equipment.compressor import Compressor
from UI.widgets import Entrada_con_unidades


class UI_equipment(UI_equip):
    """Compressor equipment edition dialog"""
    Equipment = Compressor()

    def __init__(self, equipment=None, parent=None):
        """
        equipment: Initial equipment instance to model
        """
        super().__init__(Compressor, entrada=False, salida=False,
                         parent=parent)

        # Calculate tab
        lyt_Calc = QtWidgets.QGridLayout(self.tabCalculo)
        lyt_Calc.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Method:")), 1, 1)
        self.metodo = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_METODO:
            self.metodo.addItem(txt)
        # self.metodo.addItem(QtWidgets.QApplication.translate(
            # "pychemqt", "Especificar curva de funcionamiento"))
        self.metodo.currentIndexChanged.connect(
            self.on_tipoCalculo_currentIndexChanged)
        lyt_Calc.addWidget(self.metodo, 1, 2, 1, 2)
        lyt_Calc.addItem(QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            2, 0, 1, 4)
        lyt_Calc.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Thermodynamic:")), 3, 1)
        self.termodinamica = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_TERMODINAMICA:
            self.termodinamica.addItem(txt)
        self.termodinamica.currentIndexChanged.connect(
            partial(self.changeParams, "termodinamica"))
        lyt_Calc.addWidget(self.termodinamica, 3, 2, 1, 2)
        lyt_Calc.addItem(QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            4, 0, 1, 4)
        lyt_Calc.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Out Pressure")), 5, 1)
        self.Pout = Entrada_con_unidades(Pressure)
        self.Pout.valueChanged.connect(partial(self.changeParams, "Pout"))
        lyt_Calc.addWidget(self.Pout, 5, 2)
        lyt_Calc.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Pressure ratio")), 6, 1)
        self.razon = Entrada_con_unidades(float)
        self.razon.valueChanged.connect(partial(self.changeParams, "razon"))
        lyt_Calc.addWidget(self.razon, 6, 2)
        lyt_Calc.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Efficiency")), 7, 1)
        self.rendimiento = Entrada_con_unidades(float)
        self.rendimiento.valueChanged.connect(
            partial(self.changeParams, "rendimiento"))
        lyt_Calc.addWidget(self.rendimiento, 7, 2)
        lyt_Calc.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Actual Power")), 8, 1)
        self.trabajo = Entrada_con_unidades(Power)
        self.trabajo.valueChanged.connect(
            partial(self.changeParams, "trabajo"))
        lyt_Calc.addWidget(self.trabajo, 8, 2)
        lyt_Calc.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Stages")), 9, 1)
        self.etapas = Entrada_con_unidades(
            int, spinbox=True, min=1, value=1, step=1)
        self.etapas.valueChanged.connect(partial(self.changeParams, "etapas"))
        lyt_Calc.addWidget(self.etapas, 9, 2)
        lyt_Calc.setRowStretch(10, 1)

        group = QtWidgets.QGroupBox()
        group.setTitle(QtWidgets.QApplication.translate("pychemqt", "Results"))
        lyt_Calc.addWidget(group, 12, 1, 1, 2)
        lyt = QtWidgets.QGridLayout(group)
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Power")), 1, 1)
        self.power = Entrada_con_unidades(Power, retornar=False, readOnly=True)
        lyt.addWidget(self.power, 1, 2)
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Cp/Cv ratio")), 2, 1)
        self.cp_cv = Entrada_con_unidades(float, retornar=False, readOnly=True)
        lyt.addWidget(self.cp_cv, 2, 2)
        lyt.setColumnStretch(3, 1)
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Pressure ratio")), 1, 4)
        self.razonCalculada = Entrada_con_unidades(float, readOnly=True)
        lyt.addWidget(self.razonCalculada, 1, 5)
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Efficiency")), 2, 4)
        self.rendimientoCalculado = Entrada_con_unidades(float, readOnly=True)
        lyt.addWidget(self.rendimientoCalculado, 2, 5)

        # Cost tab
        lyt_Cost = QtWidgets.QGridLayout(self.tabCostos)
        lyt_Cost.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Compressor:")), 1, 1)
        self.compresor = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_COMPRESOR:
            self.compresor.addItem(txt)
        self.compresor.currentIndexChanged.connect(
            partial(self.changeParamsCoste, "compresor"))
        lyt_Cost.addWidget(self.compresor, 1, 2)
        lyt_Cost.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Driver:")), 2, 1)
        self.transmision = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_TRANSMISION:
            self.transmision.addItem(txt)
        self.transmision.currentIndexChanged.connect(
            partial(self.changeParamsCoste, "transmision"))
        lyt_Cost.addWidget(self.transmision, 2, 2)
        lyt_Cost.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Motor:")), 3, 1)
        self.motor = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_MOTOR:
            self.motor.addItem(txt)
        self.motor.currentIndexChanged.connect(
            partial(self.changeParamsCoste, "motor"))
        lyt_Cost.addWidget(self.motor, 3, 2)
        lyt_Cost.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "RPM:")), 4, 1)
        self.rpm = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_RPM:
            self.rpm.addItem(txt)
        self.rpm.currentIndexChanged.connect(
            partial(self.changeParamsCoste, "rpm"))
        lyt_Cost.addWidget(self.rpm, 4, 2)
        lyt_Cost.addItem(QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Fixed,
            QtWidgets.QSizePolicy.Maximum), 5, 1, 1, 6)

        self.Costos = CostData(self.Equipment)
        self.Costos.valueChanged.connect(self.changeParamsCoste)
        lyt_Cost.addWidget(self.Costos, 6, 1, 1, 3)

        lyt_Cost.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 8, 1, 1, 6)
        lyt_Cost.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 10, 1, 1, 6)
        group = QtWidgets.QGroupBox()
        group.setTitle(
            QtWidgets.QApplication.translate("pychemqt", "Stimated Costs"))
        lyt_Cost.addWidget(group, 9, 1, 1, 5)
        lyt = QtWidgets.QGridLayout(group)
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Compressor")), 0, 0)
        self.C_comp = Entrada_con_unidades(Currency, retornar=False)
        self.C_comp.setReadOnly(True)
        lyt.addWidget(self.C_comp, 0, 1)
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Drive")), 1, 0)
        self.C_trans = Entrada_con_unidades(Currency, retornar=False)
        self.C_trans.setReadOnly(True)
        lyt.addWidget(self.C_trans, 1, 1)
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Motor")), 2, 0)
        self.C_motor = Entrada_con_unidades(Currency, retornar=False)
        self.C_motor.setReadOnly(True)
        lyt.addWidget(self.C_motor, 2, 1)
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Purchase cost")), 0, 4)
        self.C_adq = Entrada_con_unidades(Currency, retornar=False)
        self.C_adq.setReadOnly(True)
        lyt.addWidget(self.C_adq, 0, 5)
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Installed cost")), 1, 4)
        self.C_inst = Entrada_con_unidades(Currency, retornar=False)
        self.C_inst.setReadOnly(True)
        lyt.addWidget(self.C_inst, 1, 5)

        self.on_tipoCalculo_currentIndexChanged(0)
        if equipment:
            self.setEquipment(equipment)

    def on_tipoCalculo_currentIndexChanged(self, int):
        """Enabled or disabled widget for data entry to calculate"""
        if int == 0:
            self.trabajo.setReadOnly(True)
            self.razon.setReadOnly(True)
            self.Pout.setReadOnly(False)
            self.rendimiento.setReadOnly(False)
            self.Pout.setResaltado(True)
            self.rendimiento.setResaltado(True)
            self.razon.setResaltado(False)
            self.trabajo.setResaltado(False)
        elif int == 1:
            self.Pout.setReadOnly(True)
            self.razon.setReadOnly(False)
            self.trabajo.setReadOnly(True)
            self.rendimiento.setReadOnly(False)
            self.Pout.setResaltado(False)
            self.rendimiento.setResaltado(True)
            self.razon.setResaltado(True)
            self.trabajo.setResaltado(False)
        elif int == 2:
            self.trabajo.setReadOnly(False)
            self.Pout.setReadOnly(True)
            self.razon.setReadOnly(True)
            self.rendimiento.setReadOnly(False)
            self.Pout.setResaltado(False)
            self.rendimiento.setResaltado(True)
            self.razon.setResaltado(False)
            self.trabajo.setResaltado(True)
        elif int == 3:
            self.rendimiento.setReadOnly(True)
            self.razon.setReadOnly(True)
            self.trabajo.setReadOnly(False)
            self.Pout.setReadOnly(False)
            self.Pout.setResaltado(True)
            self.rendimiento.setResaltado(False)
            self.razon.setResaltado(False)
            self.trabajo.setResaltado(True)
        elif int == 4:
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
        self.changeParams("metodo", int)

    def rellenar(self):
        UI_equip.rellenar(self)
        if self.Equipment.status == 1 and self.metodo.currentIndex() == 5:
            self.entrada.setCorriente(self.Equipment.entrada)


if __name__ == "__main__":
    import sys
    from lib.corriente import Corriente
    app = QtWidgets.QApplication(sys.argv)
    agua = Corriente(T=400, P=101325., caudalMasico=1, fraccionMasica=[1.])
    compresor = Compressor(entrada=agua, Pout=5*101325., rendimiento=0.75,
                           compresor=1, Current_index=1000)
    dialogo = UI_equipment(compresor)
    dialogo.show()
    sys.exit(app.exec_())
