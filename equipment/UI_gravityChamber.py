#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2016, Juan José Gómez Romera <jjgomera@gmail.com>

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
# gravity chamber equipment dialog
###############################################################################

from functools import partial

from PyQt5 import QtWidgets


from lib.unidades import Length, Speed, VolFlow, DeltaP
from .gas_solid import GravityChamber
from equipment.parents import UI_equip
from UI.widgets import Entrada_con_unidades


class UI_equipment(UI_equip):
    """Gravity chamber equipment edition dialog"""
    Equipment = GravityChamber()

    def __init__(self, equipment=None, parent=None):
        """
        equipment: Initial equipment instance to model
        """
        super(UI_equipment, self).__init__(GravityChamber, entrada=False,
                                           parent=parent)

        # Calculate tab
        lyt_Calc = QtWidgets.QGridLayout(self.tabCalculo)
        lyt_Calc.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Mode")), 1, 1, 1, 1)
        self.metodo = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_TIPO:
            self.metodo.addItem(txt)
        self.metodo.currentIndexChanged.connect(self.tipoCalculoCambiado)
        lyt_Calc.addWidget(self.metodo, 1, 2, 1, 4)
        lyt_Calc.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Model")), 2, 1, 1, 1)
        self.modelo = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_MODEL:
            self.modelo.addItem(txt)
        self.modelo.currentIndexChanged.connect(partial(self.changeParams,
                                                        "modelo"))
        lyt_Calc.addWidget(self.modelo, 2, 2, 1, 1)
        lyt_Calc.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            3, 1, 1, 6)
        lyt_Calc.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Width")), 4, 1, 1, 1)
        self.W = Entrada_con_unidades(Length)
        self.W.valueChanged.connect(partial(self.changeParams, "W"))
        lyt_Calc.addWidget(self.W, 4, 2, 1, 1)
        lyt_Calc.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Height")), 5, 1, 1, 1)
        self.H = Entrada_con_unidades(Length)
        self.H.valueChanged.connect(partial(self.changeParams, "H"))
        lyt_Calc.addWidget(self.H, 5, 2, 1, 1)
        lyt_Calc.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Length")), 6, 1, 1, 1)
        self.L = Entrada_con_unidades(Length)
        self.L.valueChanged.connect(partial(self.changeParams, "L"))
        lyt_Calc.addWidget(self.L, 6, 2, 1, 1)
        lyt_Calc.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Allowable efficiency")), 7, 1, 1, 1)
        self.rendimientoAdmisible = Entrada_con_unidades(float, spinbox=True, max=1)
        self.rendimientoAdmisible.valueChanged.connect(
            partial(self.changeParams, "rendimientoAdmisible"))
        lyt_Calc.addWidget(self.rendimientoAdmisible, 7, 2, 1, 1)
        lyt_Calc.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Allowable speed")), 8, 1, 1, 1)
        self.velocidadAdmisible = Entrada_con_unidades(Speed)
        self.velocidadAdmisible.valueChanged.connect(
            partial(self.changeParams, "velocidadAdmisible"))
        lyt_Calc.addWidget(self.velocidadAdmisible, 8, 2, 1, 1)
        lyt_Calc.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Pressure loss")), 9, 1, 1, 1)
        self.deltaP = Entrada_con_unidades(DeltaP)
        self.deltaP.valueChanged.connect(partial(self.changeParams, "deltaP"))
        lyt_Calc.addWidget(self.deltaP, 9, 2, 1, 1)

        lyt_Calc.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding),
            10, 1, 1, 6)
        group_Calc = QtWidgets.QGroupBox(QtWidgets.QApplication.translate(
            "pychemqt", "Results"))
        lyt_Calc.addWidget(group_Calc, 11, 1, 1, 5)
        lyt = QtWidgets.QGridLayout(group_Calc)
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Flow")), 0, 1)
        self.Q = Entrada_con_unidades(VolFlow, "QGas", retornar=False)
        self.Q.setReadOnly(True)
        lyt.addWidget(self.Q, 0, 2)
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "V<sub>gas</sub>")), 1, 1)
        self.Vgas = Entrada_con_unidades(Speed, retornar=False, readOnly=True)
        lyt.addWidget(self.Vgas, 1, 2)
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Efficiency")), 2, 1)
        self.rendimiento = Entrada_con_unidades(float, retornar=False)
        self.rendimiento.setReadOnly(True)
        lyt.addWidget(self.rendimiento, 2, 2)
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Height")), 0, 4)
        self.HCalc = Entrada_con_unidades(Length, retornar=False)
        self.HCalc.setReadOnly(True)
        lyt.addWidget(self.HCalc, 0, 5)
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Width")), 1, 4)
        self.WCalc = Entrada_con_unidades(Length, retornar=False)
        self.WCalc.setReadOnly(True)
        lyt.addWidget(self.WCalc, 1, 5)
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Length")), 2, 4)
        self.LCalc = Entrada_con_unidades(Length, retornar=False)
        self.LCalc.setReadOnly(True)
        lyt.addWidget(self.LCalc, 2, 5)
        lyt_Calc.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding),
            12, 1, 1, 6)

        # Output tab
        self.addSalida(QtWidgets.QApplication.translate("pychemqt", "Filtered gas"))
        self.addSalida(
            QtWidgets.QApplication.translate("pychemqt", "Collected solids"))

        self.tipoCalculoCambiado(0)
        if equipment:
            self.setEquipment(equipment)

    def tipoCalculoCambiado(self, int):
        self.W.setReadOnly(int)
        self.W.setRetornar(not int)
        self.W.setResaltado(not int)
        self.H.setResaltado(not int)
        self.L.setReadOnly(int)
        self.L.setRetornar(not int)
        self.L.setResaltado(not int)
        self.rendimientoAdmisible.setReadOnly(not int)
        self.rendimientoAdmisible.setResaltado(int)
        self.velocidadAdmisible.setReadOnly(not int)
        self.velocidadAdmisible.setResaltado(False)
        self.changeParams("metodo", int)


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
    corriente = Corriente(T=300, P=101325, caudalMasico=1., fraccionMolar=[1.],
                          solido=solido)
    camara = GravityChamber(entrada=corriente, metodo=1, modelo=1, H=2.,
                            rendimientoAdmisible=0.95)
    dialogo = UI_equipment(camara)
    dialogo.show()
    sys.exit(app.exec_())
