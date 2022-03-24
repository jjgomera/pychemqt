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
# Scrubber equipment definition dialog
###############################################################################


from functools import partial

from PyQt5 import QtWidgets

from lib.unidades import Length, DeltaP
from UI import UI_corriente
from UI.widgets import Entrada_con_unidades
from equipment.gas_solid_liquid import Scrubber
from equipment.parents import UI_equip


class UI_equipment(UI_equip):
    Equipment = Scrubber()
    """Scrubber equipment definition dialog"""
    def __init__(self, equipment=None, parent=None):
        """
        equipment: Initial equipment instance to model
        """
        super().__init__(Scrubber, parent=parent)

        # Input tab
        self.entradaGas = UI_corriente.Ui_corriente(psychro=True)
        self.entradaGas.Changed.connect(
            partial(self.changeParams, "entradaGas"))
        self.Entrada.addTab(
            self.entradaGas,
            QtWidgets.QApplication.translate("equipment", "Gas"))
        self.entradaLiquido = UI_corriente.Ui_corriente()
        self.entradaLiquido.Changed.connect(
            partial(self.changeParams, "entradaLiquido"))
        self.Entrada.addTab(
            self.entradaLiquido,
            QtWidgets.QApplication.translate("pychemqt", "Liquid"))

        # Calculate tab
        gridLayout_Calculo = QtWidgets.QGridLayout(self.tabCalculo)
        gridLayout_Calculo.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Mode")), 1, 1)
        self.tipo_calculo = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_TIPO:
            self.tipo_calculo.addItem(txt)
        self.tipo_calculo.currentIndexChanged.connect(
            self.on_tipoCalculo_currentIndexChanged)
        gridLayout_Calculo.addWidget(self.tipo_calculo, 1, 2, 1, 5)
        gridLayout_Calculo.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Method")), 2, 1)
        self.modelo_rendimiento = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_MODEL:
            self.modelo_rendimiento.addItem(txt)
        self.modelo_rendimiento.currentIndexChanged.connect(
            self.on_modeloRendimiento_currentIndexChanged)
        gridLayout_Calculo.addWidget(self.modelo_rendimiento, 2, 2, 1, 5)
        gridLayout_Calculo.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "ΔP method")), 3, 1)
        self.modelo_DeltaP = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_MODEL_DELTAP:
            self.modelo_DeltaP.addItem(txt)
        self.modelo_DeltaP.currentIndexChanged.connect(
            self.on_modeloDeltaP_currentIndexChanged)
        gridLayout_Calculo.addWidget(self.modelo_DeltaP, 3, 2, 1, 5)

        gridLayout_Calculo.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Fixed,
            QtWidgets.QSizePolicy.Fixed), 4, 1, 1, 6)
        gridLayout_Calculo.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Diameter")), 5, 1)
        self.diametro = Entrada_con_unidades(Length)
        self.diametro.valueChanged.connect(
            partial(self.changeParams, "diametro"))
        gridLayout_Calculo.addWidget(self.diametro, 5, 2)
        gridLayout_Calculo.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Efficiency")), 6, 1)
        self.rendimientoAdmisible = Entrada_con_unidades(float, spinbox=True)
        self.rendimientoAdmisible.valueChanged.connect(
            partial(self.changeParams, "rendimientoAdmisible"))
        gridLayout_Calculo.addWidget(self.rendimientoAdmisible, 6, 2)

        self.groupJohnstone = QtWidgets.QWidget()
        gridLayout_Calculo.addWidget(self.groupJohnstone, 7, 1, 1, 2)
        JohnstoneLayout = QtWidgets.QHBoxLayout(self.groupJohnstone)
        JohnstoneLayout.setSpacing(0)
        JohnstoneLayout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Ventury Constant")))
        self.k = Entrada_con_unidades(float, spinbox=True)
        self.k.valueChanged.connect(partial(self.changeParams, "k"))
        JohnstoneLayout.addWidget(self.k)

        self.groupCalvert = QtWidgets.QWidget()
        gridLayout_Calculo.addWidget(self.groupCalvert, 7, 1, 1, 2)
        CalvertLayout = QtWidgets.QHBoxLayout(self.groupCalvert)
        CalvertLayout.setSpacing(0)
        CalvertLayout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "f parameter")))
        self.f = Entrada_con_unidades(float, spinbox=True)
        self.f.valueChanged.connect(partial(self.changeParams, "f"))
        CalvertLayout.addWidget(self.f)

        self.groupLt = QtWidgets.QWidget()
        gridLayout_Calculo.addWidget(self.groupLt, 5, 4, 1, 2)
        LtLayout = QtWidgets.QHBoxLayout(self.groupLt)
        LtLayout.setSpacing(0)
        LtLayout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Length throat")))
        self.Lt = Entrada_con_unidades(Length)
        self.Lt.valueChanged.connect(partial(self.changeParams, "Lt"))
        LtLayout.addWidget(self.Lt)

        gridLayout_Calculo.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 8, 1, 1, 6)
        self.groupBox_Calculo = QtWidgets.QGroupBox(
            QtWidgets.QApplication.translate("pychemqt", "Results"))
        gridLayout_Calculo.addWidget(self.groupBox_Calculo, 9, 1, 1, 5)
        gridLayout_1 = QtWidgets.QGridLayout(self.groupBox_Calculo)
        gridLayout_1.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Efficiency")), 1, 1)
        self.rendimiento = Entrada_con_unidades(
            float, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.rendimiento, 1, 2)
        gridLayout_1.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "DeltaP")), 2, 1)
        self.deltaP = Entrada_con_unidades(
            DeltaP, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.deltaP, 2, 2)

        gridLayout_Calculo.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 11, 1, 1, 6)

        # Output Tab
        self.SalidaGas = UI_corriente.Ui_corriente(readOnly=True, psychro=True)
        self.Salida.addTab(
            self.SalidaGas,
            QtWidgets.QApplication.translate("pychemqt", "Clean Gas"))
        self.SalidaLiquido = UI_corriente.Ui_corriente(readOnly=True)
        self.Salida.addTab(
            self.SalidaLiquido,
            QtWidgets.QApplication.translate("pychemqt", "Liquid"))

        self.on_tipoCalculo_currentIndexChanged(0)
        self.on_modeloRendimiento_currentIndexChanged(0)
        self.on_modeloDeltaP_currentIndexChanged(0)
        if equipment:
            self.setEquipment(equipment)

    def on_tipoCalculo_currentIndexChanged(self, modelo):
        self.rendimientoAdmisible.setEnabled(modelo)
        self.rendimientoAdmisible.setReadOnly(not modelo)
        self.diametro.setEnabled(not modelo)
        self.diametro.setReadOnly(modelo)
        self.changeParams("tipo_calculo", modelo)

    def on_modeloRendimiento_currentIndexChanged(self, modelo):
        self.groupJohnstone.setVisible(False)
        self.groupCalvert.setVisible(False)
        if modelo == 0:
            self.groupJohnstone.setVisible(True)
        elif modelo == 1:
            self.groupCalvert.setVisible(True)
        self.changeParams("modelo_rendimiento", modelo)

    def on_modeloDeltaP_currentIndexChanged(self, modelo):
        self.groupLt.setVisible(False)
        if modelo in (3, 4):
            self.groupLt.setVisible(True)
        self.changeParams("modelo_DeltaP", modelo)

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    from lib.corriente import Corriente, Solid

    diametros = [17.5e-6, 22.4e-6, 26.2e-6, 31.8e-6, 37e-6, 42.4e-6, 48e-6,
                 54e-6, 60e-6, 69e-6, 81.3e-6, 96.5e-6, 109e-6, 127e-6]
    fracciones = [0.02, 0.03, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                  0.05, 0.03, 0.02]
    solido = Solid(T=300, caudalSolido=[1/3600.],
                   distribucion_diametro=diametros,
                   distribucion_fraccion=fracciones, solids=[638])
    kw = {"fraccionMolar": [1.], "MEoS": True}
    aire = Corriente(T=350, P=101325, caudalMasico=0.01, ids=[475],
                     solido=solido, **kw)
    agua = Corriente(T=300, P=101325, caudalMasico=0.1, ids=[62], **kw)
    scrubber = Scrubber(entradaGas=aire, entradaLiquido=agua, diametro=0.25,
                        modelo_rendimiento=0, modelo_DeltaP=1, k=1000)
    dialogo = UI_equipment(scrubber)

#    dialogo = UI_equipment()

    dialogo.show()
    sys.exit(app.exec_())
