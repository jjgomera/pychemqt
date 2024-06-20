#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2023, Juan José Gómez Romera <jjgomera@gmail.com>

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
# Cyclone equipment dialog
###############################################################################


from functools import partial
import os

from lib.unidades import Length, Pressure, DeltaP, Speed, VolFlow, Currency
from equipment.gas_solid import Ciclon
from equipment.parents import UI_equip
from tools.costIndex import CostData
from tools.qt import QtGui, QtWidgets
from UI.widgets import Entrada_con_unidades


class UI_equipment(UI_equip):
    """Cyclone equipment edition dialog"""
    Equipment = Ciclon()

    def __init__(self, equipment=None, parent=None):
        """
        equipment: Initial equipment instance to model
        """
        super().__init__(Ciclon, entrada=False, parent=parent)

        # Calculate tab
        lyt_Calc = QtWidgets.QGridLayout(self.tabCalculo)
        lyt_Calc.addWidget(QtWidgets.QLabel(self.tr("Mode")), 0, 1, 1, 2)
        self.tipo_calculo = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_TIPO:
            self.tipo_calculo.addItem(txt)
        self.tipo_calculo.currentIndexChanged.connect(
            self.on_tipoCalculo_currentIndexChanged)
        lyt_Calc.addWidget(self.tipo_calculo, 0, 3, 1, 4)
        lyt_Calc.addWidget(QtWidgets.QLabel(self.tr("Method")), 1, 1, 1, 2)
        self.modelo_rendimiento = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_MODEL:
            self.modelo_rendimiento.addItem(txt)
        self.modelo_rendimiento.currentIndexChanged.connect(
            partial(self.changeParams, "modelo_rendimiento"))
        lyt_Calc.addWidget(self.modelo_rendimiento, 1, 3, 1, 4)
        lyt_Calc.addWidget(QtWidgets.QLabel(self.tr("ΔP method", None)),
            2, 1, 1, 2)
        self.modelo_DeltaP = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_MODEL_DELTAP:
            self.modelo_DeltaP.addItem(txt)
        self.modelo_DeltaP.currentIndexChanged.connect(
            partial(self.changeParams, "modelo_DeltaP"))
        lyt_Calc.addWidget(self.modelo_DeltaP, 2, 3, 1, 4)
        lyt_Calc.addWidget(QtWidgets.QLabel(
            self.tr("Design model")), 3, 1, 1, 2)
        self.modelo_ciclon = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_MODEL_CICLON:
            self.modelo_ciclon.addItem(txt)
        self.modelo_ciclon.currentIndexChanged.connect(
            self.modeloEficiencia_Changed)
        lyt_Calc.addWidget(self.modelo_ciclon, 3, 3, 1, 4)
        lyt_Calc.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 4, 0, 1, 5)

        lyt_Calc.addWidget(QtWidgets.QLabel(self.tr("Diameter")), 5, 1)
        self.Dc = Entrada_con_unidades(Length)
        self.Dc.valueChanged.connect(partial(self.changeParams, "Dc"))
        lyt_Calc.addWidget(self.Dc, 5, 2, 1, 2)
        lyt_Calc.addWidget(QtWidgets.QLabel(self.tr("Efficiency")), 6, 1)
        self.rendimientoAdmisible = Entrada_con_unidades(float, spinbox=True)
        self.rendimientoAdmisible.valueChanged.connect(
            partial(self.changeParams, "rendimientoAdmisible"))
        lyt_Calc.addWidget(self.rendimientoAdmisible, 6, 2, 1, 2)
        lyt_Calc.addWidget(QtWidgets.QLabel(self.tr("Allowable ΔP", None)),
            7, 1)
        self.DeltaPAdmisible = Entrada_con_unidades(Pressure)
        self.DeltaPAdmisible.valueChanged.connect(
            partial(self.changeParams, "DeltaPAdmisible"))
        lyt_Calc.addWidget(self.DeltaPAdmisible, 7, 2, 1, 2)
        lyt_Calc.addWidget(QtWidgets.QLabel(self.tr("Allowable speed")), 8, 1)
        self.velocidadAdmisible = Entrada_con_unidades(Speed)
        self.velocidadAdmisible.valueChanged.connect(
            partial(self.changeParams, "velocidadAdmisible"))
        lyt_Calc.addWidget(self.velocidadAdmisible, 8, 2, 1, 2)
        lyt_Calc.addWidget(QtWidgets.QLabel(self.tr("No. of ciclones")), 9, 1)
        self.num_ciclones = Entrada_con_unidades(int, spinbox=True, step=1,
                                                 decimales=0, min=1)
        self.num_ciclones.valueChanged.connect(
            partial(self.changeParams, "num_ciclones"))
        lyt_Calc.addWidget(self.num_ciclones, 9, 2, 1, 2)
        lyt_Calc.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 10, 0, 1, 4)

        group = QtWidgets.QGroupBox(self.tr("Results"))
        lyt_Calc.addWidget(group, 11, 1, 1, 3)
        lyt = QtWidgets.QGridLayout(group)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Efficiency")), 1, 1)
        self.rendimiento = Entrada_con_unidades(float, readOnly=True)
        lyt.addWidget(self.rendimiento, 1, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Pressure drop:")), 2, 1)
        self.deltaP = Entrada_con_unidades(DeltaP, readOnly=True)
        lyt.addWidget(self.deltaP, 2, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Admission speed")), 3, 1)
        self.V = Entrada_con_unidades(Speed, readOnly=True)
        lyt.addWidget(self.V, 3, 2)

        lyt_Calc.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 5, 4, 7, 1)
        group2 = QtWidgets.QGroupBox(self.tr("Geometry"))
        lyt_Calc.addWidget(group2, 5, 5, 7, 1)
        lyt = QtWidgets.QGridLayout(group2)
        lyt.addWidget(QtWidgets.QLabel("D<sub>C</sub>"), 1, 1)
        self.Dcc = Entrada_con_unidades(Length, boton=False)
        lyt.addWidget(self.Dcc, 1, 2)
        lyt.addWidget(QtWidgets.QLabel("B<sub>C</sub>"), 2, 1)
        self.Bc = Entrada_con_unidades(Length, boton=False)
        lyt.addWidget(self.Bc, 2, 2)
        lyt.addWidget(QtWidgets.QLabel("H<sub>C</sub>"), 3, 1)
        self.Hc = Entrada_con_unidades(Length, boton=False)
        lyt.addWidget(self.Hc, 3, 2)
        lyt.addWidget(QtWidgets.QLabel("J<sub>C</sub>"), 4, 1)
        self.Jc = Entrada_con_unidades(Length, boton=False)
        lyt.addWidget(self.Jc, 4, 2)
        lyt.addWidget(QtWidgets.QLabel("L<sub>C</sub>"), 5, 1)
        self.Lc = Entrada_con_unidades(Length, boton=False)
        lyt.addWidget(self.Lc, 5, 2)
        lyt.addWidget(QtWidgets.QLabel("Z<sub>C</sub>"), 6, 1)
        self.Zc = Entrada_con_unidades(Length, boton=False)
        lyt.addWidget(self.Zc, 6, 2)
        lyt.addWidget(QtWidgets.QLabel("D<sub>e</sub>"), 7, 1)
        self.De = Entrada_con_unidades(Length, boton=False)
        lyt.addWidget(self.De, 7, 2)
        lyt.addWidget(QtWidgets.QLabel("S<sub>C</sub>"), 8, 1)
        self.Sc = Entrada_con_unidades(Length, boton=False)
        lyt.addWidget(self.Sc, 8, 2)
        lyt.addWidget(QtWidgets.QLabel("N"), 9, 1)
        self.NCalc = Entrada_con_unidades(int, readOnly=True)
        lyt.addWidget(self.NCalc, 9, 2)
        lyt_Calc.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 0, 7, 10, 1)
        lyt_Calc.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 5, 6, 10, 1)

        image = QtWidgets.QLabel()
        path = os.environ["pychemqt"]+"/images/equip/ciclon.gif"
        image.setPixmap(QtGui.QPixmap(path))
        image.setScaledContents(True)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding)
        image.setSizePolicy(sizePolicy)
        lyt_Calc.addWidget(image, 0, 8, 12, 1)

        # Cost tab
        lyt_Cost = QtWidgets.QGridLayout(self.tabCostos)
        lyt_Cost.addWidget(QtWidgets.QLabel(self.tr("Model")), 1, 1)
        self.tipo_costo = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_COST:
            self.tipo_costo.addItem(txt)
        self.tipo_costo.currentIndexChanged.connect(
            partial(self.changeParamsCoste, "tipo_costo"))
        lyt_Cost.addWidget(self.tipo_costo, 1, 2)
        lyt_Cost.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed), 1, 3)

        self.Costos = CostData(self.Equipment)
        self.Costos.valueChanged.connect(self.changeParamsCoste)
        lyt_Cost.addWidget(self.Costos, 3, 1, 1, 3)
        lyt_Cost.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 4, 1, 1, 3)

        group = QtWidgets.QGroupBox(self.tr("Stimated Costs"))
        lyt_Cost.addWidget(group, 5, 1, 1, 3)
        lyt = QtWidgets.QGridLayout(group)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Number")), 0, 0)
        self.num_ciclonesCoste = Entrada_con_unidades(int, readOnly=True)
        lyt.addWidget(self.num_ciclonesCoste, 0, 1, 1, 3)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Volumetric Flow")), 1, 0)
        self.Q = Entrada_con_unidades(VolFlow, "QGas", retornar=False)
        self.Q.setReadOnly(True)
        lyt.addWidget(self.Q, 1, 1)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Purchase cost")), 0, 3)
        self.C_adq = Entrada_con_unidades(Currency, retornar=False)
        self.C_adq.setReadOnly(True)
        lyt.addWidget(self.C_adq, 0, 4)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Installed cost")), 1, 3)
        self.C_inst = Entrada_con_unidades(Currency, retornar=False)
        self.C_inst.setReadOnly(True)
        lyt.addWidget(self.C_inst, 1, 4)

        # Output tab
        self.addSalida(self.tr("Filtered gas"))
        self.addSalida(self.tr("Collected solids"))

        self.on_tipoCalculo_currentIndexChanged(0)
        self.modeloEficiencia_Changed(0)
        if equipment:
            self.setEquipment(equipment)

    def on_tipoCalculo_currentIndexChanged(self, int):
        """Habilita o desabilita los datos requeridos para el cálculo"""
        if int and self.modelo_ciclon.count() == 9:
            self.modelo_ciclon.removeItem(8)
        elif not int and self.modelo_ciclon.count() == 8:
            self.modelo_ciclon.addItem(self.tr("Custom"))

        self.Dc.setReadOnly(int)
        self.Dc.setResaltado(not int)
        self.num_ciclones.setReadOnly(int)
        self.num_ciclones.setResaltado(not int)
        self.rendimientoAdmisible.setReadOnly(not int)
        self.rendimientoAdmisible.setResaltado(int)
        self.velocidadAdmisible.setReadOnly(not int)
        self.velocidadAdmisible.setResaltado(False)
        self.DeltaPAdmisible.setReadOnly(not int)
        self.DeltaPAdmisible.setResaltado(False)
        self.changeParams("tipo_calculo", int)

    def modeloEficiencia_Changed(self, modelo):
        if modelo == 8:
            # Customized model, let user edit dimensions values
            self.Hc.setReadOnly(False)
            self.Hc.setResaltado(True)
            self.Bc.setReadOnly(False)
            self.Bc.setResaltado(True)
            self.Jc.setReadOnly(False)
            self.Jc.setResaltado(False)
            self.Lc.setReadOnly(False)
            self.Lc.setResaltado(False)
            self.Zc.setReadOnly(False)
            self.Zc.setResaltado(False)
            self.De.setReadOnly(False)
            self.De.setResaltado(True)
            self.Sc.setReadOnly(False)
            self.Sc.setResaltado(False)
        else:
            self.Hc.setReadOnly(True)
            self.Hc.setResaltado(False)
            self.Bc.setReadOnly(True)
            self.Bc.setResaltado(False)
            self.Jc.setReadOnly(True)
            self.Jc.setResaltado(False)
            self.Lc.setReadOnly(True)
            self.Lc.setResaltado(False)
            self.Zc.setReadOnly(True)
            self.Zc.setResaltado(False)
            self.De.setReadOnly(True)
            self.De.setResaltado(False)
            self.Sc.setReadOnly(True)
            self.Sc.setResaltado(False)
        self.Dcc.setReadOnly(True)
        self.Dcc.setResaltado(False)
        self.changeParams("modelo_ciclon", modelo)


if __name__ == "__main__":
    import sys
    from lib.corriente import Corriente, Solid
    # app = QtWidgets.QApplication(sys.argv)
    # diametros = [17.5e-6, 22.4e-6, 26.2e-6, 31.8e-6, 37e-6, 42.4e-6, 48e-6,
                 # 54e-6, 60e-6, 69e-6, 81.3e-6, 96.5e-6, 109e-6, 127e-6]
    # fracciones = [0.02, 0.03, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                  # 0.05, 0.03, 0.02]
    # solido = Solid(caudalSolido=[0.1], distribucion_diametro=diametros,
                   # distribucion_fraccion=fracciones)
    # corriente = Corriente(T=300, P=101325, ids=[475], solids=[638],
                          # caudalMasico=1., fraccionMolar=[1.], solido=solido)

    # ciclon = Ciclon(entrada=corriente, tipo_calculo=1, velocidadAdmisible=5)
    # print(ciclon.status)
    # ciclon.Pdf()
    dm = [17.5e-6, 22.4e-6, 26.2e-6, 31.8e-6, 37e-6, 42.4e-6, 48e-6, 54e-6,
          60e-6, 69e-6, 81.3e-6, 96.5e-6, 109e-6, 127e-6]
    fracciones = [0.02, 0.03, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                  0.05, 0.03, 0.02]
    sol = Solid(caudalSolido=[0.1], distribucion_diametro=dm,
                distribucion_fraccion=fracciones)
    kw = {"ids": [475], "solids": [638], "fraccionMolar": [1.], "MEoS": True}
    entrada = Corriente(T=300, P=1e5, caudalMasico=1, solido=sol, **kw)
    ciclon = Ciclon(entrada=entrada, tipo_calculo=1, rendimientoAdmisible=0.95,
                    velocidadAdmisible=5)
    print("%0.2f %0.2f" % (ciclon.C_instTotal, ciclon.C_adqTotal))
    # 7597.86 5427.04


    # Ciclon = UI_equipment(ciclon)
    # Ciclon.show()
    # sys.exit(app.exec())
