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
# Baghouse equipment dialog
###############################################################################


from functools import partial

from tools.qt import QtWidgets, tr

from numpy import any

from lib.unidades import Time, Pressure, Length, Area, Speed
from equipment.gas_solid import Baghouse
from equipment.parents import UI_equip
from UI.widgets import Entrada_con_unidades, Tabla


class UI_equipment(UI_equip):
    """Baghouse equipment edition dialog"""
    Equipment = Baghouse()

    def __init__(self, equipment=None, parent=None):
        """
        equipment: Initial equipment instance to model
        """
        super().__init__(Baghouse, entrada=False, parent=parent)

        # Efficiency tab
        title = [tr("pychemqt", "Diameter") +
                 ", " + Length.text("ParticleDiameter"),
                 tr("pychemqt", "Efficiency")]
        self.efic = Tabla(2, horizontalHeader=title, filas=1, stretch=False)
        self.efic.setColumnReadOnly(0, True)
        self.efic.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.efic.editingFinished.connect(self.cambiarRendimientos)
        self.tabWidget.insertTab(
            1, self.efic,
            tr("pychemqt", "Efficiencies"))

        # Calculate tab
        lyt_Calc = QtWidgets.QGridLayout(self.tabCalculo)
        lyt_Calc.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Mode")), 1, 1)
        self.metodo = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_TIPO:
            self.metodo.addItem(txt)
        self.metodo.currentIndexChanged.connect(self.tipoCalculoCambiado)
        lyt_Calc.addWidget(self.metodo, 1, 2, 1, 3)
        lyt_Calc.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed, QtWidgets.QSizePolicy.Policy.Fixed),
            2, 1, 1, 6)
        lyt_Calc.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "No cells")), 3, 1)
        self.num_filtros = Entrada_con_unidades(
            int, spinbox=True, step=1, min=1, width=50, resaltado=True,
            start=1)
        self.num_filtros.valueChanged.connect(
            partial(self.changeParams, "num_filtros"))
        lyt_Calc.addWidget(self.num_filtros, 3, 2)
        lyt_Calc.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Time")), 4, 1)
        self.tiempo = Entrada_con_unidades(Time, resaltado=True)
        self.tiempo.valueChanged.connect(partial(self.changeParams, "tiempo"))
        lyt_Calc.addWidget(self.tiempo, 4, 2)
        lyt_Calc.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Pressure drop")), 5, 1)
        self.deltaP = Entrada_con_unidades(Pressure, retornar=False)
        self.deltaP.setReadOnly(True)
        self.deltaP.valueChanged.connect(partial(self.changeParams, "deltaP"))
        lyt_Calc.addWidget(self.deltaP, 5, 2)
        lyt_Calc.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed, QtWidgets.QSizePolicy.Policy.Fixed),
            6, 1, 1, 6)
        lyt_Calc.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Bags per cell")), 7, 1)
        self.membranasFiltro = Entrada_con_unidades(int, spinbox=True, step=1,
                                                    min=1)
        self.membranasFiltro.valueChanged.connect(
            partial(self.changeParams, "membranasFiltro"))
        lyt_Calc.addWidget(self.membranasFiltro, 7, 2)
        lyt_Calc.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Bag diameter")), 8, 1)
        self.diametroMembrana = Entrada_con_unidades(Length)
        self.diametroMembrana.valueChanged.connect(
            partial(self.changeParams, "diametroMembrana"))
        lyt_Calc.addWidget(self.diametroMembrana, 8, 2)
        lyt_Calc.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Area per bag")), 9, 1)
        self.areaMembrana = Entrada_con_unidades(Area)
        self.areaMembrana.valueChanged.connect(
            partial(self.changeParams, "areaMembrana"))
        lyt_Calc.addWidget(self.areaMembrana, 9, 2)
        lyt_Calc.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Cloth resistence")), 7, 4)
        self.resistenciaFiltro = Entrada_con_unidades(float)
        self.resistenciaFiltro.valueChanged.connect(
            partial(self.changeParams, "resistenciaFiltro"))
        lyt_Calc.addWidget(self.resistenciaFiltro, 7, 5)
        lyt_Calc.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Cake resistence")), 8, 4)
        self.resistenciaTorta = Entrada_con_unidades(float)
        self.resistenciaTorta.valueChanged.connect(
            partial(self.changeParams, "resistenciaTorta"))
        lyt_Calc.addWidget(self.resistenciaTorta, 8, 5)
        lyt_Calc.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Cells cleaned")), 9, 4)
        self.limpieza = Entrada_con_unidades(int, spinbox=True, step=1, min=0)
        self.limpieza.valueChanged.connect(
            partial(self.changeParams, "limpieza"))
        lyt_Calc.addWidget(self.limpieza, 9, 5)
        lyt_Calc.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 10, 1, 1, 6)

        groupbox = QtWidgets.QGroupBox(
            tr("pychemqt", "Results"))
        lyt_Calc.addWidget(groupbox, 11, 1, 1, 5)
        lyt = QtWidgets.QGridLayout(groupbox)

        lyt.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "No cells")), 1, 1)
        self.num_filtrosCalc = Entrada_con_unidades(int, readOnly=True)
        lyt.addWidget(self.num_filtrosCalc, 1, 2)
        lyt.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Time")), 2, 1)
        self.tiempoCalc = Entrada_con_unidades(Time, readOnly=True)
        lyt.addWidget(self.tiempoCalc, 2, 2)
        lyt.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Pressure drop")), 3, 1)
        self.deltaPCalc = Entrada_con_unidades(Pressure, readOnly=True)
        lyt.addWidget(self.deltaPCalc, 3, 2)

        lyt.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Gas velocity")), 1, 4)
        self.Vgas = Entrada_con_unidades(Speed, retornar=False, readOnly=True)
        lyt.addWidget(self.Vgas, 1, 5)
        lyt.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Efficiency")), 2, 4)
        self.rendimiento = Entrada_con_unidades(float, readOnly=True)
        lyt.addWidget(self.rendimiento, 2, 5)
        lyt.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Area")), 3, 4)
        self.floorArea = Entrada_con_unidades(Area, readOnly=True)
        lyt.addWidget(self.floorArea, 3, 5)
        lyt_Calc.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 12, 1, 1, 6)

        # Output tab
        self.addSalida(
            tr("pychemqt", "Filtered gas"))
        self.addSalida(
            tr("pychemqt", "Collected solids"))

        if equipment:
            self.setEquipment(equipment)

    def cambiarRendimientos(self):
        self.changeParams("rendimientos", self.efic.getColumn(1))

    def tipoCalculoCambiado(self, tipo_calculo):
        if tipo_calculo == 0:
            self.num_filtros.setReadOnly(False)
            self.num_filtros.setResaltado(True)
            self.tiempo.setReadOnly(False)
            self.tiempo.setResaltado(True)
            self.deltaP.setReadOnly(True)
            self.deltaP.setResaltado(False)
        elif tipo_calculo == 1:
            self.num_filtros.setReadOnly(False)
            self.num_filtros.setResaltado(True)
            self.tiempo.setReadOnly(True)
            self.tiempo.setResaltado(False)
            self.deltaP.setReadOnly(False)
            self.deltaP.setResaltado(True)
        else:
            self.num_filtros.setReadOnly(True)
            self.num_filtros.setResaltado(False)
            self.tiempo.setReadOnly(False)
            self.tiempo.setResaltado(True)
            self.deltaP.setReadOnly(False)
            self.deltaP.setResaltado(True)
        self.changeParams("metodo", tipo_calculo)

    def rellenarInput(self):
        UI_equip.rellenarInput(self)
        if self.Equipment.kwargs["entrada"].solido:
            diametros = []
            for d in self.Equipment.kwargs["entrada"].solido.diametros:
                diametros.append(d.config("ParticleDiameter"))
            self.efic.setColumn(0, diametros)
        if any(self.Equipment.kwargs["rendimientos"]):
            self.efic.setColumn(1, self.Equipment.kwargs["rendimientos"])


if __name__ == "__main__":
    import sys
    from lib.corriente import Corriente, Solid
    app = QtWidgets.QApplication(sys.argv)
    diametros = [17.5e-6, 22.4e-6, 26.2e-6, 31.8e-6, 37e-6, 42.4e-6, 48e-6,
                 54e-6, 60e-6, 69e-6, 81.3e-6, 96.5e-6, 109e-6, 127e-6]
    fracciones = [0.02, 0.03, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                  0.05, 0.03, 0.02]
    # solido = Solid(caudalSolido=[0.1], distribucion_diametro=diametros,
                   # distribucion_fraccion=fracciones)
    # corriente = Corriente(T=300, P=101325, caudalMasico=1.,
                          # fraccionMolar=[1.], solido=solido)
    # filtro = Baghouse(entrada=corriente, metodo=0, num_filtros=4, tiempo=3600)
    dialogo = UI_equipment()
    dialogo.show()
    sys.exit(app.exec())
