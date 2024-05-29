#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


###############################################################################
# Grinder equipment dialog
###############################################################################


from functools import partial


from equipment.parents import UI_equip
from equipment.solids import Grinder
from lib import unidades
from tools import costIndex
from tools.qt import QtWidgets, tr
from UI.widgets import Entrada_con_unidades


BondIndex = {
    'Mineral de uranio': 17.93,
    'Escoria': 15.76,
    'Ferrocromo': 8.87,
    'Grafito': 45.03,
    'Magnesita': 16.8,
    'Mineral de plata': 17.3,
    'Molibdeno': 12.97,
    'Ferromanganeso': 7.77,
    'Arenisca': 11.53,
    'Arcilla': 7.1,
    'Mineral de níquel': 11.88,
    'Mineral de estaño': 10.81,
    'Mineral de titanio': 11.88,
    'Silicato sódico': 13.0,
    'Granito': 14.39,
    'Coque': 20.7,
    'Taconita': 14.87,
    'Hematita especular': 15.4,
    'Arena de silicato': 16.46,
    'Coque de petróleo': 73.8,
    'Gneiss': 20.13,
    'Carburo de silicio': 26.17,
    'Mineral de zinc': 12.42,
    'Granate': 12.37,
    'Caliza': 11.61,
    'Basalto': 20.41,
    'Carbón': 11.37,
    'Gabro': 18.45,
    'Dolomita': 11.31,
    'Coque de petróleo líquido': 38.6,
    'Mineral de plomo-zinc': 11.35,
    'Sal potásica': 8.23,
    'Andesita': 22.13,
    'Arcilla calcinada': 1.43,
    'Ilmenita': 13.11,
    'Mineral de hierro': 15.44,
    'Mica': 134.5,
    'Hematita': 12.68,
    'Fosfato fertilizante': 13.03,
    'Cemento en bruto': 10.57,
    'Bauxita': 9.45,
    'Mineral de plomo': 11.4,
    'Trapp': 21.1,
    'Cristal': 3.08,
    'Sienita': 14.9,
    'Coral': 10.16,
    'Roca fosfática': 10.13,
    'Caliza para cemento': 10.18,
    'Silicato': 13.53,
    'Aljez': 8.16,
    'Mineral de cromo': 9.6,
    'Feldespato': 11.67,
    'Mineral de cobre': 13.13,
    'Pizarra bituminosa': 18.1,
    'Cerámica': 15.53,
    'Pirita': 8.9,
    'Mineral de manganeso': 12.46,
    'Pirrotina': 9.57,
    'Cianita': 18.87,
    'Grava': 25.17,
    'Ferrosilicio': 12.83,
    'Sílex': 26.16,
    'Pizarra, mineral': 16.4,
    'Limanita': 8.45,
    'Barita': 6.24,
    'Esmeril': 58.18,
    'Escoria de hornos de hierro': 12.16,
    'Mineral de oro': 14.83,
    'Pumita': 11.93,
    'Rutilo': 12.12,
    'Espodumena': 13.7,
    'Fluorita': 9.76,
    'Clinker de cemento': 13.49,
    'Sinterizado': 8.77,
    'Galena': 10.19,
    'Magnetita': 10.21,
    'Cuarcita': 12.18,
    'Oolita': 11.33,
    'Pizarra, roca': 13.83,
    'Mineral de potasio': 8.88,
    'Diorita': 19.4,
    'Cuarzo': 12.77}


class UI_equipment(UI_equip):
    """Diálogo de definición de molinos trituradores de sólidos"""
    Equipment = Grinder()

    def __init__(self, equipment=None, parent=None):
        """equipment: Initial equipment instance to model"""
        super().__init__(Grinder, entrada=False, salida=False, parent=parent)

        # Calculate tab
        lyt = QtWidgets.QGridLayout(self.tabCalculo)
        lyt.addWidget(QtWidgets.QLabel(
            tr("pychemqt", "Bond work index")), 1, 0, 1, 1)
        self.Material = QtWidgets.QComboBox()
        self.Material.addItem(tr("pychemqt", "User defined"))
        for key in sorted(BondIndex.keys()):
            self.Material.addItem(key)
        self.Material.currentIndexChanged.connect(self.cambiarBondWordIndex)
        lyt.addWidget(self.Material, 1, 1, 1, 1)
        self.BondIndex = Entrada_con_unidades(float)
        self.BondIndex.valueChanged.connect(
            partial(self.changeParams, "BondIndex"))
        lyt.addWidget(self.BondIndex, 1, 2, 1, 1)
        lyt.addWidget(QtWidgets.QLabel(
            tr("pychemqt", "Exponent")), 2, 0, 1, 1)
        self.exponent = Entrada_con_unidades(float)
        self.exponent.valueChanged.connect(
            partial(self.changeParams, "exponent"))
        lyt.addWidget(self.exponent, 2, 1, 1, 2)
        lyt.addWidget(QtWidgets.QLabel(
            tr("pychemqt", "D80")), 3, 0, 1, 1)
        self.D80 = Entrada_con_unidades(unidades.Length, "ParticleDiameter")
        self.D80.valueChanged.connect(
            partial(self.changeParams, "D80"))
        lyt.addWidget(self.D80, 3, 1, 1, 2)

        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 4, 0, 1, 5)

        group = QtWidgets.QGroupBox(
            tr("pychemqt", "Results"))
        layout = QtWidgets.QGridLayout(group)
        layout.addWidget(QtWidgets.QLabel(
            tr("pychemqt", "Power")), 1, 0)
        self.power = Entrada_con_unidades(
            unidades.Power, retornar=False, readOnly=True)
        layout.addWidget(self.power, 1, 1)
        layout.addWidget(QtWidgets.QLabel(
            tr("pychemqt", "Solid Mass Flow")), 2, 0)
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
            QtWidgets.QLabel(tr("pychemqt", "Type:")), 1, 1, 1, 1)
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
        groupBox_Costos = QtWidgets.QGroupBox(tr("pychemqt", "Stimated costs"))
        lyt.addWidget(groupBox_Costos, 7, 1, 1, 6)
        lytgroup = QtWidgets.QGridLayout(groupBox_Costos)
        lytgroup.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Purchase cost")), 0, 1)
        self.C_adq = Entrada_con_unidades(
            unidades.Currency, retornar=False, readOnly=True)
        lytgroup.addWidget(self.C_adq, 0, 2)
        lytgroup.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Installed cost")), 1, 1)
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
