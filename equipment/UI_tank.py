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
# Vessel equipment dialog
###############################################################################


from functools import partial

from tools.qt import QtCore, QtWidgets

from equipment.tank import Tank
from UI import UI_corriente
from equipment import parents
from equipment.widget import helical
from lib.unidades import Currency, Density, Mass, Length, Volume
from tools.costIndex import CostData
from UI.widgets import Entrada_con_unidades


class UI_equipment(parents.UI_equip):
    """Tank dialog"""
    Equipment = Tank()

    def __init__(self, equipment=None, parent=None):
        """
        equipment: Initial equipment instance to model
        """
        super().__init__(Tank, entrada=False, salida=False, parent=parent)

        #Calculate tab
        lytCalc = QtWidgets.QGridLayout(self.tabCalculo)
        lytCalc.addWidget(QtWidgets.QLabel(self.tr("Internal diameter")), 1, 1)
        self.Di = Entrada_con_unidades(Length)
        self.Di.valueChanged.connect(partial(self.changeParams, "Di"))
        lytCalc.addWidget(self.Di, 1, 2)
        lytCalc.addWidget(QtWidgets.QLabel(self.tr("Length")), 2, 1)
        self.L = Entrada_con_unidades(Length)
        self.L.valueChanged.connect(partial(self.changeParams, "L"))
        lytCalc.addWidget(self.L, 2, 2)
        lytCalc.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 3, 3)

        # Helical Coil definition tab
        self.helicalCoil = helical.UI_Helical()
        self.helicalCoil.toggled.connect(
            partial(self.changeParams, "hasHelicalCoil"))
        self.helicalCoil.valueChanged.connect(
            partial(self.changeParams, "helicalCoil"))
        self.tabWidget.insertTab(2, self.helicalCoil, self.tr("Helical Coil"))

        # Cost tab
        lytCost = QtWidgets.QGridLayout(self.tabCostos)
        lytCost.addWidget(QtWidgets.QLabel(self.tr("Material")), 1, 1)
        self.material = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_MATERIAL:
            self.material.addItem(txt)
        self.material.currentIndexChanged.connect(
            partial(self.changeParamsCoste, "material"))
        lytCost.addWidget(self.material, 1, 2, 1, 4)
        lytCost.addWidget(QtWidgets.QLabel(self.tr("Density")), 2, 4)
        self.Densidad = Entrada_con_unidades(Density, "DenLiq")
        lytCost.addWidget(self.Densidad, 2, 5)
        lytCost.addWidget(QtWidgets.QLabel(self.tr("Diameter")), 2, 1)
        self.Diametro = Entrada_con_unidades(Length)
        lytCost.addWidget(self.Diametro, 2, 2)
        lytCost.addWidget(QtWidgets.QLabel(self.tr("Length")), 3, 1)
        self.Longitud = Entrada_con_unidades(Length)
        lytCost.addWidget(self.Longitud, 3, 2)
        lytCost.addWidget(QtWidgets.QLabel(self.tr("Thickness")), 4, 1)
        self.Espesor = Entrada_con_unidades(Length, "Thickness")
        lytCost.addWidget(self.Espesor, 4, 2)
        lytCost.addWidget(QtWidgets.QLabel(self.tr("Head")), 5, 1)
        self.Cabeza = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_HEAD:
            self.Cabeza.addItem(txt)
        self.Cabeza.currentIndexChanged.connect(
            partial(self.changeParamsCoste, "head"))
        lytCost.addWidget(self.Cabeza, 5, 2)
        lytCost.addWidget(QtWidgets.QLabel(self.tr("Thickness (Head)")), 6, 1)
        self.EspesorCabeza = Entrada_con_unidades(Length, "Thickness")
        lytCost.addWidget(self.EspesorCabeza, 6, 2)
        lytCost.addWidget(QtWidgets.QLabel(
            self.tr("Longitud reborde recto")), 7, 1)
        self.LongitudReborde = Entrada_con_unidades(Length)
        lytCost.addWidget(self.LongitudReborde, 7, 2)
        lytCost.addWidget(QtWidgets.QLabel(self.tr("Volume")), 6, 4)
        self.Volumen = Entrada_con_unidades(Volume, "VolLiq", readOnly=True)
        lytCost.addWidget(self.Volumen, 6, 5)
        lytCost.addWidget(QtWidgets.QLabel(self.tr("Weight")), 7, 4)
        self.Peso = Entrada_con_unidades(Mass, readOnly=True)
        lytCost.addWidget(self.Peso, 7, 5)
        lytCost.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 2, 3, 6, 1)
        lytCost.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 8, 0, 1, 6)

        self.Costos = CostData(self.Equipment)
        self.Costos.valueChanged.connect(self.changeParamsCoste)
        lytCost.addWidget(self.Costos,9,1,2,5)

        lytCost.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 11, 0, 1, 6)
        group = QtWidgets.QGroupBox(self.tr("Stimated Costs"))
        lytCost.addWidget(group, 12, 1, 1, 5)
        lyt = QtWidgets.QGridLayout(group)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Purchase Cost")), 0, 1)
        self.C_adq = Entrada_con_unidades(Currency, retornar=False)
        self.C_adq.setReadOnly(True)
        lyt.addWidget(self.C_adq, 0, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Installed Cost")), 1, 1)
        self.C_inst = Entrada_con_unidades(Currency, retornar=False)
        self.C_inst.setReadOnly(True)
        self.C_inst.entrada.setReadOnly(True)
        lyt.addWidget(self.C_inst, 1, 2)

        if equipment:
            self.setEquipment(equipment)


if __name__ == "__main__":
    import sys, os
    from lib.corriente import Corriente, Mezcla
    app = QtWidgets.QApplication(sys.argv)

    locale = QtCore.QLocale.system().name()
    myTranslator = QtCore.QTranslator()
    if myTranslator.load("pychemqt_" + locale, os.environ["pychemqt"] + "i18n"):
        app.installTranslator(myTranslator)

    # agua=Corriente(300, 1, 3600, Mezcla([62], [1]))
    agua = Corriente(T=140+273.15, P=361540., caudalMasico=1.36, ids=[62],
                         fraccionMolar=[1.])
    tanque = Tank(entrada=agua, Di=1, L=2)
    dialogo = UI_equipment(tanque)
    dialogo.show()
    sys.exit(app.exec())
