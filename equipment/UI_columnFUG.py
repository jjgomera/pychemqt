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
# Simplified distillation column equipment dialog
# Use the method Fenske-Underwood-Gilliland
###############################################################################


from functools import partial

from equipment.distillation import ColumnFUG
from equipment.parents import UI_equip
from lib.config import getComponents
from lib.unidades import Pressure, Volume, Length, Power, Density, Currency
from tools.costIndex import CostData
from tools.qt import QtWidgets
from UI.widgets import Entrada_con_unidades


class UI_equipment(UI_equip):
    """ColumnFUG equipment edition dialog"""
    Equipment = ColumnFUG()

    def __init__(self, equipment=None, parent=None):
        """
        equipment: Initial equipment instance to model
        """
        super().__init__(ColumnFUG, entrada=False, parent=parent)

        # Calculate tab
        lyt = QtWidgets.QGridLayout(self.tabCalculo)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Feed tray")), 2, 0)
        self.feed = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_FEED:
            self.feed.addItem(txt)
        self.feed.currentIndexChanged.connect(
            partial(self.changeParams, "feed"))
        lyt.addWidget(self.feed, 2, 1)

        lyt.addWidget(QtWidgets.QLabel(self.tr("Condenser")), 3, 0)
        self.condenser = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_CONDENSER:
            self.condenser.addItem(txt)
        self.condenser.currentIndexChanged.connect(
            partial(self.changeParams, "condenser"))
        lyt.addWidget(self.condenser, 3, 1)
        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 4, 0, 1, 5)

        group = QtWidgets.QGroupBox(self.tr("Key Components specification"))
        lyt.addWidget(group, 5, 0, 1, 5)
        layout = QtWidgets.QGridLayout(group)
        layout.addWidget(QtWidgets.QLabel(self.tr("Light")), 1, 1)
        self.LK = QtWidgets.QComboBox()
        layout.addWidget(self.LK, 1, 2)
        layout.addWidget(QtWidgets.QLabel(self.tr("Split in destilate")), 1, 4)
        self.LKsplit = Entrada_con_unidades(float, spinbox=True, max=1.)
        self.LKsplit.valueChanged.connect(
            partial(self.changeParams, "LKsplit"))
        layout.addWidget(self.LKsplit, 1, 5)

        layout.addWidget(QtWidgets.QLabel(self.tr("Heavy")), 2, 1)
        self.HK = QtWidgets.QComboBox()
        layout.addWidget(self.HK, 2, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            40, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 2, 3)
        layout.addWidget(QtWidgets.QLabel(self.tr("Spit in residue")), 2, 4)
        self.HKsplit = Entrada_con_unidades(float, spinbox=True, max=1.)
        self.HKsplit.valueChanged.connect(
            partial(self.changeParams, "HKsplit"))
        layout.addWidget(self.HKsplit, 2, 5)

        indices, nombres, M = getComponents()
        for i, nombre in enumerate(nombres):
            self.LK.addItem("%i - %s" % (i+1, nombre))
            self.HK.addItem("%i - %s" % (i+1, nombre))
        self.LK.setCurrentIndex(-1)
        self.HK.setCurrentIndex(-1)
        self.LK.currentIndexChanged.connect(partial(self.changeParams, "LK"))
        self.HK.currentIndexChanged.connect(partial(self.changeParams, "HK"))
        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 6, 0, 1, 5)

        lyt.addWidget(QtWidgets.QLabel(self.tr("Reflux ratio")), 7, 0)
        self.R = Entrada_con_unidades(float)
        self.R.valueChanged.connect(partial(self.changeParams, "R"))
        lyt.addWidget(self.R, 7, 1)
        lyt.addWidget(QtWidgets.QLabel("R/Rmin"), 8, 0)
        self.R_Rmin = Entrada_con_unidades(float)
        self.R_Rmin.valueChanged.connect(partial(self.changeParams, "R_Rmin"))
        lyt.addWidget(self.R_Rmin, 8, 1)

        lyt.addWidget(QtWidgets.QLabel(self.tr("Design Pressure")), 7, 3)
        self.Pd = Entrada_con_unidades(Pressure)
        self.Pd.valueChanged.connect(partial(self.changeParams, "Pd"))
        lyt.addWidget(self.Pd, 7, 4)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Pressure loss")), 8, 3)
        self.DeltaP = Entrada_con_unidades(Pressure)
        self.DeltaP.valueChanged.connect(partial(self.changeParams, "DeltaP"))
        lyt.addWidget(self.DeltaP, 8, 4)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 9, 0, 1, 5)
        self.buttonMcCabe = QtWidgets.QPushButton(self.tr("McCabe-Thiele"))
        self.buttonMcCabe.clicked.connect(self.mcCabe)
        lyt.addWidget(self.buttonMcCabe, 10, 0)

        groupBox_Calculo = QtWidgets.QGroupBox(self.tr("Results"))
        lyt.addWidget(groupBox_Calculo, 11, 0, 1, 5)
        layout = QtWidgets.QGridLayout(groupBox_Calculo)
        layout.addWidget(QtWidgets.QLabel(self.tr("Condenser Duty")), 0, 1)
        self.DutyCondenser = Entrada_con_unidades(Power, retornar=False)
        self.DutyCondenser.setReadOnly(True)
        layout.addWidget(self.DutyCondenser, 0, 2)
        layout.addWidget(QtWidgets.QLabel(self.tr("Reboiler Duty")), 1, 1)
        self.DutyReboiler = Entrada_con_unidades(Power, retornar=False)
        self.DutyReboiler.setReadOnly(True)
        layout.addWidget(self.DutyReboiler, 1, 2)
        layout.addWidget(QtWidgets.QLabel("Rmin"), 2, 1)
        self.Rmin = Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.Rmin, 2, 2)
        layout.addWidget(QtWidgets.QLabel(self.tr("Reflux ratio")), 3, 1)
        self.RCalculada = Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.RCalculada, 3, 2)

        layout.addWidget(QtWidgets.QLabel("Nmin"), 0, 4)
        self.Nmin = Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.Nmin, 0, 5)
        layout.addWidget(QtWidgets.QLabel(self.tr("Stages")), 1, 4)
        self.NTray = Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.NTray, 1, 5)
        layout.addWidget(QtWidgets.QLabel(self.tr("Feed stage")), 2, 4)
        self.N_feed = Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.N_feed, 2, 5)

        # Cost tab
        lyt = QtWidgets.QGridLayout(self.tabCostos)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Process")), 1, 1)
        self.proceso = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_PROCESS:
            self.proceso.addItem(txt)
        self.proceso.currentIndexChanged.connect(
            partial(self.changeParamsCoste, "proceso"))
        lyt.addWidget(self.proceso, 1, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Column tipe")), 2, 1)
        self.tipo = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_COLUMN:
            self.tipo.addItem(txt)
        self.tipo.currentIndexChanged.connect(self.mostrarSubclasificacion)
        lyt.addWidget(self.tipo, 2, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Material")), 3, 1)
        self.material = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_MATERIAL:
            self.material.addItem(txt)
        self.material.currentIndexChanged.connect(
            partial(self.changeParamsCoste, "material_columna"))
        lyt.addWidget(self.material, 3, 2)

        lyt.addItem(QtWidgets.QSpacerItem(
            30, 30, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 1, 3, 5, 1)

        self.groupBox_Pisos = QtWidgets.QGroupBox(self.tr("Tray column"))
        lyt.addWidget(self.groupBox_Pisos, 1, 4, 4, 2)
        layout = QtWidgets.QGridLayout(self.groupBox_Pisos)
        layout.addWidget(QtWidgets.QLabel(self.tr("Tray type")), 1, 1)
        self.tipoPisos = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_TRAY:
            self.tipoPisos.addItem(txt)
        self.tipoPisos.currentIndexChanged.connect(
            partial(self.changeParamsCoste, "tipo_pisos"))
        layout.addWidget(self.tipoPisos, 1, 2)
        layout.addWidget(QtWidgets.QLabel(self.tr("Material")), 2, 1)
        self.materialPisos = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_MATERIAL:
            self.materialPisos.addItem(txt)
        self.materialPisos.currentIndexChanged.connect(
            partial(self.changeParamsCoste, "material_pisos"))
        layout.addWidget(self.materialPisos, 2, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 3, 1, 1, 2)
        layout.addWidget(QtWidgets.QLabel(self.tr("Diameter")), 4, 1)
        self.diametroPisos = Entrada_con_unidades(Length)
        layout.addWidget(self.diametroPisos, 4, 2)
        layout.addWidget(QtWidgets.QLabel(self.tr("Stages")), 5, 1)
        self.NumeroPisos = Entrada_con_unidades(
            int, spinbox=True, min=1, step=1, width=50)
        layout.addWidget(self.NumeroPisos, 5, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 6, 1, 1, 2)

        self.groupBox_relleno = QtWidgets.QGroupBox(self.tr("Packed column"))
        lyt.addWidget(self.groupBox_relleno, 1, 4, 4, 2)
        layout = QtWidgets.QGridLayout(self.groupBox_relleno)
        layout.addWidget(QtWidgets.QLabel(self.tr("Volume")), 1, 1)
        self.VolumenRelleno = Entrada_con_unidades(Volume, "VolLiq")
        layout.addWidget(self.VolumenRelleno, 1, 2)
        layout.addWidget(QtWidgets.QLabel(self.tr("Unit Cost")), 2, 1)
        texto = "%s / %s" % (Currency(None).text(), Volume(None).text("VolLiq"))
        self.C_unit_relleno = Entrada_con_unidades(
            Currency, retornar=False, textounidad=texto)
        layout.addWidget(self.C_unit_relleno, 2, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 3, 1, 1, 2)

        lyt.addWidget(QtWidgets.QLabel(self.tr("Diameter")), 5, 1)
        self.Dc = Entrada_con_unidades(Length)
        lyt.addWidget(self.Dc, 5, 2, 1, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Height")), 6, 1)
        self.Hc = Entrada_con_unidades(Length)
        lyt.addWidget(self.Hc, 6, 2, 1, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Thickness (top)")), 6, 4)
        self.EspesorSuperior = Entrada_con_unidades(Length, "Thickness")
        lyt.addWidget(self.EspesorSuperior, 6, 5, 1, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Thickness (bottom)")), 7, 4)
        self.EspesorInferior = Entrada_con_unidades(Length, "Thickness")
        lyt.addWidget(self.EspesorInferior, 7, 5, 1, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Density")), 7, 1)
        self.EspesorInferior = Entrada_con_unidades(Density, "DenLiq")
        lyt.addWidget(self.EspesorInferior, 7, 2, 1, 2)

        self.Costos = CostData(self.Equipment)
        self.Costos.valueChanged.connect(self.changeParamsCoste)
        lyt.addWidget(self.Costos, 10, 1, 2, 5)

        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 12, 1, 1, 6)

        self.groupBox_Costos = QtWidgets.QGroupBox(self.tr("Stimated costs"))
        lyt.addWidget(self.groupBox_Costos, 13, 1, 1, 5)
        layout = QtWidgets.QGridLayout(self.groupBox_Costos)
        layout.addWidget(QtWidgets.QLabel(self.tr("Tray cost")), 0, 1)
        self.C_pisos = Entrada_con_unidades(Currency, retornar=False)
        self.C_pisos.setReadOnly(True)
        layout.addWidget(self.C_pisos, 0, 2)
        layout.addWidget(QtWidgets.QLabel(self.tr("Shell cost")), 1, 1)
        self.C_carcasa = Entrada_con_unidades(Currency, retornar=False)
        self.C_carcasa.setReadOnly(True)
        layout.addWidget(self.C_carcasa, 1, 2)
        layout.addWidget(QtWidgets.QLabel(
            self.tr("Platform and ladder")), 2, 1)
        self.C_accesorios = Entrada_con_unidades(Currency, retornar=False)
        self.C_accesorios.setReadOnly(True)
        layout.addWidget(self.C_accesorios, 2, 2)
        layout.addWidget(QtWidgets.QLabel(self.tr("Column cost")), 0, 4)
        self.C_columna = Entrada_con_unidades(Currency, retornar=False)
        self.C_columna.setReadOnly(True)
        layout.addWidget(self.C_columna, 0, 5)
        layout.addWidget(QtWidgets.QLabel(self.tr("Purchase costs")), 1, 4)
        self.C_adq = Entrada_con_unidades(Currency, retornar=False)
        self.C_adq.setReadOnly(True)
        layout.addWidget(self.C_adq, 1, 5)
        layout.addWidget(QtWidgets.QLabel(self.tr("Installed costs")), 2, 4)
        self.C_inst = Entrada_con_unidades(Currency, retornar=False)
        self.C_inst.setReadOnly(True)
        layout.addWidget(self.C_inst, 2, 5)

        # Output tab
        self.addSalida(self.tr("Destilate"))
        self.addSalida(self.tr("Residue"))

        self.mostrarSubclasificacion(0)
        if equipment:
            self.setEquipment(equipment)

    def mcCabe(self):
        self.Equipment.McCabe()

    def mostrarSubclasificacion(self, ind):
        self.groupBox_Pisos.setVisible(not ind)
        self.groupBox_relleno.setVisible(ind)
        self.changeParamsCoste("tipo", ind)

    def rellenar(self):
        self.buttonMcCabe.setEnabled(self.Equipment.statusMcCabe)
        UI_equip.rellenar(self)


if __name__ == "__main__":
    import sys
    from lib.corriente import Corriente
    app = QtWidgets.QApplication(sys.argv)
    blend = Corriente(T=340, P=101325, caudalMasico=1, ids=[11, 12],
                      fraccionMolar=[0.6, 0.4])
    columna = ColumnFUG(entrada=blend, LK=0, LKsplit=0.96666, HK=1,
                        HKsplit=0.95, feed=0)
    dialogo = UI_equipment(columna)
    dialogo.show()
    sys.exit(app.exec())
