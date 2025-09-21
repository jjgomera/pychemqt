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
# Double pipe heat exchanger equipment dialog
###############################################################################


from functools import partial

from lib.unidades import (Temperature, Pressure, DeltaP, Power, Length, Area,
                          ThermalConductivity, HeatTransfCoef, Currency)
from UI.widgets import Entrada_con_unidades
from equipment.heatExchanger import Hairpin
from equipment.UI_pipe import PipeCatalogDialog
from equipment.parents import UI_equip
from equipment.widget import FoulingWidget, Dialog_Finned
from tools.costIndex import CostData
from tools.qt import QtWidgets


class UI_equipment(UI_equip):
    """Double pipe heat exchanger edition dialog"""
    Equipment = Hairpin()

    def __init__(self, equipment=None, parent=None):
        """
        equipment: Initial equipment instance to model
        """
        super().__init__(Hairpin, parent=parent)

        # Input tab
        self.addEntrada(self.tr("Tube"), "entradaTubo")
        self.addEntrada(self.tr("Annulli"), "entradaExterior")

        # Pipe catalog tab
        tabCatalogo = QtWidgets.QWidget()
        self.tabWidget.insertTab(1, tabCatalogo, self.tr("Catalog"))
        lyt = QtWidgets.QGridLayout(tabCatalogo)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Tube length")), 4, 1)
        self.LTube = Entrada_con_unidades(Length)
        self.LTube.valueChanged.connect(partial(self.changeParams, "LTube"))
        lyt.addWidget(self.LTube, 4, 2)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 5, 1)
        lyt.addWidget(QtWidgets.QLabel(
            self.tr("Tube internal diameter")), 6, 1)
        self.DiTube = Entrada_con_unidades(Length, "pipeDiameter")
        self.DiTube.valueChanged.connect(partial(self.changeParams, "DiTube"))
        lyt.addWidget(self.DiTube, 6, 2)
        lyt.addWidget(QtWidgets.QLabel(
            self.tr("Tube external diameter")), 7, 1)
        self.DeTube = Entrada_con_unidades(Length, "pipeDiameter")
        self.DeTube.valueChanged.connect(partial(self.changeParams, "DeTube"))
        lyt.addWidget(self.DeTube, 7, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Tube thickness")), 8, 1)
        self.wTube = Entrada_con_unidades(Length, "Thickness")
        self.wTube.valueChanged.connect(partial(self.changeParams, "wTube"))
        lyt.addWidget(self.wTube, 8, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Tube roughness")), 9, 1)
        self.rTube = Entrada_con_unidades(Length, "Thickness")
        self.rTube.valueChanged.connect(partial(self.changeParams, "rTube"))
        lyt.addWidget(self.rTube, 9, 2)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 10, 1)
        lyt.addWidget(QtWidgets.QLabel(
            self.tr("Annulli external diameter")), 11, 1)
        self.DeeTube = Entrada_con_unidades(Length, "pipeDiameter")
        self.DeeTube.valueChanged.connect(
            partial(self.changeParams, "DeeTube"))
        lyt.addWidget(self.DeeTube, 11, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Thermal conductivity")), 12, 1)
        self.kTube = Entrada_con_unidades(ThermalConductivity)
        self.kTube.valueChanged.connect(partial(self.changeParams, "kTube"))
        lyt.addWidget(self.kTube, 12, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Tube Count")), 13, 1)
        self.nTube = Entrada_con_unidades(int)
        self.nTube.valueChanged.connect(partial(self.changeParams, "nTube"))
        lyt.addWidget(self.nTube, 13, 2)

        buttonPipe = QtWidgets.QPushButton(self.tr("Pipe Database"))
        buttonPipe.clicked.connect(self.showMaterial)
        lyt.addWidget(buttonPipe, 6, 3, 4, 1)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 14, 1)
        self.tubeFinned = QtWidgets.QCheckBox(self.tr("Finned Tube"))
        lyt.addWidget(self.tubeFinned, 15, 1, 1, 4)
        self.buttonFin = QtWidgets.QPushButton(
            self.tr("Finned Pipe Database"))
        self.buttonFin.setEnabled(False)
        self.buttonFin.clicked.connect(self.showFinTube)
        lyt.addWidget(self.buttonFin, 15, 3)
        self.tubeFinned.toggled.connect(
            partial(self.changeParams, "tubeFinned"))
        self.tubeFinned.toggled.connect(self.buttonFin.setEnabled)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Inside Fouling")), 16, 1)
        self.tubeFouling = FoulingWidget()
        self.tubeFouling.valueChanged.connect(
            partial(self.changeParams, "tubeFouling"))
        lyt.addWidget(self.tubeFouling, 16, 2, 1, 5)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Outside Fouling")), 17, 1)
        self.annulliFouling = FoulingWidget()
        self.annulliFouling.valueChanged.connect(
            partial(self.changeParams, "annulliFouling"))
        lyt.addWidget(self.annulliFouling, 17, 2, 1, 5)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 20, 1, 1, 6)

        # Calculate tab
        lyt = QtWidgets.QGridLayout(self.tabCalculo)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Mode")), 1, 1)
        self.modo = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_MODO:
            self.modo.addItem(txt)
        self.modo.currentIndexChanged.connect(
            partial(self.changeParams, "modo"))
        lyt.addWidget(self.modo, 1, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Flujo")), 2, 1)
        self.flujo = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_FLUJO:
            self.flujo.addItem(txt)
        self.flujo.currentIndexChanged.connect(
            partial(self.changeParams, "flujo"))
        lyt.addWidget(self.flujo, 2, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Layout")), 3, 1)
        self.orientacion = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_ORIENTACION:
            self.orientacion.addItem(txt)
        self.orientacion.currentIndexChanged.connect(
            partial(self.changeParams, "orientacion"))
        lyt.addWidget(self.orientacion, 3, 2)

        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 4, 1)

        lyt.addWidget(QtWidgets.QLabel(
            self.tr("Output inside temperature")), 5, 1)
        self.tubeTout = Entrada_con_unidades(Temperature)
        self.tubeTout.valueChanged.connect(
            partial(self.changeParams, "tubeTout"))
        lyt.addWidget(self.tubeTout, 5, 2)
        lyt.addWidget(QtWidgets.QLabel(
            self.tr("Output annulli temperature")), 6, 1)
        self.annulliTout = Entrada_con_unidades(Temperature)
        self.annulliTout.valueChanged.connect(
            partial(self.changeParams, "annulliTout"))
        lyt.addWidget(self.annulliTout, 6, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Output inside quality")), 5, 4)
        self.tubeXout = Entrada_con_unidades(float)
        self.tubeXout.valueChanged.connect(
            partial(self.changeParams, "tubeXout"))
        lyt.addWidget(self.tubeXout, 5, 5)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Output annulli quality")), 6, 4)
        self.annulliXout = Entrada_con_unidades(float)
        self.annulliXout.valueChanged.connect(
            partial(self.changeParams, "annulliXout"))
        lyt.addWidget(self.annulliXout, 6, 5)

        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 15, 1, 1, 6)

        group = QtWidgets.QGroupBox(self.tr("Results"))
        lyt.addWidget(group, 16, 1, 1, 6)
        layout = QtWidgets.QGridLayout(group)
        layout.addWidget(QtWidgets.QLabel(self.tr("Heat Duty")), 0, 1)
        self.Q = Entrada_con_unidades(Power, retornar=False, readOnly=True)
        layout.addWidget(self.Q, 0, 2)
        layout.addWidget(QtWidgets.QLabel(self.tr("Tout Tube")), 1, 1)
        self.ToutTube = Entrada_con_unidades(Temperature, retornar=False)
        self.ToutTube.setReadOnly(True)
        layout.addWidget(self.ToutTube, 1, 2)
        layout.addWidget(QtWidgets.QLabel(self.tr("Tout Tube")), 2, 1)
        self.ToutAnnulli = Entrada_con_unidades(Temperature, retornar=False)
        self.ToutAnnulli.setReadOnly(True)
        layout.addWidget(self.ToutAnnulli, 2, 2)
        layout.addWidget(QtWidgets.QLabel(self.tr("U")), 0, 4)
        self.U = Entrada_con_unidades(HeatTransfCoef, retornar=False)
        self.U.setReadOnly(True)
        layout.addWidget(self.U, 0, 5)
        layout.addWidget(QtWidgets.QLabel(self.tr("Area")), 1, 4)
        self.A = Entrada_con_unidades(Area, retornar=False, readOnly=True)
        layout.addWidget(self.A, 1, 5)
        layout.addWidget(QtWidgets.QLabel(self.tr("Lenght")), 2, 4)
        self.L = Entrada_con_unidades(Length, retornar=False, readOnly=True)
        layout.addWidget(self.L, 2, 5)
        layout.addWidget(QtWidgets.QLabel(self.tr("DeltaP Tube")), 0, 7)
        self.deltaPTube = Entrada_con_unidades(DeltaP, retornar=False)
        self.deltaPTube.setReadOnly(True)
        layout.addWidget(self.deltaPTube, 0, 8)
        layout.addWidget(QtWidgets.QLabel(self.tr("DeltaP Annulli")), 1, 7)
        self.deltaPAnnulli = Entrada_con_unidades(DeltaP, retornar=False)
        self.deltaPAnnulli.setReadOnly(True)
        layout.addWidget(self.deltaPAnnulli, 1, 8)
        layout.addWidget(QtWidgets.QLabel(self.tr("CF")), 2, 7)
        self.CF = Entrada_con_unidades(float, retornar=False, readOnly=True)
        layout.addWidget(self.CF, 2, 8)

        lyt.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 17, 1, 1, 6)

        # Cost tab
        lyt = QtWidgets.QGridLayout(self.tabCostos)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Material")), 2, 1)
        self.material = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_MATERIAL:
            self.material.addItem(txt)
        self.material.currentIndexChanged.connect(
            partial(self.changeParamsCoste, "material"))
        lyt.addWidget(self.material, 2, 2)
        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 3, 0, 1, 6)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Design Pressure")), 4, 1)
        self.P_dis = Entrada_con_unidades(Pressure)
        self.P_dis.valueChanged.connect(
            partial(self.changeParamsCoste, "P_dis"))
        lyt.addWidget(self.P_dis, 4, 2, 1, 1)
        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 5, 0, 1, 6)

        self.Costos = CostData(self.Equipment)
        self.Costos.valueChanged.connect(self.changeParamsCoste)
        lyt.addWidget(self.Costos, 6, 1, 2, 5)

        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 8, 0, 1, 6)
        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 10, 0, 1, 6)
        group = QtWidgets.QGroupBox(self.tr("Stimated Costs"))
        lyt.addWidget(group, 9, 1, 1, 5)
        layout = QtWidgets.QGridLayout(group)
        layout.addWidget(QtWidgets.QLabel(self.tr("Purchase Cost")), 0, 1)
        self.C_adq = Entrada_con_unidades(Currency, retornar=False)
        self.C_adq.setReadOnly(True)
        layout.addWidget(self.C_adq, 0, 2)
        layout.addWidget(QtWidgets.QLabel(self.tr("Installed Cost")), 1, 1)
        self.C_inst = Entrada_con_unidades(Currency, retornar=False)
        self.C_inst.setReadOnly(True)
        self.C_inst.entrada.setReadOnly(True)
        layout.addWidget(self.C_inst, 1, 2)

        # Output Tab
        self.addSalida(self.tr("Tube"))
        self.addSalida(self.tr("Annulli"))

        if equipment:
            self.setEquipment(equipment)

    def showMaterial(self):
        dialogo = PipeCatalogDialog()
        if dialogo.exec():
            material = dialogo.getMaterial()
            if material:
                self.rTube.setValue(material[2])
                self.DiTube.setValue(material[4])
                self.wTube.setValue(material[5])
                self.DeTube.setValue(material[6])

    def showFinTube(self):
        dialogo = Dialog_Finned(self.Equipment.kwargs)
        if dialogo.exec():
            kwarg = dialogo.kwarg()
            self.calculo(**kwarg)


if __name__ == "__main__":
    import sys
    from lib.corriente import Corriente
    app = QtWidgets.QApplication(sys.argv)
    caliente = Corriente(T=140+273.15, P=361540., caudalMasico=1.36, ids=[62],
                         fraccionMolar=[1.])
    fria = Corriente(T=20+273.15, P=101325., caudalMasico=5000/3600., ids=[62],
                     fraccionMolar=[1.])
    Cambiador = Hairpin(
        entradaTubo=caliente, entradaExterior=fria,
        modo=1, DiTube=0.0525,
        DeTube=0.0603, DeeTube=0.0779, kTube=54, rTube=0.0459994e-3,
        annulliFouling=0.000352, tubeFouling=0.000176, LTube=2.5)
    dialogo = UI_equipment(Cambiador)
    dialogo.show()
    sys.exit(app.exec())
