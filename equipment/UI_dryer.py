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
# Solid dryer equipment dialog
###############################################################################


from functools import partial

from tools.qt import QtWidgets

from equipment.gas_solid_liquid import Dryer
from equipment.parents import UI_equip
from lib import unidades
from UI.widgets import Entrada_con_unidades


class UI_equipment(UI_equip):
    """Solid dryer equipment edition dialog"""
    Equipment = Dryer()

    def __init__(self, equipment=None, parent=None):
        """
        equipment: Initial equipment instance to model
        """
        super(UI_equipment, self).__init__(Dryer, parent=parent)

        # Input tab
        self.addEntrada(self.tr("Humid Solid"), "entradaSolido")
        self.addEntrada(self.tr("Air"), "entradaAire", psychro=True)

        # Calculate tab
        lyt = QtWidgets.QGridLayout(self.tabCalculo)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Mode")), 1, 1)
        self.mode = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_MODE:
            self.mode.addItem(txt)
        self.mode.currentIndexChanged.connect(
            partial(self.changeParams, "mode"))
        lyt.addWidget(self.mode, 1, 2, 1, 4)
        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 2, 1, 1, 6)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Air Relative Humidity")), 3, 1)
        self.HumedadAire = Entrada_con_unidades(
            float, max=1, spinbox=True, step=0.01)
        self.HumedadAire.valueChanged.connect(partial(self.changeParams, "HR"))
        lyt.addWidget(self.HumedadAire, 3, 2)
        lyt.addWidget(QtWidgets.QLabel(
            self.tr("Product moisture fraction")), 4, 1)
        self.HumedadSolido = Entrada_con_unidades(
            float, max=1., spinbox=True, step=0.01,
            textounidad=unidades.Mass(None).text()+"/"+unidades.Mass(None).text())
        self.HumedadSolido.valueChanged.connect(
            partial(self.changeParams, "HumedadResidual"))
        lyt.addWidget(self.HumedadSolido, 4, 2)
        lyt.addWidget(QtWidgets.QLabel(
            self.tr("Output Solid Temperature")), 5, 1)
        self.temperatura = Entrada_con_unidades(unidades.Temperature)
        self.temperatura.valueChanged.connect(
            partial(self.changeParams, "TemperaturaSolid"))
        lyt.addWidget(self.temperatura, 5, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Heat Duty")), 6, 1)
        self.Heat = Entrada_con_unidades(unidades.Power)
        self.Heat.valueChanged.connect(partial(self.changeParams, "Heat"))
        lyt.addWidget(self.Heat, 6, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Pressure Drop")), 7, 1)
        self.DeltaP = Entrada_con_unidades(unidades.Pressure)
        self.DeltaP.valueChanged.connect(partial(self.changeParams, "DeltaP"))
        lyt.addWidget(self.DeltaP, 7, 2)

        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 8, 1, 1, 6)
        group = QtWidgets.QGroupBox(self.tr("Results"))
        lyt.addWidget(group, 9, 1, 1, 5)
        layout = QtWidgets.QGridLayout(group)
        layout.addWidget(QtWidgets.QLabel(self.tr("Output Temperature")), 1, 1)
        self.temperaturaCalculada = Entrada_con_unidades(
            unidades.Temperature, retornar=False, readOnly=True)
        layout.addWidget(self.temperaturaCalculada, 1, 2)
        layout.addWidget(QtWidgets.QLabel(self.tr("Air Flow")), 2, 1)
        self.caudalVolumetrico = Entrada_con_unidades(
            unidades.VolFlow, "QGas", retornar=False, readOnly=True)
        layout.addWidget(self.caudalVolumetrico, 2, 2)
        layout.addWidget(QtWidgets.QLabel(
            self.tr("Output Air Relative Humidity")), 3, 1)
        self.HumedadCalculada = Entrada_con_unidades(
            float, readOnly=True, textounidad="%")
        layout.addWidget(self.HumedadCalculada, 3, 2)

        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 11, 1, 1, 6)

        # Output Tab
        self.addSalida(self.tr("Air"), psychro=True)
        self.addSalida(self.tr("Dry solid"))

        if equipment:
            self.setEquipment(equipment)

#    def rellenar(self):
#        self.EntradaAire.setCorriente(self.Equipment.kwargs["entradaAire"])
#        self.EntradaSolido.setCorriente(self.Equipment.kwargs["entradaSolido"])
#        if self.Equipment.status:
#            self.temperaturaCalculada.setValue(self.Equipment.SalidaSolido.T)
#            self.caudalVolumetrico.setValue(self.Equipment.entradaAire.corriente.Q)
#            self.HumedadCalculada.setValue(self.Equipment.SalidaAire.Xw*100)
#            self.SalidaAire.setCorriente(self.Equipment.SalidaAire)
#            self.SalidaSolido.setCorriente(self.Equipment.SalidaSolido)
#            if self.Equipment.kwargs["mode"]==1:
#                self.EntradaAire.setCorriente(self.Equipment.entradaAire)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
#    from lib.corriente import Mezcla, Corriente, Solid, PsyStream
#    from lib.psycrometry import PsychroState
#    diametros=[96.5, 105, 110, 118, 125, 130, 140, 150, 170]
#    fraccion=[0.02, 0.05, 0.1, 0.15, 0.25, 0.2, 0.15, 0.05, 0.03]
#    solido=Solid(caudalSolido=[5000], distribucion_fraccion=fraccion, distribucion_diametro=diametros)
#    Solido=Corriente(T=300, P=101325., caudalMasico=50, ids=[62], fraccionMolar=[1], solido=solido)
#    Aire=PsyStream(caudalMasico=100, tdb=300, HR=50)
#    secador=Dryer(entradaSolido=Solido, entradaAire=Aire)
    dialogo = UI_equipment()
    dialogo.show()
    sys.exit(app.exec())
