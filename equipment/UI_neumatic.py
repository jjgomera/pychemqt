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
# Neumatic conveying equipment dialog
###############################################################################


from functools import partial

from equipment.neumatic import Neumatic
from equipment.parents import UI_equip
from equipment.UI_pipe import Catalogo_Materiales_Dialog, Catalogo_Accesorios_Dialog
from lib.unidades import DeltaP, Speed, Length
from tools.qt import QtWidgets
from UI.widgets import Entrada_con_unidades


class UI_equipment(UI_equip):
    """Neumatic conveying equipment dialog"""
    Equipment = Neumatic()

    def __init__(self, equipment=None, parent=None):
        """
        equipment: Initial equipment instance to model
        """
        super().__init__(Neumatic, entrada=False, salida=False, parent=parent)

        # Calculate tab
        lyt_Calc = QtWidgets.QGridLayout(self.tabCalculo)

        lyt = QtWidgets.QHBoxLayout()
        lyt.addWidget(QtWidgets.QLabel(self.tr("Saltation velocity method:")))
        self.saltation = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_SALTATION:
            self.saltation.addItem(txt)
        self.saltation.currentIndexChanged.connect(
            partial(self.changeParams, "saltation"))
        lyt.addWidget(self.saltation)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed))
        lyt_Calc.addLayout(lyt, 1, 1, 1, 4)

        lyt_Calc.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 2, 1, 1, 4)

        lyt_Calc.addWidget(QtWidgets.QLabel(self.tr("Pipe Diameter")), 3, 1)
        self.D = Entrada_con_unidades(Length, "PipeDiameter")
        self.D.valueChanged.connect(partial(self.changeParams, "D"))
        lyt_Calc.addWidget(self.D, 3, 2)
        buttonPipe = QtWidgets.QPushButton(self.tr("Pipe Database"))
        buttonPipe.clicked.connect(self.showMaterial)
        lyt_Calc.addWidget(buttonPipe, 3, 3)

        lyt_Calc.addWidget(QtWidgets.QLabel(self.tr("Pipe Roughness")), 4, 1)
        self.eD = Entrada_con_unidades(Length, "ParticleDiameter")
        self.eD.valueChanged.connect(partial(self.changeParams, "eD"))
        lyt_Calc.addWidget(self.eD, 4, 2)
        lyt_Calc.addWidget(QtWidgets.QLabel(self.tr("Pipe Length")), 5, 1)
        self.L = Entrada_con_unidades(Length)
        self.L.valueChanged.connect(partial(self.changeParams, "L"))
        lyt_Calc.addWidget(self.L, 5, 2)

        lyt_Calc.addWidget(QtWidgets.QLabel("K"), 6, 1)
        self.K = Entrada_con_unidades(float)
        self.K.valueChanged.connect(partial(self.changeParams, "K"))
        lyt_Calc.addWidget(self.K, 6, 2)
        buttonfitting = QtWidgets.QPushButton(self.tr("Pipe Fittings"))
        buttonfitting.clicked.connect(self.showFitting)
        lyt_Calc.addWidget(buttonfitting, 6, 3)
        lyt_Calc.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 10, 5)

        group = QtWidgets.QGroupBox()
        group.setTitle(self.tr("Results"))
        lyt_Calc.addWidget(group, 12, 1, 1, 5)
        lyt = QtWidgets.QGridLayout(group)
        label = QtWidgets.QLabel(self.tr("SLR"))
        label.setToolTip(self.tr("Solid Loading Ratio"))
        lyt.addWidget(label, 1, 1)
        self.SLR = Entrada_con_unidades(float, retornar=False, readOnly=True)
        lyt.addWidget(self.SLR, 1, 2)
        label = QtWidgets.QLabel(self.tr("DR"))
        label.setToolTip(self.tr("Density Ratio"))
        lyt.addWidget(label, 2, 1)
        self.DR = Entrada_con_unidades(float, retornar=False, readOnly=True)
        lyt.addWidget(self.DR, 2, 2)
        label = QtWidgets.QLabel(self.tr("VF"))
        label.setToolTip(self.tr("Volume Fraction"))
        lyt.addWidget(label, 3, 1)
        self.VF = Entrada_con_unidades(float, retornar=False, readOnly=True)
        lyt.addWidget(self.VF, 3, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Saltation velocity")), 4, 1)
        self.vs = Entrada_con_unidades(Speed, retornar=False, readOnly=True)
        lyt.addWidget(self.vs, 4, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Real velocity")), 5, 1)
        self.V = Entrada_con_unidades(Speed, retornar=False, readOnly=True)
        lyt.addWidget(self.V, 5, 2)

        lyt.addWidget(QtWidgets.QLabel("Re"), 1, 4)
        self.Re = Entrada_con_unidades(float, readOnly=True)
        lyt.addWidget(self.Re, 1, 5)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Friction factor")), 2, 4)
        self.f = Entrada_con_unidades(float, readOnly=True)
        lyt.addWidget(self.f, 2, 5)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Pressure loss")), 3, 4)
        self.DeltaP = Entrada_con_unidades(DeltaP, readOnly=True)
        lyt.addWidget(self.DeltaP, 3, 5)

        if equipment:
            self.setEquipment(equipment)

    def showMaterial(self):
        """Show pipe catalog to select it to get a nominal diameter"""
        dialogo = Catalogo_Materiales_Dialog()
        if dialogo.exec():
            material = dialogo.getMaterial()
            if material:
                # Convert internal diameter from mm to m
                diameter = material[4]/1000
                eD = material[2]/1000
                self.D.setValue(diameter)
                self.eD.setValue(eD)
                self.calculo(**{"D": diameter, "eD": eD})

    def showFitting(self):
        """Show pipe fitting dialog"""
        dialogo = Catalogo_Accesorios_Dialog()
        if dialogo.exec():
            self.K.setValue(dialogo.K)
            self.changeParams("K", dialogo.K)


if __name__ == "__main__":
    import sys
    from lib.solids import Solid
    from lib.corriente import Corriente

    app = QtWidgets.QApplication(sys.argv)
    dm = [17.5e-6, 22.4e-6, 26.2e-6, 31.8e-6, 37e-6, 42.4e-6, 48e-6, 54e-6,
          60e-6, 69e-6, 81.3e-6, 96.5e-6, 109e-6, 127e-6]
    fracciones = [0.02, 0.03, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                  0.05, 0.03, 0.02]
    sol = Solid(caudalSolido=[0.1], distribucion_diametro=dm,
                distribucion_fraccion=fracciones, solids=[638])
    kw = {"ids": [475], "fraccionMolar": [1.], "MEoS": True}
    entrada = Corriente(T=300, P=1e5, caudalMasico=1, solido=sol, **kw)
    neumatic = Neumatic(entrada=entrada)
    dialogo = UI_equipment(neumatic)
    dialogo.show()
    sys.exit(app.exec())
