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


from math import pi

from tools.qt import QtWidgets, translate

from lib.unidades import Dimensionless, Length, ThermalConductivity
from lib.pipeDatabase import finnedTube_database
from UI.widgets import Entrada_con_unidades
from equipment.accesories.gui import ToolGui


__doi__ = {
    # 1:
    #     {"autor": "",
    #      "title": "",
    #      "ref": "",
    #      "doi": ""},
        }

class FinnedPipe():
    """Finned pipe"""

    def __init__(self, H, D, delta):
        """
        Definition of twisted tape inserts

        Parameters
        ----------
        H : float
            Tape pitch for twist of π radians (180º), [m]
        D : float
            Internal diameter of tube, [m]
        delta : float
            Tape thickness, [m]
        """
        # Area tube without tape
        self.A = pi*D**2/4

        # Tape twist parameter
        self.y = Dimensionless(H/D)


class UI_FinnedPipe(ToolGui):
    """Finned pipe dialog"""

    title = translate("equipment", "Use finned tube")

    def loadUI(self):
        """Add widget"""
        lyt = self.layout()

        lytH = QtWidgets.QHBoxLayout()
        lytH.addWidget(QtWidgets.QLabel(self.tr("Pipe from database")))
        self.listTube = QtWidgets.QComboBox()
        self.listTube.addItem("")
        lytH.addWidget(self.listTube)
        lyt.addLayout(lytH, 2, 1, 1, 2)

        lyt.addWidget(QtWidgets.QLabel(self.tr("Material")), 3, 1)
        self.listMaterial = QtWidgets.QComboBox()
        self.listMaterial.addItem("")
        self.listMaterial.addItem(self.tr("Carbon Steel"))
        lyt.addWidget(self.listMaterial, 3, 2)
        lyt.addWidget(QtWidgets.QLabel(
            self.tr("Thermal Conductivity")), 4, 1)
        self.kFin = Entrada_con_unidades(ThermalConductivity)
        lyt.addWidget(self.kFin, 4, 2)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 5, 1, 1, 2)

        lyt.addWidget(QtWidgets.QLabel(self.tr("Root diameter")), 6, 1)
        self.RootD = Entrada_con_unidades(Length, "PipeDiameter")
        lyt.addWidget(self.RootD, 6, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Fin Height")), 7, 1)
        self.hFin = Entrada_con_unidades(Length, "Thickness")
        lyt.addWidget(self.hFin, 7, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Base Fin Thickness")), 8, 1)
        self.BaseThickness = Entrada_con_unidades(Length, "Thickness")
        lyt.addWidget(self.BaseThickness, 8, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Top Fin Thickness")), 9, 1)
        self.TopThickness = Entrada_con_unidades(Length, "Thickness")
        lyt.addWidget(self.TopThickness, 9, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Number of fins")), 10, 1)
        self.Nfin = Entrada_con_unidades(float, textounidad="fins/m")
        lyt.addWidget(self.Nfin, 10, 2)

        for tuberia in finnedTube_database:
            self.listTube.addItem("%s %s" % (tuberia[0], tuberia[1]))
        self.listTube.currentIndexChanged.connect(self.rellenarData)

    def rellenarData(self, idx):
        self.setReadOnly(idx)
        tuberia = finnedTube_database[idx-1]
        if tuberia[0] == "HPT":
            self.Nfin.setValue(int(tuberia[1][:2]))
            self.BaseThickness.setValue(tuberia[12]/1000.)
            self.TopThickness.setValue(tuberia[12]/1000.)
            self.RootD.setValue(tuberia[6]/1000.)
            self.hFin.setValue(tuberia[13]/1000.)

    def setReadOnly(self, boolean):
        """Set readOnly state for widget it values are from database"""
        self.RootD.setReadOnly(boolean)
        self.BaseThickness.setReadOnly(boolean)
        self.TopThickness.setReadOnly(boolean)
        self.Nfin.setReadOnly(boolean)
        self.hFin.setReadOnly(boolean)

    def kwarg(self):
        kwarg = {"hFin": self.hFin.value,
                 "thicknessBaseFin": self.BaseThickness.value,
                 "thicknessTopFin": self.TopThickness.value,
                 "kFin": self.kFin.value,
                 "nFin": self.Nfin.value,
                 "rootDoFin": self.RootD.value}
        return kwarg


class Dialog(QtWidgets.QDialog):
    """Finned pipe dialog"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(self.tr("Finned pipe"))
        lyt = QtWidgets.QVBoxLayout(self)
        self.datos = UI_FinnedPipe()
        lyt.addWidget(self.datos)
        buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Cancel
            | QtWidgets.QDialogButtonBox.StandardButton.Ok)
        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)
        lyt.addWidget(buttonBox)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = Dialog()
    Dialog.show()
    sys.exit(app.exec())
