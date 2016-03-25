#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2016, Juan José Gómez Romera <jjgomera@gmail.com>

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
# Library for equipment common widget
#   * FoulingWidget: pipe fouling input data
#   * Dialog_Finned: Finned tube definition
###############################################################################

from PyQt5 import QtCore, QtWidgets


from lib.unidades import Fouling, Length, ThermalConductivity
from lib.pipeDatabase import finnedTube_database
from UI.widgets import Entrada_con_unidades


class FoulingWidget(QtWidgets.QWidget):
    """Widget con los parametros de fouling de tuberias"""
    valueChanged = QtCore.pyqtSignal(float)
    Fouling_Factor = {
        "Industrial": {
            "Fuel oil no.2": 0.000352,
            "Fuel oil no.6": 0.000881,
            "Transformer oil": 0.000173,
            "Engine Lube oil": 0.000173,
            "Quench oil": 0.000705,
            "Manufactured gas": 0.001761,
            "Engine exhaust gas": 0.001761,
            "Steam (nonoil bearing)": 0.000088,
            "Exhaust steam (oil bearing)": 0.0003,
            "Refrigerant vapors (Oil bearing)": 0.000352,
            "Compressed air": 0.000176,
            "Ammonia vapor": 0.000176,
            "CO2 vapor": 0.000176,
            "Chlorine vapor": 0.000352,
            "Coal flue gas": 0.001761,
            "Natural gas flue gas": 0.000881,
            "Molten heat transfer salts": 0.000088,
            "Refrigerant liquids": 0.000176,
            "Hydraulic fluid": 0.000176,
            "Industrial organic heat transfer media": 0.000352,
            "Ammonia liquid": 0.000176,
            "Ammonia liquid (oil bearing)": 0.000528,
            "Calcium chloride solutions": 0.000528,
            "Sodium chloride solutions": 0.000528,
            "CO2 liquid": 0.000176,
            "Chlorine liquid": 0.000352,
            "Methanol solutions": 0.000352,
            "Ethanol solutions": 0.000352,
            "Ethilene glycol solutions": 0.000352},
        "Chemical": {
            "Acid gases": 0.00044,
            "Solvent vapors": 0.000176,
            "Stable overhead products": 0.000176,
            "MEA and DEA solutions": 0.000352,
            "DEG and TEG solutions": 0.000352,
            "Stable side draw and bottom product": 0.00026,
            "Caustic solutions": 0.000352,
            "Vegetable oils": 0.000528},
        "Natural Gas-Gasoline": {
            "Natural gas": 0.00026,
            "Overhead products": 0.00026,
            "Lean oil": 0.000352,
            "Rich oil": 0.00026,
            "Natural gasoline": 0.00026,
            "Liquified petroleum gases": 0.00026},
        "Water (T<50C, v<0.9)": {
            "Seawater": 0.000088,
            "Brackish water": 0.000352,
            "Cooling tower (treated)": 0.000176,
            "Cooling tower (untreaterd)": 0.000528,
            "City or well water": 0.000176,
            "River water minimum": 0.000352,
            "River water average": 0.000528,
            "Muddy or silty": 0.000528,
            "Hard (>15 grains/gal": 0.000528,
            "Engine jacket": 0.000176,
            "Distilled, condensate": 0.000088,
            "Distilled, boiler blowdown": 0.000352,
            "Distilled, treated boiler feedwater": 0.000176},
        "Water (T<50C, v>0.9)": {
            "Seawater": 0.000088,
            "Brackish water": 0.000176,
            "Cooling tower (treated)": 0.000176,
            "Cooling tower (untreaterd)": 0.000528,
            "City or well water": 0.000176,
            "River water minimum": 0.000176,
            "River water average": 0.000352,
            "Muddy or silty": 0.000352,
            "Hard (>15 grains/gal": 0.000528,
            "Engine jacket": 0.000176,
            "Distilled, condensate": 0.000088,
            "Distilled, boiler blowdown": 0.000352,
            "Distilled, treated boiler feedwater": 0.000088},
        "Water (T>50C, v<0.9)": {
            "Seawater": 0.000176,
            "Brackish water": 0.000528,
            "Cooling tower (treated)": 0.000352,
            "Cooling tower (untreaterd)": 0.000881,
            "City or well water": 0.000352,
            "River water minimum": 0.000528,
            "River water average": 0.000705,
            "Muddy or silty": 0.000705,
            "Hard (>15 grains/gal": 0.000881,
            "Engine jacket": 0.000176,
            "Distilled, condensate": 0.000088,
            "Distilled, boiler blowdown": 0.000352,
            "Distilled, treated boiler feedwater": 0.000176},
        "Water (T>50C, v>0.9)": {
            "Seawater": 0.000176,
            "Brackish water": 0.000352,
            "Cooling tower (treated)": 0.000352,
            "Cooling tower (untreaterd)": 0.000705,
            "City or well water": 0.000352,
            "River water minimum": 0.000352,
            "River water average": 0.000528,
            "Muddy or silty": 0.000528,
            "Hard (>15 grains/gal": 0.000881,
            "Engine jacket": 0.000176,
            "Distilled, condensate": 0.000088,
            "Distilled, boiler blowdown": 0.000352,
            "Distilled, treated boiler feedwater": 0.000176},
        "Refinery vapors": {
            "Atmospheric tower overhead vapors": 0.000176,
            "Light naphthas": 0.000176,
            "Vacuum overhead vapors": 0.000352},
        "Refinery liq.": {
            "Crude oil dry T<120C, v<0.6": 0.000528,
            "Crude oil salt T<120C, v<0.6": 0.000528,
            "Crude oil dry T<120C, 0.6<v<1.2": 0.000352,
            "Crude oil salt T<120C, 0.6<v<1.2": 0.000352,
            "Crude oil dry T<120C, 1.2<v": 0.000352,
            "Crude oil salt T<120C, 1.2<v": 0.000352,
            "Crude oil dry 120C<T<175C, v<0.6": 0.000528,
            "Crude oil salt 120C<T<175C, v<0.6": 0.000881,
            "Crude oil dry 120C<T<175C, 0.6<v<1.2": 0.000352,
            "Crude oil salt 120C<T<175C, 0.6<v<1.2": 0.000705,
            "Crude oil dry 120C<T<175C, 1.2<v": 0.000352,
            "Crude oil salt 120C<T<175C, 1.2<v": 0.000705,
            "Crude oil dry 175C<T<230C, v<0.6": 0.000705,
            "Crude oil salt 175C<T<230C, v<0.6": 0.001057,
            "Crude oil dry 175C<T<230C, 0.6<v<1.2": 0.000528,
            "Crude oil salt 175C<T<230C, 0.6<v<1.2": 0.000881,
            "Crude oil dry 175C<T<230C, 1.2<v": 0.000528,
            "Crude oil salt 175C<T<230C, 1.2<v": 0.000881,
            "Crude oil dry T>230C, v<0.6": 0.000881,
            "Crude oil salt T>230C, v<0.6": 0.001233,
            "Crude oil dry T>230C, 0.6<v<1.2": 0.000705,
            "Crude oil salt T>230C, 0.6<v<1.2": 0.001057,
            "Crude oil dry T>230C, 1.2<v": 0.000705,
            "Crude oil salt T>230C, 1.2<v": 0.001057,
            "Gasoline": 0.000352,
            "Naphtha and light distillates": 0.00044,
            "Kerosene": 0.00044,
            "Light gas oil": 0.00044,
            "Heavy gas oil": 0.00067,
            "Heavy fuel oils": 0.00105},
        "Refinery Asphalt": {
            "Vacuum tower bottoms": 0.001761,
            "Atmosphere tower bottoms": 0.001233},
        "Refinery Cracking and caking": {
            "Overhead vapors": 0.000352,
            "Light cycle oil": 0.00044,
            "Heavy cycle oil": 0.00061,
            "Light coker gas oil": 0.00061,
            "Heavy coker gas oil": 0.00079,
            "Bottoms slurry oil": 0.000528,
            "Light liquid products": 0.000176},
        "Refinery Reforming": {
            "Reformer charge": 0.000264,
            "Reformer effluent": 0.000264,
            "Hydrocracker charge and effluent": 0.000352,
            "Recycle gas": 0.000176,
            "Overhead vapors": 0.000176,
            "Liquid product >50 API": 0.000176,
            "Liquid product 30-50 API": 0.000352},
        "Refinery Light Ends": {
            "Overhead vapors and gases": 0.000176,
            "Liquid products": 0.000176,
            "Absorption oils": 0.00044,
            "Alkylation trace acid streams": 0.000352,
            "Reboiler streams": 0.00044},
        "Refinery Lube oil": {
            "Feed stock": 0.000352,
            "Solvent feed mix": 0.000352,
            "Solvent": 0.000176,
            "Extract": 0.000528,
            "Rafftnate": 0.000176,
            "Asphalt": 0.000881,
            "Wax slurries": 0.000528,
            "Refined lube oil": 0.000176},
        "Refinery Visbreaker": {
            "Overhead vapor": 0.000528,
            "Visbreaker bottoms": 0.001761},
        "Refinery Naphtha Hydrotreater": {
            "Feed": 0.000528,
            "Effluent": 0.000352,
            "Naphfthas": 0.000352,
            "Overhead vapors": 0.000264},
        "Refinery Catalytic": {
            "Charge": 0.00079,
            "Effluent": 0.000352,
            "H.T. separator": 0.000352,
            "Stripper charge": 0.000528,
            "Liquid products": 0.000352},
        "Refinery HF Alky": {
            "Alkylate, deprop. bottons, main fract": 0.000528,
            "Other": 0.000352}}

    def __init__(self, parent=None):
        super(FoulingWidget, self).__init__(parent)
        layout = QtWidgets.QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.list = QtWidgets.QComboBox()
        self.list.addItem("")
        layout.addWidget(self.list)
        self.value = Entrada_con_unidades(Fouling, decimales=6)
        self.value.valueChanged.connect(self.valueChanged.emit)
        layout.addWidget(self.value)

        for tipo in sorted(self.Fouling_Factor):
            self.list.insertSeparator(self.list.count()+1)
            for componente in sorted(self.Fouling_Factor[tipo]):
                self.list.addItem(" - ".join([tipo, componente]))
        self.list.currentIndexChanged["QString"].connect(self.rellenar)

    def setValue(self, value):
        self.value.setValue(value)

    def rellenar(self, txt):
        if txt:
            tipo, componente = txt.split(" - ")
            value = self.Fouling_Factor[str(tipo)][str(componente)]
            self.value.setReadOnly(True)
            self.value.setValue(value)
            self.valueChanged.emit(value)
        else:
            self.value.setReadOnly(False)


class Dialog_Finned(QtWidgets.QDialog):
    """Dialog to define finned tube properties"""
    def __init__(self, kwarg=None, parent=None):
        super(Dialog_Finned, self).__init__(parent=parent)
        self.setWindowTitle(QtWidgets.QApplication.translate(
            "pychemqt", "Specify tube finned characteristics"))
        layout = QtWidgets.QGridLayout(self)
        self.listTube = QtWidgets.QComboBox()
        self.listTube.addItem("")
        layout.addWidget(self.listTube, 0, 1, 1, 2)

        layout.addItem(QtWidgets.QSpacerItem(10, 10, QtWidgets.QSizePolicy.Fixed,
                                         QtWidgets.QSizePolicy.Fixed), 1, 1, 1, 2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Material")), 2, 1)
        self.listMaterial = QtWidgets.QComboBox()
        self.listMaterial.addItem("")
        self.listMaterial.addItem(QtWidgets.QApplication.translate(
            "pychemqt", "Carbon Steel"))
        layout.addWidget(self.listMaterial, 2, 2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Thermal Conductivity")), 3, 1)
        self.kFin = Entrada_con_unidades(ThermalConductivity)
        layout.addWidget(self.kFin, 3, 2)
        layout.addItem(QtWidgets.QSpacerItem(10, 10, QtWidgets.QSizePolicy.Fixed,
                                         QtWidgets.QSizePolicy.Fixed), 4, 1, 1, 2)

        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Root diameter")), 5, 1)
        self.RootD = Entrada_con_unidades(Length, "PipeDiameter")
        layout.addWidget(self.RootD, 5, 2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Fin Height")), 6, 1)
        self.hFin = Entrada_con_unidades(Length, "Thickness")
        layout.addWidget(self.hFin, 6, 2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Base Fin Thickness")), 7, 1)
        self.BaseThickness = Entrada_con_unidades(Length, "Thickness")
        layout.addWidget(self.BaseThickness, 7, 2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Top Fin Thickness")), 8, 1)
        self.TopThickness = Entrada_con_unidades(Length, "Thickness")
        layout.addWidget(self.TopThickness, 8, 2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Number of fins")), 9, 1)
        self.Nfin = Entrada_con_unidades(float, textounidad="fins/m")
        layout.addWidget(self.Nfin, 9, 2)

        self.buttonBox = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Cancel |
                                                QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox, 10, 1, 1, 2)

        for tuberia in finnedTube_database:
            self.listTube.addItem("%s %s" % (tuberia[0], tuberia[1]))
        self.listTube.currentIndexChanged.connect(self.rellenarData)
        self.listTube.currentIndexChanged.connect(self.setDisabled)

        if kwarg:
            self.hFin.setValue(kwarg["hFin"])
            self.BaseThickness.setValue(kwarg["thicknessBaseFin"])
            self.TopThickness.setValue(kwarg["thicknessTopFin"])
            self.kFin.setValue(kwarg["kFin"])
            self.Nfin.setValue(kwarg["nFin"])
            self.RootD.setValue(kwarg["rootDoFin"])

    def rellenarData(self, ind):
        tuberia = finnedTube_database[ind-1]
        if tuberia[0] == "HPT":
            self.Nfin.setValue(int(tuberia[1][:2]))
            self.BaseThickness.setValue(tuberia[12]/1000.)
            self.TopThickness.setValue(tuberia[12]/1000.)
            self.RootD.setValue(tuberia[6]/1000.)
            self.hFin.setValue(tuberia[13]/1000.)

    def setDisabled(self, bool):
        self.RootD.setReadOnly(bool)
        self.BaseThickness.setReadOnly(bool)
        self.TopThickness.setReadOnly(bool)
        self.Nfin.setReadOnly(bool)
        self.hFin.setReadOnly(bool)

    def kwarg(self):
        kwarg = {"hFin": self.hFin.value,
                 "thicknessBaseFin": self.BaseThickness.value,
                 "thicknessTopFin": self.TopThickness.value,
                 "kFin": self.kFin.value,
                 "nFin": self.Nfin.value,
                 "rootDoFin": self.RootD.value}
        return kwarg


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    #dialogo = Dialog_Finned()
    dialogo = FoulingWidget()
    dialogo.show()
    sys.exit(app.exec_())
