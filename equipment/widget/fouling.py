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
# Library for equipment common widget
#   * FoulingWidget: pipe fouling input data
###############################################################################


from tools.qt import QtCore, QtWidgets, translate

from lib.unidades import Fouling
from UI.widgets import Entrada_con_unidades
from equipment.widget.gui import ToolGui


__doi__ = {
    # 1:
    #     {"autor": "",
    #      "title": "",
    #      "ref": "",
    #      "doi": ""},
    }


Material_Fouling_Factor = {
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


class UI_Fouling(ToolGui):
    """Fouling factor for pipes dialog"""

    title = translate("equipment", "Use fouling factor")

    def loadUI(self):
        """Add widget"""
        lyt = self.layout()

        self.fouling = FoulingWidget()
        lyt.addWidget(self.fouling, 2, 1, 1, 2)


class FoulingWidget(QtWidgets.QWidget):
    """Widget with fouling factor for pipes"""
    valueChanged = QtCore.pyqtSignal(float)

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QtWidgets.QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.list = QtWidgets.QComboBox()
        self.list.addItem("")
        layout.addWidget(self.list)
        self.value = Entrada_con_unidades(Fouling, decimales=6)
        self.value.valueChanged.connect(self.valueChanged.emit)
        layout.addWidget(self.value)

        for tipo in sorted(Material_Fouling_Factor):
            self.list.insertSeparator(self.list.count()+1)
            for componente in sorted(Material_Fouling_Factor[tipo]):
                self.list.addItem(" - ".join([tipo, componente]))
        self.list.currentTextChanged.connect(self.rellenar)

    def setValue(self, value):
        self.value.setValue(value)

    def rellenar(self, txt):
        if txt:
            tipo, componente = txt.split(" - ")
            value = Material_Fouling_Factor[str(tipo)][str(componente)]
            self.value.setReadOnly(True)
            self.value.setValue(value)
            self.valueChanged.emit(value)
        else:
            self.value.setReadOnly(False)


class Dialog(QtWidgets.QDialog):
    """Finned pipe dialog"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(self.tr("Finned pipe"))
        lyt = QtWidgets.QVBoxLayout(self)
        self.datos = UI_Fouling()
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
