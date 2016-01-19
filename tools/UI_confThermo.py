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
# Thermal method calculation config
###############################################################################

import os

from PyQt5 import QtWidgets


from lib.EoS import K, H


class UI_confThermo_widget(QtWidgets.QWidget):
    """Widget to config thermal method calculations"""
    def __init__(self, config=None, parent=None):
        """Constructor, opcional config parameter with proyect configuration"""
        super(UI_confThermo_widget, self).__init__(parent)
        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "K values:")), 0, 0, 1, 2)
        self.K = QtWidgets.QComboBox()
        for eq in K:
            self.K.addItem(eq.__title__)
        layout.addWidget(self.K, 0, 2)

        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt   ", "Alfa function:")), 1, 0, 1, 2)
        self.alfa = QtWidgets.QComboBox()
        self.alfa.addItem(QtWidgets.QApplication.translate("pychemqt", "Original"))
        self.alfa.addItem("Boston-Mathias")
        self.alfa.addItem("Twu")
        self.alfa.addItem("Doridon")
        layout.addWidget(self.alfa, 1, 2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Mix rules:")), 2, 0, 1, 2)
        self.mixing_rule = QtWidgets.QComboBox()
        self.mixing_rule.addItem("van der Waals")
        self.mixing_rule.addItem("Stryjek-Vera")
        self.mixing_rule.addItem("Panagiotopoulos")
        self.mixing_rule.addItem("Melhem")
        layout.addWidget(self.mixing_rule, 2, 2)
        layout.addItem(QtWidgets.QSpacerItem(20, 20, QtWidgets.QSizePolicy.Fixed,
                                         QtWidgets.QSizePolicy.Fixed), 3, 0, 1, 4)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Enthalpy:")), 4, 0, 1, 2)
        self.H = QtWidgets.QComboBox()
        for h in H:
            self.H.addItem(h.__title__)
        layout.addWidget(self.H, 4, 2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Ideal heat capacity:")), 5, 0, 1, 2)
        self.Cp_ideal = QtWidgets.QComboBox()
        self.Cp_ideal.addItem(QtWidgets.QApplication.translate("pychemqt", "Ideal"))
        self.Cp_ideal.addItem("DIPPR")
        layout.addWidget(self.Cp_ideal, 5, 2)
        layout.addItem(QtWidgets.QSpacerItem(20, 20, QtWidgets.QSizePolicy.Fixed,
                                         QtWidgets.QSizePolicy.Fixed), 6, 0, 1, 4)
        self.MEoS = QtWidgets.QCheckBox(QtWidgets.QApplication.translate(
            "pychemqt", "Use MEoS for single compounds if it's available"))
        layout.addWidget(self.MEoS, 7, 0, 1, 3)
        self.coolProp = QtWidgets.QCheckBox(QtWidgets.QApplication.translate(
            "pychemqt", "Use external library coolProp (faster)"))
        self.coolProp.setEnabled(False)
        layout.addWidget(self.coolProp, 8, 1, 1, 2)
        self.refprop = QtWidgets.QCheckBox(QtWidgets.QApplication.translate(
            "pychemqt", "Use external library refprop (fastest)"))
        self.refprop.setEnabled(False)
        layout.addWidget(self.refprop, 9, 1, 1, 2)

        self.iapws = QtWidgets.QCheckBox(QtWidgets.QApplication.translate(
            "pychemqt", "Use IAPWS97 for water"))
        layout.addWidget(self.iapws, 10, 0, 1, 3)
        self.freesteam = QtWidgets.QCheckBox(QtWidgets.QApplication.translate(
            "pychemqt", "Use freesteam library (faster)"))
        self.freesteam.setEnabled(False)
        layout.addWidget(self.freesteam, 11, 1, 1, 2)
        self.GERG = QtWidgets.QCheckBox(QtWidgets.QApplication.translate(
            "pychemqt", "Use GERG EoS for mix if it's posible"))
        layout.addWidget(self.GERG, 12, 0, 1, 3)

        if os.environ["freesteam"] == "True":
            self.iapws.toggled.connect(self.freesteam.setEnabled)
        if os.environ["CoolProp"] == "True":
            self.MEoS.toggled.connect(self.coolProp.setEnabled)
        if os.environ["refprop"] == "True":
            self.MEoS.toggled.connect(self.refprop.setEnabled)

        if config:
            self.setConfig(config)

    def setConfig(self, config):
        if config.has_section("Thermo"):
            self.K.setCurrentIndex(config.getint("Thermo", "K"))
            self.alfa.setCurrentIndex(config.getint("Thermo", "Alfa"))
            self.mixing_rule.setCurrentIndex(config.getint("Thermo", "Mixing"))
            self.H.setCurrentIndex(config.getint("Thermo", "H"))
            self.Cp_ideal.setCurrentIndex(config.getint("Thermo", "Cp_ideal"))
            self.MEoS.setChecked(config.getboolean("Thermo", "MEoS"))
            self.iapws.setChecked(config.getboolean("Thermo", "iapws"))
            self.GERG.setChecked(config.getboolean("Thermo", "GERG"))
            self.freesteam.setChecked(config.getboolean("Thermo", "freesteam"))
            self.coolProp.setChecked(config.getboolean("Thermo", "coolProp"))
            self.refprop.setChecked(config.getboolean("Thermo", "refprop"))

    def value(self, config):
        """Function result for wizard"""
        if not config.has_section("Thermo"):
            config.add_section("Thermo")
        config.set("Thermo", "K", str(self.K.currentIndex()))
        config.set("Thermo", "Alfa", str(self.alfa.currentIndex()))
        config.set("Thermo", "Mixing", str(self.mixing_rule.currentIndex()))
        config.set("Thermo", "H", str(self.H.currentIndex()))
        config.set("Thermo", "Cp_ideal", str(self.Cp_ideal.currentIndex()))
        config.set("Thermo", "MEoS", str(self.MEoS.isChecked()))
        config.set("Thermo", "iapws", str(self.iapws.isChecked()))
        config.set("Thermo", "GERG", str(self.GERG.isChecked()))
        config.set("Thermo", "freesteam", str(self.freesteam.isChecked()))
        config.set("Thermo", "coolProp", str(self.coolProp.isChecked()))
        config.set("Thermo", "refprop", str(self.refprop.isChecked()))
        return config

    @classmethod
    def default(cls, config):
        config.add_section("Thermo")
        config.set("Thermo", "K", "0")
        config.set("Thermo", "Alfa", "0")
        config.set("Thermo", "Mixing", "0")
        config.set("Thermo", "H", "0")
        config.set("Thermo", "Cp_ideal", "0")
        config.set("Thermo", "MEoS", "False")
        config.set("Thermo", "iapws", "False")
        config.set("Thermo", "GERG", "False")
        config.set("Thermo", "freesteam", "False")
        config.set("Thermo", "coolProp", "False")
        config.set("Thermo", "refprop", "False")
        return config


class Dialog(QtWidgets.QDialog):
    """Dialog to config thermal method calculations"""
    def __init__(self, config=None, parent=None):
        super(Dialog, self).__init__(parent)
        self.setWindowTitle(QtWidgets.QApplication.translate(
            "pychemqt", "Define project thermodynamic methods"))
        layout = QtWidgets.QVBoxLayout(self)
        self.datos = UI_confThermo_widget(config)
        layout.addWidget(self.datos)
        self.buttonBox = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Cancel |
                                                QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox)

    def value(self, config):
        """Function result for wizard"""
        config = self.datos.value(config)
        return config


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = Dialog()
    Dialog.show()
    sys.exit(app.exec_())
