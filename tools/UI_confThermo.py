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
# Thermal method calculation config
###############################################################################


from ast import literal_eval
import os

from tools.qt import QtWidgets

from lib.EoS import K, H, alfa, mix, cp_ideal, K_name, H_name, K_status
from lib.bip import EoSBIP
from lib.config import getMainWindowConfig
from lib.corriente import Corriente
from UI import BIP


class UI_confThermo_widget(QtWidgets.QWidget):
    """Widget to config thermal method calculations"""
    def __init__(self, config=None, parent=None):
        """Constructor, opcional config parameter with proyect configuration"""
        super().__init__(parent)
        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(QtWidgets.QLabel(self.tr("K values:")), 0, 0, 1, 2)
        self.K = QtWidgets.QComboBox()
        for eq in K:
            self.K.addItem(eq.__title__)
        self.K.currentIndexChanged.connect(self.updateBIP)
        layout.addWidget(self.K, 0, 2)
        self.bipButton = QtWidgets.QPushButton(self.tr("BIP"))
        self.bipButton.clicked.connect(self.showBIP)
        layout.addWidget(self.bipButton, 0, 3)

        layout.addWidget(QtWidgets.QLabel(
            self.tr("Alfa function:")), 1, 0, 1, 2)
        self.alfa = QtWidgets.QComboBox()
        for a in alfa:
            self.alfa.addItem(a)
        layout.addWidget(self.alfa, 1, 2, 1, 2)
        layout.addWidget(QtWidgets.QLabel(self.tr("Mix rules:")), 2, 0, 1, 2)
        self.mixing_rule = QtWidgets.QComboBox()
        for m in mix:
            self.mixing_rule.addItem(m)
        layout.addWidget(self.mixing_rule, 2, 2, 1, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 3, 0, 1, 4)
        layout.addWidget(QtWidgets.QLabel(self.tr("Enthalpy:")), 4, 0, 1, 2)
        self.H = QtWidgets.QComboBox()
        for h in H:
            self.H.addItem(h.__title__)
        layout.addWidget(self.H, 4, 2, 1, 2)
        layout.addWidget(QtWidgets.QLabel(
            self.tr("Ideal heat capacity:")), 5, 0, 1, 2)
        self.Cp_ideal = QtWidgets.QComboBox()
        for cp in cp_ideal:
            self.Cp_ideal.addItem(cp)
        layout.addWidget(self.Cp_ideal, 5, 2, 1, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 6, 0, 1, 4)
        self.MEoS = QtWidgets.QCheckBox(
            self.tr("Use MEoS for single compounds if it's available"))
        layout.addWidget(self.MEoS, 7, 0, 1, 4)
        self.coolProp = QtWidgets.QCheckBox(
            self.tr("Use external library coolProp (faster)"))
        self.coolProp.setEnabled(False)
        layout.addWidget(self.coolProp, 8, 1, 1, 3)
        self.refprop = QtWidgets.QCheckBox(
            self.tr("Use external library refprop (fastest)"))
        self.refprop.setEnabled(False)
        layout.addWidget(self.refprop, 9, 1, 1, 3)

        self.iapws = QtWidgets.QCheckBox(self.tr("Use IAPWS97 for water"))
        layout.addWidget(self.iapws, 10, 0, 1, 4)
        self.freesteam = QtWidgets.QCheckBox(
            self.tr("Use freesteam library (faster)"))
        self.freesteam.setEnabled(False)
        layout.addWidget(self.freesteam, 11, 1, 1, 3)
        self.GERG = QtWidgets.QCheckBox(
            self.tr("Use GERG EoS for mix if it's posible"))
        layout.addWidget(self.GERG, 12, 0, 1, 4)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 13, 0, 1, 5)

        if os.environ["freesteam"] == "True":
            self.iapws.toggled.connect(self.freesteam.setEnabled)
        if os.environ["CoolProp"] == "True":
            self.MEoS.toggled.connect(self.coolProp.setEnabled)
        if os.environ["refprop"] == "True":
            self.MEoS.toggled.connect(self.refprop.setEnabled)

        if config:
            self.setConfig(config)
        self.updateBIP(0)

    def setConfig(self, config):
        """Set values from currentConfig"""
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

    def setKwargs(self, kwarg):
        """Set values from kwargs dict"""
        config = getMainWindowConfig()
        self.setConfig(config)
        for key in ["MEoS", "iapws", "GERG", "freesteam", "coolProp",
                    "refprop"]:
            if kwarg[key] != Corriente.kwargs[key]:
                getattr(self, key).setChecked(kwarg[key])

        if kwarg["K"] != Corriente.kwargs["K"]:
            index = K_name.index(kwarg["K"])
            self.K.setCurrentIndex(index)
        if kwarg["H"] != Corriente.kwargs["H"]:
            index = H_name.index(kwarg["H"])
            self.H.setCurrentIndex(index)
        if kwarg["mix"] != Corriente.kwargs["mix"]:
            index = mix.index(kwarg["mix"])
            self.mixing_rule.setCurrentIndex(index)
        if kwarg["Cp_ideal"] != Corriente.kwargs["Cp_ideal"]:
            self.Cp_ideal.setCurrentIndex(kwarg["Cp_ideal"])
        if kwarg["alfa"] != Corriente.kwargs["alfa"]:
            try:
                index = alfa.index(kwarg["alfa"])
                self.alfa.setCurrentIndex(index)
            except ValueError:
                pass

    @property
    def kwargs(self):
        """Return dict with kwarg"""
        kw = {}
        kw["K"] = self.K.currentText().split(" (")[0]
        kw["alfa"] = self.alfa.currentText()
        kw["mix"] = self.mixing_rule.currentText()
        kw["H"] = self.H.currentText().split(" (")[0]
        kw["Cp_ideal"] = self.Cp_ideal.currentIndex()
        kw["MEoS"] = self.MEoS.isChecked()
        kw["iapws"] = self.iapws.isChecked()
        kw["GERG"] = self.GERG.isChecked()
        kw["freesteam"] = self.freesteam.isChecked()
        kw["coolProp"] = self.coolProp.isChecked()
        kw["refprop"] = self.refprop.isChecked()
        return kw

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
        """Default configuration"""
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

    def updateBIP(self, index):
        """Update bipButton enabled status only for EoS with bip available"""
        if K_status[index] in EoSBIP:
            self.bipButton.setEnabled(True)
        else:
            self.bipButton.setEnabled(False)

    def showBIP(self):
        """Show BIP dialog from component for selected EoS"""
        config = getMainWindowConfig()
        ids = literal_eval(config.get("Components", "components"))
        index = EoSBIP.index(K_status[self.K.currentIndex()])
        dlg = BIP.Ui_BIP(ids, index)
        dlg.exec()


class Dialog(QtWidgets.QDialog):
    """Dialog to config thermal method calculations"""
    def __init__(self, config=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle(self.tr("Define project thermodynamic methods"))
        layout = QtWidgets.QVBoxLayout(self)
        self.datos = UI_confThermo_widget(config)
        layout.addWidget(self.datos)
        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Cancel
            | QtWidgets.QDialogButtonBox.StandardButton.Ok)
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
    sys.exit(app.exec())
