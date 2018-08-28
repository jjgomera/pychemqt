#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2017, Juan José Gómez Romera <jjgomera@gmail.com>

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
# Library to configure the psychrometric chart
#
#   - Widget: Psychrometric chart configuration
#   - Dialog: Dialog tool for standalone use
###############################################################################

import os

from PyQt5 import QtWidgets

from lib import unidades
from UI.widgets import LineConfig
from UI.prefMEOS import Isolinea


class Widget(QtWidgets.QWidget):
    """Phychrometric chart configuration"""
    lineas = [
        ("IsoTdb", unidades.Temperature,
         QtWidgets.QApplication.translate(
             "pychemqt", "Iso dry bulb temperature")),
        ("IsoW", float,
         QtWidgets.QApplication.translate(
             "pychemqt", "Iso absolute humidity")),
        ("IsoHR", float,
         QtWidgets.QApplication.translate(
             "pychemqt", "Iso relative humidity")),
        ("IsoTwb", unidades.Temperature,
         QtWidgets.QApplication.translate(
             "pychemqt", "Iso wet bulb temperature")),
        ("Isochor", unidades.SpecificVolume,
         QtWidgets.QApplication.translate("pychemqt", "Isochor"))]

    def __init__(self, config, parent=None):
        """constructor, config optional parameter to input project config"""
        super(Widget, self).__init__(parent)
        lyt = QtWidgets.QGridLayout(self)
        lyt.setContentsMargins(0, 0, 0, 0)
        scroll = QtWidgets.QScrollArea()
        scroll.setFrameStyle(QtWidgets.QFrame.NoFrame)
        lyt.addWidget(scroll)

        dlg = QtWidgets.QWidget()
        layout = QtWidgets.QGridLayout(dlg)

        groupType = QtWidgets.QGroupBox(QtWidgets.QApplication.translate(
            "pychemqt", "Chart type"))
        groupLayout = QtWidgets.QVBoxLayout(groupType)
        self.checkASHRAE = QtWidgets.QRadioButton(
            QtWidgets.QApplication.translate(
                "pychemqt", "ASHRAE Chart, W vs Tdb"))
        groupLayout.addWidget(self.checkASHRAE)
        self.checkMollier = QtWidgets.QRadioButton(
            QtWidgets.QApplication.translate("pychemqt", "Mollier Chart ix"))
        groupLayout.addWidget(self.checkMollier)
        layout.addWidget(groupType, 0, 1, 1, 2)

        self.virial = QtWidgets.QCheckBox(QtWidgets.QApplication.translate(
            "pychemqt", "Use virial equation of state"))
        layout.addWidget(self.virial, 1, 1, 1, 2)
        self.coolProp = QtWidgets.QCheckBox(QtWidgets.QApplication.translate(
            "pychemqt", "Use external library coolProp (faster)"))
        self.coolProp.setEnabled(False)
        layout.addWidget(self.coolProp, 2, 2)
        self.refprop = QtWidgets.QCheckBox(QtWidgets.QApplication.translate(
            "pychemqt", "Use external library refprop (fastest)"))
        self.refprop.setEnabled(False)
        layout.addWidget(self.refprop, 3, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Fixed,
            QtWidgets.QSizePolicy.Fixed), 4, 1)

        self.satlineconfig = LineConfig(
            "saturation", QtWidgets.QApplication.translate(
                "pychemqt", "Saturation Line Style"))
        layout.addWidget(self.satlineconfig, 5, 1, 1, 2)
        self.cruxlineconfig = LineConfig(
            "crux", QtWidgets.QApplication.translate(
                "pychemqt", "Crux Line Style"))
        layout.addWidget(self.cruxlineconfig, 6, 1, 1, 2)
        group = QtWidgets.QGroupBox(
            QtWidgets.QApplication.translate("pychemqt", "Isolines"))
        layout.addWidget(group, 7, 1, 1, 2)
        layoutgroup = QtWidgets.QGridLayout(group)
        self.comboIsolineas = QtWidgets.QComboBox()
        layoutgroup.addWidget(self.comboIsolineas, 1, 1)
        self.Isolineas = QtWidgets.QStackedWidget()
        self.comboIsolineas.currentIndexChanged.connect(
            self.Isolineas.setCurrentIndex)
        layoutgroup.addWidget(self.Isolineas, 2, 1)
        for name, unit, text in self.lineas:
            self.comboIsolineas.addItem(text)
            self.Isolineas.addWidget(Isolinea(unit, name, config, "Psychr"))
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 10, 2)

        scroll.setWidget(dlg)

        if os.environ["CoolProp"]:
            self.virial.toggled.connect(self.coolProp.setEnabled)
        if os.environ["refprop"]:
            self.virial.toggled.connect(self.refprop.setEnabled)

        if config.has_section("Psychr"):
            if config.getboolean("Psychr", 'chart'):
                self.checkASHRAE.setChecked(True)
            else:
                self.checkMollier.setChecked(True)
            self.virial.setChecked(config.getboolean("Psychr", 'virial'))
            self.coolProp.setChecked(config.getboolean("Psychr", 'coolprop'))
            self.refprop.setChecked(config.getboolean("Psychr", 'refprop'))
            self.satlineconfig.setConfig(config, "Psychr")
            self.cruxlineconfig.setConfig(config, "Psychr")

    def value(self, config):
        """Return value for main dialog"""
        if not config.has_section("Psychr"):
            config.add_section("Psychr")

        config.set("Psychr", "chart", str(self.checkASHRAE.isChecked()))
        config.set("Psychr", "virial", str(self.virial.isChecked()))
        config.set("Psychr", "coolprop", str(self.coolProp.isChecked()))
        config.set("Psychr", "refprop", str(self.refprop.isChecked()))
        config = self.satlineconfig.value(config, "Psychr")
        config = self.cruxlineconfig.value(config, "Psychr")

        for indice in range(self.Isolineas.count()):
            config = self.Isolineas.widget(indice).value(config)
        return config


class Dialog(QtWidgets.QDialog):
    """Dialog to config thermal method calculations"""
    def __init__(self, config=None, parent=None):
        super(Dialog, self).__init__(parent)
        self.setWindowTitle(QtWidgets.QApplication.translate(
            "pychemqt", "Psychrometric chart configuration"))
        layout = QtWidgets.QVBoxLayout(self)
        self.widget = Widget(config)
        layout.addWidget(self.widget)
        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Cancel | QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox)

    def value(self, config):
        """Function result for wizard"""
        config = self.widget.value(config)
        return config


if __name__ == "__main__":
    import sys
    from configparser import ConfigParser
    app = QtWidgets.QApplication(sys.argv)

    conf_dir = os.path.expanduser('~') + "/.pychemqt/"
    config = ConfigParser()
    config.read(conf_dir+"pychemqtrc")

    Dialog = Dialog(config)
    Dialog.show()
    sys.exit(app.exec_())
