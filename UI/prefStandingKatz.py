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
# Library to configure the Standing-Katz chart
#
#   - Widget: Standing-Katz chart configuration
#   - ConfigDialog: Dialog tool for standalone use
###############################################################################


import re

from PyQt5 import QtWidgets

from lib.crude import Z_list
from UI.widgets import LineConfig, Entrada_con_unidades


class Widget(QtWidgets.QWidget):
    """Standing-Katz chart configuration"""
    def __init__(self, config=None, parent=None):
        super(Widget, self).__init__(parent)
        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Method:")), 1, 1)
        self.method = QtWidgets.QComboBox()
        for Z in Z_list:
            name = Z.__name__[2:].replace("_", "-")
            year = re.search("((\d+))", Z.__doc__).group(0)
            doc = "%s (%s)" % (name, year)
            self.method.addItem(doc)
        layout.addWidget(self.method, 1, 2)

        layout.addWidget(QtWidgets.QLabel("Pr min:"), 2, 1)
        self.Prmin = Entrada_con_unidades(float, width=60, decimales=1)
        layout.addWidget(self.Prmin, 2, 2)
        layout.addWidget(QtWidgets.QLabel("Pr max:"), 3, 1)
        self.Prmax = Entrada_con_unidades(float, width=60, decimales=1)
        layout.addWidget(self.Prmax, 3, 2)

        layout.addWidget(QtWidgets.QLabel("Tr:"), 4, 1)
        self.Tr = QtWidgets.QLineEdit()
        layout.addWidget(self.Tr, 4, 2)
        self.lineconfig = LineConfig(
            "line", QtWidgets.QApplication.translate(
                "pychemqt", "Reduced temperature style line"))
        layout.addWidget(self.lineconfig, 5, 1, 1, 2)

        self.cruxconfig = LineConfig(
            "crux", QtWidgets.QApplication.translate(
                "pychemqt", "Crux style line"))
        layout.addWidget(self.cruxconfig, 6, 1, 1, 2)

        layout.addItem(QtWidgets.QSpacerItem(
            10, 0, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 10, 1, 1, 3)

        if config and config.has_section("Standing_Katz"):
            self.method.setCurrentIndex(config.getint(
                "Standing_Katz", "method"))
            self.Prmin.setValue(config.getfloat("Standing_Katz", "Prmin"))
            self.Prmax.setValue(config.getfloat("Standing_Katz", "Prmax"))
            self.Tr.setText(config.get("Standing_Katz", "Tr"))
            self.lineconfig.setConfig(config, "Standing_Katz")
            self.cruxconfig.setConfig(config, "Standing_Katz")

    def value(self, config):
        if not config.has_section("Standing_Katz"):
            config.add_section("Standing_Katz")
        config.set("Standing_Katz", "method", str(self.method.currentIndex()))
        config.set("Standing_Katz", "Prmin", str(self.Prmin.value))
        config.set("Standing_Katz", "Prmax", str(self.Prmax.value))
        config.set("Standing_Katz", "Tr", self.Tr.text())
        config = self.lineconfig.value(config, "Standing_Katz")
        config = self.cruxconfig.value(config, "Standing_Katz")
        return config


class ConfigDialog(QtWidgets.QDialog):
    """Dialog to config thermal method calculations"""
    def __init__(self, config=None, parent=None):
        super(ConfigDialog, self).__init__(parent)
        self.setWindowTitle(QtWidgets.QApplication.translate(
            "pychemqt", "Moody diagram configuration"))
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
