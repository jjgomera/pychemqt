#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2022, Juan José Gómez Romera <jjgomera@gmail.com>

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
# Library to configure qtelemental tool
#
#   - Widget: Moody chart configuration
#   - Dialog: Dialog tool for standalone use
###############################################################################


from qt import QtWidgets


class Widget(QtWidgets.QWidget):
    """External applications configuration"""
    def __init__(self, config=None, parent=None):
        super(Widget, self).__init__(parent)
        layout = QtWidgets.QGridLayout(self)

        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Color by element:")),
            1, 1)
        colorby = ["Element", "serie", "group_element", "period", "block", "phase",
                   "lattice_type", "space_group", "density_Solid",
                   "density_Liq", "density_Gas", "date", "atomic_mass",
                   "atomic_volume", "atomic_radius", "covalent_radius",
                   "vanderWaals_radius", "electronegativity",
                   "electron_affinity", "first_ionization", "Tf", "Tb",
                   "Heat_f", "Heat_b", "Cp", "k", "T_debye"]
        self.ElementalColorby = QtWidgets.QComboBox()
        for c in colorby:
            self.ElementalColorby.addItem(c)
        layout.addWidget(self.ElementalColorby, 1, 2)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Color definition")),
            2, 1)
        self.ElementalDefinition = QtWidgets.QSpinBox()
        self.ElementalDefinition.setMaximumWidth(50)
        self.ElementalDefinition.setMinimum(5)
        layout.addWidget(self.ElementalDefinition, 2, 2)
        self.ElementalLog = QtWidgets.QCheckBox(
            QtWidgets.QApplication.translate("pychemqt", "Logarithmic scale"))
        layout.addWidget(self.ElementalLog, 3, 1, 1, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 4, 3)

        if config.has_section("Applications"):
            self.ElementalColorby.setCurrentText(
                config.get("Applications", 'elementalColorby'))
            self.ElementalDefinition.setValue(
                config.getint("Applications", 'elementalDefinition'))
            self.ElementalLog.setChecked(
                config.getboolean("Applications", 'elementalLog'))

    def value(self, config):
        """Return value for main dialog"""
        if not config.has_section("Applications"):
            config.add_section("Applications")
        config.set("Applications", "elementalColorby",
                   self.ElementalColorby.currentText())
        config.set("Applications", "elementalDefinition",
                   str(self.ElementalDefinition.value()))
        config.set("Applications", "elementalLog",
                   str(self.ElementalLog.isChecked()))
        return config


class Dialog(QtWidgets.QDialog):
    """Dialog to config thermal method calculations"""
    def __init__(self, config=None, parent=None):
        super(Dialog, self).__init__(parent)
        self.setWindowTitle(QtWidgets.QApplication.translate(
            "pychemqt", "Qtelemental configuration"))
        layout = QtWidgets.QVBoxLayout(self)
        self.widget = Widget(config)
        layout.addWidget(self.widget)
        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Cancel | QtWidgets.QDialogButtonBox.StandardButton.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox)

    def value(self, config):
        """Function result for wizard"""
        config = self.widget.value(config)
        return config
