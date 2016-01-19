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
# PFD resolution tools
###############################################################################

from configparser import ConfigParser

from PyQt5 import QtWidgets


from UI.widgets import Entrada_con_unidades
from lib.config import conf_dir


class UI_confResolution_widget(QtWidgets.QWidget):
    """PFD resolution widget"""
    def __init__(self, config=None, parent=None):
        self.standards = [(600, 400), (640, 480), (720, 400), (800, 600),
                          (832, 624), (1024, 768), (1152, 864), (1280, 1024),
                          (1700, 1250), (1900, 1425), (2400, 1800), (4000, 3000)]
        super(UI_confResolution_widget, self).__init__(parent)
        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Use default resolution:")), 0, 0)
        self.standard = QtWidgets.QComboBox()
        self.standard.addItem("")
        for resolucion in self.standards:
            self.standard.addItem("%ix%i" % resolucion)
        self.standard.currentIndexChanged.connect(self.changeResolution)
        layout.addWidget(self.standard, 0, 1)

        self.checkCustom = QtWidgets.QCheckBox(QtWidgets.QApplication.translate(
            "pychemqt", "Use Custom resolution"))
        layout.addWidget(self.checkCustom, 1, 0, 1, 2)
        label = QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Width:"))
        label.setIndent(50)
        layout.addWidget(label, 2, 0)
        self.x = Entrada_con_unidades(int, width=60, spinbox=True, step=1)
        layout.addWidget(self.x, 2, 1)
        label = QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Height:"))
        label.setIndent(50)
        layout.addWidget(label, 3, 0)
        self.y = Entrada_con_unidades(int, width=60, spinbox=True, step=1)
        layout.addWidget(self.y, 3, 1)

        self.checkCustom.toggled.connect(self.x.setEnabled)
        self.checkCustom.toggled.connect(self.y.setEnabled)

        if config and config.has_section("PFD"):
            x = config.getint("PFD", "x")
            y = config.getint("PFD", "y")
            self.x.setValue(x)
            self.y.setValue(y)
            if (x, y) in self.standards:
                self.standard.setCurrentIndex(self.standards.index((x, y))+1)
                self.checkCustom.setChecked(False)
                self.x.setEnabled(False)
                self.y.setEnabled(False)
            else:
                self.standard.setCurrentIndex(0)
                self.checkCustom.setChecked(True)

    def changeResolution(self):
        """Change resolution with value of current opction selected"""
        x, y = self.standard.currentText().split("x")
        self.x.setValue(int(x))
        self.y.setValue(int(y))

    def value(self, config):
        """Function result to wizard"""
        if not config.has_section("PFD"):
            config.add_section("PFD")
        config.set("PFD", "x", str(self.x.value))
        config.set("PFD", "y", str(self.y.value))
        return config

    @classmethod
    def default(cls, config):
        config.add_section("PFD")
        Preferences = ConfigParser()
        Preferences.read(conf_dir+"pychemqtrc")
        config.set("PFD", "x", Preferences.get("PFD", "x"))
        config.set("PFD", "y", Preferences.get("PFD", "y"))
        return config


class Dialog(QtWidgets.QDialog):
    """PFD resolution dialog"""
    def __init__(self, config=None, parent=None):
        super(Dialog, self).__init__(parent)
        self.setWindowTitle(
            QtWidgets.QApplication.translate("pychemqt", "Define PFD resolution"))
        layout = QtWidgets.QVBoxLayout(self)
        self.datos = UI_confResolution_widget(config)
        layout.addWidget(self.datos)
        self.buttonBox = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Cancel |
                                                QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox)

    def value(self, config):
        """Function result to wizard"""
        config = self.datos.value(config)
        return config


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = Dialog()
    Dialog.show()
    sys.exit(app.exec_())
