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
# Library to configure the hypotethical pseudocomponents definition
#
#   - Widget: Petro pseudocompoent configuration
#   - ConfigDialog: Dialog tool for standalone use
###############################################################################


from PyQt5 import QtWidgets

from lib.petro import Petroleo


class Widget(QtWidgets.QWidget):
    """Petro new component configuration"""
    def __init__(self, config=None, parent=None):
        super(Widget, self).__init__(parent)

        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Molecular weight:")), 1, 1)
        Pm = ["Riazi Daubert", "Riazi Daubert extended", "Lee Kesler",
              "Sim Daubert", "API", "ASTM", "Goossens", "TWu"]
        self.Peso_molecular = QtWidgets.QComboBox()
        for p in Pm:
            self.Peso_molecular.addItem(p)
        layout.addWidget(self.Peso_molecular, 1, 2)

        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Critics properties:")), 2, 1)
        Critical = ["Riazi Daubert", "Riazi Daubert extended", "Riazi Adwani",
                    "Lee Kesler", "Cavett", "Sim Daubert",
                    "Watansiri Owens Starling", "Edmister", "Magoulas", "Twu",
                    "Tsonopoulos"]
        self.critical = QtWidgets.QComboBox()
        for c in Critical:
            self.critical.addItem(c)
        layout.addWidget(self.critical, 2, 2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Critic volume:")), 3, 1)
        vc = ["Riazi Daubert", "Riazi Daubert extended", "Riazi Adwani",
              "Watansiri Owens Starling", "Twu", "Tsonopoulos",
              "Hall Yarborough", "API"]
        self.vc = QtWidgets.QComboBox()
        for v in vc:
            self.vc.addItem(v)
        layout.addWidget(self.vc, 3, 2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Acentric factor:")), 4, 1)
        self.factor_acentrico = QtWidgets.QComboBox()
        self.factor_acentrico.addItem("Edmister")
        self.factor_acentrico.addItem("Lee Kesler")
        self.factor_acentrico.addItem("Watansiri Owens Starling")
        self.factor_acentrico.addItem("Magoulas")
        layout.addWidget(self.factor_acentrico, 4, 2)
        layout.addWidget(QtWidgets.QLabel("Z<sub>c</sub>:"), 5, 1)
        self.Zc = QtWidgets.QComboBox()
        self.Zc.addItem("Lee Kesler")
        self.Zc.addItem("Haugen")
        self.Zc.addItem("Reid")
        self.Zc.addItem("Salerno")
        self.Zc.addItem("Nath")
        layout.addWidget(self.Zc, 5, 2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "T boiling:")), 6, 1)
        self.t_ebull = QtWidgets.QComboBox()
        self.t_ebull.addItem("Riazi Daubert extended")
        self.t_ebull.addItem("Riazi Adwani")
        self.t_ebull.addItem("Edmister")
        self.t_ebull.addItem("Soreide")
        layout.addWidget(self.t_ebull, 6, 2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "PNA descomposition:")), 7, 1)
        self.PNA = QtWidgets.QComboBox()
        for method in Petroleo.METHODS_PNA:
            self.PNA.addItem(method)
        layout.addWidget(self.PNA, 7, 2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Destilate curve conversion:")), 12, 1)
        self.Curvas = QtWidgets.QComboBox()
        self.Curvas.addItem("Riazi")
        self.Curvas.addItem("Daubert")
        layout.addWidget(self.Curvas, 12, 2)

        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "% Hydrogen:")), 13, 1)
        self.Hidrogeno = QtWidgets.QComboBox()
        for method in Petroleo.METHODS_H:
            self.Hidrogeno.addItem(method)
        layout.addWidget(self.Hidrogeno, 13, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 0, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 15, 1, 1, 3)

        if config.has_section("petro"):
            self.Peso_molecular.setCurrentIndex(
                config.getint("petro", "molecular_weight"))
            self.critical.setCurrentIndex(config.getint("petro", "critical"))
            self.vc.setCurrentIndex(config.getint("petro", "vc"))
            self.factor_acentrico.setCurrentIndex(
                config.getint("petro", "f_acent"))
            self.t_ebull.setCurrentIndex(config.getint("petro", "t_ebull"))
            self.Zc.setCurrentIndex(config.getint("petro", "Zc"))
            self.PNA.setCurrentIndex(config.getint("petro", "PNA"))
            self.Hidrogeno.setCurrentIndex(config.getint("petro", "H"))
            self.Curvas.setCurrentIndex(config.getint("petro", "curva"))

    def value(self, config):
        if not config.has_section("petro"):
            config.add_section("petro")
        config.set("petro", "molecular_weight",
                   str(self.Peso_molecular.currentIndex()))
        config.set("petro", "critical", str(self.critical.currentIndex()))
        config.set("petro", "vc", str(self.vc.currentIndex()))
        config.set("petro", "f_acent",
                   str(self.factor_acentrico.currentIndex()))
        config.set("petro", "t_ebull", str(self.t_ebull.currentIndex()))
        config.set("petro", "Zc", str(self.Zc.currentIndex()))
        config.set("petro", "PNA", str(self.PNA.currentIndex()))
        config.set("petro", "H", str(self.Hidrogeno.currentIndex()))
        config.set("petro", "curva", str(self.Curvas.currentIndex()))
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


if __name__ == "__main__":
    import os
    import sys
    from configparser import ConfigParser
    app = QtWidgets.QApplication(sys.argv)

    conf_dir = os.path.expanduser('~') + "/.pychemqt/"
    config = ConfigParser()
    config.read(conf_dir+"pychemqtrc")

    Dialog = ConfigDialog(config)
    Dialog.show()
    sys.exit(app.exec_())
