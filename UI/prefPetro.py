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
along with this program.  If not, see <http://www.gnu.org/licenses/>.


Library to configure the hypotethical pseudocomponents definition

* :class:`Widget`: Petro pseudocompoent configuration
* :class:`ConfigDialog`: Dialog tool for standalone use
'''


from tools.qt import QtWidgets

from lib.petro import Petroleo


class Widget(QtWidgets.QWidget):
    """Petro new component configuration"""
    def __init__(self, config=None, parent=None):
        super().__init__(parent)

        lyt = QtWidgets.QGridLayout(self)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Molecular weight")), 1, 1)
        self.M = QtWidgets.QComboBox()
        for p in Petroleo.METHODS_M:
            self.M.addItem(p)
        lyt.addWidget(self.M, 1, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Critic properties")), 2, 1)
        self.critical = QtWidgets.QComboBox()
        for c in Petroleo.METHODS_crit:
            self.critical.addItem(c)
        lyt.addWidget(self.critical, 2, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Critic volume")), 3, 1)
        self.vc = QtWidgets.QComboBox()
        for v in Petroleo.METHODS_Vc:
            self.vc.addItem(v)
        lyt.addWidget(self.vc, 3, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Acentric factor")), 4, 1)
        self.f_acent = QtWidgets.QComboBox()
        for w in Petroleo.METHODS_w:
            self.f_acent.addItem(w)
        lyt.addWidget(self.f_acent, 4, 2)
        lyt.addWidget(QtWidgets.QLabel("Z<sub>c</sub>"), 5, 1)
        self.Zc = QtWidgets.QComboBox()
        for method in Petroleo.METHODS_Zc:
            self.Zc.addItem(method)
        lyt.addWidget(self.Zc, 5, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Boiling Temperature")), 6, 1)
        self.Tb = QtWidgets.QComboBox()
        for tb in Petroleo.METHODS_Tb:
            self.Tb.addItem(tb)
        lyt.addWidget(self.Tb, 6, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Specific Gravity")), 7, 1)
        self.SG = QtWidgets.QComboBox()
        for sg in Petroleo.METHODS_SG:
            self.SG.addItem(sg)
        lyt.addWidget(self.SG, 7, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Refractive Index")), 8, 1)
        self.n = QtWidgets.QComboBox()
        for n in Petroleo.METHODS_n:
            self.n.addItem(n)
        lyt.addWidget(self.n, 8, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("PNA composition")), 9, 1)
        self.PNA = QtWidgets.QComboBox()
        for method in Petroleo.METHODS_PNA:
            self.PNA.addItem(method)
        lyt.addWidget(self.PNA, 9, 2)
        lyt.addWidget(QtWidgets.QLabel(
            self.tr("Destilate curve conversion")), 10, 1)
        self.curves = QtWidgets.QComboBox()
        self.curves.addItem("Riazi")
        self.curves.addItem("Daubert")
        lyt.addWidget(self.curves, 10, 2)

        lyt.addWidget(QtWidgets.QLabel(self.tr("Hydrogen %")), 11, 1)
        self.H = QtWidgets.QComboBox()
        for method in Petroleo.METHODS_H:
            self.H.addItem(method)
        lyt.addWidget(self.H, 11, 2)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 15, 1, 1, 3)

        if config.has_section("petro"):
            self.M.setCurrentIndex(
                config.getint("petro", "M"))
            self.critical.setCurrentIndex(config.getint("petro", "critical"))
            self.vc.setCurrentIndex(config.getint("petro", "vc"))
            self.f_acent.setCurrentIndex(
                config.getint("petro", "f_acent"))
            self.Tb.setCurrentIndex(config.getint("petro", "Tb"))
            self.SG.setCurrentIndex(config.getint("petro", "SG"))
            self.n.setCurrentIndex(config.getint("petro", "n"))
            self.Zc.setCurrentIndex(config.getint("petro", "Zc"))
            self.PNA.setCurrentIndex(config.getint("petro", "PNA"))
            self.H.setCurrentIndex(config.getint("petro", "H"))
            self.curves.setCurrentIndex(config.getint("petro", "curve"))

    def value(self, config):
        """Update configparser instance with the preferences defined in this
        widget"""
        if not config.has_section("petro"):
            config.add_section("petro")
        config.set("petro", "M",
                   str(self.M.currentIndex()))
        config.set("petro", "critical", str(self.critical.currentIndex()))
        config.set("petro", "vc", str(self.vc.currentIndex()))
        config.set("petro", "f_acent",
                   str(self.f_acent.currentIndex()))
        config.set("petro", "Tb", str(self.Tb.currentIndex()))
        config.set("petro", "SG", str(self.SG.currentIndex()))
        config.set("petro", "n", str(self.n.currentIndex()))
        config.set("petro", "Zc", str(self.Zc.currentIndex()))
        config.set("petro", "PNA", str(self.PNA.currentIndex()))
        config.set("petro", "H", str(self.H.currentIndex()))
        config.set("petro", "curve", str(self.curves.currentIndex()))
        return config


class ConfigDialog(QtWidgets.QDialog):
    """Dialog to config thermal method calculations"""
    def __init__(self, config=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle(self.tr("Petrol assay definition configuration"))
        lyt = QtWidgets.QVBoxLayout(self)
        self.widget = Widget(config)
        lyt.addWidget(self.widget)
        buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Cancel
            | QtWidgets.QDialogButtonBox.StandardButton.Ok)
        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)
        lyt.addWidget(buttonBox)

    def value(self, config):
        """Function result for wizard"""
        config = self.widget.value(config)
        return config
