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
# Transport properties config section
###############################################################################


from qt import QtCore, QtWidgets

from lib.compuestos import Componente
from lib.mezcla import Mezcla


class UI_confTransport_widget(QtWidgets.QWidget):
    """Transport properties widget, tu use in dialog, wizard..."""
    def __init__(self, config=None, parent=None):
        """Constructor, opcional config parameter with project config"""
        super(UI_confTransport_widget, self).__init__(parent)
        layout = QtWidgets.QGridLayout(self)

        line = QtWidgets.QFrame()
        line.setFrameShape(QtWidgets.QFrame.Shape.HLine)
        line.setFrameShadow(QtWidgets.QFrame.Shadow.Sunken)
        layout.addWidget(line, 1, 1, 1, 7)
        line = QtWidgets.QFrame()
        line.setFrameShape(QtWidgets.QFrame.Shape.HLine)
        line.setFrameShadow(QtWidgets.QFrame.Shadow.Sunken)
        layout.addWidget(line, 3, 1, 1, 7)
        lbl_Pure = QtWidgets.QLabel()
        lbl_Pure.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        lbl_Pure.setText(QtWidgets.QApplication.translate(
            "pychemqt", "Pure Fluid Correlations"))
        layout.addWidget(lbl_Pure, 0, 1, 1, 3)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Liquid Density:")), 4, 0)
        self.RhoL = QtWidgets.QComboBox()
        for method in Componente.METHODS_RhoL:
            self.RhoL.addItem(method)
        layout.addWidget(self.RhoL, 4, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Liquid Viscosity:")), 5, 0)
        self.MuL = QtWidgets.QComboBox()
        for method in Componente.METHODS_MuL:
            self.MuL.addItem(method)
        layout.addWidget(self.MuL, 5, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Gas Viscosity:")), 6, 0)
        self.MuG = QtWidgets.QComboBox()
        for method in Componente.METHODS_MuG:
            self.MuG.addItem(method)
        layout.addWidget(self.MuG, 6, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Surface Tension:")), 7, 0)
        self.Tension = QtWidgets.QComboBox()
        for method in Componente.METHODS_Tension:
            self.Tension.addItem(method)
        layout.addWidget(self.Tension, 7, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Liquid Thermal Conductivity:")), 8, 0)
        self.ThCondL = QtWidgets.QComboBox()
        for method in Componente.METHODS_ThL:
            self.ThCondL.addItem(method)
        layout.addWidget(self.ThCondL, 8, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Gas Thermal Conductivity:")), 9, 0)
        self.ThCondG = QtWidgets.QComboBox()
        for method in Componente.METHODS_ThG:
            self.ThCondG.addItem(method)
        layout.addWidget(self.ThCondG, 9, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Vapor Pressure:")), 10, 0)
        self.Pv = QtWidgets.QComboBox()
        for method in Componente.METHODS_Pv:
            self.Pv.addItem(method)
        layout.addWidget(self.Pv, 10, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Acentric factor:")), 11, 0)
        self.w = QtWidgets.QComboBox()
        for method in Componente.METHODS_facent:
            self.w.addItem(method)
        layout.addWidget(self.w, 11, 1)

        line = QtWidgets.QFrame()
        line.setFrameShape(QtWidgets.QFrame.Shape.VLine)
        line.setFrameShadow(QtWidgets.QFrame.Shadow.Sunken)
        layout.addWidget(line, 1, 2, 11, 1)
        lbl_hP = QtWidgets.QLabel()
        lbl_hP.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        lbl_hP.setText(QtWidgets.QApplication.translate(
            "pychemqt", "High Pressure Corrections"))
        layout.addWidget(lbl_hP, 2, 3)
        self.Corr_RhoL = QtWidgets.QComboBox()
        for method in Componente.METHODS_RhoLP:
            self.Corr_RhoL.addItem(method)
        layout.addWidget(self.Corr_RhoL, 4, 3)
        self.Corr_MuL = QtWidgets.QComboBox()
        for method in Componente.METHODS_MuLP:
            self.Corr_MuL.addItem(method)
        layout.addWidget(self.Corr_MuL, 5, 3)
        self.Corr_MuG = QtWidgets.QComboBox()
        for method in Componente.METHODS_MuGP:
            self.Corr_MuG.addItem(method)
        layout.addWidget(self.Corr_MuG, 6, 3)
        self.Corr_ThCondL = QtWidgets.QComboBox()
        for method in Componente.METHODS_ThLP:
            self.Corr_ThCondL.addItem(method)
        layout.addWidget(self.Corr_ThCondL, 8, 3)
        self.Corr_ThCondG = QtWidgets.QComboBox()
        for method in Componente.METHODS_ThGP:
            self.Corr_ThCondG.addItem(method)
        layout.addWidget(self.Corr_ThCondG, 9, 3)

        line = QtWidgets.QFrame()
        line.setFrameShape(QtWidgets.QFrame.Shape.VLine)
        line.setFrameShadow(QtWidgets.QFrame.Shadow.Sunken)
        layout.addWidget(line, 0, 4, 12, 1)
        lbl_Mix = QtWidgets.QLabel()
        lbl_Mix.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        lbl_Mix.setText(QtWidgets.QApplication.translate(
            "pychemqt", "Mixture Fluid Correlations"))
        layout.addWidget(lbl_Mix, 0, 5, 1, 3)
        self.RhoLMix = QtWidgets.QComboBox()
        for method in Mezcla.METHODS_RhoL:
            self.RhoLMix.addItem(method)
        layout.addWidget(self.RhoLMix, 4, 5)
        self.MuLMix = QtWidgets.QComboBox()
        for method in Mezcla.METHODS_MuL:
            self.MuLMix.addItem(method)
        layout.addWidget(self.MuLMix, 5, 5)
        self.MuGMix = QtWidgets.QComboBox()
        for method in Mezcla.METHODS_MuG:
            self.MuGMix.addItem(method)
        layout.addWidget(self.MuGMix, 6, 5)
        self.ThLMix = QtWidgets.QComboBox()
        for method in Mezcla.METHODS_ThL:
            self.ThLMix.addItem(method)
        layout.addWidget(self.ThLMix, 8, 5)
        self.ThGMix = QtWidgets.QComboBox()
        for method in Mezcla.METHODS_ThG:
            self.ThGMix.addItem(method)
        layout.addWidget(self.ThGMix, 9, 5)

        line = QtWidgets.QFrame()
        line.setFrameShape(QtWidgets.QFrame.Shape.VLine)
        line.setFrameShadow(QtWidgets.QFrame.Shadow.Sunken)
        layout.addWidget(line, 1, 6, 11, 1)
        lbl_hPMix = QtWidgets.QLabel()
        lbl_hPMix.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        lbl_hPMix.setText(QtWidgets.QApplication.translate(
            "pychemqt", "High Pressure Corrections"))
        layout.addWidget(lbl_hPMix, 2, 7)
        self.Corr_RhoLMix = QtWidgets.QComboBox()
        for method in Mezcla.METHODS_RhoLP:
            self.Corr_RhoLMix.addItem(method)
        layout.addWidget(self.Corr_RhoLMix, 4, 7)
        self.Corr_MuGMix = QtWidgets.QComboBox()
        for method in Mezcla.METHODS_MuGP:
            self.Corr_MuGMix.addItem(method)
        layout.addWidget(self.Corr_MuGMix, 6, 7)
        self.Corr_ThGMix = QtWidgets.QComboBox()
        for method in Mezcla.METHODS_ThGP:
            self.Corr_ThGMix.addItem(method)
        layout.addWidget(self.Corr_ThGMix, 9, 7)

        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 12, 0)
        self.rhoLEoS = QtWidgets.QCheckBox(QtWidgets.QApplication.translate(
            "pychemqt", "Use liquid density from EoS if available"))
        layout.addWidget(self.rhoLEoS, 13, 0, 1, 8)

        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 15, 8)

        if config and config.has_section("Transport"):
            self.RhoL.setCurrentIndex(config.getint("Transport", "RhoL"))
            self.Corr_RhoL.setCurrentIndex(
                config.getint("Transport", "Corr_RhoL"))
            self.MuL.setCurrentIndex(config.getint("Transport", "MuL"))
            self.Corr_MuL.setCurrentIndex(
                config.getint("Transport", "Corr_MuL"))
            self.Corr_MuG.setCurrentIndex(
                config.getint("Transport", "Corr_MuG"))
            self.MuG.setCurrentIndex(config.getint("Transport", "MuG"))
            self.Tension.setCurrentIndex(config.getint("Transport", "Tension"))
            self.ThCondL.setCurrentIndex(config.getint("Transport", "ThCondL"))
            self.Corr_ThCondL.setCurrentIndex(
                config.getint("Transport", "Corr_ThCondL"))
            self.Corr_ThCondG.setCurrentIndex(
                config.getint("Transport", "Corr_ThCondG"))
            self.ThCondG.setCurrentIndex(config.getint("Transport", "ThCondG"))
            self.Pv.setCurrentIndex(config.getint("Transport", "Pv"))
            self.w.setCurrentIndex(config.getint("Transport", "f_acent"))

            self.RhoLMix.setCurrentIndex(config.getint("Transport", "RhoLMix"))
            self.Corr_RhoLMix.setCurrentIndex(
                config.getint("Transport", "Corr_RhoLMix"))
            self.MuLMix.setCurrentIndex(config.getint("Transport", "MuLMix"))
            self.MuGMix.setCurrentIndex(config.getint("Transport", "MuGMix"))
            self.Corr_MuGMix.setCurrentIndex(
                config.getint("Transport", "Corr_MuGMix"))
            self.ThLMix.setCurrentIndex(
                config.getint("Transport", "ThCondLMix"))
            self.ThGMix.setCurrentIndex(
                config.getint("Transport", "ThCondGMix"))
            self.Corr_ThGMix.setCurrentIndex(
                config.getint("Transport", "Corr_ThCondGMix"))

            self.rhoLEoS.setChecked(config.getboolean("Transport", "RhoLEoS"))

    def value(self, config):
        """Function to wizard result"""
        if not config.has_section("Transport"):
            config.add_section("Transport")
        config.set("Transport", "RhoL", str(self.RhoL.currentIndex()))
        config.set("Transport", "Corr_RhoL",
                   str(self.Corr_RhoL.currentIndex()))
        config.set("Transport", "MuL", str(self.MuL.currentIndex()))
        config.set("Transport", "Corr_MuL", str(self.Corr_MuL.currentIndex()))
        config.set("Transport", "MuG", str(self.MuG.currentIndex()))
        config.set("Transport", "Corr_MuG", str(self.Corr_MuG.currentIndex()))
        config.set("Transport", "Tension", str(self.Tension.currentIndex()))
        config.set("Transport", "ThCondL", str(self.ThCondL.currentIndex()))
        config.set("Transport", "Corr_ThCondL",
                   str(self.Corr_ThCondL.currentIndex()))
        config.set("Transport", "ThCondG", str(self.ThCondG.currentIndex()))
        config.set("Transport", "Corr_ThCondG",
                   str(self.Corr_ThCondG.currentIndex()))
        config.set("Transport", "Pv", str(self.Pv.currentIndex()))
        config.set("Transport", "f_acent", str(self.w.currentIndex()))

        config.set("Transport", "RhoLMix", str(self.RhoLMix.currentIndex()))
        config.set("Transport", "Corr_RhoLMix",
                   str(self.Corr_RhoLMix.currentIndex()))
        config.set("Transport", "MuGMix", str(self.MuGMix.currentIndex()))
        config.set("Transport", "Corr_MuGMix",
                   str(self.Corr_MuGMix.currentIndex()))
        config.set("Transport", "MuLMix", str(self.MuLMix.currentIndex()))
        config.set("Transport", "ThCondLMix", str(self.ThLMix.currentIndex()))
        config.set("Transport", "ThCondGMix", str(self.ThGMix.currentIndex()))
        config.set("Transport", "Corr_ThCondGMix",
                   str(self.Corr_ThGMix.currentIndex()))
        config.set("Transport", "RhoLEoS", str(self.rhoLEoS.isChecked()))
        return config

    @classmethod
    def default(cls, config):
        config.add_section("Transport")
        config.set("Transport", "RhoL", "0")
        config.set("Transport", "Corr_RhoL", "0")
        config.set("Transport", "MuL", "0")
        config.set("Transport", "Corr_MuL", "0")
        config.set("Transport", "MuG", "0")
        config.set("Transport", "Tension", "0")
        config.set("Transport", "ThCondL", "0")
        config.set("Transport", "Corr_ThCondL", "0")
        config.set("Transport", "ThCondG", "0")
        config.set("Transport", "Pv", "0")
        config.set("Transport", "f_acent", "0")
        config.set("Transport", "RhoLMix", "0")
        config.set("Transport", "Corr_RhoLMix", "0")
        config.set("Transport", "MuLMix", "0")
        config.set("Transport", "MuGMix", "0")
        config.set("Transport", "Corr_MuGMix", "0")
        config.set("Transport", "ThCondLMix", "0")
        config.set("Transport", "ThCondGMix", "0")
        config.set("Transport", "Corr_ThCondGMix", "0")
        config.set("Transport", "RhoLEoS", "False")
        return config


class Dialog(QtWidgets.QDialog):
    """Transport properties dialog"""
    def __init__(self, config=None, parent=None):
        super(Dialog, self).__init__(parent)
        self.setWindowTitle(QtWidgets.QApplication.translate(
            "pychemqt", "Transport Properties Methods"))
        layout = QtWidgets.QVBoxLayout(self)
        self.datos = UI_confTransport_widget(config)
        layout.addWidget(self.datos)
        btnBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Cancel | QtWidgets.QDialogButtonBox.StandardButton.Ok)
        btnBox.accepted.connect(self.accept)
        btnBox.rejected.connect(self.reject)
        layout.addWidget(btnBox)

    def value(self, config):
        """Function to wizard result"""
        config = self.datos.value(config)
        return config


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = Dialog()
    Dialog.show()
    sys.exit(app.exec())
