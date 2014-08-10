#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Transport properties config section
###############################################################################

from PyQt4 import QtCore, QtGui


class UI_confTransport_widget(QtGui.QWidget):
    """Transport properties widget, tu use in dialog, wizard..."""
    def __init__(self, config=None, parent=None):
        """Constructor, opcional config parameter with project config"""
        super(UI_confTransport_widget, self).__init__(parent)
        layout = QtGui.QGridLayout(self)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Liquid Density:")), 1, 0)
        self.RhoL = QtGui.QComboBox()
        self.RhoL.addItem("DIPPR")
        self.RhoL.addItem("Rackett")
        self.RhoL.addItem("Cavett")
        self.RhoL.addItem("Costald")
        layout.addWidget(self.RhoL, 1, 1)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Liquid Viscosity:")), 2, 0)
        self.MuL = QtGui.QComboBox()
        self.MuL.addItem("DIPPR")
        self.MuL.addItem("Parametric")
        self.MuL.addItem("Letsou & Steil")
        layout.addWidget(self.MuL, 2, 1)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Gas Viscosity:")), 3, 0)
        self.MuG = QtGui.QComboBox()
        self.MuG.addItem("DIPPR")
        self.MuG.addItem("Chapman & Enskog")
        self.MuG.addItem("Thodos")
        layout.addWidget(self.MuG, 3, 1)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Surface Tension:")), 4, 0)
        self.Tension = QtGui.QComboBox()
        self.Tension.addItem("DIPPR")
        self.Tension.addItem("Parametric")
        self.Tension.addItem("Parachor")
        self.Tension.addItem("Miller")
        self.Tension.addItem("Hakim")
        self.Tension.addItem("Hydrocarbon")
        layout.addWidget(self.Tension, 4, 1)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Liquid Thermal Conductivity:")), 5, 0)
        self.ThCondL = QtGui.QComboBox()
        self.ThCondL.addItem("DIPPR")
        self.ThCondL.addItem("Pachaiyappan")
        layout.addWidget(self.ThCondL, 5, 1)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Gas Thermal Conductivity:")), 6, 0)
        self.ThCondG = QtGui.QComboBox()
        self.ThCondG.addItem("DIPPR")
        self.ThCondG.addItem("Misic & Thodos")
        layout.addWidget(self.ThCondG, 6, 1)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Vapor Pressure:")), 7, 0)
        self.Pv = QtGui.QComboBox()
        self.Pv.addItem("DIPPR")
        self.Pv.addItem("Antoine")
        self.Pv.addItem("Lee-Kesler")
        self.Pv.addItem("Maxwell-Bonnel")
        self.Pv.addItem("Wagner")
        layout.addWidget(self.Pv, 7, 1)
        label_7 = QtGui.QLabel()
        label_7.setAlignment(QtCore.Qt.AlignCenter)
        label_7.setText(QtGui.QApplication.translate(
            "pychemqt", "High Pressure Corrections"))
        layout.addWidget(label_7, 0, 3)
        self.Corr_RhoL = QtGui.QComboBox()
        self.Corr_RhoL.addItem("Thomson, Brobst & Hankinson")
        self.Corr_RhoL.addItem("API")
        layout.addWidget(self.Corr_RhoL, 1, 3)
        self.Corr_MuL = QtGui.QComboBox()
        self.Corr_MuL.addItem("Graboski & Braun")
        self.Corr_MuL.addItem("Kouzel")
        self.Corr_MuL.addItem("Lucas")
        layout.addWidget(self.Corr_MuL, 2, 3)
        self.Corr_ThCondL = QtGui.QComboBox()
        self.Corr_ThCondL.addItem("Lenoir")
        self.Corr_ThCondL.addItem("Kanitkar & Thodos")
        layout.addWidget(self.Corr_ThCondL, 5, 3)
        layout.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding,
                                         QtGui.QSizePolicy.Expanding), 8, 0, 1, 10)

        if config and config.has_section("Transport"):
            self.RhoL.setCurrentIndex(config.getint("Transport", "RhoL"))
            self.Corr_RhoL.setCurrentIndex(config.getint("Transport", "Corr_RhoL"))
            self.MuL.setCurrentIndex(config.getint("Transport", "MuL"))
            self.Corr_MuL.setCurrentIndex(config.getint("Transport", "Corr_MuL"))
            self.MuG.setCurrentIndex(config.getint("Transport", "MuG"))
            self.Tension.setCurrentIndex(config.getint("Transport", "Tension"))
            self.ThCondL.setCurrentIndex(config.getint("Transport", "ThCondL"))
            self.Corr_ThCondL.setCurrentIndex(config.getint("Transport", "Corr_ThCondL"))
            self.ThCondG.setCurrentIndex(config.getint("Transport", "ThCondG"))
            self.Pv.setCurrentIndex(config.getint("Transport", "Pv"))

    def value(self, config):
        """Function to wizard result"""
        if not config.has_section("Transport"):
            config.add_section("Transport")
        config.set("Transport", "RhoL", str(self.RhoL.currentIndex()))
        config.set("Transport", "Corr_RhoL", str(self.Corr_RhoL.currentIndex()))
        config.set("Transport", "MuL", str(self.MuL.currentIndex()))
        config.set("Transport", "Corr_MuL", str(self.Corr_MuL.currentIndex()))
        config.set("Transport", "MuG", str(self.MuG.currentIndex()))
        config.set("Transport", "Tension", str(self.Tension.currentIndex()))
        config.set("Transport", "ThCondL", str(self.ThCondL.currentIndex()))
        config.set("Transport", "Corr_ThCondL", str(self.Corr_ThCondL.currentIndex()))
        config.set("Transport", "ThCondG", str(self.ThCondG.currentIndex()))
        config.set("Transport", "Pv", str(self.Pv.currentIndex()))
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
        return config


class Dialog(QtGui.QDialog):
    """Transport properties dialog"""
    def __init__(self, config, parent=None):
        super(Dialog, self).__init__(parent)
        self.setWindowTitle(QtGui.QApplication.translate(
            "pychemqt", "Transport Properties Methods"))
        layout = QtGui.QVBoxLayout(self)
        self.datos = UI_confTransport_widget(config)
        layout.addWidget(self.datos)
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Cancel |
                                                QtGui.QDialogButtonBox.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox)

    def value(self, config):
        """Function to wizard result"""
        config = self.datos.value(config)
        return config

if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    Dialog = Dialog()
    Dialog.show()
    sys.exit(app.exec_())
