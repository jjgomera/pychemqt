#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Transport properties config section
###############################################################################

from PyQt5 import QtCore, QtWidgets



class UI_confTransport_widget(QtWidgets.QWidget):
    """Transport properties widget, tu use in dialog, wizard..."""
    def __init__(self, config=None, parent=None):
        """Constructor, opcional config parameter with project config"""
        super(UI_confTransport_widget, self).__init__(parent)
        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Liquid Density:")), 1, 0)
        self.RhoL = QtWidgets.QComboBox()
        self.RhoL.addItem("DIPPR")
        self.RhoL.addItem("Rackett")
        self.RhoL.addItem("Cavett")
        self.RhoL.addItem("Costald")
        layout.addWidget(self.RhoL, 1, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Liquid Viscosity:")), 2, 0)
        self.MuL = QtWidgets.QComboBox()
        self.MuL.addItem("DIPPR")
        self.MuL.addItem("Parametric")
        self.MuL.addItem("Letsou & Steil")
        layout.addWidget(self.MuL, 2, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Gas Viscosity:")), 3, 0)
        self.MuG = QtWidgets.QComboBox()
        self.MuG.addItem("DIPPR")
        self.MuG.addItem("Chapman & Enskog")
        self.MuG.addItem("Thodos")
        layout.addWidget(self.MuG, 3, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Surface Tension:")), 4, 0)
        self.Tension = QtWidgets.QComboBox()
        self.Tension.addItem("DIPPR")
        self.Tension.addItem("Parametric")
        self.Tension.addItem("Parachor")
        self.Tension.addItem("Miller")
        self.Tension.addItem("Hakim")
        self.Tension.addItem("Hydrocarbon")
        layout.addWidget(self.Tension, 4, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Liquid Thermal Conductivity:")), 5, 0)
        self.ThCondL = QtWidgets.QComboBox()
        self.ThCondL.addItem("DIPPR")
        self.ThCondL.addItem("Pachaiyappan")
        layout.addWidget(self.ThCondL, 5, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Gas Thermal Conductivity:")), 6, 0)
        self.ThCondG = QtWidgets.QComboBox()
        self.ThCondG.addItem("DIPPR")
        self.ThCondG.addItem("Misic & Thodos")
        layout.addWidget(self.ThCondG, 6, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Vapor Pressure:")), 7, 0)
        self.Pv = QtWidgets.QComboBox()
        self.Pv.addItem("DIPPR")
        self.Pv.addItem("Antoine")
        self.Pv.addItem("Lee-Kesler")
        self.Pv.addItem("Maxwell-Bonnel")
        self.Pv.addItem("Wagner")
        layout.addWidget(self.Pv, 7, 1)
        label_7 = QtWidgets.QLabel()
        label_7.setAlignment(QtCore.Qt.AlignCenter)
        label_7.setText(QtWidgets.QApplication.translate(
            "pychemqt", "High Pressure Corrections"))
        layout.addWidget(label_7, 0, 3)
        self.Corr_RhoL = QtWidgets.QComboBox()
        self.Corr_RhoL.addItem("Thomson, Brobst & Hankinson")
        self.Corr_RhoL.addItem("API")
        layout.addWidget(self.Corr_RhoL, 1, 3)
        self.Corr_MuL = QtWidgets.QComboBox()
        self.Corr_MuL.addItem("Graboski & Braun")
        self.Corr_MuL.addItem("Kouzel")
        self.Corr_MuL.addItem("Lucas")
        layout.addWidget(self.Corr_MuL, 2, 3)
        self.Corr_ThCondL = QtWidgets.QComboBox()
        self.Corr_ThCondL.addItem("Lenoir")
        self.Corr_ThCondL.addItem("Kanitkar & Thodos")
        layout.addWidget(self.Corr_ThCondL, 5, 3)
        layout.addItem(QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding,
                                         QtWidgets.QSizePolicy.Expanding), 8, 0, 1, 10)

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


class Dialog(QtWidgets.QDialog):
    """Transport properties dialog"""
    def __init__(self, config, parent=None):
        super(Dialog, self).__init__(parent)
        self.setWindowTitle(QtWidgets.QApplication.translate(
            "pychemqt", "Transport Properties Methods"))
        layout = QtWidgets.QVBoxLayout(self)
        self.datos = UI_confTransport_widget(config)
        layout.addWidget(self.datos)
        self.buttonBox = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Cancel |
                                                QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox)

    def value(self, config):
        """Function to wizard result"""
        config = self.datos.value(config)
        return config

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = Dialog()
    Dialog.show()
    sys.exit(app.exec_())
