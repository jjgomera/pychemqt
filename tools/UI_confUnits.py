#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2023, Juan José Gómez Romera <jjgomera@gmail.com>

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
# Units config tools
###############################################################################

import os
from functools import partial
from configparser import ConfigParser

from tools.qt import QtGui, QtWidgets


from lib import unidades
from lib.config import conf_dir


class UI_confUnits_widget(QtWidgets.QWidget):
    """Units widget, to use in whatever need, dialog, wizard..."""
    def __init__(self, config=None, parent=None):
        """Constructor, opcional config paramater with project config"""
        super(UI_confUnits_widget, self).__init__(parent)

        layout = QtWidgets.QGridLayout(self)
        systems = QtWidgets.QGroupBox(
            QtWidgets.QApplication.translate("pychemqt", "Systems of measurement"))
        systems.setSizePolicy(QtWidgets.QSizePolicy.Policy.Preferred,
                              QtWidgets.QSizePolicy.Policy.Fixed)
        layout.addWidget(systems, 1, 1, 1, 2)
        lytSystems = QtWidgets.QHBoxLayout(systems)
        self.SI = QtWidgets.QRadioButton(
            QtWidgets.QApplication.translate("pychemqt", "SI"))
        self.SI.toggled.connect(partial(self.load, "si"))
        lytSystems.addWidget(self.SI)
        self.AltSI = QtWidgets.QRadioButton(
            QtWidgets.QApplication.translate("pychemqt", "Alt SI"))
        self.AltSI.toggled.connect(partial(self.load, "altsi"))
        lytSystems.addWidget(self.AltSI)
        self.English = QtWidgets.QRadioButton(
            QtWidgets.QApplication.translate("pychemqt", "English"))
        self.English.toggled.connect(partial(self.load, "english"))
        lytSystems.addWidget(self.English)
        self.Metric = QtWidgets.QRadioButton(
            QtWidgets.QApplication.translate("pychemqt", "Metric"))
        self.Metric.toggled.connect(partial(self.load, "metric"))
        lytSystems.addWidget(self.Metric)
        self.CGS = QtWidgets.QRadioButton(
            QtWidgets.QApplication.translate("pychemqt", "CGS"))
        self.CGS.toggled.connect(partial(self.load, "cgs"))
        lytSystems.addWidget(self.CGS)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed, QtWidgets.QSizePolicy.Policy.Fixed), 2, 1, 1, 2)

        self.tabla = QtWidgets.QTableWidget()
        layout.addWidget(self.tabla, 3, 1, 5, 1)
        self.tabla.setRowCount(len(unidades.MAGNITUDES)-1)
        self.tabla.setColumnCount(1)
        self.tabla.horizontalHeader().setVisible(False)

        self.combos = []
        for i, magnitud in enumerate(unidades.MAGNITUDES[:-1]):
            self.tabla.setVerticalHeaderItem(
                i, QtWidgets.QTableWidgetItem(magnitud[1]))
            self.tabla.setRowHeight(i, 24)
            combo = QtWidgets.QComboBox()
            if magnitud[0] == "Currency":
                for texto, unidad in zip(magnitud[2].__text__,
                                         magnitud[2].__units__):
                    combo.addItem(QtGui.QIcon(QtGui.QPixmap(
                        os.environ["pychemqt"]+"/images/flag/%s.gif" % unidad)),
                        texto + " - " + unidad)
            else:
                for unidad in magnitud[2].__text__:
                    combo.addItem(unidad)
            self.combos.append(combo)
            self.tabla.setCellWidget(i, 0, combo)
        self.tabla.resizeColumnToContents(0)

        self.tabla.setFixedWidth(self.tabla.verticalHeader().sizeHint().width()
                                 + self.tabla.columnWidth(0) + 20)

        self.nombre = QtWidgets.QLineEdit()
        self.nombre.textChanged.connect(self.nameChanged)
        layout.addWidget(self.nombre, 3, 2)
        self.Guardar = QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"]+"/images/button/fileSave.png")),
            QtWidgets.QApplication.translate("pychemqt", "Save profile"))
        self.Guardar.setEnabled(False)
        self.Guardar.clicked.connect(self.guardar)
        layout.addWidget(self.Guardar, 4, 2)
        self.perfiles = QtWidgets.QComboBox()
        layout.addWidget(self.perfiles, 6, 2)
        self.Cargar = QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"]+"/images/button/fileOpen.png")),
            QtWidgets.QApplication.translate("pychemqt", "Load profile"))
        self.Cargar.clicked.connect(self.load)
        layout.addWidget(self.Cargar, 7, 2)

        self.actualizar_lista_perfiles()
        if config and config.has_section("Units"):
            if config.getint("Units", "System") == 0:
                self.SI.setChecked(True)
                self.load("si")
            elif config.getint("Units", "System") == 1:
                self.Metric.setChecked(True)
                self.load("metric")
            elif config.getint("Units", "System") == 2:
                self.CGS.setChecked(True)
                self.load("cgs")
            elif config.getint("Units", "System") == 3:
                self.English.setChecked(True)
                self.load("english")
            elif config.getint("Units", "System") == 4:
                self.AltSI.setChecked(True)
                self.load("altsi")

            for combo, magnitud in zip(self.combos, unidades.MAGNITUDES[:-1]):
                combo.setCurrentIndex(config.getint("Units", magnitud[0]))

    def actualizar_lista_perfiles(self):
        """Update custom units profile list"""
        self.perfiles.clear()
        if os.path.isfile(conf_dir+"unitrc"):
            Config = ConfigParser()
            Config.read(conf_dir+"unitrc")
            for i in Config.options("units"):
                self.perfiles.addItem(QtWidgets.QApplication.translate("pychemqt", i))
            self.perfiles.setEnabled(True)
            self.Cargar.setEnabled(True)
        else:
            self.perfiles.setEnabled(False)
            self.Cargar.setEnabled(False)

    def guardar(self):
        """Save units profile to a file"""
        lista = []
        for combo in self.combos:
            lista.append(combo.currentIndex())
        Config = ConfigParser()
        if os.path.isfile(conf_dir+"unitrc"):
            Config.read(conf_dir+"unitrc")
        else:
            Config.add_section("units")
        Config.set("units", str(self.nombre.text()), lista)
        Config.write(open(conf_dir+"unitrc", "w"))
        self.actualizar_lista_perfiles()

    def nameChanged(self, nombre):
        if nombre:
            self.Guardar.setEnabled(True)
        else:
            self.Guardar.setEnabled(False)

    def load(self, set=None):
        """Load unit set, opcional set parameter can be any of predefined:
            si, altsi, metric, cgs, english"""
        if set:
            lista = unidades.units_set[set]
        else:
            Config = ConfigParser()
            Config.read(conf_dir+"unitrc")
            lista = eval(Config.get("units", str(self.perfiles.currentText())))
        for combo, indice in zip(self.combos, lista):
            combo.setCurrentIndex(indice)

    def value(self, config):
        """Function to wizard result"""
        if not config.has_section("Units"):
            config.add_section("Units")

        if self.SI.isChecked():
            config.set("Units", "System", "0")
        elif self.Metric.isChecked():
            config.set("Units", "System", "1")
        elif self.CGS.isChecked():
            config.set("Units", "System", "2")
        elif self.English.isChecked():
            config.set("Units", "System", "3")
        elif self.AltSI.isChecked():
            config.set("Units", "System", "4")
        else:
            config.set("Units", "System", "5")

        for combo, magnitud in zip(self.combos, unidades.MAGNITUDES[:-1]):
            config.set("Units", magnitud[0], str(combo.currentIndex()))
        return config

    @classmethod
    def default(cls, config):
        config.add_section("Units")
        config.set("Units", "System", "0")
        for magnitud in unidades.MAGNITUDES[:-1]:
            config.set("Units", magnitud[0], "0")
        return config


class Dialog(QtWidgets.QDialog):
    """Units configuration dialog"""
    def __init__(self, config=None, parent=None):
        super(Dialog, self).__init__(parent)
        self.setWindowTitle(QtWidgets.QApplication.translate("pychemqt", "Units"))
        layout = QtWidgets.QVBoxLayout(self)
        self.datos = UI_confUnits_widget(config)
        layout.addWidget(self.datos)
        self.buttonBox = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Cancel |
                                                QtWidgets.QDialogButtonBox.StandardButton.Ok)
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
    sys.exit(app.exec())
