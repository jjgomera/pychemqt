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
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


###############################################################################
# Units config tools
###############################################################################


from ast import literal_eval
from configparser import ConfigParser
from functools import partial
import os

from tools.qt import QtGui, QtWidgets

from lib import unidades
from lib.config import conf_dir, IMAGE_PATH


class UI_confUnits_widget(QtWidgets.QWidget):
    """Units widget, to use in whatever need, dialog, wizard..."""
    def __init__(self, config=None, parent=None):
        """Constructor, opcional config paramater with project config"""
        super().__init__(parent)

        lyt = QtWidgets.QGridLayout(self)
        systems = QtWidgets.QGroupBox(self.tr("Systems of measurement"))
        systems.setSizePolicy(QtWidgets.QSizePolicy.Policy.Preferred,
                              QtWidgets.QSizePolicy.Policy.Fixed)
        lyt.addWidget(systems, 1, 1, 1, 2)
        lytSystems = QtWidgets.QHBoxLayout(systems)
        self.SI = QtWidgets.QRadioButton(self.tr("SI"))
        self.SI.toggled.connect(partial(self.load, "si"))
        lytSystems.addWidget(self.SI)
        self.AltSI = QtWidgets.QRadioButton(self.tr("Alt SI"))
        self.AltSI.toggled.connect(partial(self.load, "altsi"))
        lytSystems.addWidget(self.AltSI)
        self.English = QtWidgets.QRadioButton(self.tr("English"))
        self.English.toggled.connect(partial(self.load, "english"))
        lytSystems.addWidget(self.English)
        self.Metric = QtWidgets.QRadioButton(self.tr("Metric"))
        self.Metric.toggled.connect(partial(self.load, "metric"))
        lytSystems.addWidget(self.Metric)
        self.CGS = QtWidgets.QRadioButton(self.tr("CGS"))
        self.CGS.toggled.connect(partial(self.load, "cgs"))
        lytSystems.addWidget(self.CGS)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 2, 1, 1, 2)

        table = QtWidgets.QTableWidget()
        lyt.addWidget(table, 3, 1, 5, 1)
        table.setRowCount(len(unidades.MAGNITUDES)-1)
        table.setColumnCount(1)
        table.horizontalHeader().setVisible(False)

        self.combos = []
        for i, magnitud in enumerate(unidades.MAGNITUDES[:-1]):
            table.setVerticalHeaderItem(
                i, QtWidgets.QTableWidgetItem(magnitud[1]))
            table.setRowHeight(i, 24)
            combo = QtWidgets.QComboBox()
            if magnitud[0] == "Currency":
                for texto, unidad in zip(magnitud[2].__text__,
                                         magnitud[2].__units__):
                    combo.addItem(QtGui.QIcon(QtGui.QPixmap(os.path.join(
                        IMAGE_PATH, "flag", f"{unidad}.gif"))),
                        texto + " - " + unidad)
            else:
                for unidad in magnitud[2].__text__:
                    combo.addItem(unidad)
            self.combos.append(combo)
            table.setCellWidget(i, 0, combo)
        table.resizeColumnToContents(0)

        table.setFixedWidth(table.verticalHeader().sizeHint().width()
                            + table.columnWidth(0) + 20)

        self.name = QtWidgets.QLineEdit()
        self.name.textChanged.connect(self.nameChanged)
        lyt.addWidget(self.name, 3, 2)
        self.Guardar = QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(
            os.path.join(IMAGE_PATH, "button", "fileSave.png"))),
            self.tr("Save profile"))
        self.Guardar.setEnabled(False)
        self.Guardar.clicked.connect(self.save)
        lyt.addWidget(self.Guardar, 4, 2)
        self.profiles = QtWidgets.QComboBox()
        lyt.addWidget(self.profiles, 6, 2)
        self.loadButton = QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(
            os.path.join(IMAGE_PATH, "button", "fileOpen.png"))),
            self.tr("Load profile"))
        self.loadButton.clicked.connect(self.load)
        lyt.addWidget(self.loadButton, 7, 2)

        self.updateProfileList()
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

    def updateProfileList(self):
        """Update custom units profile list"""
        self.profiles.clear()
        if os.path.isfile(conf_dir+"unitrc"):
            Config = ConfigParser()
            Config.read(conf_dir+"unitrc")
            for i in Config.options("units"):
                self.profiles.addItem(self.tr(i))
            self.profiles.setEnabled(True)
            self.loadButton.setEnabled(True)
        else:
            self.profiles.setEnabled(False)
            self.loadButton.setEnabled(False)

    def save(self):
        """Save units profile to a file"""
        lista = []
        for combo in self.combos:
            lista.append(combo.currentIndex())
        Config = ConfigParser()
        if os.path.isfile(conf_dir+"unitrc"):
            Config.read(conf_dir+"unitrc")
        else:
            Config.add_section("units")
        Config.set("units", str(self.name.text()), str(lista))
        with open(conf_dir+"unitrc", "w") as unit_file:
            Config.write(unit_file)
        self.updateProfileList()

    def nameChanged(self, name):
        """Update state of boton Guardar"""
        self.Guardar.setEnabled(bool(name))

    def load(self, unitset=None):
        """Load unit set, opcional set parameter can be any of predefined:
            si, altsi, metric, cgs, english"""
        if unitset:
            lista = unidades.units_set[unitset]
        else:
            Config = ConfigParser()
            Config.read(conf_dir+"unitrc")
            lista = literal_eval(
                Config.get("units", str(self.profiles.currentText())))
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
        """Default configuration"""
        config.add_section("Units")
        config.set("Units", "System", "0")
        for magnitud in unidades.MAGNITUDES[:-1]:
            config.set("Units", magnitud[0], "0")
        return config


class Dialog(QtWidgets.QDialog):
    """Units configuration dialog"""
    def __init__(self, config=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle(self.tr("Units"))
        lyt = QtWidgets.QVBoxLayout(self)
        self.datos = UI_confUnits_widget(config)
        lyt.addWidget(self.datos)
        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Cancel
            | QtWidgets.QDialogButtonBox.StandardButton.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        lyt.addWidget(self.buttonBox)

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
