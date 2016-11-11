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
# Library to configure MEOS tool
#
#   - ConfmEoS: mEoS parameter configuration dialog
#   - Isolinea: widget to configure isolines for mEoS
#   - ConfLine: Composite widget with line format configuration tools
###############################################################################


from PyQt5 import QtCore, QtWidgets

from lib import unidades
from lib.utilities import representacion
from UI.widgets import Entrada_con_unidades, LineConfig


class Isolinea(QtWidgets.QDialog):
    """Widget for isoline configuration for mEoS plot tools"""
    def __init__(self, unit, ConfSection, config, section="MEOS", parent=None):
        """Constructor
            unit: subclass of unidad to define the isoline type
            ConfSection: title of isoline
            config: config of pychemqt project
        """
        super(Isolinea, self).__init__(parent)
        self.ConfSection = ConfSection
        self.magnitud = unit.__name__
        self.unidad = unit
        self.section = section
        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Start")), 1, 1)
        self.inicio = Entrada_con_unidades(unit)
        layout.addWidget(self.inicio, 1, 2, 1, 3)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Fin")), 2, 1)
        self.fin = Entrada_con_unidades(unit)
        layout.addWidget(self.fin, 2, 2, 1, 3)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Intervalo")), 3, 1)
        if unit.__name__ == "Temperature":
            self.intervalo = Entrada_con_unidades(unidades.DeltaT)
        elif unit.__name__ == "Pressure":
            self.intervalo = Entrada_con_unidades(unidades.DeltaP)
        else:
            self.intervalo = Entrada_con_unidades(unit)
        layout.addWidget(self.intervalo, 3, 2, 1, 3)
        self.Personalizar = QtWidgets.QCheckBox(
            QtWidgets.QApplication.translate("pychemqt", "Customize"))
        layout.addWidget(self.Personalizar, 4, 1, 1, 4)
        self.Lista = QtWidgets.QLineEdit()
        layout.addWidget(self.Lista, 5, 1, 1, 4)
        self.Personalizar.toggled.connect(self.inicio.setDisabled)
        self.Personalizar.toggled.connect(self.fin.setDisabled)
        self.Personalizar.toggled.connect(self.intervalo.setDisabled)
        self.Personalizar.toggled.connect(self.Lista.setEnabled)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            6, 1, 1, 4)
        if unit.__name__ != "float" and section != "Psychr":
            self.Critica = QtWidgets.QCheckBox(
                QtWidgets.QApplication.translate(
                    "pychemqt", "Include critic point line"))
            layout.addWidget(self.Critica, 7, 1, 1, 4)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            8, 1, 1, 4)

        self.lineconfig = LineConfig(
            ConfSection,
            QtWidgets.QApplication.translate("pychemqt", "Line Style"))
        layout.addWidget(self.lineconfig, 9, 1, 1, 4)

        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            10, 1)
        self.label = QtWidgets.QCheckBox(
            QtWidgets.QApplication.translate("pychemqt", "Label"))
        layout.addWidget(self.label, 11, 1)
        self.variable = QtWidgets.QCheckBox(
            QtWidgets.QApplication.translate("pychemqt", "Variable in Label"))
        layout.addWidget(self.variable, 12, 1, 1, 4)
        self.unit = QtWidgets.QCheckBox(
            QtWidgets.QApplication.translate("pychemqt", "Units in Label"))
        layout.addWidget(self.unit, 13, 1, 1, 4)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Position")), 14, 1)
        self.label5 = Entrada_con_unidades(int, value=0, width=25, frame=False,
                                           readOnly=True)
        self.label5.setFixedWidth(30)
        layout.addWidget(self.label5, 14, 2)
        self.Posicion = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.Posicion.valueChanged.connect(self.label5.setValue)
        layout.addWidget(self.Posicion, 14, 3, 1, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 15, 4)

        if config.has_section(section):
            self.inicio.setValue(config.getfloat(section, ConfSection+'Start'))
            self.fin.setValue(config.getfloat(section, ConfSection+'End'))
            self.intervalo.setValue(
                config.getfloat(section, ConfSection+'Step'))
            self.Personalizar.setChecked(
                config.getboolean(section, ConfSection+'Custom'))
            if config.get(section, ConfSection+'List') != "":
                T = []
                for i in config.get(section, ConfSection+'List').split(','):
                    if unit.__name__ == "float":
                        T.append(representacion(float(i)))
                    else:
                        T.append(representacion(unit(float(i)).config()))
                self.Lista.setText(",".join(T))
            if unit.__name__ != "float" and section != "Psychr":
                self.Critica.setChecked(
                    config.getboolean(section, ConfSection+'Critic'))
            self.inicio.setDisabled(self.Personalizar.isChecked())
            self.fin.setDisabled(self.Personalizar.isChecked())
            self.intervalo.setDisabled(self.Personalizar.isChecked())
            self.Lista.setEnabled(self.Personalizar.isChecked())
            self.label.setChecked(
                config.getboolean(section, ConfSection+'Label'))
            self.variable.setChecked(
                config.getboolean(section, ConfSection+'Variable'))
            self.unit.setChecked(
                config.getboolean(section, ConfSection+'Units'))
            self.Posicion.setValue(
                config.getint(section, ConfSection+'Position'))
            self.lineconfig.setConfig(config, section)

    def value(self, config):
        config.set(self.section, self.ConfSection+"Start",
                   str(self.inicio.value))
        config.set(self.section, self.ConfSection+"End", str(self.fin.value))
        config.set(self.section, self.ConfSection+"Step",
                   str(self.intervalo.value))
        config.set(self.section, self.ConfSection+"Custom",
                   str(self.Personalizar.isChecked()))
        T = []
        if self.Lista.text():
            T1 = self.Lista.text().split(',')
            for i in T1:
                if self.unidad.__name__ == "float":
                    T.append(str(float(i)))
                else:
                    T.append(str(self.unidad(float(i), "conf")))
        config.set(self.section, self.ConfSection+"List", ", ".join(T))
        if self.unidad.__name__ != "float" and self.section != "Psychr":
            config.set(self.section, self.ConfSection+"Critic",
                       str(self.Critica.isChecked()))
        config = self.lineconfig.value(config, self.section)

        config.set(self.section, self.ConfSection+"Label",
                   str(self.label.isChecked()))
        config.set(self.section, self.ConfSection+"Variable",
                   str(self.variable.isChecked()))
        config.set(self.section, self.ConfSection+"Units",
                   str(self.unit.isChecked()))
        config.set(self.section, self.ConfSection+"Position",
                   str(self.Posicion.value()))
        return config

