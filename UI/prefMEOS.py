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
# Library to configure MEOS tool
#
#   - Isolinea: widget to configure isolines for mEoS
#   - Widget: mEoS parameter configuration dialog
#   - Dialog: Dialog tool for standalone use
###############################################################################


import os

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
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed, QtWidgets.QSizePolicy.Policy.Fixed),
            6, 1, 1, 4)
        if unit.__name__ != "float" and section != "Psychr":
            self.Critica = QtWidgets.QCheckBox(
                QtWidgets.QApplication.translate(
                    "pychemqt", "Include critic point line"))
            layout.addWidget(self.Critica, 7, 1, 1, 4)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed, QtWidgets.QSizePolicy.Policy.Fixed),
            8, 1, 1, 4)

        self.lineconfig = LineConfig(
            ConfSection,
            QtWidgets.QApplication.translate("pychemqt", "Line Style"))
        layout.addWidget(self.lineconfig, 9, 1, 1, 4)

        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed, QtWidgets.QSizePolicy.Policy.Fixed),
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
        self.Posicion = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.Posicion.valueChanged.connect(self.label5.setValue)
        layout.addWidget(self.Posicion, 14, 3, 1, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 15, 4)

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


class Widget(QtWidgets.QDialog):
    """Config mEoS parameter dialog"""
    lineas = [
        ("Isotherm", unidades.Temperature,
         QtWidgets.QApplication.translate("pychemqt", "Isotherm")),
        ("Isobar", unidades.Pressure,
         QtWidgets.QApplication.translate("pychemqt", "Isobar")),
        ("Isoenthalpic", unidades.Enthalpy,
         QtWidgets.QApplication.translate("pychemqt", "Isoenthalpic")),
        ("Isoentropic", unidades.SpecificHeat,
         QtWidgets.QApplication.translate("pychemqt", "Isoentropic")),
        ("Isochor", unidades.SpecificVolume,
         QtWidgets.QApplication.translate("pychemqt", "Isochor")),
        ("Isoquality", float,
         QtWidgets.QApplication.translate("pychemqt", "Isoquality"))]

    def __init__(self, config, parent=None):
        """constructor, config optional parameter to input project config"""
        super(Widget, self).__init__(parent)

        lyt = QtWidgets.QGridLayout(self)
        lyt.setContentsMargins(0, 0, 0, 0)
        scroll = QtWidgets.QScrollArea()
        scroll.setFrameStyle(QtWidgets.QFrame.Shape.NoFrame)
        lyt.addWidget(scroll)

        dlg = QtWidgets.QWidget()
        layout = QtWidgets.QGridLayout(dlg)

        self.coolProp = QtWidgets.QCheckBox(QtWidgets.QApplication.translate(
            "pychemqt", "Use external library coolProp (faster)"))
        self.coolProp.setEnabled(False)
        layout.addWidget(self.coolProp, 3, 1, 1, 2)
        self.refprop = QtWidgets.QCheckBox(QtWidgets.QApplication.translate(
            "pychemqt", "Use external library refprop (fastest)"))
        self.refprop.setEnabled(False)
        layout.addWidget(self.refprop, 4, 1, 1, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 4, 1)
        self.lineconfig = LineConfig(
            "saturation", QtWidgets.QApplication.translate(
                "pychemqt", "Saturation Line Style"))
        layout.addWidget(self.lineconfig, 5, 1, 1, 2)
        group = QtWidgets.QGroupBox(
            QtWidgets.QApplication.translate("pychemqt", "Isolines"))
        layout.addWidget(group, 6, 1, 1, 2)
        layoutgroup = QtWidgets.QGridLayout(group)
        self.comboIsolineas = QtWidgets.QComboBox()
        layoutgroup.addWidget(self.comboIsolineas, 1, 1)
        self.Isolineas = QtWidgets.QStackedWidget()
        self.comboIsolineas.currentIndexChanged.connect(
            self.Isolineas.setCurrentIndex)
        layoutgroup.addWidget(self.Isolineas, 2, 1)
        for nombre, unidad, text in self.lineas:
            self.comboIsolineas.addItem(text)
            self.Isolineas.addWidget(Isolinea(unidad, nombre, config))
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Plot Definition")), 7, 1)
        quality = [QtWidgets.QApplication.translate("pychemqt", "Very Low"),
                   QtWidgets.QApplication.translate("pychemqt", "Low"),
                   QtWidgets.QApplication.translate("pychemqt", "Medium"),
                   QtWidgets.QApplication.translate("pychemqt", "High"),
                   QtWidgets.QApplication.translate("pychemqt", "Ultra High")]
        self.definition = QtWidgets.QComboBox()
        for q in quality:
            self.definition.addItem(q)
        layout.addWidget(self.definition, 7, 2)
        self.grid = QtWidgets.QCheckBox(
            QtWidgets.QApplication.translate("pychemqt", "Draw grid"))
        layout.addWidget(self.grid, 9, 1, 1, 2)

        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 10, 2)

        scroll.setWidget(dlg)

        if os.environ["CoolProp"]:
            self.coolProp.setEnabled(True)
        if os.environ["refprop"]:
            self.refprop.setEnabled(True)

        if config.has_section("MEOS"):
            self.coolProp.setChecked(config.getboolean("MEOS", 'coolprop'))
            self.refprop.setChecked(config.getboolean("MEOS", 'refprop'))
            self.grid.setChecked(config.getboolean("MEOS", 'grid'))
            self.definition.setCurrentIndex(
                config.getint("MEOS", 'definition'))
            self.lineconfig.setConfig(config)

    def value(self, config):
        """Return value for main dialog"""
        if not config.has_section("MEOS"):
            config.add_section("MEOS")

        config.set("MEOS", "coolprop", str(self.coolProp.isChecked()))
        config.set("MEOS", "refprop", str(self.refprop.isChecked()))
        config = self.lineconfig.value(config)
        config.set("MEOS", "grid", str(self.grid.isChecked()))
        config.set("MEOS", "definition", str(self.definition.currentIndex()))

        for indice in range(self.Isolineas.count()):
            config = self.Isolineas.widget(indice).value(config)
        return config


class Dialog(QtWidgets.QDialog):
    """Dialog to config thermal method calculations"""
    def __init__(self, config=None, parent=None):
        super(Dialog, self).__init__(parent)
        self.setWindowTitle(QtWidgets.QApplication.translate(
            "pychemqt", "Multiparameter equation of state configuration"))
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


if __name__ == "__main__":
    import sys
    from configparser import ConfigParser
    app = QtWidgets.QApplication(sys.argv)

    conf_dir = os.path.expanduser('~') + "/.pychemqt/"
    config = ConfigParser()
    config.read(conf_dir+"pychemqtrc")

    Dialog = Dialog(config)
    Dialog.show()
    sys.exit(app.exec())
