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
# Library to show/configure pychemqt general preferences
#
#   Preferences: Preferences main dialog
#   - ConfGeneral: General configuration options
#   - ConfTooltipUnit: Tooltip with unit alternate value configuration
#   - ConfFormat: Numeric format configuration
#   - ConfApplications: External applications configuration
#   - ConfTooltipEntity: Entity properties in popup window configuration
###############################################################################


from configparser import ConfigParser
from functools import partial
from math import pi
import os
import sys

from PyQt5 import QtCore, QtGui, QtWidgets

from lib import unidades, corriente
from lib.utilities import representacion
from equipment import equipments
from tools.firstrun import which
from UI import prefElemental, prefMEOS, prefPFD, prefPsychrometric, prefPetro
from plots import drag, moody, standing
from UI.delegate import CheckEditor
from UI.widgets import ColorSelector, NumericFactor, PathConfig


class ConfGeneral(QtWidgets.QDialog):
    """General configuration options"""
    def __init__(self, config=None, parent=None):
        super(ConfGeneral, self).__init__(parent)

        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Highlight color:")), 1, 1)
        self.ColorButtonResaltado = ColorSelector()
        layout.addWidget(self.ColorButtonResaltado, 1, 2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Readonly color:")), 2, 1)
        self.ColorButtonReadOnly = ColorSelector()
        layout.addWidget(self.ColorButtonReadOnly, 2, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 0, QtWidgets.QSizePolicy.Fixed,
            QtWidgets.QSizePolicy.Fixed), 3, 1)
        group = QtWidgets.QGroupBox(
            QtWidgets.QApplication.translate("pychemqt", "Recent Files"))
        layout.addWidget(group, 4, 1, 1, 4)
        lyt = QtWidgets.QHBoxLayout(group)
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Number of recent files:")))
        self.recentFiles = QtWidgets.QSpinBox()
        self.recentFiles.setRange(1, 20)
        lyt.addWidget(self.recentFiles)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 0, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Fixed))
        self.loadLastProject = QtWidgets.QCheckBox(
            QtWidgets.QApplication.translate(
                "pychemqt", "Load last session project"))
        layout.addWidget(self.loadLastProject, 5, 1)
        self.showTrayIcon = QtWidgets.QCheckBox(
            QtWidgets.QApplication.translate("pychemqt", "Show tray icon"))
        layout.addWidget(self.showTrayIcon, 6, 1)

        layout.addItem(QtWidgets.QSpacerItem(
            10, 0, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 14, 1, 1, 4)

        if config and config.has_section("General"):
            self.ColorButtonResaltado.setColor(
                config.get("General", 'Color_Resaltado'))
            self.ColorButtonReadOnly.setColor(
                config.get("General", 'Color_ReadOnly'))
            self.recentFiles.setValue(config.getint("General", 'Recent_Files'))
            self.loadLastProject.setChecked(
                config.getboolean("General", 'Load_Last_Project'))
            self.showTrayIcon.setChecked(config.getboolean("General", 'Tray'))

    def value(self, config):
        if not config.has_section("General"):
            config.add_section("General")
        config.set("General", "Color_Resaltado",
                   self.ColorButtonResaltado.color.name())
        config.set("General", "Color_ReadOnly",
                   self.ColorButtonReadOnly.color.name())
        config.set("General", "Recent_Files", str(self.recentFiles.value()))
        config.set("General", "Load_Last_Project",
                   str(self.loadLastProject.isChecked()))
        config.set("General", "Tray", str(self.showTrayIcon.isChecked()))
        return config


class ConfTooltipUnit(QtWidgets.QDialog):
    """Tooltip with unit alternate value configuration"""
    def __init__(self, config, parent=None):
        super(ConfTooltipUnit, self).__init__(parent)

        layout = QtWidgets.QVBoxLayout(self)
        self.checkShow = QtWidgets.QCheckBox(
            QtWidgets.QApplication.translate("pychemqt", "Show Tool Tips"))
        self.checkShow.toggled.connect(self.checkShow_Toggled)
        layout.addWidget(self.checkShow)

        self.groupsystems = QtWidgets.QGroupBox(
            QtWidgets.QApplication.translate(
                "pychemqt", "Systems of measurement"))
        layout.addWidget(self.groupsystems)
        lytSystems = QtWidgets.QHBoxLayout(self.groupsystems)
        self.SI = QtWidgets.QCheckBox(
            QtWidgets.QApplication.translate("pychemqt", "SI"))
        self.SI.toggled.connect(partial(self.systems, "si"))
        lytSystems.addWidget(self.SI)
        self.AltSI = QtWidgets.QCheckBox(
            QtWidgets.QApplication.translate("pychemqt", "Alt SI"))
        self.AltSI.toggled.connect(partial(self.systems, "altsi"))
        lytSystems.addWidget(self.AltSI)
        self.English = QtWidgets.QCheckBox(
            QtWidgets.QApplication.translate("pychemqt", "English"))
        self.English.toggled.connect(partial(self.systems, "english"))
        lytSystems.addWidget(self.English)
        self.Metric = QtWidgets.QCheckBox(
            QtWidgets.QApplication.translate("pychemqt", "Metric"))
        self.Metric.toggled.connect(partial(self.systems, "metric"))
        lytSystems.addWidget(self.Metric)
        self.CGS = QtWidgets.QCheckBox(
            QtWidgets.QApplication.translate("pychemqt", "CGS"))
        self.CGS.toggled.connect(partial(self.systems, "cgs"))
        lytSystems.addWidget(self.CGS)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed))

        self.eleccion = QtWidgets.QComboBox()
        layout.addWidget(self.eleccion)
        self.stacked = QtWidgets.QStackedWidget()
        self.eleccion.currentIndexChanged.connect(self.stacked.setCurrentIndex)
        layout.addWidget(self.stacked)

        self.tabla = []
        for i, magnitud in enumerate(unidades._magnitudes[:-1]):
            textos = magnitud[2].__text__
            self.tabla.append(QtWidgets.QTableWidget())
            self.stacked.addWidget(self.tabla[i])

            self.tabla[i].setRowCount(len(textos))
            self.tabla[i].setColumnCount(1)
            self.tabla[i].setColumnWidth(0, 16)
            self.tabla[i].setItemDelegateForColumn(0, CheckEditor(self))
            self.tabla[i].horizontalHeader().setVisible(False)
            for j in range(len(textos)):
                item = QtWidgets.QTableWidgetItem(textos[j])
                self.tabla[i].setVerticalHeaderItem(j, item)
                self.tabla[i].setRowHeight(j, 24)
                self.tabla[i].setItem(j, 0, QtWidgets.QTableWidgetItem(""))
                self.tabla[i].item(j, 0).setTextAlignment(
                    QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
                self.tabla[i].openPersistentEditor(self.tabla[i].item(j, 0))
            self.rellenar(magnitud[0], i, config)
            self.eleccion.addItem(magnitud[1])

        if config.has_section("Tooltip"):
            self.checkShow.setChecked(config.getboolean("Tooltip", "Show"))
            self.SI.setChecked(config.getboolean("Tooltip", "SI"))
            self.AltSI.setChecked(config.getboolean("Tooltip", "AltSI"))
            self.English.setChecked(config.getboolean("Tooltip", "English"))
            self.Metric.setChecked(config.getboolean("Tooltip", "Metric"))
            self.CGS.setChecked(config.getboolean("Tooltip", "CGS"))

    def rellenar(self, magnitud, tabla, config):
        if config.has_section("Tooltip"):
            lista = eval(config.get("Tooltip", magnitud))
            for i in lista:
                self.tabla[tabla].item(i, 0).setText("true")

    def checkShow_Toggled(self, bool):
        self.eleccion.setEnabled(bool)
        self.groupsystems.setEnabled(bool)
        for tabla in self.tabla:
            tabla.setEnabled(bool)

    def systems(self, set, bool):
        if bool:
            txt = "true"
        else:
            txt = ""
        for tabla, value in enumerate(unidades.units_set[set][:-1]):
            self.tabla[tabla].item(value, 0).setText(txt)

    def value(self, config):
        if not config.has_section("Tooltip"):
            config.add_section("Tooltip")
        config.set("Tooltip", "Show", str(self.checkShow.isChecked()))
        config.set("Tooltip", "SI", str(self.SI.isChecked()))
        config.set("Tooltip", "CGS", str(self.CGS.isChecked()))
        config.set("Tooltip", "AltSI", str(self.AltSI.isChecked()))
        config.set("Tooltip", "English", str(self.English.isChecked()))
        config.set("Tooltip", "Metric", str(self.Metric.isChecked()))

        for i, tabla in enumerate(self.tabla):
            lista = []
            for j in range(tabla.rowCount()):
                if tabla.item(j, 0).text() == "true":
                    lista.append(j)
            config.set("Tooltip", unidades._magnitudes[i][0], str(lista))
        return config


class ConfTooltipEntity(QtWidgets.QDialog):
    """Entity properties in popup window configuration"""
    def __init__(self, config, parent=None):
        super(ConfTooltipEntity, self).__init__(parent)

        layout = QtWidgets.QVBoxLayout(self)
        self.eleccion = QtWidgets.QComboBox()
        layout.addWidget(self.eleccion)
        self.stacked = QtWidgets.QStackedWidget()
        self.eleccion.currentIndexChanged.connect(self.stacked.setCurrentIndex)
        layout.addWidget(self.stacked)

        self.tabla = [QtWidgets.QTableWidget()]
        self.tabla[0].setRowCount(len(corriente.Corriente.propertiesNames()))
        self.tabla[0].setColumnCount(1)
        self.tabla[0].setColumnWidth(0, 16)
        self.tabla[0].setItemDelegateForColumn(0, CheckEditor(self))
        self.tabla[0].horizontalHeader().setVisible(False)
        self.stacked.addWidget(self.tabla[0])
        self.eleccion.addItem(
            QtWidgets.QApplication.translate("pychemqt", "Stream"))

        for i, propiedad in enumerate(corriente.Corriente.propertiesNames()):
            item = QtWidgets.QTableWidgetItem(propiedad[0])
            self.tabla[0].setVerticalHeaderItem(i, item)
            self.tabla[0].setRowHeight(i, 24)
            self.tabla[0].setItem(i, 0, QtWidgets.QTableWidgetItem(""))
            self.tabla[0].item(i, 0).setTextAlignment(
                QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
            self.tabla[0].openPersistentEditor(self.tabla[0].item(i, 0))

        if config.has_option("TooltipEntity", "Corriente"):
            lista = eval(config.get("TooltipEntity", "Corriente"))
            for i in lista:
                self.tabla[0].item(i, 0).setText("true")

        for i, equipo in enumerate(equipments):
            propiedades = [prop[0] for prop in equipo.propertiesNames()]
            self.tabla.append(QtWidgets.QTableWidget())
            self.stacked.addWidget(self.tabla[-1])
            self.tabla[-1].setRowCount(len(propiedades))
            self.tabla[-1].setColumnCount(1)
            self.tabla[-1].setColumnWidth(0, 16)
            self.tabla[-1].setItemDelegateForColumn(0, CheckEditor(self))
            self.tabla[-1].horizontalHeader().setVisible(False)
            for j, propiedad in enumerate(propiedades):
                item = QtWidgets.QTableWidgetItem(propiedad)
                self.tabla[-1].setVerticalHeaderItem(j, item)
                self.tabla[-1].setRowHeight(j, 24)
                self.tabla[-1].setItem(j, 0, QtWidgets.QTableWidgetItem(""))
                self.tabla[-1].item(j, 0).setTextAlignment(
                    QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
                self.tabla[-1].openPersistentEditor(self.tabla[-1].item(j, 0))
            self.rellenar(equipo.__name__, i+1, config)
            self.eleccion.addItem(equipo.title)

    def rellenar(self, equipo, tabla, config):
        if config.has_section("TooltipEntity"):
            lista = eval(config.get("TooltipEntity", equipo))
            for i in lista:
                self.tabla[tabla].item(i, 0).setText("true")

    def value(self, config):
        if not config.has_section("TooltipEntity"):
            config.add_section("TooltipEntity")
        lista = []
        for j in range(self.tabla[0].rowCount()):
            if self.tabla[0].item(j, 0).text() == "true":
                lista.append(j)
        config.set("TooltipEntity", "Corriente", str(lista))

        for i, tabla in enumerate(self.tabla[1:]):
            lista = []
            for j in range(tabla.rowCount()):
                if tabla.item(j, 0).text() == "true":
                    lista.append(j)
            config.set("TooltipEntity", equipments[i].__name__, str(lista))

        return config


class ConfFormat(QtWidgets.QTableWidget):
    """Numeric format configuration"""
    def __init__(self, config=None, parent=None):
        super(ConfFormat, self).__init__(parent)
        self.setColumnCount(2)
        self.setRowCount(len(unidades._magnitudes))
        labels = [QtWidgets.QApplication.translate("pychemqt", "Format"),
                  QtWidgets.QApplication.translate("pychemqt", "Sample")]
        self.setHorizontalHeaderLabels(labels)
        self.config = []
        for i in range(self.rowCount()):
            item = QtWidgets.QTableWidgetItem(unidades._magnitudes[i][1])
            self.setVerticalHeaderItem(i, item)
            self.setRowHeight(i, 22)
            self.setItem(i, 0, QtWidgets.QTableWidgetItem(""))
            self.item(i, 0).setTextAlignment(
                QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
            self.setItem(i, 1, QtWidgets.QTableWidgetItem(""))
            self.item(i, 1).setTextAlignment(
                QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)

        if config.has_section("NumericFormat"):
            for i, magnitud in enumerate(unidades._magnitudes):
                formato = eval(config.get("NumericFormat", magnitud[0]))
                self.config.append(formato)
                self.item(i, 0).setText(self.txt(formato))
                self.item(i, 1).setText(representacion(pi, **formato))

        self.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.cellDoubleClicked.connect(self.showConfDialog)

    def showConfDialog(self, fila, columna):
        dialog = NumericFactor(self.config[fila], parent=self)
        if dialog.exec_():
            config = dialog.args()
            self.config[fila] = config
            self.item(fila, 0).setText(self.txt(config))
            self.item(fila, 1).setText(representacion(pi, **config))

    def txt(self, formato):
        if formato["signo"]:
            txt = "+"
        else:
            txt = ""
        if formato["fmt"] == 0:
            txt += "{total}.{decimales} fixed".format(**formato)
        elif formato["fmt"] == 1:
            txt += "{decimales} sign".format(**formato)
        elif formato["fmt"] == 2:
            txt += "{decimales} exp".format(**formato)
        if formato.get("exp", False):
            txt += " ({tol} exp)".format(**formato)
        return txt

    def value(self, config):
        if not config.has_section("NumericFormat"):
            config.add_section("NumericFormat")
        for i, magnitud in enumerate(unidades._magnitudes):
            config.set("NumericFormat", magnitud[0], str(self.config[i]))
        return config


class ConfApplications(QtWidgets.QDialog):
    """External applications configuration"""
    def __init__(self, config=None, parent=None):
        super(ConfApplications, self).__init__(parent)
        layout = QtWidgets.QGridLayout(self)
        l = QtWidgets.QApplication.translate("pychemqt", "External Calculator")
        msg = QtWidgets.QApplication.translate(
            "pychemqt", "Select External Calculator Application")
        self.calculadora = PathConfig(l + ":", msg=msg, patron="exe")
        layout.addWidget(self.calculadora, 1, 1)
        l = QtWidgets.QApplication.translate("pychemqt", "Text viewer")
        msg = QtWidgets.QApplication.translate(
            "pychemqt", "Select External Text Viewer Application")
        self.textViewer = PathConfig(l + ":", msg=msg, patron="exe")
        layout.addWidget(self.textViewer, 2, 1)

        terminal = QtWidgets.QGroupBox()
        layout.addWidget(terminal, 3, 1)
        layoutTerminal = QtWidgets.QGridLayout(terminal)
        msg = QtWidgets.QApplication.translate(
            "pychemqt", "Select External shell")
        self.terminal = PathConfig("", msg=msg, patron="exe")
        layoutTerminal.addWidget(self.terminal, 1, 1, 1, 3)
        layoutTerminal.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Foreground color:")),
            2, 1)
        self.ForegroundColor = ColorSelector()
        layoutTerminal.addWidget(self.ForegroundColor, 2, 2, 1, 2)
        layoutTerminal.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Background color:")),
            3, 1)
        self.BackgroundColor = ColorSelector()
        layoutTerminal.addWidget(self.BackgroundColor, 3, 2, 1, 2)
        self.ipython = QtWidgets.QCheckBox(QtWidgets.QApplication.translate(
            "pychemqt", "Use ipython if its available"))
        layoutTerminal.addWidget(self.ipython, 4, 1, 1, 3)
        self.maximized = QtWidgets.QCheckBox(
            QtWidgets.QApplication.translate("pychemqt", "Show maximized"))
        layoutTerminal.addWidget(self.maximized, 5, 1, 1, 3)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 0, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 10, 1)

        terminalTitle = QtWidgets.QApplication.translate("pychemqt", "Shell")
        if sys.platform == "win32":
            terminal.setEnabled(False)
            terminalTitle += " (" + QtWidgets.QApplication.translate(
                "pychemqt", "Only Available on linux") + ")"
        terminal.setTitle(terminalTitle)

        if config.has_section("Applications"):
            self.calculadora.setText(config.get("Applications", 'Calculator'))
            self.textViewer.setText(config.get("Applications", 'TextViewer'))
            self.terminal.setText(config.get("Applications", 'Shell'))
            self.ipython.setChecked(
                config.getboolean("Applications", 'ipython'))
            self.maximized.setChecked(
                config.getboolean("Applications", 'maximized'))
            self.ForegroundColor.setColor(
                config.get("Applications", 'foregroundColor'))
            self.BackgroundColor.setColor(
                config.get("Applications", 'backgroundColor'))

        self.ipython.setEnabled(bool(which("ipython3")))

        # TODO: Habilitar cuando añada soporte para otras terminales
        self.terminal.setEnabled(False)

    def value(self, config):
        """Return value for main dialog"""
        if not config.has_section("Applications"):
            config.add_section("Applications")
        config.set("Applications", "Calculator", self.calculadora.text())
        config.set("Applications", "TextViewer", self.textViewer.text())
        config.set("Applications", "Shell",  self.terminal.text())
        config.set("Applications", "ipython", str(self.ipython.isChecked()))
        config.set("Applications", "maximized",
                   str(self.maximized.isChecked()))
        config.set("Applications", "foregroundColor",
                   self.ForegroundColor.color.name())
        config.set("Applications", "backgroundColor",
                   self.BackgroundColor.color.name())
        return config


class ConfBabel(QtWidgets.QDialog):
    """Openbabel image generation configuration options"""
    def __init__(self, config=None, parent=None):
        super(ConfBabel, self).__init__(parent)

        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Bond color:")), 1, 1)
        self.BondColor = ColorSelector()
        layout.addWidget(self.BondColor, 1, 2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Foreground color:")), 2, 1)
        self.BackgroundColor = ColorSelector(isAlpha=True)
        layout.addWidget(self.BackgroundColor, 2, 2)
        self.checkColor = QtWidgets.QCheckBox(QtWidgets.QApplication.translate(
            "pychemqt", "Heteroatom in color"))
        layout.addWidget(self.checkColor, 3, 1, 1, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            4, 1, 1, 2)

        group = QtWidgets.QGroupBox(
            QtWidgets.QApplication.translate("pychemqt", "Atom details"))
        layout.addWidget(group, 5, 1, 1, 2)
        lyt = QtWidgets.QVBoxLayout(group)
        self.radioAll = QtWidgets.QRadioButton(
            QtWidgets.QApplication.translate("pychemqt", "Show all atoms"))
        lyt.addWidget(self.radioAll)
        self.radioEnd = QtWidgets.QRadioButton(
            QtWidgets.QApplication.translate(
                "pychemqt", "Show only terminal atoms"))
        lyt.addWidget(self.radioEnd)
        self.radioNone = QtWidgets.QRadioButton(
            QtWidgets.QApplication.translate("pychemqt", "Do not show atoms"))
        lyt.addWidget(self.radioNone)
        layout.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            6, 1, 1, 2)
        self.checkTighBond = QtWidgets.QCheckBox(
            QtWidgets.QApplication.translate("pychemqt", "Thicker bond lines"))
        layout.addWidget(self.checkTighBond, 7, 1, 1, 2)

        layout.addItem(QtWidgets.QSpacerItem(
            10, 0, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 14, 1, 1, 4)

        if config and config.has_section("Openbabel"):
            self.BondColor.setColor(config.get("Openbabel", 'BondColor'))
            alpha = config.getint("Openbabel", "BackColorAlpha")
            self.BackgroundColor.setColor(
                config.get("Openbabel", 'BackColor'), alpha)
            self.checkColor.setChecked(
                config.getboolean("Openbabel", 'AtomsColor'))
            self.radioAll.setChecked(
                config.getboolean("Openbabel", 'AtomsAll'))
            self.radioEnd.setChecked(
                config.getboolean("Openbabel", 'AtomsEnd'))
            self.radioNone.setChecked(
                config.getboolean("Openbabel", 'AtomsNone'))
            self.checkTighBond.setChecked(
                config.getboolean("Openbabel", 'TighBond'))

    def value(self, config):
        if not config.has_section("Openbabel"):
            config.add_section("Openbabel")
        config.set("Openbabel", "BondColor", self.BondColor.color.name())
        config.set("Openbabel", "BackColor",
                   str(self.BackgroundColor.color.name()))
        config.set("Openbabel", "BackColorAlpha",
                   str(self.BackgroundColor.color.alpha()))
        config.set("Openbabel", "AtomsColor", str(self.checkColor.isChecked()))
        config.set("Openbabel", "AtomsAll", str(self.radioAll.isChecked()))
        config.set("Openbabel", "AtomsEnd", str(self.radioEnd.isChecked()))
        config.set("Openbabel", "AtomsNone", str(self.radioNone.isChecked()))
        config.set("Openbabel", "TighBond",
                   str(self.checkTighBond.isChecked()))
        return config


class Preferences(QtWidgets.QDialog):
    """Preferences main dialog"""
    classes = [
        ("pychemqt.png", ConfGeneral,
         QtWidgets.QApplication.translate("pychemqt", "General")),
        ("button/PFD.png", prefPFD.Widget,
         QtWidgets.QApplication.translate("pychemqt", "PFD")),
        ("button/tooltip.png", ConfTooltipEntity,
         QtWidgets.QApplication.translate("pychemqt", "Tooltips in PFD")),
        ("button/tooltip.png", ConfTooltipUnit,
         QtWidgets.QApplication.translate("pychemqt", "Tooltips in units")),
        ("button/format_numeric.png", ConfFormat,
         QtWidgets.QApplication.translate("pychemqt", "Numeric format")),
        ("button/oil.png", prefPetro.Widget,
         QtWidgets.QApplication.translate("pychemqt", "Pseudocomponents")),
        ("button/applications.png", ConfApplications,
         QtWidgets.QApplication.translate("pychemqt", "Applications")),
        ("button/periodicTable.png", prefElemental.Widget,
         QtWidgets.QApplication.translate("pychemqt", "Elemental table")),
        ("button/tables.png", prefMEOS.Widget,
         QtWidgets.QApplication.translate("pychemqt", "mEoS")),
        ("button/psychrometric.png", prefPsychrometric.Widget,
         QtWidgets.QApplication.translate("pychemqt", "Psychrometric chart")),
        ("button/moody.png", moody.Config,
         QtWidgets.QApplication.translate("pychemqt", "Moody chart")),
        ("button/moody.png", drag.Config,
         QtWidgets.QApplication.translate("pychemqt", "Drag Sphere chart")),
        ("button/moody.png", standing.Config,
         QtWidgets.QApplication.translate("pychemqt", "Standing-Katz chart")),
        ("chemistry/grupo18.png", ConfBabel,
         QtWidgets.QApplication.translate("pychemqt", "Openbabel"))]

    def __init__(self, config, parent=None):
        """Constructor, opcional config parameter to input project config"""
        super(Preferences, self).__init__(parent)
        self.config = config
        self.setWindowTitle(
            QtWidgets.QApplication.translate("pychemqt", "Preferences"))

        layout = QtWidgets.QGridLayout(self)
        self.lista = QtWidgets.QListWidget()
        self.lista.setIconSize(QtCore.QSize(20, 20))
        self.lista.setSizePolicy(
            QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Preferred)
        layout.addWidget(self.lista, 1, 1)
        self.stacked = QtWidgets.QStackedWidget()
        layout.addWidget(self.stacked, 1, 2)
        for icon, dialog, title in self.classes:
            self.stacked.addWidget(dialog(config))
            icon = QtGui.QIcon(QtGui.QPixmap(
                os.environ["pychemqt"]+"/images/%s" % icon))
            self.lista.addItem(QtWidgets.QListWidgetItem(icon, title))

        self.lista.currentRowChanged.connect(self.stacked.setCurrentIndex)
        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Cancel | QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox, 2, 2)

    def value(self):
        """Return value for wizard"""
        config = self.config
        for indice in range(self.stacked.count()):
            config = self.stacked.widget(indice).value(config)
        return config


if __name__ == "__main__":
    conf_dir = os.path.expanduser('~') + "/.pychemqt/"
    pychemqt_dir = os.environ["PWD"] + "/"
    app = QtWidgets.QApplication(sys.argv)

#    config={"format": 0, "decimales": 4, "signo": False}
#    dialogo=NumericFactor(config)

    config = ConfigParser()
    config.read(conf_dir+"pychemqtrc")
    dialogo = Preferences(config)

    dialogo.show()
    sys.exit(app.exec_())
