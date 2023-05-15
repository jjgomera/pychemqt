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
along with this program.  If not, see <http://www.gnu.org/licenses/>.


Library to show/configure pychemqt general preferences

    * :class:`ConfGeneral`: General configuration options
    * :class:`ConfTooltipUnit`: Tooltip with unit alternate value configuration
    * :class:`ConfFormat`: Numeric format configuration
    * :class:`ConfApplications`: External applications configuration
    * :class:`ConfTooltipEntity`: Entity properties in popup configuration
'''


from ast import literal_eval
from configparser import ConfigParser
from functools import partial
from math import pi
import os
import sys

from equipment import equipments
from lib import unidades, corriente
from lib.openbabel import ConfBabel
from lib.plot import ConfPlot
from lib.utilities import representacion
import plots
from tools.firstrun import which
from tools.qt import QtCore, QtGui, QtWidgets, tr
from tools.qtelemental import Config as ConfigElemental
from tools.UI_psychrometry import Config as ConfigPsychrometry
from tools.UI_Tables import prefMEOS
from UI import prefPFD, prefPetro
from UI.delegate import CheckEditor
from UI.widgets import ColorSelector, NumericFactor, PathConfig


class ConfGeneral(QtWidgets.QDialog):
    """General configuration options"""
    def __init__(self, config=None, parent=None):
        super().__init__(parent)

        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Highlight color:")), 1, 1)
        self.ColorButtonResaltado = ColorSelector()
        layout.addWidget(self.ColorButtonResaltado, 1, 2)
        layout.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Readonly color:")), 2, 1)
        self.ColorButtonReadOnly = ColorSelector()
        layout.addWidget(self.ColorButtonReadOnly, 2, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 0, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 3, 1)
        group = QtWidgets.QGroupBox(
            tr("pychemqt", "Recent Files"))
        layout.addWidget(group, 4, 1, 1, 4)
        lyt = QtWidgets.QHBoxLayout(group)
        lyt.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Number of recent files:")))
        self.recentFiles = QtWidgets.QSpinBox()
        self.recentFiles.setRange(1, 20)
        lyt.addWidget(self.recentFiles)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed))
        self.loadLastProject = QtWidgets.QCheckBox(
            tr(
                "pychemqt", "Load last session project"))
        layout.addWidget(self.loadLastProject, 5, 1)
        self.showTrayIcon = QtWidgets.QCheckBox(
            tr("pychemqt", "Show tray icon"))
        layout.addWidget(self.showTrayIcon, 6, 1)

        layout.addItem(QtWidgets.QSpacerItem(
            10, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 14, 1, 1, 4)

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
        """Update ConfigParser instance with the config"""
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
        super().__init__(parent)

        layout = QtWidgets.QVBoxLayout(self)
        self.checkShow = QtWidgets.QCheckBox(
            tr("pychemqt", "Show Tool Tips"))
        self.checkShow.toggled.connect(self.checkShow_Toggled)
        layout.addWidget(self.checkShow)

        self.groupsystems = QtWidgets.QGroupBox(
            tr(
                "pychemqt", "Systems of measurement"))
        layout.addWidget(self.groupsystems)
        lytSystems = QtWidgets.QHBoxLayout(self.groupsystems)
        self.SI = QtWidgets.QCheckBox(
            tr("pychemqt", "SI"))
        self.SI.toggled.connect(partial(self.systems, "si"))
        lytSystems.addWidget(self.SI)
        self.AltSI = QtWidgets.QCheckBox(
            tr("pychemqt", "Alt SI"))
        self.AltSI.toggled.connect(partial(self.systems, "altsi"))
        lytSystems.addWidget(self.AltSI)
        self.English = QtWidgets.QCheckBox(
            tr("pychemqt", "English"))
        self.English.toggled.connect(partial(self.systems, "english"))
        lytSystems.addWidget(self.English)
        self.Metric = QtWidgets.QCheckBox(
            tr("pychemqt", "Metric"))
        self.Metric.toggled.connect(partial(self.systems, "metric"))
        lytSystems.addWidget(self.Metric)
        self.CGS = QtWidgets.QCheckBox(
            tr("pychemqt", "CGS"))
        self.CGS.toggled.connect(partial(self.systems, "cgs"))
        lytSystems.addWidget(self.CGS)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed))

        self.eleccion = QtWidgets.QComboBox()
        layout.addWidget(self.eleccion)
        self.stacked = QtWidgets.QStackedWidget()
        self.eleccion.currentIndexChanged.connect(self.stacked.setCurrentIndex)
        layout.addWidget(self.stacked)

        self.tabla = []
        for magnitud in unidades.MAGNITUDES[:-1]:
            # Discard last magnitud, dimensionless
            textos = magnitud[2].__text__
            tabla = QtWidgets.QTableWidget()
            self.stacked.addWidget(tabla)

            tabla.setRowCount(len(textos))
            tabla.setColumnCount(1)
            tabla.setColumnWidth(0, 16)
            tabla.setItemDelegateForColumn(0, CheckEditor(self))
            tabla.horizontalHeader().setVisible(False)

            for j, txt in enumerate(textos):
                tabla.setRowHeight(j, 24)
                tabla.setItem(j, 0, QtWidgets.QTableWidgetItem(""))
                tabla.item(j, 0).setTextAlignment(
                    QtCore.Qt.AlignmentFlag.AlignRight
                    | QtCore.Qt.AlignmentFlag.AlignVCenter)
                tabla.openPersistentEditor(tabla.item(j, 0))

                header = QtWidgets.QTableWidgetItem(txt)
                if magnitud[0] == "Currency":
                    header.setIcon(QtGui.QIcon(QtGui.QPixmap(os.path.join(
                        os.environ["pychemqt"], "images", "flag",
                        "%s.png" % unidades.Currency.__units__[j]))))
                    header.setToolTip(unidades.Currency.__tooltip__[j])

                tabla.setVerticalHeaderItem(j, header)

            self.fill(magnitud[0], tabla, config)
            self.tabla.append(tabla)
            self.eleccion.addItem(magnitud[1])

        if config.has_section("Tooltip"):
            self.checkShow.setChecked(config.getboolean("Tooltip", "Show"))
            self.SI.setChecked(config.getboolean("Tooltip", "SI"))
            self.AltSI.setChecked(config.getboolean("Tooltip", "AltSI"))
            self.English.setChecked(config.getboolean("Tooltip", "English"))
            self.Metric.setChecked(config.getboolean("Tooltip", "Metric"))
            self.CGS.setChecked(config.getboolean("Tooltip", "CGS"))

    @staticmethod
    def fill(magnitud, tabla, config):
        """Fill the the table with the data"""
        if config.has_section("Tooltip"):
            lst = map(int, config.get("Tooltip", magnitud).split(","))
            for i in lst:
                tabla.item(i, 0).setText("true")

    def checkShow_Toggled(self, boolean):
        "Enable/Disable widgets"""
        self.eleccion.setEnabled(boolean)
        self.groupsystems.setEnabled(boolean)
        for tabla in self.tabla:
            tabla.setEnabled(boolean)

    def systems(self, unitSet, boolean):
        """Check the units defined in the unit set"""
        if boolean:
            txt = "true"
        else:
            txt = ""
        for tabla, value in enumerate(unidades.units_set[unitSet][:-1]):
            self.tabla[tabla].item(value, 0).setText(txt)

    def value(self, config):
        """Update ConfigParser instance with the config"""
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
            config.set("Tooltip", unidades.MAGNITUDES[i][0],
                       ", ".join(map(str, lista)))
        return config


class ConfTooltipEntity(QtWidgets.QDialog):
    """Entity properties in popup window configuration"""
    def __init__(self, config, parent=None):
        super().__init__(parent)

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
            tr("pychemqt", "Stream"))

        for i, propiedad in enumerate(corriente.Corriente.propertiesNames()):
            item = QtWidgets.QTableWidgetItem(propiedad[0])
            self.tabla[0].setVerticalHeaderItem(i, item)
            self.tabla[0].setRowHeight(i, 24)
            self.tabla[0].setItem(i, 0, QtWidgets.QTableWidgetItem(""))
            self.tabla[0].item(i, 0).setTextAlignment(
                QtCore.Qt.AlignmentFlag.AlignRight
                | QtCore.Qt.AlignmentFlag.AlignVCenter)
            self.tabla[0].openPersistentEditor(self.tabla[0].item(i, 0))

        if config.has_option("TooltipEntity", "Corriente"):
            lst = map(int, config.get("TooltipEntity", "Corriente").split(","))
            for i in lst:
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
                    QtCore.Qt.AlignmentFlag.AlignRight
                    | QtCore.Qt.AlignmentFlag.AlignVCenter)
                self.tabla[-1].openPersistentEditor(self.tabla[-1].item(j, 0))
            self.fill(equipo.__name__, i+1, config)
            self.eleccion.addItem(equipo.title)

    def fill(self, equipo, tabla, config):
        """Fill the the table with the data"""
        if config.has_section("TooltipEntity"):
            lst = map(int, config.get("TooltipEntity", equipo).split(","))
            for i in lst:
                self.tabla[tabla].item(i, 0).setText("true")

    def value(self, config):
        """Update ConfigParser instance with the config"""
        if not config.has_section("TooltipEntity"):
            config.add_section("TooltipEntity")
        lista = []
        for j in range(self.tabla[0].rowCount()):
            if self.tabla[0].item(j, 0).text() == "true":
                lista.append(j)
        config.set("TooltipEntity", "Corriente", ", ".join(map(str, lista)))

        for i, tabla in enumerate(self.tabla[1:]):
            lista = []
            for j in range(tabla.rowCount()):
                if tabla.item(j, 0).text() == "true":
                    lista.append(j)
            config.set("TooltipEntity", equipments[i].__name__,
                       ", ".join(map(str, lista)))

        return config


class ConfFormat(QtWidgets.QTableWidget):
    """Numeric format configuration"""
    def __init__(self, config=None, parent=None):
        super().__init__(parent)
        self.setColumnCount(2)
        self.setRowCount(len(unidades.MAGNITUDES))
        labels = [tr("pychemqt", "Format"),
                  tr("pychemqt", "Sample")]
        self.setHorizontalHeaderLabels(labels)
        self.config = []
        for i in range(self.rowCount()):
            item = QtWidgets.QTableWidgetItem(unidades.MAGNITUDES[i][1])
            self.setVerticalHeaderItem(i, item)
            self.setRowHeight(i, 22)
            self.setItem(i, 0, QtWidgets.QTableWidgetItem(""))
            self.item(i, 0).setTextAlignment(
                QtCore.Qt.AlignmentFlag.AlignRight
                | QtCore.Qt.AlignmentFlag.AlignVCenter)
            self.setItem(i, 1, QtWidgets.QTableWidgetItem(""))
            self.item(i, 1).setTextAlignment(
                QtCore.Qt.AlignmentFlag.AlignRight
                | QtCore.Qt.AlignmentFlag.AlignVCenter)

        if config.has_section("NumericFormat"):
            for i, magnitud in enumerate(unidades.MAGNITUDES):
                kw = literal_eval(config.get("NumericFormat", magnitud[0]))
                self.config.append(kw)
                self.item(i, 0).setText(self.txt(kw))
                self.item(i, 1).setText(representacion(pi, **kw))

        self.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.cellDoubleClicked.connect(self.showConfDialog)

    def showConfDialog(self, fila):
        """Show dialog with numeric format configuration"""
        dialog = NumericFactor(self.config[fila], parent=self)
        if dialog.exec():
            config = dialog.args()
            self.config[fila] = config
            self.item(fila, 0).setText(self.txt(config))
            self.item(fila, 1).setText(representacion(pi, **config))

    @staticmethod
    def txt(formato):
        """Get the string format text"""
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
        """Update ConfigParser instance with the config"""
        if not config.has_section("NumericFormat"):
            config.add_section("NumericFormat")
        for i, magnitud in enumerate(unidades.MAGNITUDES):
            config.set("NumericFormat", magnitud[0], str(self.config[i]))
        return config


class ConfApplications(QtWidgets.QDialog):
    """External applications configuration"""
    def __init__(self, config=None, parent=None):
        super().__init__(parent)
        layout = QtWidgets.QGridLayout(self)
        t = tr("pychemqt", "External Calculator")
        msg = tr(
            "pychemqt", "Select External Calculator Application")
        self.calculadora = PathConfig(t + ":", msg=msg, patron="exe")
        layout.addWidget(self.calculadora, 1, 1)
        t = tr("pychemqt", "Text viewer")
        msg = tr(
            "pychemqt", "Select External Text Viewer Application")
        self.textViewer = PathConfig(t + ":", msg=msg, patron="exe")
        layout.addWidget(self.textViewer, 2, 1)

        terminal = QtWidgets.QGroupBox()
        layout.addWidget(terminal, 3, 1)
        layoutTerminal = QtWidgets.QGridLayout(terminal)
        msg = tr(
            "pychemqt", "Select External shell")
        self.terminal = PathConfig("", msg=msg, patron="exe")
        layoutTerminal.addWidget(self.terminal, 1, 1, 1, 3)
        layoutTerminal.addWidget(QtWidgets.QLabel(
            tr("pychemqt", "Foreground color:")),
            2, 1)
        self.ForegroundColor = ColorSelector()
        layoutTerminal.addWidget(self.ForegroundColor, 2, 2, 1, 2)
        layoutTerminal.addWidget(QtWidgets.QLabel(
            tr("pychemqt", "Background color:")),
            3, 1)
        self.BackgroundColor = ColorSelector()
        layoutTerminal.addWidget(self.BackgroundColor, 3, 2, 1, 2)
        self.ipython = QtWidgets.QCheckBox(tr(
            "pychemqt", "Use ipython if its available"))
        layoutTerminal.addWidget(self.ipython, 4, 1, 1, 3)
        self.maximized = QtWidgets.QCheckBox(
            tr("pychemqt", "Show maximized"))
        layoutTerminal.addWidget(self.maximized, 5, 1, 1, 3)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 10, 1)

        terminalTitle = tr("pychemqt", "Shell")
        if sys.platform == "win32":
            terminal.setEnabled(False)
            terminalTitle += " (" + tr(
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

        # TODO: Enable option when add support for other terminals
        self.terminal.setEnabled(False)

    def value(self, config):
        """Update ConfigParser instance with the config"""
        if not config.has_section("Applications"):
            config.add_section("Applications")
        config.set("Applications", "Calculator", self.calculadora.text())
        config.set("Applications", "TextViewer", self.textViewer.text())
        config.set("Applications", "Shell", self.terminal.text())
        config.set("Applications", "ipython", str(self.ipython.isChecked()))
        config.set("Applications", "maximized",
                   str(self.maximized.isChecked()))
        config.set("Applications", "foregroundColor",
                   self.ForegroundColor.color.name())
        config.set("Applications", "backgroundColor",
                   self.BackgroundColor.color.name())
        return config


class Preferences(QtWidgets.QDialog):
    """Preferences main dialog"""
    classes = [
        ("pychemqt.png", ConfGeneral,
         tr("pychemqt", "General")),
        ("button/PFD.png", prefPFD.Widget,
         tr("pychemqt", "PFD")),
        ("button/tooltip.png", ConfTooltipEntity,
         tr("pychemqt", "Tooltips in PFD")),
        ("button/tooltip.png", ConfTooltipUnit,
         tr("pychemqt", "Tooltips in units")),
        ("button/format_numeric.png", ConfFormat,
         tr("pychemqt", "Numeric format")),
        ("button/oil.png", prefPetro.Widget,
         tr("pychemqt", "Pseudocomponents")),
        ("button/applications.png", ConfApplications,
         tr("pychemqt", "Applications")),
        ("button/plot.png", ConfPlot,
         tr("pychemqt", "Plotting")),
        ("button/periodicTable.png", ConfigElemental,
         tr("pychemqt", "Elemental table")),
        ("button/tables.png", prefMEOS.Widget,
         tr("pychemqt", "mEoS")),
        ("button/psychrometric.png", ConfigPsychrometry,
         tr("pychemqt", "Psychrometric chart")),

        ("button/moody.png", plots.Pref,
         tr("pychemqt", "Chart")),

        ("chemistry/grupo18.png", ConfBabel,
         tr("pychemqt", "Openbabel"))]

    def __init__(self, config=None, parent=None):
        """Constructor, opcional config parameter to input project config"""
        super().__init__(parent)

        if config is None:
            config = ConfigParser()
        self.config = config

        self.setWindowTitle(
            tr("pychemqt", "Preferences"))

        layout = QtWidgets.QGridLayout(self)
        self.lista = QtWidgets.QTreeWidget()
        self.lista.setIconSize(QtCore.QSize(30, 30))
        self.lista.setHeaderHidden(True)
        self.lista.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Maximum,
            QtWidgets.QSizePolicy.Policy.Preferred)
        layout.addWidget(self.lista, 1, 1)
        self.stacked = QtWidgets.QStackedWidget()
        layout.addWidget(self.stacked, 1, 2)

        count = 0
        for icon, dialog, title in self.classes:
            item = QtWidgets.QTreeWidgetItem([title])
            icon = QtGui.QIcon(QtGui.QPixmap(
                os.environ["pychemqt"]+"/images/%s" % icon))
            item.setIcon(0, icon)
            item.setData(0, QtCore.Qt.ItemDataRole.UserRole, count)
            count += 1
            self.lista.addTopLevelItem(item)

            if isinstance(dialog, tuple):
                self.stacked.addWidget(QtWidgets.QWidget())
                for dlg in dialog:
                    child = QtWidgets.QTreeWidgetItem([dlg.TITLE])
                    child.setIcon(0, icon)
                    child.setData(0, QtCore.Qt.ItemDataRole.UserRole, count)
                    item.addChild(child)
                    count += 1
                    self.stacked.addWidget(dlg(config))
            else:
                self.stacked.addWidget(dialog(config))

        self.lista.currentItemChanged.connect(self.getIndex)
        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Cancel
            | QtWidgets.QDialogButtonBox.StandardButton.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox, 2, 2)

    def getIndex(self, item):
        """Get index of item"""
        self.stacked.setCurrentIndex(
            item.data(0, QtCore.Qt.ItemDataRole.UserRole))

    def value(self):
        """Return value for wizard"""
        config = self.config
        for indice in range(self.stacked.count()):
            try:
                config = self.stacked.widget(indice).value(config)
            except AttributeError:
                pass
        return config


if __name__ == "__main__":
    conf_dir = os.path.expanduser('~') + "/.pychemqt/"
    pychemqt_dir = os.environ["PWD"] + "/"
    app = QtWidgets.QApplication(sys.argv)

    conf = ConfigParser()
    conf.read(conf_dir+"pychemqtrc")
    dialogo = Preferences(conf)

    dialogo.show()
    sys.exit(app.exec())
