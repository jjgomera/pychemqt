#!/usr/bin/python3
# -*- coding: utf-8 -*-
# pylint: disable=too-many-locals, protected-access

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
# Table functionality for plugin:
#   - TablaMEoS: Tabla subclass to show meos data, add context menu options
#   - Ui_Saturation: Dialog to define input for a two-phase table calculation
#   - Ui_Isoproperty: Dialog to define input for isoproperty table calculations
#   - AddPoint: Dialog to add new point to line2D
#   - createTabla: create TablaMEoS
###############################################################################


from configparser import ConfigParser
from functools import partial
from math import ceil, floor
import os

from numpy import delete, insert
from tools.qt import QtCore, QtGui, QtWidgets, translate

from lib import meos, mEoS, coolProp, unidades, config
from lib.thermo import ThermoAdvanced
from lib.utilities import representacion, exportTable
from UI.widgets import (Entrada_con_unidades, createAction, Status, Tabla,
                        NumericFactor, ClickableLabel)

from .chooseFluid import Dialog_InfoFluid
from .library import getClassFluid, getMethod, get_propiedades, _getData
from .plot import PlotMEoS


# Table data
def createTabla(conf, title, fluidos=None, parent=None):
    """Create TablaMEoS to add to mainwindow

    Parameters
    ----------
    conf : ConfigParser
        configparser instance with project configuration
    title : str
        title for the table
    fluidos : lib.MEoS
        optional array with meos instances to fill de table
    parent : object
        mainwindow pointer
    """
    propiedades, keys, units = get_propiedades(conf)

    # Add the unit suffix to properties title
    for i, unit in enumerate(units):
        sufx = unit.text()
        if not sufx:
            sufx = "[-]"
        propiedades[i] = propiedades[i]+os.linesep+sufx

    # Add two phases properties if requested
    if conf.getboolean("MEoS", "phase"):
        for i in range(len(propiedades)-1, -1, -1):
            if keys[i] in ThermoAdvanced.propertiesPhase():
                txt = [propiedades[i]]
                prefix = translate("UI_Tables", "Liquid")
                txt.append(prefix+os.linesep+propiedades[i])
                prefix = translate("UI_Tables", "Vapour")
                txt.append(prefix+os.linesep+propiedades[i])
                propiedades[i:i+1] = txt
                units[i:i+1] = [units[i]]*3

    # Define common argument for TableMEoS
    kw = {}
    kw["horizontalHeader"] = propiedades
    kw["stretch"] = False
    kw["units"] = units
    kw["parent"] = parent
    kw["keys"] = keys

    if fluidos:
        # Generate a readOnly table filled of data
        tabla = TablaMEoS(len(propiedades), readOnly=True, **kw)
        data = []
        for fluido in fluidos:
            fila = _getData(fluido, keys, conf.getboolean("MEoS", "phase"))
            data.append(fila)
        tabla.setData(data)

    else:
        # Generate a dinamic table empty
        columnInput = []
        for key in keys:
            if key in ["P", "T", "x", "rho", "v", "h", "s"]:
                columnInput.append(False)
            else:
                columnInput.append(True)
            if conf.getboolean("MEoS", "phase") and \
                    key in ThermoAdvanced.propertiesPhase():
                columnInput.append(True)
                columnInput.append(True)
        kw["columnReadOnly"] = columnInput
        kw["dinamica"] = True

        # Discard the keys from single phase state as input values
        if conf.getboolean("MEoS", "phase"):
            for i in range(len(keys)-1, -1, -1):
                if keys[i] in ThermoAdvanced.propertiesPhase():
                    keys[i:i+1] = [keys[i], "", ""]

        tabla = TablaMEoS(len(propiedades), filas=1, **kw)

    prefix = translate("UI_Tables", "Table")
    tabla.setWindowTitle(prefix+": "+title)
    tabla.resizeColumnsToContents()

    state = {}
    state["method"] = getMethod()
    state["fluid"] = conf.getint("MEoS", "fluid")
    tabla.changeStatusThermo(state)

    return tabla


class TablaMEoS(Tabla):
    """Tabla customize to show meos data, add context menu options, save and
    load support in project"""
    Plot = None
    icon = os.path.join(config.IMAGE_PATH, "button", "table.png")
    Point = None

    def __init__(self, *args, **kwargs):
        """Constructor with additional kwargs don't recognize in Tabla
        keys: array with keys properties
        units: array of unidades subclasses
        orderUnit: array of index of unit magnitud to show
        format: array of dict with numeric format
        """
        # Manage special parameter dont recognize in Tabla
        self.parent = kwargs.get("parent", None)
        if "keys" in kwargs:
            self.keys = kwargs["keys"]
            del kwargs["keys"]
        else:
            self.keys = []
        self.units = kwargs["units"]
        del kwargs["units"]
        if "orderUnit" in kwargs:
            self.orderUnit = kwargs["orderUnit"]
            del kwargs["orderUnit"]
        else:
            self.orderUnit = []
            for unit in self.units:
                if unit == unidades.Dimensionless:
                    self.orderUnit.append(0)
                else:
                    conf = self.parent.currentConfig
                    self.orderUnit.append(conf.getint('Units', unit.__name__))

        if "format" in kwargs:
            self.format = kwargs["format"]
            del kwargs["format"]
        else:
            self.format = [
                {"fmt": 1, "decimales": 6, "signo": False}]*args[0]

        super().__init__(*args, **kwargs)
        self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(self.icon)))
        self.horizontalHeader().setContextMenuPolicy(
            QtCore.Qt.ContextMenuPolicy.CustomContextMenu)
        self.horizontalHeader().customContextMenuRequested.connect(
            self.hHeaderClicked)
        self.verticalHeader().setContextMenuPolicy(
            QtCore.Qt.ContextMenuPolicy.CustomContextMenu)
        self.verticalHeader().customContextMenuRequested.connect(
            self.vHeaderClicked)
        self.itemSelectionChanged.connect(self.selectPoint)
        self.data = []
        if not self.readOnly:
            self.cellChanged.connect(self.calculatePoint)

        # Widgets to show in the statusbar of mainwindow
        self.statusWidget = []
        self.statusThermo = ClickableLabel()
        self.statusThermo.setFrameShape(QtWidgets.QFrame.Shape.WinPanel)
        self.statusThermo.setFrameShadow(QtWidgets.QFrame.Shadow.Sunken)
        self.statusThermo.clicked.connect(self.showFluid)
        self.statusWidget.append(self.statusThermo)

    def showFluid(self):
        """Show fluid info dialog"""
        dlg = Dialog_InfoFluid(self.Point.__class__)
        dlg.exec()

    def changeStatusThermo(self, conf):
        """Change text show in status of mainwindow"""
        fluid = getClassFluid(conf["method"], conf["fluid"])
        txt = f"{fluid.name} ({conf['method']})"
        self.statusThermo.setText(txt)

    def closeEvent(self, event):
        """Force project changes to save at exit"""
        self.parent.dirty[self.parent.idTab] = True
        self.parent.saveControl()

    def _getPlot(self):
        """Return plot if it's loaded"""
        # FIXME: This procedure detect the first PlotMeos window, correct or
        # incorrect
        if not self.Plot:
            wdws = self.parent.centralWidget().currentWidget().subWindowList()
            for window in wdws:
                widget = window.widget()
                if isinstance(widget, PlotMEoS):
                    self.Plot = widget
                    break
        return self.Plot

    def hHeaderClicked(self, event):
        """Show dialog to config format and unit"""
        col = self.horizontalHeader().logicalIndexAt(event)
        unit = self.units[col]
        dialog = NumericFactor(self.format[col], unit, self.orderUnit[col])
        if dialog.exec():
            # Check unit change
            if unit != unidades.Dimensionless and \
                    dialog.unit.currentIndex() != self.orderUnit[col]:
                for i, fila in enumerate(self.data):
                    conf = unit.__units__[self.orderUnit[col]]
                    key = unit.__units__[dialog.unit.currentIndex()]
                    value = getattr(unit(fila[col], conf), key)
                    self.data[i][col] = value
                self.orderUnit[col] = dialog.unit.currentIndex()
                txt = self.horizontalHeaderItem(
                    col).text().split(os.linesep)[0]
                txt += os.linesep+unit.__text__[dialog.unit.currentIndex()]
                self.setHorizontalHeaderItem(
                    col, QtWidgets.QTableWidgetItem(txt))

            # Check format change
            self.format[col] = dialog.args()
            self.setStr()
            self.resizeColumnToContents(col)
        rango = QtWidgets.QTableWidgetSelectionRange(
            0, col, self.rowCount()-1, col)
        self.setRangeSelected(rango, True)

    def vHeaderClicked(self, position):
        """Show dialog to manage item in table"""
        row = self.verticalHeader().logicalIndexAt(position)
        rows = []
        for item in self.selectedItems():
            if item.row() not in rows:
                rows.append(item.row())
        rows.sort(reverse=True)

        actionCopy = createAction(
            self.tr("&Copy"),
            slot=self.copy, shortcut=QtGui.QKeySequence.StandardKey.Copy,
            icon=os.path.join("button", "editCopy.png"),
            parent=self)
        if not self.selectedItems():
            actionCopy.setEnabled(False)

        actionDelete = createAction(
            self.tr("Delete Point"),
            icon=os.path.join("button", "editDelete.png"),
            slot=partial(self.delete, rows), parent=self)
        if not rows:
            actionDelete.setEnabled(False)

        actionInsert = createAction(
            self.tr("Insert Point"),
            icon=os.path.join("button", "add.png"),
            slot=partial(self.add, row), parent=self)

        menu = QtWidgets.QMenu()
        menu.addAction(actionCopy)
        menu.addSeparator()
        menu.addAction(actionDelete)
        menu.addAction(actionInsert)
        menu.exec(self.mapToGlobal(position))

    def delete(self, rows):
        """Delete rows from table and for saved data"""
        self.parent.statusBar().showMessage(self.tr("Deleting point..."))
        QtWidgets.QApplication.processEvents()

        # Delete point from table
        for row in rows:
            self.removeRow(row)
            delete(self.data, row)

        # Update verticalHeader
        for row in range(self.rowCount()):
            self.setVHeader(row)

        # Delete point from data plot
        plot = self._getPlot()
        if plot:
            data = plot._getData()
            pref = self.tr("Table from")
            title = self.windowTitle().split(pref)[1][1:]
            for row in rows:
                if title == self.tr("Melting Line"):
                    for x in ThermoAdvanced.propertiesKey():
                        del data["melting"][x][row]
                elif title == self.tr("Sublimation Line"):
                    for x in ThermoAdvanced.propertiesKey():
                        del data["sublimation"][x][row]
                elif title == self.tr("Saturation Line") or \
                        title == self.tr("Liquid Saturation Line"):
                    for x in ThermoAdvanced.propertiesKey():
                        del data["saturation_0"][x][row]
                elif title == self.tr("Vapor Saturation Line"):
                    for x in ThermoAdvanced.propertiesKey():
                        del data["saturation_1"][x][row]
                else:
                    units = {"P": unidades.Pressure,
                             "T": unidades.Temperature,
                             "h": unidades.Enthalpy,
                             "s": unidades.Enthalpy,
                             "v": unidades.SpecificVolume,
                             "rho": unidades.Density}
                    var = str(title.split(" = ")[0])
                    txt = title.split(" = ")[1]
                    unit = units[var]
                    value = float(txt.split(" ")[0])
                    stdValue = str(unit(value, "conf"))
                    for x in ThermoAdvanced.propertiesKey():
                        del data[var][stdValue][x][row]
            plot._saveData(data)

            # Delete point from data
            for line in plot.plot.ax.lines:
                if str(line.get_label()) == str(title):
                    xdata = line._x
                    ydata = line._y
                    for row in rows:
                        xdata = delete(xdata, row)
                        ydata = delete(ydata, row)
                    line.set_xdata(xdata)
                    line.set_ydata(ydata)
                    plot.plot.draw()
                    break
        self.parent.statusBar().clearMessage()

    def add(self, row):
        """Add point to a table and to saved file"""
        pref = self.tr("Table from ")
        if pref in self.windowTitle():
            title = self.windowTitle().split(pref)[1]
            melting = title == self.tr("Melting Line")
        else:
            melting = False

        dlg = AddPoint(self.Point._new(), melting, self.parent)
        if dlg.exec():
            self.blockSignals(True)
            if dlg.checkBelow.isChecked():
                row += 1

            plot = self.Plot
            if plot is None:
                plot = self._getPlot()

            if plot is None:
                # If table has no associated plot, define as normal point
                units = []
                for ui, order in zip(self.units, self.orderUnit):
                    if ui is unidades.Dimensionless:
                        units.append("")
                    else:
                        units.append(ui.__units__[order])
                phase = self.parent.currentConfig.getboolean("MEoS", "phase")
                datatoTable = _getData(dlg.fluid, self.keys, phase, units)
            else:
                # If table has a associated plot, use the values of that
                datatoTable = []
                datatoTable.append(getattr(dlg.fluid, plot.x).config())
                datatoTable.append(getattr(dlg.fluid, plot.y).config())

            # Add point to table
            self.addRow(index=row)
            self.setRow(row, datatoTable)

            # Update verticalHeader
            for title in range(self.rowCount()):
                self.setVHeader(title)

            # Add point to data plot
            if plot is None:
                return

            data = plot._getData()
            if title == self.tr("Melting Line"):
                for x in ThermoAdvanced.propertiesKey():
                    data["melting"][x].insert(
                        row, getattr(dlg.fluid, x))
            elif title == self.tr("Sublimation Line"):
                for x in ThermoAdvanced.propertiesKey():
                    data["sublimation"][x].insert(
                        row, getattr(dlg.fluid, x))
            elif title == self.tr("Saturation Line") or \
                    title == self.tr("Liquid Saturation Line"):
                for x in ThermoAdvanced.propertiesKey():
                    data["saturation_0"][x].insert(
                        row, getattr(dlg.fluid, x))
            elif title == self.tr("Vapor Saturation Line"):
                for x in ThermoAdvanced.propertiesKey():
                    data["saturation_1"][x].insert(
                        row, getattr(dlg.fluid, x))
            else:
                units = {"P": unidades.Pressure,
                         "T": unidades.Temperature,
                         "h": unidades.Enthalpy,
                         "s": unidades.Enthalpy,
                         "v": unidades.SpecificVolume,
                         "rho": unidades.Density}
                var = str(title.split(" = ")[0])
                txt = title.split(" = ")[1]
                unit = units[var]
                value = float(txt.split(" ")[0])
                stdValue = str(unit(value, "conf"))

                for x in ThermoAdvanced.propertiesKey():
                    data[var][stdValue][x].insert(
                        row, getattr(dlg.fluid, x))
            plot._saveData(data)

            # Add point to data
            for line in plot.plot.ax.lines:
                if str(line.get_label()) == str(title):
                    xdata = insert(line._x, row, datatoTable[0])
                    ydata = insert(line._y, row, datatoTable[1])
                    line.set_xdata(xdata)
                    line.set_ydata(ydata)
                    plot.plot.draw()
                    break

            self.blockSignals(False)

    def selectPoint(self):
        """Show selected point in table in asociated plot if exist"""
        if self.dinamica:
            return

        plot = self._getPlot()
        if plot:
            # Remove old selected point if exist
            for i, line in enumerate(plot.plot.ax.lines):
                if line.get_label() == self.tr("Selected Point"):
                    del line
                    del plot.plot.ax.lines[i]

            # Add new selected points
            x = []
            y = []
            for item in self.selectedItems():
                if item.column():
                    y.append(float(item.text()))
                else:
                    x.append(float(item.text()))
            label = self.tr("Selected Point")
            plot.plot.ax.plot(x, y, 'ro', label=label)
            plot.plot.draw()

    def calculatePoint(self, row, column):
        """Add new value to kwargs for point, and show properties if it is
        calculable
        row, column: index for modified cell in table"""
        txt = self.item(row, column).text()
        if not txt:
            return

        key = self.keys[column]
        unit = self.units[column]
        if unit is unidades.Dimensionless:
            value = float(self.item(row, column).text())
        else:
            data = float(self.item(row, column).text())
            value = unit(data, unit.__units__[self.orderUnit[column]])
        self.Point(**{key: value})

        # If the Point is calculated, get data
        if self.Point.status:
            units = []
            for ui, order in zip(self.units, self.orderUnit):
                if ui is unidades.Dimensionless:
                    units.append("")
                else:
                    units.append(ui.__units__[order])
            phase = self.parent.currentConfig.getboolean("MEoS", "phase")
            data = _getData(self.Point, self.keys, phase, units)
            self.setRow(row, data)
            self.Point = self.Point._new()

            self.addRow()
            self.setCurrentCell(row+1, column)

    def setData(self, data):
        """Override Tabla method to adapt functionality"""
        if self.readOnly:
            self.data = data
            self.setStr()
        else:
            for i, row in enumerate(data):
                self.setRow(i, row)
        self.resizeColumnsToContents()

    def setStr(self):
        """Add data as string to cell table"""
        for fila, array in enumerate(self.data):
            if fila >= self.rowCount():
                self.addRow()
            for columna, data in enumerate(array):
                if isinstance(data, str):
                    txt = data
                else:
                    txt = representacion(data, **self.format[columna])
                self.setValue(fila, columna, txt)

    def setRow(self, row, data):
        """Add data to a row"""
        self.blockSignals(True)
        self.data.insert(row, data)
        for column, dat in enumerate(data):
            if isinstance(dat, str):
                txt = dat
            else:
                txt = representacion(dat, **self.format[column])
            self.setValue(row, column, txt)
        self.resizeColumnsToContents()

        # Set calculate point readOnly
        if not self.readOnly:
            color = config.Preferences.get("General", 'Color_ReadOnly')
            for i, boolean in enumerate(self.columnReadOnly):
                if not boolean:
                    self.item(row, i).setFlags(
                        QtCore.Qt.ItemFlag.ItemIsEnabled
                        | QtCore.Qt.ItemFlag.ItemIsSelectable)
                    self.item(row, i).setBackground(QtGui.QColor(color))
        self.blockSignals(False)

    def contextMenuEvent(self, event):
        """Show context menu over cell"""
        menu = QtWidgets.QMenu()
        actionCopy = createAction(
            self.tr("&Copy"),
            slot=partial(self.copy, event),
            shortcut=QtGui.QKeySequence.StandardKey.Copy,
            icon=os.path.join("button", "editCopy.png"),
            parent=self)
        export = createAction(
            self.tr("E&xport to csv"),
            self.exportCSV,
            icon=os.path.join("button", "export.png"),
            tip=self.tr("Export table to file"),
            parent=self)
        menu.addAction(actionCopy)
        menu.addSeparator()
        menu.addAction(export)
        menu.exec(event.globalPos())

    def copy(self, event=None):
        """Copy selected values to clipboard"""
        txt = [w.text() for w in self.selectedItems()]
        QtWidgets.QApplication.clipboard().setText(" ".join(txt))

    def exportCSV(self):
        """Export data table as a csv file"""
        if self.parent.currentFilename:
            folder = os.path.dirname(str(self.parent.currentFilename))
        else:
            folder = "."

        pat = []
        pat.append(self.tr("CSV files") + " (*.csv)")
        if os.environ["ezodf"] == "True":
            pat.append(self.tr("Libreoffice spreadsheet files") + " (*.ods)")
        if os.environ["xlwt"] == "True":
            pat.append(self.tr("Microsoft Excel 97/2000/XP/2003 XML")
                       + " (*.xls)")
        if os.environ["openpyxl"] == "True":
            pat.append(self.tr("Microsoft Excel 2007/2010 XML") + " (*.xlsx)")
        patron = ";;".join(pat)

        fname, ext = QtWidgets.QFileDialog.getSaveFileName(
            self, self.tr("Export table to file"), folder, patron)
        if fname and ext:
            ext = ext.split(".")[-1][:-1]
            exportTable(self.data, fname, ext, self.horizontalHeaderLabel)

    def writeToJSON(self, data):
        """Write instance parameter to file"""
        data["column"] = self.columnCount()

        # Save titles
        data["title"] = self.windowTitle()
        data["htitle"] = []
        for column in range(data["column"]):
            data["htitle"].append(self.horizontalHeaderItem(column).text())

        # Save units as index
        units = unidades._all
        units.append(unidades.Dimensionless)
        data["unit"] = [units.index(unit) for unit in self.units]

        # Save method calculation options
        if isinstance(self.Point, meos.MEoS):
            data["method"] = "meos"
            data["fluid"] = mEoS.__all__.index(self.Point.__class__)
            data["external_dependences"] = ""
        elif isinstance(self.Point, coolProp.CoolProp):
            data["method"] = "coolprop"
            data["fluid"] = mEoS.id_mEoS.index(self.Point.kwargs["ids"][0])
            data["external_dependences"] = "CoolProp"
        else:
            data["method"] = "refprop"
            data["fluid"] = mEoS.id_mEoS.index(self.Point.kwargs["ids"][0])
            data["external_dependences"] = "refprop"

        # Save keys
        data["keys"] = self.keys
        data["readOnly"] = self.readOnly
        data["dinamica"] = self.dinamica
        if not self.readOnly:
            data["columnReadOnly"] = self.columnReadOnly

        # Save order unit
        data["orderUnit"] = self.orderUnit

        # Save format
        data["format"] = self.format

        # Save data
        data["data"] = list(self.data)

    @classmethod
    def readFromJSON(cls, data, parent):
        """Load data table from saved file"""

        # Get units
        units = unidades._all
        units.append(unidades.Dimensionless)
        units = [units[i] for i in data["unit"]]

        # Create Tabla
        kwargs = {}
        kwargs["horizontalHeader"] = data["htitle"]
        kwargs["format"] = data["format"]
        kwargs["stretch"] = False
        kwargs["parent"] = parent
        kwargs["units"] = units
        kwargs["orderUnit"] = data["orderUnit"]
        kwargs["keys"] = data["keys"]
        kwargs["dinamica"] = data["dinamica"]

        if data["readOnly"]:
            kwargs["readOnly"] = True
        else:
            kwargs["filas"] = len(data["data"])+1
            kwargs["columnReadOnly"] = data["columnReadOnly"]

        tabla = TablaMEoS(data["column"], **kwargs)
        tabla.setWindowTitle(data["title"])
        tabla.setData(data["data"])

        fluid = getClassFluid(data["method"], data["fluid"])
        tabla.Point = fluid
        tabla.changeStatusThermo(data)

        return tabla


class Ui_Saturation(QtWidgets.QDialog):
    """Dialog to define input for a two-phase saturation table calculation"""
    def __init__(self, method=None, fluid=None, parent=None):
        """
        Parameters
        ----------
        method: str
            name of method of calculation, meos, coolprop or refprop
        fluid: int
            Index of fluid in list
        """
        super().__init__(parent)
        self.setWindowTitle(self.tr("Saturation Table"))
        layout = QtWidgets.QGridLayout(self)

        gboxType = QtWidgets.QGroupBox(self.tr("Interphase"))
        layout.addWidget(gboxType, 1, 1, 1, 2)
        layoutg1 = QtWidgets.QGridLayout(gboxType)
        self.VL = QtWidgets.QRadioButton(
            self.tr("Vapor-Liquid (boiling line)"))
        layoutg1.addWidget(self.VL, 1, 1)
        self.SL = QtWidgets.QRadioButton(self.tr("Solid-Liquid (melting line"))
        layoutg1.addWidget(self.SL, 2, 1)
        self.SV = QtWidgets.QRadioButton(
            self.tr("Solid-Vapor (Sublimation line)"))
        layoutg1.addWidget(self.SV, 3, 1)

        groupboxVariar = QtWidgets.QGroupBox(self.tr("Change"))
        layout.addWidget(groupboxVariar, 1, 3, 1, 2)
        layoutg2 = QtWidgets.QGridLayout(groupboxVariar)
        self.VariarTemperatura = QtWidgets.QRadioButton(self.tr("Temperature"))
        self.VariarTemperatura.toggled.connect(self.updateVar)
        layoutg2.addWidget(self.VariarTemperatura, 1, 1)
        self.VariarPresion = QtWidgets.QRadioButton(self.tr("Pressure"))
        self.VariarPresion.toggled.connect(self.updateVar)
        layoutg2.addWidget(self.VariarPresion, 2, 1)
        self.VariarXconT = QtWidgets.QRadioButton(
            self.tr("Quality at fixed temperature"))
        self.VariarXconT.toggled.connect(self.updateVar)
        layoutg2.addWidget(self.VariarXconT, 3, 1)
        self.VariarXconP = QtWidgets.QRadioButton(
            self.tr("Quality at fixed pressure"))
        self.VariarXconP.toggled.connect(self.updateVar)
        layoutg2.addWidget(self.VariarXconP, 4, 1)

        line = QtWidgets.QFrame()
        line.setFrameShape(QtWidgets.QFrame.Shape.HLine)
        line.setFrameShadow(QtWidgets.QFrame.Shadow.Sunken)
        layout.addWidget(line, 2, 1, 1, 4)

        self.labelFix = QtWidgets.QLabel()
        layout.addWidget(self.labelFix, 4, 3)
        self.variableFix = Entrada_con_unidades(float)
        layout.addWidget(self.variableFix, 4, 4)
        self.labelinicial = QtWidgets.QLabel(self.tr("Initial"))
        layout.addWidget(self.labelinicial, 4, 1)
        self.Inicial = Entrada_con_unidades(float)
        layout.addWidget(self.Inicial, 4, 2)
        self.labelfinal = QtWidgets.QLabel(self.tr("Final"))
        layout.addWidget(self.labelfinal, 5, 1)
        self.Final = Entrada_con_unidades(float)
        layout.addWidget(self.Final, 5, 2)
        self.labelincremento = QtWidgets.QLabel(self.tr("Increment"))
        layout.addWidget(self.labelincremento, 6, 1)
        self.Incremento = Entrada_con_unidades(float)
        layout.addWidget(self.Incremento, 6, 2)

        buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok
            | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)
        layout.addWidget(buttonBox, 10, 1, 1, 4)

        if method:
            self.fluido = getClassFluid(method, fluid)
            if isinstance(self.fluido, meos.MEoS) and (
                self.fluido._Melting_Pressure != meos.MEoS._Melting_Pressure or
                    self.fluido._melting):
                self.SL.setEnabled(True)
            else:
                self.SL.setEnabled(False)

            if isinstance(self.fluido, meos.MEoS) and (
                self.fluido._sublimation or
                self.fluido._Sublimation_Pressure !=
                    meos.MEoS._Sublimation_Pressure):
                self.SV.setEnabled(True)
            else:
                self.SV.setEnabled(False)

        self.VL.setChecked(True)
        self.VariarTemperatura.setChecked(True)
        self.updateVary()
        self.VL.toggled.connect(self.updateVary)

    def updateVary(self):
        """Update state for option to choose for properties to change"""
        self.VariarXconP.setEnabled(self.VL.isChecked())
        self.VariarXconT.setEnabled(self.VL.isChecked())
        self.VariarTemperatura.setChecked(not self.VL.isChecked())

    def updateVar(self, boolean):
        """Update input values units and text"""
        if boolean:
            # Select initial values
            fix, inicial, final, step = 0, 0, 0, 0
            if self.VL.isChecked():
                if self.sender() == self.VariarXconT:
                    fix = ceil((self.fluido.Tc-self.fluido.Tt)/2)
                    inicial = 0
                    final = 1
                    step = 0.1
                elif self.sender() == self.VariarXconP:
                    fix = ceil(self.fluido.Pc/2)
                    inicial = 0
                    final = 1
                    step = 0.1
                elif self.sender() == self.VariarTemperatura:
                    inicial = ceil(self.fluido.Tt)
                    final = floor(self.fluido.Tc)
                    step = 1.

            self.Inicial.deleteLater()
            self.Final.deleteLater()
            self.Incremento.deleteLater()
            if self.sender() == self.VariarXconT:
                self.labelFix.setVisible(True)
                self.labelFix.setText(unidades.Temperature.__title__)
                self.variableFix.deleteLater()
                self.variableFix = Entrada_con_unidades(
                    unidades.Temperature, value=fix)
                self.layout().addWidget(self.variableFix, 4, 4)
                unidadVariable = float
                self.labelinicial.setText(self.tr("Initial quality"))
                self.labelfinal.setText(self.tr("Final quality"))

            elif self.sender() == self.VariarXconP:
                self.labelFix.setVisible(True)
                self.labelFix.setText(unidades.Pressure.__title__)
                self.variableFix.deleteLater()
                self.variableFix = Entrada_con_unidades(
                    unidades.Pressure, value=fix)
                self.layout().addWidget(self.variableFix, 4, 4)
                unidadVariable = float
                self.labelinicial.setText(self.tr("Initial quality"))
                self.labelfinal.setText(self.tr("Final quality"))

            elif self.sender() == self.VariarTemperatura:
                self.labelFix.setVisible(False)
                self.variableFix.setVisible(False)
                unidadVariable = unidades.Temperature
                self.labelinicial.setText(self.tr("Initial temperature"))
                self.labelfinal.setText(self.tr("Final temperature"))

            else:
                self.labelFix.setVisible(False)
                self.variableFix.setVisible(False)
                unidadVariable = unidades.Pressure
                self.labelinicial.setText(self.tr("Initial pressure"))
                self.labelfinal.setText(self.tr("Final pressure"))

            self.Inicial = Entrada_con_unidades(unidadVariable, value=inicial)
            self.Final = Entrada_con_unidades(unidadVariable, value=final)
            if unidadVariable == unidades.Temperature:
                unidadDelta = unidades.DeltaT
            elif unidadVariable == unidades.Pressure:
                unidadDelta = unidades.DeltaP
            else:
                unidadDelta = unidadVariable

            self.Incremento = Entrada_con_unidades(unidadDelta, value=step)
            self.layout().addWidget(self.Inicial, 4, 2)
            self.layout().addWidget(self.Final, 5, 2)
            self.layout().addWidget(self.Incremento, 6, 2)


class Ui_Isoproperty(QtWidgets.QDialog):
    """Dialog to define input for isoproperty table calculations"""
    propiedades = [
        translate("UI_Tables", "Temperature"),
        translate("UI_Tables", "Pressure"),
        translate("UI_Tables", "Density"),
        translate("UI_Tables", "Volume"),
        translate("UI_Tables", "Enthalpy"),
        translate("UI_Tables", "Entropy"),
        translate("UI_Tables", "Internal Energy")]
    unidades = [unidades.Temperature, unidades.Pressure, unidades.Density,
                unidades.SpecificVolume, unidades.Enthalpy,
                unidades.SpecificHeat, unidades.Enthalpy, float]
    keys = ["T", "P", "rho", "v", "h", "s", "u", "x"]

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(self.tr("Specify Isoproperty Table"))
        layout = QtWidgets.QGridLayout(self)

        layout.addWidget(QtWidgets.QLabel(self.tr("Hold constant")), 1, 1)
        self.fix = QtWidgets.QComboBox()
        for propiedad in self.propiedades:
            self.fix.addItem(propiedad)
        self.fix.currentIndexChanged.connect(self.actualizarUI)
        layout.addWidget(self.fix, 1, 2)
        layout.addWidget(QtWidgets.QLabel(self.tr("Vary")), 2, 1)
        self.vary = QtWidgets.QComboBox()
        self.vary.currentIndexChanged.connect(self.actualizarVariable)
        layout.addWidget(self.vary, 2, 2)

        line = QtWidgets.QFrame()
        line.setFrameShape(QtWidgets.QFrame.Shape.HLine)
        line.setFrameShadow(QtWidgets.QFrame.Shadow.Sunken)
        layout.addWidget(line, 3, 1, 1, 2)

        self.labelFix = QtWidgets.QLabel()
        layout.addWidget(self.labelFix, 4, 1)
        self.variableFix = Entrada_con_unidades(float)
        layout.addWidget(self.variableFix, 4, 2)
        self.labelinicial = QtWidgets.QLabel(self.tr("Initial"))
        layout.addWidget(self.labelinicial, 5, 1)
        self.Inicial = Entrada_con_unidades(float)
        layout.addWidget(self.Inicial, 5, 2)
        self.labelfinal = QtWidgets.QLabel(self.tr("Final"))
        layout.addWidget(self.labelfinal, 6, 1)
        self.Final = Entrada_con_unidades(float)
        layout.addWidget(self.Final, 6, 2)
        self.labelincremento = QtWidgets.QLabel(self.tr("Increment"))
        layout.addWidget(self.labelincremento, 7, 1)
        self.Incremento = Entrada_con_unidades(float)
        layout.addWidget(self.Incremento, 7, 2)

        buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok
            | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)
        layout.addWidget(buttonBox, 10, 1, 1, 2)

        self.actualizarUI(0)

    def actualizarUI(self, indice):
        """Update UI"""
        self.vary.clear()
        propiedades = self.propiedades[:]
        if indice <= 1:
            propiedades.append(self.tr("Quality"))
        del propiedades[indice]
        for propiedad in propiedades:
            self.vary.addItem(propiedad)
        self.labelFix.setText(self.propiedades[indice])
        self.variableFix.deleteLater()
        self.variableFix = Entrada_con_unidades(self.unidades[indice])
        self.layout().addWidget(self.variableFix, 4, 2)

    def actualizarVariable(self, indice):
        """Update UI variables"""
        self.Inicial.deleteLater()
        self.Final.deleteLater()
        self.Incremento.deleteLater()
        if indice >= self.fix.currentIndex():
            indice += 1
        self.Inicial = Entrada_con_unidades(self.unidades[indice])
        self.Final = Entrada_con_unidades(self.unidades[indice])
        if self.unidades[indice] == unidades.Temperature:
            self.Incremento = Entrada_con_unidades(unidades.DeltaT)
        elif self.unidades[indice] == unidades.Pressure:
            self.Incremento = Entrada_con_unidades(unidades.DeltaP)
        else:
            self.Incremento = Entrada_con_unidades(self.unidades[indice])
        self.layout().addWidget(self.Inicial, 5, 2)
        self.layout().addWidget(self.Final, 6, 2)
        self.layout().addWidget(self.Incremento, 7, 2)


class AddPoint(QtWidgets.QDialog):
    """Dialog to add new point to line2D"""
    keys = ["T", "P", "x", "rho", "v", "h", "s", "u"]

    def __init__(self, fluid, melting=False, parent=None):
        """
        fluid: initial fluid instance
        melting: boolean to add melting line calculation
        """
        super().__init__(parent)
        self.setWindowTitle(self.tr("Add Point to line"))
        layout = QtWidgets.QGridLayout(self)
        self.fluid = fluid

        self.Inputs = []
        for i, (title, key, unit) in enumerate(meos.inputData):
            layout.addWidget(QtWidgets.QLabel(title), i, 1)
            if unit is unidades.Dimensionless:
                entrada = Entrada_con_unidades(float)
            else:
                entrada = Entrada_con_unidades(unit)
            entrada.valueChanged.connect(partial(self.update, key))
            self.Inputs.append(entrada)
            layout.addWidget(entrada, i, 2)

        self.status = Status(self.fluid.status, self.fluid.msg)
        row = len(meos.inputData)
        layout.addWidget(self.status, row+1, 1, 1, 2)

        if isinstance(fluid, meos.MEoS) and fluid._melting:
            self.checkMelting = QtWidgets.QRadioButton(
                self.tr("Melting Point"))
            self.checkMelting.setChecked(melting)
            layout.addWidget(self.checkMelting, row+2, 1, 1, 2)
            row += 1
        layout.addWidget(QtWidgets.QLabel(self.tr("To")), row+2, 1)
        self.To = Entrada_con_unidades(unidades.Temperature)
        self.To.valueChanged.connect(partial(self.update, "To"))
        layout.addWidget(self.To, row+2, 2)
        layout.addWidget(QtWidgets.QLabel(self.tr("rhoo")), row+3, 1)
        self.rhoo = Entrada_con_unidades(unidades.Density)
        self.rhoo.valueChanged.connect(partial(self.update, "rhoo"))
        layout.addWidget(self.rhoo, row+3, 2)

        self.checkBelow = QtWidgets.QCheckBox(
            self.tr("Add below selected point"))
        layout.addWidget(self.checkBelow, row+4, 1, 1, 2)

        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Reset
            | QtWidgets.QDialogButtonBox.StandardButton.Ok
            | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        self.buttonBox.clicked.connect(self.click)
        layout.addWidget(self.buttonBox, row+5, 1, 1, 2)

    def click(self, button):
        """Manage mouse click event over buttonbox"""
        if QtWidgets.QDialogButtonBox.StandardButton.Reset == \
                self.buttonBox.standardButton(button):
            self.reset()
        elif QtWidgets.QDialogButtonBox.StandardButton.Ok == \
                self.buttonBox.standardButton(button):
            self.accept()
        elif QtWidgets.QDialogButtonBox.StandardButton.Cancel == \
                self.buttonBox.standardButton(button):
            self.reject()

    def update(self, key, value):
        """Update fluid instance with new parameter key with value"""
        self.status.setState(4)
        QtWidgets.QApplication.processEvents()
        if isinstance(self.fluid, meos.MEoS) and self.fluid._melting and \
                self.checkMelting.isChecked() and key == "T":
            P = self.fluid._Melting_Pressure(value)
            self.fluid(**{key: value, "P": P})
        else:
            self.fluid(**{key: value})
        if self.fluid.status in (1, 3):
            self.fill(self.fluid)
        self.status.setState(self.fluid.status, self.fluid.msg)

    def fill(self, fluid):
        """Fill dialog widget with fluid properties values"""
        self.blockSignals(True)
        Config = ConfigParser()
        Config.read(config.conf_dir + "pychemqtrc")
        for key, data in zip(self.keys, self.Inputs):
            data.setValue(getattr(fluid, key))
            if key in fluid.kwargs and \
                    fluid.kwargs[key] != fluid.__class__.kwargs[key]:
                data.setResaltado(True)
            else:
                data.setResaltado(False)
        self.blockSignals(False)

    def reset(self):
        """Reset dialog widgets to initial clear status"""
        self.fluid = self.fluid.__class__()
        self.status.setState(self.fluid.status, self.fluid.msg)
        self.rhoo.clear()
        self.To.clear()
        for inputData in self.Inputs:
            inputData.clear()
            inputData.setResaltado(False)


if __name__ == "__main__":
    import sys
    from lib.mEoS import H2O

    app = QtWidgets.QApplication(sys.argv)
    cmp = H2O()
    SteamTables = AddPoint(cmp)
    SteamTables.show()
    sys.exit(app.exec())
