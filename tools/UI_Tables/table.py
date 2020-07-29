#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2017, Juan José Gómez Romera <jjgomera@gmail.com>

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

from PyQt5 import QtCore, QtGui, QtWidgets
from numpy import delete, insert

from lib import meos, mEoS, coolProp, refProp, unidades, config
from lib.thermo import ThermoAdvanced
from lib.utilities import representacion, exportTable
from UI.widgets import (Entrada_con_unidades, createAction, Status, Tabla,
                        NumericFactor)

from tools.UI_Tables.plot import PlotMEoS
from tools.UI_Tables.library import get_propiedades, _getData, getClassFluid


# Table data
def createTabla(config, title, fluidos=None, parent=None):
    """Create TablaMEoS to add to mainwindow
        config: configparser instance with project configuration
        title: title for the table
        fluidos: optional array with meos instances to fill de table
        parent: mainwindow pointer
        """
    propiedades, keys, units = get_propiedades(config)

    # Add the unit suffix to properties title
    for i, unit in enumerate(units):
        sufx = unit.text()
        if not sufx:
            sufx = "[-]"
        propiedades[i] = propiedades[i]+os.linesep+sufx

    # Add two phases properties if requested
    if config.getboolean("MEoS", "phase"):
        for i in range(len(propiedades)-1, -1, -1):
            if keys[i] in ThermoAdvanced.propertiesPhase():
                txt = [propiedades[i]]
                prefix = QtWidgets.QApplication.translate("pychemqt", "Liquid")
                txt.append(prefix+os.linesep+propiedades[i])
                prefix = QtWidgets.QApplication.translate("pychemqt", "Vapour")
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
            fila = _getData(fluido, keys, config.getboolean("MEoS", "phase"))
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
            if config.getboolean("MEoS", "phase") and \
                    key in ThermoAdvanced.propertiesPhase():
                columnInput.append(True)
                columnInput.append(True)
        kw["columnReadOnly"] = columnInput
        kw["dinamica"] = True

        # Discard the keys from single phase state as input values
        if config.getboolean("MEoS", "phase"):
            for i in range(len(keys)-1, -1, -1):
                if keys[i] in ThermoAdvanced.propertiesPhase():
                    keys[i:i+1] = [keys[i], "", ""]

        tabla = TablaMEoS(len(propiedades), filas=1, **kw)

    prefix = QtWidgets.QApplication.translate("pychemqt", "Table")
    tabla.setWindowTitle(prefix+": "+title)
    tabla.resizeColumnsToContents()
    return tabla


class TablaMEoS(Tabla):
    """Tabla customize to show meos data, add context menu options, save and
    load support in project"""
    Plot = None
    icon = os.path.join(config.IMAGE_PATH, "button", "table.png")

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
                {"format": 1, "decimales": 6, "signo": False}]*args[0]

        super(TablaMEoS, self).__init__(*args, **kwargs)
        self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(self.icon)))
        self.horizontalHeader().setContextMenuPolicy(
            QtCore.Qt.CustomContextMenu)
        self.horizontalHeader().customContextMenuRequested.connect(
            self.hHeaderClicked)
        self.verticalHeader().setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.verticalHeader().customContextMenuRequested.connect(
            self.vHeaderClicked)
        self.itemSelectionChanged.connect(self.selectPoint)
        self.data = []
        if not self.readOnly:
            self.cellChanged.connect(self.calculatePoint)

    def closeEvent(self, event):
        self.parent.dirty[self.parent.idTab] = True
        self.parent.saveControl()

    def _getPlot(self):
        """Return plot if it's loaded"""
        # FIXME: This procedure detect the first PlotMeos window, correct or
        # incorrect
        if not self.Plot:
            windows = self.parent.centralwidget.currentWidget().subWindowList()
            for window in windows:
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
        if dialog.exec_():
            # Check unit change
            if unit != unidades.Dimensionless and \
                    dialog.unit.currentIndex() != self.orderUnit[col]:
                for i, fila in enumerate(self.data):
                    conf = unit.__units__[self.orderUnit[col]]
                    key = unit.__units__[dialog.unit.currentIndex()]
                    value = unit(fila[col], conf).__getattribute__(key)
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
        range = QtWidgets.QTableWidgetSelectionRange(
            0, col, self.rowCount()-1, col)
        self.setRangeSelected(range, True)

    def vHeaderClicked(self, position):
        """Show dialog to manage item in table"""
        row = self.verticalHeader().logicalIndexAt(position)
        rows = []
        for item in self.selectedItems():
            if item.row() not in rows:
                rows.append(item.row())
        rows.sort(reverse=True)

        actionCopy = createAction(
            QtWidgets.QApplication.translate("pychemqt", "&Copy"),
            slot=self.copy, shortcut=QtGui.QKeySequence.Copy,
            icon=os.environ["pychemqt"] +
            os.path.join("images", "button", "editCopy"),
            parent=self)
        if not self.selectedItems():
            actionCopy.setEnabled(False)

        actionDelete = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Delete Point"),
            icon=os.environ["pychemqt"]+"/images/button/editDelete",
            slot=partial(self.delete, rows), parent=self)
        if not rows:
            actionDelete.setEnabled(False)

        actionInsert = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Insert Point"),
            icon=os.environ["pychemqt"]+"/images/button/add",
            slot=partial(self.add, row), parent=self)

        menu = QtWidgets.QMenu()
        menu.addAction(actionCopy)
        menu.addSeparator()
        menu.addAction(actionDelete)
        menu.addAction(actionInsert)
        menu.exec_(self.mapToGlobal(position))

    def delete(self, rows):
        """Delete rows from table and for saved data"""
        self.parent.statusbar.showMessage(QtWidgets.QApplication.translate(
            "pychemqt", "Deleting point..."))
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
            pref = QtWidgets.QApplication.translate("pychemqt", "Table from")
            title = self.windowTitle().split(pref)[1][1:]
            for row in rows:
                if title == QtWidgets.QApplication.translate(
                        "pychemqt", "Melting Line"):
                    for x in ThermoAdvanced.propertiesKey():
                        del data["melting"][x][row]
                elif title == QtWidgets.QApplication.translate(
                        "pychemqt", "Sublimation Line"):
                    for x in ThermoAdvanced.propertiesKey():
                        del data["sublimation"][x][row]
                elif title == QtWidgets.QApplication.translate(
                        "pychemqt", "Saturation Line") or \
                        title == QtWidgets.QApplication.translate(
                            "pychemqt", "Liquid Saturation Line"):
                    for x in ThermoAdvanced.propertiesKey():
                        del data["saturation_0"][x][row]
                elif title == QtWidgets.QApplication.translate(
                        "pychemqt", "Vapor Saturation Line"):
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
        self.parent.statusbar.clearMessage()

    def add(self, row):
        """Add point to a table and to saved file"""
        pref = QtWidgets.QApplication.translate("pychemqt", "Table from ")
        if pref in self.windowTitle():
            title = self.windowTitle().split(pref)[1]
            melting = title == QtWidgets.QApplication.translate(
                "pychemqt", "Melting Line")
        else:
            melting = False

        dlg = AddPoint(self.Point._new(), melting, self.parent)
        if dlg.exec_():
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
                datatoTable.append(dlg.fluid.__getattribute__(plot.x).config())
                datatoTable.append(dlg.fluid.__getattribute__(plot.y).config())

            # Add point to table
            self.addRow(index=row)
            self.setRow(row, datatoTable)

            # Update verticalHeader
            for row in range(self.rowCount()):
                self.setVHeader(row)

            # Add point to data plot
            if plot is None:
                return

            data = plot._getData()
            if title == QtWidgets.QApplication.translate(
                    "pychemqt", "Melting Line"):
                for x in ThermoAdvanced.propertiesKey():
                    data["melting"][x].insert(
                        row, dlg.fluid.__getattribute__(x))
            elif title == QtWidgets.QApplication.translate(
                    "pychemqt", "Sublimation Line"):
                for x in ThermoAdvanced.propertiesKey():
                    data["sublimation"][x].insert(
                        row, dlg.fluid.__getattribute__(x))
            elif title == QtWidgets.QApplication.translate(
                    "pychemqt", "Saturation Line") or \
                    title == QtWidgets.QApplication.translate(
                        "pychemqt", "Liquid Saturation Line"):
                for x in ThermoAdvanced.propertiesKey():
                    data["saturation_0"][x].insert(
                        row, dlg.fluid.__getattribute__(x))
            elif title == QtWidgets.QApplication.translate(
                    "pychemqt", "Vapor Saturation Line"):
                for x in ThermoAdvanced.propertiesKey():
                    data["saturation_1"][x].insert(
                        row, dlg.fluid.__getattribute__(x))
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
                        row, dlg.fluid.__getattribute__(x))
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
                if line.get_label() == QtWidgets.QApplication.translate(
                        "pychemqt", "Selected Point"):
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
            label = QtWidgets.QApplication.translate(
                "pychemqt", "Selected Point")
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
        for column, data in enumerate(data):
            if isinstance(data, str):
                txt = data
            else:
                txt = representacion(data, **self.format[column])
            self.setValue(row, column, txt)
        self.resizeColumnsToContents()

        # Set calculate point readOnly
        if not self.readOnly:
            flags = QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable
            color = config.Preferences.get("General", 'Color_ReadOnly')
            for i, bool in enumerate(self.columnReadOnly):
                if not bool:
                    self.item(row, i).setFlags(flags)
                    self.item(row, i).setBackground(QtGui.QColor(color))
        self.blockSignals(False)

    def contextMenuEvent(self, event):
        """Show context menu over cell"""
        menu = QtWidgets.QMenu()
        actionCopy = createAction(
            QtWidgets.QApplication.translate("pychemqt", "&Copy"),
            slot=partial(self.copy, event), shortcut=QtGui.QKeySequence.Copy,
            icon=os.environ["pychemqt"] +
            os.path.join("images", "button", "editCopy"),
            parent=self)
        export = createAction(
            QtWidgets.QApplication.translate("pychemqt", "E&xport to csv"),
            self.exportCSV,
            icon=os.environ["pychemqt"] +
            os.path.join("images", "button", "export"),
            tip=QtWidgets.QApplication.translate(
                "pychemqt", "Export table to file"),
            parent=self)
        menu.addAction(actionCopy)
        menu.addSeparator()
        menu.addAction(export)
        menu.exec_(event.globalPos())

    def copy(self, event=None):
        """Copy selected values to clipboard"""
        txt = [w.text() for w in self.selectedItems()]
        QtWidgets.QApplication.clipboard().setText(" ".join(txt))

    def exportCSV(self):
        """Export data table as a csv file"""
        if self.parent.currentFilename:
            dir = os.path.dirname(str(self.parent.currentFilename))
        else:
            dir = "."

        pat = []
        pat.append(QtWidgets.QApplication.translate(
            "pychemqt", "CSV files") + " (*.csv)")
        if os.environ["ezodf"] == "True":
            pat.append(QtWidgets.QApplication.translate(
                "pychemqt", "Libreoffice spreadsheet files") + " (*.ods)")
        if os.environ["xlwt"] == "True":
            pat.append(QtWidgets.QApplication.translate(
                "pychemqt",
                "Microsoft Excel 97/2000/XP/2003 XML") + " (*.xls)")
        if os.environ["openpyxl"] == "True":
            pat.append(QtWidgets.QApplication.translate(
                "pychemqt", "Microsoft Excel 2007/2010 XML") + " (*.xlsx)")
        patron = ";;".join(pat)

        fname, ext = QtWidgets.QFileDialog.getSaveFileName(
            self, QtWidgets.QApplication.translate(
                "pychemqt", "Export table to file"), dir, patron)
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
        all = unidades._all
        all.append(unidades.Dimensionless)
        data["unit"] = [all.index(unit) for unit in self.units]

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
        all = unidades._all
        all.append(unidades.Dimensionless)
        units = [all[i] for i in data["unit"]]

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
        super(Ui_Saturation, self).__init__(parent)
        self.setWindowTitle(
            QtWidgets.QApplication.translate("pychemqt", "Saturation Table"))
        layout = QtWidgets.QGridLayout(self)

        gboxType = QtWidgets.QGroupBox(
            QtWidgets.QApplication.translate("pychemqt", "Interphase"))
        layout.addWidget(gboxType, 1, 1, 1, 2)
        layoutg1 = QtWidgets.QGridLayout(gboxType)
        self.VL = QtWidgets.QRadioButton(QtWidgets.QApplication.translate(
            "pychemqt", "Vapor-Liquid (boiling line)"))
        layoutg1.addWidget(self.VL, 1, 1)
        self.SL = QtWidgets.QRadioButton(QtWidgets.QApplication.translate(
            "pychemqt", "Solid-Liquid (melting line"))
        layoutg1.addWidget(self.SL, 2, 1)
        self.SV = QtWidgets.QRadioButton(QtWidgets.QApplication.translate(
            "pychemqt", "Solid-Vapor (Sublimation line)"))
        layoutg1.addWidget(self.SV, 3, 1)

        groupboxVariar = QtWidgets.QGroupBox(
            QtWidgets.QApplication.translate("pychemqt", "Change"))
        layout.addWidget(groupboxVariar, 1, 3, 1, 2)
        layoutg2 = QtWidgets.QGridLayout(groupboxVariar)
        self.VariarTemperatura = QtWidgets.QRadioButton(
            QtWidgets.QApplication.translate("pychemqt", "Temperature"))
        self.VariarTemperatura.toggled.connect(self.updateVar)
        layoutg2.addWidget(self.VariarTemperatura, 1, 1)
        self.VariarPresion = QtWidgets.QRadioButton(
            QtWidgets.QApplication.translate("pychemqt", "Pressure"))
        self.VariarPresion.toggled.connect(self.updateVar)
        layoutg2.addWidget(self.VariarPresion, 2, 1)
        self.VariarXconT = QtWidgets.QRadioButton(
            QtWidgets.QApplication.translate(
                "pychemqt", "Quality at fixed temperature"))
        self.VariarXconT.toggled.connect(self.updateVar)
        layoutg2.addWidget(self.VariarXconT, 3, 1)
        self.VariarXconP = QtWidgets.QRadioButton(
            QtWidgets.QApplication.translate(
                "pychemqt", "Quality at fixed pressure"))
        self.VariarXconP.toggled.connect(self.updateVar)
        layoutg2.addWidget(self.VariarXconP, 4, 1)

        line = QtWidgets.QFrame()
        line.setFrameShape(QtWidgets.QFrame.HLine)
        line.setFrameShadow(QtWidgets.QFrame.Sunken)
        layout.addWidget(line, 2, 1, 1, 4)

        self.labelFix = QtWidgets.QLabel()
        layout.addWidget(self.labelFix, 4, 3)
        self.variableFix = Entrada_con_unidades(float)
        layout.addWidget(self.variableFix, 4, 4)
        self.labelinicial = QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Initial"))
        layout.addWidget(self.labelinicial, 4, 1)
        self.Inicial = Entrada_con_unidades(float)
        layout.addWidget(self.Inicial, 4, 2)
        self.labelfinal = QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Final"))
        layout.addWidget(self.labelfinal, 5, 1)
        self.Final = Entrada_con_unidades(float)
        layout.addWidget(self.Final, 5, 2)
        self.labelincremento = QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Increment"))
        layout.addWidget(self.labelincremento, 6, 1)
        self.Incremento = Entrada_con_unidades(float)
        layout.addWidget(self.Incremento, 6, 2)

        buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
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

    def updateVar(self, bool):
        """Update input values units and text"""
        if bool:
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
                self.labelinicial.setText(QtWidgets.QApplication.translate(
                    "pychemqt", "Initial quality"))
                self.labelfinal.setText(QtWidgets.QApplication.translate(
                    "pychemqt", "Final quality"))

            elif self.sender() == self.VariarXconP:
                self.labelFix.setVisible(True)
                self.labelFix.setText(unidades.Pressure.__title__)
                self.variableFix.deleteLater()
                self.variableFix = Entrada_con_unidades(
                    unidades.Pressure, value=fix)
                self.layout().addWidget(self.variableFix, 4, 4)
                unidadVariable = float
                self.labelinicial.setText(QtWidgets.QApplication.translate(
                    "pychemqt", "Initial quality"))
                self.labelfinal.setText(QtWidgets.QApplication.translate(
                    "pychemqt", "Final quality"))

            elif self.sender() == self.VariarTemperatura:
                self.labelFix.setVisible(False)
                self.variableFix.setVisible(False)
                unidadVariable = unidades.Temperature
                self.labelinicial.setText(QtWidgets.QApplication.translate(
                    "pychemqt", "Initial temperature"))
                self.labelfinal.setText(QtWidgets.QApplication.translate(
                    "pychemqt", "Final temperature"))

            else:
                self.labelFix.setVisible(False)
                self.variableFix.setVisible(False)
                unidadVariable = unidades.Pressure
                self.labelinicial.setText(QtWidgets.QApplication.translate(
                    "pychemqt", "Initial pressure"))
                self.labelfinal.setText(QtWidgets.QApplication.translate(
                    "pychemqt", "Final pressure"))

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
        QtWidgets.QApplication.translate("pychemqt", "Temperature"),
        QtWidgets.QApplication.translate("pychemqt", "Pressure"),
        QtWidgets.QApplication.translate("pychemqt", "Density"),
        QtWidgets.QApplication.translate("pychemqt", "Volume"),
        QtWidgets.QApplication.translate("pychemqt", "Enthalpy"),
        QtWidgets.QApplication.translate("pychemqt", "Entropy"),
        QtWidgets.QApplication.translate("pychemqt", "Internal Energy")]
    unidades = [unidades.Temperature, unidades.Pressure, unidades.Density,
                unidades.SpecificVolume, unidades.Enthalpy,
                unidades.SpecificHeat, unidades.Enthalpy, float]
    keys = ["T", "P", "rho", "v", "h", "s", "u", "x"]

    def __init__(self, parent=None):
        super(Ui_Isoproperty, self).__init__(parent)
        self.setWindowTitle(QtWidgets.QApplication.translate(
            "pychemqt", "Specify Isoproperty Table"))
        layout = QtWidgets.QGridLayout(self)

        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Hold constant")), 1, 1)
        self.fix = QtWidgets.QComboBox()
        for propiedad in self.propiedades:
            self.fix.addItem(propiedad)
        self.fix.currentIndexChanged.connect(self.actualizarUI)
        layout.addWidget(self.fix, 1, 2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Vary")), 2, 1)
        self.vary = QtWidgets.QComboBox()
        self.vary.currentIndexChanged.connect(self.actualizarVariable)
        layout.addWidget(self.vary, 2, 2)

        line = QtWidgets.QFrame()
        line.setFrameShape(QtWidgets.QFrame.HLine)
        line.setFrameShadow(QtWidgets.QFrame.Sunken)
        layout.addWidget(line, 3, 1, 1, 2)

        self.labelFix = QtWidgets.QLabel()
        layout.addWidget(self.labelFix, 4, 1)
        self.variableFix = Entrada_con_unidades(float)
        layout.addWidget(self.variableFix, 4, 2)
        self.labelinicial = QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Initial"))
        layout.addWidget(self.labelinicial, 5, 1)
        self.Inicial = Entrada_con_unidades(float)
        layout.addWidget(self.Inicial, 5, 2)
        self.labelfinal = QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Final"))
        layout.addWidget(self.labelfinal, 6, 1)
        self.Final = Entrada_con_unidades(float)
        layout.addWidget(self.Final, 6, 2)
        self.labelincremento = QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Increment"))
        layout.addWidget(self.labelincremento, 7, 1)
        self.Incremento = Entrada_con_unidades(float)
        layout.addWidget(self.Incremento, 7, 2)

        buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)
        layout.addWidget(buttonBox, 10, 1, 1, 2)

        self.actualizarUI(0)

    def actualizarUI(self, indice):
        self.vary.clear()
        propiedades = self.propiedades[:]
        if indice <= 1:
            propiedades.append(QtWidgets.QApplication.translate(
                "pychemqt", "Quality"))
        del propiedades[indice]
        for propiedad in propiedades:
            self.vary.addItem(propiedad)
        self.labelFix.setText(self.propiedades[indice])
        self.variableFix.deleteLater()
        self.variableFix = Entrada_con_unidades(self.unidades[indice])
        self.layout().addWidget(self.variableFix, 4, 2)

    def actualizarVariable(self, indice):
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
        super(AddPoint, self).__init__(parent)
        self.setWindowTitle(
            QtWidgets.QApplication.translate("pychemqt", "Add Point to line"))
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
        layout.addWidget(self.status, i+1, 1, 1, 2)

        if isinstance(fluid, meos.MEoS) and fluid._melting:
            self.checkMelting = QtWidgets.QRadioButton(
                QtWidgets.QApplication.translate("pychemqt", "Melting Point"))
            self.checkMelting.setChecked(melting)
            layout.addWidget(self.checkMelting, i+2, 1, 1, 2)
            i += 1
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "To")), i+2, 1)
        self.To = Entrada_con_unidades(unidades.Temperature)
        self.To.valueChanged.connect(partial(self.update, "To"))
        layout.addWidget(self.To, i+2, 2)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "rhoo")), i+3, 1)
        self.rhoo = Entrada_con_unidades(unidades.Density)
        self.rhoo.valueChanged.connect(partial(self.update, "rhoo"))
        layout.addWidget(self.rhoo, i+3, 2)

        self.checkBelow = QtWidgets.QCheckBox(QtWidgets.QApplication.translate(
            "pychemqt", "Add below selected point"))
        layout.addWidget(self.checkBelow, i+4, 1, 1, 2)

        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Reset | QtWidgets.QDialogButtonBox.Ok |
            QtWidgets.QDialogButtonBox.Cancel)
        self.buttonBox.clicked.connect(self.click)
        layout.addWidget(self.buttonBox, i+5, 1, 1, 2)

    def click(self, button):
        """Manage mouse click event over buttonbox"""
        if QtWidgets.QDialogButtonBox.Reset == \
                self.buttonBox.standardButton(button):
            self.reset()
        elif QtWidgets.QDialogButtonBox.Ok == \
                self.buttonBox.standardButton(button):
            self.accept()
        elif QtWidgets.QDialogButtonBox.Cancel == \
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
        for key, input in zip(self.keys, self.Inputs):
            input.setValue(fluid.__getattribute__(key))
            if key in fluid.kwargs and \
                    fluid.kwargs[key] != fluid.__class__.kwargs[key]:
                input.setResaltado(True)
            else:
                input.setResaltado(False)
        self.blockSignals(False)

    def reset(self):
        """Reset dialog widgets to initial clear status"""
        self.fluid = self.fluid.__class__()
        self.status.setState(self.fluid.status, self.fluid.msg)
        self.rhoo.clear()
        self.To.clear()
        for input in self.Inputs:
            input.clear()
            input.setResaltado(False)


if __name__ == "__main__":
    import sys
    from lib.refProp import RefProp
    from lib.mEoS import H2O

    app = QtWidgets.QApplication(sys.argv)
    conf = config.getMainWindowConfig()

    fluid = RefProp(ids=[62])
    # fluid = H2O()

    SteamTables = AddPoint(fluid)

    SteamTables.show()
    sys.exit(app.exec_())
