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
# Plot functionality for plugin:
#   - PlotMEoS: Plot widget to show meos plot data, add context menu options
#   - Plot2D: Dialog for select a special 2D plot
#   - Plot3D: Dialog for define a 3D plot
#   - EditPlot: Dialog to edit plot
#   - AddLine: Dialog to add new isoline to plot
#   - EditAxis: Dialog to configure axes plot properties
#   - AxisWidget: Dialog to configure axes plot properties
#   - calcIsoline: Isoline calculation procedure
#   - get_points: Get point number to plot lines from Preferences
#   - getLineFormat: get matplotlib line format from preferences
#   - plotIsoline: plot isoline procedure
#   - plot2D3D: general procedure for plotting 2D and 3D
#   - _getunitTransform
###############################################################################


from functools import partial
import gzip
from math import log10, atan, pi
import os
import pickle

from PyQt5 import QtCore, QtGui, QtWidgets
from numpy import concatenate, linspace, logspace, transpose, log, nan
from matplotlib.font_manager import FontProperties

from lib import meos, mEoS, unidades, plot, config
from lib.thermo import ThermoAdvanced
from lib.utilities import formatLine
from UI.widgets import (Entrada_con_unidades, createAction, LineStyleCombo,
                        MarkerCombo, ColorSelector, InputFont)

# from .table import TablaMEoS
from .library import calcPoint


class PlotMEoS(QtWidgets.QWidget):
    """Plot widget to show meos plot data, add context menu options"""
    def __init__(self, dim, toolbar=False, filename="", parent=None):
        """constructor
        Input:
            dim: dimension of plot, | 2 | 3
            toolbar: boolean to add the matplotlib toolbar
            filename: filename for data
        """
        super(PlotMEoS, self).__init__(parent)
        self.parent = parent
        self.dim = dim
        self.filename = filename
        self.notes = []

        layout = QtWidgets.QVBoxLayout(self)
        self.plot = plot.matplotlib(dim)

        self.plot.lx = self.plot.ax.axhline(c="#888888", ls=":")  # horiz line
        self.plot.ly = self.plot.ax.axvline(c="#888888", ls=":")  # vert line

        self.plot.lx.set_visible(False)
        self.plot.ly.set_visible(False)

        layout.addWidget(self.plot)
        self.toolbar = plot.NavigationToolbar2QT(self.plot, self.plot)
        self.toolbar.setVisible(toolbar)
        layout.addWidget(self.toolbar)

        self.editAxesAction = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Edit &Axis"),
            icon=os.environ["pychemqt"]+"/images/button/editor",
            slot=self.editAxis, parent=self)
        self.editAction = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Edit &Plot"),
            slot=self.edit,
            icon=os.environ["pychemqt"]+"/images/button/fit",
            parent=self)
        self.editMarginAction = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Edit &Margins"),
            slot=self.toolbar.configure_subplots, parent=self)
        self.saveAction = createAction(
            QtWidgets.QApplication.translate("pychemqt", "&Save Plot"),
            slot=self.toolbar.save_figure,
            icon=os.environ["pychemqt"]+"/images/button/fileSave", parent=self)
        self.toolbarVisibleAction = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Toggle &Toolbar"),
            self.toolbar.setVisible, checkable=True, parent=self)
        self.gridToggleAction = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Toggle &Grid"),
            self.grid, checkable=True, parent=self)
        grid = config.Preferences.getboolean("MEOS", "grid")
        self.gridToggleAction.setChecked(grid)

        if dim == 2:
            self.plot.fig.canvas.mpl_connect('button_press_event', self.click)
        else:
            self.editMarginAction.setEnabled(False)

    def contextMenuEvent(self, event):
        """Create context menu"""
        menuTable = QtWidgets.QMenu(
            QtWidgets.QApplication.translate("pychemqt", "Tabulated data"))
        menuTable.setIcon(
            QtGui.QIcon(os.environ["pychemqt"]+"/images/button/table"))
        for linea in self.plot.ax.lines:
            action = createAction(linea.get_label(),
                                  slot=partial(self.table, linea), parent=self)
            menuTable.addAction(action)

        menu = QtWidgets.QMenu()
        menu.addAction(self.editAxesAction)
        menu.addAction(self.editAction)
        menu.addAction(self.editMarginAction)
        menu.addSeparator()
        menu.addAction(self.saveAction)
        menu.addAction(menuTable.menuAction())
        menu.addSeparator()
        menu.addAction(self.toolbarVisibleAction)
        menu.addAction(self.gridToggleAction)
        menu.exec_(event.globalPos())

        if self.plot.ax._gridOn:
            self.gridToggleAction.setChecked(True)

    def grid(self, bool):
        self.plot.ax.grid(bool)
        self.plot.ax._gridOn = bool
        self.plot.draw()

    def edit(self):
        dialog = EditPlot(self, self.parent)
        dialog.exec_()

    def editAxis(self):
        dialog = EditAxis(self.plot, self.parent)
        dialog.exec_()

    def table(self, obj):
        """Export plot data to table
        Input:
            obj: object (Line2D instance) to show data
        """
        xtxt = meos.propiedades[meos.keys.index(self.x)]
        ytxt = meos.propiedades[meos.keys.index(self.y)]
        xunit = meos.units[meos.keys.index(self.x)]
        yunit = meos.units[meos.keys.index(self.y)]
        HHeader = [xtxt+os.linesep+xunit.text(), ytxt+os.linesep+yunit.text()]
        units = [xunit, yunit]
        if self.dim == 3:
            ztxt = meos.propiedades[meos.keys.index(self.z)]
            zunit = meos.units[meos.keys.index(self.z)]
            HHeader.append(ztxt+os.linesep+zunit.text())
            units.append(zunit)
            data = obj._verts3d
        else:
            data = obj.get_data(orig=True)

        from .table import TablaMEoS
        tabla = TablaMEoS(self.dim, horizontalHeader=HHeader, units=units,
                          stretch=False, readOnly=True, parent=self.parent)
        tabla.Point = mEoS.__all__[self.config["fluid"]]
        tabla.setData(transpose(data))
        tabla.verticalHeader().setContextMenuPolicy(
            QtCore.Qt.CustomContextMenu)

        title = QtWidgets.QApplication.translate("pychemqt", "Table from") + \
            " " + obj.get_label()
        tabla.setWindowTitle(title)
        self.parent.centralwidget.currentWidget().addSubWindow(tabla)
        tabla.show()

    def _getData(self):
        """Get data from file"""
        filenameHard = os.environ["pychemqt"]+"dat"+os.sep+"mEoS" + \
            os.sep + self.filename+".gz"
        filenameSoft = config.conf_dir+self.filename
        if os.path.isfile(filenameSoft):
            print(filenameSoft)
            with open(filenameSoft, "rb") as archivo:
                data = pickle.load(archivo, fix_imports=False, errors="strict")
            return data
        elif os.path.isfile(filenameHard):
            print(filenameHard)
            with gzip.GzipFile(filenameHard, 'rb') as archivo:
                data = pickle.load(archivo, encoding="latin1")
            self._saveData(data)
            return data

    def _saveData(self, data):
        """Save changes in data to file"""
        with open(config.conf_dir+self.filename, 'wb') as file:
            pickle.dump(data, file)

    def click(self, event):
        """Update input and graph annotate when mouse click over chart"""
        # Accept only left click
        print(event, self.x, self.y)
        if event.button != 1:
            return
        units = {"x": unidades.Dimensionless,
                 "T": unidades.Temperature,
                 "P": unidades.Pressure,
                 "h": unidades.Enthalpy,
                 "u": unidades.Enthalpy,
                 "s": unidades.SpecificHeat,
                 "v": unidades.SpecificVolume,
                 "rho": unidades.Density}
        if self.x in units and self.y in units:
            x = units[self.x](event.xdata, "conf")
            y = units[self.y](event.ydata, "conf")

            fluid = mEoS.__all__[self.config["fluid"]]
            kwargs = {self.x: x, self.y: y}
            print(fluid, self.config, kwargs)
            fluido = calcPoint(fluid, self.config, **kwargs)
            if fluido and fluido.status and \
                    fluido._constants["Tmin"] <= fluido.T and \
                    fluido.T <= fluido._constants["Tmax"] and \
                    0 < fluido.P.kPa and \
                    fluido.P.kPa <= fluido._constants["Pmax"]:
                self.plot.lx.set_ydata(event.ydata)
                self.plot.ly.set_xdata(event.xdata)
                self.plot.lx.set_visible(True)
                self.plot.ly.set_visible(True)
                self.showPointData(fluido)
            else:
                self.plot.lx.set_visible(False)
                self.plot.ly.set_visible(False)
                self.clearPointData()

    def showPointData(self, state):
        self.clearPointData()
        yi = 0.98
        for key in ("T", "P", "x", "v", "rho", "h", "s", "u"):
            self.notes.append(self.plot.ax.annotate(
                "%s: %s" % (key, state.__getattribute__(key).str), (0.01, yi),
                xycoords='axes fraction', size="small", va="center"))
            yi -= 0.025
        self.plot.draw()

    def clearPointData(self):
        while self.notes:
            anotation = self.notes.pop()
            anotation.remove()
        self.plot.draw()

    def writeToJSON(self, data):
        """Write instance parameter to file"""
        data["filename"] = self.filename
        data["windowTitle"] = self.windowTitle()
        data["x"] = self.x
        data["y"] = self.y
        data["z"] = self.z

        # TODO: Add support for save font properties
        # Title format
        title = {}
        title["txt"] = self.plot.ax.get_title()
        title["color"] = self.plot.ax.title.get_color()
        title["family"] = self.plot.ax.title.get_fontfamily()
        title["style"] = self.plot.ax.title.get_style()
        title["weight"] = self.plot.ax.title.get_weight()
        title["stretch"] = self.plot.ax.title.get_stretch()
        title["size"] = self.plot.ax.title.get_size()
        data["title"] = title

        # xlabel format
        xlabel = {}
        xlabel["txt"] = self.plot.ax.get_xlabel()
        xlabel["color"] = self.plot.ax.xaxis.get_label().get_color()
        xlabel["family"] = self.plot.ax.xaxis.get_label().get_fontfamily()
        xlabel["style"] = self.plot.ax.xaxis.get_label().get_style()
        xlabel["weight"] = self.plot.ax.xaxis.get_label().get_weight()
        xlabel["stretch"] = self.plot.ax.xaxis.get_label().get_stretch()
        xlabel["size"] = self.plot.ax.xaxis.get_label().get_size()
        data["xlabel"] = xlabel

        # ylable format
        ylabel = {}
        ylabel["txt"] = self.plot.ax.get_ylabel()
        ylabel["color"] = self.plot.ax.yaxis.get_label().get_color()
        ylabel["family"] = self.plot.ax.yaxis.get_label().get_fontfamily()
        ylabel["style"] = self.plot.ax.yaxis.get_label().get_style()
        ylabel["weight"] = self.plot.ax.yaxis.get_label().get_weight()
        ylabel["stretch"] = self.plot.ax.yaxis.get_label().get_stretch()
        ylabel["size"] = self.plot.ax.yaxis.get_label().get_size()
        data["ylabel"] = ylabel

        # zlable format
        zlabel = {}
        if self.z:
            zlabel["txt"] = self.plot.ax.get_zlabel()
            zlabel["color"] = self.plot.ax.zaxis.get_label().get_color()
            zlabel["family"] = self.plot.ax.zaxis.get_label().get_fontfamily()
            zlabel["style"] = self.plot.ax.zaxis.get_label().get_style()
            zlabel["weight"] = self.plot.ax.zaxis.get_label().get_weight()
            zlabel["stretch"] = self.plot.ax.zaxis.get_label().get_stretch()
            zlabel["size"] = self.plot.ax.zaxis.get_label().get_size()
        data["zlabel"] = zlabel

        data["grid"] = self.plot.ax._gridOn
        data["xscale"] = self.plot.ax.get_xscale()
        data["yscale"] = self.plot.ax.get_yscale()
        xmin, xmax = self.plot.ax.get_xlim()
        data["xmin"] = xmin
        data["xmax"] = xmax
        ymin, ymax = self.plot.ax.get_ylim()
        data["ymin"] = ymin
        data["ymax"] = ymax
        if self.z:
            zmin, zmax = self.plot.ax.get_zlim()
            data["zmin"] = zmin
            data["zmax"] = zmax
        else:
            data["zmin"] = None
            data["zmax"] = None

        data["marginleft"] = self.plot.fig.subplotpars.left
        data["marginbottom"] = self.plot.fig.subplotpars.bottom
        data["marginright"] = self.plot.fig.subplotpars.right
        data["margintop"] = self.plot.fig.subplotpars.top

        # Config
        data["fluid"] = self.config["fluid"]
        data["eq"] = self.config["eq"]
        data["visco"] = self.config["visco"]
        data["thermal"] = self.config["thermal"]

        # data
        lines = {}
        for line in self.plot.ax.lines[2:]:
            dat = {}
            dat["x"] = list(line.get_xdata())
            dat["y"] = list(line.get_ydata())
            dat["label"] = line.get_label()

            # line style
            dat["lw"] = line.get_lw()
            dat["ls"] = line.get_ls()
            dat["marker"] = line.get_marker()
            dat["color"] = line.get_color()
            dat["ms"] = line.get_ms()
            dat["mfc"] = line.get_mfc()
            dat["mew"] = line.get_mew()
            dat["mec"] = line.get_mec()
            dat["visible"] = line.get_visible()
            dat["antialiased"] = line.get_antialiased()

            # line text
            # saturation and melting line dont define it at plot creation
            try:
                text = {}
                text["visible"] = line.text.get_visible()
                text["txt"] = line.text.get_text()
                text["rot"] = line.text.get_rotation()
                text["pos"] = line.text.pos
                text["family"] = line.text.get_fontfamily()
                text["style"] = line.text.get_style()
                text["weight"] = line.text.get_weight()
                text["stretch"] = line.text.get_stretch()
                text["size"] = line.text.get_size()
                text["va"] = line.text.get_va()
            except AttributeError:
                text = {"visible": False, "txt": "", "pos": 50, "rot": 0,
                        "family": "sans-serif", "style": "normal",
                        "weight": "normal", "stretch": "normal",
                        "size": "small", "va": "center"}
            dat["annotation"] = text

            lines[line._label] = dat
        data["lines"] = lines

    @classmethod
    def readFromJSON(cls, data, parent):
        filename = data["filename"]
        title = data["windowTitle"]
        x = data["x"]
        y = data["y"]
        z = data["z"]
        if z:
            dim = 3
        else:
            dim = 2
        grafico = PlotMEoS(dim=dim, parent=parent, filename=filename)
        grafico.x = x
        grafico.y = y
        grafico.z = z
        grafico.setWindowTitle(title)

        title = data["title"]["txt"]
        if title:
            grafico.plot.ax.set_title(title)
            grafico.plot.ax.title.set_color(data["title"]["color"])
            grafico.plot.ax.title.set_family(data["title"]["family"])
            grafico.plot.ax.title.set_style(data["title"]["style"])
            grafico.plot.ax.title.set_weight(data["title"]["weight"])
            grafico.plot.ax.title.set_stretch(data["title"]["stretch"])
            grafico.plot.ax.title.set_size(data["title"]["size"])

        xlabel = data["xlabel"]["txt"]
        if xlabel:
            grafico.plot.ax.set_xlabel(xlabel)
            label = grafico.plot.ax.xaxis.get_label()
            label.set_color(data["xlabel"]["color"])
            label.set_family(data["xlabel"]["family"])
            label.set_style(data["xlabel"]["style"])
            label.set_weight(data["xlabel"]["weight"])
            label.set_stretch(data["xlabel"]["stretch"])
            label.set_size(data["xlabel"]["size"])

        ylabel = data["ylabel"]["txt"]
        if ylabel:
            grafico.plot.ax.set_ylabel(ylabel)
            label = grafico.plot.ax.yaxis.get_label()
            label.set_color(data["ylabel"]["color"])
            label.set_family(data["ylabel"]["family"])
            label.set_style(data["ylabel"]["style"])
            label.set_weight(data["ylabel"]["weight"])
            label.set_stretch(data["ylabel"]["stretch"])
            label.set_size(data["ylabel"]["size"])

        if z:
            zlabel = data["zlabel"]["txt"]
            if zlabel:
                grafico.plot.ax.set_zlabel(zlabel)
                label = grafico.plot.ax.zaxis.get_label()
                label.set_color(data["zlabel"]["color"])
                label.set_family(data["zlabel"]["family"])
                label.set_style(data["zlabel"]["style"])
                label.set_weight(data["zlabel"]["weight"])
                label.set_stretch(data["zlabel"]["stretch"])
                label.set_size(data["zlabel"]["size"])

        grafico.plot.ax._gridOn = data["grid"]
        grafico.plot.ax.grid(data["grid"])

        grafico.plot.ax.set_xlim(data["xmin"], data["xmax"])
        grafico.plot.ax.set_ylim(data["ymin"], data["ymax"])
        if z:
            grafico.plot.ax.set_zlim(data["zmin"], data["zmax"])

        for label, line in data["lines"].items():
            x = line["x"]
            y = line["y"]

            format = {}
            format["lw"] = line["lw"]
            format["ls"] = line["ls"]
            format["marker"] = line["marker"]
            format["color"] = line["color"]
            format["ms"] = line["ms"]
            format["mfc"] = line["mfc"]
            format["mew"] = line["mew"]
            format["mec"] = line["mec"]

            ln, = grafico.plot.ax.plot(x, y, label=label, **format)
            ln.set_visible(line["visible"])
            ln.set_antialiased(line["antialiased"])

            txt = line["annotation"]["txt"]
            rot = line["annotation"]["rot"]
            pos = line["annotation"]["pos"]
            i = int(len(x)*pos/100)
            kw = {}
            kw["ha"] = "center"
            kw["rotation_mode"] = "anchor"
            for key in ("va", "visible", "family", "style", "weight",
                        "stretch", "size"):
                kw[key] = line["annotation"][key]

            if i >= len(x):
                i = len(x)-1
            text = grafico.plot.ax.text(x[i], y[i], txt, rotation=rot, **kw)

            # We creating a link between line and its annotation text
            ln.text = text
            # We save position value in % unit to avoid index find
            ln.text.pos = pos

        grafico.plot.ax.set_xscale(data["xscale"])
        grafico.plot.ax.set_yscale(data["yscale"])

        # Load margins
        left = data["marginleft"]
        bottom = data["marginbottom"]
        right = data["marginright"]
        top = data["margintop"]
        grafico.plot.fig.subplots_adjust(left, bottom, right, top)

        # Load config
        conf = {}
        conf["fluid"] = data["fluid"]
        conf["eq"] = data["eq"]
        conf["visco"] = data["visco"]
        conf["thermal"] = data["thermal"]
        grafico.config = conf

        return grafico


class Plot2D(QtWidgets.QDialog):
    """Dialog for select a special 2D plot"""
    def __init__(self, parent=None):
        super(Plot2D, self).__init__(parent)
        self.setWindowTitle(
            QtWidgets.QApplication.translate("pychemqt", "Setup 2D Plot"))
        layout = QtWidgets.QVBoxLayout(self)
        group_Ejex = QtWidgets.QGroupBox(
            QtWidgets.QApplication.translate("pychemqt", "Axis X"))
        layout.addWidget(group_Ejex)
        layout_GroupX = QtWidgets.QVBoxLayout(group_Ejex)
        self.ejeX = QtWidgets.QComboBox()
        layout_GroupX.addWidget(self.ejeX)
        self.Xscale = QtWidgets.QCheckBox(
            QtWidgets.QApplication.translate("pychemqt", "Logarithmic scale"))
        layout_GroupX.addWidget(self.Xscale)
        for prop in ThermoAdvanced.propertiesName():
            self.ejeX.addItem(prop)

        group_Ejey = QtWidgets.QGroupBox(
            QtWidgets.QApplication.translate("pychemqt", "Axis Y"))
        layout.addWidget(group_Ejey)
        layout_GroupY = QtWidgets.QVBoxLayout(group_Ejey)
        self.ejeY = QtWidgets.QComboBox()
        layout_GroupY.addWidget(self.ejeY)
        self.Yscale = QtWidgets.QCheckBox(
            QtWidgets.QApplication.translate("pychemqt", "Logarithmic scale"))
        layout_GroupY.addWidget(self.Yscale)

        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox)

        self.ejeXChanged(0)
        self.ejeX.currentIndexChanged.connect(self.ejeXChanged)

    def ejeXChanged(self, index):
        """Fill variables available in ejeY, all except the active in ejeX"""
        # Save current status to restore
        current = self.ejeY.currentIndex()
        if current == -1:
            current = 0

        # Refill ejeY combo
        self.ejeY.clear()
        props = ThermoAdvanced.propertiesName()
        del props[index]
        for prop in props:
            self.ejeY.addItem(prop)

        # Restore inicial state
        if index == 0 and current == 0:
            self.ejeY.setCurrentIndex(0)
        elif index <= current:
            self.ejeY.setCurrentIndex(current)
        else:
            self.ejeY.setCurrentIndex(current+1)


class Plot3D(QtWidgets.QDialog):
    """Dialog for configure a 3D plot"""

    def __init__(self, parent=None):
        super(Plot3D, self).__init__(parent)
        self.setWindowTitle(
            QtWidgets.QApplication.translate("pychemqt", "Setup 3D Plot"))
        layout = QtWidgets.QGridLayout(self)

        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Axis X")), 1, 1)
        self.ejeX = QtWidgets.QComboBox()
        for prop in ThermoAdvanced.propertiesName():
            self.ejeX.addItem(prop)
        layout.addWidget(self.ejeX, 1, 2)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Axis Y")), 2, 1)
        self.ejeY = QtWidgets.QComboBox()
        layout.addWidget(self.ejeY, 2, 2)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Axis Z")), 3, 1)
        self.ejeZ = QtWidgets.QComboBox()
        layout.addWidget(self.ejeZ, 3, 2)

        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox, 4, 1, 1, 2)

        self.ejeX.currentIndexChanged.connect(self.ejeXChanged)
        self.ejeY.currentIndexChanged.connect(self.ejeYChanged)
        self.ejeXChanged(0)

    def ejeXChanged(self, index):
        """Fill variables available in ejeY, all except the active in ejeX"""
        # Save current status to restore
        current = self.ejeY.currentIndex()
        if current == -1:
            current = 0

        # Refill ejeY combo
        self.ejeY.clear()
        props = ThermoAdvanced.propertiesName()
        del props[index]
        for prop in props:
            self.ejeY.addItem(prop)

        # Restore inicial state
        if index == 0 and current == 0:
            self.ejeY.setCurrentIndex(0)
        elif index <= current:
            self.ejeY.setCurrentIndex(current)
        else:
            self.ejeY.setCurrentIndex(current+1)

    def ejeYChanged(self, indY):
        """Fill variables available in ejeZ, all except the actives in other"""
        # Save current status to restore
        current = self.ejeZ.currentIndex()
        if current == -1:
            current = 0

        # Refill ejeY combo
        self.ejeZ.clear()
        prop2 = ThermoAdvanced.propertiesName()[:]
        indX = self.ejeX.currentIndex()
        del prop2[indX]
        del prop2[indY]
        for prop in prop2:
            self.ejeZ.addItem(prop)

        # Restore inicial state
        if indX == 0 and indY == 0 and current == 0:
            self.ejeZ.setCurrentIndex(0)
        elif indY <= current or indX <= current:
            self.ejeZ.setCurrentIndex(current)
        else:
            self.ejeZ.setCurrentIndex(current+1)


class EditPlot(QtWidgets.QDialog):
    """Dialog to edit plot. This dialog let user change plot p"""
    def __init__(self, plotMEoS, parent=None):
        super(EditPlot, self).__init__(parent)
        self.setWindowTitle(
            QtWidgets.QApplication.translate("pychemqt", "Edit Plot"))
        layout = QtWidgets.QGridLayout(self)
        self.plotMEoS = plotMEoS
        self.fig = plotMEoS.plot
        self.parent = parent
        self.semaforo = QtCore.QSemaphore(1)

        self.lista = QtWidgets.QListWidget()
        layout.addWidget(self.lista, 0, 1, 1, 3)

        lytTitle = QtWidgets.QHBoxLayout()
        label = QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Label"))
        lytTitle.addWidget(label)
        self.label = QtWidgets.QLineEdit()
        lytTitle.addWidget(self.label)
        layout.addLayout(lytTitle, 1, 1, 1, 3)

        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Line Width")), 2, 1)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Line Style")), 2, 2)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Color")), 2, 3)
        self.Grosor = QtWidgets.QDoubleSpinBox()
        self.Grosor.setAlignment(QtCore.Qt.AlignRight)
        self.Grosor.setRange(0.1, 5)
        self.Grosor.setDecimals(1)
        self.Grosor.setSingleStep(0.1)
        layout.addWidget(self.Grosor, 3, 1)
        self.Linea = LineStyleCombo()
        layout.addWidget(self.Linea, 3, 2)
        self.ColorButton = ColorSelector()
        layout.addWidget(self.ColorButton, 3, 3)

        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Marker")), 4, 1)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Marker Size")), 4, 2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Marker Color")), 4, 3)
        self.Marca = MarkerCombo()
        layout.addWidget(self.Marca, 5, 1)
        self.markerSize = QtWidgets.QDoubleSpinBox()
        self.markerSize.setAlignment(QtCore.Qt.AlignRight)
        self.markerSize.setDecimals(1)
        self.markerSize.setSingleStep(0.1)
        layout.addWidget(self.markerSize, 5, 2)
        self.markerfacecolor = ColorSelector()
        layout.addWidget(self.markerfacecolor, 5, 3)

        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Marker edge")), 7, 1)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Width")), 6, 2)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Color")), 6, 3)
        self.markerEdgeSize = QtWidgets.QDoubleSpinBox()
        self.markerEdgeSize.setAlignment(QtCore.Qt.AlignRight)
        self.markerEdgeSize.setDecimals(1)
        self.markerEdgeSize.setSingleStep(0.1)
        layout.addWidget(self.markerEdgeSize, 7, 2)
        self.markeredgecolor = ColorSelector()
        layout.addWidget(self.markeredgecolor, 7, 3)

        grpAnnotate = QtWidgets.QGroupBox(
            QtWidgets.QApplication.translate("pychemqt", "Annotation"))
        layout.addWidget(grpAnnotate, 8, 1, 1, 3)
        lytAnnotation = QtWidgets.QGridLayout(grpAnnotate)
        self.annotationVisible = QtWidgets.QCheckBox(
            QtWidgets.QApplication.translate("pychemqt", "Visible"))
        lytAnnotation.addWidget(self.annotationVisible, 1, 1, 1, 3)

        lytTitle = QtWidgets.QHBoxLayout()
        label = QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Label"))
        lytTitle.addWidget(label)
        # self.annotationLabel = QtWidgets.QLineEdit()
        self.annotationLabel = InputFont()
        lytTitle.addWidget(self.annotationLabel)
        lytAnnotation.addLayout(lytTitle, 2, 1, 1, 3)

        lytPosition = QtWidgets.QHBoxLayout()
        lytPosition.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Location")))
        self.labelAnnotationPos = Entrada_con_unidades(
            int, value=50, width=40, frame=False, readOnly=True, suffix="%",
            showNull=True)
        self.labelAnnotationPos.setFixedWidth(40)
        lytPosition.addWidget(self.labelAnnotationPos)
        self.annotationPos = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.annotationPos.setRange(0, 100)
        self.annotationPos.setValue(50)
        self.annotationPos.valueChanged.connect(
            partial(self._updateLabel, self.labelAnnotationPos))
        lytPosition.addWidget(self.annotationPos)
        lytAnnotation.addLayout(lytPosition, 3, 1, 1, 3)

        lytAngle = QtWidgets.QHBoxLayout()
        lytAngle.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Rotation")))
        self.labelAnnotationRot = Entrada_con_unidades(
            int, value=50, width=40, frame=False, readOnly=True, suffix="º",
            showNull=True)
        self.labelAnnotationRot.setFixedWidth(40)
        lytAngle.addWidget(self.labelAnnotationRot)
        self.annotationRot = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.annotationRot.setRange(0, 360)
        self.annotationRot.setValue(0)
        self.annotationRot.valueChanged.connect(
            partial(self._updateLabel, self.labelAnnotationRot))
        lytAngle.addWidget(self.annotationRot)
        lytAnnotation.addLayout(lytAngle, 4, 1, 1, 3)

        lytVA = QtWidgets.QHBoxLayout()
        lytVA.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Aligment")))
        self.annotationVA = QtWidgets.QComboBox()
        alignment = [
            QtWidgets.QApplication.translate("pychemqt", "Center"),
            QtWidgets.QApplication.translate("pychemqt", "Top"),
            QtWidgets.QApplication.translate("pychemqt", "Bottom"),
            QtWidgets.QApplication.translate("pychemqt", "Baseline"),
            QtWidgets.QApplication.translate("pychemqt", "Center baseline")]
        for alig in alignment:
            self.annotationVA.addItem(alig)
        lytVA.addWidget(self.annotationVA)
        lytVA.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding))
        lytAnnotation.addLayout(lytVA, 5, 1, 1, 3)

        self.annotationVisible.stateChanged.connect(
            self.annotationLabel.setEnabled)
        self.annotationVisible.stateChanged.connect(
            self.annotationPos.setEnabled)
        self.annotationVisible.stateChanged.connect(
            self.annotationRot.setEnabled)

        self.visible = QtWidgets.QCheckBox(
            QtWidgets.QApplication.translate("pychemqt", "Visible"))
        layout.addWidget(self.visible, 13, 1, 1, 3)
        self.antialiases = QtWidgets.QCheckBox(
            QtWidgets.QApplication.translate("pychemqt", "Antialiases"))
        layout.addWidget(self.antialiases, 14, 1, 1, 3)

        layoutButton = QtWidgets.QHBoxLayout()
        layout.addLayout(layoutButton, 15, 1, 1, 3)
        self.botonAdd = QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"] + "/images/button/add.png")), "")
        self.botonAdd.clicked.connect(self.add)
        layoutButton.addWidget(self.botonAdd)
        self.botonRemove = QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"] + "/images/button/remove.png")), "")
        self.botonRemove.clicked.connect(self.remove)
        layoutButton.addWidget(self.botonRemove)
        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Close)
        self.buttonBox.rejected.connect(self.close)
        layoutButton.addWidget(self.buttonBox)

        for linea in self.fig.ax.lines[2:]:
            self.lista.addItem(linea._label)

        self.lista.currentRowChanged.connect(self.update)
        self.lista.setCurrentRow(0)
        self.label.textChanged.connect(partial(self.changeValue, "label"))
        self.Grosor.valueChanged.connect(partial(self.changeValue, "lw"))
        self.Linea.valueChanged.connect(partial(self.changeValue, "ls"))
        self.Linea.currentIndexChanged.connect(self.ColorButton.setEnabled)
        self.ColorButton.valueChanged.connect(
            partial(self.changeValue, "color"))
        self.Marca.valueChanged.connect(partial(self.changeValue, "marker"))
        self.Marca.currentIndexChanged.connect(self.markerSize.setEnabled)
        self.Marca.currentIndexChanged.connect(self.markerfacecolor.setEnabled)
        self.Marca.currentIndexChanged.connect(self.markerEdgeSize.setEnabled)
        self.Marca.currentIndexChanged.connect(self.markeredgecolor.setEnabled)
        self.markerSize.valueChanged.connect(partial(self.changeValue, "ms"))
        self.markerfacecolor.valueChanged.connect(
            partial(self.changeValue, "mfc"))
        self.markerEdgeSize.valueChanged.connect(
            partial(self.changeValue, "mew"))
        self.markeredgecolor.valueChanged.connect(
            partial(self.changeValue, "mec"))
        self.visible.toggled.connect(partial(self.changeValue, "visible"))
        self.antialiases.toggled.connect(
            partial(self.changeValue, "antialiases"))

        self.annotationVisible.toggled.connect(
            partial(self.changeValue, "textVisible"))
        self.annotationLabel.textChanged.connect(
            partial(self.changeValue, "textLabel"))
        self.annotationLabel.colorChanged.connect(
                partial(self.changeValue, "textcolor"))
        self.annotationLabel.fontChanged.connect(
            partial(self.changeValue, "textfont"))
        self.annotationPos.valueChanged.connect(
            partial(self.changeValue, "textPos"))
        self.annotationRot.valueChanged.connect(
            partial(self.changeValue, "textRot"))
        self.annotationVA.currentIndexChanged.connect(
            partial(self.changeValue, "textVA"))

    def _updateLabel(self, label, value):
        label.setValue(value)

    def update(self, i):
        """Fill format widget with value of selected line"""
        if self.semaforo.available() > 0:
            self.semaforo.acquire(1)
            linea = self.fig.ax.lines[i+2]
            self.label.setText(linea.get_label())
            self.Grosor.setValue(linea.get_lw())
            self.Linea.setCurrentValue(linea.get_ls())
            self.ColorButton.setColor(linea.get_color())
            self.Marca.setCurrentValue(linea.get_marker())
            self.markerSize.setValue(linea.get_ms())
            self.markerfacecolor.setColor(linea.get_mfc())
            self.markerEdgeSize.setValue(linea.get_mew())
            self.markeredgecolor.setColor(linea.get_mec())
            self.visible.setChecked(linea.get_visible())
            self.antialiases.setChecked(linea.get_antialiased())

            try:
                self.annotationVisible.setChecked(linea.text.get_visible())
                self.annotationLabel.setText(linea.text.get_text())
                self.annotationPos.setValue(linea.text.pos)
                self.annotationRot.setValue(linea.text.get_rotation())
                va = ["center", "top", "bottom", "baseline", "center_baseline"]
                self.annotationVA.setCurrentIndex(va.index(linea.text.get_va()))
            except AttributeError:
                self.annotationVisible.setChecked(False)
            self.semaforo.release(1)

    def changeValue(self, key, value):
        """Update plot data"""
        if self.semaforo.available() > 0:
            self.semaforo.acquire(1)
            linea = self.fig.ax.lines[self.lista.currentRow()+2]
            func = {"label": linea.set_label,
                    "lw": linea.set_lw,
                    "ls": linea.set_ls,
                    "marker": linea.set_marker,
                    "color": linea.set_color,
                    "ms": linea.set_ms,
                    "mfc": linea.set_mfc,
                    "mew": linea.set_mew,
                    "mec": linea.set_mec,
                    "visible": linea.set_visible,
                    "antialiases": linea.set_antialiased,
                    "textVisible": linea.text.set_visible,
                    "textLabel": linea.text.set_text,
                    "textcolor": linea.text.set_color,
                    "textfont": linea.text.set_fontproperties,
                    "textPos": linea.text.set_position,
                    "textRot": linea.text.set_rotation,
                    "textVA": linea.text.set_va}

            if key == "textPos":
                linea.text.pos = value
                xi = linea.get_xdata()
                yi = linea.get_ydata()
                i = int(len(xi)*value/100)
                if i >= len(xi):
                    i = len(yi)-1
                value = xi[i], yi[i]
            elif key == "textVA":
                va = ["center", "top", "bottom", "baseline", "center_baseline"]
                value = va[value]
            elif key == "textfont":
                value = convertFont(value)
            elif key in ("ls", "marker", "color", "mfc", "mec"):
                value = str(value)
            func[key](value)
            if key == "label":
                self.lista.currentItem().setText(value)
            else:
                self.fig.draw()
                self.parent.dirty[self.parent.idTab] = True
                self.parent.saveControl()
            self.semaforo.release(1)

    def add(self):
        """Add a isoline to plot"""
        dialog = AddLine()
        if dialog.exec_():
            points = get_points(config.Preferences)
            self.parent.progressBar.setVisible(True)
            index = self.parent.currentConfig.getint("MEoS", "fluid")
            # fluid = getClassFluid(self.config)
            fluid = mEoS.__all__[index]
            prop = dialog.tipo.currentIndex()
            value = dialog.input[prop].value

            eq = fluid.eq[self.parent.currentConfig.getint("MEoS", "eq")]
            T = list(concatenate([
                linspace(eq["Tmin"], 0.9*fluid.Tc, points),
                linspace(0.9*fluid.Tc, 0.99*fluid.Tc, points),
                linspace(0.99*fluid.Tc, fluid.Tc, points),
                linspace(fluid.Tc, 1.01*fluid.Tc, points),
                linspace(1.01*fluid.Tc, 1.1*fluid.Tc, points),
                linspace(1.1*fluid.Tc, eq["Tmax"], points)]))

            Pmin = fluid(T=eq["Tmin"], x=0).P
            Pmax = eq["Pmax"]*1000
            P = list(concatenate([
                logspace(log10(Pmin), log10(0.9*fluid.Pc), points),
                linspace(0.9*fluid.Pc, 0.99*fluid.Pc, points),
                linspace(0.99*fluid.Pc, fluid.Pc, points),
                linspace(fluid.Pc, 1.01*fluid.Pc, points),
                linspace(1.01*fluid.Pc, 1.1*fluid.Pc, points),
                logspace(log10(1.1*fluid.Pc), log10(Pmax), points)]))
            for i in range(5, 0, -1):
                del T[points*i]
                del P[points*i]

            if prop == 0:
                # Calcualte isotherm line
                self.parent.statusbar.showMessage(
                    QtWidgets.QApplication.translate(
                        "pychemqt", "Adding isotherm line..."))
                fluidos = calcIsoline(
                    fluid, self.parent.currentConfig, "P", "T", P, value,
                    0, 0, 100, 1, self.parent.progressBar)
                var = "T"
                name = "Isotherm"
                unit = unidades.Temperature
            elif prop == 1:
                # Calculate isobar line
                self.parent.statusbar.showMessage(
                    QtWidgets.QApplication.translate(
                        "pychemqt", "Adding isobar line..."))
                fluidos = calcIsoline(
                    fluid, self.parent.currentConfig, "T", "P", T, value,
                    0, 0, 100, 1, self.parent.progressBar)
                var = "P"
                name = "Isobar"
                unit = unidades.Pressure
            elif prop == 2:
                # Calculate isoenthalpic line
                self.parent.statusbar.showMessage(
                    QtWidgets.QApplication.translate(
                        "pychemqt", "Adding isoenthalpic line..."))
                fluidos = calcIsoline(
                    fluid, self.parent.currentConfig, "P", "h", P, value,
                    0, 0, 100, 1, self.parent.progressBar)
                var = "h"
                name = "Isoenthalpic"
                unit = unidades.Enthalpy
            elif prop == 3:
                # Calculate isoentropic line
                self.parent.statusbar.showMessage(
                    QtWidgets.QApplication.translate(
                        "pychemqt", "Adding isoentropic line..."))
                fluidos = calcIsoline(
                    fluid, self.parent.currentConfig, "T", "s", T, value,
                    0, 0, 100, 1, self.parent.progressBar)
                var = "s"
                name = "Isoentropic"
                unit = unidades.SpecificHeat
            elif prop == 4:
                # Calculate isochor line
                self.parent.statusbar.showMessage(
                    QtWidgets.QApplication.translate(
                        "pychemqt", "Adding isochor line..."))
                fluidos = calcIsoline(
                    fluid, self.parent.currentConfig, "T", "v", T, value,
                    0, 0, 100, 1, self.parent.progressBar)
                var = "v"
                name = "Isochor"
                unit = unidades.SpecificVolume
            elif prop == 5:
                # Calculate isodensity line
                self.parent.statusbar.showMessage(
                    QtWidgets.QApplication.translate(
                        "pychemqt", "Adding isodensity line..."))
                fluidos = calcIsoline(
                    fluid, self.parent.currentConfig, "T", "rho", T, value,
                    0, 0, 100, 1, self.parent.progressBar)
                var = "rho"
                name = "Isochor"
                unit = unidades.Density
            elif prop == 6:
                # Calculate isoquality line
                self.parent.statusbar.showMessage(
                    QtWidgets.QApplication.translate(
                        "pychemqt", "Adding isoquality line..."))
                T = T[:3*points-2]
                fluidos = calcIsoline(
                    fluid, self.parent.currentConfig, "T", "x", T, value,
                    0, 0, 100, 1, self.parent.progressBar)
                var = "x"
                name = "Isoquality"
                unit = unidades.Dimensionless

            line = {value: {}}
            for x in ThermoAdvanced.propertiesKey():
                dat_propiedad = []
                for fluido in fluidos:
                    num = fluido.__getattribute__(x)
                    if isinstance(num, str):
                        dat_propiedad.append(num)
                    elif x in ("f", "fi"):
                        dat_propiedad.append(num[0])
                    elif num is not None:
                        dat_propiedad.append(num._data)
                    else:
                        dat_propiedad.append(None)
                line[value][x] = dat_propiedad

            style = getLineFormat(config.Preferences, name)
            functionx = _getunitTransform(self.plotMEoS.x)
            functiony = _getunitTransform(self.plotMEoS.y)
            functionz = _getunitTransform(self.plotMEoS.z)
            transform = (functionx, functiony, functionz)
            ax = self.plotMEoS.x, self.plotMEoS.y, self.plotMEoS.z
            plotIsoline(line, ax, var, unit, self.plotMEoS, transform, **style)

            self.plotMEoS.plot.draw()
            self.parent.progressBar.setVisible(False)
            self.parent.dirty[self.parent.idTab] = True
            self.parent.saveControl()
            self.lista.addItem(self.fig.ax.lines[-1].get_label())
            self.lista.setCurrentRow(self.lista.count()-1)

            # Save new line to file
            data = self.plotMEoS._getData()
            if var not in data:
                data[var] = {}
            data[var][value] = line[value]
            self.plotMEoS._saveData(data)

    def remove(self):
        """Remove a line from plot"""
        self.parent.statusbar.showMessage(QtWidgets.QApplication.translate(
            "pychemqt", "Deleting line..."))
        QtWidgets.QApplication.processEvents()

        # Remove data from file
        data = self.plotMEoS._getData()
        txt = self.lista.currentItem().text().split()
        var = txt[0]
        units = {"T": unidades.Temperature,
                 "P": unidades.Pressure,
                 "v": unidades.SpecificVolume,
                 "rho": unidades.Density,
                 "h": unidades.Enthalpy,
                 "s": unidades.SpecificHeat,
                 "x": unidades.Dimensionless}
        if var in units:
            unit = units[var]
            for key in data[var]:
                str = unit(key).str
                if str[1:] == " ".join(txt[2:]):
                    del data[var][key]
                    self.plotMEoS._saveData(data)
                    break

        # Remove line to plot and update list element
        index = self.lista.currentRow()
        del self.fig.ax.lines[index+2]
        if index == 0:
            self.lista.setCurrentRow(1)
        else:
            self.lista.setCurrentRow(index-1)
        self.lista.takeItem(index)
        self.fig.draw()
        self.parent.statusbar.clearMessage()
        self.parent.dirty[self.parent.idTab] = True
        self.parent.saveControl()


class AddLine(QtWidgets.QDialog):
    """Dialog to add new isoline to plot"""
    lineas = [(QtWidgets.QApplication.translate("pychemqt", "Isotherm"),
               unidades.Temperature, None),
              (QtWidgets.QApplication.translate("pychemqt", "Isobar"),
               unidades.Pressure, None),
              (QtWidgets.QApplication.translate("pychemqt", "Isoenthalpic"),
               unidades.Enthalpy, None),
              (QtWidgets.QApplication.translate("pychemqt", "Isoentropic"),
              unidades.SpecificHeat, "SpecificEntropy"),
              (QtWidgets.QApplication.translate("pychemqt", "Isochor"),
              unidades.SpecificVolume, None),
              (QtWidgets.QApplication.translate("pychemqt", "Isodensity"),
              unidades.Density, None),
              (QtWidgets.QApplication.translate("pychemqt", "Isoquality"),
              float, None)]

    def __init__(self, parent=None):
        super(AddLine, self).__init__(parent)
        self.setWindowTitle(
            QtWidgets.QApplication.translate("pychemqt", "Add Line to Plot"))
        layout = QtWidgets.QGridLayout(self)

        self.tipo = QtWidgets.QComboBox()
        layout.addWidget(self.tipo, 1, 1, 1, 2)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Value")), 2, 1)

        self.input = []
        for title, unidad, magnitud in self.lineas:
            self.input.append(Entrada_con_unidades(unidad, magnitud))
            layout.addWidget(self.input[-1], 2, 2)
            self.tipo.addItem(title)

        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox, 10, 1, 1, 2)

        self.isolineaChanged(0)
        self.tipo.currentIndexChanged.connect(self.isolineaChanged)

    def isolineaChanged(self, int):
        """Let show only the active inputs"""
        for i in self.input:
            i.setVisible(False)
        self.input[int].setVisible(True)


class EditAxis(QtWidgets.QDialog):
    """Dialog to configure axes plot properties, label, margins, scales"""
    def __init__(self, fig=None, parent=None):
        super(EditAxis, self).__init__(parent)
        self.setWindowTitle(
            QtWidgets.QApplication.translate("pychemqt", "Edit Axis"))
        layout = QtWidgets.QGridLayout(self)
        self.fig = fig

        lytTitle = QtWidgets.QHBoxLayout()
        lb = QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Title"))
        lb.setSizePolicy(
            QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Maximum)
        lytTitle.addWidget(lb)
        self.title = InputFont()
        lytTitle.addWidget(self.title)
        layout.addLayout(lytTitle, 1, 1, 1, self.fig.dim)

        self.axisX = AxisWidget("x", self)
        layout.addWidget(self.axisX, 2, 1)
        self.axisY = AxisWidget("y", self)
        layout.addWidget(self.axisY, 2, 2)

        if self.fig.dim == 3:
            self.axisZ = AxisWidget("z", self)
            layout.addWidget(self.axisZ, 2, 3)
            self.axisX.scale.setEnabled(False)
            self.axisY.scale.setEnabled(False)
            self.axisZ.scale.setEnabled(False)

        self.gridCheckbox = QtWidgets.QCheckBox(
            QtWidgets.QApplication.translate("pychemqt", "Show Grid"))
        layout.addWidget(self.gridCheckbox, 3, 1, 1, self.fig.dim)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 5, 1, 1, self.fig.dim)
        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Close)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox, 10, 1, 1, self.fig.dim)

        if fig:
            self.populate()

        self.title.textChanged.connect(partial(self.update, "title"))
        self.title.colorChanged.connect(partial(self.update, "titlecolor"))
        self.title.fontChanged.connect(partial(self.update, "titlefont"))
        self.axisX.label.textChanged.connect(partial(self.update, "xlabel"))
        self.axisX.label.colorChanged.connect(
            partial(self.update, "xlabelcolor"))
        self.axisX.label.fontChanged.connect(
            partial(self.update, "xlabelfont"))
        self.axisY.label.textChanged.connect(partial(self.update, "ylabel"))
        self.axisY.label.colorChanged.connect(
            partial(self.update, "ylabelcolor"))
        self.axisY.label.fontChanged.connect(
            partial(self.update, "ylabelfont"))
        self.gridCheckbox.toggled.connect(partial(self.update, "grid"))
        self.axisX.scale.toggled.connect(partial(self.update, "xscale"))
        self.axisY.scale.toggled.connect(partial(self.update, "yscale"))
        self.axisX.min.valueChanged.connect(partial(self.update, "xmin"))
        self.axisY.min.valueChanged.connect(partial(self.update, "ymin"))
        self.axisX.max.valueChanged.connect(partial(self.update, "xmax"))
        self.axisY.max.valueChanged.connect(partial(self.update, "ymax"))
        if self.fig.dim == 3:
            self.axisZ.label.textChanged.connect(
                partial(self.update, "zlabel"))
            self.axisZ.label.colorChanged.connect(
                partial(self.update, "zlabelcolor"))
            self.axisZ.label.fontChanged.connect(
                partial(self.update, "zlabelfont"))
            self.axisZ.min.valueChanged.connect(partial(self.update, "zmin"))
            self.axisZ.max.valueChanged.connect(partial(self.update, "zmax"))

    def populate(self):
        """Fill widget with plot parameters"""
        self.title.setText(self.fig.ax.get_title())
        self.title.setColor(QtGui.QColor(self.fig.ax.title.get_color()))
        self.axisX.label.setText(self.fig.ax.get_xlabel())
        xcolor = self.fig.ax.xaxis.get_label().get_color()
        self.axisX.label.setColor(QtGui.QColor(xcolor))
        self.axisY.label.setText(self.fig.ax.get_ylabel())
        ycolor = self.fig.ax.yaxis.get_label().get_color()
        self.axisY.label.setColor(QtGui.QColor(ycolor))
        self.gridCheckbox.setChecked(self.fig.ax._gridOn)
        self.axisX.scale.setChecked(self.fig.ax.get_xscale() == "log")
        self.axisY.scale.setChecked(self.fig.ax.get_yscale() == "log")
        xmin, xmax = self.fig.ax.get_xlim()
        self.axisX.min.setValue(xmin)
        self.axisX.max.setValue(xmax)
        ymin, ymax = self.fig.ax.get_ylim()
        self.axisY.min.setValue(ymin)
        self.axisY.max.setValue(ymax)
        if self.fig.dim == 3:
            self.axisZ.label.setText(self.fig.ax.get_zlabel())
            zcolor = self.fig.ax.zaxis.get_label().get_color()
            self.axisZ.label.setColor(QtGui.QColor(zcolor))
            zmin, zmax = self.fig.ax.get_zlim()
            self.axisZ.min.setValue(zmin)
            self.axisZ.max.setValue(zmax)

    def update(self, key, value):
        """Update plot
        Input:
            key: plot parameter key to update
            value: new value for key
        """
        f = {"xlabel": self.fig.ax.set_xlabel,
             "xlabelcolor": self.fig.ax.xaxis.get_label().set_color,
             "xlabelfont": self.fig.ax.xaxis.get_label().set_fontproperties,
             "ylabel": self.fig.ax.set_ylabel,
             "ylabelcolor": self.fig.ax.yaxis.get_label().set_color,
             "ylabelfont": self.fig.ax.yaxis.get_label().set_fontproperties,
             "title": self.fig.ax.set_title,
             "titlecolor": self.fig.ax.title.set_color,
             "titlefont": self.fig.ax.title.set_fontproperties,
             "xscale": self.fig.ax.set_xscale,
             "yscale": self.fig.ax.set_yscale,
             "grid": self.fig.ax.grid}

        if self.fig.dim == 3:
            f["zlabel"] = self.fig.ax.set_zlabel
            f["zlabelcolor"] = self.fig.ax.zaxis.get_label().set_color
            f["zlabelfont"] = self.fig.ax.zaxis.get_label().set_fontproperties

        if key in ("xscale", "yscale"):
            if value:
                value = "log"
            else:
                value = "linear"
        if key == "grid":
            self.fig.ax._gridOn = value
        if key in ("titlecolor", "xlabelcolor", "ylabelcolor"):
            value = str(value)
        if key in ("titlefont", "xlabelfont", "ylabelfont"):
            value = convertFont(value)

        if key in ("xmin", "xmax"):
            xmin = self.axisX.min.value
            xmax = self.axisX.max.value
            self.fig.ax.set_xlim(xmin, xmax)
        elif key in ("ymin", "ymax"):
            ymin = self.axisY.min.value
            ymax = self.axisY.max.value
            self.fig.ax.set_ylim(ymin, ymax)
        elif key in ("zmin", "zmax"):
            ymin = self.axisZ.min.value
            ymax = self.axisZ.max.value
            self.fig.ax.set_zlim(ymin, ymax)
        else:
            f[key](value)
        self.parent().dirty[self.parent().idTab] = True
        self.parent().saveControl()
        self.fig.draw()


def convertFont(qfont):
    """Convert qt QFont class properties to FontProperties to use in
    matplotlib

    Parameters
    ----------
    qfont : QFont
        QFont with properties to extract

    Returns
    -------
    font : FontProperties
        FontProperties instance to use in any matplotlib text instance
    """
    family = str(qfont.family())

    # Matplotlib use 0-1000 scale, qt only 0-100 scale
    weight = 10*qfont.weight()

    if qfont.style() == 0:
        style = "normal"
    elif qfont.style() == 1:
        style = "italic"
    elif qfont.style() == 2:
        style = "oblique"
    else:
        style = None
    # print(family, style, qfont.stretch(), weight, qfont.pointSize())
    font = FontProperties(family, style, None, qfont.stretch(),
                          weight, qfont.pointSize())

    return font


class AxisWidget(QtWidgets.QGroupBox):
    """Dialog to configure axes plot properties"""
    def __init__(self, name, parent=None):
        title = name+" "+QtWidgets.QApplication.translate("pychemqt", "Axis")
        super(AxisWidget, self).__init__(title, parent)
        lyt = QtWidgets.QGridLayout(self)
        lyt.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Label")), 1, 1)
        self.label = InputFont()
        lyt.addWidget(self.label, 1, 2)
        self.scale = QtWidgets.QCheckBox(
            QtWidgets.QApplication.translate("pychemqt", "Logarithmic scale"))
        lyt.addWidget(self.scale, 2, 1, 1, 2)
        lyt.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "from")), 3, 1)
        self.min = Entrada_con_unidades(float, min=float("-inf"))
        lyt.addWidget(self.min, 3, 2)
        lyt.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "to")), 4, 1)
        self.max = Entrada_con_unidades(float, min=float("-inf"))
        lyt.addWidget(self.max, 4, 2)


def calcIsoline(f, conf, var, fix, vvar, vfix, ini, step, end, total, bar):
    """Procedure to calculate isoline. In isotherm and isobar add to calculate
    point the saturated states in two-phases region"""
    fluidos = []
    fail = 0
    N_points = get_points(config.Preferences)
    fase = None
    rhoo = 0
    To = 0
    for Ti in vvar:
        kwargs = {var: Ti, fix: vfix, "rho0": rhoo, "T0": To}
        print(kwargs)
        fluido = calcPoint(f, conf, **kwargs)

        avance = ini + end*step/total + \
                end/total*(len(fluidos)+fail)/(len(vvar)+N_points)
        bar.setValue(avance)
        QtWidgets.QApplication.processEvents()

        if fluido and fluido.status and (fluido.rho != rhoo or fluido.T != To):
            fluidos.append(fluido)


            # Save values of last point as initial guess for next calculation
            if var not in ("T", "P") or fix not in ("T", "P"):
                rhoo = fluido.rho
                To = fluido.T

            if var in ("T", "P") and fix in ("T", "P"):
                if fase is None:
                    fase = fluido.x

                if fase == fluido.x:
                    continue

                print("Calculating two phase additional point")
                if fluido.P < f.Pc and fluido.T < f.Tc:
                    if fase != fluido.x and fase <= 0:
                        xi = linspace(0, 1, N_points)
                    elif fase != fluido.x and fase >= 1:
                        xi = linspace(1, 0, N_points)

                    for x in xi:
                        print({fix: vfix, "x": x})
                        fluido_x = calcPoint(f, conf, **{fix: vfix, "x": x})
                        fluidos.insert(-1, fluido_x)
                        avance = ini + end*step/total + end/total * \
                            (len(fluidos)+fail)/(len(vvar)+N_points)
                        bar.setValue(avance)

                fase = fluido.x
        else:
            fail += 1

    return fluidos


def get_points(Preferences):
    """Get point number to plot lines from Preferences"""
    definition = Preferences.getint("MEOS", "definition")
    if definition == 1:
        points = 10
    elif definition == 2:
        points = 25
    elif definition == 3:
        points = 50
    elif definition == 4:
        points = 100
    else:
        points = 5
    return points


def getLineFormat(Preferences, name):
    """get matplotlib line format from preferences
        Preferences: configparser instance with pycheqmt preferences
        name: name of isoline"""
    format = formatLine(Preferences, "MEOS", name)

    # Anotation
    if name != "saturation":
        format["annotate"] = Preferences.getboolean("MEOS", name+"label")
        format["pos"] = Preferences.getint("MEOS", name+"position")
        format["unit"] = Preferences.getboolean("MEOS", name+"units")
        format["variable"] = Preferences.getboolean("MEOS", name+"variable")

    return format


def plotIsoline(data, axis, title, unidad, grafico, transform, **format):
    """Procedure to plot any isoline
    Input:
        data: section of property isoline of matrix data
        axis: array with keys of three axis, z None in 2D plot
        title: key of isoline type
        unidad: unidades subclass with isoline unit
        grafico: PlotMEoS instance to plot data
        transform: unit transform function for use configurated units in plots
        format: any matplotlib plot kwargs
    """
    x, y, z = axis
    fx, fy, fz = transform
    xscale = grafico.plot.ax.get_xscale()
    yscale = grafico.plot.ax.get_yscale()
    annotate = format.pop("annotate")
    pos = format.pop("pos")
    unit = format.pop("unit")
    variable = format.pop("variable")
    for key in sorted(data.keys()):
        xi = list(map(fx, data[key][x]))
        yi = list(map(fy, data[key][y]))
        label = "%s =%s" % (title, unidad(key).str)
        if z:
            zi = list(map(fz, data[key][z]))
            line, = grafico.plot.ax.plot(xi, yi, zi, label=label, **format)
        else:
            line, = grafico.plot.ax.plot(xi, yi, label=label, **format)

        # Add annotate for isolines
        # if annotate and not z:
        if variable and unit:
            txt = label
        elif variable:
            txt = "%s =%s" % (title, unidad(key).config())
        elif unit:
            txt = unidad(key).str
        else:
            txt = unidad(key).config()

        xmin, xmax = grafico.plot.ax.get_xlim()
        ymin, ymax = grafico.plot.ax.get_ylim()

        i = int(len(xi)*pos/100)
        if i >= len(xi):
            i = len(yi)-2

        if pos > 50:
            j = i-1
        else:
            j = i+1
        if xscale == "log":
            f_x = (log(xi[i])-log(xi[j]))/(log(xmax)-log(xmin))
        else:
            f_x = (xi[i]-xi[j])/(xmax-xmin)
        if yscale == "log":
            f_y = (log(yi[i])-log(yi[j]))/(log(ymax)-log(ymin))
        else:
            f_y = (yi[i]-yi[j])/(ymax-ymin)

        rot = atan(f_y/f_x)*360/2/pi

        kw = {}
        kw["ha"] = "center"
        kw["va"] = "center_baseline"
        kw["rotation_mode"] = "anchor"
        kw["rotation"] = rot
        kw["size"] = "small"
        text = grafico.plot.ax.text(xi[i], yi[i], txt, **kw)

        line.text = text
        line.text.pos = pos
        if not annotate:
            text.set_visible(False)


def plot2D3D(grafico, data, Preferences, x, y, z=None):
    """Plot procedure
    Parameters:
        grafico: plot
        data: data to plot
        Preferences: ConfigParser instance from mainwindow preferencesChanged
        x: Key for x axis
        y: Key for y axis
        z: Key for z axis Optional for 3D plot"""

    functionx = _getunitTransform(x)
    functiony = _getunitTransform(y)
    functionz = _getunitTransform(z)
    transform = (functionx, functiony, functionz)

    # Plot saturation lines
    format = getLineFormat(Preferences, "saturation")
    if x == "P" and y == "T":
        satLines = QtWidgets.QApplication.translate(
            "pychemqt", "Saturation Line"),
    else:
        satLines = [
            QtWidgets.QApplication.translate(
                "pychemqt", "Liquid Saturation Line"),
            QtWidgets.QApplication.translate(
                "pychemqt", "Vapor Saturation Line")]
    for fase, label in enumerate(satLines):
        xsat = list(map(functionx, data["saturation_%i" % fase][x]))
        ysat = list(map(functiony, data["saturation_%i" % fase][y]))
        if z:
            zsat = list(map(functionz, data["saturation_%i" % fase][z]))
            grafico.plot.ax.plot(xsat, ysat, zsat, label=label, **format)
        else:
            grafico.plot.ax.plot(xsat, ysat, label=label, **format)

    # Plot melting and sublimation lines
    if "melting" in data:
        label = QtWidgets.QApplication.translate("pychemqt", "Melting Line")
        xmel = list(map(functionx, data["melting"][x]))
        ymel = list(map(functiony, data["melting"][y]))
        if z:
            zmel = list(map(functionz, data["melting"][z]))
            grafico.plot.ax.plot(xmel, ymel, zmel, label=label, **format)
        else:
            grafico.plot.ax.plot(xmel, ymel, label=label, **format)
    if "sublimation" in data:
        xsub = list(map(functionx, data["sublimation"][x]))
        ysub = list(map(functiony, data["sublimation"][y]))
        label = QtWidgets.QApplication.translate(
            "pychemqt", "Sublimation Line")
        if z:
            zmel = list(map(functionz, data["melting"][z]))
            grafico.plot.ax.plot(xmel, ymel, zmel, label=label, **format)
        else:
            grafico.plot.ax.plot(xsub, ysub, label=label, **format)

    # Plot quality isolines
    if x not in ["P", "T"] or y not in ["P", "T"] or z:
        format = getLineFormat(Preferences, "Isoquality")
        plotIsoline(data["x"], (x, y, z), "x", unidades.Dimensionless, grafico,
                    transform, **format)

    # Plot isotherm lines
    if x != "T" and y != "T" or z:
        format = getLineFormat(Preferences, "Isotherm")
        plotIsoline(data["T"], (x, y, z), "T", unidades.Temperature, grafico,
                    transform, **format)

    # Plot isobar lines
    if x != "P" and y != "P" or z:
        format = getLineFormat(Preferences, "Isobar")
        plotIsoline(data["P"], (x, y, z), "P", unidades.Pressure, grafico,
                    transform, **format)

    # Plot isochor lines
    if x not in ["rho", "v"] and y not in ["rho", "v"] or z:
        format = getLineFormat(Preferences, "Isochor")
        plotIsoline(data["v"], (x, y, z), "v", unidades.SpecificVolume,
                    grafico, transform, **format)
        # Plot isodensity lines
        if "rho" in data:
            plotIsoline(data["rho"], (x, y, z), "rho", unidades.Density,
                        grafico, transform, **format)

    # Plot isoenthalpic lines
    if x != "h" and y != "h" or z:
        format = getLineFormat(Preferences, "Isoenthalpic")
        plotIsoline(data["h"], (x, y, z), "h", unidades.Enthalpy, grafico,
                    transform, **format)

    # Plot isoentropic lines
    if x != "s" and y != "s" or z:
        format = getLineFormat(Preferences, "Isoentropic")
        plotIsoline(data["s"], (x, y, z), "s", unidades.SpecificHeat, grafico,
                    transform, **format)


def _getunitTransform(eje):
    """Return the axis unit transform function to map data to configurated unit
        Parameters:
            seq: list with axis property keys
    """
    if not eje:
        return None
    elif eje == "T":
        index = config.getMainWindowConfig().getint("Units", "Temperature")
        func = [float, unidades.K2C, unidades.K2R, unidades.K2F, unidades.K2Re]
        return func[index]
    else:
        unit = meos.units[meos.keys.index(eje)]
        factor = unit(1.).config()
        return lambda val: val*factor if val is not None else nan


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)

    conf = config.getMainWindowConfig()

    # SteamTables = AddPoint(conf)
    # SteamTables=AddLine(None)
    # SteamTables = Dialog(conf)
    SteamTables = Plot3D()

    SteamTables.show()
    sys.exit(app.exec_())
