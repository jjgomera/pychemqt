#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>."""


from configparser import ConfigParser
from functools import partial
import json
import os
import platform
import subprocess
import sys
import time


import equipment
from lib import config
from lib.config import conf_dir, setMainWindowConfig, IMAGE_PATH, Preferences
from lib.project import Project
import plots as charts
from tools.qt import QtCore, QtGui, QtWidgets
from tools import (UI_confComponents, UI_Preferences, UI_confTransport,
                   UI_confThermo, UI_confUnits, UI_databank, UI_Tables,
                   UI_unitConverter, UI_psychrometry, costIndex, doi,
                   dependences, terminal, qtelemental, wizard)
from UI import newComponent, flujo, plots, viewComponents
from UI.prefPFD import BrushCombo
from UI.petro import Definicion_Petro
from UI.widgets import createAction


with open("VERSION") as version_file:
    __version__ = version_file.read().strip()

year = config.__doc__.split()[7][:-1]

other_window = (plots.Binary_distillation, UI_Tables.table.TablaMEoS,
                UI_Tables.PlotMEoS)
other_window_names = [cl.__name__ for cl in other_window]


class TabWidget(QtWidgets.QTabWidget):
    """Custom QTableWidget to populate a empty page with a welcome and
    information message"""

    def paintEvent(self, event):
        """Special paint event"""
        txt = self.tr(
            "Welcome to pychemqt, \n"
            "a software for simulating Chemical Engineering units operations\n"
            "\n"
            f"Copyright © {year} Juan José Gómez Romera (jjgomera)\n"
            "Licenced with GPL.v3\n"
            "This software is distributed in the hope that it will be useful,\n"
            "it is provided 'as is' without warranty of any kind\n"
            "\n"
            "You can start by creating a new project or opening an existing.")

        if self.count():
            QtWidgets.QTabWidget.paintEvent(self, event)
        else:
            painter = QtGui.QPainter(self)
            rect = event.rect()
            image = QtGui.QImage("images/pychemqt.png")
            rectImage = QtCore.QRect(25, rect.center().y()-50, 100, 100)
            painter.drawImage(rectImage, image)
            rect.setLeft(150)

            painter.drawText(rect, QtCore.Qt.AlignmentFlag.AlignVCenter, txt)


class TreeEquipment(QtWidgets.QTreeWidget):
    """Custom QTreWidget to list the equipment and stream from the active
    project"""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setIconSize(QtCore.QSize(30, 30))
        self.headerItem().setHidden(True)
        self.setContextMenuPolicy(QtCore.Qt.ContextMenuPolicy.CustomContextMenu)
        self.Stream = QtWidgets.QTreeWidgetItem(self, 0)
        self.Stream.setText(0, self.tr("Streams"))
        self.Stream.setExpanded(True)
        self.Equipment = QtWidgets.QTreeWidgetItem(self, 0)
        self.Equipment.setText(0, self.tr("Equipments"))
        self.Equipment.setExpanded(True)

    def updateList(self, items):
        """Update list with project items"""
        # self.clear()

        self.Stream.takeChildren()
        self.Equipment.takeChildren()

        ins = []
        outs = []
        for stream in items["in"]:
            for stream in items["in"][stream].down:
                ins.append(stream.id)
        for stream in items["out"]:
            for stream in items["out"][stream].up:
                outs.append(stream.id)

        for key in sorted(items["stream"].keys()):
            idx = items["stream"][key].id
            if idx in ins:
                item = QtWidgets.QTreeWidgetItem(self.Stream, 1)
                item.setText(0, str(idx))
                item.setIcon(0, QtGui.QIcon(QtGui.QPixmap(
                    os.path.join(IMAGE_PATH, "equipment", "in.svg"))))
            elif idx in outs:
                item = QtWidgets.QTreeWidgetItem(self.Stream, 2)
                item.setText(0, str(idx))
                item.setIcon(0, QtGui.QIcon(QtGui.QPixmap(
                    os.path.join(IMAGE_PATH, "equipment", "out.svg"))))
            else:
                item = QtWidgets.QTreeWidgetItem(self.Stream, 3)
                item.setText(0, str(idx))
                item.setIcon(0, QtGui.QIcon(QtGui.QPixmap(
                    os.path.join(IMAGE_PATH, "equipment", "stream.png"))))

        for equip in items["equip"]:
            item = QtWidgets.QTreeWidgetItem(self.Equipment, 4)
            eq = items["equip"][equip]
            item.setText(0, f'{eq.id} - {eq.name}')
            item.setIcon(0, QtGui.QIcon(QtGui.QPixmap(
                items["equip"][equip].imagen)))


class FlowLayout(QtWidgets.QLayout):
    """Dynamic layout to reorder object using all available space"""
    def __init__(self, margin=0, spacing=0, parent=None):
        super().__init__(parent)

        if parent is not None:
            self.setMargin(margin)

        self.setSpacing(spacing)

        self.itemList = []

    def __del__(self):
        item = self.takeAt(0)
        while item:
            item = self.takeAt(0)

    def addItem(self, item):
        """Add item to layout"""
        self.itemList.append(item)

    def count(self):
        """Return count of item"""
        return len(self.itemList)

    def itemAt(self, index):
        """Return item at position index"""
        if 0 <= index < len(self.itemList):
            return self.itemList[index]

        return None

    def takeAt(self, index):
        """Return item at position index and remove from layout"""
        if 0 <= index < len(self.itemList):
            return self.itemList.pop(index)

        return None

    def expandingDirections(self):
        """Limit expanding direction to horizontal orientation"""
        return QtCore.Qt.Orientation.Horizontal

    def hasHeightForWidth(self):
        """Change default value"""
        return True

    def heightForWidth(self, width):
        """Get height from layout"""
        height = self.doLayout(QtCore.QRect(0, 0, width, 0), True)
        return height

    def setGeometry(self, rect):
        """Special geomtry from layout"""
        super().setGeometry(rect)
        self.doLayout(rect, False)

    def sizeHint(self):
        """Set size hint to minimum size"""
        return self.minimumSize()

    def minimumSize(self):
        """Get minimum size of layout"""
        size = QtCore.QSize()

        for item in self.itemList:
            size = size.expandedTo(item.minimumSize())

        size = QtCore.QSize(
            size.width() + self.contentsMargins().left()
            + self.contentsMargins().right(),
            size.height() + self.contentsMargins().bottom()
            + self.contentsMargins().top())
        return size

    def doLayout(self, rect, testOnly):
        """Dynamic generation of layout"""
        x = rect.x()
        y = rect.y()
        lineHeight = 0

        for item in self.itemList:
            wid = item.widget()
            spaceX = self.spacing() + wid.style().layoutSpacing(
                QtWidgets.QSizePolicy.ControlType.PushButton,
                QtWidgets.QSizePolicy.ControlType.PushButton,
                QtCore.Qt.Orientation.Horizontal)*self.spacing()
            spaceY = self.spacing() + wid.style().layoutSpacing(
                QtWidgets.QSizePolicy.ControlType.PushButton,
                QtWidgets.QSizePolicy.ControlType.PushButton,
                QtCore.Qt.Orientation.Vertical)*self.spacing()
            nextX = x + item.sizeHint().width() + spaceX
            if nextX - spaceX > rect.right() and lineHeight > 0:
                x = rect.x()
                y = y + lineHeight + spaceY
                nextX = x + item.sizeHint().width() + spaceX
                lineHeight = 0

            if not testOnly:
                item.setGeometry(QtCore.QRect(
                    QtCore.QPoint(x, y), item.sizeHint()))

            x = nextX
            lineHeight = max(lineHeight, item.sizeHint().height())

        return y + lineHeight - rect.y()


class UI_pychemqt(QtWidgets.QMainWindow):
    """Main window UI definition"""
    idNew = 0  # Count for new projects created in session, for autoname
    config = []
    dirty = []
    filename = []
    pfd = []

    def __init__(self, parent=None):
        super().__init__(parent)
        self.wdg = []
        self.setWindowTitle("pychemqt")
        centralwidget = TabWidget()
        centralwidget.currentChanged.connect(self.currentTabChanged)
        centralwidget.setTabsClosable(True)
        self.setCentralWidget(centralwidget)
        self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(
            IMAGE_PATH + "pychemqt.png")))

        # Actions
        fileNewAction = createAction(
            self.tr("&New"),
            slot=self.fileNew,
            shortcut=QtGui.QKeySequence.StandardKey.New,
            icon=os.path.join("button", "fileNew.png"),
            tip=self.tr("Start new project"),
            parent=self)
        fileOpenAction = createAction(
            self.tr("&Open"),
            slot=self.fileOpen,
            shortcut=QtGui.QKeySequence.StandardKey.Open,
            icon=os.path.join("button", "fileOpen.png"),
            tip=self.tr("Open project"),
            parent=self)
        self.fileSaveAction = createAction(
            self.tr("&Save"),
            slot=self.fileSave,
            shortcut=QtGui.QKeySequence.StandardKey.Save,
            icon=os.path.join("button", "fileSave.png"),
            tip=self.tr("Save project"),
            parent=self)
        self.fileSaveAsAction = createAction(
            self.tr("Save &as"),
            slot=self.fileSaveAs,
            shortcut=QtGui.QKeySequence.StandardKey.SaveAs,
            icon=os.path.join("button", "fileSaveAs.png"),
            tip=self.tr("Save project as"),
            parent=self)
        self.fileSaveAllAction = createAction(
            self.tr("Save A&ll"),
            slot=self.fileSaveAll,
            icon=os.path.join("button", "fileSaveAll.png"),
            tip=self.tr("Save all open project"),
            parent=self)
        self.fileCloseAction = createAction(
            self.tr("&Close"),
            slot=self.fileClose,
            shortcut=QtGui.QKeySequence.StandardKey.Close,
            icon=os.path.join("button", "fileClose.png"),
            tip=self.tr("Close project"),
            parent=self)
        ExitAction = createAction(
            self.tr("&Exit"),
            slot=self.closeEvent,
            shortcut=QtGui.QKeySequence.StandardKey.Quit,
            icon=os.path.join("button", "exit.png"),
            tip=self.tr("Quit pychemqt"),
            parent=self)

        self.actionWizard = createAction(
            self.tr("Wizard"),
            icon=os.path.join("button", "wizard.png"),
            slot=self.wizard,
            tip=self.tr("Launch configuration wizard"),
            parent=self)
        self.actionComponentList = createAction(
            self.tr("Components list"),
            slot=partial(self.dialogConfig, UI_confComponents),
            tip=self.tr("Defining componente list dialog"),
            parent=self)
        self.actionThermo = createAction(
            self.tr("Thermodynamic properties"),
            slot=partial(self.dialogConfig, UI_confThermo),
            tip=self.tr("Defining thermodynamic properties methods"),
            parent=self)
        self.actionTransporte = createAction(
            self.tr("Transport properties"),
            slot=partial(self.dialogConfig, UI_confTransport),
            tip=self.tr("Defining transport properties methods"),
            parent=self)
        self.actionUnidades = createAction(
            self.tr("Units"),
            slot=partial(self.dialogConfig, UI_confUnits),
            tip=self.tr("Defining preferred units"),
            parent=self)
        self.actioncostIndex = createAction(
            self.tr("&Cost Index"),
            slot=partial(self.launch, costIndex.Ui_CostIndex),
            tip=self.tr("Defining cost index"),
            parent=self)
        self.actionPreferencias = createAction(
            self.tr("&Preferences"),
            slot=self.Preferences,
            icon=os.path.join("button", "configure.png"),
            shortcut=QtGui.QKeySequence.StandardKey.Preferences,
            tip=self.tr("Defining general preferences"),
            parent=self)

        self.actionZoomIn = createAction(
            self.tr("Zoom in"),
            slot=partial(self.zoom, "+"),
            icon=os.path.join("button", "zoomIn.png"),
            shortcut=QtGui.QKeySequence.StandardKey.ZoomIn, parent=self)
        self.actionZoom = createAction(
            self.tr("Zoom"),
            slot=partial(self.zoom, "Dialog"),
            icon=os.path.join("button", "zoomIn.png"),
            parent=self)
        self.actionZoomOut = createAction(
            self.tr("Zoom out"),
            slot=partial(self.zoom, "-"),
            icon=os.path.join("button", "zoomOut.png"),
            shortcut=QtGui.QKeySequence.StandardKey.ZoomOut,
            parent=self)
        self.actionOverviewWindow = createAction(
            self.tr("Overview window"),
            slot=self.overview, checkable=True, parent=self)
        actionVerStatus = createAction(
            self.tr("Status"),
            shortcut="Ctrl+Alt+S",
            tip=self.tr("Show/Hide status toolbar"),
            checkable=True, parent=self)
        self.actionVerToolbar = createAction(
            self.tr("Palette"),
            shortcut="Ctrl+Alt+T",
            icon=os.path.join("button", "palette.png"),
            tip=self.tr("Show/Hide equipment palette"),
            checkable=True, parent=self)
        self.actionVerItem = createAction(
            self.tr("Item"),
            icon=os.path.join("button", "list.png"),
            tip=self.tr("Show/Hide item list"),
            checkable=True, parent=self)

        self.calculatorAction = createAction(
            self.tr("&Calculator"),
            slot=self.calculator,
            icon=os.path.join("button", "calculator.png"),
            shortcut="F2",
            tip=self.tr("Open system calculator"),
            parent=self)

        terminalAction = createAction(
            self.tr("Python Shell"),
            shortcut="F3",
            slot=self.terminal,
            icon=os.path.join("button", "terminal.png"),
            tip=self.tr("Open system terminal"),
            parent=self)
        if sys.platform == "win32":
            terminalAction.setEnabled(False)

        conversorUnidadesAction = createAction(
            self.tr("&Units converter"),
            slot=partial(self.launch, UI_unitConverter.UI_unitConverter),
            shortcut="F4",
            icon=os.path.join("button", "unitConverter.png"),
            tip=self.tr("Open Units converter"),
            parent=self)
        currencyAction = createAction(
            self.tr("&Currency converter"),
            slot=partial(self.launch, UI_unitConverter.moneda),
            shortcut="F5",
            icon=os.path.join("button", "currency.png"),
            tip=self.tr("Open Currency converter"),
            parent=self)
        TablaPeriodicaAction = createAction(
            self.tr("&Periodic Table"),
            slot=partial(self.launch, qtelemental.qtelemental),
            shortcut="F6",
            icon=os.path.join("button", "periodicTable.png"),
            tip=self.tr("Show a basic Mendeleiev periodic table"),
            parent=self)
        TablesAction = createAction(
            self.tr("MEOS"),
            slot=self.meos,
            shortcut="F7",
            icon=os.path.join("button", "tables.png"),
            tip=self.tr("Open a advanced thermodynamic properties application"),
            parent=self)
        psychrometricChartAction = createAction(
            self.tr("&Psicrometric Chart"),
            slot=partial(self.launch, UI_psychrometry.UI_Psychrometry),
            shortcut="F9",
            icon=os.path.join("button", "psychrometric.png"),
            tip=self.tr("Open a humid-air application"),
            parent=self)
        externalProgramAction = createAction(
            self.tr("External Programs"),
            slot=partial(self.launch, dependences.ShowDependences),
            icon=os.path.join("button", "showPrograms.png"),
            tip=self.tr("Show External Programs Status"),
            parent=self)
        self.saveAsImage = createAction(
            self.tr("Save PFD as image"),
            slot=self.savePFDImage,
            icon=os.path.join("button", "image.png"),
            parent=self)
        actionAyuda = createAction(
            self.tr("Help"),
            slot=self.help,
            icon=os.path.join("button", "help.png"),
            parent=self)
        actionDocum = createAction(
            self.tr("References"),
            slot=self.documentation,
            parent=self)
        actionLog = createAction(
            self.tr("View Log"),
            slot=self.log,
            parent=self)
        actionAcerca_de = createAction(
            self.tr("About pychemqt"),
            slot=self.acerca,
            icon=os.path.join("button/helpAbout.png"),
            parent=self)
        actionAcerca_deQt = createAction(
            self.tr("About Qt"),
            slot=self.acercaQt,
            icon=os.path.join("button", "AboutQt.png"),
            parent=self)

        self.zoomValue = QtWidgets.QSpinBox()
        self.zoomValue.setSuffix("%")
        self.zoomValue.setRange(5, 1000)
        self.zoomValue.setValue(100)
        self.zoomValue.setSingleStep(5)
        self.zoomValue.valueChanged.connect(self.zoom)

        # Toolbars
        self.BarraArchivo = QtWidgets.QToolBar(
            self.tr("File"), self)
        self.BarraArchivo.setObjectName("BarraArchivo")
        self.BarraArchivo.setIconSize(QtCore.QSize(16, 16))
        self.BarraArchivo.addAction(fileNewAction)
        self.BarraArchivo.addAction(fileOpenAction)
        self.BarraArchivo.addAction(self.fileCloseAction)
        self.BarraArchivo.addAction(self.fileSaveAction)
        self.BarraArchivo.addAction(self.fileSaveAsAction)
        self.BarraArchivo.addAction(self.fileSaveAllAction)
        self.BarraArchivo.addAction(ExitAction)
        self.addToolBar(QtCore.Qt.ToolBarArea.TopToolBarArea, self.BarraArchivo)

        self.BarraVer = QtWidgets.QToolBar(
            self.tr("View"), self)
        self.BarraVer.setObjectName("BarraVer")
        self.BarraVer.setIconSize(QtCore.QSize(16, 16))
        self.BarraVer.addAction(self.actionZoomOut)
        self.BarraVer.addWidget(self.zoomValue)
        self.BarraVer.addAction(self.actionZoomIn)
        self.BarraVer.addSeparator()
        self.BarraVer.addAction(self.actionVerToolbar)
        self.BarraVer.addAction(self.actionVerItem)
        self.addToolBar(QtCore.Qt.ToolBarArea.TopToolBarArea, self.BarraVer)

        self.BarraHerramientas = QtWidgets.QToolBar(
            self.tr("Tools"), self)
        self.BarraHerramientas.setObjectName("BarraHerramientas")
        self.BarraHerramientas.setIconSize(QtCore.QSize(16, 16))
        self.BarraHerramientas.addAction(self.calculatorAction)
        self.BarraHerramientas.addAction(terminalAction)
        self.BarraHerramientas.addAction(conversorUnidadesAction)
        self.BarraHerramientas.addAction(currencyAction)
        self.BarraHerramientas.addAction(TablaPeriodicaAction)
        self.BarraHerramientas.addAction(TablesAction)
        self.BarraHerramientas.addAction(psychrometricChartAction)
        self.BarraHerramientas.addAction(externalProgramAction)
        self.BarraHerramientas.addSeparator()
        self.BarraHerramientas.addAction(actionAyuda)
        self.addToolBar(
            QtCore.Qt.ToolBarArea.TopToolBarArea, self.BarraHerramientas)

        # Paleta Toolbox
        self.toolboxPalette = QtWidgets.QDockWidget(
            self.tr("Equipment"))
        self.toolboxPalette.setObjectName("toolbox")
        toolboxContenido = QtWidgets.QWidget()
        self.toolboxPalette.setWidget(toolboxContenido)
        self.addDockWidget(QtCore.Qt.DockWidgetArea(1), self.toolboxPalette)
        self.toolboxPalette.setFeatures(
            QtWidgets.QDockWidget.DockWidgetFeature.DockWidgetClosable
            | QtWidgets.QDockWidget.DockWidgetFeature.DockWidgetMovable
            | QtWidgets.QDockWidget.DockWidgetFeature.DockWidgetFloatable)
        self.toolboxPalette.setAllowedAreas(
            QtCore.Qt.DockWidgetArea.LeftDockWidgetArea
            | QtCore.Qt.DockWidgetArea.RightDockWidgetArea)

        self.toolboxPalette.visibilityChanged.connect(
            self.actionVerToolbar.setChecked)
        self.actionVerToolbar.triggered.connect(self.toolboxPalette.setVisible)
        layouttoolbox = QtWidgets.QVBoxLayout(toolboxContenido)
        layouttoolbox.setContentsMargins(5, 5, 5, 5)

        txt = QtWidgets.QLabel(
            self.tr("Plot"))
        txt.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        layouttoolbox.addWidget(txt)
        l1 = FlowLayout(spacing=10)
        actionTexto, botonTexto = createAction(
            self.tr("Insert text"),
            icon=os.path.join("equipment", "text.png"),
            slot=self.addText,
            button=True, parent=toolboxContenido)
        l1.addWidget(botonTexto)
        actionCuadrado, botonCuadrado = createAction(
            self.tr("Draw square"),
            icon=os.path.join("equipment", "square.png"),
            slot=partial(self.addItem, "square"),
            button=True, parent=toolboxContenido)
        l1.addWidget(botonCuadrado)
        actionCircle, botonCircle = createAction(
            self.tr("Draw circle"),
            icon=os.path.join("equipment", "circle.png"),
            slot=partial(self.addItem, "ellipse"),
            button=True, parent=toolboxContenido)
        l1.addWidget(botonCircle)
        layouttoolbox.addItem(l1)
        layouttoolbox.addStretch(1)

        txt = QtWidgets.QLabel(self.tr("Flux"))
        txt.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        layouttoolbox.addWidget(txt)
        l2 = FlowLayout(spacing=10)
        actionEntrada, botonEntrada = createAction(
            self.tr("Input"),
            icon=os.path.join("equipment", "in.png"),
            slot=partial(self.addItem, "in"),
            button=True, parent=toolboxContenido)
        l2.addWidget(botonEntrada)
        actionCorriente, self.botonCorriente = createAction(
            self.tr("Stream"),
            icon=os.path.join("equipment", "stream.png"),
            slot=partial(self.addItem, "stream"),
            button=True, checkable=True, parent=toolboxContenido)
        l2.addWidget(self.botonCorriente)
        actionSalida, botonSalida = createAction(
            self.tr("Output"),
            icon=os.path.join("equipment", "out.png"),
            slot=partial(self.addItem, "out"),
            button=True, parent=toolboxContenido)
        l2.addWidget(botonSalida)
        actionDivider, botonDivider = createAction(
            self.tr("Divider"),
            icon=os.path.join("equipment", "divider.png"),
            slot=partial(self.addEquipment, equipment.UI_divider),
            button=True, parent=toolboxContenido)
        l2.addWidget(botonDivider)
        actionValve, botonValve = createAction(
            self.tr("Valve"),
            icon=os.path.join("equipment", "valve.png"),
            slot=partial(self.addEquipment, equipment.UI_valve),
            button=True, parent=toolboxContenido)
        l2.addWidget(botonValve)
        actionMixer, botonMixer = createAction(
            self.tr("Mixer"),
            icon=os.path.join("equipment", "mixer.png"),
            slot=partial(self.addEquipment, equipment.UI_mixer),
            button=True, parent=toolboxContenido)
        l2.addWidget(botonMixer)
        actionCompresor, botonCompresor = createAction(
            self.tr("Compressor"),
            icon=os.path.join("equipment", "compressor.png"),
            slot=partial(self.addEquipment, equipment.UI_compressor),
            button=True, parent=toolboxContenido)
        l2.addWidget(botonCompresor)
        actionTurbine, botonTurbine = createAction(
            self.tr("Turbine"),
            icon=os.path.join("equipment", "turbine.png"),
            slot=partial(self.addEquipment, equipment.UI_turbine),
            button=True, parent=toolboxContenido)
        l2.addWidget(botonTurbine)
        actionPump, botonPump = createAction(
            self.tr("Pump"),
            icon=os.path.join("equipment", "pump.png"),
            slot=partial(self.addEquipment, equipment.UI_pump),
            button=True, parent=toolboxContenido)
        l2.addWidget(botonPump)
        actionPipe, botonPipe = createAction(
            self.tr("Pipe"),
            icon=os.path.join("equipment", "pipe.png"),
            slot=partial(self.addEquipment, equipment.UI_pipe),
            button=True, parent=toolboxContenido)
        l2.addWidget(botonPipe)
        layouttoolbox.addItem(l2)
        layouttoolbox.addStretch(1)

        txt = QtWidgets.QLabel(
            self.tr("Basics"))
        txt.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        layouttoolbox.addWidget(txt)
        l3 = FlowLayout(spacing=10)
        actionTorreFUG, botonTorreFUG = createAction(
            self.tr("Distillation tower (method FUG)"),
            icon=os.path.join("equipment", "columnFUG.png"),
            slot=partial(self.addEquipment, equipment.UI_columnFUG),
            button=True, parent=toolboxContenido)
        l3.addWidget(botonTorreFUG)
        actionFlash, botonFlash = createAction(
            self.tr("Flash"),
            icon=os.path.join("equipment", "flash.png"),
            slot=partial(self.addEquipment, equipment.UI_flash),
            button=True, parent=toolboxContenido)
        l3.addWidget(botonFlash)
        actionTorre, botonTorre = createAction(
            self.tr("Distillation tower (exact method)"),
            icon=os.path.join("equipment", "tower.png"),
            slot=partial(self.addEquipment, equipment.UI_tower),
            button=True, parent=toolboxContenido)
        l3.addWidget(botonTorre)
        botonTorre.setEnabled(False)
        actionheatExchanger, botonheatExchanger = createAction(
            self.tr("Generic heat exchanger"),
            icon=os.path.join("equipment", "heatExchanger"),
            slot=partial(self.addEquipment, equipment.UI_heatExchanger),
            button=True, parent=toolboxContenido)
        l3.addWidget(botonheatExchanger)
        actionhairpin, botonhairpin = createAction(
            self.tr("Hairpin heat exchanger"),
            icon=os.path.join("equipment", "hairpin.png"),
            slot=partial(self.addEquipment, equipment.UI_hairpin),
            button=True, parent=toolboxContenido)
        l3.addWidget(botonhairpin)
        actionShellTube, botonShellTube = createAction(
            self.tr("Shell and tube heat exchanger"),
            icon=os.path.join("equipment", "shellTube.png"),
            slot=partial(self.addEquipment, equipment.UI_shellTube),
            button=True, parent=toolboxContenido)
        l3.addWidget(botonShellTube)
        actionFireHeater, botonFireHeater = createAction(
            self.tr("Fired Heater heat exchanger"),
            icon=os.path.join("equipment", "fireHeater.png"),
            slot=partial(self.addEquipment, equipment.UI_fireHeater),
            button=True, parent=toolboxContenido)
        l3.addWidget(botonFireHeater)
        actionReactor, botonReactor = createAction(
            self.tr("Reactor"),
            icon=os.path.join("equipment", "reactor.png"),
            slot=partial(self.addEquipment, equipment.UI_reactor),
            button=True, parent=toolboxContenido)
        l3.addWidget(botonReactor)
        layouttoolbox.addItem(l3)
        layouttoolbox.addStretch(1)

        txt = QtWidgets.QLabel(
            self.tr("Solids"))
        txt.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        layouttoolbox.addWidget(txt)
        l4 = FlowLayout(spacing=10)
        actionBaghouse, botonBaghouse = createAction(
            self.tr("Baghouse"),
            icon=os.path.join("equipment", "baghouse.png"),
            slot=partial(self.addEquipment, equipment.UI_baghouse),
            button=True, parent=toolboxContenido)
        l4.addWidget(botonBaghouse)
        actionCentrifuge, botonCentrifuge = createAction(
            self.tr("Centrifuge"),
            icon=os.path.join("equipment", "centrifuge.png"),
            slot=partial(self.addEquipment, equipment.UI_centrifuge),
            button=True, parent=toolboxContenido)
        l4.addWidget(botonCentrifuge)
        botonCentrifuge.setEnabled(False)
        actionCiclon, botonCiclon = createAction(
            self.tr("Cyclone"),
            icon=os.path.join("equipment", "ciclon.png"),
            slot=partial(self.addEquipment, equipment.UI_ciclon),
            button=True, parent=toolboxContenido)
        l4.addWidget(botonCiclon)
        actionElectroPrecipitator, botonElectroPrecipitator = createAction(
            self.tr("Electric precipitator"),
            icon=os.path.join("equipment", "electricPrecipitator"),
            slot=partial(self.addEquipment, equipment.UI_electricPrecipitator),
            button=True, parent=toolboxContenido)
        l4.addWidget(botonElectroPrecipitator)
        actionGrinder, botonGrinder = createAction(
            self.tr("Grinder"),
            icon=os.path.join("equipment", "grinder.png"),
            slot=partial(self.addEquipment, equipment.UI_grinder),
            button=True, parent=toolboxContenido)
        l4.addWidget(botonGrinder)
        actionDryer, botonDryer = createAction(
            self.tr("Solids dryer"),
            icon=os.path.join("equipment", "dryer.png"),
            slot=partial(self.addEquipment, equipment.UI_dryer),
            button=True, parent=toolboxContenido)
        l4.addWidget(botonDryer)
        actionWasher, botonWasher = createAction(
            self.tr("Solid washer"),
            icon=os.path.join("equipment", "solidWasher.png"),
            slot=partial(self.addEquipment, equipment.UI_solidWasher),
            button=True, parent=toolboxContenido)
        l4.addWidget(botonWasher)
        botonWasher.setEnabled(False)
        actionVacuum, botonVacuum = createAction(
            self.tr("Vacuum filter"),
            icon=os.path.join("equipment", "vacuumfilter.png"),
            slot=partial(self.addEquipment, equipment.UI_vacuum),
            button=True, parent=toolboxContenido)
        l4.addWidget(botonVacuum)
        botonVacuum.setEnabled(False)
        actionScrubber, botonScrubber = createAction(
            self.tr("Scrubber"),
            icon=os.path.join("equipment", "scrubber.png"),
            slot=partial(self.addEquipment, equipment.UI_scrubber),
            button=True, parent=toolboxContenido)
        l4.addWidget(botonScrubber)
        actionNeumatic, botonNeumatic = createAction(
            self.tr("Neumatic"),
            icon=os.path.join("equipment", "neumatic.png"),
            slot=partial(self.addEquipment, equipment.UI_neumatic),
            button=True, parent=toolboxContenido)
        l4.addWidget(botonNeumatic)
        actionGravityChandler, botonGravityChandler = createAction(
            self.tr("Gravity settling chamber"),
            icon=os.path.join("equipment", "gravityChamber.png"),
            slot=partial(self.addEquipment, equipment.UI_gravityChamber),
            button=True, parent=toolboxContenido)
        l4.addWidget(botonGravityChandler)
        layouttoolbox.addItem(l4)
        layouttoolbox.addStretch(1)

        txt = QtWidgets.QLabel(self.tr("Tools"))
        txt.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        layouttoolbox.addWidget(txt)

        l5 = FlowLayout()
        actionControler, botonControler = createAction(
            self.tr("PID controller"),
            icon=os.path.join("equipment", "controller.png"),
            slot=partial(self.addEquipment, equipment.UI_solidWasher),
            button=True, parent=toolboxContenido)
        l5.addWidget(botonControler)
        botonControler.setEnabled(False)
        actionControlValve, botonControlValve = createAction(
            self.tr("Control valve"),
            icon=os.path.join("equipment", "controlvalve.png"),
            slot=partial(self.addEquipment, equipment.UI_vacuum),
            button=True, parent=toolboxContenido)
        l5.addWidget(botonControlValve)
        botonControlValve.setEnabled(False)
        actionSpreadsheet, botonSpreadsheet = createAction(
            self.tr("External spreadsheet module"),
            icon=os.path.join("equipment", "spreadsheet.png"),
            slot=partial(self.addEquipment, equipment.UI_spreadsheet),
            button=True,
            parent=toolboxContenido)
        if os.environ["ezodf"] != "True" and os.environ["openpyxl"] != "True" \
                and os.environ["xlwt"] != "True":
            actionSpreadsheet.setEnabled(False)
            botonSpreadsheet.setEnabled(False)
        l5.addWidget(botonSpreadsheet)
        layouttoolbox.addItem(l5)
        layouttoolbox.addStretch(10)

        # Menus
        self.menubar = QtWidgets.QMenuBar()
        self.setMenuBar(self.menubar)

        self.menuArchivo = QtWidgets.QMenu(self.tr("&File"))
        self.menuRecentFiles = QtWidgets.QMenu(
            self.tr("Open Recent Files"),
            self.menuArchivo)
        self.menuRecentFiles.aboutToShow.connect(
            self.aboutToShow_MenuRecentFiles)

        self.menuArchivo.addAction(fileNewAction)
        self.menuArchivo.addAction(fileOpenAction)
        self.menuArchivo.addAction(self.fileSaveAction)
        self.menuArchivo.addAction(self.fileSaveAsAction)
        self.menuArchivo.addAction(self.fileSaveAllAction)
        self.menuArchivo.addAction(self.fileCloseAction)
        self.menuArchivo.addSeparator()
        self.menuArchivo.addAction(self.menuRecentFiles.menuAction())
        self.menuArchivo.addSeparator()
        self.menuArchivo.addAction(ExitAction)
        self.menubar.addAction(self.menuArchivo.menuAction())

        self.menuEditar = QtWidgets.QMenu(self.tr("&Edit"))
        self.menuEditar.aboutToShow.connect(self.aboutToShow_MenuEdit)
        self.menubar.addAction(self.menuEditar.menuAction())

        self.menuVer = QtWidgets.QMenu(self.tr("&View"))
        self.menuVer.addAction(self.actionZoomOut)
        self.menuVer.addAction(self.actionZoom)
        self.menuVer.addAction(self.actionZoomIn)
        self.menuVer.addSeparator()
        self.menuVer.addAction(self.actionOverviewWindow)
        self.menuVer.addAction(actionVerStatus)
        self.menuVer.addAction(self.actionVerToolbar)
        self.menuVer.addAction(self.actionVerItem)
        self.menubar.addAction(self.menuVer.menuAction())

        self.menuObjetos = QtWidgets.QMenu(self.tr("Insert component"))
        self.menuObjetosGraficos = QtWidgets.QMenu(self.tr("Plot"))
        self.menuObjetosGraficos.addAction(actionTexto)
        self.menuObjetosGraficos.addAction(actionCuadrado)
        self.menuObjetosGraficos.addAction(actionCircle)
        self.menuObjetos.addAction(self.menuObjetosGraficos.menuAction())

        self.menuObjetosFlujo = QtWidgets.QMenu(self.tr("Flux"))
        self.menuObjetosFlujo.addAction(actionEntrada)
        self.menuObjetosFlujo.addAction(actionCorriente)
        self.menuObjetosFlujo.addAction(actionSalida)
        self.menuObjetosFlujo.addAction(actionDivider)
        self.menuObjetosFlujo.addAction(actionPipe)
        self.menuObjetosFlujo.addAction(actionMixer)
        self.menuObjetosFlujo.addAction(actionCompresor)
        self.menuObjetosFlujo.addAction(actionTurbine)
        self.menuObjetosFlujo.addAction(actionPump)
        self.menuObjetosFlujo.addAction(actionValve)
        self.menuObjetos.addAction(self.menuObjetosFlujo.menuAction())

        self.menuObjetosBasics = QtWidgets.QMenu(self.tr("Basics"))
        self.menuObjetosBasics.addAction(actionTorreFUG)
        self.menuObjetosBasics.addAction(actionFlash)
        self.menuObjetosBasics.addAction(actionTorre)
        self.menuObjetosBasics.addAction(actionheatExchanger)
        self.menuObjetosBasics.addAction(actionhairpin)
        self.menuObjetosBasics.addAction(actionShellTube)
        self.menuObjetosBasics.addAction(actionFireHeater)
        self.menuObjetosBasics.addAction(actionReactor)
        self.menuObjetos.addAction(self.menuObjetosBasics.menuAction())

        self.menuObjetosSolids = QtWidgets.QMenu(self.tr("Solids"))
        self.menuObjetosSolids.addAction(actionBaghouse)
        self.menuObjetosSolids.addAction(actionCentrifuge)
        self.menuObjetosSolids.addAction(actionCiclon)
        self.menuObjetosSolids.addAction(actionElectroPrecipitator)
        self.menuObjetosSolids.addAction(actionGrinder)
        self.menuObjetosSolids.addAction(actionDryer)
        self.menuObjetosSolids.addAction(actionWasher)
        self.menuObjetosSolids.addAction(actionVacuum)
        self.menuObjetosSolids.addAction(actionScrubber)
        self.menuObjetosSolids.addAction(actionNeumatic)
        self.menuObjetosSolids.addAction(actionGravityChandler)
        self.menuObjetos.addAction(self.menuObjetosSolids.menuAction())

        self.menuObjetosTools = QtWidgets.QMenu(self.tr("Tools"))
        self.menuObjetosTools.addAction(actionControler)
        self.menuObjetosTools.addAction(actionControlValve)
        self.menuObjetosTools.addAction(actionSpreadsheet)
        self.menuObjetos.addAction(self.menuObjetosTools.menuAction())

        self.menuPFD = QtWidgets.QMenu(self.tr("&PFD"))
        self.menuPFD.aboutToShow.connect(self.aboutToShow_MenuPFD)
        self.menubar.addAction(self.menuPFD.menuAction())

        # self.menuPlot = QtWidgets.QMenu(self.tr("Pl&ot"))
        # for indice, grafico in enumerate(plots._all):
            # self.menuPlot.addAction(grafico.title, partial(self.plot, indice))
        # self.menubar.addAction(self.menuPlot.menuAction())

        self.menuCharts = QtWidgets.QMenu(self.tr("Charts"), self)
        for titulo, lista in charts._all.items():
            menu = QtWidgets.QMenu(titulo, self)
            for grafico in lista:
                menu.addAction(grafico.title, partial(self.launch, grafico))
            self.menuCharts.addAction(menu.menuAction())

        self.menuAddComponent = QtWidgets.QMenu(self.tr("New Component"))
        self.menuAddComponent.addAction(
            self.tr("Component"),
            partial(self.launch, viewComponents.View_Component))
        self.menuAddComponent.addAction(
            self.tr("Pseudocomponent"),
            partial(self.launch, Definicion_Petro))
        self.menuAddComponent.addSeparator()

        # Group contribution new component menu
        self.menuNewComponent = QtWidgets.QMenu(self.tr("Group contribution"))
        self.menuAddComponent.addAction(self.menuNewComponent.menuAction())
        for f in newComponent._methods:
            self.menuNewComponent.addAction(
                f.__title__,
                partial(self.launch, newComponent.Ui_Contribution, f.__name__))

        self.menuHerramientas = QtWidgets.QMenu(self.tr("&Tools"))
        self.menuHerramientas.addAction(
            self.tr("Component database"),
            partial(self.launch, UI_databank.UI_databank))
        self.menuHerramientas.addAction(self.menuAddComponent.menuAction())
        self.menuHerramientas.addSeparator()
        self.menuHerramientas.addAction(self.calculatorAction)
        self.menuHerramientas.addAction(terminalAction)
        self.menuHerramientas.addAction(conversorUnidadesAction)
        self.menuHerramientas.addAction(currencyAction)
        self.menuHerramientas.addAction(TablaPeriodicaAction)
        self.menuMEoS = UI_Tables.Menu(parent=self)
        self.menuHerramientas.addAction(self.menuMEoS.menuAction())
        self.menuHerramientas.addAction(psychrometricChartAction)
        self.menuHerramientas.addSeparator()
        self.menuHerramientas.addAction(self.menuCharts.menuAction())
        self.menubar.addAction(self.menuHerramientas.menuAction())

        self.menuVentana = QtWidgets.QMenu(self.tr("&Window"))
        self.menuVentana.aboutToShow.connect(self.aboutToShow_MenuWindow)
        self.menubar.addAction(self.menuVentana.menuAction())

        self.menuAyuda = QtWidgets.QMenu(self.tr("&Help"))
        self.menuAyuda.addAction(actionAyuda)
        self.menuAyuda.addAction(actionLog)
        self.menuAyuda.addSeparator()
        self.menuAyuda.addAction(actionDocum)
        self.menuAyuda.addAction(externalProgramAction)
        self.menuAyuda.addSeparator()
        self.menuAyuda.addAction(actionAcerca_de)
        self.menuAyuda.addAction(actionAcerca_deQt)
        self.menubar.addAction(self.menuAyuda.menuAction())

        # Toolbox ListEquipment
        self.toolboxItem = QtWidgets.QDockWidget(self.tr("Item"))
        self.toolboxItem.setObjectName("item")
        self.list = TreeEquipment()
        self.list.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection)
        self.list.itemSelectionChanged.connect(self.selectionChanged)
        self.list.customContextMenuRequested.connect(self.contextListMenu)
        self.toolboxItem.setWidget(self.list)
        self.addDockWidget(QtCore.Qt.DockWidgetArea(1), self.toolboxItem)
        self.toolboxItem.setFeatures(
            QtWidgets.QDockWidget.DockWidgetFeature.DockWidgetClosable
            | QtWidgets.QDockWidget.DockWidgetFeature.DockWidgetMovable
            | QtWidgets.QDockWidget.DockWidgetFeature.DockWidgetFloatable)
        self.toolboxItem.setAllowedAreas(
            QtCore.Qt.DockWidgetArea.LeftDockWidgetArea
            | QtCore.Qt.DockWidgetArea.RightDockWidgetArea)
        self.toolboxItem.visibilityChanged.connect(
            self.actionVerItem.setChecked)
        self.actionVerItem.triggered.connect(self.toolboxItem.setVisible)

        # Toolbox Status
        toolbox = QtWidgets.QDockWidget(self.tr("Status"))
        toolbox.setObjectName("status")
        self.status = QtWidgets.QTextEdit()
        self.status.setReadOnly(True)
        toolbox.setWidget(self.status)
        self.addDockWidget(QtCore.Qt.DockWidgetArea(8), toolbox)
        toolbox.setFeatures(
            QtWidgets.QDockWidget.DockWidgetFeature.NoDockWidgetFeatures)
        toolbox.setAllowedAreas(QtCore.Qt.DockWidgetArea.BottomDockWidgetArea)
        toolbox.visibilityChanged.connect(actionVerStatus.setChecked)
        actionVerStatus.triggered.connect(toolbox.setVisible)

        # Statusbar
        statusbar = QtWidgets.QStatusBar(self)
        statusbar.setMaximumHeight(20)
        statusbar.setSizeGripEnabled(True)
        self.progressBar = QtWidgets.QProgressBar()
        self.progressBar.setVisible(False)
        self.progressBar.setFixedWidth(80)
        statusbar.addPermanentWidget(self.progressBar)
        self.setStatusBar(statusbar)

        # TrayIcon
        self.systemtray = QtWidgets.QSystemTrayIcon(
            QtGui.QIcon(QtGui.QPixmap(IMAGE_PATH + "pychemqt.png")), self)
        self.systemtray.setToolTip("pychemqt")
        self.systemtray.setContextMenu(self.menuHerramientas)

        # Load settings
        settings = QtCore.QSettings()
        self.recentFiles = settings.value("RecentFiles")
        self.lastFile = settings.value("LastFile")
        if settings.value("Geometry"):
            self.restoreGeometry(settings.value("Geometry"))
        else:
            self.showMaximized()
        if settings.value("MainWindow/State"):
            self.restoreState(settings.value("MainWindow/State"))
        if self.recentFiles is None:
            self.recentFiles = []
        self.menuRecentFiles.setEnabled(bool(self.recentFiles))

        self.updateStatus("Loaded pychemqt")
        self.activeControl(False)
        self.changePreferenceLive()
        self.centralWidget().tabCloseRequested.connect(self.fileClose)

    @property
    def currentScene(self):
        """Get current graphicscene"""
        if not self.centralWidget().count():
            return None
        return self.currentView.scene()

    @property
    def currentView(self):
        """Get current graphicsview"""
        if not self.centralWidget().count():
            return None
        return self.centralWidget().currentWidget().subWindowList()[0].widget()

    @property
    def currentMdi(self):
        """Get mdi area of current project"""
        if not self.centralWidget().count():
            return None
        return self.centralWidget().currentWidget()

    @property
    def currentConfig(self):
        """Get config of current project"""
        if not self.centralWidget().count():
            return None
        return self.config[self.idTab]

    @property
    def currentFilename(self):
        """Get file name of current project"""
        if not self.centralWidget().count():
            return None
        return self.filename[self.idTab]

    def getScene(self, index):
        """Get graphicscene of tab with selected index"""
        if index > self.centralWidget().count():
            return None
        return self.getView(index).scene()

    def getView(self, index):
        """Get graphics view of tab with selected index"""
        if index > self.centralWidget().count():
            return None
        return self.centralWidget().widget(index).subWindowList()[0].widget()

    @property
    def idTab(self):
        """Get id of current tab"""
        if not self.centralWidget().count():
            return None
        return self.centralWidget().currentIndex()

    def closeEvent(self, event=None):
        """Catch main close event to do actions like avoid lost changes or
        save main program settings"""
        if self.okToContinue():
            for tab in range(self.centralWidget().count()):
                centralWidget = self.centralWidget().widget(tab)
                scene = centralWidget.subWindowList()[0].widget().scene()
                scene.clearSelection()

            settings = QtCore.QSettings()
            if self.filename:
                filename = QtCore.QVariant(self.filename)
            else:
                filename = QtCore.QVariant()
            settings.setValue("LastFile", filename)
            if self.recentFiles:
                recentFiles = QtCore.QVariant(self.recentFiles)
            else:
                recentFiles = QtCore.QVariant()
            settings.setValue("RecentFiles", recentFiles)
            settings.setValue("Geometry", QtCore.QVariant(self.saveGeometry()))
            settings.setValue("MainWindow/State",
                              QtCore.QVariant(self.saveState()))
            self.close()
        else:
            event.ignore()

    def okToContinue(self, ind=-1):
        """Ask user to save project when close to avoid lost changes"""
        if not self.dirty:
            return True

        if ind != -1:
            ind = list(range(self.centralWidget().count()))
        else:
            ind = [ind]
        dirty = False
        for tab in ind:
            if self.dirty[tab]:
                dirty = True
                break

        if dirty:
            dialog = QtWidgets.QMessageBox.question(
                self,
                self.tr("Unsaved changes"),
                self.tr("Save unsaved changes?"),
                QtWidgets.QMessageBox.StandardButton.Yes
                | QtWidgets.QMessageBox.StandardButton.No
                | QtWidgets.QMessageBox.StandardButton.Cancel,
                QtWidgets.QMessageBox.StandardButton.Yes)

            if dialog == QtWidgets.QMessageBox.StandardButton.Cancel:
                return False
            if dialog == QtWidgets.QMessageBox.StandardButton.No:
                return True
            if dialog == QtWidgets.QMessageBox.StandardButton.Yes:
                self.fileSaveAll()
                return True
        else:
            return True

    # Menus configuration
    def aboutToShow_MenuEdit(self):
        """Dynamic generation of edit menu"""
        self.menuEditar.clear()
        self.menuEditar.addAction(self.actionWizard)
        self.menuEditar.addSeparator()
        self.menuEditar.addAction(self.actionComponentList)
        self.menuEditar.addAction(self.actionThermo)
        self.menuEditar.addAction(self.actionTransporte)
        self.menuEditar.addAction(self.actionUnidades)
        self.menuEditar.addSeparator()
        self.menuEditar.addAction(self.actioncostIndex)
        self.menuEditar.addSeparator()
        self.menuEditar.addAction(self.actionPreferencias)

    def aboutToShow_MenuPFD(self):
        """Dynamic generation of PFD menu"""
        self.menuPFD.clear()
        if self.currentScene:
            self.currentScene.addActions(self.menuPFD)

        self.menuPFD.addAction(self.menuObjetos.menuAction())
        self.menuPFD.addSeparator()
        self.menuPFD.addAction(self.saveAsImage)

    def aboutToShow_MenuWindow(self):
        """Dynamic generation of window menu"""
        self.menuVentana.clear()

        # Add subwindow options
        actionLeft = createAction(
            self.tr("&Previous"),
            icon=os.path.join("button", "arrow-left.png"),
            slot=self.currentMdi.activatePreviousSubWindow,
            shortcut=QtGui.QKeySequence.StandardKey.PreviousChild, parent=self)
        self.menuVentana.addAction(actionLeft)
        actionRight = createAction(
            self.tr("&Next"),
            icon=os.path.join("button", "arrow-right.png"),
            slot=self.currentMdi.activateNextSubWindow,
            shortcut=QtGui.QKeySequence.StandardKey.NextChild, parent=self)
        self.menuVentana.addAction(actionRight)
        self.menuVentana.addSeparator()

        actionTile = createAction(
            self.tr("&Tile"),
            icon=os.path.join("button", "tile.png"),
            slot=self.currentMdi.tileSubWindows, parent=self)
        self.menuVentana.addAction(actionTile)
        actionCascade = createAction(
            self.tr("&Cascade"),
            icon=os.path.join("button", "cascade.png"),
            slot=self.currentMdi.cascadeSubWindows, parent=self)
        self.menuVentana.addAction(actionCascade)
        actionRestore = createAction(
            self.tr("&Restore All"),
            slot=self.windowRestoreAll, parent=self)
        self.menuVentana.addAction(actionRestore)
        actionIconize = createAction(
            self.tr("&Iconize All"),
            slot=self.windowMinimizeAll, parent=self)
        self.menuVentana.addAction(actionIconize)
        self.menuVentana.addSeparator()

        # Add subwindow list
        active = self.currentMdi.activeSubWindow()
        for i, window in enumerate(self.currentMdi.subWindowList()):
            if window is active:
                iconPath = IMAGE_PATH + "button/ok.png"
            else:
                iconPath = ""
            self.menuVentana.addAction(
                QtGui.QIcon(iconPath),
                f"&{i+1} {window.windowTitle()}",
                partial(self.windowSelect, i))
        self.menuVentana.addSeparator()

        actionClose = createAction(
            self.tr("&Close window"),
            icon=os.path.join("button", "fileClose.png"),
            slot=self.windowClose, parent=self)
        mdiwindow = self.currentMdi.activeSubWindow()
        if mdiwindow.windowTitle() in (
                self.tr("Overview Window"), self.tr("Flow Diagram")):
            actionClose.setEnabled(False)
        self.menuVentana.addAction(actionClose)

    def windowClose(self):
        """Close active subwindow"""
        self.currentMdi.closeActiveSubWindow()
        self.dirty[self.idTab] = True
        self.saveControl()

    def windowSelect(self, index):
        """Show the selected subwindow"""
        window = self.currentMdi.subWindowList()[index]
        self.currentMdi.setActiveSubWindow(window)

    def windowRestoreAll(self):
        """Restore all subwindows to last window size and position"""
        for window in self.currentMdi.subWindowList():
            window.showNormal()

    def windowMinimizeAll(self):
        """Minimize all subwindows"""
        for window in self.currentMdi.subWindowList():
            window.showMinimized()

    def contextListMenu(self, event):
        """Generation of context menu for list toolbox"""
        contextMenu = QtWidgets.QMenu()
        self.currentScene.addActions(contextMenu, event)
        contextMenu.exec(event)

    def aboutToShow_MenuRecentFiles(self):
        """Dynamic generation of recent files menu"""
        self.menuRecentFiles.clear()
        recentFiles = []
        for fname in self.recentFiles:
            if fname not in self.filename and QtCore.QFile.exists(fname):
                recentFiles.append(fname)
        if recentFiles:
            self.menuRecentFiles.addSeparator()
            for i, fname in enumerate(recentFiles):
                action = QtGui.QAction(f"&{i+1} {fname}", self)
                action.setData(QtCore.QVariant(fname))
                action.triggered.connect(self.loadFile)
                self.menuRecentFiles.addAction(action)
        self.menuRecentFiles.addSeparator()
        self.menuRecentFiles.addAction(
            QtGui.QIcon(os.path.join(IMAGE_PATH, "button", "clear.png")),
            self.tr("Clear"),
            self.clearRecentFiles)

        # Disable clear option if menu is empty
        if not recentFiles:
            self.menuRecentFiles.actions()[-1].setEnabled(False)

    # File Manipulation
    def clearRecentFiles(self):
        """Clear recent open files list"""
        self.recentFiles = []
        self.menuRecentFiles.setEnabled(False)

    def addRecentFile(self, fname):
        """Populate recent file menu"""
        if fname and fname not in self.recentFiles:
            self.recentFiles.insert(0, fname)
            if len(self.recentFiles) > 9:
                self.recentFiles = self.recentFiles[:9]
        self.menuRecentFiles.setEnabled(len(self.recentFiles))

    def fileNew(self):
        """Actions to be done when create a new project"""
        UI_pychemqt.idNew += 1
        self.dirty.append(True)
        self.filename.append("")
        conf = ConfigParser()
        conf.add_section("PFD")
        conf.set("PFD", "x", Preferences.get("PFD", "x"))
        conf.set("PFD", "y", Preferences.get("PFD", "y"))
        self.config.append(conf)
        mdiArea = QtWidgets.QMdiArea()
        self.centralWidget().addTab(
            mdiArea,
            self.tr("New Project") + f" {UI_pychemqt.idNew}")
        self.centralWidget().setCurrentIndex(self.centralWidget().count()-1)
        self.loadPFD(mdiArea)
        self.updateStatus(self.tr("New project created"))
        self.activeControl(True)
        self.wizard()

    def fileSave(self, indice=None):
        """Save project to file"""
        if indice is None:
            indice = self.idTab

        if not self.filename[indice]:
            self.fileSaveAs()
        else:
            with open(self.filename[indice], "w") as file:
                data = {}
                self.getScene(indice).project.writeToJSON(data)

                PFD = {}
                win = self.centralWidget().currentWidget().subWindowList()[0]

                PFD["minimized"] = win.isMinimized()
                PFD["maximized"] = win.isMaximized()
                if win.isMinimized():
                    win.showNormal()
                    PFD["x"] = win.pos().x()
                    PFD["y"] = win.pos().y()
                    PFD["height"] = win.size().height()
                    PFD["width"] = win.size().width()
                    win.showMinimized()
                elif win.isMaximized():
                    win.showNormal()
                    PFD["x"] = win.pos().x()
                    PFD["y"] = win.pos().y()
                    PFD["height"] = win.size().height()
                    PFD["width"] = win.size().width()
                    win.showMaximized()
                else:
                    PFD["x"] = win.pos().x()
                    PFD["y"] = win.pos().y()
                    PFD["height"] = win.size().height()
                    PFD["width"] = win.size().width()

                self.currentScene.writeToJSON(PFD)
                data["PFD"] = PFD

                other = {}
                ventanas = self.centralWidget().currentWidget().subWindowList()
                for ind, win in enumerate(ventanas[1:]):
                    ventana = {}
                    ventana["class"] = win.widget().__class__.__name__
                    ventana["minimized"] = win.isMinimized()
                    ventana["maximized"] = win.isMaximized()
                    if win.isMinimized():
                        win.showNormal()
                        ventana["x"] = win.pos().x()
                        ventana["y"] = win.pos().y()
                        ventana["height"] = win.size().height()
                        ventana["width"] = win.size().width()
                        win.showMinimized()
                    elif win.isMaximized():
                        win.showNormal()
                        ventana["x"] = win.pos().x()
                        ventana["y"] = win.pos().y()
                        ventana["height"] = win.size().height()
                        ventana["width"] = win.size().width()
                        win.showMaximized()
                    else:
                        ventana["x"] = win.pos().x()
                        ventana["y"] = win.pos().y()
                        ventana["height"] = win.size().height()
                        ventana["width"] = win.size().width()

                    widget = {}
                    win.widget().writeToJSON(widget)
                    ventana["window"] = widget
                    other[ind] = ventana

                    # Add dependences from other windows
                    if widget.get("external_dependences", None):
                        data["external_dependences"].add(
                            widget["external_dependences"])

                data["other"] = other

                # python set are not serializable so convert to list
                data["external_dependences"] = list(
                    data["external_dependences"])

                json.dump(data, file, indent=4)

            self.dirty[self.idTab] = False
            self.updateStatus(self.tr("Saved as")+" "+self.filename[indice])
            self.dirty[indice] = False
            self.saveControl()

    def fileSaveAs(self, indice=None):
        """Save current project with a new name"""
        if indice is None:
            indice = self.idTab
        folder = self.filename[indice] if self.filename[indice] else "."
        fname = QtWidgets.QFileDialog.getSaveFileName(
            self,
            self.tr("Save project"),
            folder, "pychemqt project file (*.pcq)")

        if fname[0]:
            name = fname[0]
            if name.split(".")[-1] != "pcq":
                name += ".pcq"
            self.addRecentFile(name)
            self.filename[indice] = name
            self.fileSave(indice)
            self.centralWidget().setTabText(
                indice, os.path.splitext(os.path.basename(name))[0])

    def fileSaveAll(self):
        """Save all open project with changes"""
        for tab in range(self.centralWidget().count()):
            if self.dirty[tab]:
                self.fileSave(tab)

    def fileOpen(self, fname=None):
        """Open project"""
        if not fname:
            if self.filename:
                folder = os.path.dirname(str(self.filename[-1]))
            else:
                folder = "."
            fname = QtWidgets.QFileDialog.getOpenFileName(
                self,
                self.tr("Open project"),
                folder, self.tr("pychemqt project file") + " (*.pcq)")[0]

        if fname:
            try:
                self.loadFile(fname)
            except ImportError as e:
                QtWidgets.QMessageBox.warning(self, self.tr("Error"), e.msg)
            except Exception as error:
                QtWidgets.QMessageBox.critical(
                    self,
                    self.tr("Error"),
                    self.tr("Failed to load file") + "\n" + fname)
                raise error
            else:
                self.activeControl(True)

    def loadFile(self, fname=None):
        """Actions to be done when load a project file"""
        # Retrive fname if sender is actions in recent files list
        if not fname:
            action = self.sender()
            if isinstance(action, QtGui.QAction):
                fname = str(action.data())
            else:
                return

        with open(fname, "r") as file:
            data = json.load(file)

        # Check availability of optional dependences necessary for the file
        if "external_dependences" in data:
            available = True
            for dep in data["external_dependences"]:
                if os.environ[dep] != "True":
                    available = False
                    break

            if not available:
                msg = self.tr("Failed to load")
                msg += " " + fname + os.linesep
                msg += self.tr("This project require")
                msg += ": " + ", ".join(data["external_dependences"])
                raise ImportError(msg)

        self.dirty.append(False)
        self.filename.append(fname)
        self.addRecentFile(fname)
        self.clearWindow()

        project = Project()
        project.readFromJSON(data)
        self.config.append(project.config)

        mdiArea = QtWidgets.QMdiArea()
        self.loadPFD(mdiArea)

        x = data["PFD"]["x"]
        y = data["PFD"]["y"]
        width = data["PFD"]["width"]
        height = data["PFD"]["height"]
        mdiArea.subWindowList()[0].move(QtCore.QPoint(x, y))
        mdiArea.subWindowList()[0].resize(QtCore.QSize(width, height))

        if data["PFD"]["maximized"]:
            mdiArea.subWindowList()[-1].showMaximized()
        elif data["PFD"]["minimized"]:
            mdiArea.subWindowList()[-1].showMinimized()

        mdiArea.subWindowList()[0].widget().scene().readFromJSON(data)
        self.list.updateList(
            mdiArea.subWindowList()[0].widget().scene().objects)

        for ventana in data["other"].values():
            indice = other_window_names.index(ventana["class"])
            widget = other_window[indice]
            instance = widget.readFromJSON(ventana["window"], self)
            mdiArea.addSubWindow(instance)
            x = ventana["x"]
            y = ventana["y"]
            h = ventana["height"]
            w = ventana["width"]
            mdiArea.subWindowList()[-1].move(QtCore.QPoint(x, y))
            mdiArea.subWindowList()[-1].resize(QtCore.QSize(w, h))

            # FIXME: This line raise error in matplotlib as plot has no
            # window yet, disabled minimized and maximized state at start
            # meanwhile
            # if ventana["maximized"]:
                # mdiArea.subWindowList()[-1].showMaximized()
            # elif ventana["minimized"]:
                # mdiArea.subWindowList()[-1].showMinimized()

        self.centralWidget().addTab(
            mdiArea, os.path.splitext(os.path.basename(str(fname)))[0])
        self.centralWidget().setCurrentIndex(self.centralWidget().count()-1)
        mdiArea.subWindowActivated.connect(self.changeWindow)

        self.currentScene.project = project
        self.updateStatus(self.tr("Load") + " " + fname, True)

        self.activeControl(True)
        self.changeWindow(mdiArea.subWindowList()[0])
        self.changeStatusThermo(self.config[self.idTab])

    def loadPFD(self, mdiarea):
        """Actions to be done to load a PFD"""
        PFD = flujo.GraphicsView(True, self)
        PFD.setWindowTitle(self.tr("Flow Diagram"))
        scene = flujo.GraphicsScene(self)
        scene.selectionChanged.connect(self.selectionChanged)
        x = self.config[-1].getint("PFD", "x")
        y = self.config[-1].getint("PFD", "y")
        scene.setSceneRect(0, 0, x, y)
        PFD.setScene(scene)
        scene.defineShortcut()
        PFD.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(
            os.path.join(IMAGE_PATH, "button", "PFD.png"))))
        PFD.zoomChanged.connect(self.zoomValue.setValue)
        mdiarea.addSubWindow(PFD)
        mdiarea.subWindowList()[-1].setWindowFlags(
            QtCore.Qt.WindowType.CustomizeWindowHint
            | QtCore.Qt.WindowType.WindowTitleHint
            | QtCore.Qt.WindowType.WindowMinMaxButtonsHint)
        PFD.show()

    def changeWindow(self, window):
        """Update status info when change subwindow"""
        # Avoid remove message for updated widget when calculation is going on
        if self.statusBar().currentMessage():
            return

        try:
            wdgs = window.widget().statusWidget
        except AttributeError:
            return
        self.clearWindow()

        self.wdg = wdgs
        for wdg in wdgs:
            self.statusBar().addWidget(wdg)
            wdg.show()

    def clearWindow(self):
        """Actions to be done when close a window in a tab mdi"""
        for wdg in self.wdg:
            self.statusBar().removeWidget(wdg)
        self.statusBar().reformat()

    def fileClose(self, index):
        """Actions to be done when close a project in mainwindow"""
        if self.okToContinue(index):
            self.updateStatus(self.tr( "Closed") + " " + self.currentFilename)
            self.centralWidget().removeTab(index)
            del self.dirty[index]
            del self.config[index]
            del self.filename[index]
            if self.centralWidget().count():
                self.activeControl(True)
            else:
                self.activeControl(False)
                self.list.clear()
                flujo.StreamItem.id = 0
                flujo.EquipmentItem.id = 0
            self.clearWindow()

    # Configuration
    def wizard(self):
        """Show a wizard dialog for project configuration"""
        dialog = wizard.Wizard(self.config[self.idTab])
        if dialog.exec():
            self.updateConfig(dialog.value)
            self.updateStatus(self.tr("Project configuration"), True)
        else:
            self.updateConfig(wizard.Wizard.default())
            self.updateStatus(self.tr("Project configuration"), False)

    def updateConfig(self, conf):
        """Do changes in mainwindow at any config change"""
        self.config[self.idTab] = conf
        self.currentScene.project.setConfig(conf)
        self.dirty[self.idTab] = True
        x = conf.getint("PFD", "x")
        y = conf.getint("PFD", "y")
        self.currentScene.setSceneRect(0, 0, x, y)
        self.changeStatusThermo(conf)
        setMainWindowConfig(conf)

        # TODO: Delete this when its not necessary to run library isolated
        with open(conf_dir+"pychemqtrc_temporal", "w") as conf_file:
            conf.write(conf_file)

    def updateStatus(self, text, success=True):
        """Función que añade entradas al cuadro de status
        text: texto a mostrar
        success: boolean que indica si va todo bien"""
        if success:
            txt = self.tr("Success")
            color = "#00aa00"
        else:
            txt = self.tr("Failure")
            color = "#ff0000"
        self.status.append(
            '<b>' + time.strftime("%H:%M:%S", time.localtime()) + '</b> - '
            + text + f' [<font color="{color}">{txt}</font>]')
        QtWidgets.QApplication.processEvents()

    def changeStatusThermo(self, conf):
        """Update status thermo for pfd subwindow"""
        self.currentView.changeStatusThermo(conf)

    def dialogConfig(self, UIconfig):
        """Generic procedure for launch a config dialog"""
        Dialog = UIconfig.Dialog(self.config[self.idTab])
        if Dialog.exec():
            conf = Dialog.value(self.config[self.idTab])
            self.updateConfig(conf)
            self.saveControl()

    def Preferences(self):
        """Show preferences dialog"""
        dialog = UI_Preferences.Preferences(Preferences)
        if dialog.exec():
            preferences = dialog.value()
            self.updatePreferences(preferences)
            self.updateStatus(self.tr("pychemqt configuration change"), True)
        else:
            self.updateStatus(self.tr("pychemqt configuration change"), False)

    def updatePreferences(self, preferences):
        """Update preferences in program and save to file"""

        # Save preferences to file
        with open(conf_dir+"pychemqtrc", "w") as pref_file:
            preferences.write(pref_file)

        # Load new preferences in main program preferences variable
        global Preferences
        Preferences = ConfigParser()
        Preferences.read(conf_dir+"pychemqtrc")
        config.Preferences = Preferences
        self.changePreferenceLive()

    def changePreferenceLive(self):
        """Upgrade UI when change preferences"""

        # Upgrade tray icon status
        if Preferences.getboolean("General", 'Tray'):
            self.systemtray.show()
        else:
            self.systemtray.hide()

        # Change PFD brush
        if self.currentMdi:
            for subwindow in self.currentMdi.subWindowList():
                if isinstance(subwindow.widget(), flujo.GraphicsView):
                    brushColor = Preferences.get("PFD", "brushColor")
                    stl = BrushCombo.BRUSH[Preferences.getint("PFD", "brush")]
                    subwindow.widget().setBackgroundBrush(
                        QtGui.QBrush(QtGui.QColor(brushColor), stl))

        # Check availability of system calculator
        calculator = str(Preferences.get("Applications", 'Calculator'))
        self.calculatorAction.setEnabled(bool(calculator))

    def activeControl(self, boolean):
        """Enable/Disable actions if any project is loaded"""
        self.fileSaveAsAction.setEnabled(boolean)
        self.fileSaveAllAction.setEnabled(boolean)
        self.fileCloseAction.setEnabled(boolean)
        self.actionWizard.setEnabled(boolean)
        self.actionComponentList.setEnabled(boolean)
        self.actionThermo.setEnabled(boolean)
        self.actionTransporte.setEnabled(boolean)
        self.actionUnidades.setEnabled(boolean)
        self.actionZoom.setEnabled(boolean)
        self.actionZoomIn.setEnabled(boolean)
        self.actionZoomOut.setEnabled(boolean)
        self.zoomValue.setEnabled(boolean)
        self.menuPFD.setEnabled(boolean)
        # self.menuPlot.setEnabled(boolean)
        self.menuVentana.setEnabled(boolean)
        self.menuMEoS.setEnabled(boolean)
        self.toolboxItem.setEnabled(boolean)
        self.toolboxPalette.setEnabled(boolean)
        if boolean:
            self.saveControl()
        else:
            self.fileSaveAction.setEnabled(False)

    def saveControl(self):
        """Enable/Disable save actions is current tab is modified"""
        self.fileSaveAction.setEnabled(self.dirty[self.idTab])
        self.tabModified(self.idTab)

    def tabModified(self, index):
        """Set icon as modified for tab with index specified"""
        if self.dirty[index]:
            icon = QtGui.QIcon(IMAGE_PATH + "button/editor.png")
        else:
            icon = QtGui.QIcon()
        self.centralWidget().setTabIcon(index, icon)

    def currentTabChanged(self, index):
        """Actions to be done when change current tab in main window"""
        if index == -1:
            self.list.clear()
            self.activeControl(False)
        elif self.filename[-1]:
            scene = self.currentScene
            self.list.updateList(scene.objects)
            self.centralWidget().currentWidget().subWindowList()[0].widget()
            self.list.updateList(self.currentScene.objects)
            self.changeStatusThermo(self.config[self.idTab])
        else:
            self.list.clear()

    # Tools
    def calculator(self):
        """Show external calculator application"""
        command = str(Preferences.get("Applications", 'Calculator'))
        os.system(command)

    def terminal(self):
        """Show instrospection console"""
        terminal.XTerm(Preferences)

    def meos(self):
        """Launch main dialog of UI_Tables tool"""
        dialog = UI_Tables.Dialog(self.currentConfig, self)
        dialog.exec()

    def launch(self, dlg, *args):
        """Generic method to launch external application"""
        # Remove checked args send by a qaction
        if dlg != newComponent.Ui_Contribution:
            dialog = dlg(*args[1:])
        else:
            dialog = dlg(*args)

        # Add text to statusbar
        txt = dialog.__dict__.get("__TEXTSTATUS__", None)
        if txt:
            self.updateStatus(txt)

        dialog.exec()

    # plots
    # def plot(self, indice, x=None, y=None):
        # grafico = plots._all[indice]()
        # if grafico.exec():
            # self.currentMdi.addSubWindow(grafico.plot)
            # grafico.plot.show()

        # indices, nombres, M=getComponents(config=self.config[self.idTab])
        # dialog=grafico(indices, nombres, x, y)
        # self.currentMdi.addSubWindow(dialog)
        # dialog.show()

    # Help
    def help(self):
        """Open pychemqt documentation"""
        # First search for a local documentation build
        # By default in docs/_build/html/
        path = os.path.join(os.environ["pychemqt"], "docs") + os.sep
        if os.path.isdir(path):
            indexpath = os.path.join(path, "_build", "html", "index.html")
            url = QtCore.QUrl(indexpath)
        else:
            url = QtCore.QUrl("http://pychemqt.readthedocs.io/")
        QtGui.QDesktopServices.openUrl(url)

    def documentation(self):
        """Open tools/doi app with scientific references"""
        dialog = doi.ShowReference()
        dialog.exec()

    def log(self):
        """Open log file"""
        command = Preferences.get("Applications", 'TextViewer')
        logfile = os.path.join(conf_dir, "pychemqt.log")
        subprocess.Popen([command, logfile])

    def acerca(self):
        """Show basic information about pychemqt including version of
        library used"""
        txt = self.tr(
            "Software for simulate units operations in Chemical Engineering")
        QtWidgets.QMessageBox.about(
            self,
            self.tr("About pychemqt"),
            f"<b>pychemqt</b> v {__version__}"
            f"<p>Copyright &copy; {year} Juan José Gómez Romera (jjgomera)<br>"
            "Licenced with GPL.v3"
            f"<p>{txt}"
            f"<p>Python {platform.python_version()}"
            f" - Qt {QtCore.QT_VERSION_STR} - PyQt {QtCore.PYQT_VERSION_STR}"
            f" on {platform.system()}")

    def acercaQt(self):
        """Show info about qt library"""
        QtWidgets.QMessageBox.aboutQt(self, self.tr("About Qt"))

    def savePFDImage(self):
        """Save PFD screenshot in a file"""
        if self.filename[self.idTab]:
            folder = os.path.dirname(str(self.filename[self.idTab]))
        else:
            folder = "."
        fname = QtWidgets.QFileDialog.getSaveFileName(
            None,
            self.tr("Save PFD as image"),
            folder, "Portable Network Graphics (*.png)")[0]
        if fname:
            if fname.split(".")[-1] != "png":
                fname += ".png"
            rect = self.currentScene.sceneRect()
            img = QtGui.QImage(
                int(rect.width()), int(rect.height()),
                QtGui.QImage.Format.Format_ARGB32_Premultiplied)
            p = QtGui.QPainter(img)
            self.currentScene.render(p)
            p.end()
            img.save(fname)

    def zoom(self, value=None):
        """Change zoom value of PFD

        Parameters
        ----------
        Value : float
            Zoom value or any of this str
            + : Zoom in view
            - : Zoom out view
            Dialog : Show dialog to set custom zoom value
        """
        # Zoom in and zoon out buttom step
        step = 5

        if value == "+":
            self.zoomValue.setValue(self.zoomValue.value() + step)
        elif value == "-":
            self.zoomValue.setValue(self.zoomValue.value() - step)
        elif value == "Dialog":
            value, check = QtWidgets.QInputDialog.getInt(
                self,
                self.tr("Zoom"),
                self.tr("Zoom factor:"),
                self.zoomValue.value())
            if check:
                self.zoomValue.setValue(value)
        else:
            self.currentView.zoom(value)

    def overview(self):
        """Show a overview window of PFD or hide if exist"""

        # Check if overview window is loaded
        for window in self.currentMdi.subWindowList():
            if window.windowTitle() == self.tr("Overview Window"):
                window.close()
                return

        PFD = flujo.GraphicsView(False, self)
        PFD.setWindowTitle(self.tr("Overview Window"))
        PFD.setScene(self.currentScene)

        # Use a custom 20% zoom out
        PFD.zoom(20)

        self.currentMdi.addSubWindow(PFD)

        # Disable close button in subwindow titlebar
        self.currentMdi.subWindowList()[-1].setWindowFlags(
            QtCore.Qt.WindowType.CustomizeWindowHint
            | QtCore.Qt.WindowType.WindowTitleHint
            | QtCore.Qt.WindowType.WindowMinMaxButtonsHint)
        PFD.show()

    def selectionChanged(self):
        """Selection changed from both PFD and item treewidget to let both
        stay syncronized with same selected items"""
        sender = self.sender()
        if isinstance(sender, QtWidgets.QTreeWidget):
            self.currentScene.blockSignals(True)
            self.currentScene.clearSelection()
            for element in sender.selectedItems():
                if element.parent() == self.list.Equipment:
                    obj = self.currentScene.getObject(
                        "e", int(element.text(0).split(" ")[0]))
                    obj.setSelected(True)
                if element.parent() == self.list.Stream:
                    obj = self.currentScene.getObject(
                        "s", int(element.text(0)))
                    obj.setSelected(True)
            self.currentScene.blockSignals(False)
        elif isinstance(sender, flujo.GraphicsScene):
            self.list.blockSignals(True)
            self.list.clearSelection()
            for element in sender.selectedItems():
                if isinstance(element, flujo.StreamItem):
                    element = self.list.Stream.child(element.id-1)
                    if element:
                        element.setSelected(True)
                elif isinstance(element, flujo.EquipmentItem) and \
                        element.tipo == "e":
                    element = self.list.Equipment.child(element.id-1)
                    if element:
                        element.setSelected(True)
            self.list.blockSignals(False)

    def addText(self):
        """Add text item to PFD of current project"""
        dialog = flujo.TextItemDlg()
        if dialog.exec():
            txtItem = flujo.TextItem(dialog.editor.texto)
            self.currentScene.waitClick(1, "txt", txtItem)

    def addItem(self, tipo, check=True):
        """Add item to PFD of current project
        tipo: type of item
          * square: square item, defining 2 points of diagonal
          * ellipse: ellipse item, enclosed in a square defined with 2 points
          * Input stream to PFD, defined only by one point
          * Output stream to PFD, defined only by one point
          * Strem: Stream to join several equipment"""
        if tipo == "square":
            obj = flujo.RectItem()
            num = 2
        elif tipo == "ellipse":
            obj = flujo.EllipseItem()
            num = 2
        elif tipo == "in":
            obj = flujo.EquipmentItem("in", None)
            num = 1
        elif tipo == "out":
            obj = flujo.EquipmentItem("out", None)
            num = 1
        elif tipo == "stream":
            if check:
                obj = flujo.StreamItem()
                num = 2
            else:
                self.currentScene.clickCollector.quit()
                return
        else:
            return

        self.currentScene.waitClick(num, tipo, obj)

    def addEquipment(self, equipo):
        """Add equipment to PFD of current project"""
        equip = equipment.UI_equipments.index(equipo)
        obj = flujo.EquipmentItem(equipo.__name__.split("_")[-1], equip)

        # Wait click in PFD to set position
        self.currentScene.waitClick(1, "equip", obj)
