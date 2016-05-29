#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>."""


from configparser import ConfigParser
from functools import partial
import json
import os
import platform
import sys
import time

from PyQt5 import QtCore, QtGui, QtWidgets

from UI import newComponent, flujo, wizard, plots, viewComponents
import plots as charts
from UI.widgets import createAction, ClickableLabel, TreeEquipment, FlowLayout
from lib.config import (conf_dir, QTSETTING_FILE, setMainWindowConfig,
                        IMAGE_PATH)
from lib.project import Project
from lib.EoS import K, H
from equipment import *  # noqa
from tools import (UI_confComponents, UI_Preferences, UI_confTransport,
                   UI_confThermo, UI_confUnits, UI_confResolution, UI_databank,
                   UI_Tables, UI_unitConverter, UI_steamTables,
                   UI_psychrometry, costIndex, doi, dependences)
from UI.conversor_unidades import moneda

__version__ = "0.1.0"

other_window = (plots.Binary_distillation, UI_Tables.TablaMEoS,
                UI_Tables.PlotMEoS)
other_window_names = [cl.__name__ for cl in other_window]


class TabWidget(QtWidgets.QTabWidget):
    """Customize an empty tabwidget to show a welcome and information message"""
    def paintEvent(self, event):
        if self.count():
            QtWidgets.QTabWidget.paintEvent(self, event)
        else:
            painter = QtGui.QPainter(self)
            rect = event.rect()
            image = QtGui.QImage("images/pychemqt.png")
            rectImage = QtCore.QRect(25, rect.center().y()-50, 100, 100)
            painter.drawImage(rectImage, image)
            txt = QtWidgets.QApplication.translate(
                "pychemqt", """Welcome to pychemqt,
a software for simulating Chemical Engineering units operations,

Copyright © 2012 Juan José Gómez Romera (jjgomera)
Licenced with GPL.v3
This software is distributed in the hope that it will be useful,
but without any warranty, it is provided "as is" without warranty of any kind

You can start by creating a new project or opening an existing project.""",
                None)
            rect.setLeft(150)
            painter.drawText(rect, QtCore.Qt.AlignVCenter, txt)


class UI_pychemqt(QtWidgets.QMainWindow):
    """Main window definition"""
    idNew = 0
    config = []
    dirty = []
    filename = []
    pfd = []

    def __init__(self, parent=None):
        super(UI_pychemqt, self).__init__(parent)
        self.setWindowTitle("pychemqt")
        self.Preferences = ConfigParser()
        self.Preferences.read(conf_dir+"pychemqtrc")
        self.centralwidget = TabWidget()
        self.centralwidget.currentChanged.connect(self.currentTabChanged)
        self.centralwidget.setTabsClosable(True)
        self.setCentralWidget(self.centralwidget)
        icon = IMAGE_PATH + "pychemqt.png"
        self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(icon)))

        # Acciones
        fileNewAction = createAction(
            QtWidgets.QApplication.translate("pychemqt", "&New"),
            slot=self.fileNew,
            shortcut=QtGui.QKeySequence.New,
            icon="button/fileNew",
            tip=QtWidgets.QApplication.translate(
                "pychemqt", "Start new project"),
            parent=self)
        fileOpenAction = createAction(
            QtWidgets.QApplication.translate("pychemqt", "&Open"),
            slot=self.fileOpen,
            shortcut=QtGui.QKeySequence.Open,
            icon="button/fileOpen",
            tip=QtWidgets.QApplication.translate("pychemqt", "Open project"),
            parent=self)
        self.fileSaveAction = createAction(
            QtWidgets.QApplication.translate("pychemqt", "&Save"),
            slot=self.fileSave,
            shortcut=QtGui.QKeySequence.Save,
            icon="button/fileSave",
            tip=QtWidgets.QApplication.translate("pychemqt", "Save project"),
            parent=self)
        self.fileSaveAsAction = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Save &as"),
            slot=self.fileSaveAs,
            shortcut=QtGui.QKeySequence.SaveAs,
            icon="button/fileSaveAs",
            tip=QtWidgets.QApplication.translate(
                "pychemqt", "Save project as"),
            parent=self)
        self.fileSaveAllAction = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Save A&ll"),
            slot=self.fileSaveAll,
            icon="button/fileSaveAll",
            tip=QtWidgets.QApplication.translate(
                "pychemqt", "Save all open project"),
            parent=self)
        self.fileCloseAction = createAction(
            QtWidgets.QApplication.translate("pychemqt", "&Close"),
            slot=self.fileClose,
            shortcut=QtGui.QKeySequence.Close,
            icon="button/fileClose",
            tip=QtWidgets.QApplication.translate("pychemqt", "Close project"),
            parent=self)
        ExitAction = createAction(
            QtWidgets.QApplication.translate("pychemqt", "&Exit"),
            slot=self.closeEvent,
            shortcut=QtGui.QKeySequence.Quit,
            icon="button/exit",
            tip=QtWidgets.QApplication.translate("pychemqt", "Quit pychemqt"),
            parent=self)

        self.actionWizard = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Wizard"),
            icon="button/wizard",
            slot=self.wizard,
            tip=QtWidgets.QApplication.translate(
                "pychemqt", "Launch configuration wizard"),
            parent=self)
        self.actionComponentList = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Components list"),
            slot=partial(self.dialogConfig, UI_confComponents),
            tip=QtWidgets.QApplication.translate(
                "pychemqt", "Defining componente list dialog"),
            parent=self)
        self.actionThermo = createAction(
            QtWidgets.QApplication.translate(
                "pychemqt", "Thermodynamic properties"),
            slot=partial(self.dialogConfig, UI_confThermo),
            tip=QtWidgets.QApplication.translate(
                "pychemqt", "Defining thermodynamic properties methods"),
            parent=self)
        self.actionTransporte = createAction(
            QtWidgets.QApplication.translate(
                "pychemqt", "Transport properties"),
            slot=partial(self.dialogConfig, UI_confTransport),
            tip=QtWidgets.QApplication.translate(
                "pychemqt", "Defining transport properties methods"),
            parent=self)
        self.actionUnidades = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Units"),
            slot=partial(self.dialogConfig, UI_confUnits),
            tip=QtWidgets.QApplication.translate(
                "pychemqt", "Defining preferred units"),
            parent=self)
        self.actioncostIndex = createAction(
            QtWidgets.QApplication.translate("pychemqt", "&Cost Index"),
            slot=self.costos,
            tip=QtWidgets.QApplication.translate(
                "pychemqt", "Defining cost index"),
            parent=self)
        self.actionPreferencias = createAction(
            QtWidgets.QApplication.translate("pychemqt", "&Preferences"),
            slot=self.Preferencias,
            icon="button/configure",
            shortcut=QtGui.QKeySequence.Preferences,
            tip=QtWidgets.QApplication.translate(
                "pychemqt", "Defining general preferences"),
            parent=self)

        self.actionZoomIn = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Zoom in"),
            slot=partial(self.zoom, "+"),
            icon="button/zoomIn",
            shortcut=QtGui.QKeySequence.ZoomIn, parent=self)
        self.actionZoom = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Zoom"),
            slot=partial(self.zoom, "Dialog"),
            icon="button/zoomIn",
            parent=self)
        self.actionZoomOut = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Zoom out"),
            slot=partial(self.zoom, "-"),
            icon="button/zoomOut",
            shortcut=QtGui.QKeySequence.ZoomOut,
            parent=self)
        actionOverviewWindow = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Overview window"),
            slot=self.overview, parent=self)
        actionVerStatus = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Status"),
            shortcut="Ctrl+Alt+S",
            tip=QtWidgets.QApplication.translate(
                "pychemqt", "Show/Hide status toolbar"),
            checkable=True, parent=self)
        self.actionVerToolbar = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Palette"),
            shortcut="Ctrl+Alt+T",
            icon="button/palette",
            tip=QtWidgets.QApplication.translate(
                "pychemqt", "Show/Hide equipment palette"),
            checkable=True, parent=self)
        self.actionVerItem = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Item"),
            icon="button/list",
            tip=QtWidgets.QApplication.translate(
                "pychemqt", "Show/Hide item list"),
            checkable=True, parent=self)

        calculatorAction = createAction(
            QtWidgets.QApplication.translate("pychemqt", "&Calculator"),
            slot=self.calculator,
            icon="button/calculator",
            shortcut="F2",
            tip=QtWidgets.QApplication.translate(
                "pychemqt", "Open system calculator"),
            parent=self)
        terminalAction = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Python Shell"),
            shortcut="F3",
            slot=self.terminal,
            icon="button/terminal",
            tip=QtWidgets.QApplication.translate(
                "pychemqt", "Open system terminal"),
            parent=self)
        if sys.platform == "win32":
            terminalAction.setEnabled(False)

        conversorUnidadesAction = createAction(
            QtWidgets.QApplication.translate("pychemqt", "&Units converter"),
            slot=self.conversor_unidades,
            shortcut="F4",
            icon="button/unitConverter",
            tip=QtWidgets.QApplication.translate(
                "pychemqt", "Open Units converter"),
            parent=self)
        currencyAction = createAction(
            QtWidgets.QApplication.translate(
                "pychemqt", "&Currency converter"),
            slot=self.conversor_moneda,
            shortcut="F5",
            icon="button/currency",
            tip=QtWidgets.QApplication.translate(
                "pychemqt", "Open Currency converter"),
            parent=self)
        TablaPeriodicaAction = createAction(
            QtWidgets.QApplication.translate("pychemqt", "&Periodic Table"),
            slot=self.tablaPeriodica,
            shortcut="F6",
            icon="button/periodicTable",
            tip=QtWidgets.QApplication.translate(
                "pychemqt", "Show a basic Mendeleiev periodic table"),
            parent=self)
        steamTablesAction = createAction(
            QtWidgets.QApplication.translate("pychemqt", "&Steam Tables"),
            slot=self.tablasVapor,
            shortcut="F7",
            icon="button/steamTables",
            tip=QtWidgets.QApplication.translate(
                "pychemqt",
                "Open a water-steam table and graphic application"),
            parent=self)
        psychrometricChartAction = createAction(
            QtWidgets.QApplication.translate(
                "pychemqt", "&Psicrometric Chart"),
            slot=self.diagramaPsicrometrico,
            shortcut="F9",
            icon="button/psychrometric",
            tip=QtWidgets.QApplication.translate(
                "pychemqt", "Open a humid-air application"),
            parent=self)
        externalProgramAction = createAction(
            QtWidgets.QApplication.translate("pychemqt", "External Programs"),
            slot=self.externalPrograms,
            icon="button/showPrograms",
            tip=QtWidgets.QApplication.translate(
                "pychemqt", "Show External Programs Status"),
            parent=self)
        saveAsImage = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Save PFD as image"),
            slot=self.savePFDImage,
            icon="button/image",
            parent=self)
        actionAyuda = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Help"),
            slot=self.help,
            icon="button/help",
            parent=self)
        actionAcerca_de = createAction(
            QtWidgets.QApplication.translate("pychemqt", "About pychemqt"),
            slot=self.acerca,
            icon="button/helpAbout",
            parent=self)
        actionAcerca_deQt = createAction(
            QtWidgets.QApplication.translate("pychemqt", "About Qt"),
            slot=self.acercaQt,
            icon="button/AboutQt",
            parent=self)

        self.zoomValue = QtWidgets.QSpinBox()
        self.zoomValue.setSuffix("%")
        self.zoomValue.setRange(5, 1000)
        self.zoomValue.setValue(100)
        self.zoomValue.setSingleStep(5)
        self.zoomValue.valueChanged.connect(self.zoom)

        # Toolbar
        self.BarraArchivo = QtWidgets.QToolBar(
            QtWidgets.QApplication.translate("pychemqt", "File"), self)
        self.BarraArchivo.setObjectName("BarraArchivo")
        self.BarraArchivo.setIconSize(QtCore.QSize(16, 16))
        self.BarraArchivo.addAction(fileNewAction)
        self.BarraArchivo.addAction(fileOpenAction)
        self.BarraArchivo.addAction(self.fileCloseAction)
        self.BarraArchivo.addAction(self.fileSaveAction)
        self.BarraArchivo.addAction(self.fileSaveAsAction)
        self.BarraArchivo.addAction(self.fileSaveAllAction)
        self.BarraArchivo.addAction(ExitAction)
        self.addToolBar(QtCore.Qt.TopToolBarArea, self.BarraArchivo)

        self.BarraVer = QtWidgets.QToolBar(
            QtWidgets.QApplication.translate("pychemqt", "View"), self)
        self.BarraVer.setObjectName("BarraVer")
        self.BarraVer.setIconSize(QtCore.QSize(16, 16))
        self.BarraVer.addAction(self.actionZoomOut)
        self.BarraVer.addWidget(self.zoomValue)
        self.BarraVer.addAction(self.actionZoomIn)
        self.BarraVer.addSeparator()
        self.BarraVer.addAction(self.actionVerToolbar)
        self.BarraVer.addAction(self.actionVerItem)
        self.addToolBar(QtCore.Qt.TopToolBarArea, self.BarraVer)

        self.BarraHerramientas = QtWidgets.QToolBar(
            QtWidgets.QApplication.translate("pychemqt", "Tools"), self)
        self.BarraHerramientas.setObjectName("BarraHerramientas")
        self.BarraHerramientas.setIconSize(QtCore.QSize(16, 16))
        self.BarraHerramientas.addAction(calculatorAction)
        self.BarraHerramientas.addAction(terminalAction)
        self.BarraHerramientas.addAction(conversorUnidadesAction)
        self.BarraHerramientas.addAction(currencyAction)
        self.BarraHerramientas.addAction(TablaPeriodicaAction)
        self.BarraHerramientas.addAction(steamTablesAction)
        self.BarraHerramientas.addAction(psychrometricChartAction)
        self.BarraHerramientas.addAction(externalProgramAction)
        self.BarraHerramientas.addSeparator()
        self.BarraHerramientas.addAction(actionAyuda)
        self.addToolBar(QtCore.Qt.TopToolBarArea, self.BarraHerramientas)

        # Paleta Toolbox
        self.toolboxPalette = QtWidgets.QDockWidget(
            QtWidgets.QApplication.translate("pychemqt", "Equipos"))
        self.toolboxPalette.setObjectName("toolbox")
        toolboxContenido = QtWidgets.QWidget()
        self.toolboxPalette.setWidget(toolboxContenido)
        self.addDockWidget(QtCore.Qt.DockWidgetArea(1), self.toolboxPalette)
        self.toolboxPalette.setFeatures(
            QtWidgets.QDockWidget.AllDockWidgetFeatures)
        self.toolboxPalette.setAllowedAreas(
            QtCore.Qt.LeftDockWidgetArea | QtCore.Qt.RightDockWidgetArea)

        self.toolboxPalette.visibilityChanged.connect(
            self.actionVerToolbar.setChecked)
        self.actionVerToolbar.triggered.connect(self.toolboxPalette.setVisible)
        layouttoolbox = QtWidgets.QVBoxLayout(toolboxContenido)
        layouttoolbox.setContentsMargins(5, 5, 5, 5)

        txt = QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Plot"))
        txt.setAlignment(QtCore.Qt.AlignCenter)
        layouttoolbox.addWidget(txt)
        l1 = FlowLayout(spacing=10)
        actionTexto, botonTexto = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Insert text"),
            icon="equipment/text",
            slot=self.addText,
            button=True, parent=toolboxContenido)
        l1.addWidget(botonTexto)
        actionCuadrado, botonCuadrado = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Draw square"),
            icon="equipment/square",
            slot=partial(self.addItem, "square"),
            button=True, parent=toolboxContenido)
        l1.addWidget(botonCuadrado)
        actionCircle, botonCircle = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Draw circle"),
            icon="equipment/circle",
            slot=partial(self.addItem, "ellipse"),
            button=True, parent=toolboxContenido)
        l1.addWidget(botonCircle)
        layouttoolbox.addItem(l1)
        layouttoolbox.addStretch(1)

        txt = QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Flux"))
        txt.setAlignment(QtCore.Qt.AlignCenter)
        layouttoolbox.addWidget(txt)
        l2 = FlowLayout(spacing=10)
        actionEntrada, botonEntrada = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Input"),
            icon="equipment/in",
            slot=partial(self.addItem, "in"),
            button=True, parent=toolboxContenido)
        l2.addWidget(botonEntrada)
        actionCorriente, botonCorriente = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Stream"),
            icon="equipment/stream",
            slot=partial(self.addItem, "stream"),
            button=True, checkable=True, parent=toolboxContenido)
        l2.addWidget(botonCorriente)
        actionSalida, botonSalida = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Output"),
            icon="equipment/out",
            slot=partial(self.addItem, "out"),
            button=True, parent=toolboxContenido)
        l2.addWidget(botonSalida)
        actionDivider, botonDivider = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Divider"),
            icon="equipment/divider",
            slot=partial(self.addEquipment, UI_divider),
            button=True, parent=toolboxContenido)
        l2.addWidget(botonDivider)
        actionValve, botonValve = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Valve"),
            icon="equipment/valve",
            slot=partial(self.addEquipment, UI_valve),
            button=True, parent=toolboxContenido)
        l2.addWidget(botonValve)
        actionMixer, botonMixer = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Mixer"),
            icon="equipment/mixer",
            slot=partial(self.addEquipment, UI_mixer),
            button=True, parent=toolboxContenido)
        l2.addWidget(botonMixer)
        actionCompresor, botonCompresor = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Compressor"),
            icon="equipment/compressor",
            slot=partial(self.addEquipment, UI_compressor),
            button=True, parent=toolboxContenido)
        l2.addWidget(botonCompresor)
        actionTurbine, botonTurbine = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Turbine"),
            icon="equipment/turbine",
            slot=partial(self.addEquipment, UI_turbine),
            button=True, parent=toolboxContenido)
        l2.addWidget(botonTurbine)
        actionPump, botonPump = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Pump"),
            icon="equipment/pump",
            slot=partial(self.addEquipment, UI_pump),
            button=True, parent=toolboxContenido)
        l2.addWidget(botonPump)
        actionPipe, botonPipe = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Pipe"),
            icon="equipment/pipe",
            slot=partial(self.addEquipment, UI_pipe),
            button=True, parent=toolboxContenido)
        l2.addWidget(botonPipe)
        layouttoolbox.addItem(l2)
        layouttoolbox.addStretch(1)

        txt = QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Basics"))
        txt.setAlignment(QtCore.Qt.AlignCenter)
        layouttoolbox.addWidget(txt)
        l3 = FlowLayout(spacing=10)
        actionTorreFUG, botonTorreFUG = createAction(
            QtWidgets.QApplication.translate(
                "pychemqt", "Distillation tower (method FUG)"),
            icon="equipment/columnFUG",
            slot=partial(self.addEquipment, UI_columnFUG),
            button=True, parent=toolboxContenido)
        l3.addWidget(botonTorreFUG)
        actionFlash, botonFlash = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Flash"),
            icon="equipment/flash",
            slot=partial(self.addEquipment, UI_flash),
            button=True, parent=toolboxContenido)
        l3.addWidget(botonFlash)
        actionTorre, botonTorre = createAction(
            QtWidgets.QApplication.translate(
                "pychemqt", "Distillation tower (exact method)"),
            icon="equipment/tower",
            slot=partial(self.addEquipment, UI_tower),
            button=True, parent=toolboxContenido)
        l3.addWidget(botonTorre)
        botonTorre.setEnabled(False)
        actionheatExchanger, botonheatExchanger = createAction(
            QtWidgets.QApplication.translate(
                "pychemqt", "Generic heat exchanger"),
            icon="equipment/heatExchanger",
            slot=partial(self.addEquipment, UI_heatExchanger),
            button=True, parent=toolboxContenido)
        l3.addWidget(botonheatExchanger)
        actionhairpin, botonhairpin = createAction(
            QtWidgets.QApplication.translate(
                "pychemqt", "Hairpin heat exchanger"),
            icon="equipment/hairpin",
            slot=partial(self.addEquipment, UI_hairpin),
            button=True, parent=toolboxContenido)
        l3.addWidget(botonhairpin)
        actionShellTube, botonShellTube = createAction(
            QtWidgets.QApplication.translate(
                "pychemqt", "Shell and tube heat exchanger"),
            icon="equipment/shellTube",
            slot=partial(self.addEquipment, UI_shellTube),
            button=True, parent=toolboxContenido)
        l3.addWidget(botonShellTube)
        actionFireHeater, botonFireHeater = createAction(
            QtWidgets.QApplication.translate(
                "pychemqt", "Fired Heater heat exchanger"),
            icon="equipment/fireHeater",
            slot=partial(self.addEquipment, UI_fireHeater),
            button=True, parent=toolboxContenido)
        l3.addWidget(botonFireHeater)
        actionReactor, botonReactor = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Reactor"),
            icon="equipment/reactor",
            slot=partial(self.addEquipment, UI_reactor),
            button=True, parent=toolboxContenido)
        l3.addWidget(botonReactor)
        layouttoolbox.addItem(l3)
        layouttoolbox.addStretch(1)

        txt = QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Solids"))
        txt.setAlignment(QtCore.Qt.AlignCenter)
        layouttoolbox.addWidget(txt)
        l4 = FlowLayout(spacing=10)
        actionBaghouse, botonBaghouse = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Baghouse"),
            icon="equipment/baghouse",
            slot=partial(self.addEquipment, UI_baghouse),
            button=True, parent=toolboxContenido)
        l4.addWidget(botonBaghouse)
        actionCentrifuge, botonCentrifuge = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Centrifuge"),
            icon="equipment/centrifuge",
            slot=partial(self.addEquipment, UI_centrifuge),
            button=True, parent=toolboxContenido)
        l4.addWidget(botonCentrifuge)
        botonCentrifuge.setEnabled(False)
        actionCiclon, botonCiclon = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Cyclone"),
            icon="equipment/ciclon",
            slot=partial(self.addEquipment, UI_ciclon),
            button=True, parent=toolboxContenido)
        l4.addWidget(botonCiclon)
        actionElectroPrecipitator, botonElectroPrecipitator = createAction(
            QtWidgets.QApplication.translate(
                "pychemqt", "Electric precipitator"),
            icon="equipment/electricPrecipitator",
            slot=partial(self.addEquipment, UI_electricPrecipitator),
            button=True, parent=toolboxContenido)
        l4.addWidget(botonElectroPrecipitator)
        actionGrinder, botonGrinder = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Grinder"),
            icon="equipment/grinder",
            slot=partial(self.addEquipment, UI_grinder),
            button=True, parent=toolboxContenido)
        l4.addWidget(botonGrinder)
        botonGrinder.setEnabled(False)
        actionDryer, botonDryer = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Solids dryer"),
            icon="equipment/dryer",
            slot=partial(self.addEquipment, UI_dryer),
            button=True, parent=toolboxContenido)
        l4.addWidget(botonDryer)
        actionWasher, botonWasher = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Solid washer"),
            icon="equipment/solidWasher",
            slot=partial(self.addEquipment, UI_solidWasher),
            button=True, parent=toolboxContenido)
        l4.addWidget(botonWasher)
        botonWasher.setEnabled(False)
        actionVacuum, botonVacuum = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Vacuum filter"),
            icon="equipment/vacuumfilter",
            slot=partial(self.addEquipment, UI_vacuum),
            button=True, parent=toolboxContenido)
        l4.addWidget(botonVacuum)
        botonVacuum.setEnabled(False)
        actionScrubber, botonScrubber = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Scrubber"),
            icon="equipment/scrubber",
            slot=partial(self.addEquipment, UI_scrubber),
            button=True, parent=toolboxContenido)
        l4.addWidget(botonScrubber)
        actionGravityChandler, botonGravityChandler = createAction(
            QtWidgets.QApplication.translate(
                "pychemqt", "Gravity settling chamber"),
            icon="equipment/gravityChamber",
            slot=partial(self.addEquipment, UI_gravityChamber),
            button=True, parent=toolboxContenido)
        l4.addWidget(botonGravityChandler)
        layouttoolbox.addItem(l4)
        layouttoolbox.addStretch(1)

        txt = QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Tools"))
        txt.setAlignment(QtCore.Qt.AlignCenter)
        layouttoolbox.addWidget(txt)
        l5 = FlowLayout()
        actionControler, botonControler = createAction(
            QtWidgets.QApplication.translate("pychemqt", "PID controller"),
            icon="equipment/controller",
            slot=partial(self.addEquipment, UI_solidWasher),
            button=True, parent=toolboxContenido)
        l5.addWidget(botonControler)
        botonControler.setEnabled(False)
        actionControlValve, botonControlValve = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Control valve"),
            icon="equipment/controlvalve",
            slot=partial(self.addEquipment, UI_vacuum),
            button=True, parent=toolboxContenido)
        l5.addWidget(botonControlValve)
        botonControlValve.setEnabled(False)
        actionSpreadsheet, botonSpreadsheet = createAction(
            QtWidgets.QApplication.translate(
                "pychemqt", "External spreadsheet module"),
            icon="equipment/spreadsheet",
            slot=partial(self.addEquipment, UI_spreadsheet),
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

        self.menuArchivo = QtWidgets.QMenu(
            QtWidgets.QApplication.translate("pychemqt", "&File"))
        self.menuRecentFiles = QtWidgets.QMenu(
            QtWidgets.QApplication.translate("pychemqt", "Open Recent Files"),
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

        self.menuEditar = QtWidgets.QMenu(
            QtWidgets.QApplication.translate("pychemqt", "&Edit"))
        self.menuEditar.aboutToShow.connect(self.aboutToShow_MenuEdit)
        self.menubar.addAction(self.menuEditar.menuAction())

        self.menuVer = QtWidgets.QMenu(
            QtWidgets.QApplication.translate("pychemqt", "&View"))
        self.menuVer.addAction(self.actionZoomOut)
        self.menuVer.addAction(self.actionZoom)
        self.menuVer.addAction(self.actionZoomIn)
        self.menuVer.addAction(actionOverviewWindow)
        self.menuVer.addSeparator()
        self.menuVer.addAction(actionVerStatus)
        self.menuVer.addAction(self.actionVerToolbar)
        self.menuVer.addAction(self.actionVerItem)
        self.menubar.addAction(self.menuVer.menuAction())

        self.menuObjetosGraficos = QtWidgets.QMenu(
            QtWidgets.QApplication.translate("pychemqt", "Plot"))
        self.menuObjetosGraficos.addAction(actionTexto)
        self.menuObjetosGraficos.addAction(actionCuadrado)
        self.menuObjetosGraficos.addAction(actionCircle)

        self.menuObjetosFlujo = QtWidgets.QMenu(
            QtWidgets.QApplication.translate("pychemqt", "Flux"))
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

        self.menuObjetosBasics = QtWidgets.QMenu(
            QtWidgets.QApplication.translate("pychemqt", "Basics"))
        self.menuObjetosBasics.addAction(actionTorreFUG)
        self.menuObjetosBasics.addAction(actionFlash)
        self.menuObjetosBasics.addAction(actionTorre)
        self.menuObjetosBasics.addAction(actionheatExchanger)
        self.menuObjetosBasics.addAction(actionhairpin)
        self.menuObjetosBasics.addAction(actionShellTube)
        self.menuObjetosBasics.addAction(actionFireHeater)
        self.menuObjetosBasics.addAction(actionReactor)

        self.menuObjetosSolids = QtWidgets.QMenu(
            QtWidgets.QApplication.translate("pychemqt", "Solids"))
        self.menuObjetosSolids.addAction(actionBaghouse)
        self.menuObjetosSolids.addAction(actionCentrifuge)
        self.menuObjetosSolids.addAction(actionCiclon)
        self.menuObjetosSolids.addAction(actionElectroPrecipitator)
        self.menuObjetosSolids.addAction(actionGrinder)
        self.menuObjetosSolids.addAction(actionDryer)
        self.menuObjetosSolids.addAction(actionWasher)
        self.menuObjetosSolids.addAction(actionVacuum)
        self.menuObjetosSolids.addAction(actionScrubber)
        self.menuObjetosSolids.addAction(actionGravityChandler)

        self.menuObjetosTools = QtWidgets.QMenu(
            QtWidgets.QApplication.translate("pychemqt", "Tools"))
        self.menuObjetosTools.addAction(actionControler)
        self.menuObjetosTools.addAction(actionControlValve)
        self.menuObjetosTools.addAction(actionSpreadsheet)

        self.menuPFD = QtWidgets.QMenu(
            QtWidgets.QApplication.translate("pychemqt", "&PFD"))
        self.actionResolution = createAction(
            QtWidgets.QApplication.translate("pychemqt", "Resolution"),
            slot=partial(self.dialogConfig, UI_confResolution),
            tip=QtWidgets.QApplication.translate(
                "pychemqt", "Defining PFD resolution dialog"),
            parent=self)
        self.menuPFD.addAction(self.actionResolution)
        self.menuPFD.addSeparator()
        self.menuPFD.addAction(self.menuObjetosGraficos.menuAction())
        self.menuPFD.addAction(self.menuObjetosFlujo.menuAction())
        self.menuPFD.addAction(self.menuObjetosBasics.menuAction())
        self.menuPFD.addAction(self.menuObjetosSolids.menuAction())
        self.menuPFD.addAction(self.menuObjetosTools.menuAction())
        self.menuPFD.addSeparator()
        self.menuPFD.addAction(saveAsImage)
        self.menubar.addAction(self.menuPFD.menuAction())

        self.menuPlot = QtWidgets.QMenu(
            QtWidgets.QApplication.translate("pychemqt", "Pl&ot"))
        for indice, grafico in enumerate(plots.__all__):
            self.menuPlot.addAction(grafico.title, partial(self.plot, indice))
        self.menubar.addAction(self.menuPlot.menuAction())

        self.menuCharts = QtWidgets.QMenu(
            QtWidgets.QApplication.translate("pychemqt", "Charts"), self)
        for titulo, lista in charts._all.items():
            menu = QtWidgets.QMenu(titulo, self)
            for grafico in lista:
                menu.addAction(grafico.title, partial(self.chart, grafico))
            self.menuCharts.addAction(menu.menuAction())

        self.menuAddComponent = QtWidgets.QMenu(
            QtWidgets.QApplication.translate("pychemqt", "New Component"))
        self.menuAddComponent.addAction(
            QtWidgets.QApplication.translate("pychemqt", "Component"),
            self.newcomponente)
        self.menuAddComponent.addAction(
            QtWidgets.QApplication.translate("pychemqt", "Pseudocomponent"),
            self.pseudocomponente)
        self.menuAddComponent.addSeparator()
        self.menuAddComponent.addAction(
            "Joback", partial(self.newComponent_Contribution, "Joback"))
        self.menuAddComponent.addAction(
            "Constantinou-Gani",
            partial(self.newComponent_Contribution, "Constantinou"))
        self.menuAddComponent.addAction(
            "Wilson-Jasperson",
            partial(self.newComponent_Contribution, "Wilson"))
        self.menuAddComponent.addAction(
            "Marrero-Pardillo",
            partial(self.newComponent_Contribution, "Marrero"))
        self.menuAddComponent.addAction(
            "Elliott", partial(self.newComponent_Contribution, "Elliott"))
        self.menuAddComponent.addAction(
            "Ambrose", partial(self.newComponent_Contribution, "Ambrose"))

        self.menuHerramientas = QtWidgets.QMenu(
            QtWidgets.QApplication.translate("pychemqt", "&Tools"))
        self.menuHerramientas.addAction(
            QtWidgets.QApplication.translate("pychemqt", "Component database"),
            self.verComponentes)
        self.menuHerramientas.addAction(self.menuAddComponent.menuAction())
        self.menuHerramientas.addSeparator()
        self.menuHerramientas.addAction(calculatorAction)
        self.menuHerramientas.addAction(terminalAction)
        self.menuHerramientas.addAction(conversorUnidadesAction)
        self.menuHerramientas.addAction(currencyAction)
        self.menuHerramientas.addAction(TablaPeriodicaAction)
        self.menuHerramientas.addAction(steamTablesAction)
        self.menuMEoS = UI_Tables.plugin(parent=self)
        self.menuHerramientas.addAction(self.menuMEoS.menuAction())
        self.menuHerramientas.addAction(psychrometricChartAction)
        self.menuHerramientas.addSeparator()
        self.menuHerramientas.addAction(self.menuCharts.menuAction())
        self.menuHerramientas.addSeparator()
        self.menuHerramientas.addAction(externalProgramAction)
        self.menubar.addAction(self.menuHerramientas.menuAction())

        self.menuVentana = QtWidgets.QMenu(
            QtWidgets.QApplication.translate("pychemqt", "&Window"))
        self.menuVentana.aboutToShow.connect(self.aboutToShow_MenuWindow)
        self.menubar.addAction(self.menuVentana.menuAction())

        self.menuAyuda = QtWidgets.QMenu(
            QtWidgets.QApplication.translate("pychemqt", "&Help"))
        self.menuAyuda.addAction(actionAyuda)
        self.menuAyuda.addSeparator()
        self.menuAyuda.addAction(actionAcerca_de)
        self.menuAyuda.addAction(actionAcerca_deQt)
        self.menubar.addAction(self.menuAyuda.menuAction())

        # Toolbox ListEquipment
        self.toolboxItem = QtWidgets.QDockWidget(
            QtWidgets.QApplication.translate("pychemqt", "Item"))
        self.toolboxItem.setObjectName("item")
        self.list = TreeEquipment()
        self.list.setSelectionMode(
            QtWidgets.QAbstractItemView.ExtendedSelection)
        self.list.itemSelectionChanged.connect(self.selectionChanged)
        self.list.customContextMenuRequested.connect(self.contextListMenu)
        self.toolboxItem.setWidget(self.list)
        self.addDockWidget(QtCore.Qt.DockWidgetArea(1), self.toolboxItem)
        self.toolboxItem.setFeatures(
            QtWidgets.QDockWidget.AllDockWidgetFeatures)
        self.toolboxItem.setAllowedAreas(
            QtCore.Qt.LeftDockWidgetArea | QtCore.Qt.RightDockWidgetArea)
        self.toolboxItem.visibilityChanged.connect(
            self.actionVerItem.setChecked)
        self.actionVerItem.triggered.connect(self.toolboxItem.setVisible)

        # Toolbox Status
        toolbox = QtWidgets.QDockWidget(
            QtWidgets.QApplication.translate("pychemqt", "Status"))
        toolbox.setObjectName("status")
        self.status = QtWidgets.QTextEdit()
        self.status.setReadOnly(True)
        toolbox.setWidget(self.status)
        self.addDockWidget(QtCore.Qt.DockWidgetArea(8), toolbox)
        toolbox.setFeatures(QtWidgets.QDockWidget.NoDockWidgetFeatures)
        toolbox.setAllowedAreas(QtCore.Qt.BottomDockWidgetArea)
        toolbox.visibilityChanged.connect(actionVerStatus.setChecked)
        actionVerStatus.triggered.connect(toolbox.setVisible)

        # Statusbar
        self.statusbar = QtWidgets.QStatusBar(self)
        self.setStatusBar(self.statusbar)
        self.statusbar.setMaximumHeight(20)
        self.progressBar = QtWidgets.QProgressBar()
        self.progressBar.setVisible(False)
        self.progressBar.setFixedWidth(80)
        self.statusbar.addPermanentWidget(self.progressBar)
        self.statusPosition = QtWidgets.QLabel(self)
        self.statusbar.addPermanentWidget(self.statusPosition)
        self.statusResolution = ClickableLabel(self)
        self.statusResolution.clicked.connect(
            partial(self.dialogConfig, UI_confResolution))
        self.statusbar.addPermanentWidget(self.statusResolution)
        self.statusThermo = ClickableLabel(self)
        self.statusThermo.clicked.connect(
            partial(self.dialogConfig, UI_confThermo))
        self.statusbar.addPermanentWidget(self.statusThermo)

        # TrayIcon
        self.systemtray = QtWidgets.QSystemTrayIcon(
            QtGui.QIcon(QtGui.QPixmap(IMAGE_PATH + "pychemqt.png")), self)
        self.systemtray.setToolTip("pychemqt")
        self.systemtray.setContextMenu(self.menuHerramientas)

        # Iniciar valores
        if os.path.isfile(QTSETTING_FILE):
            settings = QtCore.QSettings()
            self.recentFiles = settings.value("RecentFiles")
            self.lastFile = settings.value("LastFile")
            self.restoreGeometry(settings.value("Geometry"))
            self.restoreState(settings.value("MainWindow/State"))
            if self.recentFiles is None:
                self.recentFiles = []
            self.menuRecentFiles.setEnabled(bool(self.recentFiles))
        else:
            self.recentFiles = []
            self.menuRecentFiles.setEnabled(False)
            self.lastFile = []

        self.updateStatus("Loaded pychemqt")
        self.activeControl(False)
        self.changePreferenceLive()
        self.centralwidget.tabCloseRequested.connect(self.fileClose)

    @property
    def currentScene(self):
        if self.centralwidget.count():
            return self.currentView.scene()

    @property
    def currentView(self):
        if not self.centralwidget.count():
            return False
        return self.centralwidget.currentWidget().subWindowList()[0].widget()

    @property
    def currentMdi(self):
        if self.centralwidget.count():
            return self.centralwidget.currentWidget()

    @property
    def currentConfig(self):
        if self.centralwidget.count():
            return self.config[self.idTab]

    @property
    def currentFilename(self):
        if self.centralwidget.count():
            return self.filename[self.idTab]

    def getScene(self, indice):
        return self.getView(indice).scene()

    def getView(self, indice):
        return self.centralwidget.widget(indice).subWindowList()[0].widget()

    @property
    def idTab(self):
        if self.centralwidget.count():
            return self.centralwidget.currentIndex()

    def closeEvent(self, event=None):
        if self.okToContinue():
            for tab in range(self.centralwidget.count()):
                centralwidget = self.centralwidget.widget(tab)
                scene = centralwidget.subWindowList()[0].widget().scene()
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
        if not self.dirty:
            return True
        if ind != -1:
            ind = list(range(self.centralwidget.count()))
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
                QtWidgets.QApplication.translate(
                    "pychemqt", "Unsaved changes"),
                QtWidgets.QApplication.translate(
                    "pychemqt", "Save unsaved changes?"),
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No |
                QtWidgets.QMessageBox.Cancel,
                QtWidgets.QMessageBox.Yes)
            if dialog == QtWidgets.QMessageBox.Cancel:
                return False
            elif dialog == QtWidgets.QMessageBox.No:
                return True
            elif dialog == QtWidgets.QMessageBox.Yes:
                self.fileSaveAll()
                return True
        else:
            return True

# Menus configuration
    def aboutToShow_MenuEdit(self):
        self.menuEditar.clear()
        if self.currentScene:
            self.currentScene.addActions(self.menuEditar)
        self.menuEditar.addSeparator()
        self.menuEditar.addAction(self.actionWizard)
        self.menuEditar.addSeparator()
        self.menuEditar.addAction(self.actionComponentList)
        self.menuEditar.addAction(self.actionThermo)
        self.menuEditar.addAction(self.actionTransporte)
        self.menuEditar.addAction(self.actionUnidades)
        self.menuEditar.addAction(self.actioncostIndex)
        self.menuEditar.addSeparator()
        self.menuEditar.addAction(self.actionPreferencias)

    def aboutToShow_MenuWindow(self):
        self.menuVentana.clear()
        self.menuVentana.addAction(
            QtGui.QIcon(IMAGE_PATH + "button/arrow-left.png"),
            QtWidgets.QApplication.translate("pychemqt", "&Previous"),
            self.currentMdi.activatePreviousSubWindow,
            QtGui.QKeySequence.PreviousChild)
        self.menuVentana.addAction(
            QtGui.QIcon(IMAGE_PATH + "button/arrow-right.png"),
            QtWidgets.QApplication.translate("pychemqt", "&Next"),
            self.currentMdi.activateNextSubWindow,
            QtGui.QKeySequence.NextChild)
        self.menuVentana.addSeparator()
        self.menuVentana.addAction(
            QtGui.QIcon(IMAGE_PATH + "button/tile.png"),
            QtWidgets.QApplication.translate("pychemqt", "&Tile"),
            self.currentMdi.tileSubWindows)
        self.menuVentana.addAction(
            QtGui.QIcon(IMAGE_PATH + "button/cascade.png"),
            QtWidgets.QApplication.translate("pychemqt", "&Cascade"),
            self.currentMdi.cascadeSubWindows)
        self.menuVentana.addAction(
            QtWidgets.QApplication.translate("pychemqt", "&Restore All"),
            self.windowRestoreAll)
        self.menuVentana.addAction(
            QtWidgets.QApplication.translate("pychemqt", "&Iconize All"),
            self.windowMinimizeAll)
        self.menuVentana.addSeparator()
        for i, window in enumerate(self.currentMdi.subWindowList()):
            self.menuVentana.addAction("&%i %s" % (i+1, window.windowTitle()))
        self.menuVentana.addSeparator()
        self.menuVentana.addAction(
            QtGui.QIcon(IMAGE_PATH + "button/fileClose.png"),
            QtWidgets.QApplication.translate("pychemqt", "&Close window"),
            self.currentMdi.closeActiveSubWindow)

    def windowRestoreAll(self):
        for window in self.currentMdi.subWindowList():
            window.showNormal()

    def windowMinimizeAll(self):
        for window in self.currentMdi.subWindowList():
            window.showMinimized()

    def contextListMenu(self, event):
        contextMenu = QtWidgets.QMenu()
        self.currentScene.addActions(contextMenu, event)
        contextMenu.exec_(event)

    def aboutToShow_MenuRecentFiles(self):
        self.menuRecentFiles.clear()
        recentFiles = []
        for fname in self.recentFiles:
            if fname not in self.filename and QtCore.QFile.exists(fname):
                recentFiles.append(fname)
        if recentFiles:
            self.menuRecentFiles.addSeparator()
            for i, fname in enumerate(recentFiles):
                action = QtWidgets.QAction("&%d %s" % (i + 1, fname), self)
                action.setData(QtCore.QVariant(fname))
                action.triggered.connect(self.loadFile)
                self.menuRecentFiles.addAction(action)
        self.menuRecentFiles.addSeparator()
        self.menuRecentFiles.addAction(
            QtGui.QIcon(IMAGE_PATH + "button/clear.png"),
            QtWidgets.QApplication.translate("pychemqt", "Clear"),
            self.clearRecentFiles)

# File Manipulation
    def clearRecentFiles(self):
        self.recentFiles = []
        self.menuRecentFiles.setEnabled(False)

    def addRecentFile(self, fname):
        if fname and fname not in self.recentFiles:
            self.recentFiles.insert(0, fname)
            if len(self.recentFiles) > 9:
                self.recentFiles = self.recentFiles[:9]
        self.menuRecentFiles.setEnabled(len(self.recentFiles))

    def loadPFD(self, mdiarea):
        PFD = flujo.GraphicsView()
        PFD.setWindowTitle(
            QtWidgets.QApplication.translate("pychemqt", "Flow Diagram"))
        PFD.mouseMove.connect(self.updatePosition)
        scene = flujo.GraphicsScene(self)
        scene.selectionChanged.connect(self.selectionChanged)
        x = self.config[-1].getint("PFD", "x")
        y = self.config[-1].getint("PFD", "y")
        scene.setSceneRect(0, 0, x, y)
        PFD.setScene(scene)
        mdiarea.addSubWindow(PFD)
        PFD.show()

    def fileNew(self):
        UI_pychemqt.idNew += 1
        self.dirty.append(True)
        self.filename.append("")
        config = ConfigParser()
        config.add_section("PFD")
        config.set("PFD", "x", self.Preferences.get("PFD", "x"))
        config.set("PFD", "y", self.Preferences.get("PFD", "y"))
        self.config.append(config)
        mdiArea = QtWidgets.QMdiArea()
#        style=StyleCustom()
#        mdiArea.setStyle(style)
        self.centralwidget.addTab(
            mdiArea,
            QtWidgets.QApplication.translate("pychemqt", "New Project") +
            " %i" % UI_pychemqt.idNew)
        self.centralwidget.setCurrentIndex(self.centralwidget.count()-1)
        self.loadPFD(mdiArea)
        self.updateStatus(QtWidgets.QApplication.translate(
            "pychemqt", "New project created"))
        self.activeControl(True)
        self.wizard()

    def fileSave(self, indice=None):
        if indice is None:
            indice = self.idTab
        if not self.filename[indice]:
            self.fileSaveAs()
        else:
            with open(self.filename[indice], "w") as file:
                data = {}
                self.getScene(indice).project.writeToJSON(data)

                PFD = {}
                win = self.centralwidget.currentWidget().subWindowList()[0]
                PFD["x"] = win.pos().x()
                PFD["y"] = win.pos().y()
                PFD["height"] = win.size().height()
                PFD["width"] = win.size().width()
                self.currentScene.writeToJSON(PFD)
                data["PFD"] = PFD

                other = {}
                ventanas = self.centralwidget.currentWidget().subWindowList()
                for ind, win in enumerate(ventanas[1:]):
                    ventana = {}
                    ventana["class"] = window.widget().__class__
                    ventana["x"] = win.pos().x()
                    ventana["y"] = win.pos().y()
                    ventana["height"] = win.size().height()
                    vetnana["width"] = win.size().width()

                    widget = {}
                    win.widget().writeToJSON(widget)
                    ventana["window"] = widget
                    other["ind"] = ventana
                data["other"] = other

                json.dump(data, file, indent=4)

            self.dirty[self.idTab] = False
            self.updateStatus(
                QtWidgets.QApplication.translate("pychemqt", "Saved as") +
                " %s" % self.filename[indice])
            self.dirty[indice] = False
            self.saveControl()

    def fileSaveAs(self, indice=None):
        if indice is None:
            indice = self.idTab
        dir = self.filename[indice] if self.filename[indice] else "."
        fname = QtWidgets.QFileDialog.getSaveFileName(
            self,
            QtWidgets.QApplication.translate("pychemqt", "Save project"),
            dir, "pychemqt project file (*.pcq)")
        if fname:
            name = fname[0]
            if name.split(".")[-1] != "pcq":
                name += ".pcq"
            self.addRecentFile(name)
            self.filename[indice] = name
            self.fileSave(indice)
            self.centralwidget.setTabText(
                indice, os.path.splitext(os.path.basename(name))[0])

    def fileSaveAll(self):
        for tab in range(self.centralwidget.count()):
            if self.dirty[tab]:
                self.fileSave(tab)

    def fileOpen(self, fname=None):
        if not fname:
            if self.filename:
                dir = os.path.dirname(str(self.filename[-1]))
            else:
                dir = "."
            fname = QtWidgets.QFileDialog.getOpenFileName(
                self,
                QtWidgets.QApplication.translate("pychemqt", "Open project"),
                dir, QtWidgets.QApplication.translate(
                    "pychemqt", "pychemqt project file") + " (*.pcq)")[0]
        if fname:
            try:
                self.loadFile(fname)
                self.activeControl(True)
            except Exception as error:
                QtWidgets.QMessageBox.critical(
                    self,
                    QtWidgets.QApplication.translate("pychemqt", "Error"),
                    QtWidgets.QApplication.translate(
                        "pychemqt", "Failed to load file") + "\n" + fname)
                self.activeControl(False)
                del self.filename[-1]
                raise error
#                print(error)

    def loadFile(self, fname=None):
        if not fname:
            action = self.sender()
            if isinstance(action, QtWidgets.QAction):
                fname = str(action.data())
            else:
                return

        if fname:
            self.dirty.append(False)
            self.filename.append(fname)
            self.addRecentFile(fname)

            with open(fname, "r") as file:
                data = json.load(file)

            project = Project()
            project.readFromJSON(data)
            self.config.append(project.config)

            mdiArea = QtWidgets.QMdiArea()
            self.loadPFD(mdiArea)

            x = data["PFD"]["x"]
            y = data["PFD"]["y"]
            pos = QtCore.QPoint(x, y)
            width = data["PFD"]["width"]
            height = data["PFD"]["height"]
            size = QtCore.QSize(width, height)
            mdiArea.subWindowList()[0].move(pos)
            mdiArea.subWindowList()[0].resize(size)

            mdiArea.subWindowList()[0].widget().scene().readFromJSON(data)
            self.list.updateList(
                mdiArea.subWindowList()[0].widget().scene().objects)

            for ventana in data["other"]:
                name = ventana["class"]
                indice = other_window_names.index(name)
                widget = other_window[indice]()
                widget.readFromJSON(ventana)
                mdiArea.addSubWindow(widget)
                x = ventana["x"]
                y = ventana["y"]
                pos = QtCore.QPoint(x, y)
                h = ventana["height"]
                w = ventana["width"]
                size = QtCore.QSize(w, h)
                mdiArea.subWindowList()[ventana+1].move(pos)
                mdiArea.subWindowList()[ventana+1].resize(size)

            self.centralwidget.addTab(
                mdiArea, os.path.splitext(os.path.basename(str(fname)))[0])
            self.centralwidget.setCurrentIndex(self.centralwidget.count()-1)

            self.currentScene.project = project
            self.updateStatus(QtWidgets.QApplication.translate(
                "pychemqt", "Load") + " " + fname, True)

            self.activeControl(True)
            self.changeStatusThermo(self.config[self.idTab])

    def fileClose(self, int):
        if self.okToContinue(int):
            self.updateStatus(QtWidgets.QApplication.translate(
                "pychemqt", "Closed") + " " + self.currentFilename)
            self.centralwidget.removeTab(int)
            del self.dirty[int]
            del self.config[int]
            del self.filename[int]
            if self.centralwidget.count():
                self.activeControl(True)
            else:
                self.activeControl(False)
                self.list.clear()
                flujo.StreamItem.id = 0
                flujo.EquipmentItem.id = 0

# Configuration
    def wizard(self):
        dialog = wizard.Wizard(self.config[self.idTab])
        if dialog.exec_():
            self.updateConfig(dialog.value)
            self.updateStatus(QtWidgets.QApplication.translate(
                "pychemqt", "Project configuration"), True)
        else:
            self.updateConfig(wizard.Wizard.default())
            self.updateStatus(QtWidgets.QApplication.translate(
                "pychemqt", "Project configuration"), False)

    def updateConfig(self, config):
        self.config[self.idTab] = config
        self.currentScene.project.setConfig(config)
        self.dirty[self.idTab] = True
        x = config.getint("PFD", "x")
        y = config.getint("PFD", "y")
        self.currentScene.setSceneRect(0, 0, x, y)
        self.changeStatusThermo(config)
        setMainWindowConfig(config)

        # TODO: Delete this when its not necessary to run library isolated
        config.write(open(conf_dir+"pychemqtrc_temporal", "w"))

    def updateStatus(self, text, success=True):
        """Función que añade entradas al cuadro de status
        text: texto a mostrar
        success: boolean que indica si va todo bien"""
        if success:
            txt = QtWidgets.QApplication.translate("pychemqt", "Success")
            color = "#00aa00"
        else:
            txt = QtWidgets.QApplication.translate("pychemqt", "Failure")
            color = "#ff0000"
        self.status.append(
            '<b>' + time.strftime("%H:%M:%S", time.localtime()) +
            '</b> - ' + text + ' [<font color="%s">%s</font>]' % (color, txt))
        QtWidgets.QApplication.processEvents()

    def updatePosition(self, point):
        self.statusPosition.setText("(%i, %i)" % (point.x(), point.y()))

    def changeStatusThermo(self, config):
            if config.has_section("Thermo") and \
                    config.has_section("Components"):
                components = config.get("Components", "components")
                if config.getboolean("Thermo", "iapws") and \
                        config.getboolean("Thermo", "freesteam") and \
                        len(components) == 1 and components[0] == 62:
                    txt = "Freesteam"
                elif config.getboolean("Thermo", "iapws") and \
                        len(components) == 1 and components[0] == 62:
                    txt = "IAPWS97"
                elif config.getboolean("Thermo", "meos") and \
                        config.getboolean("Thermo", "refprop"):
                    txt = "Refprop"
                elif config.getboolean("Thermo", "meos") and \
                        config.getboolean("Thermo", "coolprop"):
                    txt = "CoolProp"
                elif config.getboolean("Thermo", "meos"):
                    txt = "MEoS"
                else:
                    txt = "K: %s  H: %s" % (
                        K[config.getint("Thermo", "K")].__status__,
                        H[config.getint("Thermo", "H")].__status__)

                self.statusThermo.setText(txt)
            if config.has_section("PFD"):
                x = config.getint("PFD", "x")
                y = config.getint("PFD", "y")
                self.statusResolution.setText("%i, %i" % (x, y))

    def dialogConfig(self, UIconfig):
        Dialog = UIconfig.Dialog(self.config[self.idTab])
        if Dialog.exec_():
            config = Dialog.value(self.config[self.idTab])
            self.updateConfig(config)
            self.saveControl()

    def costos(self):
        dialog = costIndex.Ui_CostIndex()
        dialog.exec_()

    def Preferencias(self):
        dialog = UI_Preferences.Preferences(self.Preferences)
        if dialog.exec_():
            preferences = dialog.value()
            preferences.write(open(conf_dir+"pychemqtrc", "w"))
            self.Preferences = ConfigParser()
            self.Preferences.read(conf_dir+"pychemqtrc")
            self.updateStatus(QtWidgets.QApplication.translate(
                "pychemqt", "pychemqt configuration change"), True)
            self.changePreferenceLive()
        else:
            self.updateStatus(QtWidgets.QApplication.translate(
                "pychemqt", "pychemqt configuration change"), False)

    def changePreferenceLive(self):
        if self.Preferences.getboolean("General", 'Tray'):
            self.systemtray.show()
        else:
            self.systemtray.hide()

    def activeControl(self, boolean):
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
        self.menuPlot.setEnabled(boolean)
        self.menuVentana.setEnabled(boolean)
        self.menuMEoS.setEnabled(boolean)
        self.toolboxItem.setEnabled(boolean)
        self.toolboxPalette.setEnabled(boolean)
        if boolean:
            self.saveControl()
        else:
            self.fileSaveAction.setEnabled(False)
            self.statusPosition.clear()
            self.statusResolution.clear()
            self.statusThermo.clear()

    def saveControl(self):
        self.fileSaveAction.setEnabled(self.dirty[self.idTab])
        self.tabModified(self.idTab)

    def tabModified(self, indice):
        if self.dirty[indice]:
            icon = QtGui.QIcon(IMAGE_PATH + "button/editor.png")
        else:
            icon = QtGui.QIcon()
        self.centralwidget.setTabIcon(indice, icon)

    def currentTabChanged(self, indice):
        if indice == -1:
            self.list.clear()
            self.activeControl(False)
        elif self.filename[-1]:
            scene = self.currentScene
            self.list.updateList(scene.objects)
            self.centralwidget.currentWidget().subWindowList()[0].widget()
            self.list.updateList(self.currentScene.objects)
            self.changeStatusThermo(self.config[self.idTab])
        else:
            self.list.clear()
            self.statusThermo.clear()
            self.statusResolution.clear()

    # Help
    def help(self):
        dialog = doi.ShowReference()
        dialog.exec_()

    def acerca(self):
        txt = QtWidgets.QApplication.translate(
            "pychemqt",
            "Software for simulate units operations in Chemical Engineering")
        QtWidgets.QMessageBox.about(
            self,
            QtWidgets.QApplication.translate("pychemqt", "About pychemqt"),
            """<b>pychemqt</b> v %s
            <p>Copyright &copy; 2012 Juan José Gómez Romera (jjgomera)<br>
            Licenced with GPL.v3
            <p>%s
            <p>Python %s - Qt %s - PyQt %s on %s""" % (
                __version__, txt, platform.python_version(),
                QtCore.QT_VERSION_STR, QtCore.PYQT_VERSION_STR,
                platform.system()))

    def acercaQt(self):
        QtWidgets.QMessageBox.aboutQt(
            self, QtWidgets.QApplication.translate("pychemqt", "About Qt"))

    # Tools
    def calculator(self):
        command = str(self.Preferences.get("Applications", 'Calculator'))
        os.system(command)

    def terminal(self):
        from tools import terminal
        terminal.XTerm(self.Preferences)

    def tablaPeriodica(self):
        from tools import qtelemental
        Tabla_Periodica = qtelemental.qtelemental()
        self.updateStatus(QtWidgets.QApplication.translate(
            "pychemqt", "Launched periodic table aplication"))
        Tabla_Periodica.exec_()

    def tablasVapor(self):
        SteamTables = UI_steamTables.Ui_SteamTables()
        self.updateStatus(QtWidgets.QApplication.translate(
            "pychemqt", "Launched steam-water properties aplication"))
        SteamTables.exec_()

    def diagramaPsicrometrico(self):
        Psychrometry = UI_psychrometry.UI_Psychrometry()
        self.updateStatus(QtWidgets.QApplication.translate(
            "pychemqt", "Launched humid air properties aplication"))
        Psychrometry.exec_()

    def externalPrograms(self):
        dialog = dependences.ShowDependences()
        dialog.exec_()

    def conversor_unidades(self):
        Conversor = UI_unitConverter.UI_unitConverter()
        self.updateStatus(QtWidgets.QApplication.translate(
            "pychemqt", "Launched unit converter aplication"))
        Conversor.exec_()

    def conversor_moneda(self):
        Conversor = moneda()
        self.updateStatus(QtWidgets.QApplication.translate(
            "pychemqt", "Launched currency converter aplication"))
        Conversor.exec_()

    def chart(self, grafico):
        dialog = grafico()
        self.updateStatus(QtWidgets.QApplication.translate(
            "pychemqt", "Show") + " " + grafico.title)
        dialog.exec_()

    def verComponentes(self):
        Base_datos = UI_databank.UI_databank()
        Base_datos.exec_()

    def newcomponente(self):
        dialog = viewComponents.View_Component()
        dialog.exec_()

    def pseudocomponente(self):
        Dialog = newComponent.Definicion_Petro()
        Dialog.exec_()

    def newComponent_Contribution(self, name):
        Dialog = newComponent.Ui_Contribution(name)
        Dialog.exec_()

    # PFD
    def plot(self, indice, x=None, y=None):
        grafico = plots.__all__[indice]()
        if grafico.exec_():
            self.currentMdi.addSubWindow(grafico.plot)
            grafico.plot.show()

#        indices, nombres, M=getComponents(config=self.config[self.idTab])
#        dialog=grafico(indices, nombres, x, y)
#        self.currentMdi.addSubWindow(dialog)
#        dialog.show()

    def savePFDImage(self):
        if self.filename[self.idTab]:
            dir = os.path.dirname(str(self.filename[self.idTab]))
        else:
            dir = "."
        fname = QtWidgets.QFileDialog.getSaveFileName(
            None,
            QtWidgets.QApplication.translate("pychemqt", "Save PFD as image"),
            dir, "Portable Network Graphics (*.png)")[0]
        if fname:
            rect = self.currentScene.sceneRect()
            img = QtGui.QImage(
                rect.width(), rect.height(),
                QtGui.QImage.Format_ARGB32_Premultiplied)
            p = QtGui.QPainter(img)
            self.currentScene.render(p)
            p.end()
            img.save(fname)

    def zoom(self, value=None):
        if value == "+":
            self.zoomValue.setValue(self.zoomValue.value()+5)
        elif value == "-":
            self.zoomValue.setValue(self.zoomValue.value()-5)
        elif value == "Dialog":
            value, bool = QtWidgets.QInputDialog.getInt(
                self,
                QtWidgets.QApplication.translate("pychemqt", "Zoom"),
                QtWidgets.QApplication.translate("pychemqt", "Zoom factor:"),
                self.zoomValue.value())
            if bool:
                self.zoomValue.setValue(value)
        else:
            self.currentView.zoom(value)

    def overview(self):
        PFD = flujo.GraphicsView(False)
        PFD.zoom(20)
        PFD.setWindowTitle(QtWidgets.QApplication.translate(
            "pychemqt", "Overview Window"))
        PFD.setScene(self.currentScene)
        self.centralwidget.currentWidget().addSubWindow(PFD)
        PFD.show()

    def selectionChanged(self):
        sender = self.sender()
        if isinstance(sender, QtWidgets.QTreeWidget):
            self.currentScene.blockSignals(True)
            self.currentScene.clearSelection()
            for element in self.list.selectedItems():
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
                    self.list.Stream.child(element.id-1).setSelected(True)
                elif isinstance(element, flujo.EquipmentItem) and \
                        element.tipo == "e":
                    self.list.Equipment.child(element.id-1).setSelected(True)
            self.list.blockSignals(False)

    def addText(self):
        dialog = flujo.TextItemDlg()
        if dialog.exec_():
            self.currentScene.waitClick(
                1, "txt", flujo.TextItem(dialog.editor.texto))

    def addItem(self, type, bool=True):
        if type == "square":
            obj = flujo.RectItem()
            num = 2
        elif type == "ellipse":
            obj = flujo.EllipseItem()
            num = 2
        elif type == "in":
            obj = flujo.EquipmentItem("in", None)
            num = 1
        elif type == "out":
            obj = flujo.EquipmentItem("out", None)
            num = 1
        elif type == "stream":
            if bool:
                obj = flujo.StreamItem()
                num = 2
            else:
                self.currentScene.clickCollector.quit()
                return
        self.currentScene.waitClick(num, type, obj)

    def addEquipment(self, equipo):
        equip = UI_equipments.index(equipo)
        object = flujo.EquipmentItem(equipo.__name__.split("_")[-1], equip)
        self.currentScene.waitClick(1, "equip", object)
