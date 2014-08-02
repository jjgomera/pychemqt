#!/usr/bin/python
# -*- coding: utf-8 -*-

from ConfigParser import ConfigParser
import os, time, platform, sys
from functools import partial

from PyQt4 import QtCore, QtGui
from scipy import arange
from scipy.optimize import fsolve

from UI import texteditor, newComponent, flujo, wizard, charts, plots, viewComponents
from UI.widgets import createAction, ClickableLabel, TreeEquipment, FlowLayout, Tabla, TabWidget
from lib.config import conf_dir, getComponents
from lib.project import Project
from lib.EoS import K, H
from lib import unidades, mEoS
from equipment import *
from tools import UI_confComponents, UI_Preferences, UI_confTransport, UI_confThermo, UI_confUnits, UI_confResolution, UI_databank, UI_Tables, UI_unitConverter, UI_steamTables, UI_psychrometry, costIndex, doi, dependences
from UI.conversor_unidades import moneda

__version__ = "0.1.0"


class UI_pychemqt(QtGui.QMainWindow):

    idNew=0
    config= []
    dirty = []
    filename=QtCore.QStringList()
    pfd=[]

    def __init__(self, parent=None):
        super(UI_pychemqt, self).__init__(parent)
        self.setWindowTitle("pychemqt")
        self.Preferences=ConfigParser()
        self.Preferences.read(conf_dir+"pychemqtrc")
        self.centralwidget = TabWidget()
        self.centralwidget.setTabsClosable(True)
        self.setCentralWidget(self.centralwidget)
        self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/pychemqt.png")))

        #Acciones
        fileNewAction = createAction(QtGui.QApplication.translate("pychemqt", "&New"), self.fileNew, QtGui.QKeySequence.New, os.environ["pychemqt"]+"/images/button/fileNew", QtGui.QApplication.translate("pychemqt", "Start new project"), parent=self)
        fileOpenAction = createAction(QtGui.QApplication.translate("pychemqt", "&Open"), self.fileOpen, QtGui.QKeySequence.Open, os.environ["pychemqt"]+"/images/button/fileOpen", QtGui.QApplication.translate("pychemqt", "Open project"), parent=self)
        self.fileSaveAction = createAction(QtGui.QApplication.translate("pychemqt", "&Save"), self.fileSave, QtGui.QKeySequence.Save, os.environ["pychemqt"]+"/images/button/fileSave", QtGui.QApplication.translate("pychemqt", "Save project"), parent=self)
        self.fileSaveAsAction = createAction(QtGui.QApplication.translate("pychemqt", "Save &as"), self.fileSaveAs, QtGui.QKeySequence.SaveAs, os.environ["pychemqt"]+"/images/button/fileSaveAs", QtGui.QApplication.translate("pychemqt", "Save project as"), parent=self)
        self.fileSaveAllAction = createAction(QtGui.QApplication.translate("pychemqt", "Save A&ll"), self.fileSaveAll, "", os.environ["pychemqt"]+"/images/button/fileSaveAll", QtGui.QApplication.translate("pychemqt", "Save all open project"), parent=self)
        self.fileCloseAction = createAction(QtGui.QApplication.translate("pychemqt", "&Close"), self.fileClose, QtGui.QKeySequence.Close, os.environ["pychemqt"]+"/images/button/fileClose", QtGui.QApplication.translate("pychemqt", "Close project"), parent=self)
        ExitAction = createAction(QtGui.QApplication.translate("pychemqt", "&Exit"), self.closeEvent, QtGui.QKeySequence.Quit, os.environ["pychemqt"]+"/images/button/exit", QtGui.QApplication.translate("pychemqt", "Salir de pychemqt"), parent=self)

        self.actionWizard = createAction(QtGui.QApplication.translate("pychemqt", "Wizard"), icon=os.environ["pychemqt"]+"/images/button/wizard", slot=self.wizard, tip=QtGui.QApplication.translate("pychemqt", "Launch configuration wizard"), parent=self)
        self.actionComponentList = createAction(QtGui.QApplication.translate("pychemqt", "Components list"), slot=partial(self.dialogConfig, UI_confComponents), tip=QtGui.QApplication.translate("pychemqt", "Defining componente list dialog"), parent=self)
        self.actionThermo = createAction(QtGui.QApplication.translate("pychemqt", "Thermodynamic properties"), slot=partial(self.dialogConfig, UI_confThermo), tip=QtGui.QApplication.translate("pychemqt", "Defining thermodynamic properties methods"), parent=self)
        self.actionTransporte = createAction(QtGui.QApplication.translate("pychemqt", "Transport properties"), slot=partial(self.dialogConfig, UI_confTransport), tip=QtGui.QApplication.translate("pychemqt", "Defining transport properties methods"), parent=self)
        self.actionUnidades = createAction(QtGui.QApplication.translate("pychemqt", "Units"), slot=partial(self.dialogConfig, UI_confUnits), tip=QtGui.QApplication.translate("pychemqt", "Defining preferred units"), parent=self)
        self.actioncostIndex = createAction(QtGui.QApplication.translate("pychemqt", "&Cost Index"), slot=self.costos, tip=QtGui.QApplication.translate("pychemqt", "Defining cost index"), parent=self)
        self.actionPreferencias = createAction(QtGui.QApplication.translate("pychemqt", "&Preferences"), slot=self.Preferencias, icon=os.environ["pychemqt"]+"/images/button/configure", shortcut=QtGui.QKeySequence.Preferences, tip=QtGui.QApplication.translate("pychemqt", "Defining general preferences"), parent=self)

        self.actionZoomIn = createAction(QtGui.QApplication.translate("pychemqt", "Zoom in"), slot=partial(self.zoom, "+"), icon=os.environ["pychemqt"]+"/images/button/zoomIn", shortcut=QtGui.QKeySequence.ZoomIn, parent=self)
        self.actionZoom = createAction(QtGui.QApplication.translate("pychemqt", "Zoom"), slot=partial(self.zoom, "Dialog"), icon=os.environ["pychemqt"]+"/images/button/zoomIn", parent=self)
        self.actionZoomOut = createAction(QtGui.QApplication.translate("pychemqt", "Zoom out"), slot=partial(self.zoom, "-"), icon=os.environ["pychemqt"]+"/images/button/zoomOut", shortcut=QtGui.QKeySequence.ZoomOut, parent=self)
        actionOverviewWindow = createAction(QtGui.QApplication.translate("pychemqt", "Overview window"), slot=self.overview, parent=self)
        actionVerStatus = createAction(QtGui.QApplication.translate("pychemqt", "Status"), shortcut="Ctrl+Alt+S", tip=QtGui.QApplication.translate("pychemqt", "Show/Hide status toolbar"), checkable=True, parent=self)
        self.actionVerToolbar = createAction(QtGui.QApplication.translate("pychemqt", "Palette"), shortcut="Ctrl+Alt+T", icon=os.environ["pychemqt"]+"/images/button/palette", tip=QtGui.QApplication.translate("pychemqt", "Show/Hide equipment palette"), checkable=True, parent=self)
        self.actionVerItem = createAction(QtGui.QApplication.translate("pychemqt", "Item"), icon=os.environ["pychemqt"]+"/images/button/list", tip=QtGui.QApplication.translate("pychemqt", "Show/Hide item list"), checkable=True, parent=self)

        calculatorAction = createAction(QtGui.QApplication.translate("pychemqt", "&Calculator"), slot=self.calculator, icon=os.environ["pychemqt"]+"/images/button/calculator", shortcut="F2", tip=QtGui.QApplication.translate("pychemqt", "Open system calculator"), parent=self)
        terminalAction = createAction(QtGui.QApplication.translate("pychemqt", "Python Shell"), shortcut="F3", slot=self.terminal, icon=os.environ["pychemqt"]+"/images/button/terminal", tip=QtGui.QApplication.translate("pychemqt", "Open system terminal"), parent=self)
        if sys.platform=="win32":
            terminalAction.setEnabled(False)
        conversorUnidadesAction = createAction(QtGui.QApplication.translate("pychemqt", "&Units converter"), slot=self.conversor_unidades, shortcut="F4", icon=os.environ["pychemqt"]+"/images/button/unitConverter", tip=QtGui.QApplication.translate("pychemqt", "Open Units converter"), parent=self)
        currencyAction = createAction(QtGui.QApplication.translate("pychemqt", "&Currency converter"), slot=self.conversor_moneda, shortcut="F5", icon=os.environ["pychemqt"]+"/images/button/currency", tip=QtGui.QApplication.translate("pychemqt", "Open Currency converter"), parent=self)
        TablaPeriodicaAction = createAction(QtGui.QApplication.translate("pychemqt", "&Periodic Table"), slot=self.tablaPeriodica, shortcut="F6", icon=os.environ["pychemqt"]+"/images/button/periodicTable", tip=QtGui.QApplication.translate("pychemqt", "Show a basic Mendeleiev periodic table"), parent=self)
        if not os.environ["Elemental"]:
            TablaPeriodicaAction.setEnabled(False)
        steamTablesAction = createAction(QtGui.QApplication.translate("pychemqt", "&Steam Tables"), slot=self.tablasVapor, shortcut="F7", icon=os.environ["pychemqt"]+"/images/button/steamTables", tip=QtGui.QApplication.translate("pychemqt", "Open a water-steam table and graphic application"), parent=self)
        psychrometricChartAction = createAction(QtGui.QApplication.translate("pychemqt", "&Psicrometric Chart"), slot=self.diagramaPsicrometrico, shortcut="F9", icon=os.environ["pychemqt"]+"/images/button/psychrometric", tip=QtGui.QApplication.translate("pychemqt", "Open a humid-air application"), parent=self)
        externalProgramAction = createAction(QtGui.QApplication.translate("pychemqt", "External Programs"), slot=self.externalPrograms, icon=os.environ["pychemqt"]+"/images/button/showPrograms", tip=QtGui.QApplication.translate("pychemqt", "Show External Programs Status"), parent=self)

        saveAsImage = createAction(QtGui.QApplication.translate("pychemqt", "Save PFD as image"), slot=self.savePFDImage, icon=os.environ["pychemqt"]+"/images/button/image", parent=self)

        actionAyuda = createAction(QtGui.QApplication.translate("pychemqt", "Help"), slot=self.help, icon=os.environ["pychemqt"]+"/images/button/help", parent=self)
        actionAcerca_de = createAction(QtGui.QApplication.translate("pychemqt", "About pychemqt"), slot=self.acerca, icon=os.environ["pychemqt"]+"/images/button/helpAbout", parent=self)
        actionAcerca_deQt = createAction(QtGui.QApplication.translate("pychemqt", "About Qt"), slot=self.acercaQt, icon=os.environ["pychemqt"]+"/images/button/AboutQt", parent=self)

        self.zoomValue=QtGui.QSpinBox()
        self.zoomValue.setSuffix("%")
        self.zoomValue.setRange(5, 1000)
        self.zoomValue.setValue(100)
        self.zoomValue.setSingleStep(5)
        self.zoomValue.valueChanged.connect(self.zoom)

        #Toolbar
        self.BarraArchivo = QtGui.QToolBar(QtGui.QApplication.translate("pychemqt", "File"), self)
        self.BarraArchivo.setObjectName("BarraArchivo")
        self.BarraArchivo.setIconSize(QtCore.QSize(16,16))
        self.BarraArchivo.addAction(fileNewAction)
        self.BarraArchivo.addAction(fileOpenAction)
        self.BarraArchivo.addAction(self.fileCloseAction)
        self.BarraArchivo.addAction(self.fileSaveAction)
        self.BarraArchivo.addAction(self.fileSaveAsAction)
        self.BarraArchivo.addAction(self.fileSaveAllAction)
        self.BarraArchivo.addAction(ExitAction)
        self.addToolBar(QtCore.Qt.TopToolBarArea,self.BarraArchivo)


        self.BarraVer = QtGui.QToolBar(QtGui.QApplication.translate("pychemqt", "View"), self)
        self.BarraVer.setObjectName("BarraVer")
        self.BarraVer.setIconSize(QtCore.QSize(16,16))
        self.BarraVer.addAction(self.actionZoomOut)
        self.BarraVer.addWidget(self.zoomValue)
        self.BarraVer.addAction(self.actionZoomIn)
        self.BarraVer.addSeparator()
        self.BarraVer.addAction(self.actionVerToolbar)
        self.BarraVer.addAction(self.actionVerItem)
        self.addToolBar(QtCore.Qt.TopToolBarArea,self.BarraVer)

        self.BarraHerramientas = QtGui.QToolBar(QtGui.QApplication.translate("pychemqt", "Tools"), self)
        self.BarraHerramientas.setObjectName("BarraHerramientas")
        self.BarraHerramientas.setIconSize(QtCore.QSize(16,16))
        self.BarraHerramientas.addAction(calculatorAction)
        self.BarraHerramientas.addAction(terminalAction)
        self.BarraHerramientas.addAction(conversorUnidadesAction)
        self.BarraHerramientas.addAction(currencyAction)
        self.BarraHerramientas.addAction(TablaPeriodicaAction)
        self.BarraHerramientas.addAction(steamTablesAction)
        self.BarraHerramientas.addAction(psychrometricChartAction)
        self.BarraHerramientas.addAction(externalProgramAction)
        self.addToolBar(QtCore.Qt.TopToolBarArea,self.BarraHerramientas)


        #Paleta Toolbox
        self.toolboxPalette = QtGui.QDockWidget(QtGui.QApplication.translate("pychemqt", "Equipos"))
        self.toolboxPalette.setObjectName("toolbox")
        toolboxContenido=QtGui.QWidget()
        self.toolboxPalette.setWidget(toolboxContenido)
        self.addDockWidget(QtCore.Qt.DockWidgetArea(1), self.toolboxPalette)
        self.toolboxPalette.setFeatures(QtGui.QDockWidget.AllDockWidgetFeatures)
        self.toolboxPalette.setAllowedAreas(QtCore.Qt.LeftDockWidgetArea|QtCore.Qt.RightDockWidgetArea)

        self.toolboxPalette.visibilityChanged.connect(self.actionVerToolbar.setChecked)
        self.actionVerToolbar.triggered.connect(self.toolboxPalette.setVisible)
        layouttoolbox=QtGui.QVBoxLayout(toolboxContenido)
        layouttoolbox.setContentsMargins(5, 5, 5, 5)

        txt=QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Plot"))
        txt.setAlignment(QtCore.Qt.AlignCenter)
        layouttoolbox.addWidget(txt)
        l1=FlowLayout()
        actionTexto, botonTexto=createAction(icon=os.environ["pychemqt"]+"/images/equipment/text", text=QtGui.QApplication.translate("pychemqt", "Insert text"), slot=self.addText, button=True, parent=toolboxContenido)
        l1.addWidget(botonTexto)
        actionCuadrado, botonCuadrado=createAction(icon=os.environ["pychemqt"]+"/images/equipment/square", text=QtGui.QApplication.translate("pychemqt", "Draw square"), slot=partial(self.addItem, "square"), button=True, parent=toolboxContenido)
        l1.addWidget(botonCuadrado)
        actionCircle, botonCircle=createAction(icon=os.environ["pychemqt"]+"/images/equipment/circle", text=QtGui.QApplication.translate("pychemqt", "Draw circle"), slot=partial(self.addItem, "ellipse"), button=True, parent=toolboxContenido)
        l1.addWidget(botonCircle)
        layouttoolbox.addItem(l1)
        layouttoolbox.addStretch(1)

        txt=QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Flux"))
        txt.setAlignment(QtCore.Qt.AlignCenter)
        layouttoolbox.addWidget(txt)
        l2=FlowLayout()
        actionEntrada, botonEntrada=createAction(icon=os.environ["pychemqt"]+"/images/equipment/in", text=QtGui.QApplication.translate("pychemqt", "Input"), slot=partial(self.addItem, "in"), button=True, parent=toolboxContenido)
        l2.addWidget(botonEntrada)
        actionCorriente, botonCorriente=createAction(icon=os.environ["pychemqt"]+"/images/equipment/stream", text=QtGui.QApplication.translate("pychemqt", "Stream"), slot=partial(self.addItem, "stream"), button=True, checkable=True, parent=toolboxContenido)
        l2.addWidget(botonCorriente)
        actionSalida, botonSalida=createAction(icon=os.environ["pychemqt"]+"/images/equipment/out", text=QtGui.QApplication.translate("pychemqt", "Output"), slot=partial(self.addItem, "out"), button=True, parent=toolboxContenido)
        l2.addWidget(botonSalida)
        actionDivider, botonDivider=createAction(icon=os.environ["pychemqt"]+"/images/equipment/divider", text=QtGui.QApplication.translate("pychemqt", "Divider"), slot=partial(self.addEquipment, UI_divider), button=True, parent=toolboxContenido)
        l2.addWidget(botonDivider)
        actionValve, botonValve=createAction(icon=os.environ["pychemqt"]+"/images/equipment/valve", text=QtGui.QApplication.translate("pychemqt", "Valve"), slot=partial(self.addEquipment, UI_valve), button=True, parent=toolboxContenido)
        l2.addWidget(botonValve)
        actionMixer, botonMixer=createAction(icon=os.environ["pychemqt"]+"/images/equipment/mixer", text=QtGui.QApplication.translate("pychemqt", "Mixer"), slot=partial(self.addEquipment, UI_mixer), button=True, parent=toolboxContenido)
        l2.addWidget(botonMixer)
        actionCompresor, botonCompresor=createAction(icon=os.environ["pychemqt"]+"/images/equipment/compressor", text=QtGui.QApplication.translate("pychemqt", "Compressor"), slot=partial(self.addEquipment, UI_compressor), button=True, parent=toolboxContenido)
        l2.addWidget(botonCompresor)
        actionTurbine, botonTurbine=createAction(icon=os.environ["pychemqt"]+"/images/equipment/turbine", text=QtGui.QApplication.translate("pychemqt", "Turbine"), slot=partial(self.addEquipment, UI_turbine), button=True, parent=toolboxContenido)
        l2.addWidget(botonTurbine)
        actionPump, botonPump=createAction(icon=os.environ["pychemqt"]+"/images/equipment/pump", text=QtGui.QApplication.translate("pychemqt", "Pump"), slot=partial(self.addEquipment, UI_pump), button=True, parent=toolboxContenido)
        l2.addWidget(botonPump)
        actionPipe, botonPipe=createAction(icon=os.environ["pychemqt"]+"/images/equipment/pipe", text=QtGui.QApplication.translate("pychemqt", "Pipe"), slot=partial(self.addEquipment, UI_pipe), button=True, parent=toolboxContenido)
        l2.addWidget(botonPipe)
        layouttoolbox.addItem(l2)
        layouttoolbox.addStretch(1)

        txt=QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Basics"))
        txt.setAlignment(QtCore.Qt.AlignCenter)
        layouttoolbox.addWidget(txt)
        l3=FlowLayout()
        actionTorreFUG, botonTorreFUG=createAction(icon=os.environ["pychemqt"]+"/images/equipment/columnFUG", text=QtGui.QApplication.translate("pychemqt", "Distillation tower (method FUG)"), slot=partial(self.addEquipment, UI_columnFUG), button=True, parent=toolboxContenido)
        l3.addWidget(botonTorreFUG)
        actionFlash, botonFlash=createAction(icon=os.environ["pychemqt"]+"/images/equipment/flash", text=QtGui.QApplication.translate("pychemqt", "Flash"), slot=partial(self.addEquipment, UI_flash), button=True, parent=toolboxContenido)
        l3.addWidget(botonFlash)
        actionTorre, botonTorre=createAction(icon=os.environ["pychemqt"]+"/images/equipment/tower", text=QtGui.QApplication.translate("pychemqt", "Distillation tower (exact method)"), slot=partial(self.addEquipment, UI_tower), button=True, parent=toolboxContenido)
        l3.addWidget(botonTorre)
        botonTorre.setEnabled(False)
        actionheatExchanger, botonheatExchanger=createAction(icon=os.environ["pychemqt"]+"/images/equipment/heatExchanger", text=QtGui.QApplication.translate("pychemqt", "Generic heat exchanger"), slot=partial(self.addEquipment, UI_heatExchanger), button=True, parent=toolboxContenido)
        l3.addWidget(botonheatExchanger)
        actionhairpin, botonhairpin=createAction(icon=os.environ["pychemqt"]+"/images/equipment/hairpin", text=QtGui.QApplication.translate("pychemqt", "Hairpin heat exchanger"), slot=partial(self.addEquipment, UI_hairpin), button=True, parent=toolboxContenido)
        l3.addWidget(botonhairpin)
        actionShellTube, botonShellTube=createAction(icon=os.environ["pychemqt"]+"/images/equipment/shellTube", text=QtGui.QApplication.translate("pychemqt", "Shell and tube heat exchanger"), slot=partial(self.addEquipment, UI_shellTube), button=True, parent=toolboxContenido)
        l3.addWidget(botonShellTube)
        actionFireHeater, botonFireHeater=createAction(icon=os.environ["pychemqt"]+"/images/equipment/fireHeater", text=QtGui.QApplication.translate("pychemqt", "Fired Heater heat exchanger"), slot=partial(self.addEquipment, UI_fireHeater), button=True, parent=toolboxContenido)
        l3.addWidget(botonFireHeater)
        actionReactor, botonReactor=createAction(icon=os.environ["pychemqt"]+"/images/equipment/reactor", text=QtGui.QApplication.translate("pychemqt", "Reactor"), slot=partial(self.addEquipment, UI_reactor), button=True, parent=toolboxContenido)
        l3.addWidget(botonReactor)
        botonReactor.setEnabled(False)
        layouttoolbox.addItem(l3)
        layouttoolbox.addStretch(1)

        txt=QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Solids"))
        txt.setAlignment(QtCore.Qt.AlignCenter)
        layouttoolbox.addWidget(txt)
        l4=FlowLayout()
        actionBaghouse, botonBaghouse=createAction(icon=os.environ["pychemqt"]+"/images/equipment/baghouse", text=QtGui.QApplication.translate("pychemqt", "Baghouse"), slot=partial(self.addEquipment, UI_baghouse), button=True, parent=toolboxContenido)
        l4.addWidget(botonBaghouse)
        actionCentrifuge, botonCentrifuge=createAction(icon=os.environ["pychemqt"]+"/images/equipment/centrifuge", text=QtGui.QApplication.translate("pychemqt", "Centrifuge"), slot=partial(self.addEquipment, UI_centrifuge), button=True, parent=toolboxContenido)
        l4.addWidget(botonCentrifuge)
        botonCentrifuge.setEnabled(False)
        actionCiclon, botonCiclon=createAction(icon=os.environ["pychemqt"]+"/images/equipment/ciclon", text=QtGui.QApplication.translate("pychemqt", "Cyclone"), slot=partial(self.addEquipment, UI_ciclon), button=True, parent=toolboxContenido)
        l4.addWidget(botonCiclon)
        actionElectroPrecipitator, botonElectroPrecipitator=createAction(icon=os.environ["pychemqt"]+"/images/equipment/electricPrecipitator", text=QtGui.QApplication.translate("pychemqt", "Electric precipitator"), slot=partial(self.addEquipment, UI_electricPrecipitator), button=True, parent=toolboxContenido)
        l4.addWidget(botonElectroPrecipitator)
        actionGrinder, botonGrinder=createAction(icon=os.environ["pychemqt"]+"/images/equipment/grinder", text=QtGui.QApplication.translate("pychemqt", "Grinder"), slot=partial(self.addEquipment, UI_grinder), button=True, parent=toolboxContenido)
        l4.addWidget(botonGrinder)
        botonGrinder.setEnabled(False)
        actionDryer, botonDryer=createAction(icon=os.environ["pychemqt"]+"/images/equipment/dryer", text=QtGui.QApplication.translate("pychemqt", "Solids dryer"), slot=partial(self.addEquipment, UI_dryer), button=True, parent=toolboxContenido)
        l4.addWidget(botonDryer)
        actionWasher, botonWasher=createAction(icon=os.environ["pychemqt"]+"/images/equipment/solidWasher", text=QtGui.QApplication.translate("pychemqt", "Solid washer"), slot=partial(self.addEquipment, UI_solidWasher), button=True, parent=toolboxContenido)
        l4.addWidget(botonWasher)
        botonWasher.setEnabled(False)
        actionVacuum, botonVacuum=createAction(icon=os.environ["pychemqt"]+"/images/equipment/vacuumfilter", text=QtGui.QApplication.translate("pychemqt", "Vacuum filter"), slot=partial(self.addEquipment, UI_vacuum), button=True, parent=toolboxContenido)
        l4.addWidget(botonVacuum)
        botonVacuum.setEnabled(False)
        actionVenturi, botonVenturi=createAction(icon=os.environ["pychemqt"]+"/images/equipment/venturi", text=QtGui.QApplication.translate("pychemqt", "Venturi"), slot=partial(self.addEquipment, UI_scrubber), button=True, parent=toolboxContenido)
        l4.addWidget(botonVenturi)
        botonVenturi.setEnabled(False)
        actionGravityChandler, botonGravityChandler=createAction(icon=os.environ["pychemqt"]+"/images/equipment/gravityChamber", text=QtGui.QApplication.translate("pychemqt", "Gravity chandler"), slot=partial(self.addEquipment, UI_gravityChamber), button=True, parent=toolboxContenido)
        l4.addWidget(botonGravityChandler)
        layouttoolbox.addItem(l4)
        layouttoolbox.addStretch(1)

        txt=QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Tools"))
        txt.setAlignment(QtCore.Qt.AlignCenter)
        layouttoolbox.addWidget(txt)
        l5=FlowLayout()
        actionControler, botonControler=createAction(icon=os.environ["pychemqt"]+"/images/equipment/controller", text=QtGui.QApplication.translate("pychemqt", "PID controller"), slot=partial(self.addEquipment, UI_solidWasher), button=True, parent=toolboxContenido)
        l5.addWidget(botonControler)
        botonControler.setEnabled(False)
        actionControlValve, botonControlValve=createAction(icon=os.environ["pychemqt"]+"/images/equipment/controlvalve", text=QtGui.QApplication.translate("pychemqt", "Control valve"), slot=partial(self.addEquipment, UI_vacuum), button=True, parent=toolboxContenido)
        l5.addWidget(botonControlValve)
        botonControlValve.setEnabled(False)
        actionSpreadsheet, botonSpreadsheet=createAction(icon=os.environ["pychemqt"]+"/images/equipment/spreadsheet", text=QtGui.QApplication.translate("pychemqt", "External spreadsheet module"), slot=partial(self.addEquipment, UI_spreadsheet), button=True, parent=toolboxContenido)
        if not os.environ["ezodf"]:
            actionSpreadsheet.setEnabled(False)
            botonSpreadsheet.setEnabled(False)
        l5.addWidget(botonSpreadsheet)
        layouttoolbox.addItem(l5)
        layouttoolbox.addStretch(10)


        #Menus
        self.menubar = QtGui.QMenuBar()
        self.setMenuBar(self.menubar)

        self.menuArchivo = QtGui.QMenu(QtGui.QApplication.translate("pychemqt", "&File"))
        self.menuRecentFiles = QtGui.QMenu(QtGui.QApplication.translate("pychemqt", "Open Recent Files"), self.menuArchivo)
        self.menuRecentFiles.aboutToShow.connect(self.aboutToShow_MenuRecentFiles)

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

        self.menuEditar = QtGui.QMenu(QtGui.QApplication.translate("pychemqt", "&Edit"))
        self.menuEditar.aboutToShow.connect(self.aboutToShow_MenuEdit)
        self.menubar.addAction(self.menuEditar.menuAction())

        self.menuVer = QtGui.QMenu(QtGui.QApplication.translate("pychemqt", "&View"))
        self.menuVer.addAction(self.actionZoomOut)
        self.menuVer.addAction(self.actionZoom)
        self.menuVer.addAction(self.actionZoomIn)
        self.menuVer.addAction(actionOverviewWindow)
        self.menuVer.addSeparator()
        self.menuVer.addAction(actionVerStatus)
        self.menuVer.addAction(self.actionVerToolbar)
        self.menuVer.addAction(self.actionVerItem)
        self.menubar.addAction(self.menuVer.menuAction())

        self.menuObjetosGraficos = QtGui.QMenu(QtGui.QApplication.translate("pychemqt", "Plot"))
        self.menuObjetosGraficos.addAction(actionTexto)
        self.menuObjetosGraficos.addAction(actionCuadrado)
        self.menuObjetosGraficos.addAction(actionCircle)

        self.menuObjetosFlujo = QtGui.QMenu(QtGui.QApplication.translate("pychemqt", "Flux"))
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

        self.menuObjetosBasics = QtGui.QMenu(QtGui.QApplication.translate("pychemqt", "Basics"))
        self.menuObjetosBasics.addAction(actionTorreFUG)
        self.menuObjetosBasics.addAction(actionFlash)
        self.menuObjetosBasics.addAction(actionTorre)
        self.menuObjetosBasics.addAction(actionheatExchanger)
        self.menuObjetosBasics.addAction(actionhairpin)
        self.menuObjetosBasics.addAction(actionShellTube)
        self.menuObjetosBasics.addAction(actionFireHeater)
        self.menuObjetosBasics.addAction(actionReactor)

        self.menuObjetosSolids = QtGui.QMenu(QtGui.QApplication.translate("pychemqt", "Solids"))
        self.menuObjetosSolids.addAction(actionBaghouse)
        self.menuObjetosSolids.addAction(actionCentrifuge)
        self.menuObjetosSolids.addAction(actionCiclon)
        self.menuObjetosSolids.addAction(actionElectroPrecipitator)
        self.menuObjetosSolids.addAction(actionGrinder)
        self.menuObjetosSolids.addAction(actionDryer)
        self.menuObjetosSolids.addAction(actionWasher)
        self.menuObjetosSolids.addAction(actionVacuum)
        self.menuObjetosSolids.addAction(actionVenturi)
        self.menuObjetosSolids.addAction(actionGravityChandler)

        self.menuObjetosTools = QtGui.QMenu(QtGui.QApplication.translate("pychemqt", "Tools"))
        self.menuObjetosTools.addAction(actionControler)
        self.menuObjetosTools.addAction(actionControlValve)
        self.menuObjetosTools.addAction(actionSpreadsheet)

        self.menuPFD = QtGui.QMenu(QtGui.QApplication.translate("pychemqt", "&PFD"))
        self.actionResolution = createAction(QtGui.QApplication.translate("pychemqt", "Resolution"), slot=partial(self.dialogConfig, UI_confResolution), tip=QtGui.QApplication.translate("pychemqt", "Defining PFD resolution dialog"), parent=self)
        self.menuPFD.addAction( self.actionResolution)
        self.menuPFD.addSeparator()
        self.menuPFD.addAction(self.menuObjetosGraficos.menuAction())
        self.menuPFD.addAction(self.menuObjetosFlujo.menuAction())
        self.menuPFD.addAction(self.menuObjetosBasics.menuAction())
        self.menuPFD.addAction(self.menuObjetosSolids.menuAction())
        self.menuPFD.addAction(self.menuObjetosTools.menuAction())
        self.menuPFD.addSeparator()
        self.menuPFD.addAction(saveAsImage)
        self.menubar.addAction(self.menuPFD.menuAction())

        self.menuPlot = QtGui.QMenu(QtGui.QApplication.translate("pychemqt", "Pl&ot"))
        for indice, grafico in enumerate(plots.__all__):
            self.menuPlot.addAction(grafico.title, partial(self.plot, indice))
        self.menubar.addAction(self.menuPlot.menuAction())

        self.menuCharts = QtGui.QMenu(QtGui.QApplication.translate("pychemqt", "Charts"), self)
        for titulo, lista in charts.__all__.iteritems():
            menu= QtGui.QMenu(titulo, self)
            for grafico in lista:
                menu.addAction(grafico.title, partial(self.chart, grafico))
            self.menuCharts.addAction(menu.menuAction())

        self.menuAddComponent = QtGui.QMenu(QtGui.QApplication.translate("pychemqt", "New Component"))
        self.menuAddComponent.addAction(QtGui.QApplication.translate("pychemqt", "Component"), self.newcomponente)
        self.menuAddComponent.addAction(QtGui.QApplication.translate("pychemqt", "Pseudocomponent"), self.pseudocomponente)
        self.menuAddComponent.addSeparator()
        self.menuAddComponent.addAction("Joback", partial(self.newComponent_Contribution, "Joback"))
        self.menuAddComponent.addAction("Constantinou-Gani", partial(self.newComponent_Contribution, "Constantinou"))
        self.menuAddComponent.addAction("Wilson-Jasperson", partial(self.newComponent_Contribution, "Wilson"))
        self.menuAddComponent.addAction("Marrero-Pardillo", partial(self.newComponent_Contribution, "Marrero"))
        self.menuAddComponent.addAction("Elliott", partial(self.newComponent_Contribution, "Elliott"))
        self.menuAddComponent.addAction("Ambrose", partial(self.newComponent_Contribution, "Ambrose"))

        self.menuHerramientas = QtGui.QMenu(QtGui.QApplication.translate("pychemqt", "&Tools"))
        self.menuHerramientas.addAction(QtGui.QApplication.translate("pychemqt", "Component database"), self.verComponentes)
        self.menuHerramientas.addAction(self.menuAddComponent.menuAction())
        self.menuHerramientas.addSeparator()
        self.menuHerramientas.addAction(calculatorAction)
        self.menuHerramientas.addAction(terminalAction)
        self.menuHerramientas.addAction(conversorUnidadesAction)
        self.menuHerramientas.addAction(currencyAction)
        self.menuHerramientas.addAction(TablaPeriodicaAction)
        self.menuHerramientas.addAction(steamTablesAction)
        self.menuMEoS=UI_Tables.plugin(QtGui.QApplication.translate("pychemqt", "MEoS properties"), parent=self)
        self.menuHerramientas.addAction(self.menuMEoS.menuAction())
        self.menuHerramientas.addAction(psychrometricChartAction)
        self.menuHerramientas.addSeparator()
        self.menuHerramientas.addAction(self.menuCharts.menuAction())
        self.menuHerramientas.addSeparator()
        self.menuHerramientas.addAction(externalProgramAction)
        self.menubar.addAction(self.menuHerramientas.menuAction())

        self.menuVentana = QtGui.QMenu(QtGui.QApplication.translate("pychemqt", "&Window"))
        self.menuVentana.aboutToShow.connect(self.aboutToShow_MenuWindow)
        self.menubar.addAction(self.menuVentana.menuAction())

        self.menuAyuda = QtGui.QMenu(QtGui.QApplication.translate("pychemqt", "&Help"))
        self.menuAyuda.addAction(actionAyuda)
        self.menuAyuda.addSeparator()
        self.menuAyuda.addAction(actionAcerca_de)
        self.menuAyuda.addAction(actionAcerca_deQt)
        self.menubar.addAction(self.menuAyuda.menuAction())

        #Toolbox ListEquipment
        self.toolboxItem = QtGui.QDockWidget(QtGui.QApplication.translate("pychemqt", "Item"))
        self.toolboxItem.setObjectName("item")
        self.list = TreeEquipment()
        self.list.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
        self.list.itemSelectionChanged.connect(self.selectionChanged)
        self.list.customContextMenuRequested.connect(self.contextListMenu)
        self.toolboxItem.setWidget(self.list)
        self.addDockWidget(QtCore.Qt.DockWidgetArea(1), self.toolboxItem)
        self.toolboxItem.setFeatures(QtGui.QDockWidget.AllDockWidgetFeatures)
        self.toolboxItem.setAllowedAreas(QtCore.Qt.LeftDockWidgetArea|QtCore.Qt.RightDockWidgetArea)
        self.toolboxItem.visibilityChanged.connect(self.actionVerItem.setChecked)
        self.actionVerItem.triggered.connect(self.toolboxItem.setVisible)

        #Toolbox Status
        toolbox = QtGui.QDockWidget(QtGui.QApplication.translate("pychemqt", "Status"))
        toolbox.setObjectName("status")
        self.status = QtGui.QTextEdit()
        self.status.setReadOnly(True)
        toolbox.setWidget(self.status)
        self.addDockWidget(QtCore.Qt.DockWidgetArea(8), toolbox)
        toolbox.setFeatures(QtGui.QDockWidget.NoDockWidgetFeatures)
        toolbox.setAllowedAreas(QtCore.Qt.BottomDockWidgetArea)
        toolbox.visibilityChanged.connect(actionVerStatus.setChecked)
        actionVerStatus.triggered.connect(toolbox.setVisible)

        #Statusbar
        self.statusbar = QtGui.QStatusBar(self)
        self.setStatusBar(self.statusbar)
        self.statusbar.setMaximumHeight(20)
        self.progressBar=QtGui.QProgressBar()
        self.progressBar.setVisible(False)
        self.progressBar.setFixedWidth(80)
        self.statusbar.addPermanentWidget(self.progressBar)
        self.statusPosition=QtGui.QLabel(self)
        self.statusbar.addPermanentWidget(self.statusPosition)
        self.statusResolution=ClickableLabel(self)
        self.statusResolution.clicked.connect(partial(self.dialogConfig, UI_confResolution))
        self.statusbar.addPermanentWidget(self.statusResolution)
        self.statusThermo=ClickableLabel(self)
        self.statusThermo.clicked.connect(partial(self.dialogConfig, UI_confThermo))
        self.statusbar.addPermanentWidget(self.statusThermo)


        #TrayIcon
        self.systemtray=QtGui.QSystemTrayIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/pychemqt.png")), self)
        self.systemtray.setToolTip("pychemqt")
        self.systemtray.setContextMenu(self.menuHerramientas)

        #Iniciar valores
        self.settings = QtCore.QSettings()
        self.recentFiles = self.settings.value("RecentFiles").toStringList()
        self.restoreGeometry(self.settings.value("Geometry").toByteArray())
        self.restoreState(self.settings.value("MainWindow/State").toByteArray())
        self.menuRecentFiles.setEnabled(len(self.recentFiles))

        self.updateStatus("Loaded pychemqt")
        self.activeControl(False)
        self.changePreferenceLive()
        self.centralwidget.tabCloseRequested.connect(self.fileClose)
        self.centralwidget.currentChanged.connect(self.currentTabChanged)

    @property
    def currentScene(self):
        if self.centralwidget.count():
            return self.currentView.scene()

    @property
    def currentView(self):
        if self.centralwidget.count():
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
                self.centralwidget.widget(tab).subWindowList()[0].widget().scene().clearSelection()
            settings = QtCore.QSettings()
            filename = QtCore.QVariant(self.filename) if self.filename else QtCore.QVariant()
            settings.setValue("LastFile", filename)
            recentFiles = QtCore.QVariant(self.recentFiles) if self.recentFiles else QtCore.QVariant()
            settings.setValue("RecentFiles", recentFiles)
            settings.setValue("Geometry", QtCore.QVariant(self.saveGeometry()))
            settings.setValue("MainWindow/State", QtCore.QVariant(self.saveState()))
            self.close()
        else:
            event.ignore()


    def okToContinue(self, ind=-1):
        if not self.dirty:
            return True
        if ind!=-1:
            ind=range(self.centralwidget.count())
        else:
            ind=[ind]
        dirty=False
        for tab in ind:
            if self.dirty[tab]:
                dirty=True
                break
        if dirty:
            dialog=QtGui.QMessageBox.question(self,
                        QtGui.QApplication.translate("pychemqt", "Unsaved changes"),
                        QtGui.QApplication.translate("pychemqt", "Save unsaved changes?"),
                        QtGui.QMessageBox.Yes| QtGui.QMessageBox.No|QtGui.QMessageBox.Cancel,
                        QtGui.QMessageBox.Yes)
            if dialog == QtGui.QMessageBox.Cancel:
                return False
            elif dialog == QtGui.QMessageBox.No:
                return True
            elif dialog == QtGui.QMessageBox.Yes:
                self.fileSaveAll()
                return True
        else: return True


#Menus configuration
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
        self.menuVentana.addAction(QtGui.QIcon(os.environ["pychemqt"]+"/images/button/arrow-left.png"), QtGui.QApplication.translate("pychemqt", "&Previous"), self.currentMdi.activatePreviousSubWindow, QtGui.QKeySequence.PreviousChild)
        self.menuVentana.addAction(QtGui.QIcon(os.environ["pychemqt"]+"/images/button/arrow-right.png"), QtGui.QApplication.translate("pychemqt", "&Next"), self.currentMdi.activateNextSubWindow, QtGui.QKeySequence.NextChild)
        self.menuVentana.addSeparator()
        self.menuVentana.addAction(QtGui.QIcon(os.environ["pychemqt"]+"/images/button/tile.png"), QtGui.QApplication.translate("pychemqt", "&Tile"), self.currentMdi.tileSubWindows)
        self.menuVentana.addAction(QtGui.QIcon(os.environ["pychemqt"]+"/images/button/cascade.png"), QtGui.QApplication.translate("pychemqt", "&Cascade"), self.currentMdi.cascadeSubWindows)
        self.menuVentana.addAction(QtGui.QApplication.translate("pychemqt", "&Restore All"), self.windowRestoreAll)
        self.menuVentana.addAction(QtGui.QApplication.translate("pychemqt", "&Iconize All"), self.windowMinimizeAll)
        self.menuVentana.addSeparator()
        for i, window in enumerate(self.currentMdi.subWindowList()):
            self.menuVentana.addAction("&%i " % (i+1)+str(window.windowTitle()))
        self.menuVentana.addSeparator()
        self.menuVentana.addAction(QtGui.QIcon(os.environ["pychemqt"]+"/images/button/fileClose.png"), QtGui.QApplication.translate("pychemqt", "&Close window"), self.currentMdi.closeActiveSubWindow)

    def windowRestoreAll(self):
        for window in self.currentMdi.subWindowList():
            window.showNormal()

    def windowMinimizeAll(self):
        for window in self.currentMdi.subWindowList():
            window.showMinimized()


    def contextListMenu(self, event):
        contextMenu= QtGui.QMenu()
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
                action = QtGui.QAction("&%d %s" % (i + 1, fname), self)
                action.setData(QtCore.QVariant(fname))
                action.triggered.connect(self.loadFile)
                self.menuRecentFiles.addAction(action)
        self.menuRecentFiles.addSeparator()
        self.menuRecentFiles.addAction(QtGui.QIcon(os.environ["pychemqt"]+"/images/button/clear.png"), QtGui.QApplication.translate("pychemqt", "Clear"), self.clearRecentFiles)


#File Manipulation
    def clearRecentFiles(self):
        self.recentFiles=QtCore.QStringList()
        self.menuRecentFiles.setEnabled(False)

    def addRecentFile(self, fname):
        if fname and fname not in self.recentFiles:
            self.recentFiles.prepend(QtCore.QString(fname))
            while self.recentFiles.count() > 9:
                self.recentFiles.takeLast()
        self.menuRecentFiles.setEnabled(len(self.recentFiles))

    def loadPFD(self, mdiarea):
        PFD = flujo.GraphicsView()
        PFD.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Flow Diagram"))
        PFD.mouseMove.connect(self.updatePosition)
        scene = flujo.GraphicsScene(self)
        scene.selectionChanged.connect(self.selectionChanged)
        scene.setSceneRect(0, 0, self.config[-1].getint("PFD", "x"), self.config[-1].getint("PFD", "y"))
        PFD.setScene(scene)
        mdiarea.addSubWindow(PFD)
        PFD.show()


    def fileNew(self):
        UI_pychemqt.idNew+=1
        self.dirty.append(True)
        self.filename.append("")
        config=ConfigParser()
        config.add_section("PFD")
        config.set("PFD", "x", self.Preferences.get("PFD", "x"))
        config.set("PFD", "y", self.Preferences.get("PFD", "y"))
        self.config.append(config)
        mdiArea = QtGui.QMdiArea()
#        style=StyleCustom()
#        mdiArea.setStyle(style)
        self.centralwidget.addTab(mdiArea, QtGui.QApplication.translate("pychemqt", "New Project")+" %i" % UI_pychemqt.idNew)
        self.centralwidget.setCurrentIndex(self.centralwidget.count()-1)
        self.loadPFD(mdiArea)
        self.updateStatus(QtGui.QApplication.translate("pychemqt", "Created new project"))
        self.activeControl(True)
        self.wizard()


    def fileSave(self, indice=None):
        if indice is None:
            indice=self.idTab
        if not self.filename[indice]:
            self.fileSaveAs()
        else:

            fh = QtCore.QFile(self.filename[indice])
            if not fh.open(QtCore.QIODevice.WriteOnly):
                raise IOError, unicode(fh.errorString())
            stream = QtCore.QDataStream(fh)
            stream.writeInt32(Project.MAGIC_NUMBER)
            stream.writeInt32(Project.FILE_VERSION)
            stream.setVersion(QtCore.QDataStream.Qt_4_2)

            self.getScene(indice).project.writeToStream(stream)

            stream << self.centralwidget.currentWidget().subWindowList()[0].pos()
            stream << self.centralwidget.currentWidget().subWindowList()[0].size()

            self.currentScene.saveToFile(stream)

            otras_ventanas=len(self.centralwidget.currentWidget().subWindowList())-1
            stream.writeInt32(otras_ventanas)
            for ventana in self.centralwidget.currentWidget().subWindowList()[1:]:
                stream.writeInt32(ventana.widget().Comp1.currentIndex())
                stream.writeInt32(ventana.widget().Comp2.currentIndex())
                stream.writeInt32(len(ventana.widget().x))
                for i in ventana.widget().x:
                    stream.writeFloat(i)
                for i in ventana.widget().y:
                    stream.writeFloat(i)
                stream << ventana.pos()
                stream << ventana.size()

            fh.close()
            self.dirty[self.idTab]=False

            self.updateStatus(QtGui.QApplication.translate("pychemqt", "Saved as %s" % self.filename[indice]))
            self.dirty[indice]=False
            self.saveControl()


    def fileSaveAs(self, indice=None):
        if indice is None:
            indice=self.idTab
        dir = self.filename[indice] if self.filename[indice] else "."
        fname = unicode(QtGui.QFileDialog.getSaveFileName(self,
                            QtGui.QApplication.translate("pychemqt", "Save project"), dir,
                            "pychemqt project file (*.pcq)"))
        if fname:
            if fname.split(".")[-1]!="pcq":
                fname+=".pcq"
            self.addRecentFile(fname)
            self.filename[indice] = fname
            self.fileSave(indice)
            self.centralwidget.setTabText(indice, os.path.splitext(os.path.basename(fname))[0])


    def fileSaveAll(self):
        for tab in range(self.centralwidget.count()):
            if self.dirty[tab]:
                self.fileSave(tab)

    def fileOpen(self, fname=None):
        if not fname:
            dir = os.path.dirname(str(self.filename[-1])) if self.filename else "."
            fname = unicode(QtGui.QFileDialog.getOpenFileName(self, QtGui.QApplication.translate("pychemqt", "Open project"), dir, QtGui.QApplication.translate("pychemqt", "pychemqt project file")+" (*.pcq)"))
        if fname:
            try:
                self.loadFile(fname)
                self.activeControl(True)
            except Exception as error:
                aviso=QtGui.QMessageBox.critical(self, QtGui.QApplication.translate("pychemqt", "Error"), QtGui.QApplication.translate("pychemqt", "Failed to load file")+"\n"+fname)
                self.activeControl(False)
                del self.filename[-1]
                print error


    def loadFile(self, fname=None):
        if not fname:
            action = self.sender()
            if isinstance(action, QtGui.QAction):
                fname = unicode(action.data().toString())
            else:
                return

        if fname:
            self.dirty.append(False)
            self.filename.append(fname)
            self.addRecentFile(fname)

            fh = QtCore.QFile(fname)
            if not fh.open(QtCore.QIODevice.ReadOnly):
                raise IOError, unicode(fh.errorString())
            stream = QtCore.QDataStream(fh)

            magic = stream.readInt32()
            if magic != Project.MAGIC_NUMBER:
                raise IOError, "unrecognized file type"
            version = stream.readInt32()
            if version < Project.FILE_VERSION:
                raise IOError, "old and unreadable file format"
            elif version > Project.FILE_VERSION:
                raise IOError, "new and unreadable file format"
            stream.setVersion(QtCore.QDataStream.Qt_4_2)

            project=Project()
            project.loadFromStream(stream)

            self.config.append(project.config)

            mdiArea = QtGui.QMdiArea()

            self.loadPFD(mdiArea)

            pos=QtCore.QPoint()
            size=QtCore.QSize()
            stream >> pos >> size
            mdiArea.subWindowList()[0].move(pos)
            mdiArea.subWindowList()[0].resize(size)

            mdiArea.subWindowList()[0].widget().scene().readFromFile(stream)
            self.list.updateList(mdiArea.subWindowList()[0].widget().scene().objects)

            otras_ventanas=stream.readInt32()
            for ventana in range(otras_ventanas):
                id1=stream.readInt32()
                id2=stream.readInt32()
                len=stream.readInt32()
                x=[]
                for i in range(len):
                    x.append(stream.readFloat())
                y=[]
                for i in range(len):
                    y.append(stream.readFloat())
                self.plot(0, x, y)

                pos=QtCore.QPoint()
                size=QtCore.QSize()
                stream >> pos >> size
                mdiArea.subWindowList()[ventana+1].move(pos)
                mdiArea.subWindowList()[ventana+1].resize(size)

            self.centralwidget.addTab(mdiArea, os.path.splitext(os.path.basename(str(fname)))[0])
            self.centralwidget.setCurrentIndex(self.centralwidget.count()-1)

            self.currentScene.project=project
            self.updateStatus(QtGui.QApplication.translate("pychemqt", "Load")+" "+ fname, True)

            self.activeControl(True)
            self.changeStatusThermo(self.config[self.idTab])


    def fileClose(self, int):
        if self.okToContinue(int):
            self.updateStatus(QtGui.QApplication.translate("pychemqt", "Closed")+" "+ self.currentFilename)
            self.centralwidget.removeTab(int)
            del self.dirty[int]
            del self.config[int]
            del self.filename[int]
            if self.centralwidget.count():
                self.activeControl(True)
            else:
                self.activeControl(False)
                self.list.clear()
                flujo.StreamItem.id=0
                flujo.EquipmentItem.id=0


#Configuration
    def wizard(self):
        dialog=wizard.Wizard(self.config[self.idTab])
        if dialog.exec_():
            self.updateConfig(dialog.value)
            self.updateStatus(QtGui.QApplication.translate("pychemqt", "Project configuration"), True)
        else:
            self.updateConfig(wizard.Wizard.default())
            self.updateStatus(QtGui.QApplication.translate("pychemqt", "Project configuration"), False)

    def updateConfig(self, config):
        self.config[self.idTab]=config
        self.currentScene.project.setConfig(config)
        self.dirty[self.idTab]=True
        self.currentScene.setSceneRect(0, 0, config.getint("PFD", "x"), config.getint("PFD", "y"))
        self.changeStatusThermo(config)

#TODO: Delete this when its not necessary to run library isolated
        config.write(open(conf_dir+"pychemqtrc_temporal", "w"))

    def updateStatus(self, text, success=True):
        """Función que añade entradas al cuadro de status
        text: texto a mostrar
        success: boolean que indica si va todo bien"""
        if success:
            txt=QtGui.QApplication.translate("pychemqt", "Success")
            color="#00aa00"
        else:
            txt=QtGui.QApplication.translate("pychemqt", "Failure")
            color="#ff0000"
        self.status.append( '<b>' + time.strftime("%H:%M:%S", time.localtime()) + '</b> - ' + text + ' [<font color="%s">%s</font>]' %(color, txt))
        QtGui.QApplication.processEvents()

    def updatePosition(self, point):
        self.statusPosition.setText("(%i, %i)" %(point.x(), point.y()))

    def changeStatusThermo(self, config):
            if config.has_section("Thermo") and config.has_section("Components"):
                components=eval(config.get("Components", "components"))
                if config.getboolean("Thermo", "iapws") and config.getboolean("Thermo", "freesteam") and len(components)==1 and components[0]==62:
                    txt="Freesteam"
                elif config.getboolean("Thermo", "iapws") and len(components)==1 and components[0]==62:
                    txt="IAPWS97"
                elif config.getboolean("Thermo", "meos") and config.getboolean("Thermo", "refprop"):
                    txt="Refprop"
                elif config.getboolean("Thermo", "meos") and config.getboolean("Thermo", "coolprop"):
                    txt="CoolProp"
                elif config.getboolean("Thermo", "meos"):
                    txt="MEoS"
                else:
                    txt="K: %s  H: %s" % (K[config.getint("Thermo","K")].__status__, H[config.getint("Thermo","H")].__status__)

                self.statusThermo.setText(txt)
            if config.has_section("PFD"):
                self.statusResolution.setText("%i, %i" % (config.getint("PFD","x"), config.getint("PFD","y")))

    def dialogConfig(self, UIconfig):
        Dialog = UIconfig.Dialog(self.config[self.idTab])
        if Dialog.exec_():
            config=Dialog.value(self.config[self.idTab])
            self.updateConfig(config)
            self.saveControl()

    def costos(self):
        dialog = costIndex.Ui_CostIndex()
        dialog.exec_()

    def Preferencias(self):
        dialog = UI_Preferences.Preferences(self.Preferences)
        if dialog.exec_():
            preferences=dialog.value()
            preferences.write(open(conf_dir+"pychemqtrc", "w"))
            self.Preferences=ConfigParser()
            self.Preferences.read(conf_dir+"pychemqtrc")
            self.updateStatus(QtGui.QApplication.translate("pychemqt", "pychemqt configuration change"), True)
            self.changePreferenceLive()
        else:
            self.updateStatus(QtGui.QApplication.translate("pychemqt", "pychemqt configuration change"), False)

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
            icon=QtGui.QIcon(os.environ["pychemqt"]+"/images/button/editor.png")
        else:
            icon=QtGui.QIcon()
        self.centralwidget.setTabIcon(indice, icon)

    def currentTabChanged(self, indice):
#        flujo.StreamItem.id=0
#        flujo.EquipmentItem.id=0
        if indice==-1:
            self.list.clear()
            self.activeControl(False)
        elif self.filename[-1]:
            self.list.updateList(self.centralwidget.currentWidget().subWindowList()[0].widget().scene().objects)
#            print self.currentScene.objects
            self.centralwidget.currentWidget().subWindowList()[0].widget()
            self.list.updateList(self.currentScene.objects)
            self.changeStatusThermo(self.config[self.idTab])
        else:
            self.list.clear()
            self.statusThermo.clear()
            self.statusResolution.clear()

#Help
    def help(self):
        dialog=doi.ShowReference()
        dialog.exec_()


    def acerca(self):
        txt= QtGui.QApplication.translate("pychemqt", "Software for simulate units operations in Chemical Engineering")
        QtGui.QMessageBox.about(self, QtGui.QApplication.translate("pychemqt", "About pychemqt"),
                u"""<b>pychemqt</b> v %s
                <p>Copyright &copy; 2012 Juan José Gómez Romera (jjgomera)<br>
                Licenced with GPL.v3
                <p>%s
                <p>Python %s - Qt %s - PyQt %s on %s""" % (
                __version__, txt, platform.python_version(),
                QtCore.QT_VERSION_STR, QtCore.PYQT_VERSION_STR, platform.system()))

    def acercaQt(self):
        QtGui.QMessageBox.aboutQt(self,QtGui.QApplication.translate("pychemqt", "About Qt"))


#Tools
    def calculator(self):
        command=str(self.Preferences.get("Applications", 'Calculator'))
        os.system(command)

    def terminal(self):
        from tools import terminal
        shell = terminal.XTerm(self.Preferences)
        shell.show()

    def tablaPeriodica(self):
        from tools import qtelemental
        Tabla_Periodica = qtelemental.qtelemental()
        self.updateStatus(QtGui.QApplication.translate("pychemqt", "Launched periodic table aplication"))
        Tabla_Periodica.exec_()

    def tablasVapor(self):
        SteamTables = UI_steamTables.Ui_SteamTables()
        self.updateStatus(QtGui.QApplication.translate("pychemqt", "Launched steam-water properties aplication"))
        SteamTables.exec_()

    def diagramaPsicrometrico(self):
        Psychrometry=UI_psychrometry.UI_Psychrometry()
        self.updateStatus(QtGui.QApplication.translate("pychemqt", "Launched humid air properties aplication"))
        Psychrometry.show()

    def externalPrograms(self):
        dialog=dependences.ShowDependences()
        dialog.exec_()

    def conversor_unidades(self):
        Conversor = UI_unitConverter.UI_unitConverter()
        self.updateStatus(QtGui.QApplication.translate("pychemqt", "Launched unit converter aplication"))
        Conversor.exec_()

    def conversor_moneda(self):
        Conversor = moneda()
        self.updateStatus(QtGui.QApplication.translate("pychemqt", "Launched currency converter aplication"))
        Conversor.exec_()

    def chart(self, grafico):
        dialog=grafico()
        self.updateStatus(QtGui.QApplication.translate("pychemqt", "Show")+" "+grafico.title)
        dialog.exec_()

    def verComponentes(self):
        Base_datos = UI_databank.UI_databank()
        Base_datos.exec_()

    def newcomponente(self):
        dialog=viewComponents.View_Component()
        dialog.exec_()

    def pseudocomponente(self):
        Dialog = newComponent.Definicion_Petro()
        Dialog.exec_()

    def newComponent_Contribution(self, name):
        Dialog = newComponent.Ui_Contribution(name)
        Dialog.exec_()


#PFD
    def plot(self, indice, x=None, y=None):
        grafico=plots.__all__[indice]
        indices, nombres, M=getComponents(config=self.config[self.idTab])
        dialog=grafico(indices, nombres, x, y)
        self.currentMdi.addSubWindow(dialog)
        dialog.show()

    def savePFDImage(self):
        dir = os.path.dirname(str(self.filename[self.idTab])) if self.filename[self.idTab] else "."
        fname = unicode(QtGui.QFileDialog.getSaveFileName(None, QtGui.QApplication.translate("pychemqt", "Save PFD as image"), dir, "Portable Network Graphics (*.png)"))
        if fname:
            rect=self.currentScene.sceneRect()
            img=QtGui.QImage(rect.width(), rect.height(),QtGui.QImage.Format_ARGB32_Premultiplied)
            p=QtGui.QPainter(img)
            self.currentScene.render(p)
            p.end()
            img.save(fname)

    def zoom(self, value=None):
        if value=="+":
            self.zoomValue.setValue(self.zoomValue.value()+5)
        elif value=="-":
            self.zoomValue.setValue(self.zoomValue.value()-5)
        elif value=="Dialog":
            value, bool=QtGui.QInputDialog.getInteger(self, QtGui.QApplication.translate("pychemqt", "Zoom"), QtGui.QApplication.translate("pychemqt", "Zoom factor:"), self.zoomValue.value())
            if bool:
                self.zoomValue.setValue(value)
        else:
            self.currentView.zoom(value)

    def overview(self):
        PFD = flujo.GraphicsView(False)
        PFD.zoom(20)
        PFD.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Overview Window"))
        PFD.setScene(self.currentScene)
        self.centralwidget.currentWidget().addSubWindow(PFD)
        PFD.show()



    def selectionChanged(self):
        sender=self.sender()
        if isinstance(sender, QtGui.QTreeWidget):
            self.currentScene.blockSignals(True)
            self.currentScene.clearSelection()
            for element in self.list.selectedItems():
                if element.parent()==self.list.Equipment:
                    self.currentScene.getObject("e", int(element.text(0).split(" ")[0])).setSelected(True)
                if element.parent()==self.list.Stream:
                    self.currentScene.getObject("s", int(element.text(0))).setSelected(True)
            self.currentScene.blockSignals(False)
        elif isinstance(sender, flujo.GraphicsScene):
            self.list.blockSignals(True)
            self.list.clearSelection()
            for element in sender.selectedItems():
                if isinstance(element, flujo.StreamItem):
                    self.list.Stream.child(element.id-1).setSelected(True)
                elif isinstance(element, flujo.EquipmentItem) and element.tipo=="e":
                    self.list.Equipment.child(element.id-1).setSelected(True)
            self.list.blockSignals(False)


    def addText(self):
        dialog = flujo.TextItemDlg()
        if dialog.exec_():
            self.currentScene.waitClick(1, "txt", flujo.TextItem(dialog.editor.texto))

    def addItem(self, type, bool=True):
        if type=="square":
            object=flujo.RectItem()
            num=2
        elif type=="ellipse":
            object=flujo.EllipseItem()
            num=2
        elif type=="in":
            object=flujo.EquipmentItem("in", None)
            num=1
        elif type=="out":
            object=flujo.EquipmentItem("out", None)
            num=1
        elif type=="stream":
            if bool:
                object = flujo.StreamItem()
                num=2
            else:
                self.currentScene.clickCollector.quit()
                return
        self.currentScene.waitClick(num, type, object)

    def addEquipment(self, equipo):
        equip=UI_equipments.index(equipo)
        object=flujo.EquipmentItem(equipo.__name__.split("_")[-1], equip)
        self.currentScene.waitClick(1, "equip", object)

