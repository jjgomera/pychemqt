#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Shell and tube heat exchanger equipment dialog
###############################################################################

from functools import partial

from PyQt5 import QtWidgets


from lib.unidades import Length, ThermalConductivity, Pressure, Currency
from UI.widgets import Entrada_con_unidades
from equipment.parents import UI_equip
from equipment.widget import FoulingWidget, Dialog_Finned
from equipment.UI_pipe import Catalogo_Materiales_Dialog
from equipment.heatExchanger import Shell_Tube
from tools.costIndex import CostData


class Dialog_Methods(QtWidgets.QDialog):
    """Dialog to select methods calculations"""

    def __init__(self, equipment, parent=None):
        super(Dialog_Methods, self).__init__(parent=parent)
        self.setWindowTitle(QtCore.QCoreApplication.translate(
            "pychemqt", "Specify calculation methods"))
        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate(
            "pychemqt", "Tubeside laminar flow")), 1, 1)
        self.tubesideLaminar = QtWidgets.QComboBox()
        for txt in equipment.TEXT_METHOD_TUBE_LAMINAR:
            self.tubesideLaminar.addItem(txt)
        layout.addWidget(self.tubesideLaminar, 1, 2)

        layout.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate(
            "pychemqt", "Tubeside turbulent flow")), 2, 1)
        self.tubesideTurbulent = QtWidgets.QComboBox()
        for txt in equipment.TEXT_METHOD_TUBE_TURBULENT:
            self.tubesideTurbulent.addItem(txt)
        layout.addWidget(self.tubesideTurbulent, 2, 2)

        layout.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "ShellSide")), 3, 1)
        self.shellSide = QtWidgets.QComboBox()
        for txt in equipment.TEXT_METHOD_SHELL:
            self.shellSide.addItem(txt)
        layout.addWidget(self.shellSide, 3, 2)

        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Cancel | QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox, 10, 1, 1, 2)

        if equipment.kwargs:
            self.tubesideLaminar.setCurrentIndex(
                equipment.kwargs["tubesideLaminar"])
            self.tubesideTurbulent.setCurrentIndex(
                equipment.kwargs["tubesideTurbulent"])
            self.shellSide.setCurrentIndex(
                equipment.kwargs["shellsideSensible"])

    @property
    def kwargs(self):
        kwargs = {}
        kwargs["tubesideLaminar"] = self.tubesideLaminar.currentIndex()
        kwargs["tubesideTurbulent"] = self.tubesideTurbulent.currentIndex()
        kwargs["shellsideSensible"] = self.shellSide.currentIndex()
        return kwargs


class UI_equipment(UI_equip):
    """Shell and tube heat exchanger edition dialog"""
    Equipment = Shell_Tube()

    def __init__(self, equipment=None, parent=None):
        """
        equipment: Initial equipment instance to model
        """
        super(UI_equipment, self).__init__(Shell_Tube, parent=parent)

        # Input tab
        self.addEntrada(QtCore.QCoreApplication.translate("pychemqt", "Tubes"),
                        "entradaTubo")
        self.addEntrada(QtCore.QCoreApplication.translate("pychemqt", "Shell"),
                        "entradaCarcasa")

        # Model tab
        tab = QtWidgets.QWidget()
        self.tabWidget.insertTab(
            1, tab, QtCore.QCoreApplication.translate("pychemqt", "Model"))
        lyt = QtWidgets.QGridLayout(tab)
        lyt.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Class")), 2, 1)
        self.class_ = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_CLASS:
            self.class_.addItem(txt)
        self.class_.currentIndexChanged.connect(
            partial(self.changeParams, "class_"))
        lyt.addWidget(self.class_, 2, 2)
        lyt.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Front end head")), 3, 1)
        self.frontHead = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_FRONTHEAD:
            self.frontHead.addItem(txt)
        self.frontHead.currentIndexChanged.connect(
            partial(self.changeParams, "frontHead"))
        lyt.addWidget(self.frontHead, 3, 2, 1, 3)
        lyt.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Shell type")), 4, 1)
        self.shell = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_SHELL:
            self.shell.addItem(txt)
        self.shell.currentIndexChanged.connect(
            partial(self.changeParams, "shell"))
        lyt.addWidget(self.shell, 4, 2)
        lyt.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Rear end head")), 5, 1)
        self.rearHead = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_REARHEAD:
            self.rearHead.addItem(txt)

        self.rearHead.currentIndexChanged.connect(
            partial(self.changeParams, "rearHead"))
        lyt.addWidget(self.rearHead, 5, 2, 1, 2)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            6, 1, 1, 6)
        lyt.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Orientation")), 7, 1)
        self.orientacion = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_ORIENTATION:
            self.orientacion.addItem(txt)
        self.orientacion.currentIndexChanged.connect(
            partial(self.changeParams, "orientation"))
        lyt.addWidget(self.orientacion, 7, 2)

        botonMetodos = QtWidgets.QPushButton(
            QtCore.QCoreApplication.translate("pychemqt", "Calculation methods"))
        botonMetodos.clicked.connect(self.selectMethods)
        lyt.addWidget(botonMetodos, 9, 1)
        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding),
            10, 1, 1, 6)

        # Tubes tab
        tab = QtWidgets.QWidget()
        self.tabWidget.insertTab(
            2, tab, QtCore.QCoreApplication.translate("pychemqt", "Tubes"))
        lyt = QtWidgets.QGridLayout(tab)
        lyt.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Number of tubes")), 1, 1)
        self.NTubes = Entrada_con_unidades(int, width=60, spinbox=True, step=1)
        self.NTubes.valueChanged.connect(partial(self.changeParams, "NTube"))
        lyt.addWidget(self.NTubes, 1, 2)
        lyt.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate(
            "pychemqt", "Number of tube passes")), 2, 1)
        self.NPases = Entrada_con_unidades(int, width=60, spinbox=True, step=1)
        self.NPases.valueChanged.connect(partial(self.changeParams, "NPases"))
        lyt.addWidget(self.NPases, 2, 2)
        lyt.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Tube length")), 3, 1)
        self.LTube = Entrada_con_unidades(Length)
        self.LTube.valueChanged.connect(partial(self.changeParams, "LTube"))
        lyt.addWidget(self.LTube, 3, 2)
        lyt.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate(
            "pychemqt", "Tube external diameter")), 4, 1)
        self.DeTube = Entrada_con_unidades(Length, "pipeDiameter")
        self.DeTube.valueChanged.connect(partial(self.changeParams, "DeTube"))
        lyt.addWidget(self.DeTube, 4, 2)
        lyt.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Tube thickness")), 5, 1)
        self.wTube = Entrada_con_unidades(Length, "Thickness")
        self.wTube.valueChanged.connect(partial(self.changeParams, "wTube"))
        lyt.addWidget(self.wTube, 5, 2)
        lyt.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Tube roughness")), 6, 1)
        self.rTube = Entrada_con_unidades(Length, "Thickness")
        self.rTube.valueChanged.connect(partial(self.changeParams, "rTube"))
        lyt.addWidget(self.rTube, 6, 2)
        lyt.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate(
            "pychemqt", "Thermal conductivity")), 7, 1)
        self.kTube = Entrada_con_unidades(ThermalConductivity)
        self.kTube.valueChanged.connect(partial(self.changeParams, "kTube"))
        lyt.addWidget(self.kTube, 7, 2)
        self.buttonPipe = QtWidgets.QPushButton(
            QtCore.QCoreApplication.translate("pychemqt", "Pipe Database"))
        self.buttonPipe.clicked.connect(self.showMaterial)
        lyt.addWidget(self.buttonPipe, 4, 4, 4, 1)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            8, 1, 1, 6)

        lyt.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Tube pattern")), 9, 1)
        self.distribucionTube = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_DISTRIBUTION_TUBE:
            self.distribucionTube.addItem(txt)
        self.distribucionTube.currentIndexChanged.connect(
            partial(self.changeParams, "distribucionTube"))
        lyt.addWidget(self.distribucionTube, 9, 2)
        lyt.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Tube pitch")), 10, 1)
        self.pitch = Entrada_con_unidades(Length)
        self.pitch.valueChanged.connect(partial(self.changeParams, "pitch"))
        lyt.addWidget(self.pitch, 10, 2)
        lyt.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Fin Tube")), 11, 1)
        self.buttonFin = QtWidgets.QPushButton(
            QtCore.QCoreApplication.translate("pychemqt", "Finned Pipe Database"))
        self.buttonFin.setEnabled(False)
        self.buttonFin.clicked.connect(self.showFinTube)
        lyt.addWidget(self.buttonFin, 11, 4, 1, 1)
        self.finned = QtWidgets.QComboBox()
        self.finned.addItem(QtCore.QCoreApplication.translate("pychemqt", "Bared tube"))
        self.finned.addItem(QtCore.QCoreApplication.translate("pychemqt", "Finned tube"))
        self.finned.currentIndexChanged.connect(self.finnedChanged)
        lyt.addWidget(self.finned, 11, 2)
        lyt.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Fouling")), 12, 1)
        self.tubeFouling = FoulingWidget()
        self.tubeFouling.valueChanged.connect(
            partial(self.changeParams, "tubeFouling"))
        lyt.addWidget(self.tubeFouling, 12, 2, 1, 3)
        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding),
            15, 1, 1, 6)

        # Shell tab
        tab = QtWidgets.QWidget()
        self.tabWidget.insertTab(
            3, tab, QtCore.QCoreApplication.translate("pychemqt", "Shell"))
        lyt = QtWidgets.QGridLayout(tab)
        lyt.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate(
            "pychemqt", "Exchangers in paralell")), 1, 1)
        self.paralelo = Entrada_con_unidades(int, width=60)
        self.paralelo.valueChanged.connect(
            partial(self.changeParams, "parallel"))
        lyt.addWidget(self.paralelo, 1, 2)
        lyt.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate(
            "pychemqt", "Exchangers in serie")), 2, 1)
        self.serie = Entrada_con_unidades(int, width=60)
        self.serie.valueChanged.connect(partial(self.changeParams, "serie"))
        lyt.addWidget(self.serie, 2, 2)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            3, 1, 1, 6)
        lyt.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Shell Diameter")), 4, 1)
        self.DShell = Entrada_con_unidades(Length)
        self.DShell.valueChanged.connect(partial(self.changeParams, "DShell"))
        lyt.addWidget(self.DShell, 4, 2)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            5, 1, 1, 6)
        lyt.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Shell Material")), 6, 1)
        self.materialShell = QtWidgets.QComboBox()
        lyt.addWidget(self.materialShell, 6, 2)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            7, 1, 1, 6)
        lyt.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Fouling")), 8, 1)
        self.shellFouling = FoulingWidget()
        self.shellFouling.valueChanged.connect(
            partial(self.changeParams, "shellFouling"))
        lyt.addWidget(self.shellFouling, 8, 2, 1, 2)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            9, 1, 1, 6)
        lyt.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Sealing Strips")), 10, 1)
        self.sealingStrips = Entrada_con_unidades(float)
        self.sealingStrips.valueChanged.connect(
            partial(self.changeParams, "sealingStrips"))
        lyt.addWidget(self.sealingStrips, 10, 2)

        group = QtWidgets.QGroupBox(
            QtCore.QCoreApplication.translate("pychemqt", "Clearances"))
        lyt.addWidget(group, 11, 1, 1, 6)
        lyt = QtWidgets.QGridLayout(group)
        lyt.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate(
            "pychemqt", "Tube to baffle hole")), 1, 1)
        self.ClearanceTubeBaffle = Entrada_con_unidades(Length, "Thickness")
        self.ClearanceTubeBaffle.valueChanged.connect(
            partial(self.changeParams, "clearanceTubeBaffle"))
        lyt.addWidget(self.ClearanceTubeBaffle, 1, 2)
        lyt.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Shell to baffle")), 2, 1)
        self.ClearanceShellBaffle = Entrada_con_unidades(Length, "Thickness")
        self.ClearanceShellBaffle.valueChanged.connect(
            partial(self.changeParams, "clearanceShellBaffle"))
        lyt.addWidget(self.ClearanceShellBaffle, 2, 2)
        lyt.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Shell to bundle")), 3, 1)
        self.ClearanceShellBundle = Entrada_con_unidades(Length, "Thickness")
        self.ClearanceShellBundle.valueChanged.connect(
            partial(self.changeParams, "clearanceShellBundle"))
        lyt.addWidget(self.ClearanceShellBundle, 3, 2)
        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding),
            15, 1, 1, 6)

        # Fitting tab
        tab = QtWidgets.QWidget()
        self.tabWidget.insertTab(
            4, tab, QtCore.QCoreApplication.translate("pychemqt", "Baffle"))
        lyt = QtWidgets.QGridLayout(tab)
        lyt.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Baffle type")), 1, 1)
        self.baffleType = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_BAFFLE_TYPE:
            self.baffleType.addItem(txt)
        self.baffleType.currentIndexChanged.connect(
            partial(self.changeParams, "baffleType"))
        lyt.addWidget(self.baffleType, 1, 2)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            2, 1, 1, 6)
        lyt.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Inlet spacing")), 3, 1)
        self.baffleSpacingIn = Entrada_con_unidades(Length)
        self.baffleSpacingIn.valueChanged.connect(
            partial(self.changeParams, "baffleSpacingIn"))
        lyt.addWidget(self.baffleSpacingIn, 3, 2)
        lyt.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Outlet spacing")), 3, 4)
        self.baffleSpacingOut = Entrada_con_unidades(Length)
        self.baffleSpacingOut.valueChanged.connect(
            partial(self.changeParams, "baffleSpacingOut"))
        lyt.addWidget(self.baffleSpacingOut, 3, 5)
        lyt.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Baffle spacing")), 4, 1)
        self.baffleSpacing = Entrada_con_unidades(Length)
        self.baffleSpacing.valueChanged.connect(
            partial(self.changeParams, "baffleSpacing"))
        lyt.addWidget(self.baffleSpacing, 4, 2)
        lyt.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate(
            "pychemqt", "Baffle thickness")), 5, 1)
        self.baffleThickness = Entrada_con_unidades(Length, "Thickness")
        self.baffleThickness.valueChanged.connect(
            partial(self.changeParams, "baffleThickness"))
        lyt.addWidget(self.baffleThickness, 5, 2)
        lyt.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Orientation")), 6, 1)
        self.baffleOrientation = QtWidgets.QComboBox()
        self.baffleOrientation.addItem(QtCore.QCoreApplication.translate("pychemqt", "Horizontal"))
        self.baffleOrientation.addItem(QtCore.QCoreApplication.translate("pychemqt", "Vertical"))
        self.baffleOrientation.currentIndexChanged.connect(
            partial(self.changeParams, "baffleOrientation"))
        lyt.addWidget(self.baffleOrientation, 6, 2)
        lyt.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Cut percent")), 7, 1)
        self.baffleCut = Entrada_con_unidades(float, textounidad="%")
        self.baffleCut.valueChanged.connect(
            partial(self.changeParams, "baffleCut"))
        lyt.addWidget(self.baffleCut, 7, 2)
        lyt.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Cut base")), 7, 4)
        self.baffleCutBase = QtWidgets.QComboBox()
        self.baffleCutBase.addItem(QtCore.QCoreApplication.translate("pychemqt", "Diameter"))
        self.baffleCutBase.addItem(QtCore.QCoreApplication.translate("pychemqt", "Area"))
        self.baffleCutBase.currentIndexChanged.connect(
            partial(self.changeParams, "baffleCutBase"))
        lyt.addWidget(self.baffleCutBase, 7, 5)
        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            1, 3, 6, 1)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            8, 1, 1, 6)

        group = QtWidgets.QGroupBox(
            QtCore.QCoreApplication.translate("pychemqt", "Nozzles"))
        lyt.addWidget(group, 9, 1, 1, 6)
        layout = QtWidgets.QGridLayout(group)
        layout.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Shellside")), 0, 2)
        layout.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Tubeside")), 0, 3)
        layout.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Input diameter")), 1, 1)
        self.nozzleInShellsideDiameter = Entrada_con_unidades(
            Length, "PipeDiameter")
        self.nozzleInShellsideDiameter.valueChanged.connect(
            partial(self.changeParams, "nozzleInShellsideDiameter"))
        layout.addWidget(self.nozzleInShellsideDiameter, 1, 2)
        self.nozzleInTubesideDiameter = Entrada_con_unidades(
            Length, "PipeDiameter")
        self.nozzleInTubesideDiameter.valueChanged.connect(
            partial(self.changeParams, "nozzleInTubesideDiameter"))
        layout.addWidget(self.nozzleInTubesideDiameter, 1, 3)
        layout.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Output diameter")), 2, 1)
        self.nozzleOutShellsideDiameter = Entrada_con_unidades(
            Length, "PipeDiameter")
        self.nozzleOutShellsideDiameter.valueChanged.connect(
            partial(self.changeParams, "nozzleOutShellsideDiameter"))
        layout.addWidget(self.nozzleOutShellsideDiameter, 2, 2)
        self.nozzleOutTubesideDiameter = Entrada_con_unidades(
            Length, "PipeDiameter")
        self.nozzleOutTubesideDiameter.valueChanged.connect(
            partial(self.changeParams, "nozzleOutTubesideDiameter"))
        layout.addWidget(self.nozzleOutTubesideDiameter, 2, 3)
        layout.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed),
            1, 4, 2, 1)

        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding),
            15, 1, 1, 6)

        # Calculate tab
        lyt = QtWidgets.QGridLayout(self.tabCalculo)
        lyt.addWidget(QtWidgets.QLabel(QtCore.QCoreApplication.translate(
            "pychemqt", "Calculation Mode")), 1, 1)
        self.modo = QtWidgets.QComboBox()
        self.modo.addItem(QtCore.QCoreApplication.translate("pychemqt", "Rating"))
        self.modo.addItem(QtCore.QCoreApplication.translate("pychemqt", "Design"))
        lyt.addWidget(self.modo, 1, 2)
        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding),
            15, 1, 1, 6)

        # Cost tab
        lyt = QtWidgets.QGridLayout(self.tabCostos)
        lyt.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Type")), 1, 1)
        self.tipoCoste = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_COST_TYPE:
            self.tipoCoste.addItem(txt)
        self.tipoCoste.currentIndexChanged.connect(
            partial(self.changeParamsCoste, "tipoCoste"))
        lyt.addWidget(self.tipoCoste, 1, 2)
        lyt.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Material")), 2, 1)
        self.materialCoste = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_COST_MATERIAL:
            self.materialCoste.addItem(txt)
        self.materialCoste.currentIndexChanged.connect(
            partial(self.changeParamsCoste, "P_dis"))
        lyt.addWidget(self.materialCoste, 2, 2)
        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding),
            3, 0, 1, 6)
        lyt.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Design Pressure")), 4, 1)
        self.Pdiseno = Entrada_con_unidades(Pressure)
        lyt.addWidget(self.Pdiseno, 4, 2, 1, 1)
        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding),
            5, 0, 1, 6)

        self.Costos = CostData(self.Equipment)
        self.Costos.valueChanged.connect(self.calcularCostos)
        lyt.addWidget(self.Costos, 6, 1, 2, 5)

        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding),
            8, 0, 1, 6)
        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding),
            10, 0, 1, 6)
        group = QtWidgets.QGroupBox(
            QtCore.QCoreApplication.translate("pychemqt", "Stimated Costs"))
        lyt.addWidget(group, 9, 1, 1, 5)
        layout = QtWidgets.QGridLayout(group)
        layout.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Purchase Cost")), 0, 1)
        self.C_adq = Entrada_con_unidades(
            Currency, retornar=False, readOnly=True)
        layout.addWidget(self.C_adq, 0, 2)
        layout.addWidget(QtWidgets.QLabel(
            QtCore.QCoreApplication.translate("pychemqt", "Installed Cost")), 1, 1)
        self.C_inst = Entrada_con_unidades(
            Currency, retornar=False, readOnly=True)
        self.C_inst.entrada.setReadOnly(True)
        layout.addWidget(self.C_inst, 1, 2)

        # Output Tab
        self.addSalida(QtCore.QCoreApplication.translate("pychemqt", "Tubes"))
        self.addSalida(QtCore.QCoreApplication.translate("pychemqt", "Shell"))

        if equipment:
            self.setEquipment(equipment)

    def selectMethods(self):
        """Show dialog for select calculation methods"""
        dialogo = Dialog_Methods(self.Equipment)
        if dialogo.exec_():
            self.Equipment(**dialogo.kwargs)

    def showMaterial(self):
        dialogo = Catalogo_Materiales_Dialog()
        if dialogo.exec_():
            pass

    def rellenarFouling(self, widget, txt):
        if txt:
            widget.setReadOnly(True)
            widget.setValue(Fouling_Factor_Shell_Tube_Exchanger[str(txt)])
        else:
            widget.setReadOnly(False)

    def finnedChanged(self, ind):
        self.buttonFin.setEnabled(ind)
        self.changeParams("finned", ind)

    def showFinTube(self):
        dialogo = Dialog_Finned()
        if dialogo.exec_():
            pass

    def rellenar(self):
        pass

    def calcularCostos(self):
        if self.todos_datos():
            self.Equipment.Coste(self.factorInstalacion.value(), 0, self.tipo.currentIndex(), self.material.currentIndex())
            self.C_adq.setValue(self.Equipment.C_adq.config())
            self.C_inst.setValue(self.Equipment.C_inst.config())


if __name__ == "__main__":
    import sys
    from lib.corriente import Corriente, Mezcla
    app = QtWidgets.QApplication(sys.argv)
#    aguaFria=Corriente(300, 1, 100, Mezcla([62], [1]))
#    aguaCaliente=Corriente(370, 1, 250, Mezcla([62], [1]))
    dialogo = UI_equipment()
    dialogo.show()
    sys.exit(app.exec_())
