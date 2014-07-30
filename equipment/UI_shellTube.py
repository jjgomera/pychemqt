#!/usr/bin/python
# -*- coding: utf-8 -*-

#######################################################################
###                Diálogo de definición de cambiadores de calor de carcasa y tubos, UI_shellTube                           ###
#######################################################################

from functools import partial

from PyQt4 import QtCore, QtGui

from lib import unidades
from UI import UI_corriente
from UI.widgets import Entrada_con_unidades
from equipment.parents import UI_equip, FoulingWidget, Dialog_Finned
from equipment.UI_pipe import Catalogo_Materiales_Dialog
from equipment.heatExchanger import Shell_Tube
from tools.costIndex import CostData



class Dialog_Methods(QtGui.QDialog):
    def __init__(self, equipment, parent=None):
        super(Dialog_Methods, self).__init__(parent=parent)
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Specify calculation methods"))
        layout = QtGui.QGridLayout(self)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Tubeside laminar flow")),1,1)
        self.tubesideLaminar=QtGui.QComboBox()
        for txt in equipment.TEXT_METHOD_TUBE_LAMINAR:
            self.tubesideLaminar.addItem(txt)
        layout.addWidget(self.tubesideLaminar,1,2)

        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Tubeside turbulent flow")),2,1)
        self.tubesideTurbulent=QtGui.QComboBox()
        for txt in equipment.TEXT_METHOD_TUBE_TURBULENT:
            self.tubesideTurbulent.addItem(txt)
        layout.addWidget(self.tubesideTurbulent,2,2)

        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "ShellSide")),3,1)
        self.shellSide=QtGui.QComboBox()
        for txt in equipment.TEXT_METHOD_SHELL:
            self.shellSide.addItem(txt)
        layout.addWidget(self.shellSide,3,2)

        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox,10,1,1,2)

        if equipment.kwargs:
            self.tubesideLaminar.setCurrentIndex(equipment.kwargs["tubesideLaminar"])
            self.tubesideTurbulent.setCurrentIndex(equipment.kwargs["tubesideTurbulent"])
            self.shellSide.setCurrentIndex(equipment.kwargs["shellsideSensible"])

    @property
    def kwargs(self):
        kwargs={}
        kwargs["tubesideLaminar"]=self.tubesideLaminar.currentIndex()
        kwargs["tubesideTurbulent"]=self.tubesideTurbulent.currentIndex()
        kwargs["shellsideSensible"]=self.shellSide.currentIndex()
        return kwargs


        

class UI_equipment(UI_equip):
    """Diálogo de definición de cambiadores de calor de carcasa y tubos"""
    Equipment=Shell_Tube()
    def __init__(self, equipment=None, parent=None):
        """
        equipment: instancia de equipo inicial
        """
        super(UI_equipment, self).__init__(Shell_Tube, parent=parent)
        
        #Pestaña entrada
        self.EntradaTubo= UI_corriente.Ui_corriente()
        self.EntradaTubo.Changed.connect(partial(self.changeParams, "entradaTubo"))
        self.entrada.addTab(self.EntradaTubo,QtGui.QApplication.translate("pychemqt", "Tubes"))
        self.EntradaCarcasa= UI_corriente.Ui_corriente()
        self.EntradaCarcasa.Changed.connect(partial(self.changeParams, "entradaCarcasa"))
        self.entrada.addTab(self.EntradaCarcasa, QtGui.QApplication.translate("pychemqt", "Shell"))
        
        #Pestaña modelo
        self.tabModelo = QtGui.QWidget()
        self.tabWidget.insertTab(1, self.tabModelo,QtGui.QApplication.translate("pychemqt", "Model"))
        gridLayout_Modelo = QtGui.QGridLayout(self.tabModelo)
        gridLayout_Modelo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Class")),2,1)
        self.class_=QtGui.QComboBox()
        for txt in self.Equipment.TEXT_CLASS:
            self.class_.addItem(txt)
        self.class_.currentIndexChanged.connect(partial(self.changeParams, "class_"))
        gridLayout_Modelo.addWidget(self.class_,2,2)
        gridLayout_Modelo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Front end head")),3,1)
        self.frontHead=QtGui.QComboBox()
        for txt in self.Equipment.TEXT_FRONTHEAD:
            self.frontHead.addItem(txt)
        self.frontHead.currentIndexChanged.connect(partial(self.changeParams, "frontHead"))
        gridLayout_Modelo.addWidget(self.frontHead,3,2,1,3)
        gridLayout_Modelo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Shell type")),4,1)
        self.shell=QtGui.QComboBox()
        for txt in self.Equipment.TEXT_SHELL:
            self.shell.addItem(txt)
        self.shell.currentIndexChanged.connect(partial(self.changeParams, "shell"))
        gridLayout_Modelo.addWidget(self.shell,4,2)
        gridLayout_Modelo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Rear end head")),5,1)
        self.rearHead=QtGui.QComboBox()
        for txt in self.Equipment.TEXT_REARHEAD:
            self.rearHead.addItem(txt)
        
        self.rearHead.currentIndexChanged.connect(partial(self.changeParams, "rearHead"))
        gridLayout_Modelo.addWidget(self.rearHead,5,2,1,2)
        gridLayout_Modelo.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),6,1,1,6)
        gridLayout_Modelo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Orientation")),7,1)
        self.orientacion=QtGui.QComboBox()
        for txt in self.Equipment.TEXT_ORIENTATION:
            self.orientacion.addItem(txt)
        self.orientacion.currentIndexChanged.connect(partial(self.changeParams, "orientation"))
        gridLayout_Modelo.addWidget(self.orientacion,7,2)
        
        botonMetodos=QtGui.QPushButton(QtGui.QApplication.translate("pychemqt", "Calculation methods"))
        botonMetodos.clicked.connect(self.selectMethods)
        gridLayout_Modelo.addWidget(botonMetodos,9,1)
        gridLayout_Modelo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),10,1,1,6)
        
        
        #Pestaña Tubos
        self.tabTube = QtGui.QWidget()
        self.tabWidget.insertTab(2, self.tabTube,QtGui.QApplication.translate("pychemqt", "Tubes"))
        gridLayout_Tube = QtGui.QGridLayout(self.tabTube)
        gridLayout_Tube.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Number of tubes")),1,1)
        self.NTubes=Entrada_con_unidades(int, width=60, spinbox=True, step=1)
        self.NTubes.valueChanged.connect(partial(self.changeParams, "NTube"))
        gridLayout_Tube.addWidget(self.NTubes,1,2)
        gridLayout_Tube.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Number of tube passes")),2,1)
        self.NPases=Entrada_con_unidades(int, width=60, spinbox=True, step=1)
        self.NPases.valueChanged.connect(partial(self.changeParams, "NPases"))
        gridLayout_Tube.addWidget(self.NPases,2,2)
        gridLayout_Tube.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Tube length")),3,1)
        self.LTube=Entrada_con_unidades(unidades.Length)
        self.LTube.valueChanged.connect(partial(self.changeParams, "LTube"))
        gridLayout_Tube.addWidget(self.LTube,3,2)
        gridLayout_Tube.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Tube external diameter")),4,1)
        self.DeTube=Entrada_con_unidades(unidades.Length, "pipeDiameter")
        self.DeTube.valueChanged.connect(partial(self.changeParams, "DeTube"))
        gridLayout_Tube.addWidget(self.DeTube,4,2)
        gridLayout_Tube.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Tube thickness")),5,1)
        self.wTube=Entrada_con_unidades(unidades.Length, "Thickness")
        self.wTube.valueChanged.connect(partial(self.changeParams, "wTube"))
        gridLayout_Tube.addWidget(self.wTube,5,2)
        gridLayout_Tube.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Tube roughness")),6,1)
        self.rTube=Entrada_con_unidades(unidades.Length, "Thickness")
        self.rTube.valueChanged.connect(partial(self.changeParams, "rTube"))
        gridLayout_Tube.addWidget(self.rTube,6,2)
        gridLayout_Tube.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Thermal conductivity")),7,1)
        self.kTube=Entrada_con_unidades(unidades.ThermalConductivity)
        self.kTube.valueChanged.connect(partial(self.changeParams, "kTube"))
        gridLayout_Tube.addWidget(self.kTube,7,2)
        self.buttonPipe=QtGui.QPushButton(QtGui.QApplication.translate("pychemqt", "Pipe Database"))
        self.buttonPipe.clicked.connect(self.showMaterial)
        gridLayout_Tube.addWidget(self.buttonPipe,4,4,4,1)
        gridLayout_Tube.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),8,1,1,6)

        gridLayout_Tube.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Tube pattern")),9,1)
        self.distribucionTube=QtGui.QComboBox()
        self.distribucionTube.addItem(QtGui.QApplication.translate("pychemqt", "Triangular")+u", 30º")
        self.distribucionTube.addItem(QtGui.QApplication.translate("pychemqt", "Diamond")+u", 45º")
        self.distribucionTube.addItem(QtGui.QApplication.translate("pychemqt", "Rotated Triangular")+u", 60º")
        self.distribucionTube.addItem(QtGui.QApplication.translate("pychemqt", "Square")+u", 90º")
        self.distribucionTube.currentIndexChanged.connect(partial(self.changeParams, "distribucionTube"))
        gridLayout_Tube.addWidget(self.distribucionTube,9,2)
        gridLayout_Tube.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Tube pitch")),10,1)
        self.pitch=Entrada_con_unidades(unidades.Length)
        self.pitch.valueChanged.connect(partial(self.changeParams, "pitch"))
        gridLayout_Tube.addWidget(self.pitch,10,2)
        gridLayout_Tube.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Fin Tube")),11,1)
        self.buttonFin=QtGui.QPushButton(QtGui.QApplication.translate("pychemqt", "Finned Pipe Database"))
        self.buttonFin.setEnabled(False)
        self.buttonFin.clicked.connect(self.showFinTube)
        gridLayout_Tube.addWidget(self.buttonFin,11,4,1,1)
        self.finned=QtGui.QComboBox()
        self.finned.addItem(QtGui.QApplication.translate("pychemqt", "Bared tube"))
        self.finned.addItem(QtGui.QApplication.translate("pychemqt", "Finned tube"))
        self.finned.currentIndexChanged.connect(self.finnedChanged)
        gridLayout_Tube.addWidget(self.finned,11,2)
        gridLayout_Tube.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Fouling")),12,1)
        self.tubeFouling=FoulingWidget()
        self.tubeFouling.valueChanged.connect(partial(self.changeParams, "tubeFouling"))
        gridLayout_Tube.addWidget(self.tubeFouling,12,2,1,3)
        gridLayout_Tube.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),15,1,1,6)


        #Pestaña Carcasa
        self.tabShell = QtGui.QWidget()
        self.tabWidget.insertTab(3, self.tabShell,QtGui.QApplication.translate("pychemqt", "Shell"))
        gridLayout_Shell = QtGui.QGridLayout(self.tabShell)
        gridLayout_Shell.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Exchangers in paralell")),1,1)
        self.paralelo=Entrada_con_unidades(int, width=60)
        self.paralelo.valueChanged.connect(partial(self.changeParams, "parallel"))
        gridLayout_Shell.addWidget(self.paralelo,1,2)
        gridLayout_Shell.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Exchangers in serie")),2,1)
        self.serie=Entrada_con_unidades(int, width=60)
        self.serie.valueChanged.connect(partial(self.changeParams, "serie"))
        gridLayout_Shell.addWidget(self.serie,2,2)
        gridLayout_Shell.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),3,1,1,6)
        gridLayout_Shell.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Shell Diameter")),4,1)
        self.DShell=Entrada_con_unidades(unidades.Length)
        self.DShell.valueChanged.connect(partial(self.changeParams, "DShell"))
        gridLayout_Shell.addWidget(self.DShell,4,2)
        gridLayout_Shell.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),5,1,1,6)
        gridLayout_Shell.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Shell Material")),6,1)
        self.materialShell=QtGui.QComboBox()
        gridLayout_Shell.addWidget(self.materialShell,6,2)
        gridLayout_Shell.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),7,1,1,6)
        gridLayout_Shell.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Fouling")),8,1)
        self.shellFouling=FoulingWidget()
        self.shellFouling.valueChanged.connect(partial(self.changeParams, "shellFouling"))
        gridLayout_Shell.addWidget(self.shellFouling,8,2,1,2)
        gridLayout_Shell.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),9,1,1,6)    
        gridLayout_Shell.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Sealing Strips")),10,1)
        self.sealingStrips=Entrada_con_unidades(float)
        self.sealingStrips.valueChanged.connect(partial(self.changeParams, "sealingStrips"))
        gridLayout_Shell.addWidget(self.sealingStrips,10,2)
        
        group=QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Clearances"))
        gridLayout_Shell.addWidget(group,11,1,1,6)    
        lyt=QtGui.QGridLayout(group)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Tube to baffle hole")),1,1)
        self.ClearanceTubeBaffle=Entrada_con_unidades(unidades.Length, "Thickness")
        self.ClearanceTubeBaffle.valueChanged.connect(partial(self.changeParams, "clearanceTubeBaffle"))
        lyt.addWidget(self.ClearanceTubeBaffle,1,2)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Shell to baffle")),2,1)
        self.ClearanceShellBaffle=Entrada_con_unidades(unidades.Length, "Thickness")
        self.ClearanceShellBaffle.valueChanged.connect(partial(self.changeParams, "clearanceShellBaffle"))
        lyt.addWidget(self.ClearanceShellBaffle,2,2)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Shell to bundle")),3,1)
        self.ClearanceShellBundle=Entrada_con_unidades(unidades.Length, "Thickness")
        self.ClearanceShellBundle.valueChanged.connect(partial(self.changeParams, "clearanceShellBundle"))
        lyt.addWidget(self.ClearanceShellBundle,3,2)

        gridLayout_Shell.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),15,1,1,6)

        #Pestaña Accesories
        self.tabFitting = QtGui.QWidget()
        self.tabWidget.insertTab(4, self.tabFitting,QtGui.QApplication.translate("pychemqt", "Baffle"))
        gridLayout_Fitting = QtGui.QGridLayout(self.tabFitting)
        gridLayout_Fitting.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Baffle type")),1,1)
        self.baffleType=QtGui.QComboBox()
        self.baffleType.addItem(QtGui.QApplication.translate("pychemqt", "Single segmental"))
        self.baffleType.addItem(QtGui.QApplication.translate("pychemqt", "Double segmental"))
        self.baffleType.addItem(QtGui.QApplication.translate("pychemqt", "Triple segmental"))
        self.baffleType.addItem(QtGui.QApplication.translate("pychemqt", "No tubes in window"))
        self.baffleType.addItem(QtGui.QApplication.translate("pychemqt", "Disk & donut"))
        self.baffleType.addItem(QtGui.QApplication.translate("pychemqt", "Rod"))
        self.baffleType.currentIndexChanged.connect(partial(self.changeParams, "typeBaffle"))
        gridLayout_Fitting.addWidget(self.baffleType,1,2)     
        gridLayout_Fitting.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),2,1,1,6)
        gridLayout_Fitting.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Inlet spacing")),3,1)
        self.baffleSpacingIn=Entrada_con_unidades(unidades.Length)
        self.baffleSpacingIn.valueChanged.connect(partial(self.changeParams, "baffleSpacingIn"))
        gridLayout_Fitting.addWidget(self.baffleSpacingIn,3,2)
        gridLayout_Fitting.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Outlet spacing")),3,4)
        self.baffleSpacingOut=Entrada_con_unidades(unidades.Length)
        self.baffleSpacingOut.valueChanged.connect(partial(self.changeParams, "baffleSpacingOut"))
        gridLayout_Fitting.addWidget(self.baffleSpacingOut,3,5)
        gridLayout_Fitting.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Baffle spacing")),4,1)
        self.baffleSpacing=Entrada_con_unidades(unidades.Length)
        self.baffleSpacing.valueChanged.connect(partial(self.changeParams, "baffleSpacing"))
        gridLayout_Fitting.addWidget(self.baffleSpacing,4,2)
        gridLayout_Fitting.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Baffle thickness")),5,1)
        self.baffleThickness=Entrada_con_unidades(unidades.Length, "Thickness")
        self.baffleThickness.valueChanged.connect(partial(self.changeParams, "baffleThickness"))
        gridLayout_Fitting.addWidget(self.baffleThickness,5,2)
        gridLayout_Fitting.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Orientation")),6,1)
        self.baffleOrientation=QtGui.QComboBox()
        self.baffleOrientation.addItem(QtGui.QApplication.translate("pychemqt", "Horizontal"))
        self.baffleOrientation.addItem(QtGui.QApplication.translate("pychemqt", "Vertical"))
        self.baffleOrientation.currentIndexChanged.connect(partial(self.changeParams, "baffleOrientation"))
        gridLayout_Fitting.addWidget(self.baffleOrientation,6,2)
        gridLayout_Fitting.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Cut percent")),7,1)
        self.baffleCut=Entrada_con_unidades(float, textounidad="%")
        self.baffleCut.valueChanged.connect(partial(self.changeParams, "baffleCut"))
        gridLayout_Fitting.addWidget(self.baffleCut,7,2)
        gridLayout_Fitting.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Cut base")),7,4)
        self.baffleCutBase=QtGui.QComboBox()
        self.baffleCutBase.addItem(QtGui.QApplication.translate("pychemqt", "Diameter"))
        self.baffleCutBase.addItem(QtGui.QApplication.translate("pychemqt", "Area"))
        self.baffleCutBase.currentIndexChanged.connect(partial(self.changeParams, "baffleCutBase"))
        gridLayout_Fitting.addWidget(self.baffleCutBase,7,5)
        gridLayout_Fitting.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),1,3,6,1)
        gridLayout_Fitting.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),8,1,1,6)

        group=QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Nozzles"))
        gridLayout_Fitting.addWidget(group,9,1,1,6)    
        lyt=QtGui.QGridLayout(group)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Shellside")),0,2)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Tubeside")),0,3)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Input diameter")),1,1)
        self.nozzleInShellsideDiameter=Entrada_con_unidades(unidades.Length, "PipeDiameter")
        self.nozzleInShellsideDiameter.valueChanged.connect(partial(self.changeParams, "nozzleInShellsideDiameter"))
        lyt.addWidget(self.nozzleInShellsideDiameter,1,2)
        self.nozzleInTubesideDiameter=Entrada_con_unidades(unidades.Length, "PipeDiameter")
        self.nozzleInTubesideDiameter.valueChanged.connect(partial(self.changeParams, "nozzleInTubesideDiameter"))
        lyt.addWidget(self.nozzleInTubesideDiameter,1,3)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Output diameter")),2,1)
        self.nozzleOutShellsideDiameter=Entrada_con_unidades(unidades.Length, "PipeDiameter")
        self.nozzleOutShellsideDiameter.valueChanged.connect(partial(self.changeParams, "nozzleOutShellsideDiameter"))
        lyt.addWidget(self.nozzleOutShellsideDiameter,2,2)
        self.nozzleOutTubesideDiameter=Entrada_con_unidades(unidades.Length, "PipeDiameter")
        self.nozzleOutTubesideDiameter.valueChanged.connect(partial(self.changeParams, "nozzleOutTubesideDiameter"))
        lyt.addWidget(self.nozzleOutTubesideDiameter,2,3)
        lyt.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Fixed),1,4,2,1)

        gridLayout_Fitting.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),15,1,1,6)
        
        #Pestaña calculo
        gridLayout_Calculo = QtGui.QGridLayout(self.tabCalculo)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Calculation Mode")),1,1)
        self.modo=QtGui.QComboBox()
        self.modo.addItem(QtGui.QApplication.translate("pychemqt", "Rating"))
        self.modo.addItem(QtGui.QApplication.translate("pychemqt", "Design"))
        gridLayout_Calculo.addWidget(self.modo,1,2)

        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),15,1,1,6)


        #Pestaña costos
        gridLayout_Costos = QtGui.QGridLayout(self.tabCostos)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Type")),1,1)
        self.tipo=QtGui.QComboBox()
        self.tipo.addItem(QtGui.QApplication.translate("pychemqt", "Fixed Head"))
        self.tipo.addItem(QtGui.QApplication.translate("pychemqt", "Kettle Reboiler"))
        self.tipo.addItem(QtGui.QApplication.translate("pychemqt", "U-Tube"))
        self.tipo.currentIndexChanged.connect(self.calcularCostos)
        gridLayout_Costos.addWidget(self.tipo,1,2)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Material")),2,1)
        self.material=QtGui.QComboBox()
        self.material.addItem(QtGui.QApplication.translate("pychemqt", "Carbon Steel"))
        self.material.addItem(QtGui.QApplication.translate("pychemqt", "Stainless Steel 316"))
        self.material.addItem(QtGui.QApplication.translate("pychemqt", "Stainless Steel 304"))
        self.material.addItem(QtGui.QApplication.translate("pychemqt", "Stainless Steel 347"))
        self.material.addItem(QtGui.QApplication.translate("pychemqt", "Nickel 200"))
        self.material.addItem(QtGui.QApplication.translate("pychemqt", "Monel 400"))
        self.material.addItem(QtGui.QApplication.translate("pychemqt", "Inconel 600"))
        self.material.addItem(QtGui.QApplication.translate("pychemqt", "Incoloy 825"))
        self.material.addItem(QtGui.QApplication.translate("pychemqt", "Titanium"))
        self.material.addItem(QtGui.QApplication.translate("pychemqt", "Hastelloy"))
        self.material.currentIndexChanged.connect(self.calcularCostos)
        gridLayout_Costos.addWidget(self.material,2,2)
        gridLayout_Costos.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),3,0,1,6)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Design Pressure")),4,1)
        self.Pdiseno=Entrada_con_unidades(unidades.Pressure)
        gridLayout_Costos.addWidget(self.Pdiseno,4,2,1,1)
        gridLayout_Costos.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),5,0,1,6)

        self.Costos=CostData(self.Equipment)
        self.Costos.valueChanged.connect(self.calcularCostos)
        gridLayout_Costos.addWidget(self.Costos,6,1,2,5)

        gridLayout_Costos.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),8,0,1,6)
        gridLayout_Costos.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),10,0,1,6)
        self.groupBox_Costos = QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Stimated Costs"))
        gridLayout_Costos.addWidget(self.groupBox_Costos,9,1,1,5)
        gridLayout_5 = QtGui.QGridLayout(self.groupBox_Costos)
        gridLayout_5.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Purchase Cost")),0,1)
        self.C_adq=Entrada_con_unidades(unidades.Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_adq,0,2)
        gridLayout_5.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Installed Cost")),1,1)
        self.C_inst=Entrada_con_unidades(unidades.Currency, retornar=False, readOnly=True)
        self.C_inst.entrada.setReadOnly(True)
        gridLayout_5.addWidget(self.C_inst,1,2)

        #Pestaña salida
        self.SalidaTubo= UI_corriente.Ui_corriente(readOnly=True)
        self.Salida.addTab(self.SalidaTubo, QtGui.QApplication.translate("pychemqt", "Tubes"))
        self.SalidaCarcasa= UI_corriente.Ui_corriente(readOnly=True)
        self.Salida.addTab(self.SalidaCarcasa, QtGui.QApplication.translate("pychemqt", "Shell"))
        
        if equipment:
            self.setEquipment(equipment)


    def selectMethods(self):
        dialogo=Dialog_Methods(self.Equipment)
        if dialogo.exec_():
            self.Equipment(**dialogo.kwargs)
        
    def showMaterial(self):
        dialogo=Catalogo_Materiales_Dialog()
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
        dialogo=Dialog_Finned()
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
    app = QtGui.QApplication(sys.argv)
#    aguaFria=Corriente(300, 1, 100, Mezcla([62], [1]))
#    aguaCaliente=Corriente(370, 1, 250, Mezcla([62], [1]))
    dialogo = UI_equipment()
    dialogo.show()
    sys.exit(app.exec_())
