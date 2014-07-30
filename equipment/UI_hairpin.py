#!/usr/bin/python
# -*- coding: utf-8 -*-

#####################################################
###  Diálogo de definición de cambiadores de calor de doble tubo, UI_doublePipe  ###
#####################################################

from functools import partial

from PyQt4 import QtCore, QtGui

from lib.unidades import  Temperature, DeltaT, Pressure, DeltaP, Power, Length, Area, ThermalConductivity, HeatTransfCoef, Currency
from UI.widgets import Entrada_con_unidades
from UI import UI_corriente
from equipment.heatExchanger import Hairpin
from equipment.UI_pipe import Catalogo_Materiales_Dialog
from equipment.parents import UI_equip, FoulingWidget, Dialog_Finned
from tools.costIndex import CostData





class UI_equipment(UI_equip):
    """Diálogo de definición de cambiadores de calor generales"""
    Equipment=Hairpin()
    def __init__(self, equipment=None, parent=None):
        """
        equipment: instancia de equipo inicial
        """
        super(UI_equipment, self).__init__(Hairpin, parent=parent)
        
        #Pestaña entrada
        self.entradaTubo= UI_corriente.Ui_corriente()
        self.entradaTubo.Changed.connect(partial(self.changeParams, "entradaTubo"))
        self.entrada.addTab(self.entradaTubo,QtGui.QApplication.translate("pychemqt", "Tube"))
        self.entradaExterior= UI_corriente.Ui_corriente()
        self.entradaExterior.Changed.connect(partial(self.changeParams, "entradaExterior"))
        self.entrada.addTab(self.entradaExterior, QtGui.QApplication.translate("pychemqt", "Annulli"))

        #Pestaña Catalogo
        tabCatalogo = QtGui.QWidget()
        self.tabWidget.insertTab(1, tabCatalogo,QtGui.QApplication.translate("pychemqt", "Catalog"))
        gridLayout_Catalogo = QtGui.QGridLayout(tabCatalogo)
        gridLayout_Catalogo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Tube length")),4,1)
        self.LTube=Entrada_con_unidades(Length)
        self.LTube.valueChanged.connect(partial(self.changeParams, "LTube"))
        gridLayout_Catalogo.addWidget(self.LTube,4,2)
        gridLayout_Catalogo.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),5,1)
        gridLayout_Catalogo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Tube internal diameter")),6,1)
        self.DiTube=Entrada_con_unidades(Length, "pipeDiameter")
        self.DiTube.valueChanged.connect(partial(self.changeParams, "DiTube"))
        gridLayout_Catalogo.addWidget(self.DiTube,6,2)
        gridLayout_Catalogo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Tube external diameter")),7,1)
        self.DeTube=Entrada_con_unidades(Length, "pipeDiameter")
        self.DeTube.valueChanged.connect(partial(self.changeParams, "DeTube"))
        gridLayout_Catalogo.addWidget(self.DeTube,7,2)
        gridLayout_Catalogo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Tube thickness")),8,1)
        self.wTube=Entrada_con_unidades(Length, "Thickness")
        self.wTube.valueChanged.connect(partial(self.changeParams, "wTube"))
        gridLayout_Catalogo.addWidget(self.wTube,8,2)
        gridLayout_Catalogo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Tube roughness")),9,1)
        self.rTube=Entrada_con_unidades(Length, "Thickness")
        self.rTube.valueChanged.connect(partial(self.changeParams, "rTube"))
        gridLayout_Catalogo.addWidget(self.rTube,9,2)
        gridLayout_Catalogo.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),10,1)
        gridLayout_Catalogo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Annulli external diameter")),11,1)
        self.DeeTube=Entrada_con_unidades(Length, "pipeDiameter")
        self.DeeTube.valueChanged.connect(partial(self.changeParams, "DeeTube"))
        gridLayout_Catalogo.addWidget(self.DeeTube,11,2)
        gridLayout_Catalogo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Thermal conductivity")),12,1)
        self.kTube=Entrada_con_unidades(ThermalConductivity)
        self.kTube.valueChanged.connect(partial(self.changeParams, "kTube"))
        gridLayout_Catalogo.addWidget(self.kTube,12,2)
        gridLayout_Catalogo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Tube Count")),13,1)
        self.nTube=Entrada_con_unidades(int)
        self.nTube.valueChanged.connect(partial(self.changeParams, "nTube"))
        gridLayout_Catalogo.addWidget(self.nTube,13,2)

        self.buttonPipe=QtGui.QPushButton(QtGui.QApplication.translate("pychemqt", "Pipe Database"))
        self.buttonPipe.clicked.connect(self.showMaterial)
        gridLayout_Catalogo.addWidget(self.buttonPipe,6,3,4,1)
        gridLayout_Catalogo.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),14,1)
        self.tubeFinned=QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Finned Tube"))
        gridLayout_Catalogo.addWidget(self.tubeFinned,15,1,1,4)
        self.buttonFin=QtGui.QPushButton(QtGui.QApplication.translate("pychemqt", "Finned Pipe Database"))
        self.buttonFin.setEnabled(False)
        self.buttonFin.clicked.connect(self.showFinTube)
        gridLayout_Catalogo.addWidget(self.buttonFin,15,3,1,1)
        self.tubeFinned.toggled.connect(partial(self.changeParams, "tubeFinned"))
        self.tubeFinned.toggled.connect(self.buttonFin.setEnabled)
        gridLayout_Catalogo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Inside Fouling")),16,1)
        self.tubeFouling=FoulingWidget()
        self.tubeFouling.valueChanged.connect(partial(self.changeParams, "tubeFouling"))
        gridLayout_Catalogo.addWidget(self.tubeFouling,16,2,1,5)
        gridLayout_Catalogo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Outside Fouling")),17,1)
        self.annulliFouling=FoulingWidget()
        self.annulliFouling.valueChanged.connect(partial(self.changeParams, "annulliFouling"))
        gridLayout_Catalogo.addWidget(self.annulliFouling,17,2,1,5)
        gridLayout_Catalogo.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),20,1,1,6)
        
        #Pestaña calculo
        gridLayout_Calculo = QtGui.QGridLayout(self.tabCalculo)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Mode")),1,1)
        self.modo= QtGui.QComboBox()
        for txt in self.Equipment.TEXT_MODO:
            self.modo.addItem(txt)
        self.modo.currentIndexChanged.connect(partial(self.changeParams, "modo"))
        gridLayout_Calculo.addWidget(self.modo,1,2)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Flujo")),2,1)
        self.flujo= QtGui.QComboBox()
        for txt in self.Equipment.TEXT_FLUJO:
            self.flujo.addItem(txt)
        self.flujo.currentIndexChanged.connect(partial(self.changeParams, "flujo"))
        gridLayout_Calculo.addWidget(self.flujo,2,2)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Layout")),3,1)
        self.orientacion= QtGui.QComboBox()
        for txt in self.Equipment.TEXT_ORIENTACION:
            self.orientacion.addItem(txt)
        self.orientacion.currentIndexChanged.connect(partial(self.changeParams, "orientacion"))
        gridLayout_Calculo.addWidget(self.orientacion,3,2)
        
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),4,1)

        
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Output inside temperature")),5,1)
        self.tubeTout=Entrada_con_unidades(Temperature)
        self.tubeTout.valueChanged.connect(partial(self.changeParams, "tubeTout"))
        gridLayout_Calculo.addWidget(self.tubeTout,5,2)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Output annulli temperature")),6,1)
        self.annulliTout=Entrada_con_unidades(Temperature)
        self.annulliTout.valueChanged.connect(partial(self.changeParams, "annulliTout"))
        gridLayout_Calculo.addWidget(self.annulliTout,6,2)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Output inside quality")),5,4)
        self.tubeXout=Entrada_con_unidades(float)
        self.tubeXout.valueChanged.connect(partial(self.changeParams, "tubeXout"))
        gridLayout_Calculo.addWidget(self.tubeXout,5,5)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Output annulli quality")),6,4)
        self.annulliXout=Entrada_con_unidades(float)
        self.annulliXout.valueChanged.connect(partial(self.changeParams, "annulliXout"))
        gridLayout_Calculo.addWidget(self.annulliXout,6,5)

        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),15,1,1,6)
        
        self.groupBox_Calculo = QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Results"))
        gridLayout_Calculo.addWidget(self.groupBox_Calculo,16,1,1,6)
        gridLayout_1 = QtGui.QGridLayout(self.groupBox_Calculo)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Heat Duty")),0,1)
        self.Q=Entrada_con_unidades(Power, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.Q,0,2)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Tout Tube")),1,1)
        self.ToutTube=Entrada_con_unidades(Temperature, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.ToutTube,1,2)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Tout Tube")),2,1)
        self.ToutAnnulli=Entrada_con_unidades(Temperature, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.ToutAnnulli,2,2)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "U")),0,4)
        self.U=Entrada_con_unidades(HeatTransfCoef, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.U,0,5)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Area")),1,4)
        self.A=Entrada_con_unidades(Area, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.A,1,5)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Lenght")),2,4)
        self.L=Entrada_con_unidades(Length, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.L,2,5)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "DeltaP Tube")),0,7)
        self.deltaPTube=Entrada_con_unidades(DeltaP, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.deltaPTube,0,8)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "DeltaP Annulli")),1,7)
        self.deltaPAnnulli=Entrada_con_unidades(DeltaP, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.deltaPAnnulli,1,8)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "CF")),2,7)
        self.CF=Entrada_con_unidades(float, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.CF,2,8)

        gridLayout_Calculo.addItem(QtGui.QSpacerItem(0,0,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),17,1,1,6)


        #Pestaña costos
        gridLayout_Costos = QtGui.QGridLayout(self.tabCostos)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Material")),2,1)
        self.material=QtGui.QComboBox()
        self.material.addItem(QtGui.QApplication.translate("pychemqt", "Carbon steel/carbon steel"))
        self.material.addItem(QtGui.QApplication.translate("pychemqt", "Carbon steel/304 stainless"))
        self.material.addItem(QtGui.QApplication.translate("pychemqt", "Carbon steel/316 stainless"))
        self.material.currentIndexChanged.connect(partial(self.changeParamsCoste, "material"))
        gridLayout_Costos.addWidget(self.material,2,2)
        gridLayout_Costos.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),3,0,1,6)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Design Pressure")),4,1)
        self.P_dis=Entrada_con_unidades(Pressure)
        self.P_dis.valueChanged.connect(partial(self.changeParamsCoste, "P_dis"))
        gridLayout_Costos.addWidget(self.P_dis,4,2,1,1)
        gridLayout_Costos.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),5,0,1,6)

        self.Costos=CostData(self.Equipment)
        self.Costos.valueChanged.connect(self.changeParamsCoste)
        gridLayout_Costos.addWidget(self.Costos,6,1,2,5)

        gridLayout_Costos.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),8,0,1,6)
        gridLayout_Costos.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),10,0,1,6)
        self.groupBox_Costos = QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Stimated Costs"))
        gridLayout_Costos.addWidget(self.groupBox_Costos,9,1,1,5)
        gridLayout_5 = QtGui.QGridLayout(self.groupBox_Costos)
        gridLayout_5.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Purchase Cost")),0,1)
        self.C_adq=Entrada_con_unidades(Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_adq,0,2)
        gridLayout_5.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Installed Cost")),1,1)
        self.C_inst=Entrada_con_unidades(Currency, retornar=False, readOnly=True)
        self.C_inst.entrada.setReadOnly(True)
        gridLayout_5.addWidget(self.C_inst,1,2)

        #Pestaña salida
        self.addSalida(QtGui.QApplication.translate("pychemqt", "Tube"))
        self.addSalida(QtGui.QApplication.translate("pychemqt", "Annulli"))

        if equipment:
            self.setEquipment(equipment)


    def showMaterial(self):
        dialogo=Catalogo_Materiales_Dialog()
        if dialogo.exec_():
            material=dialogo.getMaterial()
            if material:
                self.rTube.setValue(material[2])
                self.DiTube.setValue(material[4])
                self.wTube.setValue(material[5])
                self.DeTube.setValue(material[6])

    def showFinTube(self):
        dialogo=Dialog_Finned(self.Equipment.kwargs)
        if dialogo.exec_():
            kwarg=dialogo.kwarg()
            self.calculo(**kwarg)
    
    
if __name__ == "__main__":
    import sys        
    from lib.corriente import Corriente
    app = QtGui.QApplication(sys.argv)
    caliente=Corriente(T=140+273.15, P=361540., caudalMasico=1.36, ids=[62], fraccionMolar=[1.])
    fria=Corriente(T=20+273.15, P=101325., caudalMasico=5000/3600., ids=[62], fraccionMolar=[1.])
    Cambiador=Hairpin(entradaTubo=caliente, entradaExterior=fria, modo=1,  
                      DiTube=0.0525, DeTube=0.0603, DeeTube=0.0779, kTube=54, rTube=0.0459994e-3, 
                      annulliFouling= 0.000352, tubeFouling=0.000176, LTube=2.5)
    dialogo = UI_equipment(Cambiador)
    dialogo.show()
    sys.exit(app.exec_())
