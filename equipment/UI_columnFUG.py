#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################
###   Diálogo de definición de columnas de destilación UI_columnFUG    ###
###    Usa el método aproximado de Fenske-Underwood-Gilliland            ###
###############################################

from functools import partial

from PyQt4 import QtCore, QtGui

from lib.config import getComponents
from lib.unidades import Pressure, Volume, Length, Power, Density, Currency
from tools.costIndex import CostData
from equipment.parents import UI_equip
from equipment.distillation import ColumnFUG
from UI import UI_corriente
from UI.widgets import Entrada_con_unidades



class UI_equipment(UI_equip):
    """Diálogo de definición de columnas de destilación"""
    Equipment=ColumnFUG()
    def __init__(self, equipment=None, parent=None):
        """
        equipment: instancia de equipo inicial
        """
        super(UI_equipment, self).__init__(ColumnFUG, entrada=False, parent=parent)
        
        #Pestaña calculo
        gridLayout_Calculo = QtGui.QGridLayout(self.tabCalculo)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Feed tray")),2,0)
        self.feed=QtGui.QComboBox()
        for txt in self.Equipment.TEXT_FEED:
            self.feed.addItem(txt)
        self.feed.currentIndexChanged.connect(partial(self.changeParams, "feed"))
        gridLayout_Calculo.addWidget(self.feed,2,1)

        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Condenser")),3,0)
        self.condenser=QtGui.QComboBox()
        for txt in self.Equipment.TEXT_CONDENSER:
            self.condenser.addItem(txt)
        self.condenser.currentIndexChanged.connect(partial(self.changeParams, "condenser"))
        gridLayout_Calculo.addWidget(self.condenser,3,1)
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),4,0,1,5)

        groupbox=QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Key Components specification"))
        gridLayout_Calculo.addWidget(groupbox,5,0,1,5)
        layout=QtGui.QGridLayout(groupbox)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Light")),1,1)
        self.LK=QtGui.QComboBox()
        layout.addWidget(self.LK,1,2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Split in destilate")),1,4)
        self.LKsplit=Entrada_con_unidades(float, spinbox=True, max=1.)
        self.LKsplit.valueChanged.connect(partial(self.changeParams, "LKsplit"))
        layout.addWidget(self.LKsplit,1,5)
        
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Heavy")),2,1)
        self.HK=QtGui.QComboBox()
        layout.addWidget(self.HK,2,2)
        layout.addItem(QtGui.QSpacerItem(40,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),2,3)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Spit in residue")),2,4)
        self.HKsplit=Entrada_con_unidades(float, spinbox=True, max=1.)
        self.HKsplit.valueChanged.connect(partial(self.changeParams, "HKsplit"))
        layout.addWidget(self.HKsplit,2,5)
        
        indices, nombres, M=getComponents()
        for i, nombre in enumerate(nombres):
            self.LK.addItem("%i - %s" %(i+1, nombre))
            self.HK.addItem("%i - %s" %(i+1, nombre))
        self.LK.setCurrentIndex(-1)
        self.HK.setCurrentIndex(-1)
        self.LK.currentIndexChanged.connect(partial(self.changeParams, "LK"))
        self.HK.currentIndexChanged.connect(partial(self.changeParams, "HK"))
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),6,0,1,5)

        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Reflux ratio")),7,0)
        self.R=Entrada_con_unidades(float)
        self.R.valueChanged.connect(partial(self.changeParams, "R"))
        gridLayout_Calculo.addWidget(self.R,7,1)
        gridLayout_Calculo.addWidget(QtGui.QLabel("R/Rmin"),8,0)
        self.R_Rmin=Entrada_con_unidades(float)
        self.R_Rmin.valueChanged.connect(partial(self.changeParams, "R_Rmin"))
        gridLayout_Calculo.addWidget(self.R_Rmin,8,1)

        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Design Pressure")),7,3)
        self.Pd=Entrada_con_unidades(Pressure)
        self.Pd.valueChanged.connect(partial(self.changeParams, "Pd"))
        gridLayout_Calculo.addWidget(self.Pd,7,4)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Pressure loss")),8,3)
        self.DeltaP=Entrada_con_unidades(Pressure)
        self.DeltaP.valueChanged.connect(partial(self.changeParams, "DeltaP"))
        gridLayout_Calculo.addWidget(self.DeltaP,8,4)
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),9,0,1,5)
        self.buttonMcCabe=QtGui.QPushButton(QtGui.QApplication.translate("pychemqt", "McCabe-Thiele"))
        self.buttonMcCabe.clicked.connect(self.mcCabe)
        gridLayout_Calculo.addWidget(self.buttonMcCabe,10,0)
        
        groupBox_Calculo = QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Results"))
        gridLayout_Calculo.addWidget(groupBox_Calculo,11,0,1,5)
        gridLayout_1 = QtGui.QGridLayout(groupBox_Calculo)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Condenser Duty")),0,1)
        self.DutyCondenser=Entrada_con_unidades(Power, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.DutyCondenser,0,2,1,1)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Reboiler Duty")),1,1)
        self.DutyReboiler=Entrada_con_unidades(Power, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.DutyReboiler,1,2)
        gridLayout_1.addWidget(QtGui.QLabel("Rmin"),2,1)
        self.Rmin=Entrada_con_unidades(float, readOnly=True)
        gridLayout_1.addWidget(self.Rmin,2,2)    
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Reflux ratio")),3,1)
        self.RCalculada=Entrada_con_unidades(float, readOnly=True)
        gridLayout_1.addWidget(self.RCalculada,3,2)      

        gridLayout_1.addWidget(QtGui.QLabel("Nmin"),0,4)
        self.Nmin=Entrada_con_unidades(float, readOnly=True)
        gridLayout_1.addWidget(self.Nmin,0,5)      
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Stages")),1,4)
        self.NTray=Entrada_con_unidades(float, readOnly=True)
        gridLayout_1.addWidget(self.NTray,1,5)      
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Feed stage")),2,4)
        self.N_feed=Entrada_con_unidades(float, readOnly=True)
        gridLayout_1.addWidget(self.N_feed,2,5)      

        
        #Pestaña costos
        gridLayout_Costos = QtGui.QGridLayout(self.tabCostos)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Process")),1,1)
        self.proceso=QtGui.QComboBox()
        for txt in self.Equipment.TEXT_PROCESS:
            self.proceso.addItem(txt)
        self.proceso.currentIndexChanged.connect(partial(self.changeParamsCoste, "proceso"))
        gridLayout_Costos.addWidget(self.proceso,1,2)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Column tipe")),2,1)
        self.tipo=QtGui.QComboBox()
        for txt in self.Equipment.TEXT_COLUMN:
            self.tipo.addItem(txt)
        self.tipo.currentIndexChanged.connect(self.mostrarSubclasificacion)
        gridLayout_Costos.addWidget(self.tipo,2,2)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Material")),3,1)
        self.material=QtGui.QComboBox()
        for txt in self.Equipment.TEXT_MATERIAL:
            self.material.addItem(txt)
        self.material.currentIndexChanged.connect(partial(self.changeParamsCoste, "material_columna"))
        gridLayout_Costos.addWidget(self.material,3,2)

        gridLayout_Costos.addItem(QtGui.QSpacerItem(30,30,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),1,3,5,1)

        self.groupBox_Pisos = QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Tray column"))
        gridLayout_Costos.addWidget(self.groupBox_Pisos,1,4,4,2)
        gridLayout_1 = QtGui.QGridLayout(self.groupBox_Pisos)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Tray type")),1,1)
        self.tipoPisos=QtGui.QComboBox()
        for txt in self.Equipment.TEXT_TRAY:
            self.tipoPisos.addItem(txt)
        self.tipoPisos.currentIndexChanged.connect(partial(self.changeParamsCoste, "tipo_pisos"))
        gridLayout_1.addWidget(self.tipoPisos,1,2)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Material")),2,1)
        self.materialPisos=QtGui.QComboBox()
        for txt in self.Equipment.TEXT_MATERIAL:
            self.materialPisos.addItem(txt)
        self.materialPisos.currentIndexChanged.connect(partial(self.changeParamsCoste, "material_pisos"))
        gridLayout_1.addWidget(self.materialPisos,2,2)
        gridLayout_1.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),3,1,1,2)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Diameter")),4,1)
        self.diametroPisos=Entrada_con_unidades(Length)
        gridLayout_1.addWidget(self.diametroPisos,4,2)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Stages")),5,1)
        self.NumeroPisos=Entrada_con_unidades(int, spinbox=True, min=1, step=1, width=50)
        gridLayout_1.addWidget(self.NumeroPisos,5,2)
        gridLayout_1.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),6,1,1,2)

        self.groupBox_relleno = QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Packed column"))
        gridLayout_Costos.addWidget(self.groupBox_relleno,1,4,4,2)
        gridLayout_2 = QtGui.QGridLayout(self.groupBox_relleno)
        gridLayout_2.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Volume")),1,1)
        self.VolumenRelleno=Entrada_con_unidades(Volume, "VolLiq")
        gridLayout_2.addWidget(self.VolumenRelleno,1,2)
        gridLayout_2.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Unit Cost")),2,1)
        self.C_unit_relleno=Entrada_con_unidades(Currency, retornar=False, textounidad="%s / %s" % (Currency(None).text(), Volume(None).text("VolLiq")))
        gridLayout_2.addWidget(self.C_unit_relleno,2,2)
        gridLayout_2.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),3,1,1,2)
        
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Diameter")),5,1)
        self.Dc=Entrada_con_unidades(Length)
        gridLayout_Costos.addWidget(self.Dc,5,2,1,2)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Height")),6,1)
        self.Hc=Entrada_con_unidades(Length)
        gridLayout_Costos.addWidget(self.Hc,6,2,1,2)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Thickness (top)")),6,4)
        self.EspesorSuperior=Entrada_con_unidades(Length, "Thickness")
        gridLayout_Costos.addWidget(self.EspesorSuperior,6,5,1,2)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Thickness (bottom)")),7,4)
        self.EspesorInferior=Entrada_con_unidades(Length, "Thickness")
        gridLayout_Costos.addWidget(self.EspesorInferior,7,5,1,2)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Density")),7,1)
        self.EspesorInferior=Entrada_con_unidades(Density, "DenLiq")
        gridLayout_Costos.addWidget(self.EspesorInferior,7,2,1,2)
        
        self.Costos=CostData(self.Equipment)
        self.Costos.valueChanged.connect(self.changeParamsCoste)
        gridLayout_Costos.addWidget(self.Costos,10,1,2,5)

        gridLayout_Costos.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),12,1,1,6)
        
        self.groupBox_Costos = QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Stimated costs"))
        gridLayout_Costos.addWidget(self.groupBox_Costos,13,1,1,5)
        gridLayout_5 = QtGui.QGridLayout(self.groupBox_Costos)
        gridLayout_5.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Tray cost")),0,1)
        self.C_pisos=Entrada_con_unidades(Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_pisos,0,2)
        gridLayout_5.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Shell cost")),1,1)
        self.C_carcasa=Entrada_con_unidades(Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_carcasa,1,2)
        gridLayout_5.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Platform and ladder")),2,1)
        self.C_accesorios=Entrada_con_unidades(Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_accesorios,2,2)
        gridLayout_5.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Column cost")),0,4)
        self.C_columna=Entrada_con_unidades(Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_columna,0,5)
        gridLayout_5.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Purchase costs")),1,4)
        self.C_adq=Entrada_con_unidades(Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_adq,1,5)
        gridLayout_5.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Installed costs")),2,4)
        self.C_inst=Entrada_con_unidades(Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_inst,2,5)

        #Salidas
        self.SalidaDestilado= UI_corriente.Ui_corriente(readOnly=True)
        self.Salida.addTab(self.SalidaDestilado,QtGui.QApplication.translate("pychemqt", "Destilate"))
        self.SalidaResiduo= UI_corriente.Ui_corriente(readOnly=True)
        self.Salida.addTab(self.SalidaResiduo,QtGui.QApplication.translate("pychemqt", "Residue"))

        self.mostrarSubclasificacion(0)

        if equipment:
            self.setEquipment(equipment)

    def mcCabe(self):
        self.Equipment.McCabe()

    def mostrarSubclasificacion(self, ind):
        self.groupBox_Pisos.setVisible(not ind)
        self.groupBox_relleno.setVisible(ind)
        self.changeParamsCoste("tipo", ind)

    def rellenar(self):
        self.rellenarInput()
        self.buttonMcCabe.setEnabled(self.Equipment.statusMcCabe)
        if self.Equipment.status==1:
            UI_equip.rellenar(self)
            self.SalidaDestilado.setCorriente(self.Equipment.salida[0])
            self.SalidaResiduo.setCorriente(self.Equipment.salida[1])


if __name__ == "__main__":
    import sys        
    from lib.corriente import Corriente, Mezcla, Solid
    app = QtGui.QApplication(sys.argv)
    blend=Corriente(T=340, P=101325, caudalMasico=1, ids=[11, 12], fraccionMolar=[0.6, 0.4])
    columna=ColumnFUG(entrada=blend, LK=0, LKsplit=0.96666, HK=1, HKsplit=0.95, feed=0)
    dialogo = UI_equipment(columna)
    dialogo.show()
    sys.exit(app.exec_())
