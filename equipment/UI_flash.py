#!/usr/bin/python
# -*- coding: utf-8 -*-

#######################################
###     Diálogo de definición de tambores flash, UI_flash     ###
#######################################

from functools import partial

from PyQt4 import QtCore, QtGui

from lib.unidades import Length, Mass, Volume, Density, Currency
from tools.costIndex import CostData
from equipment.parents import UI_equip
from equipment.distillation import Flash
from UI import UI_corriente
from UI.widgets import Entrada_con_unidades


class UI_equipment (UI_equip):
    """Diálogo de definición de tambores flash"""
    Equipment=Flash()
    def __init__(self, equipment=None, parent=None):
        """
        equipment: instancia de equipo inicial
        """
        super(UI_equipment , self).__init__(Flash, entrada=False, parent=parent)
        
        #Pestaña calculo
        gridLayout_Calculo = QtGui.QGridLayout(self.tabCalculo)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Method")),0,1)
        self.flash=QtGui.QComboBox()
        for txt in self.Equipment.TEXT_FLASH:
            self.flash.addItem(txt)
        self.flash.currentIndexChanged.connect(partial(self.changeParams, "metodo"))
        gridLayout_Calculo.addWidget(self.flash,0,2)
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),1,1,1,6)
        
        #Pestaña costos
        gridLayout_Costos = QtGui.QGridLayout(self.tabCostos)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Orientation")),0,1)
        self.orientacion=QtGui.QComboBox()
        for txt in self.Equipment.TEXT_ORIENTATION:
            self.orientacion.addItem(txt)
        self.orientacion.currentIndexChanged.connect(partial(self.changeParamsCoste, "orientacion"))
        gridLayout_Costos.addWidget(self.orientacion,0,2)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Material")),1,1)
        self.material=QtGui.QComboBox()
        for txt in self.Equipment.TEXT_MATERIAL:
            self.material.addItem(txt)
        self.material.currentIndexChanged.connect(partial(self.changeParamsCoste, "material"))
        gridLayout_Costos.addWidget(self.material,1,2,1,4)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Density")),2,4)
        self.Densidad=Entrada_con_unidades(Density, "DenLiq")
        self.Densidad.valueChanged.connect(partial(self.changeParamsCoste, "densidad"))
        gridLayout_Costos.addWidget(self.Densidad,2,5)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Diameter")),2,1)
        self.diametro=Entrada_con_unidades(Length)
        self.diametro.valueChanged.connect(partial(self.changeParamsCoste, "diametro"))
        gridLayout_Costos.addWidget(self.diametro,2,2)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Length")),3,1)
        self.longitud=Entrada_con_unidades(Length)
        self.longitud.valueChanged.connect(partial(self.changeParamsCoste, "longitud"))
        gridLayout_Costos.addWidget(self.longitud,3,2)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Thickness")),4,1)
        self.espesor=Entrada_con_unidades(Length, "Thickness")
        self.espesor.valueChanged.connect(partial(self.changeParamsCoste, "espesor"))
        gridLayout_Costos.addWidget(self.espesor,4,2)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Head type")),5,1)
        self.cabeza=QtGui.QComboBox()
        for txt in self.Equipment.TEXT_HEAD:
            self.cabeza.addItem(txt)
        self.cabeza.currentIndexChanged.connect(partial(self.changeParamsCoste, "cabeza"))
        gridLayout_Costos.addWidget(self.cabeza,5,2) 
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Head Thickness")),6,1)
        self.espesor_cabeza=Entrada_con_unidades(Length, "Thickness")
        self.espesor_cabeza.valueChanged.connect(partial(self.changeParamsCoste, "espesor_cabeza"))
        gridLayout_Costos.addWidget(self.espesor_cabeza,6,2)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Straight flange length")),7,1)
        self.reborde=Entrada_con_unidades(Length)
        self.reborde.valueChanged.connect(partial(self.changeParamsCoste, "reborde"))
        gridLayout_Costos.addWidget(self.reborde,7,2)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Volume")),6,4)
        self.Volumen=Entrada_con_unidades(Volume, "VolLiq", readOnly=True, retornar=False)
        gridLayout_Costos.addWidget(self.Volumen,6,5)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Weight")),7,4)
        self.Peso=Entrada_con_unidades(Mass, readOnly=True)
        gridLayout_Costos.addWidget(self.Peso,7,5,1,1)
        gridLayout_Costos.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),2,3,6,1)
        gridLayout_Costos.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),8,0,1,6)
        
        self.Costos=CostData(self.Equipment)
        self.Costos.valueChanged.connect(self.changeParamsCoste)
        gridLayout_Costos.addWidget(self.Costos,9,1,2,5)
        
        gridLayout_Costos.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),11,0,1,6)
        self.groupBox_Costos = QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Stimated Costs"))
        gridLayout_Costos.addWidget(self.groupBox_Costos,12,1,1,5)
        gridLayout_5 = QtGui.QGridLayout(self.groupBox_Costos)
        gridLayout_5.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Purchase costs")),0,1)
        self.C_adq=Entrada_con_unidades(Currency, retornar=False, readOnly=True, tolerancia=8, decimales=2)
        gridLayout_5.addWidget(self.C_adq,0,2)
        gridLayout_5.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Installed costs")),1,1)
        self.C_inst=Entrada_con_unidades(Currency, retornar=False, readOnly=True, tolerancia=8, decimales=2)
        gridLayout_5.addWidget(self.C_inst,1,2)
        gridLayout_Costos.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),13,0,1,6)
        
        #Pestaña salida
        self.SalidaVapor= UI_corriente.Ui_corriente(readOnly=True)
        self.Salida.addTab(self.SalidaVapor,QtGui.QApplication.translate("pychemqt", "Gas"))
        self.SalidaLiquido= UI_corriente.Ui_corriente(readOnly=True)
        self.Salida.addTab(self.SalidaLiquido,QtGui.QApplication.translate("pychemqt", "Liquid"))

        if equipment:
            self.setEquipment(equipment)


    def rellenar(self):
        self.rellenarInput()
        if self.Equipment.status==1:
            UI_equip.rellenar(self)
            
            self.SalidaVapor.setCorriente(self.Equipment.salida[0])
            self.SalidaLiquido.setCorriente(self.Equipment.salida[1])


if __name__ == "__main__":
    import sys        
    from lib.corriente import Corriente
    app = QtGui.QApplication(sys.argv)
    entrada=Corriente(T=340, P=101325, caudalMasico=0.01, ids=[10, 38, 22, 61], fraccionMolar=[.3, 0.5, 0.05, 0.15])
    flash=Flash(entrada=entrada)
    dialogo = UI_equipment(flash)
    dialogo.show()
    sys.exit(app.exec_())
