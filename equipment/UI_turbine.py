#!/usr/bin/python
# -*- coding: utf-8 -*-

#######################################
###      Diálogo de definición de expansores, UI_turbine      ###
#######################################

from functools import partial

from PyQt4 import QtCore, QtGui

from lib.unidades import Pressure, Power, Currency
from tools.costIndex import CostData
from equipment.parents import UI_equip
from equipment.compressor import Turbine
from UI import UI_corriente
from UI.widgets import Entrada_con_unidades


class UI_equipment(UI_equip):
    """Diálogo de definición de turbinas o expansores de gases"""
    Equipment=Turbine()
    def __init__(self, equipment=None, parent=None):
        """entrada: Parametro opcional de clase corriente que indica la corriente de entrada"""
        super(UI_equipment, self).__init__(Turbine, entrada=False, salida=False, parent=parent)
                
        #Pestaña calculo
        gridLayout_Calculo = QtGui.QGridLayout(self.tabCalculo)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Method:")),1,1)
        self.metodo=QtGui.QComboBox()
        for txt in self.Equipment.TEXT_METODO:
            self.metodo.addItem(txt)
#        self.metodo.addItem(QtGui.QApplication.translate("pychemqt", "Especificar curva de funcionamiento"))
        self.metodo.currentIndexChanged.connect(self.on_tipoCalculo_currentIndexChanged)
        gridLayout_Calculo.addWidget(self.metodo, 1, 2, 1, 2)
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(40,20,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),2,0,1,4)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Thermodynamic:")),3,1)
        self.termodinamica=QtGui.QComboBox()
        for txt in self.Equipment.TEXT_TERMODINAMICA:
            self.termodinamica.addItem(txt)
        self.termodinamica.currentIndexChanged.connect(partial(self.changeParams, "termodinamica"))
        gridLayout_Calculo.addWidget(self.termodinamica, 3, 2, 1, 2)
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(40,20,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),4,0,1,4)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Out Pressure")),5,1)
        self.Pout=Entrada_con_unidades(Pressure)
        self.Pout.valueChanged.connect(partial(self.changeParams, "Pout"))
        gridLayout_Calculo.addWidget(self.Pout,5,2)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Pressure ratio")),6,1)
        self.razon=Entrada_con_unidades(float)
        self.razon.valueChanged.connect(partial(self.changeParams, "razon"))
        gridLayout_Calculo.addWidget(self.razon,6,2)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Efficiency")),7,1)
        self.rendimiento=Entrada_con_unidades(float)
        self.rendimiento.valueChanged.connect(partial(self.changeParams, "rendimiento"))
        gridLayout_Calculo.addWidget(self.rendimiento,7,2)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Actual Power")),8,1)
        self.trabajo=Entrada_con_unidades(Power)
        self.trabajo.valueChanged.connect(partial(self.changeParams, "trabajo"))
        gridLayout_Calculo.addWidget(self.trabajo,8,2)
        gridLayout_Calculo.setRowStretch(10, 1)
        
        self.groupBox_Resultados = QtGui.QGroupBox()
        self.groupBox_Resultados.setTitle(QtGui.QApplication.translate("pychemqt", "Results"))
        gridLayout_Calculo.addWidget(self.groupBox_Resultados,12,1,1,2)
        gridLayout_3 = QtGui.QGridLayout(self.groupBox_Resultados)
        gridLayout_3.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Power")),1,1)
        self.power=Entrada_con_unidades(Power, retornar=False, readOnly=True)
        gridLayout_3.addWidget(self.power,1,2)
        gridLayout_3.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Cp/Cv")),2,1)
        self.cp_cv=Entrada_con_unidades(float, retornar=False, readOnly=True)
        gridLayout_3.addWidget(self.cp_cv,2,2)
        gridLayout_3.setColumnStretch(3, 1)
        gridLayout_3.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Pressure ratio")),1,4)
        self.razonCalculada=Entrada_con_unidades(float, readOnly=True)
        gridLayout_3.addWidget(self.razonCalculada,1,5)
        gridLayout_3.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Efficiency")),2,4)
        self.rendimientoCalculado=Entrada_con_unidades(float, readOnly=True)
        gridLayout_3.addWidget(self.rendimientoCalculado,2,5)


        #Pestaña costos
        layout_costos = QtGui.QGridLayout(self.tabCostos)
        self.Costos=CostData(self.Equipment)
        self.Costos.valueChanged.connect(self.changeParamsCoste)
        layout_costos.addWidget(self.Costos,1,0,1,2)
        
        layout_costos.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),2,0,1,2)
        self.groupBox_Costos = QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Stimated Costs"))
        layout_costos.addWidget(self.groupBox_Costos,3,0,1,2)
        gridLayout_5 = QtGui.QGridLayout(self.groupBox_Costos)
        gridLayout_5.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Purchase cost")),0,1)
        self.C_adq=Entrada_con_unidades(Currency, retornar=False, tolerancia=8, decimales=2, readOnly=True)
        gridLayout_5.addWidget(self.C_adq,0,2)
        gridLayout_5.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Installed cost")),1,1)
        self.C_inst=Entrada_con_unidades(Currency, retornar=False, tolerancia=8, decimales=2, readOnly=True)
        gridLayout_5.addWidget(self.C_inst,1,2)
        layout_costos.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),4,0,1,2)
        
        self.on_tipoCalculo_currentIndexChanged(0)
        
        if equipment:
            self.setEquipment(equipment)


    def on_tipoCalculo_currentIndexChanged(self,int):
        """Habilita o desabilita los datos requeridos para el cálculo"""
        if int == 0:
            self.trabajo.setReadOnly(True)
            self.razon.setReadOnly(True)
            self.Pout.setReadOnly(False)
            self.rendimiento.setReadOnly(False)
            self.Pout.setResaltado(True)
            self.rendimiento.setResaltado(True)
            self.razon.setResaltado(False)
            self.trabajo.setResaltado(False)
        elif int==1:
            self.Pout.setReadOnly(True)
            self.razon.setReadOnly(False)
            self.trabajo.setReadOnly(True)
            self.rendimiento.setReadOnly(False)
            self.Pout.setResaltado(False)
            self.rendimiento.setResaltado(True)
            self.razon.setResaltado(True)
            self.trabajo.setResaltado(False)
        elif int==2:
            self.trabajo.setReadOnly(False)
            self.Pout.setReadOnly(True)
            self.razon.setReadOnly(True)
            self.rendimiento.setReadOnly(False)
            self.Pout.setResaltado(False)
            self.rendimiento.setResaltado(True)
            self.razon.setResaltado(False)
            self.trabajo.setResaltado(True)
        elif int==3:
            self.rendimiento.setReadOnly(True)
            self.razon.setReadOnly(True)
            self.trabajo.setReadOnly(False)
            self.Pout.setReadOnly(False)
            self.Pout.setResaltado(True)
            self.rendimiento.setResaltado(False)
            self.razon.setResaltado(False)
            self.trabajo.setResaltado(True)
        elif int==4:
            self.rendimiento.setReadOnly(True)
            self.razon.setReadOnly(False)
            self.trabajo.setReadOnly(False)
            self.Pout.setReadOnly(True)
            self.Pout.setResaltado(False)
            self.rendimiento.setResaltado(False)
            self.razon.setResaltado(True)
            self.trabajo.setResaltado(True)
        else:
            self.rendimiento.setReadOnly(False)
            self.razon.setReadOnly(False)
            self.trabajo.setReadOnly(False)
            self.Pout.setReadOnly(False)
            self.Pout.setResaltado(True)
            self.rendimiento.setResaltado(True)
            self.razon.setResaltado(True)
            self.trabajo.setResaltado(True)
            
        self.changeParams("metodo", int)

    
    def rellenar(self):
        UI_equip.rellenar(self)
        if self.Equipment.status==1 and self.metodo.currentIndex()==5:
                self.entrada.setCorriente(self.Equipment.entrada)


if __name__ == "__main__":
    import sys        
    from lib.corriente import Corriente, Mezcla
    app = QtGui.QApplication(sys.argv)
    agua=Corriente(T=500, P=3*101325, caudalMasico=1, fraccionMasica=[1.])
    turbina=Turbine(entrada=agua, metodo=1, razon=0.3, rendimiento=1.)
    dialogo = UI_equipment(turbina)
    dialogo.show()
    sys.exit(app.exec_())
