#!/usr/bin/python
# -*- coding: utf-8 -*-

#######################################################################
###                             Diálogo de definición de tanques de almacenamiento, UI_tank                                        ###
#######################################################################

from PyQt4 import QtCore, QtGui

from equipment.tank import Tank
from UI import UI_corriente
from equipment import parents
from lib import unidades
from tools import costIndex
from UI.widgets import Entrada_con_unidades


class UI_equipment(parents.UI_equipment):
    """Diálogo de definición de tuberías"""
    def __init__(self, entrada=None, parent=None):
        """entrada: Parametro opcional de clase corriente que indica la corriente de entrada en kla tubería"""
        super(UI_equipment, self).__init__(Tank, entrada=False, salida=False, parent=parent)
        self.entrada=entrada
        
        #Pestaña entrada
        self.Entrada= UI_corriente.Ui_corriente(entrada)
        self.Entrada.Changed.connect(self.cambiar_entrada)
        self.tabWidget.insertTab(0, self.Entrada,QtGui.QApplication.translate("equipment", "Entrada", None, QtGui.QApplication.UnicodeUTF8))
        
        #Pestaña calculo
        gridLayout_Calculo = QtGui.QGridLayout(self.tabCalculo)


        #Pestaña costos
        gridLayout_Costos = QtGui.QGridLayout(self.tabCostos)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Material:", None, QtGui.QApplication.UnicodeUTF8)), 1, 1, 1, 1)
        self.material=QtGui.QComboBox()
        self.material.addItem(QtGui.QApplication.translate("equipment", "Acero al carbon", None, QtGui.QApplication.UnicodeUTF8))
        self.material.addItem(QtGui.QApplication.translate("equipment", "Acero inoxidable 316", None, QtGui.QApplication.UnicodeUTF8))
        self.material.addItem(QtGui.QApplication.translate("equipment", "Acero inoxidable 304", None, QtGui.QApplication.UnicodeUTF8))
        self.material.addItem(QtGui.QApplication.translate("equipment", "Acero inoxidable 347", None, QtGui.QApplication.UnicodeUTF8))
        self.material.addItem(QtGui.QApplication.translate("equipment", "Niquel", None, QtGui.QApplication.UnicodeUTF8))
        self.material.addItem(QtGui.QApplication.translate("equipment", "Monel", None, QtGui.QApplication.UnicodeUTF8))
        self.material.addItem(QtGui.QApplication.translate("equipment", "Inconel", None, QtGui.QApplication.UnicodeUTF8))
        self.material.addItem(QtGui.QApplication.translate("equipment", "Zirconio", None, QtGui.QApplication.UnicodeUTF8))
        self.material.addItem(QtGui.QApplication.translate("equipment", "Titanio", None, QtGui.QApplication.UnicodeUTF8))
        self.material.addItem(QtGui.QApplication.translate("equipment", "Ladrillo y caucho o ladrillo y acero recubierto de poliester", None, QtGui.QApplication.UnicodeUTF8))
        self.material.addItem(QtGui.QApplication.translate("equipment", "Caucho o acero recubierto de plomo", None, QtGui.QApplication.UnicodeUTF8))
        self.material.addItem(QtGui.QApplication.translate("equipment", "Poliester reforzado con fiberglass", None, QtGui.QApplication.UnicodeUTF8))
        self.material.addItem(QtGui.QApplication.translate("equipment", "Aluminio", None, QtGui.QApplication.UnicodeUTF8))
        self.material.addItem(QtGui.QApplication.translate("equipment", "Cobre", None, QtGui.QApplication.UnicodeUTF8))
        self.material.addItem(QtGui.QApplication.translate("equipment", "Hormigón", None, QtGui.QApplication.UnicodeUTF8))
        self.material.currentIndexChanged.connect(self.calcularCostos)
        gridLayout_Costos.addWidget(self.material, 1, 2, 1, 4)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Densidad:", None, QtGui.QApplication.UnicodeUTF8)), 2, 4, 1, 1)
        self.Densidad=Entrada_con_unidades(unidades.Density, "DenLiq")
        gridLayout_Costos.addWidget(self.Densidad,2,5,1,1)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Diametro:", None, QtGui.QApplication.UnicodeUTF8)), 2, 1, 1, 1)
        self.Diametro=Entrada_con_unidades(unidades.Length)
        gridLayout_Costos.addWidget(self.Diametro,2,2,1,1)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Longitud:", None, QtGui.QApplication.UnicodeUTF8)), 3, 1, 1, 1)
        self.Longitud=Entrada_con_unidades(unidades.Length)
        gridLayout_Costos.addWidget(self.Longitud,3,2,1,1)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Espesor:", None, QtGui.QApplication.UnicodeUTF8)), 4, 1, 1, 1)
        self.Espesor=Entrada_con_unidades(unidades.Length, "Thickness")
        gridLayout_Costos.addWidget(self.Espesor,4,2,1,1)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Cabeza:", None, QtGui.QApplication.UnicodeUTF8)), 5, 1, 1, 1)
        self.Cabeza=QtGui.QComboBox()
        self.Cabeza.addItem(QtGui.QApplication.translate("equipment", "Elipsoidal", None, QtGui.QApplication.UnicodeUTF8))
        self.Cabeza.addItem(QtGui.QApplication.translate("equipment", "Semiesférico", None, QtGui.QApplication.UnicodeUTF8))
        self.Cabeza.addItem(QtGui.QApplication.translate("equipment", "Bumped", None, QtGui.QApplication.UnicodeUTF8))
        self.Cabeza.addItem(QtGui.QApplication.translate("equipment", "Liso", None, QtGui.QApplication.UnicodeUTF8))        
        self.Cabeza.currentIndexChanged.connect(self.calcularCostos)
        gridLayout_Costos.addWidget(self.Cabeza, 5, 2, 1, 1) 
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Espesor (cabeza):", None, QtGui.QApplication.UnicodeUTF8)), 6, 1, 1, 1)
        self.EspesorCabeza=Entrada_con_unidades(unidades.Length, "Thickness")
        gridLayout_Costos.addWidget(self.EspesorCabeza,6,2,1,1)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Longitud reborde recto:", None, QtGui.QApplication.UnicodeUTF8)), 7, 1, 1, 1)
        self.LongitudReborde=Entrada_con_unidades(unidades.Length)
        gridLayout_Costos.addWidget(self.LongitudReborde,7,2,1,1)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Volumen:", None, QtGui.QApplication.UnicodeUTF8)), 6, 4, 1, 1)
        self.Volumen=Entrada_con_unidades(unidades.Volume, "VolLiq", readOnly=True)
        gridLayout_Costos.addWidget(self.Volumen,6,5,1,1)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Peso:", None, QtGui.QApplication.UnicodeUTF8)), 7, 4, 1, 1)
        self.Peso=Entrada_con_unidades(unidades.Mass, readOnly=True)
        gridLayout_Costos.addWidget(self.Peso,7,5,1,1)
        gridLayout_Costos.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),2,3,6,1)
        gridLayout_Costos.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),8,0,1,6)
        
        self.Costos=costIndex.CostData(1.7, 3)
        self.Costos.valueChanged.connect(self.calcularCostos)
        gridLayout_Costos.addWidget(self.Costos,9,1,2,5)

        gridLayout_Costos.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),11,0,1,6)
        self.groupBox_Costos = QtGui.QGroupBox(QtGui.QApplication.translate("equipment", "Costos calculados", None, QtGui.QApplication.UnicodeUTF8))
        gridLayout_Costos.addWidget(self.groupBox_Costos,12,1,1,5)
        gridLayout_5 = QtGui.QGridLayout(self.groupBox_Costos)
        gridLayout_5.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Coste Adquisición:", None, QtGui.QApplication.UnicodeUTF8)),0,1,1,1)
        self.C_adq=Entrada_con_unidades(unidades.Currency, retornar=False, readOnly=True, tolerancia=8, decimales=2)
        gridLayout_5.addWidget(self.C_adq,0,2,1,1)
        gridLayout_5.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Coste Instalación:", None, QtGui.QApplication.UnicodeUTF8)),1,1,1,1)
        self.C_inst=Entrada_con_unidades(unidades.Currency, retornar=False, readOnly=True, tolerancia=8, decimales=2)
        gridLayout_5.addWidget(self.C_inst,1,2,1,1)
        gridLayout_Costos.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),13,0,1,6)
        
        #Pestaña salida
        self.Salida= UI_corriente.Ui_corriente(readOnly=True)
        self.tabWidget.insertTab(2, self.Salida,QtGui.QApplication.translate("equipment", "Salida", None, QtGui.QApplication.UnicodeUTF8))

        self.tabWidget.setCurrentIndex(0)

        
    def cambiar_entrada(self, corriente):
        selfentrada=corriente
        self.calculo()

    def calculo(self):
        if self.todos_datos():
            
            self.rellenoSalida()
    
    def rellenoSalida(self):
        pass    
        
    def todos_datos(self):
        pass

    def calcularCostos(self):
        if self.todos_datos():
            if self.tipo.currentIndex()==0:
                self.FireHeater.Coste(self.factorInstalacion.value(), 0, self.tipobox.currentIndex(), self.material.currentIndex())
            else:
                self.FireHeater.Coste(self.factorInstalacion.value(), 1, self.tipocilindrico.currentIndex(), self.material.currentIndex())
            self.C_adq.setValue(self.FireHeater.C_adq.config())
            self.C_inst.setValue(self.FireHeater.C_inst.config())



if __name__ == "__main__":
    import sys        
    from lib.corriente import Corriente, Mezcla
    app = QtGui.QApplication(sys.argv)
    agua=Corriente(300, 1, 3600, Mezcla([62], [1]))
    dialogo = UI_equipment(agua)
    dialogo.show()
    sys.exit(app.exec_())
