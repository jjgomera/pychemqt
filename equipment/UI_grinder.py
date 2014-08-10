#!/usr/bin/python
# -*- coding: utf-8 -*-

#######################################################################
###                                              Diálogo de definición de molinos, UI_grinder                                                  ###
#######################################################################

from PyQt4 import QtCore, QtGui

from equipment.solids import Grinder
from UI import UI_corriente
from equipment import parents
from lib import unidades
from tools import costIndex
from UI.widgets import Entrada_con_unidades


BondIndex={   'Mineral de uranio': 17.93,
                        'Escoria': 15.76,
                        'Ferrocromo': 8.87,
                        'Grafito': 45.03,
                        'Magnesita': 16.8,
                        'Mineral de plata': 17.3,
                        'Molibdeno': 12.97,
                        'Ferromanganeso': 7.77,
                        'Arenisca': 11.53,
                        'Arcilla': 7.1,
                        'Mineral de níquel': 11.88,
                        'Mineral de estaño': 10.81,
                        'Mineral de titanio': 11.88,
                        'Silicato sódico': 13.0,
                        'Granito': 14.39,
                        'Coque': 20.7,
                        'Taconita': 14.87,
                        'Hematita especular': 15.4,
                        'Arena de silicato': 16.46,
                        'Coque de petróleo': 73.8,
                        'Gneiss': 20.13,
                        'Carburo de silicio': 26.17,
                        'Mineral de zinc': 12.42,
                        'Granate': 12.37,
                        'Caliza': 11.61,
                        'Basalto': 20.41,
                        'Carbón': 11.37,
                        'Gabro': 18.45,
                        'Dolomita': 11.31,
                        'Coque de petróleo líquido': 38.6,
                        'Mineral de plomo-zinc': 11.35,
                        'Sal potásica': 8.23,
                        'Andesita': 22.13,
                        'Arcilla calcinada': 1.43,
                        'Ilmenita': 13.11,
                        'Mineral de hierro': 15.44,
                        'Mica': 134.5,
                        'Hematita': 12.68,
                        'Fosfato fertilizante': 13.03,
                        'Cemento en bruto': 10.57,
                        'Bauxita': 9.45,
                        'Mineral de plomo': 11.4,
                        'Trapp': 21.1,
                        'Cristal': 3.08,
                        'Sienita': 14.9,
                        'Coral': 10.16,
                        'Roca fosfática': 10.13,
                        'Caliza para cemento': 10.18,
                        'Silicato': 13.53,
                        'Aljez': 8.16,
                        'Mineral de cromo': 9.6,
                        'Feldespato': 11.67,
                        'Mineral de cobre': 13.13,
                        'Pizarra bituminosa': 18.1,
                        'Cerámica': 15.53,
                        'Pirita': 8.9,
                        'Mineral de manganeso': 12.46,
                        'Pirrotina': 9.57,
                        'Cianita': 18.87,
                        'Grava': 25.17,
                        'Ferrosilicio': 12.83,
                        'Sílex': 26.16,
                        'Pizarra, mineral': 16.4,
                        'Limanita': 8.45,
                        'Barita': 6.24,
                        'Esmeril': 58.18,
                        'Escoria de hornos de hierro': 12.16,
                        'Mineral de oro': 14.83,
                        'Pumita': 11.93,
                        'Rutilo': 12.12,
                        'Espodumena': 13.7,
                        'Fluorita': 9.76,
                        'Clinker de cemento': 13.49,
                        'Sinterizado': 8.77,
                        'Galena': 10.19,
                        'Magnetita': 10.21,
                        'Cuarcita': 12.18,
                        'Oolita': 11.33,
                        'Pizarra, roca': 13.83,
                        'Mineral de potasio': 8.88,
                        'Diorita': 19.4,
                        'Cuarzo': 12.77}

class UI_equipment(parents.UI_equip):
    """Diálogo de definición de molinos trituradores de sólidos"""

    def __init__(self, entrada=None, parent=None):
        """entrada: Parametro opcional de clase corriente que indica la corriente de entrada en kla tubería"""
        super(UI_equipment, self).__init__(Grinder, entrada=False, salida=False, parent=parent)
        self.entrada=entrada

        #Pestaña entrada
        self.Entrada= UI_corriente.Ui_corriente(entrada)
        self.Entrada.Changed.connect(self.cambiar_entrada)
        self.tabWidget.insertTab(0, self.Entrada, QtGui.QApplication.translate("equipment", "Entrada", None, QtGui.QApplication.UnicodeUTF8))

        #Pestaña calculo
        gridLayout_Calculo = QtGui.QGridLayout(self.tabCalculo)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Índice de trabajo de bond:", None, QtGui.QApplication.UnicodeUTF8)), 1, 0, 1, 1)
        self.Material=QtGui.QComboBox()
        self.Material.addItem(QtGui.QApplication.translate("equipment", "Definido por el usuario", None, QtGui.QApplication.UnicodeUTF8))
        for key in sorted(BondIndex.keys()):
            self.Material.addItem(key)
        self.Material.currentIndexChanged[str].connect(self.cambiarBondWordIndex)
        gridLayout_Calculo.addWidget(self.Material, 1, 1, 1, 1)
        self.BondWorkIndex=Entrada_con_unidades(float)
        gridLayout_Calculo.addWidget(self.BondWorkIndex, 1, 2, 1, 1)
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),10,0,1,5)

        #Pestaña costos
        gridLayout_Costos = QtGui.QGridLayout(self.tabCostos)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Tipo:", None, QtGui.QApplication.UnicodeUTF8)), 1, 1, 1, 1)
        self.tipo=QtGui.QComboBox()
        self.tipo.addItem(QtGui.QApplication.translate("equipment", "De cono", None, QtGui.QApplication.UnicodeUTF8))
        self.tipo.addItem(QtGui.QApplication.translate("equipment", "Giratorio", None, QtGui.QApplication.UnicodeUTF8))
        self.tipo.addItem(QtGui.QApplication.translate("equipment", "Dentado", None, QtGui.QApplication.UnicodeUTF8))
        self.tipo.addItem(QtGui.QApplication.translate("equipment", "De martillo", None, QtGui.QApplication.UnicodeUTF8))
        self.tipo.addItem(QtGui.QApplication.translate("equipment", "De bolas", None, QtGui.QApplication.UnicodeUTF8))
        self.tipo.addItem(QtGui.QApplication.translate("equipment", "Pulverizador", None, QtGui.QApplication.UnicodeUTF8))
        self.tipo.currentIndexChanged.connect(self.calcularCostos)
        gridLayout_Costos.addWidget(self.tipo, 1, 2, 1, 1)
        gridLayout_Costos.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),2,1,1,2)

        self.Costos=costIndex.CostData(1.3, 2)
        self.Costos.valueChanged.connect(self.calcularCostos)
        gridLayout_Costos.addWidget(self.Costos,4,1,2,5)

        gridLayout_Costos.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),6,1,1,6)
        gridLayout_Costos.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),10,1,1,6)
        self.groupBox_Costos = QtGui.QGroupBox(QtGui.QApplication.translate("equipment", "Costos calculados", None, QtGui.QApplication.UnicodeUTF8))
        gridLayout_Costos.addWidget(self.groupBox_Costos,7,1,1,6)
        gridLayout_5 = QtGui.QGridLayout(self.groupBox_Costos)
        gridLayout_5.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Coste Adquisición:", None, QtGui.QApplication.UnicodeUTF8)),0,1,1,1)
        self.C_adq=Entrada_con_unidades(unidades.Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_adq,0,2,1,1)
        gridLayout_5.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Coste Instalación:", None, QtGui.QApplication.UnicodeUTF8)),1,1,1,1)
        self.C_inst=Entrada_con_unidades(unidades.Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_inst,1,2,1,1)

        #Pestaña salida
        self.Salida= UI_corriente.Ui_corriente(readOnly=True)
        self.tabWidget.insertTab(3, self.Salida,QtGui.QApplication.translate("equipment", "Salida", None, QtGui.QApplication.UnicodeUTF8))

        self.tabWidget.setCurrentIndex(0)


    def cambiarBondWordIndex(self, txt):
        try:
            value=BondIndex[str(txt)]
        except KeyError:
            self.BondWorkIndex.setReadOnly(False)
            self.BondWorkIndex.clear()
        else:
            self.BondWorkIndex.setValue(value)
            self.BondWorkIndex.setReadOnly(True)

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

    def calcularCostos(self, factor=None, indiceBase=None, indiceActual=None):
        if self.todos_datos():
            if not factor: factor=self.Costos.factor
            if not indiceBase: indiceBase=self.Costos.Base
            if not indiceActual: indiceActual=self.Costos.Actual
            if self.tipo.currentIndex()==0:
                self.FireHeater.Coste(factor, indiceBase, indiceActual, 0, self.tipobox.currentIndex(), self.material.currentIndex())
            else:
                self.FireHeater.Coste(factor, indiceBase, indiceActual, 1, self.tipocilindrico.currentIndex(), self.material.currentIndex())
            self.C_adq.setValue(self.FireHeater.C_adq.config())
            self.C_inst.setValue(self.FireHeater.C_inst.config())

if __name__ == "__main__":
    import sys
    from lib.corriente import Corriente, Mezcla, Solid
    app = QtGui.QApplication(sys.argv)
    agua=Corriente(300, 1, 3600, Mezcla([62], [1]))
    dialogo = UI_equipment(agua)
    dialogo.show()
    sys.exit(app.exec_())
