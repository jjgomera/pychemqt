#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2023, Juan José Gómez Romera <jjgomera@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''



#######################################################################
###                                              Diálogo de definición de molinos, UI_grinder                                                  ###
#######################################################################

from tools.qt import QtCore, QtWidgets

from equipment.parents import UI_equip
from equipment.solids import Grinder
from lib import unidades
from tools import costIndex
from UI import UI_corriente
from UI.widgets import Entrada_con_unidades


BondIndex = {
    'Mineral de uranio': 17.93,
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


class UI_equipment(UI_equip):
    """Diálogo de definición de molinos trituradores de sólidos"""
    Equipment = Grinder()

    def __init__(self, entrada=None, parent=None):
        """entrada: Parametro opcional de clase corriente que indica la corriente de entrada en kla tubería"""
        super(UI_equipment, self).__init__(Grinder, entrada=False, salida=False, parent=parent)

        #Pestaña calculo
        gridLayout_Calculo = QtWidgets.QGridLayout(self.tabCalculo)
        gridLayout_Calculo.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("equipment", "Índice de trabajo de bond:", None)), 1, 0, 1, 1)
        self.Material=QtWidgets.QComboBox()
        self.Material.addItem(QtWidgets.QApplication.translate("equipment", "Definido por el usuario", None))
        for key in sorted(BondIndex.keys()):
            self.Material.addItem(key)
        self.Material.currentIndexChanged[str].connect(self.cambiarBondWordIndex)
        gridLayout_Calculo.addWidget(self.Material, 1, 1, 1, 1)
        self.BondWorkIndex=Entrada_con_unidades(float)
        gridLayout_Calculo.addWidget(self.BondWorkIndex, 1, 2, 1, 1)
        gridLayout_Calculo.addItem(QtWidgets.QSpacerItem(10,10,QtWidgets.QSizePolicy.Policy.Expanding,QtWidgets.QSizePolicy.Policy.Expanding),10,0,1,5)

        #Pestaña costos
        gridLayout_Costos = QtWidgets.QGridLayout(self.tabCostos)
        gridLayout_Costos.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("equipment", "Tipo:", None)), 1, 1, 1, 1)
        self.tipo=QtWidgets.QComboBox()
        self.tipo.addItem(QtWidgets.QApplication.translate("equipment", "De cono", None))
        self.tipo.addItem(QtWidgets.QApplication.translate("equipment", "Giratorio", None))
        self.tipo.addItem(QtWidgets.QApplication.translate("equipment", "Dentado", None))
        self.tipo.addItem(QtWidgets.QApplication.translate("equipment", "De martillo", None))
        self.tipo.addItem(QtWidgets.QApplication.translate("equipment", "De bolas", None))
        self.tipo.addItem(QtWidgets.QApplication.translate("equipment", "Pulverizador", None))
        self.tipo.currentIndexChanged.connect(self.calcularCostos)
        gridLayout_Costos.addWidget(self.tipo, 1, 2, 1, 1)
        gridLayout_Costos.addItem(QtWidgets.QSpacerItem(10,10,QtWidgets.QSizePolicy.Policy.Fixed,QtWidgets.QSizePolicy.Policy.Fixed),2,1,1,2)

        self.Costos = costIndex.CostData(self.Equipment)
        self.Costos.valueChanged.connect(self.calcularCostos)
        gridLayout_Costos.addWidget(self.Costos,4,1,2,5)

        gridLayout_Costos.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Policy.Expanding,QtWidgets.QSizePolicy.Policy.Expanding),6,1,1,6)
        gridLayout_Costos.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Policy.Expanding,QtWidgets.QSizePolicy.Policy.Expanding),10,1,1,6)
        self.groupBox_Costos = QtWidgets.QGroupBox(QtWidgets.QApplication.translate("equipment", "Costos calculados", None))
        gridLayout_Costos.addWidget(self.groupBox_Costos,7,1,1,6)
        gridLayout_5 = QtWidgets.QGridLayout(self.groupBox_Costos)
        gridLayout_5.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("equipment", "Coste Adquisición:", None)),0,1,1,1)
        self.C_adq=Entrada_con_unidades(unidades.Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_adq,0,2,1,1)
        gridLayout_5.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("equipment", "Coste Instalación:", None)),1,1,1,1)
        self.C_inst=Entrada_con_unidades(unidades.Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_inst,1,2,1,1)


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
    app = QtWidgets.QApplication(sys.argv)
    agua=Corriente(T=300, P=101325, caudalMasico=1., ids=[62], fraccionMolar=[1.], MEoS=True)
    dialogo = UI_equipment(agua)
    dialogo.show()
    sys.exit(app.exec())
