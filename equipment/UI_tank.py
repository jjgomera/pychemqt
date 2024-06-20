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
###                             Diálogo de definición de tanques de almacenamiento, UI_tank                                        ###
#######################################################################

from tools.qt import QtCore, QtWidgets


from equipment.tank import Tank
from UI import UI_corriente
from equipment import parents
from lib import unidades
from tools import costIndex
from UI.widgets import Entrada_con_unidades


class UI_equipment(parents.UI_equip):
    """Diálogo de definición de tuberías"""
    def __init__(self, entrada=None, parent=None):
        """entrada: Parametro opcional de clase corriente que indica la corriente de entrada en kla tubería"""
        super().__init__(Tank, entrada=False, salida=False, parent=parent)
        self.entrada=entrada

        #Pestaña entrada
        self.Entrada= UI_corriente.Ui_corriente(entrada)
        self.Entrada.Changed.connect(self.cambiar_entrada)
        self.tabWidget.insertTab(0, self.Entrada,self.tr("Entrada", None))

        #Pestaña calculo
        gridLayout_Calculo = QtWidgets.QGridLayout(self.tabCalculo)


        #Pestaña costos
        gridLayout_Costos = QtWidgets.QGridLayout(self.tabCostos)
        gridLayout_Costos.addWidget(QtWidgets.QLabel(self.tr("Material:", None)), 1, 1, 1, 1)
        self.material=QtWidgets.QComboBox()
        self.material.addItem(self.tr("Acero al carbon", None))
        self.material.addItem(self.tr("Acero inoxidable 316", None))
        self.material.addItem(self.tr("Acero inoxidable 304", None))
        self.material.addItem(self.tr("Acero inoxidable 347", None))
        self.material.addItem(self.tr("Niquel", None))
        self.material.addItem(self.tr("Monel", None))
        self.material.addItem(self.tr("Inconel", None))
        self.material.addItem(self.tr("Zirconio", None))
        self.material.addItem(self.tr("Titanio", None))
        self.material.addItem(self.tr("Ladrillo y caucho o ladrillo y acero recubierto de poliester", None))
        self.material.addItem(self.tr("Caucho o acero recubierto de plomo", None))
        self.material.addItem(self.tr("Poliester reforzado con fiberglass", None))
        self.material.addItem(self.tr("Aluminio", None))
        self.material.addItem(self.tr("Cobre", None))
        self.material.addItem(self.tr("Hormigón", None))
        self.material.currentIndexChanged.connect(self.calcularCostos)
        gridLayout_Costos.addWidget(self.material, 1, 2, 1, 4)
        gridLayout_Costos.addWidget(QtWidgets.QLabel(self.tr("Densidad:", None)), 2, 4, 1, 1)
        self.Densidad=Entrada_con_unidades(unidades.Density, "DenLiq")
        gridLayout_Costos.addWidget(self.Densidad,2,5,1,1)
        gridLayout_Costos.addWidget(QtWidgets.QLabel(self.tr("Diametro:", None)), 2, 1, 1, 1)
        self.Diametro=Entrada_con_unidades(unidades.Length)
        gridLayout_Costos.addWidget(self.Diametro,2,2,1,1)
        gridLayout_Costos.addWidget(QtWidgets.QLabel(self.tr("Longitud:", None)), 3, 1, 1, 1)
        self.Longitud=Entrada_con_unidades(unidades.Length)
        gridLayout_Costos.addWidget(self.Longitud,3,2,1,1)
        gridLayout_Costos.addWidget(QtWidgets.QLabel(self.tr("Espesor:", None)), 4, 1, 1, 1)
        self.Espesor=Entrada_con_unidades(unidades.Length, "Thickness")
        gridLayout_Costos.addWidget(self.Espesor,4,2,1,1)
        gridLayout_Costos.addWidget(QtWidgets.QLabel(self.tr("Cabeza:", None)), 5, 1, 1, 1)
        self.Cabeza=QtWidgets.QComboBox()
        self.Cabeza.addItem(self.tr("Elipsoidal", None))
        self.Cabeza.addItem(self.tr("Semiesférico", None))
        self.Cabeza.addItem(self.tr("Bumped", None))
        self.Cabeza.addItem(self.tr("Liso", None))
        self.Cabeza.currentIndexChanged.connect(self.calcularCostos)
        gridLayout_Costos.addWidget(self.Cabeza, 5, 2, 1, 1)
        gridLayout_Costos.addWidget(QtWidgets.QLabel(self.tr("Espesor (cabeza):", None)), 6, 1, 1, 1)
        self.EspesorCabeza=Entrada_con_unidades(unidades.Length, "Thickness")
        gridLayout_Costos.addWidget(self.EspesorCabeza,6,2,1,1)
        gridLayout_Costos.addWidget(QtWidgets.QLabel(self.tr("Longitud reborde recto:", None)), 7, 1, 1, 1)
        self.LongitudReborde=Entrada_con_unidades(unidades.Length)
        gridLayout_Costos.addWidget(self.LongitudReborde,7,2,1,1)
        gridLayout_Costos.addWidget(QtWidgets.QLabel(self.tr("Volumen:", None)), 6, 4, 1, 1)
        self.Volumen=Entrada_con_unidades(unidades.Volume, "VolLiq", readOnly=True)
        gridLayout_Costos.addWidget(self.Volumen,6,5,1,1)
        gridLayout_Costos.addWidget(QtWidgets.QLabel(self.tr("Peso:", None)), 7, 4, 1, 1)
        self.Peso=Entrada_con_unidades(unidades.Mass, readOnly=True)
        gridLayout_Costos.addWidget(self.Peso,7,5,1,1)
        gridLayout_Costos.addItem(QtWidgets.QSpacerItem(10,10,QtWidgets.QSizePolicy.Policy.Fixed,QtWidgets.QSizePolicy.Policy.Fixed),2,3,6,1)
        gridLayout_Costos.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Policy.Expanding,QtWidgets.QSizePolicy.Policy.Expanding),8,0,1,6)

        self.Costos=costIndex.CostData(1.7, 3)
        self.Costos.valueChanged.connect(self.calcularCostos)
        gridLayout_Costos.addWidget(self.Costos,9,1,2,5)

        gridLayout_Costos.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Policy.Expanding,QtWidgets.QSizePolicy.Policy.Expanding),11,0,1,6)
        self.groupBox_Costos = QtWidgets.QGroupBox(self.tr("Costos calculados", None))
        gridLayout_Costos.addWidget(self.groupBox_Costos,12,1,1,5)
        gridLayout_5 = QtWidgets.QGridLayout(self.groupBox_Costos)
        gridLayout_5.addWidget(QtWidgets.QLabel(self.tr("Coste Adquisición:", None)),0,1,1,1)
        self.C_adq=Entrada_con_unidades(unidades.Currency, retornar=False, readOnly=True, tolerancia=8, decimales=2)
        gridLayout_5.addWidget(self.C_adq,0,2,1,1)
        gridLayout_5.addWidget(QtWidgets.QLabel(self.tr("Coste Instalación:", None)),1,1,1,1)
        self.C_inst=Entrada_con_unidades(unidades.Currency, retornar=False, readOnly=True, tolerancia=8, decimales=2)
        gridLayout_5.addWidget(self.C_inst,1,2,1,1)
        gridLayout_Costos.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Policy.Expanding,QtWidgets.QSizePolicy.Policy.Expanding),13,0,1,6)

        #Pestaña salida
        self.Salida= UI_corriente.Ui_corriente(readOnly=True)
        self.tabWidget.insertTab(2, self.Salida,self.tr("Salida", None))

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
    app = QtWidgets.QApplication(sys.argv)
    agua=Corriente(300, 1, 3600, Mezcla([62], [1]))
    dialogo = UI_equipment(agua)
    dialogo.show()
    sys.exit(app.exec())
