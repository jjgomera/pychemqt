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
###                                                  Diálogo de definición de filtros, UI_filter                                                     ###
#######################################################################

from tools.qt import QtCore, QtWidgets, tr


from equipment.liquid_solid import Filter
from UI import UI_corriente
from equipment import parents
from lib import unidades
from tools import costIndex
from UI.widgets import Entrada_con_unidades


class UI_equipment(parents.UI_equip):
    """Diálogo de definición de filtros por presión o a vación para la separación de sólidos de corrientes líquidas"""
    def __init__(self, entrada=None, parent=None):
        """entrada: Parametro opcional de clase corriente que indica la corriente de entrada en kla tubería"""
        super(UI_equipment, self).__init__(Filter, entrada=False, parent=parent)
        self.entrada=entrada

        #Pestaña entrada
        self.Entrada= UI_corriente.Ui_corriente(entrada)
        self.Entrada.Changed.connect(self.cambiar_entrada)
        self.tabWidget.insertTab(0, self.Entrada, tr("equipment", "Entrada", None))

        #Pestaña calculo
        gridLayout_Calculo = QtWidgets.QGridLayout(self.tabCalculo)

        #Pestaña costos
        gridLayout_Costos = QtWidgets.QGridLayout(self.tabCostos)
        gridLayout_Costos.addWidget(QtWidgets.QLabel(tr("equipment", "Tipo:", None)), 1, 0, 1, 1)
        self.tipo=QtWidgets.QComboBox()
        self.tipo.addItem(tr("equipment", "Rotary vacuum belt discharge", None))
        self.tipo.addItem(tr("equipment", "Rotary vacuum drum scraper discharge", None))
        self.tipo.addItem(tr("equipment", "Rotary vacuum disk", None))
        self.tipo.addItem(tr("equipment", "Horizontal vacuum belt", None))
        self.tipo.addItem(tr("equipment", "Pressure leaf", None))
        self.tipo.addItem(tr("equipment", "Plate and frame", None))
        self.tipo.currentIndexChanged.connect(self.calcularCostos)
        gridLayout_Costos.addWidget(self.tipo, 1, 1, 1, 3)
        gridLayout_Costos.addItem(QtWidgets.QSpacerItem(10,10,QtWidgets.QSizePolicy.Policy.Fixed,QtWidgets.QSizePolicy.Policy.Fixed),2,0,1,2)

        self.Costos=costIndex.CostData(1.3, 2)
        self.Costos.valueChanged.connect(self.calcularCostos)
        gridLayout_Costos.addWidget(self.Costos,4,0,2,5)

        gridLayout_Costos.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Policy.Expanding,QtWidgets.QSizePolicy.Policy.Expanding),6,0,1,6)
        gridLayout_Costos.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Policy.Expanding,QtWidgets.QSizePolicy.Policy.Expanding),10,0,1,6)
        self.groupBox_Costos = QtWidgets.QGroupBox(tr("equipment", "Costos calculados", None))
        gridLayout_Costos.addWidget(self.groupBox_Costos,7,0,1,4)
        gridLayout_5 = QtWidgets.QGridLayout(self.groupBox_Costos)
        gridLayout_5.addWidget(QtWidgets.QLabel(tr("equipment", "Coste Adquisición:", None)),1,1)
        self.C_adq=Entrada_con_unidades(unidades.Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_adq,1,2)
        gridLayout_5.addWidget(QtWidgets.QLabel(tr("equipment", "Coste Instalación:", None)),2,1)
        self.C_inst=Entrada_con_unidades(unidades.Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_inst,2,2)


        #Pestaña salida
        self.SalidaGas= UI_corriente.Ui_corriente(readOnly=True)
        self.SalidaSolido= UI_corriente.Ui_corriente(readOnly=True)
        self.Salida.addTab(self.SalidaGas,tr("equipment", "Gas filtrado", None))
        self.Salida.addTab(self.SalidaSolido,tr("equipment", "Sólidos recogidos", None))

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
    from lib.corriente import Corriente, Mezcla, Solid
    app = QtWidgets.QApplication(sys.argv)
    distribucion=[[17.5, 0.02],
                                [22.4, 0.03],
                                [26.2,  0.05],
                                [31.8,  0.1],
                                [37, 0.1],
                                [42.4, 0.1],
                                [48, 0.1],
                                [54, 0.1],
                                [60, 0.1],
                                [69, 0.1],
                                [81.3, 0.1],
                                [96.5, 0.05],
                                [109, 0.03],
                                [127, 0.02]]

    solido=Solid([64], [138718], distribucion)
    agua=Corriente(300, 1, 3600, Mezcla([62], [1]), solido)
    dialogo = UI_equipment(agua)
    dialogo.show()
    sys.exit(app.exec())
