#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2022, Juan José Gómez Romera <jjgomera@gmail.com>

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
# Diálogo de definición de tamices, UI_screen
#######################################################################

import os, sys
path=os.path.dirname("/home/jjgomera/pychemqt/")
sys.path.append(path)

from PyQt5 import QtWidgets


from lib import unidades
from tools.costIndex import CostData
from equipment.parents import UI_equip
from UI.widgets import Entrada_con_unidades
from UI import UI_corriente
from equipment.solids import Screen


class UI_equipment(UI_equip):
    """Diálogo de definición de tamices de sólidos"""
    Equipment=Screen()
    def __init__(self, equipment=None, parent=None):
        """
        equipment: instancia de equipo inicial
        """
        super(UI_equipment, self).__init__(Screen, entrada=False, parent=parent)

        #Pestaña entrada
#        self.Entrada= UI_corriente.Ui_corriente(entrada)
#        self.Entrada.Changed.connect(self.cambiar_entrada)
#        self.tabWidget.insertTab(0, self.Entrada, QtGui.QApplication.translate("equipment", "Entrada", None, QtGui.QApplication.UnicodeUTF8))

        #Pestaña calculo
        gridLayout_Calculo = QtWidgets.QGridLayout(self.tabCalculo)

        #Pestaña costos
        gridLayout_Costos = QtWidgets.QGridLayout(self.tabCostos)
        gridLayout_Costos.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("equipment", "Area:", None)), 1, 1, 1, 1)
        self.Area=Entrada_con_unidades(unidades.Area)
        self.Area.valueChanged.connect(self.calcularCostos)
        gridLayout_Costos.addWidget(self.Area, 1, 2, 1, 1)
        gridLayout_Costos.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Fixed,QtWidgets.QSizePolicy.Fixed),2,0,1,2)

        self.Costos=CostData(self.Equipment)
        self.Costos.valueChanged.connect(self.changeParamsCoste)
        gridLayout_Costos.addWidget(self.Costos,4,1,2,5)

        gridLayout_Costos.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Expanding),6,0,1,6)
        gridLayout_Costos.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Expanding),10,0,1,6)
        self.groupBox_Costos = QtWidgets.QGroupBox(QtWidgets.QApplication.translate("equipment", "Costos calculados", None))
        gridLayout_Costos.addWidget(self.groupBox_Costos,7,0,1,6)
        gridLayout_5 = QtWidgets.QGridLayout(self.groupBox_Costos)
        gridLayout_5.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("equipment", "Coste Adquisición:", None)),0,1,1,1)
        self.C_adq=Entrada_con_unidades(unidades.Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_adq,0,2,1,1)
        gridLayout_5.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("equipment", "Coste Instalación:", None)),1,1,1,1)
        self.C_inst=Entrada_con_unidades(unidades.Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_inst,1,2,1,1)

        #Pestaña salida
        self.SalidaGas= UI_corriente.Ui_corriente(readOnly=True)
        self.SalidaSolido= UI_corriente.Ui_corriente(readOnly=True)
        self.Salida.addTab(self.SalidaGas,QtWidgets.QApplication.translate("equipment", "Gas filtrado", None))
        self.Salida.addTab(self.SalidaSolido,QtWidgets.QApplication.translate("equipment", "Sólidos recogidos", None))

        self.tabWidget.setCurrentIndex(0)


    def cambiar_entrada(self, corriente):
        selfentrada=corriente
        self.calculo()

    def calculo(self):
        if self.todos_datos():

            self.rellenoSalida()

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
    agua=Corriente(T=300, P=1e5, caudal=3600)
    dialogo = UI_equipment()
    dialogo.show()
    sys.exit(app.exec_())
