#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2025, Juan José Gómez Romera <jjgomera@gmail.com>

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
###                                         Diálogo de definición de centrifugas, UI_centrifuge                                             ###
#######################################################################


from tools.qt import QtCore, QtWidgets

from equipment.liquid_solid import Centrifuge
from UI import UI_corriente
from equipment import parents
from lib.corriente import Corriente, Solid
from lib import unidades
from tools import costIndex
from UI.widgets import Entrada_con_unidades


class UI_equipment(parents.UI_equip):
    """Diálogo de definición de tamices de sólidos"""
    def __init__(self, entrada=None, parent=None):
        """entrada: Parametro opcional de clase corriente que indica la corriente de entrada en kla tubería"""
        super(UI_equipment, self).__init__(Centrifuge, entrada=False, parent=parent)
        self.entrada=entrada

        #Pestaña entrada
        self.Entrada= UI_corriente.Ui_corriente(entrada)
        self.Entrada.Changed.connect(self.cambiar_entrada)
        self.tabWidget.insertTab(0, self.Entrada,self.tr("Entrada", None))

        #Pestaña calculo
        gridLayout_Calculo = QtWidgets.QGridLayout(self.tabCalculo)

        #Pestaña costos
        gridLayout_Costos = QtWidgets.QGridLayout(self.tabCostos)
        gridLayout_Costos.addWidget(QtWidgets.QLabel(self.tr("Tipo:", None)), 1, 1, 1, 1)
        self.tipo=QtWidgets.QComboBox()
        self.tipo.addItem(self.tr("Proceso inorgánico", None))
        self.tipo.addItem(self.tr("Proceso orgánico", None))
        self.tipo.currentIndexChanged.connect(self.mostrarSubclasificacion)
        self.tipo.currentIndexChanged.connect(self.calcularCostos)
        gridLayout_Costos.addWidget(self.tipo, 1, 2, 1, 1)
        gridLayout_Costos.addWidget(QtWidgets.QLabel(self.tr("Material:", None)), 2, 1, 1, 1)
        self.materialInorganico=QtWidgets.QComboBox()
        self.materialInorganico.addItem(self.tr("Acero al carbón", None))
        self.materialInorganico.addItem(self.tr("Acero inoxidable 316", None))
        self.materialInorganico.addItem(self.tr("Monel", None))
        self.materialInorganico.addItem(self.tr("Níquel", None))
        self.materialInorganico.currentIndexChanged.connect(self.calcularCostos)
        gridLayout_Costos.addWidget(self.materialInorganico, 2, 2, 1, 1)
        gridLayout_Costos.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Policy.Fixed,QtWidgets.QSizePolicy.Policy.Fixed),3,0,1,2)
        self.materialOrganico=QtWidgets.QComboBox()
        self.materialOrganico.addItem(self.tr("Acero al carbón", None))
        self.materialOrganico.addItem(self.tr("Acero inoxidable 316", None))
        self.materialOrganico.addItem(self.tr("Monel", None))
        self.materialOrganico.addItem(self.tr("Níquel", None))
        self.materialOrganico.addItem(self.tr("Hastelloy", None))
        self.materialOrganico.currentIndexChanged.connect(self.calcularCostos)
        gridLayout_Costos.addWidget(self.materialOrganico, 2, 2, 1, 1)
        gridLayout_Costos.addItem(QtWidgets.QSpacerItem(10,10,QtWidgets.QSizePolicy.Policy.Fixed,QtWidgets.QSizePolicy.Policy.Fixed),3,0,1,6)
        gridLayout_Costos.addWidget(QtWidgets.QLabel(self.tr("Caudal calculado:", None)), 4, 1, 1, 1)
        self.caudalcalculado=Entrada_con_unidades(unidades.MassFlow, readOnly=True, retornar=False)
        gridLayout_Costos.addWidget(self.caudalcalculado,4,2,1,1)
        gridLayout_Costos.addWidget(QtWidgets.QLabel(self.tr("Caudal de diseño:", None)), 5, 1, 1, 1)
        self.caudaldiseno=Entrada_con_unidades(unidades.MassFlow)
        gridLayout_Costos.addWidget(self.caudaldiseno,5,2,1,1)
        gridLayout_Costos.addItem(QtWidgets.QSpacerItem(10,10,QtWidgets.QSizePolicy.Policy.Fixed,QtWidgets.QSizePolicy.Policy.Fixed),6,0,1,6)

        self.Costos=costIndex.CostData(1.3, 2)
        self.Costos.valueChanged.connect(self.calcularCostos)
        gridLayout_Costos.addWidget(self.Costos,7,1,2,2)

        gridLayout_Costos.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Policy.Expanding,QtWidgets.QSizePolicy.Policy.Expanding),9,0,1,6)
        gridLayout_Costos.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Policy.Expanding,QtWidgets.QSizePolicy.Policy.Expanding),11,0,1,6)
        self.groupBox_Costos = QtWidgets.QGroupBox(self.tr("Costos calculados", None))
        gridLayout_Costos.addWidget(self.groupBox_Costos,10,0,1,6)
        gridLayout_5 = QtWidgets.QGridLayout(self.groupBox_Costos)
        gridLayout_5.addWidget(QtWidgets.QLabel(self.tr("Coste Adquisición:", None)),0,1,1,1)
        self.C_adq=Entrada_con_unidades(unidades.Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_adq,0,2,1,1)
        gridLayout_5.addWidget(QtWidgets.QLabel(self.tr("Coste Instalación:", None)),1,1,1,1)
        self.C_inst=Entrada_con_unidades(unidades.Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_inst,1,2,1,1)

        self.tabWidget.setCurrentIndex(0)
        self.mostrarSubclasificacion(0)

    def mostrarSubclasificacion(self, ind):
        self.materialInorganico.setVisible(not ind)
        self.materialOrganico.setVisible(ind)

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
