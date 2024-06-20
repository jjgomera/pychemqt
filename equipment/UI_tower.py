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
###                                               Diálogo de definición de tuberías, UI_pipe                                                     ###
#######################################################################

from tools.qt import QtCore, QtWidgets


from equipment.distillation import Tower
from UI import UI_corriente
from equipment import parents
from lib import unidades, config
from tools import costIndex
from UI.widgets import Entrada_con_unidades


class UI_equipment(parents.UI_equip):
    """Diálogo de definición de tuberías"""
    def __init__(self, entrada=None, parent=None):
        """entrada: Parametro opcional de clase corriente que indica la corriente de entrada en kla tubería"""
        super(UI_equipment, self).__init__(Tower, parent=parent)
        self.entrada=entrada

#        #Pestaña entrada
#        self.Entrada= UI_corriente.Ui_corriente(entrada)
#        self.Entrada.Changed.connect(self.cambiar_entrada)
#        self.tabWidget.addTab(self.Entrada, self.tr("Entrada", None, QtGui.QApplication.UnicodeUTF8))

        #Pestaña calculo
        gridLayout_Calculo = QtWidgets.QGridLayout(self.tabCalculo)

        #Pestaña costos
        gridLayout_Costos = QtWidgets.QGridLayout(self.tabCostos)
        gridLayout_Costos.addWidget(QtWidgets.QLabel(self.tr("Proceso:", None)), 1, 1, 1, 1)
        self.proceso=QtWidgets.QComboBox()
        self.proceso.addItem(self.tr("Destilación", None))
        self.proceso.addItem(self.tr("Absorción", None))
        self.proceso.currentIndexChanged.connect(self.calcularCostos)
        gridLayout_Costos.addWidget(self.proceso, 1, 2, 1, 1)
        gridLayout_Costos.addWidget(QtWidgets.QLabel(self.tr("Tipo de columna:", None)), 2, 1, 1, 1)
        self.tipo=QtWidgets.QComboBox()
        self.tipo.addItem(self.tr("De pisos", None))
        self.tipo.addItem(self.tr("De relleno", None))
        self.tipo.currentIndexChanged.connect(self.mostrarSubclasificacion)
        self.tipo.currentIndexChanged.connect(self.calcularCostos)
        gridLayout_Costos.addWidget(self.tipo, 2, 2, 1, 1)
        gridLayout_Costos.addWidget(QtWidgets.QLabel(self.tr("Material:", None)), 3, 1, 1, 1)
        self.material=QtWidgets.QComboBox()
        self.material.addItem(self.tr("Acero al carbon", None))
        self.material.addItem(self.tr("Acero inoxidable 304", None))
        self.material.addItem(self.tr("Acero inoxidable 316", None))
        self.material.addItem(self.tr("Carpenter 20CB-3", None))
        self.material.addItem(self.tr("Niquel 200", None))
        self.material.addItem(self.tr("Monel 400", None))
        self.material.addItem(self.tr("Inconel 600", None))
        self.material.addItem(self.tr("Incoloy 825", None))
        self.material.addItem(self.tr("Titanio", None))
        self.material.currentIndexChanged.connect(self.calcularCostos)
        gridLayout_Costos.addWidget(self.material, 3, 2, 1, 1)

        gridLayout_Costos.addItem(QtWidgets.QSpacerItem(30,30,QtWidgets.QSizePolicy.Policy.Fixed,QtWidgets.QSizePolicy.Policy.Fixed),1,3,5,1)

        self.groupBox_Pisos = QtWidgets.QGroupBox(self.tr("Torre de pisos", None))
        gridLayout_Costos.addWidget(self.groupBox_Pisos,1,4,4,2)
        gridLayout_1 = QtWidgets.QGridLayout(self.groupBox_Pisos)
        gridLayout_1.addWidget(QtWidgets.QLabel(self.tr("Tipo:", None)), 1, 1, 1, 1)
        self.tipoPisos=QtWidgets.QComboBox()
        self.tipoPisos.addItem(self.tr("De válvula", None))
        self.tipoPisos.addItem(self.tr("De rejilla", None))
        self.tipoPisos.addItem(self.tr("De borboteo", None))
        self.tipoPisos.addItem(self.tr("De tamiz", None))
        self.tipoPisos.currentIndexChanged.connect(self.calcularCostos)
        gridLayout_1.addWidget(self.tipoPisos, 1, 2, 1, 1)
        gridLayout_1.addWidget(QtWidgets.QLabel(self.tr("Material:", None)), 2, 1, 1, 1)
        self.materialPisos=QtWidgets.QComboBox()
        self.materialPisos.addItem(self.tr("Acero al carbon", None))
        self.materialPisos.addItem(self.tr("Acero inoxidable 304", None))
        self.materialPisos.addItem(self.tr("Acero inoxidable 316", None))
        self.materialPisos.addItem(self.tr("Carpenter 20CB-3", None))
        self.materialPisos.addItem(self.tr("Niquel 200", None))
        self.materialPisos.addItem(self.tr("Monel 400", None))
        self.materialPisos.addItem(self.tr("Inconel 600", None))
        self.materialPisos.addItem(self.tr("Incoloy 825", None))
        self.materialPisos.addItem(self.tr("Titanio", None))
        self.materialPisos.currentIndexChanged.connect(self.calcularCostos)
        gridLayout_1.addWidget(self.materialPisos, 2, 2, 1, 1)
        gridLayout_1.addItem(QtWidgets.QSpacerItem(10,10,QtWidgets.QSizePolicy.Policy.Fixed,QtWidgets.QSizePolicy.Policy.Fixed),3,1,1,2)
        gridLayout_1.addWidget(QtWidgets.QLabel(self.tr("Diametro:", None)), 4, 1, 1, 1)
        self.diametroPisos=Entrada_con_unidades(unidades.Length)
        gridLayout_1.addWidget(self.diametroPisos,4,2,1,1)
        gridLayout_1.addWidget(QtWidgets.QLabel(self.tr("Número:", None)), 5, 1, 1, 1)
        self.NumeroPisos=Entrada_con_unidades(int, spinbox=True, min=1, step=1, width=50)
        gridLayout_1.addWidget(self.NumeroPisos,5,2,1,1)
        gridLayout_1.addItem(QtWidgets.QSpacerItem(10,10,QtWidgets.QSizePolicy.Policy.Fixed,QtWidgets.QSizePolicy.Policy.Fixed),6,1,1,2)

        self.groupBox_relleno = QtWidgets.QGroupBox(self.tr("Torre de relleno", None))
        gridLayout_Costos.addWidget(self.groupBox_relleno,1, 4, 4, 2)
        gridLayout_2 = QtWidgets.QGridLayout(self.groupBox_relleno)
        gridLayout_2.addWidget(QtWidgets.QLabel(self.tr("Volumen:", None)), 1, 1, 1, 1)
        self.VolumenRelleno=Entrada_con_unidades(unidades.Volume, "VolLiq")
        gridLayout_2.addWidget(self.VolumenRelleno,1,2,1,1)
        gridLayout_2.addWidget(QtWidgets.QLabel(self.tr("Coste unitario:", None)),2,1,1,1)
        self.C_unit_relleno=Entrada_con_unidades(unidades.Currency, retornar=False, textounidad="%s / %s" % (unidades.Currency(None).text(), unidades.Volume(None).text("VolLiq")))
        gridLayout_2.addWidget(self.C_unit_relleno,2,2,1,1)
        gridLayout_2.addItem(QtWidgets.QSpacerItem(10,10,QtWidgets.QSizePolicy.Policy.Fixed,QtWidgets.QSizePolicy.Policy.Fixed),3,1,1,2)

        gridLayout_Costos.addWidget(QtWidgets.QLabel(self.tr("Diametro:", None)),5,1,1,1)
        self.Dc=Entrada_con_unidades(unidades.Length)
        gridLayout_Costos.addWidget(self.Dc,5,2,1,2)
        gridLayout_Costos.addWidget(QtWidgets.QLabel(self.tr("Altura:", None)), 6, 1, 1, 1)
        self.Hc=Entrada_con_unidades(unidades.Length)
        gridLayout_Costos.addWidget(self.Hc,6,2,1,2)
        gridLayout_Costos.addWidget(QtWidgets.QLabel(self.tr("Espesor (Tapa):", None)), 7, 1, 1, 1)
        self.EspesorSuperior=Entrada_con_unidades(unidades.Length, "Thickness")
        gridLayout_Costos.addWidget(self.EspesorSuperior,7,2,1,2)
        gridLayout_Costos.addWidget(QtWidgets.QLabel(self.tr("Espesor (Fondo):", None)), 8, 1, 1, 1)
        self.EspesorInferior=Entrada_con_unidades(unidades.Length, "Thickness")
        gridLayout_Costos.addWidget(self.EspesorInferior,8,2,1,2)
        gridLayout_Costos.addWidget(QtWidgets.QLabel(self.tr("Densidad:", None)), 9, 1, 1, 1)
        self.EspesorInferior=Entrada_con_unidades(unidades.Density, "DenLiq")
        gridLayout_Costos.addWidget(self.EspesorInferior,9,2,1,2)

        self.Costos=costIndex.CostData(3, 3)
        self.Costos.valueChanged.connect(self.calcularCostos)
        gridLayout_Costos.addWidget(self.Costos,10,1,2,4)

        gridLayout_Costos.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Policy.Expanding,QtWidgets.QSizePolicy.Policy.Expanding),12,1,1,6)

        self.groupBox_Costos = QtWidgets.QGroupBox(self.tr("Costos calculados", None))
        gridLayout_Costos.addWidget(self.groupBox_Costos,13,1,1,5)
        gridLayout_5 = QtWidgets.QGridLayout(self.groupBox_Costos)
        gridLayout_5.addWidget(QtWidgets.QLabel(self.tr("Coste Pisos:", None)),0,1,1,1)
        self.C_pisos=Entrada_con_unidades(unidades.Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_pisos,0,2,1,1)
        gridLayout_5.addWidget(QtWidgets.QLabel(self.tr("Coste Carcasa:", None)),1,1,1,1)
        self.C_carcasa=Entrada_con_unidades(unidades.Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_carcasa,1,2,1,1)
        gridLayout_5.addWidget(QtWidgets.QLabel(self.tr("Coste Accesorios:", None)),2,1,1,1)
        self.C_accesorios=Entrada_con_unidades(unidades.Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_accesorios,2,2,1,1)
        gridLayout_5.addWidget(QtWidgets.QLabel(self.tr("Coste Columna:", None)),0,4,1,1)
        self.C_columna=Entrada_con_unidades(unidades.Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_columna,0,5,1,1)
        gridLayout_5.addWidget(QtWidgets.QLabel(self.tr("Coste Adquisición:", None)),1,4,1,1)
        self.C_adq=Entrada_con_unidades(unidades.Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_adq,1,5,1,1)
        gridLayout_5.addWidget(QtWidgets.QLabel(self.tr("Coste Instalación:", None)),2,4,1,1)
        self.C_inst=Entrada_con_unidades(unidades.Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_inst,2,5,1,1)


#        #Pestaña salida
#        self.pSalida = QtGui.QTabWidget()
#        self.tabWidget.addTab(self.pSalida,self.tr("Salida", None, QtGui.QApplication.UnicodeUTF8))


        self.mostrarSubclasificacion(0)


    def mostrarSubclasificacion(self, ind):
        self.groupBox_Pisos.setVisible(not ind)
        self.groupBox_relleno.setVisible(ind)

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
            self.ShellTube.Coste(self.factorInstalacion.value(), 0, self.tipo.currentIndex(), self.material.currentIndex())
            self.C_adq.setValue(self.ShellTube.C_adq.config())
            self.C_inst.setValue(self.ShellTube.C_inst.config())

if __name__ == "__main__":
    import sys
    from lib.corriente import Corriente, Mezcla
    app = QtWidgets.QApplication(sys.argv)
    agua=Corriente(T=300, P=101325, caudalMasico=3600, ids=[62], fraccion=[1])
    dialogo = UI_equipment(agua)
    dialogo.show()
    sys.exit(app.exec())
