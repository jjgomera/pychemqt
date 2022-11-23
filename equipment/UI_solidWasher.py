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



############################################################################
###                                       Diálogo de definición de unidades de secado de sólidos, UI_dryer                                      ###
############################################################################

from tools.qt import QtCore, QtWidgets


from equipment.gas_solid_liquid import Dryer
from lib import unidades, config
from UI import UI_corriente
from equipment import parents
from UI.widgets import Entrada_con_unidades
from tools import costIndex


class UI_equipment(parents.UI_equip):
    """Dialogo de definición de unidades de secado de sólidos"""
    def __init__(self, entradaSolido=None, entradaAire=None,  parent=None):
        """entrada: Parametro opcional de clase corriente que indica la corriente de entrada"""
        super(UI_equipment, self).__init__(Dryer, parent=parent)
        self.entradaSolido=entradaSolido
        self.entradaAire=entradaAire
        self.Equipment=Dryer()

        #Pestaña entrada
        self.EntradaSolido= UI_corriente.Ui_corriente(self.entradaSolido)
        self.EntradaSolido.Changed.connect(self.cambiar_entrada)
        self.Entrada.addTab(self.EntradaSolido,QtWidgets.QApplication.translate("equipment", "Solido", None))
        self.EntradaAire= UI_corriente.Ui_psychrometry(self.entradaAire)
        self.EntradaAire.Changed.connect(self.cambiar_aire)
        self.Entrada.addTab(self.EntradaAire,QtWidgets.QApplication.translate("equipment", "Aire", None))

        #Pestaña calculo
        gridLayout_Calculo = QtWidgets.QGridLayout(self.tabCalculo)
        gridLayout_Calculo.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("equipment", "Tipo de cálculo:", None)), 1, 1)
        self.TipoCalculo=QtWidgets.QComboBox()
        self.TipoCalculo.addItem(QtWidgets.QApplication.translate("equipment", "Cálculo, conocido el flujo de vapor, calcular las corrientes de salida", None))
        self.TipoCalculo.addItem(QtWidgets.QApplication.translate("equipment", "Diseño, calcular el flujo de aire necesario", None))
        self.TipoCalculo.currentIndexChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.TipoCalculo, 1, 2, 1, 4)
        gridLayout_Calculo.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Policy.Fixed,QtWidgets.QSizePolicy.Policy.Fixed),2,1,1,6)
        gridLayout_Calculo.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("equipment", "Humedad relativa en el aire:", None)), 3, 1, 1, 1)
        self.HumedadAire=Entrada_con_unidades(float, max=1, spinbox=True, step=0.01)
        self.HumedadAire.valueChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.HumedadAire,3,2,1,1)
        gridLayout_Calculo.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("equipment", "Humedad residual del sólido:", None)), 4, 1, 1, 1)
        self.HumedadSolido=Entrada_con_unidades(float, max=1., spinbox=True, step=0.01, textounidad=unidades.Mass(None).text()+"/"+unidades.Mass(None).text())
        self.HumedadSolido.valueChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.HumedadSolido,4,2,1,1)
        gridLayout_Calculo.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("equipment", "Temperatura del sólido a la salida:", None)), 5, 1, 1, 1)
        self.temperatura=Entrada_con_unidades(unidades.Temperature)
        self.temperatura.valueChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.temperatura,5,2,1,1)
        gridLayout_Calculo.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("equipment", "Intercambio de calor:", None)), 6, 1, 1, 1)
        self.Heat=Entrada_con_unidades(unidades.Power)
        self.Heat.valueChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.Heat,6,2,1,1)
        gridLayout_Calculo.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("equipment", "Pérdida de presión:", None)), 7, 1, 1, 1)
        self.DeltaP=Entrada_con_unidades(unidades.Pressure)
        self.DeltaP.valueChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.DeltaP,7,2,1,1)

        gridLayout_Calculo.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Policy.Expanding,QtWidgets.QSizePolicy.Policy.Expanding),8,1,1,6)
        self.groupBox_Calculo = QtWidgets.QGroupBox(QtWidgets.QApplication.translate("equipment", "Datos calculados", None))
        gridLayout_Calculo.addWidget(self.groupBox_Calculo,9,1,1,5)
        gridLayout_1 = QtWidgets.QGridLayout(self.groupBox_Calculo)
        gridLayout_1.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("equipment", "Temperatura a la salida:", None)), 1, 1, 1, 1)
        self.temperaturaCalculada=Entrada_con_unidades(unidades.Temperature, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.temperaturaCalculada,1,2,1,1)
        gridLayout_1.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("equipment", "Caudal de aire:", None)),2,1)
        self.caudalVolumetrico=Entrada_con_unidades(unidades.VolFlow, "QGas", retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.caudalVolumetrico,2,2,1,1)
        gridLayout_1.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("equipment", "Humedad del aire:", None)), 3, 1)
        self.HumedadCalculada=Entrada_con_unidades(float, readOnly=True, textounidad="%")
        gridLayout_1.addWidget(self.HumedadCalculada,3,2,1,1)

        gridLayout_Calculo.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Policy.Expanding,QtWidgets.QSizePolicy.Policy.Expanding),11,1,1,6)

        #Pestaña salida
        self.SalidaSolido= UI_corriente.Ui_corriente(readOnly=True)
        self.Salida.addTab(self.SalidaSolido,QtWidgets.QApplication.translate("equipment", "Sólido secado", None))
        self.SalidaAire= UI_corriente.Ui_psychrometry(readOnly=True)
        self.Salida.addTab(self.SalidaAire,QtWidgets.QApplication.translate("equipment", "Aire", None))


    def cambiar_entrada(self, corriente):
        self.entradaSolido=corriente
        self.calculo()

    def cambiar_aire(self, punto):
        self.entradaAire=punto
        self.calculo()

    def calculo(self):
        if self.todos_datos():
            self.Equipment(entrada=self.entrada, calculo=0, modelo=self.Modelo.currentIndex(), anchura=self.anchura.value, altura=self.altura.value, longitud=self.longitud.value)
            self.rellenoSalida(1)

    def rellenoSalida(self, estado=1, texto=""):
        self.caudalVolumetrico.setValue(self.entrada.caudal_volumetrico)
        self.velocidadGasCalculada.setValue(self.Equipment.Vgas)
        self.rendimientoCalculado.setValue(self.Equipment.rendimiento)
        self.alturaCalculada.setValue(self.Equipment.H)
        self.anchuraCalculada.setValue(self.Equipment.B)
        self.longitudCalculada.setValue(self.Equipment.L)
        self.SalidaGas.rellenar(self.Equipment.SalidaAire)
        self.SalidaSolido.rellenar(self.Equipment.SalidaSolido)
        self.status.setState(estado, texto)

    def todos_datos(self):
        return self.EntradaSolido.todos_datos() and self.EntradaAire.todos_datos()


if __name__ == "__main__":
    import sys
    from lib.corriente import Mezcla, Punto_Psicrometrico, Corriente, Solid
    app = QtWidgets.QApplication(sys.argv)
    distribucion=[[96.5, 0.02],
                        [105, 0.05],
                        [110,  0.1],
                        [118, 0.15],
                        [125, 0.25],
                        [130, 0.2],
                        [140, 0.15],
                        [150, 0.05],
                        [170, 0.03]]

    solido=Solid([638], [100], distribucion)
    entradaSolido=Corriente(300, 1, 1,  Mezcla([62], [1.]), solido)
    aire=Punto_Psicrometrico(caudal=1000, tdb=300, H=0)
    dialogo = UI_equipment(entradaSolido, aire)
    dialogo.show()
    sys.exit(app.exec())
