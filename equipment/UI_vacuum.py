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



##########################################################################
###                                              Diálogo de definición de filtros de mangas, UI_baghouse                                       ###
##########################################################################

from tools.qt import QtCore, QtWidgets


from equipment.gas_solid import Baghouse
from lib import unidades
from lib.utilities import representacion
from UI import UI_corriente
from equipment import parents
from UI.delegate import CellEditor
from UI.widgets import Entrada_con_unidades


class UI_equipment(parents.UI_equip):
    """Diálogo de definición de filtros de mangas"""
    def __init__(self, entrada=None, parent=None):
        """entrada: Parametro opcional de clase corriente que indica la corriente de entrada"""
        super(UI_equipment, self).__init__(Baghouse, entrada=False, parent=parent)
        self.entrada=entrada

        #Pestaña entrada
        self.Entrada= UI_corriente.Ui_corriente(entrada)
        self.Entrada.Changed.connect(self.cambiar_entrada)
        self.tabWidget.insertTab(0, self.Entrada,self.tr("Entrada", None))

        #Pestaña definición rendimientos
        self.Rendimientos= QtWidgets.QTableWidget(1, 2)
        self.Rendimientos.setItemDelegateForColumn(1, CellEditor(self))
        self.Rendimientos.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.Rendimientos.verticalHeader().hide()
        self.Rendimientos.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.AllEditTriggers)
        if self.entrada:
            self.rellenarTablaRendimientos()
        self.rendimientos=[]
        self.Rendimientos.cellChanged.connect(self.cambiarRendimientos)
        self.tabWidget.insertTab(1, self.Rendimientos,self.tr("Rendimientos", None))

        #Cálculo
        gridLayout_Calculo = QtWidgets.QGridLayout(self.tabCalculo)
        gridLayout_Calculo.addWidget(QtWidgets.QLabel(self.tr("Tipo de cálculo:", None)), 1, 1, 1, 1)
        self.TipoCalculo=QtWidgets.QComboBox()
        self.TipoCalculo.addItem(self.tr("Calcular caída de presión", None))
        self.TipoCalculo.addItem(self.tr("Calcular tiempo de filtración", None))
        self.TipoCalculo.addItem(self.tr("Calcular número de filtros", None))
        self.TipoCalculo.currentIndexChanged.connect(self.tipoCalculoCambiado)
        gridLayout_Calculo.addWidget(self.TipoCalculo, 1, 2, 1, 4)
        gridLayout_Calculo.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Policy.Fixed,QtWidgets.QSizePolicy.Policy.Fixed),2,1,1,6)
        gridLayout_Calculo.addWidget(QtWidgets.QLabel(self.tr("Nº de filtros:", None)), 3, 1, 1, 1)
        self.numFiltros=Entrada_con_unidades(int, spinbox=True, step=1, width=50, resaltado=True, min=1, start=1)
        self.numFiltros.valueChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.numFiltros,3,2,1,1)
        gridLayout_Calculo.addWidget(QtWidgets.QLabel(self.tr("Tiempo de filtración:", None)), 4, 1, 1, 1)
        self.tiempo=Entrada_con_unidades(unidades.Time, resaltado=True)
        self.tiempo.valueChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.tiempo,4,2,1,1)
        gridLayout_Calculo.addWidget(QtWidgets.QLabel(self.tr("Pérdida de presión:", None)), 5, 1, 1, 1)
        self.deltaP=Entrada_con_unidades(unidades.Pressure, retornar=False, readOnly=True)
        self.deltaP.valueChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.deltaP,5,2,1,1)
        gridLayout_Calculo.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Policy.Fixed,QtWidgets.QSizePolicy.Policy.Fixed),6,1,1,6)
        gridLayout_Calculo.addWidget(QtWidgets.QLabel(self.tr("Membranas por filtro:", None)), 7, 1, 1, 1)
        self.MembranaCelda=Entrada_con_unidades(int, spinbox=True, step=1, width=70, value=78, min=1)
        self.MembranaCelda.valueChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.MembranaCelda,7,2,1,1)
        gridLayout_Calculo.addWidget(QtWidgets.QLabel(self.tr("Diametro de membrana:", None)), 8, 1, 1, 1)
        self.Diametro=Entrada_con_unidades(unidades.Length, value=unidades.Length(0.5, "ft"))
        self.Diametro.valueChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.Diametro,8,2,1,1)
        gridLayout_Calculo.addWidget(QtWidgets.QLabel(self.tr("Area por membrana:", None)), 9, 1, 1, 1)
        self.Area=Entrada_con_unidades(unidades.Area, value=unidades.Area(16, "ft2"))
        self.Area.valueChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.Area,9,2,1,1)
        gridLayout_Calculo.addWidget(QtWidgets.QLabel(self.tr("Resistencia filtro:", None)), 7, 4, 1, 1)
        self.resistenciaFiltro=Entrada_con_unidades(float, spinbox=True, step=0.01, width=70, value=0.84, min=0)
        self.resistenciaFiltro.valueChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.resistenciaFiltro,7,5,1,1)
        gridLayout_Calculo.addWidget(QtWidgets.QLabel(self.tr("Resistencia torta:", None)), 8, 4, 1, 1)
        self.resistenciaTorta=Entrada_con_unidades(float, spinbox=True, step=0.01, width=70, value=0.1, min=0)
        self.resistenciaTorta.valueChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.resistenciaTorta,8,5,1,1)
        gridLayout_Calculo.addWidget(QtWidgets.QLabel(self.tr("Filtros en limpieza:", None)), 9, 4, 1, 1)
        self.Limpieza=Entrada_con_unidades(int, spinbox=True, step=1, width=70, value=1, min=0)
        self.Limpieza.valueChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.Limpieza,9,5,1,1)
        gridLayout_Calculo.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Policy.Expanding,QtWidgets.QSizePolicy.Policy.Expanding),10,1,1,6)

        self.groupBox_Calculo = QtWidgets.QGroupBox(self.tr("Datos calculados", None))
        gridLayout_Calculo.addWidget(self.groupBox_Calculo,11,1,1,5)
        self.gridLayout_1 = QtWidgets.QGridLayout(self.groupBox_Calculo)
        self.gridLayout_1.addWidget(QtWidgets.QLabel(self.tr("V<sub>gas</sub>:", None)),0,1,1,1)
        self.Vgas=Entrada_con_unidades(unidades.Speed, retornar=False, readOnly=True)
        self.gridLayout_1.addWidget(self.Vgas,0,2,1,1)
        self.gridLayout_1.addWidget(QtWidgets.QLabel(self.tr("Rendimiento:", None)),1,1,1,1)
        self.rendimientoCalculado=Entrada_con_unidades(float, readOnly=True)
        self.gridLayout_1.addWidget(self.rendimientoCalculado,1,2,1,1)
        self.gridLayout_1.addWidget(QtWidgets.QLabel(self.tr("Superficie:", None)),2,1,1,1)
        self.superficie=Entrada_con_unidades(unidades.Area, readOnly=True)
        self.gridLayout_1.addWidget(self.superficie,2,2,1,1)
        gridLayout_Calculo.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Policy.Expanding,QtWidgets.QSizePolicy.Policy.Expanding),12,1,1,6)

        #Salidas
        self.SalidaGas= UI_corriente.Ui_corriente(readOnly=True)
        self.SalidaSolido= UI_corriente.Ui_corriente(readOnly=True)
        self.Salida.addTab(self.SalidaGas,self.tr("Gas filtrado", None))
        self.Salida.addTab(self.SalidaSolido,self.tr("Sólidos recogidos", None))

        self.tabWidget.setCurrentIndex(0)


    def cambiar_entrada(self, corriente):
        self.entrada=corriente
        self.rellenarTablaRendimientos()
        self.calculo()

    def rellenarTablaRendimientos(self):
        self.Rendimientos.clearContents()
        self.Rendimientos.setRowCount(len(self.entrada.solido.distribucion))
        self.Rendimientos.setHorizontalHeaderLabels([self.tr("Diámetro, µm", None), self.tr("Rendimiento", None)])
        for i in range(len(self.entrada.solido.distribucion)):
            self.Rendimientos.setRowHeight(i, 22)
            self.Rendimientos.setItem(i, 0, QtWidgets.QTableWidgetItem(representacion(1e6*self.entrada.solido.diametros[i])))
            self.Rendimientos.item(i, 0).setTextAlignment(QtCore.Qt.AlignmentFlag.AlignRight|QtCore.Qt.AlignmentFlag.AlignVCenter)
            self.Rendimientos.item(i, 0).setFlags(QtCore.Qt.ItemFlag.ItemIsSelectable|QtCore.Qt.ItemFlag.ItemIsEnabled)
            self.Rendimientos.setItem(i, 1, QtWidgets.QTableWidgetItem(""))
            self.Rendimientos.item(i, 1).setTextAlignment(QtCore.Qt.AlignmentFlag.AlignRight|QtCore.Qt.AlignmentFlag.AlignVCenter)

    def cambiarRendimientos(self, fila, columna):
        numero=float(self.Rendimientos.item(fila, columna).text())
        if numero<0 or numero >1:
            self.Rendimientos.item(fila, columna).setText("")
        else:
            if self.rendimientos==[]:
                self.rendimientos=[0]*self.Rendimientos.rowCount()
            self.rendimientos[fila]=numero

    def todos_datos(self):
        if self.TipoCalculo.currentIndex()==0:
            todos_datos=self.numFiltros.value and self.tiempo.value
        elif self.TipoCalculo.currentIndex()==1:
            todos_datos=self.numFiltros.value and self.deltaP.value
        else:
            todos_datos=self.tiempo.value and self.deltaP.value
        return todos_datos and self.Entrada.todos_datos()

    def calculo(self):
        if self.todos_datos():
            if self.Limpieza.value==self.numFiltros.value:
                self.status.setState(5, self.tr("Todos los filtros en limpieza", None))
            else:
                self.status.setState(4)
                self.Equipment(entrada=self.entrada, metodo=self.TipoCalculo.currentIndex(), num_filtros=self.numFiltros.value, tiempo=self.tiempo.value, deltaP=self.deltaP.value.atm, resistenciaFiltro=self.resistenciaFiltro.value, resistenciaTorta=self.resistenciaTorta.value, limpieza=self.Limpieza.value, membranasFiltro=self.MembranaCelda.value, diametroMembrana=self.Diametro.value, areaMembrana=self.Area.value, rendimientos=self.rendimientos)
                self.rellenoSalida()
                if self.rendimientos==[]:
                    self.status.setState(3, self.tr("Usando rendimiento por defecto", None))
                else:
                    self.status.setState(1)

    def rellenoSalida(self):
        if self.TipoCalculo.currentIndex()==0:
            self.deltaP.setValue(self.Equipment.deltaP)
        elif self.TipoCalculo.currentIndex()==1:
            self.tiempo.setValue(self.Equipment.tiempo)
        else:
            self.numFiltros.setValue(self.Equipment.num_filtros)

        self.Vgas.setValue(self.Equipment.Vgas)
        self.rendimientoCalculado.setValue(self.Equipment.rendimiento)
        self.superficie.setValue(self.Equipment.floorArea)
        self.SalidaGas.rellenar(self.Equipment.SalidaAire)
        self.SalidaSolido.rellenar(self.Equipment.SalidaSolido)


    def tipoCalculoCambiado(self, tipo_calculo):
        if tipo_calculo==0:
            self.numFiltros.setReadOnly(False)
            self.numFiltros.setRetornar(True)
            self.numFiltros.setResaltado(True)
            self.tiempo.setReadOnly(False)
            self.tiempo.setRetornar(True)
            self.tiempo.setResaltado(True)
            self.deltaP.setReadOnly(True)
            self.deltaP.setRetornar(False)
            self.deltaP.setResaltado(False)
        elif tipo_calculo==1:
            self.numFiltros.setReadOnly(False)
            self.numFiltros.setRetornar(True)
            self.numFiltros.setResaltado(True)
            self.tiempo.setReadOnly(True)
            self.tiempo.setRetornar(False)
            self.tiempo.setResaltado(False)
            self.deltaP.setReadOnly(False)
            self.deltaP.setRetornar(True)
            self.deltaP.setResaltado(True)
        else:
            self.numFiltros.setReadOnly(True)
            self.numFiltros.setRetornar(False)
            self.numFiltros.setResaltado(False)
            self.tiempo.setReadOnly(False)
            self.tiempo.setRetornar(True)
            self.tiempo.setResaltado(True)
            self.deltaP.setReadOnly(False)
            self.deltaP.setRetornar(True)
            self.deltaP.setResaltado(True)

if __name__ == "__main__":
    import sys
    from lib.corriente import Mezcla, Corriente, Solid
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

    # solido=Solid([64], [138718], distribucion)
    # entrada=Corriente(423.15, 3, 11784,  Mezcla([475], [1]), solido)
    dialogo = UI_equipment()
    dialogo.show()
    sys.exit(app.exec())
