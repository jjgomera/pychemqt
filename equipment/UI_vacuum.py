#!/usr/bin/python
# -*- coding: utf-8 -*-

##########################################################################
###                                              Diálogo de definición de filtros de mangas, UI_baghouse                                       ###
##########################################################################

from PyQt4 import QtCore, QtGui

from gas_solid import Baghouse
from lib import unidades, config
from UI import UI_corriente
from equipment import parents
from UI.delegate import CellEditor
from UI.widgets import Entrada_con_unidades


class UI_equipment(parents.UI_equipment):
    """Diálogo de definición de filtros de mangas"""
    def __init__(self, entrada=None, parent=None):
        """entrada: Parametro opcional de clase corriente que indica la corriente de entrada"""
        super(UI_equipment, self).__init__(Baghouse, entrada=False, costos=False, parent=parent)
        self.entrada=entrada

        #Pestaña entrada
        self.Entrada= UI_corriente.Ui_corriente(entrada)
        self.Entrada.Changed.connect(self.cambiar_entrada)
        self.tabWidget.insertTab(0, self.Entrada,QtGui.QApplication.translate("equipment", "Entrada", None, QtGui.QApplication.UnicodeUTF8))

        #Pestaña definición rendimientos
        self.Rendimientos= QtGui.QTableWidget(1, 2)
        self.Rendimientos.setItemDelegateForColumn(1, CellEditor(self))
        self.Rendimientos.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
        self.Rendimientos.verticalHeader().hide()
        self.Rendimientos.setEditTriggers(QtGui.QAbstractItemView.AllEditTriggers)
        if self.entrada:
            self.rellenarTablaRendimientos()
        self.rendimientos=[]
        self.Rendimientos.cellChanged.connect(self.cambiarRendimientos)
        self.tabWidget.insertTab(1, self.Rendimientos,QtGui.QApplication.translate("equipment", "Rendimientos", None, QtGui.QApplication.UnicodeUTF8))
        
        #Cálculo
        gridLayout_Calculo = QtGui.QGridLayout(self.tabCalculo)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Tipo de cálculo:", None, QtGui.QApplication.UnicodeUTF8)), 1, 1, 1, 1)
        self.TipoCalculo=QtGui.QComboBox()
        self.TipoCalculo.addItem(QtGui.QApplication.translate("equipment", "Calcular caída de presión", None, QtGui.QApplication.UnicodeUTF8))
        self.TipoCalculo.addItem(QtGui.QApplication.translate("equipment", "Calcular tiempo de filtración", None, QtGui.QApplication.UnicodeUTF8))
        self.TipoCalculo.addItem(QtGui.QApplication.translate("equipment", "Calcular número de filtros", None, QtGui.QApplication.UnicodeUTF8))
        self.TipoCalculo.currentIndexChanged.connect(self.tipoCalculoCambiado)
        gridLayout_Calculo.addWidget(self.TipoCalculo, 1, 2, 1, 4)
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),2,1,1,6)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Nº de filtros:", None, QtGui.QApplication.UnicodeUTF8)), 3, 1, 1, 1)
        self.numFiltros=Entrada_con_unidades(int, spinbox=True, step=1, width=50, resaltado=True, min=1, start=1)
        self.numFiltros.valueChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.numFiltros,3,2,1,1)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Tiempo de filtración:", None, QtGui.QApplication.UnicodeUTF8)), 4, 1, 1, 1)
        self.tiempo=Entrada_con_unidades(unidades.Time, resaltado=True)
        self.tiempo.valueChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.tiempo,4,2,1,1)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Pérdida de presión:", None, QtGui.QApplication.UnicodeUTF8)), 5, 1, 1, 1)
        self.deltaP=Entrada_con_unidades(unidades.Pressure, retornar=False, readOnly=True)
        self.deltaP.valueChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.deltaP,5,2,1,1)
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),6,1,1,6)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Membranas por filtro:", None, QtGui.QApplication.UnicodeUTF8)), 7, 1, 1, 1)
        self.MembranaCelda=Entrada_con_unidades(int, spinbox=True, step=1, width=70, value=78, min=1)
        self.MembranaCelda.valueChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.MembranaCelda,7,2,1,1) 
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Diametro de membrana:", None, QtGui.QApplication.UnicodeUTF8)), 8, 1, 1, 1)
        self.Diametro=Entrada_con_unidades(unidades.Length, value=unidades.Length(0.5, "ft"))
        self.Diametro.valueChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.Diametro,8,2,1,1) 
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Area por membrana:", None, QtGui.QApplication.UnicodeUTF8)), 9, 1, 1, 1)
        self.Area=Entrada_con_unidades(unidades.Area, value=unidades.Area(16, "ft2"))
        self.Area.valueChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.Area,9,2,1,1) 
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Resistencia filtro:", None, QtGui.QApplication.UnicodeUTF8)), 7, 4, 1, 1)
        self.resistenciaFiltro=Entrada_con_unidades(float, spinbox=True, step=0.01, width=70, value=0.84, min=0)
        self.resistenciaFiltro.valueChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.resistenciaFiltro,7,5,1,1) 
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Resistencia torta:", None, QtGui.QApplication.UnicodeUTF8)), 8, 4, 1, 1)
        self.resistenciaTorta=Entrada_con_unidades(float, spinbox=True, step=0.01, width=70, value=0.1, min=0)
        self.resistenciaTorta.valueChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.resistenciaTorta,8,5,1,1) 
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Filtros en limpieza:", None, QtGui.QApplication.UnicodeUTF8)), 9, 4, 1, 1)
        self.Limpieza=Entrada_con_unidades(int, spinbox=True, step=1, width=70, value=1, min=0)
        self.Limpieza.valueChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.Limpieza,9,5,1,1)  
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),10,1,1,6)
        
        self.groupBox_Calculo = QtGui.QGroupBox(QtGui.QApplication.translate("equipment", "Datos calculados", None, QtGui.QApplication.UnicodeUTF8))
        gridLayout_Calculo.addWidget(self.groupBox_Calculo,11,1,1,5)
        self.gridLayout_1 = QtGui.QGridLayout(self.groupBox_Calculo)
        self.gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "V<sub>gas</sub>:", None, QtGui.QApplication.UnicodeUTF8)),0,1,1,1)
        self.Vgas=Entrada_con_unidades(unidades.Speed, retornar=False, readOnly=True)
        self.gridLayout_1.addWidget(self.Vgas,0,2,1,1)
        self.gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Rendimiento:", None, QtGui.QApplication.UnicodeUTF8)),1,1,1,1)
        self.rendimientoCalculado=Entrada_con_unidades(float, readOnly=True)
        self.gridLayout_1.addWidget(self.rendimientoCalculado,1,2,1,1)
        self.gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Superficie:", None, QtGui.QApplication.UnicodeUTF8)),2,1,1,1)
        self.superficie=Entrada_con_unidades(unidades.Area, readOnly=True)
        self.gridLayout_1.addWidget(self.superficie,2,2,1,1)        
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),12,1,1,6)

        #Salidas
        self.SalidaGas= UI_corriente.Ui_corriente(readOnly=True)
        self.SalidaSolido= UI_corriente.Ui_corriente(readOnly=True)
        self.Salida.addTab(self.SalidaGas,QtGui.QApplication.translate("equipment", "Gas filtrado", None, QtGui.QApplication.UnicodeUTF8))
        self.Salida.addTab(self.SalidaSolido,QtGui.QApplication.translate("equipment", "Sólidos recogidos", None, QtGui.QApplication.UnicodeUTF8))
        
        self.tabWidget.setCurrentIndex(0)


    def cambiar_entrada(self, corriente):
        self.entrada=corriente
        self.rellenarTablaRendimientos()
        self.calculo()
        
    def rellenarTablaRendimientos(self):
        self.Rendimientos.clearContents()
        self.Rendimientos.setRowCount(len(self.entrada.solido.distribucion))
        self.Rendimientos.setHorizontalHeaderLabels([QtGui.QApplication.translate("equipment", "Diámetro, µm", None, QtGui.QApplication.UnicodeUTF8), QtGui.QApplication.translate("equipment", "Rendimiento", None, QtGui.QApplication.UnicodeUTF8)])
        for i in range(len(self.entrada.solido.distribucion)):
            self.Rendimientos.setRowHeight(i, 22)
            self.Rendimientos.setItem(i, 0, QtGui.QTableWidgetItem(config.representacion(1e6*self.entrada.solido.diametros[i])))
            self.Rendimientos.item(i, 0).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
            self.Rendimientos.item(i, 0).setFlags(QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsEnabled)
            self.Rendimientos.setItem(i, 1, QtGui.QTableWidgetItem(""))
            self.Rendimientos.item(i, 1).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)

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
                self.status.setState(5, QtGui.QApplication.translate("equipment", "Todos los filtros en limpieza", None, QtGui.QApplication.UnicodeUTF8))
            else:
                self.status.setState(4)
                self.Equipment(entrada=self.entrada, metodo=self.TipoCalculo.currentIndex(), num_filtros=self.numFiltros.value, tiempo=self.tiempo.value, deltaP=self.deltaP.value.atm, resistenciaFiltro=self.resistenciaFiltro.value, resistenciaTorta=self.resistenciaTorta.value, limpieza=self.Limpieza.value, membranasFiltro=self.MembranaCelda.value, diametroMembrana=self.Diametro.value, areaMembrana=self.Area.value, rendimientos=self.rendimientos)
                self.rellenoSalida()
                if self.rendimientos==[]:
                    self.status.setState(3, QtGui.QApplication.translate("equipment", "Usando rendimiento por defecto", None, QtGui.QApplication.UnicodeUTF8))
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
    app = QtGui.QApplication(sys.argv)
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
    entrada=Corriente(423.15, 3, 11784,  Mezcla([475], [1]), solido)
    dialogo = UI_equipment(entrada)
    dialogo.show()
    sys.exit(app.exec_())
