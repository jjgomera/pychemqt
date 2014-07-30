#!/usr/bin/python
# -*- coding: utf-8 -*-

#########################################################################
# Diálogo de definición de filtros de mangas, UI_baghouse
#########################################################################

from functools import partial

from PyQt4 import QtGui
from numpy import any

from lib.unidades import Time, Pressure, Length, Area, Speed
from equipment.gas_solid import Baghouse, UI_equipment_Solid
from UI import UI_corriente
from UI.widgets import Entrada_con_unidades, Tabla


class UI_equipment(UI_equipment_Solid):
    """Diálogo de definición de filtros de mangas"""
    Equipment=Baghouse()
    def __init__(self, equipment=None, parent=None):
        """
        equipment: instancia de equipo inicial
        """
        super(UI_equipment, self).__init__(Baghouse, entrada=False, costos=False, parent=parent)

        #Pestaña definición rendimientos
        self.Rendimientos=Tabla(2, horizontalHeader=[QtGui.QApplication.translate("pychemqt", "Diameter")+", "+Length.text("ParticleDiameter"), QtGui.QApplication.translate("pychemqt", "Efficiency")], filas=1, stretch=False)
        self.Rendimientos.setColumnReadOnly(0, True)
        self.Rendimientos.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
        self.Rendimientos.editingFinished.connect(self.cambiarRendimientos)
        self.tabWidget.insertTab(1, self.Rendimientos, QtGui.QApplication.translate("pychemqt", "Efficiencies"))

        #Cálculo
        gridLayout_Calculo = QtGui.QGridLayout(self.tabCalculo)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Mode")),1,1)
        self.metodo=QtGui.QComboBox()
        for txt in self.Equipment.TEXT_TIPO:
            self.metodo.addItem(txt)
        self.metodo.currentIndexChanged.connect(self.tipoCalculoCambiado)
        gridLayout_Calculo.addWidget(self.metodo,1,2,1,3)
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),2,1,1,6)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "No cells")),3,1)
        self.num_filtros=Entrada_con_unidades(int, spinbox=True, step=1, width=50, resaltado=True, min=1, start=1)
        self.num_filtros.valueChanged.connect(partial(self.changeParams, "num_filtros"))
        gridLayout_Calculo.addWidget(self.num_filtros,3,2)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Time")),4,1)
        self.tiempo=Entrada_con_unidades(Time, resaltado=True)
        self.tiempo.valueChanged.connect(partial(self.changeParams, "tiempo"))
        gridLayout_Calculo.addWidget(self.tiempo,4,2)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Pressure drop")),5,1)
        self.deltaP=Entrada_con_unidades(Pressure, retornar=False, readOnly=True)
        self.deltaP.valueChanged.connect(partial(self.changeParams, "deltaP"))
        gridLayout_Calculo.addWidget(self.deltaP,5,2)
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),6,1,1,6)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Bags per cell")),7,1)
        self.membranasFiltro=Entrada_con_unidades(int, spinbox=True, step=1, min=1)
        self.membranasFiltro.valueChanged.connect(partial(self.changeParams, "membranasFiltro"))
        gridLayout_Calculo.addWidget(self.membranasFiltro,7,2)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Bag diameter")),8,1)
        self.diametroMembrana=Entrada_con_unidades(Length)
        self.diametroMembrana.valueChanged.connect(partial(self.changeParams, "diametroMembrana"))
        gridLayout_Calculo.addWidget(self.diametroMembrana,8,2)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Area per bag")),9,1)
        self.areaMembrana=Entrada_con_unidades(Area)
        self.areaMembrana.valueChanged.connect(partial(self.changeParams, "areaMembrana"))
        gridLayout_Calculo.addWidget(self.areaMembrana,9,2)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Cloth resistence")),7,4)
        self.resistenciaFiltro=Entrada_con_unidades(float)
        self.resistenciaFiltro.valueChanged.connect(partial(self.changeParams, "resistenciaFiltro"))
        gridLayout_Calculo.addWidget(self.resistenciaFiltro,7,5)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Cake resistence")),8,4)
        self.resistenciaTorta=Entrada_con_unidades(float)
        self.resistenciaTorta.valueChanged.connect(partial(self.changeParams, "resistenciaTorta"))
        gridLayout_Calculo.addWidget(self.resistenciaTorta,8,5)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Cells cleaned")),9,4)
        self.limpieza=Entrada_con_unidades(int, spinbox=True, step=1, min=0)
        self.limpieza.valueChanged.connect(partial(self.changeParams, "limpieza"))
        gridLayout_Calculo.addWidget(self.limpieza,9,5)
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),10,1,1,6)

        self.groupBox_Calculo = QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Results"))
        gridLayout_Calculo.addWidget(self.groupBox_Calculo,11,1,1,5)
        gridLayout_1 = QtGui.QGridLayout(self.groupBox_Calculo)

        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "No cells")),1,1)
        self.num_filtrosCalc=Entrada_con_unidades(int, readOnly=True)
        gridLayout_1.addWidget(self.num_filtrosCalc,1,2)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Time")),2,1)
        self.tiempoCalc=Entrada_con_unidades(Time, readOnly=True)
        gridLayout_1.addWidget(self.tiempoCalc,2,2)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Pressure drop")),3,1)
        self.deltaPCalc=Entrada_con_unidades(Pressure, readOnly=True)
        gridLayout_1.addWidget(self.deltaPCalc,3,2)

        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Gas velocity")),1,4)
        self.Vgas=Entrada_con_unidades(Speed, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.Vgas,1,5)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Efficiency")),2,4)
        self.rendimiento=Entrada_con_unidades(float, readOnly=True)
        gridLayout_1.addWidget(self.rendimiento,2,5)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Area")),3,4)
        self.floorArea=Entrada_con_unidades(Area, readOnly=True)
        gridLayout_1.addWidget(self.floorArea,3,5)
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),12,1,1,6)

        #Salidas
        self.SalidaGas= UI_corriente.Ui_corriente(readOnly=True)
        self.SalidaSolido= UI_corriente.Ui_corriente(readOnly=True)
        self.Salida.addTab(self.SalidaGas,QtGui.QApplication.translate("pychemqt", "Filtered gas"))
        self.Salida.addTab(self.SalidaSolido,QtGui.QApplication.translate("pychemqt", "Collected solids"))

        if equipment:
            self.setEquipment(equipment)

    def cambiarRendimientos(self):
        self.changeParams("rendimientos", self.Rendimientos.getColumn(1))


    def tipoCalculoCambiado(self, tipo_calculo):
        if tipo_calculo==0:
            self.num_filtros.setReadOnly(False)
            self.num_filtros.setResaltado(True)
            self.tiempo.setReadOnly(False)
            self.tiempo.setResaltado(True)
            self.deltaP.setReadOnly(True)
            self.deltaP.setResaltado(False)
        elif tipo_calculo==1:
            self.num_filtros.setReadOnly(False)
            self.num_filtros.setResaltado(True)
            self.tiempo.setReadOnly(True)
            self.tiempo.setResaltado(False)
            self.deltaP.setReadOnly(False)
            self.deltaP.setResaltado(True)
        else:
            self.num_filtros.setReadOnly(True)
            self.num_filtros.setResaltado(False)
            self.tiempo.setReadOnly(False)
            self.tiempo.setResaltado(True)
            self.deltaP.setReadOnly(False)
            self.deltaP.setResaltado(True)
        self.changeParams("metodo", tipo_calculo)


    def rellenarInput(self):
        UI_equipment_Solid.rellenarInput(self)
        if self.Equipment.kwargs["entrada"].solido:
            diametros=[d.config("ParticleDiameter") for d in self.Equipment.kwargs["entrada"].solido.diametros]
            self.Rendimientos.setColumn(0, diametros)
        if any(self.Equipment.kwargs["rendimientos"]):
            self.Rendimientos.setColumn(1, self.Equipment.kwargs["rendimientos"])


if __name__ == "__main__":
    import sys
    from lib.corriente import Corriente, Solid
    app = QtGui.QApplication(sys.argv)
    diametros=[17.5e-6, 22.4e-6, 26.2e-6, 31.8e-6, 37e-6, 42.4e-6, 48e-6, 54e-6, 60e-6, 69e-6, 81.3e-6, 96.5e-6, 109e-6, 127e-6]
    fracciones=[0.02, 0.03, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.03, 0.02]
    solido=Solid(caudalSolido=[0.1], distribucion_diametro=diametros, distribucion_fraccion=fracciones)
    corriente=Corriente(T=300, P=101325, caudalMasico=1.,  fraccionMolar=[1.], solido=solido)
    filtro=Baghouse(entrada=corriente, metodo=0, num_filtros=4, tiempo=3600)
    dialogo = UI_equipment(filtro)
    dialogo.show()
    sys.exit(app.exec_())
