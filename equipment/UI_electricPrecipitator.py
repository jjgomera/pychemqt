#!/usr/bin/python
# -*- coding: utf-8 -*-

#######################################################
###    Diálogo de definición de precipitadores electrostáticos, UI_electricPrecipitator   ###
#######################################################

from functools import partial

from PyQt4 import QtCore, QtGui

from lib.unidades import DeltaP, PotencialElectric, Area
from equipment.gas_solid import ElectricPrecipitator, UI_equipment_Solid
from UI import UI_corriente
from UI.widgets import Entrada_con_unidades


class UI_equipment(UI_equipment_Solid):
    """Diálogo de definición de precipitadores electrostáticos"""
    Equipment=ElectricPrecipitator()
    def __init__(self, equipment=None, parent=None):
        """
        equipment: instancia de equipo inicial
        """
        super(UI_equipment, self).__init__(ElectricPrecipitator, entrada=False, parent=parent)

        #Pestaña calculo
        gridLayout_Calculo = QtGui.QGridLayout(self.tabCalculo)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Mode")),1,1)
        self.metodo=QtGui.QComboBox()
        for txt in self.Equipment.TEXT_TIPO:
            self.metodo.addItem(txt)
        self.metodo.currentIndexChanged.connect(self.tipoCalculoCambiado)
        gridLayout_Calculo.addWidget(self.metodo, 1, 2, 1, 4)
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),2,1,1,6)

        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Area")),3,1)
        self.area=Entrada_con_unidades(Area, resaltado=True)
        self.area.valueChanged.connect(partial(self.changeParams, "area"))
        gridLayout_Calculo.addWidget(self.area,3,2)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Allowable efficiency")),4,1)
        self.rendimientoAdmisible=Entrada_con_unidades(float,  readOnly=True)
        self.rendimientoAdmisible.valueChanged.connect(partial(self.changeParams, "rendimientoAdmisible"))
        gridLayout_Calculo.addWidget(self.rendimientoAdmisible,4,2)
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),5,1,1,6)

        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Dielectric constant")),6,1)
        self.epsilon=Entrada_con_unidades(float)
        self.epsilon.valueChanged.connect(partial(self.changeParams, "epsilon"))
        gridLayout_Calculo.addWidget(self.epsilon,6,2)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Charging field")),7,1)
        self.potencialCarga=Entrada_con_unidades(PotencialElectric)
        self.potencialCarga.valueChanged.connect(partial(self.changeParams, "potencialCarga"))
        gridLayout_Calculo.addWidget(self.potencialCarga,7,2)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Collecting field")),8,1)
        self.potencialDescarga=Entrada_con_unidades(PotencialElectric)
        self.potencialDescarga.valueChanged.connect(partial(self.changeParams, "potencialDescarga"))
        gridLayout_Calculo.addWidget(self.potencialDescarga,8,2)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Pressure drop")),9,1)
        self.deltaP=Entrada_con_unidades(DeltaP)
        self.deltaP.valueChanged.connect(partial(self.changeParams, "deltaP"))
        gridLayout_Calculo.addWidget(self.deltaP,9,2)
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),10,1,1,6)

        self.groupBox_Calculo = QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Result"))
        gridLayout_Calculo.addWidget(self.groupBox_Calculo,11,1,1,5)
        gridLayout_1 = QtGui.QGridLayout(self.groupBox_Calculo)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Area")),0,1)
        self.areaCalculada=Entrada_con_unidades(Area, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.areaCalculada,0,2)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Efficiency")),1,1)
        self.rendimiento=Entrada_con_unidades(float, readOnly=True)
        gridLayout_1.addWidget(self.rendimiento,1,2)
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),12,1,1,6)

        #Salidas
        self.SalidaGas= UI_corriente.Ui_corriente(readOnly=True)
        self.Salida.addTab(self.SalidaGas,QtGui.QApplication.translate("pychemqt", "Filtered gas"))
        self.SalidaSolido= UI_corriente.Ui_corriente(readOnly=True)
        self.Salida.addTab(self.SalidaSolido,QtGui.QApplication.translate("pychemqt", "Collected solids"))

        if equipment:
            self.setEquipment(equipment)


    def tipoCalculoCambiado(self, tipo_calculo):
        self.area.setReadOnly(tipo_calculo)
        self.area.setResaltado(not tipo_calculo)
        self.rendimientoAdmisible.setReadOnly(not tipo_calculo)
        self.rendimientoAdmisible.setResaltado(tipo_calculo)
        self.changeParams("metodo", tipo_calculo)


if __name__ == "__main__":
    import sys
    from lib.corriente import Corriente, Solid
    app = QtGui.QApplication(sys.argv)
    diametros=[17.5e-6, 22.4e-6, 26.2e-6, 31.8e-6, 37e-6, 42.4e-6, 48e-6, 54e-6, 60e-6, 69e-6, 81.3e-6, 96.5e-6, 109e-6, 127e-6]
    fracciones=[0.02, 0.03, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.03, 0.02]
    solido=Solid(caudalSolido=[0.1], distribucion_diametro=diametros, distribucion_fraccion=fracciones)
    corriente=Corriente(T=300, P=101325, caudalMasico=1.,  fraccionMolar=[1.], solido=solido)
    precipitador=ElectricPrecipitator(entrada=corriente, metodo=1, rendimientoAdmisible=0.9)
    dialogo = UI_equipment(precipitador)
    dialogo.show()
    sys.exit(app.exec_())
