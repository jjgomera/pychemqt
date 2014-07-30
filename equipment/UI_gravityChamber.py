#!/usr/bin/python
# -*- coding: utf-8 -*-

############################################################
###  Diálogo de definición de camaras de sedimentación por gravedad, UI_gravityChambers   ###
############################################################

from functools import partial

from PyQt4 import QtCore, QtGui

from lib.unidades import Length, Speed, VolFlow, DeltaP
from gas_solid import GravityChamber, UI_equipment_Solid
from UI import UI_corriente
from UI.widgets import Entrada_con_unidades


class UI_equipment(UI_equipment_Solid):
    """Diálogo de definición de cámaras de sedimentación por gravedad"""
    Equipment=GravityChamber()
    def __init__(self, equipment=None, parent=None):
        """
        equipment: instancia de equipo inicial
        """
        super(UI_equipment, self).__init__(GravityChamber, costos=False, entrada=False, parent=parent)

        #Cálculo
        gridLayout_Calculo = QtGui.QGridLayout(self.tabCalculo)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Mode")), 1, 1, 1, 1)
        self.metodo=QtGui.QComboBox()
        for txt in self.Equipment.TEXT_TIPO:
            self.metodo.addItem(txt)
        self.metodo.currentIndexChanged.connect(self.tipoCalculoCambiado)
        gridLayout_Calculo.addWidget(self.metodo, 1, 2, 1, 4)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Model")), 2, 1, 1, 1)
        self.modelo=QtGui.QComboBox()
        for txt in self.Equipment.TEXT_MODEL:
            self.modelo.addItem(txt)
        self.modelo.currentIndexChanged.connect(partial(self.changeParams, "modelo"))
        gridLayout_Calculo.addWidget(self.modelo, 2, 2, 1, 1)
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),3,1,1,6)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Width")), 4, 1, 1, 1)
        self.W=Entrada_con_unidades(Length)
        self.W.valueChanged.connect(partial(self.changeParams, "W"))
        gridLayout_Calculo.addWidget(self.W,4,2,1,1)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Height")), 5, 1, 1, 1)
        self.H=Entrada_con_unidades(Length)
        self.H.valueChanged.connect(partial(self.changeParams, "H"))
        gridLayout_Calculo.addWidget(self.H,5,2,1,1)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Length")), 6, 1, 1, 1)
        self.L=Entrada_con_unidades(Length)
        self.L.valueChanged.connect(partial(self.changeParams, "L"))
        gridLayout_Calculo.addWidget(self.L,6,2,1,1)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Allowable efficiency")), 7, 1, 1, 1)
        self.rendimientoAdmisible=Entrada_con_unidades(float, spinbox=True, max=1)
        self.rendimientoAdmisible.valueChanged.connect(partial(self.changeParams, "rendimientoAdmisible"))
        gridLayout_Calculo.addWidget(self.rendimientoAdmisible,7,2,1,1) 
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Allowable speed")), 8, 1, 1, 1)
        self.velocidadAdmisible=Entrada_con_unidades(Speed)
        self.velocidadAdmisible.valueChanged.connect(partial(self.changeParams, "velocidadAdmisible"))
        gridLayout_Calculo.addWidget(self.velocidadAdmisible,8,2,1,1) 
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Pressure loss")), 9, 1, 1, 1)
        self.deltaP=Entrada_con_unidades(DeltaP)
        self.deltaP.valueChanged.connect(partial(self.changeParams, "deltaP"))
        gridLayout_Calculo.addWidget(self.deltaP,9,2,1,1) 
  
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),10,1,1,6)
        self.groupBox_Calculo = QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Results"))
        gridLayout_Calculo.addWidget(self.groupBox_Calculo,11,1,1,5)
        gridLayout_1 = QtGui.QGridLayout(self.groupBox_Calculo)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Flow")),0,1,1,1)
        self.Q=Entrada_con_unidades(VolFlow, "QGas", retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.Q,0,2,1,1)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "V<sub>gas</sub>")),1,1,1,1)
        self.Vgas=Entrada_con_unidades(Speed, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.Vgas,1,2,1,1)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Efficiency")),2,1,1,1)
        self.rendimiento=Entrada_con_unidades(float, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.rendimiento,2,2,1,1)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Height")),0,4,1,1)
        self.HCalc=Entrada_con_unidades(Length, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.HCalc,0,5,1,1)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Width")),1,4,1,1)
        self.WCalc=Entrada_con_unidades(Length, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.WCalc,1,5,1,1)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Length")),2,4,1,1)
        self.LCalc=Entrada_con_unidades(Length, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.LCalc,2,5,1,1)
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),12,1,1,6)

        #Salidas
        self.SalidaGas= UI_corriente.Ui_corriente(readOnly=True)
        self.Salida.addTab(self.SalidaGas, QtGui.QApplication.translate("pychemqt", "Filtered gas"))
        self.SalidaSolido= UI_corriente.Ui_corriente(readOnly=True)
        self.Salida.addTab(self.SalidaSolido, QtGui.QApplication.translate("pychemqt", "Collected solids"))
        
        self.tipoCalculoCambiado(0)
        
        if equipment:
            self.setEquipment(equipment)


    def tipoCalculoCambiado(self, int):
        self.W.setReadOnly(int)
        self.W.setRetornar(not int)
        self.W.setResaltado(not int)
        self.H.setResaltado(not int)
        self.L.setReadOnly(int)
        self.L.setRetornar(not int)
        self.L.setResaltado(not int)
        self.rendimientoAdmisible.setReadOnly(not int)
        self.rendimientoAdmisible.setResaltado(int)
        self.velocidadAdmisible.setReadOnly(not int)
        self.velocidadAdmisible.setResaltado(False)
        self.changeParams("metodo", int)


if __name__ == "__main__":
    import sys        
    from lib.corriente import Corriente, Solid
    app = QtGui.QApplication(sys.argv)
    diametros=[17.5e-6, 22.4e-6, 26.2e-6, 31.8e-6, 37e-6, 42.4e-6, 48e-6, 54e-6, 60e-6, 69e-6, 81.3e-6, 96.5e-6, 109e-6, 127e-6]
    fracciones=[0.02, 0.03, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.03, 0.02]
    solido=Solid(caudalSolido=[0.1], distribucion_diametro=diametros, distribucion_fraccion=fracciones)
    corriente=Corriente(T=300, P=101325, caudalMasico=1.,  fraccionMolar=[1.], solido=solido)
    camara=GravityChamber(entrada=corriente, metodo=1, modelo=1, H=2., rendimientoAdmisible=0.95)
    dialogo = UI_equipment(camara)
    dialogo.show()
    sys.exit(app.exec_())
