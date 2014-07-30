#!/usr/bin/python
# -*- coding: utf-8 -*-

############################################################################
###                                       Diálogo de definición de unidades de secado de sólidos, UI_dryer                                      ###
############################################################################

from functools import partial

from PyQt4 import QtCore, QtGui

from equipment.gas_solid_liquid import Dryer
from lib import unidades
from UI import UI_corriente
from equipment.parents import UI_equip
from UI.widgets import Entrada_con_unidades
from tools import costIndex


class UI_equipment(UI_equip):
    """Dialogo de definición de unidades de secado de sólidos"""
    Equipment=Dryer()
    def __init__(self, equipment=None, parent=None):
        """entrada: Parametro opcional de clase corriente que indica la corriente de entrada"""
        super(UI_equipment, self).__init__(Dryer, costos=False, parent=parent)

        #Pestaña entrada
        self.EntradaSolido= UI_corriente.Ui_corriente()
        self.EntradaSolido.Changed.connect(partial(self.changeParams, "entradaSolido"))
        self.Entrada.addTab(self.EntradaSolido,QtGui.QApplication.translate("pychemqt", "Humid Solid"))
        self.EntradaAire= UI_corriente.Ui_psychrometry()
        self.EntradaAire.Changed.connect(partial(self.changeParams, "entradaAire"))
        self.Entrada.addTab(self.EntradaAire,QtGui.QApplication.translate("pychemqt", "Air"))
        
        #Pestaña calculo
        gridLayout_Calculo = QtGui.QGridLayout(self.tabCalculo)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Mode")),1,1)
        self.TipoCalculo=QtGui.QComboBox()
        self.TipoCalculo.addItem(QtGui.QApplication.translate("pychemqt", "Rating, calculate output conditions"))
        self.TipoCalculo.addItem(QtGui.QApplication.translate("pychemqt", "Design, calculate air flow necessary"))
        self.TipoCalculo.currentIndexChanged.connect(partial(self.changeParams, "tipoCalculo"))
        gridLayout_Calculo.addWidget(self.TipoCalculo,1,2,1,4)
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),2,1,1,6)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Air Relative Humidity")),3,1)
        self.HumedadAire=Entrada_con_unidades(float, max=1, spinbox=True, step=0.01)
        self.HumedadAire.valueChanged.connect(partial(self.changeParams, "HR"))
        gridLayout_Calculo.addWidget(self.HumedadAire,3,2)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Product moisture fraction")),4,1)
        self.HumedadSolido=Entrada_con_unidades(float, max=1., spinbox=True, step=0.01, textounidad=unidades.Mass(None).text()+"/"+unidades.Mass(None).text())
        self.HumedadSolido.valueChanged.connect(partial(self.changeParams, "HumedadResidual"))
        gridLayout_Calculo.addWidget(self.HumedadSolido,4,2)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Output Solid Temperature")),5,1)
        self.temperatura=Entrada_con_unidades(unidades.Temperature)
        self.temperatura.valueChanged.connect(partial(self.changeParams, "TemperaturaSolid"))
        gridLayout_Calculo.addWidget(self.temperatura,5,2)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Heat Duty")),6,1)
        self.Heat=Entrada_con_unidades(unidades.Power)
        self.Heat.valueChanged.connect(partial(self.changeParams, "Heat"))
        gridLayout_Calculo.addWidget(self.Heat,6,2) 
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Pressure Drop")),7,1)
        self.DeltaP=Entrada_con_unidades(unidades.Pressure)
        self.DeltaP.valueChanged.connect(partial(self.changeParams, "DeltaP"))
        gridLayout_Calculo.addWidget(self.DeltaP,7,2) 
  
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),8,1,1,6)
        self.groupBox_Calculo = QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Results"))
        gridLayout_Calculo.addWidget(self.groupBox_Calculo,9,1,1,5)
        gridLayout_1 = QtGui.QGridLayout(self.groupBox_Calculo)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Output Temperature")),1,1)
        self.temperaturaCalculada=Entrada_con_unidades(unidades.Temperature, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.temperaturaCalculada,1,2)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Air Flow")),2,1)
        self.caudalVolumetrico=Entrada_con_unidades(unidades.VolFlow, "QGas", retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.caudalVolumetrico,2,2)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Output Air Relative Humidity")),3,1)
        self.HumedadCalculada=Entrada_con_unidades(float, readOnly=True, textounidad="%")
        gridLayout_1.addWidget(self.HumedadCalculada,3,2)
        
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),11,1,1,6)

        #Pestaña salida
        self.SalidaSolido= UI_corriente.Ui_corriente(readOnly=True)
        self.Salida.addTab(self.SalidaSolido,QtGui.QApplication.translate("pychemqt", "Dry solid"))
        self.SalidaAire= UI_corriente.Ui_psychrometry(readOnly=True)
        self.Salida.addTab(self.SalidaAire,QtGui.QApplication.translate("pychemqt", "Air"))
        
        if equipment:
            self.setEquipment(equipment)

    def rellenar(self):
        self.EntradaAire.setCorriente(self.Equipment.kwargs["entradaAire"])
        self.EntradaSolido.setCorriente(self.Equipment.kwargs["entradaSolido"])
        if self.Equipment.status:
            self.temperaturaCalculada.setValue(self.Equipment.SalidaSolido.T)
            self.caudalVolumetrico.setValue(self.Equipment.entradaAire.corriente.Q)
            self.HumedadCalculada.setValue(self.Equipment.SalidaAire.Xw*100)
            self.SalidaAire.setCorriente(self.Equipment.SalidaAire)
            self.SalidaSolido.setCorriente(self.Equipment.SalidaSolido)
            if self.Equipment.kwargs["tipoCalculo"]==1:
                self.EntradaAire.setCorriente(self.Equipment.entradaAire)


if __name__ == "__main__":
    import sys        
    from lib.corriente import Mezcla, Corriente, Solid
    from lib.psycrometry import Punto_Psicrometrico
    app = QtGui.QApplication(sys.argv)
    diametros=[96.5, 105, 110, 118, 125, 130, 140, 150, 170]
    fraccion=[0.02, 0.05, 0.1, 0.15, 0.25, 0.2, 0.15, 0.05, 0.03]
    solido=Solid(caudalSolido=[5000], distribucion_fraccion=fraccion, distribucion_diametro=diametros)
    Solido=Corriente(T=300, P=101325., caudalMasico=50, fraccionMolar=[1, 0], solido=solido)
    Aire=Punto_Psicrometrico(caudal=100, tdb=300, HR=50)
    secador=Dryer(entradaSolido=Solido, entradaAire=Aire)
    dialogo = UI_equipment(secador)
    dialogo.show()
    sys.exit(app.exec_())
