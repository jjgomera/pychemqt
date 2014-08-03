#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Diálogo de definición de lavadores de gases, UI_scrubber
###############################################################################


from functools import partial

from PyQt4 import QtGui

from lib.unidades import Length, DeltaP
from UI import UI_corriente
from UI.widgets import Entrada_con_unidades
from equipment.gas_solid_liquid import Scrubber
from equipment.parents import UI_equip


class UI_equipment(UI_equip):
    Equipment = Scrubber()
    """Dialogo de definición de unidades de lavado de gases"""
    def __init__(self, equipment=None, parent=None):
        """
        equipment: instancia de equipo inicial
        """
        super(UI_equipment, self).__init__(Scrubber, costos=False, parent=parent)

        # Pestaña entrada
        self.entradaGas = UI_corriente.Ui_corriente()
        self.entradaGas.Changed.connect(partial(self.changeParams,"entradaGas"))
        self.entrada.addTab(self.entradaGas, QtGui.QApplication.translate("equipment", "Gas"))
        self.entradaLiquido = UI_corriente.Ui_corriente()
        self.entradaLiquido.Changed.connect(partial(self.changeParams, "entradaLiquido"))
        self.entrada.addTab(self.entradaLiquido, QtGui.QApplication.translate("pychemqt", "Liquid"))

        # Pestaña calculo
        gridLayout_Calculo = QtGui.QGridLayout(self.tabCalculo)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Mode")),1,1)
        self.tipo_calculo = QtGui.QComboBox()
        for txt in self.Equipment.TEXT_TIPO:
            self.tipo_calculo.addItem(txt)
        self.tipo_calculo.currentIndexChanged.connect(self.on_tipoCalculo_currentIndexChanged)
        gridLayout_Calculo.addWidget(self.tipo_calculo,1,2,1,5)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Method")),2,1)
        self.modelo_rendimiento = QtGui.QComboBox()
        for txt in self.Equipment.TEXT_MODEL:
            self.modelo_rendimiento.addItem(txt)
        self.modelo_rendimiento.currentIndexChanged.connect(self.on_modeloRendimiento_currentIndexChanged)
        gridLayout_Calculo.addWidget(self.modelo_rendimiento,2,2,1,5)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "ΔP method", None, QtGui.QApplication.UnicodeUTF8)),3,1)
        self.modelo_DeltaP = QtGui.QComboBox()
        for txt in self.Equipment.TEXT_MODEL_DELTAP:
            self.modelo_DeltaP.addItem(txt)
        self.modelo_DeltaP.currentIndexChanged.connect(self.on_modeloDeltaP_currentIndexChanged)
        gridLayout_Calculo.addWidget(self.modelo_DeltaP,3,2,1,5)

        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),4,1,1,6)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Diameter")),5,1)
        self.diametro=Entrada_con_unidades(Length)
        self.diametro.valueChanged.connect(partial(self.changeParams, "diametro"))
        gridLayout_Calculo.addWidget(self.diametro,5,2)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Efficiency")),6,1)
        self.rendimiento=Entrada_con_unidades(float, spinbox=True)
        self.rendimiento.valueChanged.connect(partial(self.changeParams, "rendimiento"))
        gridLayout_Calculo.addWidget(self.rendimiento,6,2)

        self.groupJohnstone=QtGui.QGroupBox()
        self.groupJohnstone.setFlat(True)
        gridLayout_Calculo.addWidget(self.groupJohnstone,7,1,1,2)
        JohnstoneLayout=QtGui.QHBoxLayout(self.groupJohnstone)
        JohnstoneLayout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Ventury Constant")))
        self.k=Entrada_con_unidades(float, spinbox=True)
        self.k.valueChanged.connect(partial(self.changeParams, "k"))
        JohnstoneLayout.addWidget(self.k)

        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Length throat")),5,4)
        self.Lt=Entrada_con_unidades(Length)
        self.Lt.valueChanged.connect(partial(self.changeParams, "Lt"))
        gridLayout_Calculo.addWidget(self.Lt,5,5)

        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),8,1,1,6)
        self.groupBox_Calculo = QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Results"))
        gridLayout_Calculo.addWidget(self.groupBox_Calculo,9,1,1,5)
        gridLayout_1 = QtGui.QGridLayout(self.groupBox_Calculo)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Efficiency")),1,1)
        self.rendimientoCalc=Entrada_con_unidades(float, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.rendimientoCalc,1,2)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "DeltaP")),2,1)
        self.deltaP=Entrada_con_unidades(DeltaP, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.deltaP,2,2)

        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),11,1,1,6)

        #Pestaña salida
        self.SalidaGas= UI_corriente.Ui_corriente(readOnly=True)
        self.Salida.addTab(self.SalidaGas,QtGui.QApplication.translate("pychemqt", "Clean Gas"))
        self.SalidaLiquido= UI_corriente.Ui_corriente(readOnly=True)
        self.Salida.addTab(self.SalidaLiquido,QtGui.QApplication.translate("pychemqt", "Liquid"))

        self.on_tipoCalculo_currentIndexChanged(0)
        self.on_modeloRendimiento_currentIndexChanged(0)
        self.on_modeloDeltaP_currentIndexChanged(0)
        if equipment:
            self.setEquipment(equipment)

    def on_tipoCalculo_currentIndexChanged(self, modelo):
        self.rendimiento.setEnabled(modelo)
        self.rendimiento.setReadOnly(not modelo)
        self.diametro.setEnabled(not modelo)
        self.diametro.setReadOnly(modelo)
        self.changeParams("tipo_calculo", modelo)

    def on_modeloRendimiento_currentIndexChanged(self, modelo):
        self.groupJohnstone.setVisible(False)
        if modelo == 0:
            self.groupJohnstone.setVisible(True)
        self.changeParams("modelo_rendimiento", modelo)

    def on_modeloDeltaP_currentIndexChanged(self, modelo):
        if modelo==3:
            self.Lt.setEnabled(True)
            self.Lt.setReadOnly(False)
        else:
            self.Lt.setEnabled(False)
            self.Lt.setReadOnly(True)
        self.changeParams("modelo_DeltaP", modelo)

if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
#    from lib.corriente import Corriente, Solid
#    diametros=[17.5e-6, 22.4e-6, 26.2e-6, 31.8e-6, 37e-6, 42.4e-6, 48e-6, 54e-6, 60e-6, 69e-6, 81.3e-6, 96.5e-6, 109e-6, 127e-6]
#    fracciones=[0.02, 0.03, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.03, 0.02]
#    solido=Solid(caudalSolido=[0.1], distribucion_diametro=diametros, distribucion_fraccion=fracciones)
#    aire=Corriente(T=300, P=101325, caudalMasico=1., ids=[475], fraccionMolar=[1.], solido=solido)
#    agua=Corriente(T=300, P=101325, caudalMasico=1., ids=[62], fraccionMolar=[1.])
#    secador=Scrubber(entradaGas=aire, entradaLiquido=agua, diametro=0.05)

    dialogo = UI_equipment()
    dialogo.show()
    sys.exit(app.exec_())
