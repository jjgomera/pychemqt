#!/usr/bin/python
# -*- coding: utf-8 -*-

######################################
###      Diálogo de definición de mezcladores, UI_mixer    ###
######################################

from functools import partial

from PyQt4 import QtCore, QtGui

from lib.unidades import Pressure
from equipment.parents import UI_equip
from equipment.flux import Mixer
from UI import UI_corriente
from UI.widgets import Entrada_con_unidades


class UI_equipment(UI_equip):
    """Diálogo de definición de mezcladores"""
    Equipment=Mixer()
    def __init__(self, equipment=None, entradas=1, parent=None):
        """
        equipment: instancia de equipo inicial
        entradas: Numero de entradas
        """
        super(UI_equipment, self).__init__(Mixer, salida=False, parent=parent)

        #Pestaña entrada
        for i in range(entradas):
            entrada=UI_corriente.Ui_corriente()
            entrada.Changed.connect(partial(self.cambiarEntrada, i))
            self.entrada.addTab(entrada, str(i+1))

        #Pestaña calculo
        gridLayout_Calculo = QtGui.QGridLayout(self.tabCalculo)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Output Pressure Method")),1,1)
        self.criterio=QtGui.QComboBox()
        for txt in self.Equipment.TEXT_METODO:
            self.criterio.addItem(txt)
        self.criterio.currentIndexChanged.connect(self.criterio_Changed)
        gridLayout_Calculo.addWidget(self.criterio,1,2)

        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),2,1,1,3)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Output Pressure")),3,1)
        self.Pout=Entrada_con_unidades(Pressure)
        self.Pout.valueChanged.connect(partial(self.changeParams, "Pout"))
        gridLayout_Calculo.addWidget(self.Pout,3,2)
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),4,1,1,3)

        self.criterio_Changed(0)

        if equipment:
            self.setEquipment(equipment)
        else:
            self.Equipment=Mixer(entradas=entradas)

    def criterio_Changed(self, int):
        self.Pout.setEnabled(int==2)
        self.changeParams("criterio", int)

    def cambiarEntrada(self, ind, corriente):
        self.Equipment(id_entrada=ind, entrada=corriente)

    def rellenarInput(self):
        UI_equip.rellenarInput(self)
        for i, entrada in enumerate(self.Equipment.kwargs["entrada"]):
            if entrada:
                self.entrada.widget(i).setCorriente(entrada)


if __name__ == "__main__":
    import sys
    from lib.corriente import Corriente, Mezcla
    app = QtGui.QApplication(sys.argv)
    agua=Corriente(T=300, P=101325., caudalMasico=1, fraccionMasica=[1.])
    agua2=Corriente(T=300, P=101325.*2, caudalMasico=2, fraccionMasica=[1.])
    mezclador=Mixer(entrada=[agua, agua2], criterio=0)
#    mezclador=Mixer(criterio=0)
    dialogo = UI_equipment(mezclador, entradas=2)
    dialogo.show()
    sys.exit(app.exec_())
