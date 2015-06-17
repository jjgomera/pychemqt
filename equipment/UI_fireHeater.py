#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Fireheater equipment dialog
###############################################################################

from functools import partial

from PyQt4 import QtGui

from lib.unidades import Temperature, Pressure, Power, VolFlow, Currency
from tools.costIndex import CostData
from equipment.parents import UI_equip
from equipment.heatExchanger import Fired_Heater
from UI.widgets import Entrada_con_unidades


class UI_equipment(UI_equip):
    """Fireheater equipment edition dialog"""
    Equipment = Fired_Heater()

    def __init__(self, equipment=None, parent=None):
        """
        equipment: Initial equipment instance to model
        """
        super(UI_equipment, self).__init__(Fired_Heater, entrada=False,
                                           salida=False, parent=parent)

        # Calculate tab
        layout = QtGui.QGridLayout(self.tabCalculo)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Output Temperature")), 1, 1)
        self.Tout = Entrada_con_unidades(Temperature, resaltado=True)
        self.Tout.valueChanged.connect(partial(self.changeParams, "Tout"))
        layout.addWidget(self.Tout, 1, 2)
        layout.addItem(QtGui.QSpacerItem(
            10, 10, QtGui.QSizePolicy. Fixed, QtGui.QSizePolicy.Fixed),
            2, 0, 1, 6)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Pressure drop")), 3, 1)
        self.deltaP = Entrada_con_unidades(Pressure)
        self.deltaP.valueChanged.connect(partial(self.changeParams, "deltaP"))
        layout.addWidget(self.deltaP, 3, 2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Maximum heat flux")), 4, 1)
        self.Hmax = Entrada_con_unidades(Power)
        self.Hmax.valueChanged.connect(partial(self.changeParams, "Hmax"))
        layout.addWidget(self.Hmax, 4, 2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Fuel calorific value")), 5, 1)
        self.poderCalorifico = Entrada_con_unidades(float)
        self.poderCalorifico.valueChanged.connect(
            partial(self.changeParams, "poderCalorifico"))
        layout.addWidget(self.poderCalorifico, 5, 2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Efficiency")), 6, 1)
        self.eficiencia = Entrada_con_unidades(float, spinbox=True)
        self.eficiencia.valueChanged.connect(
            partial(self.changeParams, "eficiencia"))
        layout.addWidget(self.eficiencia, 6, 2)
        layout.addItem(QtGui.QSpacerItem(
            10, 10, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding),
            7, 0, 1, 6)

        group = QtGui.QGroupBox(
            QtGui.QApplication.translate("pychemqt", "Results"))
        layout.addWidget(group, 8, 1, 1, 5)
        lyt = QtGui.QGridLayout(group)
        lyt.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Heat")), 0, 1)
        self.Heat = Entrada_con_unidades(Power, retornar=False, readOnly=True)
        lyt.addWidget(self.Heat, 0, 2)
        lyt.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Fuel")), 1, 1)
        self.CombustibleRequerido = Entrada_con_unidades(
            VolFlow, "QLiq", retornar=False, readOnly=True)
        lyt.addWidget(self.CombustibleRequerido, 1, 2)
        layout.addItem(QtGui.QSpacerItem(
            20, 20, QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed),
            9, 0, 1, 6)

        # Cost tab
        lyt_Cost = QtGui.QGridLayout(self.tabCostos)
        lyt_Cost.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Type")), 1, 1)
        self.tipo = QtGui.QComboBox()
        for txt in self.Equipment.TEXT_TIPO:
            self.tipo.addItem(txt)
        self.tipo.currentIndexChanged.connect(self.mostrarSubclasificacion)
        lyt_Cost.addWidget(self.tipo, 1, 2)
        self.label = QtGui.QLabel()
        lyt_Cost.addWidget(self.label, 2, 1)
        self.subtipoBox = QtGui.QComboBox()
        for txt in self.Equipment.TEXT_SUBTIPOBOX:
            self.subtipoBox.addItem(txt)
        self.subtipoBox.currentIndexChanged.connect(
            partial(self.changeParamsCoste, "subtipoBox"))
        lyt_Cost.addWidget(self.subtipoBox, 2, 2)
        self.subtipoCylindrical = QtGui.QComboBox()
        for txt in self.Equipment.TEXT_SUBTIPOCYLINDRICAL:
            self.subtipoCylindrical.addItem(txt)
        self.subtipoCylindrical.currentIndexChanged.connect(
            partial(self.changeParamsCoste, "subtipoCylindrical"))
        lyt_Cost.addWidget(self.subtipoCylindrical, 2, 2)
        lyt_Cost.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Material")), 3, 1)
        self.material = QtGui.QComboBox()
        for txt in self.Equipment.TEXT_MATERIAL:
            self.material.addItem(txt)
        self.material.currentIndexChanged.connect(
            partial(self.changeParamsCoste, "material"))
        lyt_Cost.addWidget(self.material, 3, 2)
        lyt_Cost.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Design pressure")), 4, 1)
        self.P_dis = Entrada_con_unidades(Pressure)
        self.P_dis.valueChanged.connect(partial(self.changeParamsCoste, "P_dis"))
        lyt_Cost.addWidget(self.P_dis, 4, 2)
        lyt_Cost.addItem(QtGui.QSpacerItem(
            20, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding),
            5, 1, 1, 6)

        self.Costos = CostData(self.Equipment)
        self.Costos.valueChanged.connect(self.changeParamsCoste)
        lyt_Cost.addWidget(self.Costos, 6, 1, 2, 5)

        lyt_Cost.addItem(QtGui.QSpacerItem(
            20, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding),
            8, 1, 1, 6)
        group = QtGui.QGroupBox(
            QtGui.QApplication.translate("pychemqt", "Stimated Costs"))
        lyt_Cost.addWidget(group, 9, 1, 1, 6)
        layout = QtGui.QGridLayout(group)
        layout.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Purchase costs")), 0, 1)
        self.C_adq = Entrada_con_unidades(Currency, retornar=False)
        self.C_adq.setReadOnly(True)
        layout.addWidget(self.C_adq, 0, 2)
        layout.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Installed costs")), 1, 1)
        self.C_inst = Entrada_con_unidades(Currency, retornar=False)
        self.C_inst.setReadOnly(True)
        layout.addWidget(self.C_inst, 1, 2)
        lyt_Cost.addItem(QtGui.QSpacerItem(
            20, 20, QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed),
            10, 1, 1, 6)

        self.mostrarSubclasificacion(0)
        if equipment:
            self.setEquipment(equipment)

    def mostrarSubclasificacion(self, ind):
        if ind:
            txt = QtGui.QApplication.translate("pychemqt",
                                               "Cylindrical heater type")
        else:
            txt = QtGui.QApplication.translate("pychemqt", "Box heater type")
        self.label.setText(txt)
        self.subtipoBox.setVisible(not ind)
        self.subtipoCylindrical.setVisible(ind)
        self.changeParamsCoste("tipo", ind)


if __name__ == "__main__":
    import sys
    from lib.corriente import Corriente
    from equipment.heatExchanger import Fired_Heater
    app = QtGui.QApplication(sys.argv)
    agua = Corriente(T=300, P=101325., caudalMasico=0.01, fraccionMasica=[1.])
    fireheater = Fired_Heater(entrada=agua, Tout=450)
    dialogo = UI_equipment(fireheater)
    dialogo.show()
    sys.exit(app.exec_())
