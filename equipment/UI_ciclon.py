#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Cyclone equipment dialog
###############################################################################

from functools import partial
import os

from PyQt4 import QtGui

from lib.unidades import Length, Pressure, DeltaP, Speed, VolFlow, Currency
from tools.costIndex import CostData
from equipment.gas_solid import Ciclon, UI_equipment_Solid
from UI import UI_corriente
from UI.widgets import Entrada_con_unidades


class UI_equipment(UI_equipment_Solid):
    """Cyclone equipment edition dialog"""
    Equipment = Ciclon()

    def __init__(self, equipment=None, parent=None):
        """
        equipment: Initial equipment instance to model
        """
        super(UI_equipment, self).__init__(Ciclon, entrada=False, parent=parent)

        # Calculate tab
        lyt_Calc = QtGui.QGridLayout(self.tabCalculo)
        lyt_Calc.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Mode")), 0, 1, 1, 2)
        self.tipo_calculo = QtGui.QComboBox()
        for txt in self.Equipment.TEXT_TIPO:
            self.tipo_calculo.addItem(txt)
        self.tipo_calculo.currentIndexChanged.connect(
            self.on_tipoCalculo_currentIndexChanged)
        lyt_Calc.addWidget(self.tipo_calculo, 0, 3, 1, 4)
        lyt_Calc.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Method")), 1, 1, 1, 2)
        self.modelo_rendimiento = QtGui.QComboBox()
        for txt in self.Equipment.TEXT_MODEL:
            self.modelo_rendimiento.addItem(txt)
        self.modelo_rendimiento.currentIndexChanged.connect(
            partial(self.changeParams, "modelo_rendimiento"))
        lyt_Calc.addWidget(self.modelo_rendimiento, 1, 3, 1, 4)
        lyt_Calc.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "ΔP method", None, QtGui.QApplication.UnicodeUTF8)),
            2, 1, 1, 2)
        self.modelo_DeltaP = QtGui.QComboBox()
        for txt in self.Equipment.TEXT_MODEL_DELTAP:
            self.modelo_DeltaP.addItem(txt)
        self.modelo_DeltaP.currentIndexChanged.connect(
            partial(self.changeParams, "modelo_DeltaP"))
        lyt_Calc.addWidget(self.modelo_DeltaP, 2, 3, 1, 4)
        lyt_Calc.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Design model")), 3, 1, 1, 2)
        self.modelo_ciclon = QtGui.QComboBox()
        for txt in self.Equipment.TEXT_MODEL_CICLON:
            self.modelo_ciclon.addItem(txt)
        self.modelo_ciclon.currentIndexChanged.connect(
            self.modeloEficiencia_Changed)
        lyt_Calc.addWidget(self.modelo_ciclon, 3, 3, 1, 4)
        lyt_Calc.addItem(QtGui.QSpacerItem(
            20, 20, QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed),
            4, 0, 1, 5)

        lyt_Calc.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Diameter")), 5, 1)
        self.Dc = Entrada_con_unidades(Length)
        self.Dc.valueChanged.connect(partial(self.changeParams, "Dc"))
        lyt_Calc.addWidget(self.Dc, 5, 2, 1, 2)
        lyt_Calc.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Efficiency")), 6, 1)
        self.rendimientoAdmisible = Entrada_con_unidades(float, spinbox=True)
        self.rendimientoAdmisible.valueChanged.connect(
            partial(self.changeParams, "rendimientoAdmisible"))
        lyt_Calc.addWidget(self.rendimientoAdmisible, 6, 2, 1, 2)
        lyt_Calc.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Allowable ΔP", None, QtGui.QApplication.UnicodeUTF8)),
            7, 1)
        self.DeltaPAdmisible = Entrada_con_unidades(Pressure)
        self.DeltaPAdmisible.valueChanged.connect(
            partial(self.changeParams, "DeltaPAdmisible"))
        lyt_Calc.addWidget(self.DeltaPAdmisible, 7, 2, 1, 2)
        lyt_Calc.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Allowable speed")), 8, 1)
        self.velocidadAdmisible = Entrada_con_unidades(Speed)
        self.velocidadAdmisible.valueChanged.connect(
            partial(self.changeParams, "velocidadAdmisible"))
        lyt_Calc.addWidget(self.velocidadAdmisible, 8, 2, 1, 2)
        lyt_Calc.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "No. of ciclones")), 9, 1)
        self.num_ciclones = Entrada_con_unidades(int, spinbox=True, step=1,
                                                 decimales=0, min=1)
        self.num_ciclones.valueChanged.connect(
            partial(self.changeParams, "num_ciclones"))
        lyt_Calc.addWidget(self.num_ciclones, 9, 2, 1, 2)
        lyt_Calc.addItem(QtGui.QSpacerItem(
            20, 20, QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed),
            10, 0, 1, 4)

        group = QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt",
                                                             "Results"))
        lyt_Calc.addWidget(group, 11, 1, 1, 3)
        lyt = QtGui.QGridLayout(group)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Efficiency")), 1, 1)
        self.rendimientoCalc = Entrada_con_unidades(float, readOnly=True)
        lyt.addWidget(self.rendimientoCalc, 1, 2)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Pressure drop:")), 2, 1)
        self.deltaP = Entrada_con_unidades(DeltaP, readOnly=True)
        lyt.addWidget(self.deltaP, 2, 2)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Admission speed")), 3, 1)
        self.V = Entrada_con_unidades(Speed, readOnly=True)
        lyt.addWidget(self.V, 3, 2)

        lyt_Calc.addItem(QtGui.QSpacerItem(
            20, 20, QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed),
            5, 4, 7, 1)
        group2 = QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt",
                                                              "Geometry"))
        lyt_Calc.addWidget(group2, 5, 5, 7, 1)
        lyt = QtGui.QGridLayout(group2)
        lyt.addWidget(QtGui.QLabel("D<sub>C</sub>"), 1, 1)
        self.Dcc = Entrada_con_unidades(Length, boton=False)
        lyt.addWidget(self.Dcc, 1, 2)
        lyt.addWidget(QtGui.QLabel("B<sub>C</sub>"), 2, 1)
        self.Bc = Entrada_con_unidades(Length, boton=False)
        lyt.addWidget(self.Bc, 2, 2)
        lyt.addWidget(QtGui.QLabel("H<sub>C</sub>"), 3, 1)
        self.Hc = Entrada_con_unidades(Length, boton=False)
        lyt.addWidget(self.Hc, 3, 2)
        lyt.addWidget(QtGui.QLabel("J<sub>C</sub>"), 4, 1)
        self.Jc = Entrada_con_unidades(Length, boton=False)
        lyt.addWidget(self.Jc, 4, 2)
        lyt.addWidget(QtGui.QLabel("L<sub>C</sub>"), 5, 1)
        self.Lc = Entrada_con_unidades(Length, boton=False)
        lyt.addWidget(self.Lc, 5, 2)
        lyt.addWidget(QtGui.QLabel("Z<sub>C</sub>"), 6, 1)
        self.Zc = Entrada_con_unidades(Length, boton=False)
        lyt.addWidget(self.Zc, 6, 2)
        lyt.addWidget(QtGui.QLabel("D<sub>e</sub>"), 7, 1)
        self.De = Entrada_con_unidades(Length, boton=False)
        lyt.addWidget(self.De, 7, 2)
        lyt.addWidget(QtGui.QLabel("S<sub>C</sub>"), 8, 1)
        self.Sc = Entrada_con_unidades(Length, boton=False)
        lyt.addWidget(self.Sc, 8, 2)
        lyt.addWidget(QtGui.QLabel("N"), 9, 1)
        self.NCalc = Entrada_con_unidades(int, readOnly=True)
        lyt.addWidget(self.NCalc, 9, 2)
        lyt_Calc.addItem(QtGui.QSpacerItem(
            20, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding),
            0, 7, 10, 1)
        lyt_Calc.addItem(QtGui.QSpacerItem(
            20, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding),
            5, 6, 10, 1)

        image = QtGui.QLabel()
        path = os.environ["pychemqt"]+"/images/equip/ciclon.gif"
        image.setPixmap(QtGui.QPixmap(path))
        image.setScaledContents(True)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding,
                                       QtGui.QSizePolicy.Expanding)
        image.setSizePolicy(sizePolicy)
        lyt_Calc.addWidget(image, 0, 8, 12, 1)

        # Cost tab
        lyt_Cost = QtGui.QGridLayout(self.tabCostos)
        lyt_Cost.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Model")), 1, 1)
        self.tipo_costo = QtGui.QComboBox()
        for txt in self.Equipment.TEXT_COST:
            self.tipo_costo.addItem(txt)
        self.tipo_costo.currentIndexChanged.connect(
            partial(self.changeParamsCoste, "tipo_costo"))
        lyt_Cost.addWidget(self.tipo_costo, 1, 2)
        lyt_Cost.addItem(QtGui.QSpacerItem(
            10, 10, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Fixed), 1, 3)

        self.Costos = CostData(self.Equipment)
        self.Costos.valueChanged.connect(self.changeParamsCoste)
        lyt_Cost.addWidget(self.Costos, 3, 1, 1, 3)
        lyt_Cost.addItem(QtGui.QSpacerItem(
            20, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding),
            4, 1, 1, 3)

        group = QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt",
                                                             "Stimated Costs"))
        lyt_Cost.addWidget(group, 5, 1, 1, 3)
        lyt = QtGui.QGridLayout(group)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Number")), 0, 0)
        self.num_ciclonesCoste = Entrada_con_unidades(int, readOnly=True)
        lyt.addWidget(self.num_ciclonesCoste, 0, 1, 1, 3)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Volumetric Flow")), 1, 0)
        self.Q = Entrada_con_unidades(VolFlow, "QGas", retornar=False)
        self.Q.setReadOnly(True)
        lyt.addWidget(self.Q, 1, 1)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Purchase cost")), 0, 3)
        self.C_adq = Entrada_con_unidades(Currency, retornar=False,
                                          decimales=2, tolerancia=8)
        self.C_adq.setReadOnly(True)
        lyt.addWidget(self.C_adq, 0, 4)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Installed cost")), 1, 3)
        self.C_inst = Entrada_con_unidades(Currency, retornar=False,
                                           decimales=2, tolerancia=8)
        self.C_inst.setReadOnly(True)
        lyt.addWidget(self.C_inst, 1, 4)

        # Output tab
        self.SalidaGas = UI_corriente.Ui_corriente(readOnly=True)
        self.SalidaSolido = UI_corriente.Ui_corriente(readOnly=True)
        self.Salida.addTab(
            self.SalidaGas,
            QtGui.QApplication.translate("pychemqt", "Filtered gas"))
        self.Salida.addTab(
            self.SalidaSolido,
            QtGui.QApplication.translate("pychemqt", "Collected solids"))

        self.on_tipoCalculo_currentIndexChanged(0)
        self.modeloEficiencia_Changed(0)
        if equipment:
            self.setEquipment(equipment)

    def on_tipoCalculo_currentIndexChanged(self, int):
        """Habilita o desabilita los datos requeridos para el cálculo"""
        if int and self.modelo_ciclon.count() == 9:
            self.modelo_ciclon.removeItem(8)
        elif not int and self.modelo_ciclon.count() == 8:
            self.modelo_ciclon.addItem(QtGui.QApplication.translate("pychemqt",
                                                                    "Custom"))

        self.Dc.setReadOnly(int)
        self.Dc.setResaltado(not int)
        self.num_ciclones.setReadOnly(int)
        self.num_ciclones.setResaltado(not int)
        self.rendimientoAdmisible.setReadOnly(not int)
        self.rendimientoAdmisible.setResaltado(int)
        self.velocidadAdmisible.setReadOnly(not int)
        self.velocidadAdmisible.setResaltado(False)
        self.DeltaPAdmisible.setReadOnly(not int)
        self.DeltaPAdmisible.setResaltado(False)
        self.changeParams("tipo_calculo", int)

    def modeloEficiencia_Changed(self, modelo):
        if modelo == 8:
            # Customized model, let user edit dimensions values
            self.Hc.setReadOnly(False)
            self.Hc.setResaltado(True)
            self.Bc.setReadOnly(False)
            self.Bc.setResaltado(True)
            self.Jc.setReadOnly(False)
            self.Jc.setResaltado(False)
            self.Lc.setReadOnly(False)
            self.Lc.setResaltado(False)
            self.Zc.setReadOnly(False)
            self.Zc.setResaltado(False)
            self.De.setReadOnly(False)
            self.De.setResaltado(True)
            self.Sc.setReadOnly(False)
            self.Sc.setResaltado(False)
        else:
            self.Hc.setReadOnly(True)
            self.Hc.setResaltado(False)
            self.Bc.setReadOnly(True)
            self.Bc.setResaltado(False)
            self.Jc.setReadOnly(True)
            self.Jc.setResaltado(False)
            self.Lc.setReadOnly(True)
            self.Lc.setResaltado(False)
            self.Zc.setReadOnly(True)
            self.Zc.setResaltado(False)
            self.De.setReadOnly(True)
            self.De.setResaltado(False)
            self.Sc.setReadOnly(True)
            self.Sc.setResaltado(False)
        self.Dcc.setReadOnly(True)
        self.Dcc.setResaltado(False)
        self.changeParams("modelo_ciclon", modelo)


if __name__ == "__main__":
    import sys
    from lib.corriente import Corriente, Solid
    app = QtGui.QApplication(sys.argv)
    diametros = [17.5e-6, 22.4e-6, 26.2e-6, 31.8e-6, 37e-6, 42.4e-6, 48e-6,
                 54e-6, 60e-6, 69e-6, 81.3e-6, 96.5e-6, 109e-6, 127e-6]
    fracciones = [0.02, 0.03, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                  0.05, 0.03, 0.02]
    solido = Solid(caudalSolido=[0.1], distribucion_diametro=diametros,
                   distribucion_fraccion=fracciones)
    corriente = Corriente(T=300, P=101325, caudalMasico=1., fraccionMolar=[1.],
                          solido=solido)

    ciclon = Ciclon(entrada=corriente, tipo_calculo=1, velocidadAdmisible=5)
    Ciclon = UI_equipment(ciclon)
    Ciclon.show()
    sys.exit(app.exec_())
