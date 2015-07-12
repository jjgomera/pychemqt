#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Simplified distillation column equipment dialog
# Use the method Fenske-Underwood-Gilliland
###############################################################################

from functools import partial

from PyQt5 import QtWidgets


from lib.config import getComponents
from lib.unidades import Pressure, Volume, Length, Power, Density, Currency
from tools.costIndex import CostData
from equipment.parents import UI_equip
from equipment.distillation import ColumnFUG
from UI.widgets import Entrada_con_unidades


class UI_equipment(UI_equip):
    """ColumnFUG equipment edition dialog"""
    Equipment = ColumnFUG()

    def __init__(self, equipment=None, parent=None):
        """
        equipment: Initial equipment instance to model
        """
        super(UI_equipment, self).__init__(ColumnFUG, entrada=False,
                                           parent=parent)

        # Calculate tab
        lyt = QtWidgets.QGridLayout(self.tabCalculo)
        lyt.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Feed tray")), 2, 0)
        self.feed = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_FEED:
            self.feed.addItem(txt)
        self.feed.currentIndexChanged.connect(
            partial(self.changeParams, "feed"))
        lyt.addWidget(self.feed, 2, 1)

        lyt.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Condenser")), 3, 0)
        self.condenser = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_CONDENSER:
            self.condenser.addItem(txt)
        self.condenser.currentIndexChanged.connect(
            partial(self.changeParams, "condenser"))
        lyt.addWidget(self.condenser, 3, 1)
        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            4, 0, 1, 5)

        group = QtWidgets.QGroupBox(QtWidgets.QApplication.translate(
            "pychemqt", "Key Components specification"))
        lyt.addWidget(group, 5, 0, 1, 5)
        layout = QtWidgets.QGridLayout(group)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Light")), 1, 1)
        self.LK = QtWidgets.QComboBox()
        layout.addWidget(self.LK, 1, 2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Split in destilate")), 1, 4)
        self.LKsplit = Entrada_con_unidades(float, spinbox=True, max=1.)
        self.LKsplit.valueChanged.connect(
            partial(self.changeParams, "LKsplit"))
        layout.addWidget(self.LKsplit, 1, 5)

        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Heavy")), 2, 1)
        self.HK = QtWidgets.QComboBox()
        layout.addWidget(self.HK, 2, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            40, 10, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed), 2, 3)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Spit in residue")), 2, 4)
        self.HKsplit = Entrada_con_unidades(float, spinbox=True, max=1.)
        self.HKsplit.valueChanged.connect(
            partial(self.changeParams, "HKsplit"))
        layout.addWidget(self.HKsplit, 2, 5)

        indices, nombres, M = getComponents()
        for i, nombre in enumerate(nombres):
            self.LK.addItem("%i - %s" % (i+1, nombre))
            self.HK.addItem("%i - %s" % (i+1, nombre))
        self.LK.setCurrentIndex(-1)
        self.HK.setCurrentIndex(-1)
        self.LK.currentIndexChanged.connect(partial(self.changeParams, "LK"))
        self.HK.currentIndexChanged.connect(partial(self.changeParams, "HK"))
        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            6, 0, 1, 5)

        lyt.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Reflux ratio")), 7, 0)
        self.R = Entrada_con_unidades(float)
        self.R.valueChanged.connect(partial(self.changeParams, "R"))
        lyt.addWidget(self.R, 7, 1)
        lyt.addWidget(QtWidgets.QLabel("R/Rmin"), 8, 0)
        self.R_Rmin = Entrada_con_unidades(float)
        self.R_Rmin.valueChanged.connect(partial(self.changeParams, "R_Rmin"))
        lyt.addWidget(self.R_Rmin, 8, 1)

        lyt.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Design Pressure")), 7, 3)
        self.Pd = Entrada_con_unidades(Pressure)
        self.Pd.valueChanged.connect(partial(self.changeParams, "Pd"))
        lyt.addWidget(self.Pd, 7, 4)
        lyt.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Pressure loss")), 8, 3)
        self.DeltaP = Entrada_con_unidades(Pressure)
        self.DeltaP.valueChanged.connect(partial(self.changeParams, "DeltaP"))
        lyt.addWidget(self.DeltaP, 8, 4)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding),
            9, 0, 1, 5)
        self.buttonMcCabe = QtWidgets.QPushButton(
            QtWidgets.QApplication.translate("pychemqt", "McCabe-Thiele"))
        self.buttonMcCabe.clicked.connect(self.mcCabe)
        lyt.addWidget(self.buttonMcCabe, 10, 0)

        groupBox_Calculo = QtWidgets.QGroupBox(
            QtWidgets.QApplication.translate("pychemqt", "Results"))
        lyt.addWidget(groupBox_Calculo, 11, 0, 1, 5)
        layout = QtWidgets.QGridLayout(groupBox_Calculo)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Condenser Duty")), 0, 1)
        self.DutyCondenser = Entrada_con_unidades(Power, retornar=False)
        self.DutyCondenser.setReadOnly(True)
        layout.addWidget(self.DutyCondenser, 0, 2)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Reboiler Duty")), 1, 1)
        self.DutyReboiler = Entrada_con_unidades(Power, retornar=False)
        self.DutyReboiler.setReadOnly(True)
        layout.addWidget(self.DutyReboiler, 1, 2)
        layout.addWidget(QtWidgets.QLabel("Rmin"), 2, 1)
        self.Rmin = Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.Rmin, 2, 2)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Reflux ratio")), 3, 1)
        self.RCalculada = Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.RCalculada, 3, 2)

        layout.addWidget(QtWidgets.QLabel("Nmin"), 0, 4)
        self.Nmin = Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.Nmin, 0, 5)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Stages")), 1, 4)
        self.NTray = Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.NTray, 1, 5)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Feed stage")), 2, 4)
        self.N_feed = Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.N_feed, 2, 5)

        # Cost tab
        lyt = QtWidgets.QGridLayout(self.tabCostos)
        lyt.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Process")), 1, 1)
        self.proceso = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_PROCESS:
            self.proceso.addItem(txt)
        self.proceso.currentIndexChanged.connect(
            partial(self.changeParamsCoste, "proceso"))
        lyt.addWidget(self.proceso, 1, 2)
        lyt.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Column tipe")), 2, 1)
        self.tipo = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_COLUMN:
            self.tipo.addItem(txt)
        self.tipo.currentIndexChanged.connect(self.mostrarSubclasificacion)
        lyt.addWidget(self.tipo, 2, 2)
        lyt.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Material")), 3, 1)
        self.material = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_MATERIAL:
            self.material.addItem(txt)
        self.material.currentIndexChanged.connect(
            partial(self.changeParamsCoste, "material_columna"))
        lyt.addWidget(self.material, 3, 2)

        lyt.addItem(QtWidgets.QSpacerItem(
            30, 30, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            1, 3, 5, 1)

        self.groupBox_Pisos = QtWidgets.QGroupBox(
            QtWidgets.QApplication.translate("pychemqt", "Tray column"))
        lyt.addWidget(self.groupBox_Pisos, 1, 4, 4, 2)
        layout = QtWidgets.QGridLayout(self.groupBox_Pisos)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Tray type")), 1, 1)
        self.tipoPisos = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_TRAY:
            self.tipoPisos.addItem(txt)
        self.tipoPisos.currentIndexChanged.connect(
            partial(self.changeParamsCoste, "tipo_pisos"))
        layout.addWidget(self.tipoPisos, 1, 2)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Material")), 2, 1)
        self.materialPisos = QtWidgets.QComboBox()
        for txt in self.Equipment.TEXT_MATERIAL:
            self.materialPisos.addItem(txt)
        self.materialPisos.currentIndexChanged.connect(
            partial(self.changeParamsCoste, "material_pisos"))
        layout.addWidget(self.materialPisos, 2, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            3, 1, 1, 2)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Diameter")), 4, 1)
        self.diametroPisos = Entrada_con_unidades(Length)
        layout.addWidget(self.diametroPisos, 4, 2)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Stages")), 5, 1)
        self.NumeroPisos = Entrada_con_unidades(int, spinbox=True, min=1,
                                                step=1, width=50)
        layout.addWidget(self.NumeroPisos, 5, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            6, 1, 1, 2)

        self.groupBox_relleno = QtWidgets.QGroupBox(
            QtWidgets.QApplication.translate("pychemqt", "Packed column"))
        lyt.addWidget(self.groupBox_relleno, 1, 4, 4, 2)
        layout = QtWidgets.QGridLayout(self.groupBox_relleno)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Volume")), 1, 1)
        self.VolumenRelleno = Entrada_con_unidades(Volume, "VolLiq")
        layout.addWidget(self.VolumenRelleno, 1, 2)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Unit Cost")), 2, 1)
        texto = "%s / %s" % (Currency(None).text(), Volume(None).text("VolLiq"))
        self.C_unit_relleno = Entrada_con_unidades(Currency, retornar=False,
                                                   textounidad=texto)
        layout.addWidget(self.C_unit_relleno, 2, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            3, 1, 1, 2)

        lyt.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Diameter")), 5, 1)
        self.Dc = Entrada_con_unidades(Length)
        lyt.addWidget(self.Dc, 5, 2, 1, 2)
        lyt.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Height")), 6, 1)
        self.Hc = Entrada_con_unidades(Length)
        lyt.addWidget(self.Hc, 6, 2, 1, 2)
        lyt.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Thickness (top)")), 6, 4)
        self.EspesorSuperior = Entrada_con_unidades(Length, "Thickness")
        lyt.addWidget(self.EspesorSuperior, 6, 5, 1, 2)
        lyt.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Thickness (bottom)")), 7, 4)
        self.EspesorInferior = Entrada_con_unidades(Length, "Thickness")
        lyt.addWidget(self.EspesorInferior, 7, 5, 1, 2)
        lyt.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Density")), 7, 1)
        self.EspesorInferior = Entrada_con_unidades(Density, "DenLiq")
        lyt.addWidget(self.EspesorInferior, 7, 2, 1, 2)

        self.Costos = CostData(self.Equipment)
        self.Costos.valueChanged.connect(self.changeParamsCoste)
        lyt.addWidget(self.Costos, 10, 1, 2, 5)

        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding),
            12, 1, 1, 6)

        self.groupBox_Costos = QtWidgets.QGroupBox(
            QtWidgets.QApplication.translate("pychemqt", "Stimated costs"))
        lyt.addWidget(self.groupBox_Costos, 13, 1, 1, 5)
        layout = QtWidgets.QGridLayout(self.groupBox_Costos)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Tray cost")), 0, 1)
        self.C_pisos = Entrada_con_unidades(Currency, retornar=False)
        self.C_pisos.setReadOnly(True)
        layout.addWidget(self.C_pisos, 0, 2)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Shell cost")), 1, 1)
        self.C_carcasa = Entrada_con_unidades(Currency, retornar=False)
        self.C_carcasa.setReadOnly(True)
        layout.addWidget(self.C_carcasa, 1, 2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Platform and ladder")), 2, 1)
        self.C_accesorios = Entrada_con_unidades(Currency, retornar=False)
        self.C_accesorios.setReadOnly(True)
        layout.addWidget(self.C_accesorios, 2, 2)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Column cost")), 0, 4)
        self.C_columna = Entrada_con_unidades(Currency, retornar=False)
        self.C_columna.setReadOnly(True)
        layout.addWidget(self.C_columna, 0, 5)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Purchase costs")), 1, 4)
        self.C_adq = Entrada_con_unidades(Currency, retornar=False)
        self.C_adq.setReadOnly(True)
        layout.addWidget(self.C_adq, 1, 5)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Installed costs")), 2, 4)
        self.C_inst = Entrada_con_unidades(Currency, retornar=False)
        self.C_inst.setReadOnly(True)
        layout.addWidget(self.C_inst, 2, 5)

        # Output tab
        self.addSalida(QtWidgets.QApplication.translate("pychemqt", "Destilate"))
        self.addSalida(QtWidgets.QApplication.translate("pychemqt", "Residue"))

        self.mostrarSubclasificacion(0)
        if equipment:
            self.setEquipment(equipment)

    def mcCabe(self):
        self.Equipment.McCabe()

    def mostrarSubclasificacion(self, ind):
        self.groupBox_Pisos.setVisible(not ind)
        self.groupBox_relleno.setVisible(ind)
        self.changeParamsCoste("tipo", ind)

    def rellenar(self):
        self.buttonMcCabe.setEnabled(self.Equipment.statusMcCabe)
        UI_equip.rellenar(self)


if __name__ == "__main__":
    import sys
    from lib.corriente import Corriente
    app = QtWidgets.QApplication(sys.argv)
    blend = Corriente(T=340, P=101325, caudalMasico=1, ids=[11, 12],
                      fraccionMolar=[0.6, 0.4])
    columna = ColumnFUG(entrada=blend, LK=0, LKsplit=0.96666, HK=1,
                        HKsplit=0.95, feed=0)
    dialogo = UI_equipment(columna)
    dialogo.show()
    sys.exit(app.exec_())
