#!/usr/bin/python
# -*- coding: utf-8 -*-

#######################################################################
###                                               Diálogo de definición de tuberías, UI_pipe                                                     ###
#######################################################################

from PyQt4 import QtCore, QtGui

from equipment.distillation import Tower
from UI import UI_corriente
from equipment import parents
from lib import unidades, config
from tools import costIndex
from UI.widgets import Entrada_con_unidades


class UI_equipment(parents.UI_equip):
    """Diálogo de definición de tuberías"""
    def __init__(self, entrada=None, parent=None):
        """entrada: Parametro opcional de clase corriente que indica la corriente de entrada en kla tubería"""
        super(UI_equipment, self).__init__(Tower, parent=parent)
        self.entrada=entrada

#        #Pestaña entrada
#        self.Entrada= UI_corriente.Ui_corriente(entrada)
#        self.Entrada.Changed.connect(self.cambiar_entrada)
#        self.tabWidget.addTab(self.Entrada, QtGui.QApplication.translate("equipment", "Entrada", None, QtGui.QApplication.UnicodeUTF8))

        #Pestaña calculo
        gridLayout_Calculo = QtGui.QGridLayout(self.tabCalculo)

        #Pestaña costos
        gridLayout_Costos = QtGui.QGridLayout(self.tabCostos)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Proceso:", None, QtGui.QApplication.UnicodeUTF8)), 1, 1, 1, 1)
        self.proceso=QtGui.QComboBox()
        self.proceso.addItem(QtGui.QApplication.translate("equipment", "Destilación", None, QtGui.QApplication.UnicodeUTF8))
        self.proceso.addItem(QtGui.QApplication.translate("equipment", "Absorción", None, QtGui.QApplication.UnicodeUTF8))
        self.proceso.currentIndexChanged.connect(self.calcularCostos)
        gridLayout_Costos.addWidget(self.proceso, 1, 2, 1, 1)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Tipo de columna:", None, QtGui.QApplication.UnicodeUTF8)), 2, 1, 1, 1)
        self.tipo=QtGui.QComboBox()
        self.tipo.addItem(QtGui.QApplication.translate("equipment", "De pisos", None, QtGui.QApplication.UnicodeUTF8))
        self.tipo.addItem(QtGui.QApplication.translate("equipment", "De relleno", None, QtGui.QApplication.UnicodeUTF8))
        self.tipo.currentIndexChanged.connect(self.mostrarSubclasificacion)
        self.tipo.currentIndexChanged.connect(self.calcularCostos)
        gridLayout_Costos.addWidget(self.tipo, 2, 2, 1, 1)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Material:", None, QtGui.QApplication.UnicodeUTF8)), 3, 1, 1, 1)
        self.material=QtGui.QComboBox()
        self.material.addItem(QtGui.QApplication.translate("equipment", "Acero al carbon", None, QtGui.QApplication.UnicodeUTF8))
        self.material.addItem(QtGui.QApplication.translate("equipment", "Acero inoxidable 304", None, QtGui.QApplication.UnicodeUTF8))
        self.material.addItem(QtGui.QApplication.translate("equipment", "Acero inoxidable 316", None, QtGui.QApplication.UnicodeUTF8))
        self.material.addItem(QtGui.QApplication.translate("equipment", "Carpenter 20CB-3", None, QtGui.QApplication.UnicodeUTF8))
        self.material.addItem(QtGui.QApplication.translate("equipment", "Niquel 200", None, QtGui.QApplication.UnicodeUTF8))
        self.material.addItem(QtGui.QApplication.translate("equipment", "Monel 400", None, QtGui.QApplication.UnicodeUTF8))
        self.material.addItem(QtGui.QApplication.translate("equipment", "Inconel 600", None, QtGui.QApplication.UnicodeUTF8))
        self.material.addItem(QtGui.QApplication.translate("equipment", "Incoloy 825", None, QtGui.QApplication.UnicodeUTF8))
        self.material.addItem(QtGui.QApplication.translate("equipment", "Titanio", None, QtGui.QApplication.UnicodeUTF8))
        self.material.currentIndexChanged.connect(self.calcularCostos)
        gridLayout_Costos.addWidget(self.material, 3, 2, 1, 1)

        gridLayout_Costos.addItem(QtGui.QSpacerItem(30,30,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),1,3,5,1)

        self.groupBox_Pisos = QtGui.QGroupBox(QtGui.QApplication.translate("equipment", "Torre de pisos", None, QtGui.QApplication.UnicodeUTF8))
        gridLayout_Costos.addWidget(self.groupBox_Pisos,1,4,4,2)
        gridLayout_1 = QtGui.QGridLayout(self.groupBox_Pisos)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Tipo:", None, QtGui.QApplication.UnicodeUTF8)), 1, 1, 1, 1)
        self.tipoPisos=QtGui.QComboBox()
        self.tipoPisos.addItem(QtGui.QApplication.translate("equipment", "De válvula", None, QtGui.QApplication.UnicodeUTF8))
        self.tipoPisos.addItem(QtGui.QApplication.translate("equipment", "De rejilla", None, QtGui.QApplication.UnicodeUTF8))
        self.tipoPisos.addItem(QtGui.QApplication.translate("equipment", "De borboteo", None, QtGui.QApplication.UnicodeUTF8))
        self.tipoPisos.addItem(QtGui.QApplication.translate("equipment", "De tamiz", None, QtGui.QApplication.UnicodeUTF8))
        self.tipoPisos.currentIndexChanged.connect(self.calcularCostos)
        gridLayout_1.addWidget(self.tipoPisos, 1, 2, 1, 1)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Material:", None, QtGui.QApplication.UnicodeUTF8)), 2, 1, 1, 1)
        self.materialPisos=QtGui.QComboBox()
        self.materialPisos.addItem(QtGui.QApplication.translate("equipment", "Acero al carbon", None, QtGui.QApplication.UnicodeUTF8))
        self.materialPisos.addItem(QtGui.QApplication.translate("equipment", "Acero inoxidable 304", None, QtGui.QApplication.UnicodeUTF8))
        self.materialPisos.addItem(QtGui.QApplication.translate("equipment", "Acero inoxidable 316", None, QtGui.QApplication.UnicodeUTF8))
        self.materialPisos.addItem(QtGui.QApplication.translate("equipment", "Carpenter 20CB-3", None, QtGui.QApplication.UnicodeUTF8))
        self.materialPisos.addItem(QtGui.QApplication.translate("equipment", "Niquel 200", None, QtGui.QApplication.UnicodeUTF8))
        self.materialPisos.addItem(QtGui.QApplication.translate("equipment", "Monel 400", None, QtGui.QApplication.UnicodeUTF8))
        self.materialPisos.addItem(QtGui.QApplication.translate("equipment", "Inconel 600", None, QtGui.QApplication.UnicodeUTF8))
        self.materialPisos.addItem(QtGui.QApplication.translate("equipment", "Incoloy 825", None, QtGui.QApplication.UnicodeUTF8))
        self.materialPisos.addItem(QtGui.QApplication.translate("equipment", "Titanio", None, QtGui.QApplication.UnicodeUTF8))
        self.materialPisos.currentIndexChanged.connect(self.calcularCostos)
        gridLayout_1.addWidget(self.materialPisos, 2, 2, 1, 1)
        gridLayout_1.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),3,1,1,2)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Diametro:", None, QtGui.QApplication.UnicodeUTF8)), 4, 1, 1, 1)
        self.diametroPisos=Entrada_con_unidades(unidades.Length)
        gridLayout_1.addWidget(self.diametroPisos,4,2,1,1)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Número:", None, QtGui.QApplication.UnicodeUTF8)), 5, 1, 1, 1)
        self.NumeroPisos=Entrada_con_unidades(int, spinbox=True, min=1, step=1, width=50)
        gridLayout_1.addWidget(self.NumeroPisos,5,2,1,1)
        gridLayout_1.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),6,1,1,2)

        self.groupBox_relleno = QtGui.QGroupBox(QtGui.QApplication.translate("equipment", "Torre de relleno", None, QtGui.QApplication.UnicodeUTF8))
        gridLayout_Costos.addWidget(self.groupBox_relleno,1, 4, 4, 2)
        gridLayout_2 = QtGui.QGridLayout(self.groupBox_relleno)
        gridLayout_2.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Volumen:", None, QtGui.QApplication.UnicodeUTF8)), 1, 1, 1, 1)
        self.VolumenRelleno=Entrada_con_unidades(unidades.Volume, "VolLiq")
        gridLayout_2.addWidget(self.VolumenRelleno,1,2,1,1)
        gridLayout_2.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Coste unitario:", None, QtGui.QApplication.UnicodeUTF8)),2,1,1,1)
        self.C_unit_relleno=Entrada_con_unidades(unidades.Currency, retornar=False, textounidad="%s / %s" % (unidades.Currency(None).text(), unidades.Volume(None).text("VolLiq")))
        gridLayout_2.addWidget(self.C_unit_relleno,2,2,1,1)
        gridLayout_2.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),3,1,1,2)

        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Diametro:", None, QtGui.QApplication.UnicodeUTF8)),5,1,1,1)
        self.Dc=Entrada_con_unidades(unidades.Length)
        gridLayout_Costos.addWidget(self.Dc,5,2,1,2)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Altura:", None, QtGui.QApplication.UnicodeUTF8)), 6, 1, 1, 1)
        self.Hc=Entrada_con_unidades(unidades.Length)
        gridLayout_Costos.addWidget(self.Hc,6,2,1,2)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Espesor (Tapa):", None, QtGui.QApplication.UnicodeUTF8)), 7, 1, 1, 1)
        self.EspesorSuperior=Entrada_con_unidades(unidades.Length, "Thickness")
        gridLayout_Costos.addWidget(self.EspesorSuperior,7,2,1,2)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Espesor (Fondo):", None, QtGui.QApplication.UnicodeUTF8)), 8, 1, 1, 1)
        self.EspesorInferior=Entrada_con_unidades(unidades.Length, "Thickness")
        gridLayout_Costos.addWidget(self.EspesorInferior,8,2,1,2)
        gridLayout_Costos.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Densidad:", None, QtGui.QApplication.UnicodeUTF8)), 9, 1, 1, 1)
        self.EspesorInferior=Entrada_con_unidades(unidades.Density, "DenLiq")
        gridLayout_Costos.addWidget(self.EspesorInferior,9,2,1,2)

        self.Costos=costIndex.CostData(3, 3)
        self.Costos.valueChanged.connect(self.calcularCostos)
        gridLayout_Costos.addWidget(self.Costos,10,1,2,4)

        gridLayout_Costos.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),12,1,1,6)

        self.groupBox_Costos = QtGui.QGroupBox(QtGui.QApplication.translate("equipment", "Costos calculados", None, QtGui.QApplication.UnicodeUTF8))
        gridLayout_Costos.addWidget(self.groupBox_Costos,13,1,1,5)
        gridLayout_5 = QtGui.QGridLayout(self.groupBox_Costos)
        gridLayout_5.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Coste Pisos:", None, QtGui.QApplication.UnicodeUTF8)),0,1,1,1)
        self.C_pisos=Entrada_con_unidades(unidades.Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_pisos,0,2,1,1)
        gridLayout_5.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Coste Carcasa:", None, QtGui.QApplication.UnicodeUTF8)),1,1,1,1)
        self.C_carcasa=Entrada_con_unidades(unidades.Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_carcasa,1,2,1,1)
        gridLayout_5.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Coste Accesorios:", None, QtGui.QApplication.UnicodeUTF8)),2,1,1,1)
        self.C_accesorios=Entrada_con_unidades(unidades.Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_accesorios,2,2,1,1)
        gridLayout_5.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Coste Columna:", None, QtGui.QApplication.UnicodeUTF8)),0,4,1,1)
        self.C_columna=Entrada_con_unidades(unidades.Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_columna,0,5,1,1)
        gridLayout_5.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Coste Adquisición:", None, QtGui.QApplication.UnicodeUTF8)),1,4,1,1)
        self.C_adq=Entrada_con_unidades(unidades.Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_adq,1,5,1,1)
        gridLayout_5.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Coste Instalación:", None, QtGui.QApplication.UnicodeUTF8)),2,4,1,1)
        self.C_inst=Entrada_con_unidades(unidades.Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_inst,2,5,1,1)


#        #Pestaña salida
#        self.pSalida = QtGui.QTabWidget()
#        self.tabWidget.addTab(self.pSalida,QtGui.QApplication.translate("equipment", "Salida", None, QtGui.QApplication.UnicodeUTF8))


        self.mostrarSubclasificacion(0)


    def mostrarSubclasificacion(self, ind):
        self.groupBox_Pisos.setVisible(not ind)
        self.groupBox_relleno.setVisible(ind)

    def cambiar_entrada(self, corriente):
        selfentrada=corriente
        self.calculo()

    def calculo(self):
        if self.todos_datos():

            self.rellenoSalida()

    def rellenoSalida(self):
        pass

    def todos_datos(self):
        pass

    def calcularCostos(self):
        if self.todos_datos():
            self.ShellTube.Coste(self.factorInstalacion.value(), 0, self.tipo.currentIndex(), self.material.currentIndex())
            self.C_adq.setValue(self.ShellTube.C_adq.config())
            self.C_inst.setValue(self.ShellTube.C_inst.config())

if __name__ == "__main__":
    import sys
    from lib.corriente import Corriente, Mezcla
    app = QtGui.QApplication(sys.argv)
    agua=Corriente(T=300, P=101325, caudalMasico=3600, ids=[62], fraccion=[1])
    dialogo = UI_equipment(agua)
    dialogo.show()
    sys.exit(app.exec_())
