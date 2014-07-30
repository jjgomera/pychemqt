#!/usr/bin/python
# -*- coding: utf-8 -*-

import os

from PyQt4 import QtCore, QtGui
from numpy import loadtxt

from widgets import Tabla
from tools.HelpView import HelpView
from UI.widgets import Entrada_con_unidades
from lib.unidades import Temperature


class eqDIPPR(QtGui.QWidget):
    def __init__(self, value, parent=None):
        super(eqDIPPR, self).__init__(parent)
        layout=QtGui.QHBoxLayout(self)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Eq DIPPR")+" ", self))
        self.eqDIPPR = QtGui.QSpinBox(self)
        self.eqDIPPR.setValue(value)
        self.eqDIPPR.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.eqDIPPR.setFixedWidth(50)
        self.eqDIPPR.setRange(1, 9)
        txt=u"""
    1:     Y = A+B*T+C*T²+D*T³+E*T⁴
    2:     Y = exp(A+B*T+C*ln(T)+D*T^E)
    3:     Y = A*T^(B/(1+C*T+D*T^2))
    4:     Y = A+B*exp(-C/T^D)
    5:     Y = A + B/T + C/T³ + D/T⁸ + E/T⁹
    6:     Y = A/(B^(1+(1-T/C)^D)
    7:     Y = A*(1-Tr)^(B+C*Tr+D*Tr²+E*Tr³)
    8:     Y = A+ B*((C/T)/sinh(C/T))² + D*((E/T)/cosh(E/T))²
    9:     Y = A²/Tr + B - 2ACTr - ADTr² - C²Tr³/3 - CDTr⁴/2 - D²Tr⁵/5
        """
        var=QtGui.QApplication.translate("pychemqt", """where:
                    Y Property to fit
                    T temperature in Kelvin
                    Tr: reduced temperature T/Tc
                    A,B,C,D,E parameters""")
        self.eqDIPPR.setToolTip(QtGui.QApplication.translate("pychemqt", "Equation")+txt+var)
        layout.addWidget(self.eqDIPPR)
        layout.addStretch()

    @property
    def value(self):
        return self.eqDIPPR.value

    def setValue(self, value):
        self.eqDIPPR.setValue(value)

    def clear(self):
        self.eqDIPPR.clear()

class Entrada_Datos(QtGui.QDialog):
    """Dialogo de entrada de datos en forma de tabla"""
    def __init__(self, data=None, t=[], property=[], horizontalHeader=[], title="", help=False, helpFile="", DIPPR=False, tc=0, tcValue=None, eq=1, parent=None):
        """
        data=matriz con los datos de la tabla
        t: valores de la columna de temperatura
        property: valores de la columna de propiedades
        horizontalHeader: Lista con los títulos de las columnas
        title: titulo de la ventana
        help: boolean que indica si se activa el botón de ayuda
        helpFile: Ruta del archivo de ayuda, ya sea un archivo del disco duro o una url
        DIPPR: boolean que indica si se añade un campo de texto para el número de eq de DIPPR
        tc: boolean que indica si se muestra una entrada para la tc
        tcValue: valor de tc
        eq: valor inicial de la ecuación
        """
        super(Entrada_Datos, self).__init__(parent)
        self.setWindowTitle(title)
        self.columnas=len(horizontalHeader)
        self.horizontalHeader=horizontalHeader
        self.title=title
        self.helpFile=helpFile
        gridLayout = QtGui.QGridLayout(self)
        self.botonAbrir = QtGui.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/fileOpen.png")), QtGui.QApplication.translate("pychemqt", "Open"))
        self.botonAbrir.clicked.connect(self.Abrir)
        gridLayout.addWidget(self.botonAbrir,1,1)
        self.botonGuardar = QtGui.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/fileSave.png")), QtGui.QApplication.translate("pychemqt", "Save"))
        self.botonGuardar.clicked.connect(self.Guardar)
        gridLayout.addWidget(self.botonGuardar,1,2)
        self.botonDelete=QtGui.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/clear.png")), QtGui.QApplication.translate("pychemqt", "Clear"))
        self.botonDelete.clicked.connect(self.Borrar)
        gridLayout.addWidget(self.botonDelete,1,3)
        gridLayout.addItem(QtGui.QSpacerItem(0,0,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),1,4)

        self.tabla=Tabla(self.columnas, horizontalHeader=horizontalHeader, verticalHeader=False, stretch=False)
        self.tabla.setConnected()
        if data:
            self.tabla.setMatrix(data)
            self.tabla.addRow()
        elif t and property:
            self.tabla.setColumn(0, t)
            self.tabla.setColumn(1, property)
        gridLayout.addWidget(self.tabla,2,1,1,4)

        if DIPPR:
            self.eqDIPPR=eqDIPPR(eq)
            gridLayout.addWidget(self.eqDIPPR,3,1,1,4)
            self.eqDIPPR.eqDIPPR.valueChanged.connect(self.showTc)

        if tc:
            lyt=QtGui.QHBoxLayout()
            self.labelTc=QtGui.QLabel("Tc: ", self)
            lyt.addWidget(self.labelTc)
            self.tc=Entrada_con_unidades(Temperature, value=tcValue, parent=self)
            lyt.addWidget(self.tc)
            lyt.addItem(QtGui.QSpacerItem(0,0,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding))
            gridLayout.addItem(lyt,4,1,1,4)
            self.showTc(1)

        if help:
            botones=QtGui.QDialogButtonBox.Help|QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok
        else:
            botones=QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok
        self.boton = QtGui.QDialogButtonBox(botones)
        self.boton.accepted.connect(self.accept)
        self.boton.rejected.connect(self.reject)
        self.boton.helpRequested.connect(self.ayuda)
        gridLayout.addWidget(self.boton,5,1,1,4)


    def showTc(self, value):
        self.labelTc.setVisible(value in (7, 9))
        self.tc.setVisible(value in (7, 9))

    def Abrir(self):
        fname = unicode(QtGui.QFileDialog.getOpenFileName(self, QtGui.QApplication.translate("pychemqt", "Open text file"), "./"))
        if fname:
            data=loadtxt(fname)
            self.tabla.setMatrix(data)
            self.tabla.addRow()

    def Guardar(self):
        fname = unicode(QtGui.QFileDialog.getSaveFileName(self, QtGui.QApplication.translate("pychemqt", "Save data to file"), "./"))
        if fname:
            with open(fname, 'w') as file:
                file.write("#"+self.title+"\n")
                file.write("#")
                try:
                    for i in self.horizontalHeader:
                        file.write(i+"\t")
                except UnicodeEncodeError:
                    pass
                file.write("\n")
                data=self.data
                for fila in range(len(data)):
                    for columna in range(self.tabla.columnCount()):
                        file.write(str(data[fila][columna])+"\t")
                    file.write("\n")

    def Borrar(self):
        self.tabla.setRowCount(1)
        self.tabla.clearContents()

    def ayuda(self):
        Dialog = HelpView(self.windowTitle(), QtCore.QUrl(self.helpFile))
        Dialog.exec_()

    @property
    def data(self):
        return self.tabla.getMatrix()


if __name__ == "__main__":
    import sys, os
    app = QtGui.QApplication(sys.argv)
    titulo=u"Distribución de tamaño de sólidos"
    encabezado=[u"Diametro, μm", u"Fracción másica", "acumulado"]
    ui = Entrada_Datos(horizontalHeader=encabezado, title=titulo, help=True, DIPPR=True, helpFile=os.environ["pychemqt"]+"doc/doc/pychemqt.UI.entrada_datos.html")
    ui.show()
    sys.exit(app.exec_())
