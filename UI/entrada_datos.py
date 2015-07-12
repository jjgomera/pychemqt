#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Dialog to common data entry
#   -Entrada_Datos: Dialog to data entry
#   -eqDIPPR: Widget for select DIPPR equation
###############################################################################

import os

from PyQt5 import QtCore, QtGui, QtWidgets

from numpy import loadtxt

from .widgets import Tabla
from tools.HelpView import HelpView
from UI.widgets import Entrada_con_unidades
from lib.unidades import Temperature


class eqDIPPR(QtWidgets.QWidget):
    """Custom widget to define DIPPR equation input"""
    def __init__(self, value, parent=None):
        super(eqDIPPR, self).__init__(parent)
        layout = QtWidgets.QHBoxLayout(self)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Eq DIPPR") + " "))
        self.eqDIPPR = QtWidgets.QSpinBox()
        self.eqDIPPR.setValue(value)
        self.eqDIPPR.setRange(1, 9)
        self.eqDIPPR.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.eqDIPPR.setFixedWidth(50)
        txt = """
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
        var = QtWidgets.QApplication.translate("pychemqt", """where:
                    Y Property to fit
                    T temperature in Kelvin
                    Tr: reduced temperature T/Tc
                    A,B,C,D,E parameters""")
        self.eqDIPPR.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Equation") + txt + var)
        layout.addWidget(self.eqDIPPR)
        layout.addStretch()

    @property
    def value(self):
        return self.eqDIPPR.value

    def setValue(self, value):
        self.eqDIPPR.setValue(value)

    def clear(self):
        self.eqDIPPR.clear()


class Entrada_Datos(QtWidgets.QDialog):
    """Table data input dialog"""
    def __init__(self, data=None, t=[], property=[], horizontalHeader=[],
                 title="", help=False, helpFile="", DIPPR=False, tc=0,
                 tcValue=None, eq=1, parent=None):
        """
        title: window title
        data: mrray with original data
        t: values for x column, generally temperature
        property: values for 2...n columns
        horizontalHeader: List with column title
        help: boolean to show help button
        helpFile: Path for help file, file or url
        DIPPR: boolean to show DIPPR widget
        tc: boolean to show critical temperature (same DIPPR eq need it)
        tcValue: value for critical temperature
        eq: Value for DIPPR equation
        """
        super(Entrada_Datos, self).__init__(parent)
        self.setWindowTitle(title)
        self.columnas = len(horizontalHeader)
        self.horizontalHeader = horizontalHeader
        self.title = title
        self.helpFile = helpFile
        gridLayout = QtWidgets.QGridLayout(self)
        self.botonAbrir = QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"]+"/images/button/fileOpen.png")),
            QtWidgets.QApplication.translate("pychemqt", "Open"))
        self.botonAbrir.clicked.connect(self.Abrir)
        gridLayout.addWidget(self.botonAbrir, 1, 1)
        self.botonGuardar = QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"]+"/images/button/fileSave.png")),
            QtWidgets.QApplication.translate("pychemqt", "Save"))
        self.botonGuardar.clicked.connect(self.Guardar)
        gridLayout.addWidget(self.botonGuardar, 1, 2)
        self.botonDelete = QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"]+"/images/button/clear.png")),
            QtWidgets.QApplication.translate("pychemqt", "Clear"))
        self.botonDelete.clicked.connect(self.Borrar)
        gridLayout.addWidget(self.botonDelete, 1, 3)
        gridLayout.addItem(QtWidgets.QSpacerItem(0, 0, QtWidgets.QSizePolicy.Expanding,
                                             QtWidgets.QSizePolicy.Expanding), 1, 4)

        self.tabla = Tabla(self.columnas, horizontalHeader=horizontalHeader,
                           verticalHeader=False, stretch=False)
        self.tabla.setConnected()
        if data:
            self.tabla.setMatrix(data)
            self.tabla.addRow()
        elif t and property:
            self.tabla.setColumn(0, t)
            self.tabla.setColumn(1, property)
        gridLayout.addWidget(self.tabla, 2, 1, 1, 4)

        if DIPPR:
            self.eqDIPPR = eqDIPPR(eq)
            gridLayout.addWidget(self.eqDIPPR, 3, 1, 1, 4)
            self.eqDIPPR.eqDIPPR.valueChanged.connect(self.showTc)

        if tc:
            lyt = QtWidgets.QHBoxLayout()
            self.labelTc = QtWidgets.QLabel("Tc: ", self)
            lyt.addWidget(self.labelTc)
            self.tc = Entrada_con_unidades(Temperature, value=tcValue)
            lyt.addWidget(self.tc)
            lyt.addItem(QtWidgets.QSpacerItem(0, 0, QtWidgets.QSizePolicy.Expanding,
                                          QtWidgets.QSizePolicy.Expanding))
            gridLayout.addItem(lyt, 4, 1, 1, 4)
            self.showTc(1)

        if help:
            botones = QtWidgets.QDialogButtonBox.Help | \
                QtWidgets.QDialogButtonBox.Cancel | QtWidgets.QDialogButtonBox.Ok
        else:
            botones = QtWidgets.QDialogButtonBox.Cancel | QtWidgets.QDialogButtonBox.Ok
        self.boton = QtWidgets.QDialogButtonBox(botones)
        self.boton.accepted.connect(self.accept)
        self.boton.rejected.connect(self.reject)
        self.boton.helpRequested.connect(self.ayuda)
        gridLayout.addWidget(self.boton, 5, 1, 1, 4)

    def showTc(self, value):
        self.labelTc.setVisible(value in (7, 9))
        self.tc.setVisible(value in (7, 9))

    def Abrir(self):
        fname = str(QtWidgets.QFileDialog.getOpenFileName(
            self, QtWidgets[0].QCoreApplication.translate("pychemqt", "Open text file"), "./"))
        if fname:
            data = loadtxt(fname)
            self.tabla.setMatrix(data)
            self.tabla.addRow()

    def Guardar(self):
        fname = str(QtWidgets.QFileDialog.getSaveFileName(
            self, QtWidgets[0].QCoreApplication.translate("pychemqt", "Save data to file"), "./"))
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
                data = self.data
                for fila in range(len(data)):
                    for columna in range(self.tabla.columnCount()):
                        file.write(str(data[fila][columna])+"\t")
                    file.write("\n")

    def Borrar(self):
        """Clear table"""
        self.tabla.setRowCount(1)
        self.tabla.clearContents()

    def ayuda(self):
        """Show help file"""
        Dialog = HelpView(self.windowTitle(), QtCore.QUrl(self.helpFile))
        Dialog.exec_()

    @property
    def data(self):
        return self.tabla.getMatrix()


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    titulo = "Distribución de tamaño de sólidos"
    encabezado = ["Diametro, μm", "Fracción másica", "acumulado"]
    ui = Entrada_Datos(horizontalHeader=encabezado, title=titulo, help=True,
                       DIPPR=True, helpFile=os.environ["pychemqt"] +
                       "doc/doc/pychemqt.UI.entrada_datos.html")
    ui.show()
    sys.exit(app.exec_())
