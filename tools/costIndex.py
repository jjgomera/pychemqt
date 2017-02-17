#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2016, Juan José Gómez Romera <jjgomera@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


###############################################################################
# Library with costindex functionality
#   -UI_CostIndex: Dialog to configure costindex values
#   -CostData: Common widget to equipment with cost section
###############################################################################


from functools import partial
import os

from PyQt5 import QtCore, QtWidgets

from UI.widgets import Entrada_con_unidades
from lib import config

indiceBase = ["Jan-1982", 313.95, 336.19, 326.01, 312.03, 383.18, 297.63,
              421.1, 235.42, 338.2, 263.92, 290.13, 303.26]

# Load data from config file
indiceActual = []
with open(config.conf_dir+"CostIndex.dat", "r") as archivo:
    indiceActual.append(archivo.readline()[:-1])
    while True:
        data = archivo.readline().rstrip("\n")
        if data:
            indiceActual.append(float(data))
            if len(indiceActual) == 13:
                break


class Ui_CostIndex(QtWidgets.QDialog):
    """Dialog to show/configure costIndex"""
    def __init__(self, parent=None):
        super(Ui_CostIndex, self).__init__(parent)
        self.setWindowTitle(
            QtWidgets.QApplication.translate("pychemqt", "Cost Index"))
        self.custom = True
        layout = QtWidgets.QGridLayout(self)
        self.fecha = QtWidgets.QComboBox()
        layout.addWidget(self.fecha, 1, 1, 1, 3)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "CE INDEX")), 2, 1, 1, 2)
        self.index = Entrada_con_unidades(float, width=70, decimales=1)
        layout.addWidget(self.index, 2, 3, 1, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Equipments")), 3, 1, 1, 2)
        self.equipos = Entrada_con_unidades(float, width=70, decimales=1)
        layout.addWidget(self.equipos, 3, 3, 1, 1)
        layout.addItem(QtWidgets.QSpacerItem(
            30, 0, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            4, 1, 1, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Heat exchangers & Tanks")), 4, 2, 1, 1)
        self.cambiadores_calor = Entrada_con_unidades(
            float, width=70, decimales=1)
        layout.addWidget(self.cambiadores_calor, 4, 3, 1, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Process machinery")), 5, 2, 1, 1)
        self.maquinaria = Entrada_con_unidades(float, width=70, decimales=1)
        layout.addWidget(self.maquinaria, 5, 3, 1, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Pipe, valves & fittings")), 6, 2, 1, 1)
        self.tuberias = Entrada_con_unidades(float, width=70, decimales=1)
        layout.addWidget(self.tuberias, 6, 3, 1, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Process instruments")), 7, 2, 1, 1)
        self.instrumentos = Entrada_con_unidades(float, width=70, decimales=1)
        layout.addWidget(self.instrumentos, 7, 3, 1, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Pumps & compressors")), 8, 2, 1, 1)
        self.bombas = Entrada_con_unidades(float, width=70, decimales=1)
        layout.addWidget(self.bombas, 8, 3, 1, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Electrical equipments")), 9, 2, 1, 1)
        self.equipos_electricos = Entrada_con_unidades(
            float, width=70, decimales=1)
        layout.addWidget(self.equipos_electricos, 9, 3, 1, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Structural supports & misc")), 10, 2, 1, 1)
        self.soportes = Entrada_con_unidades(float, width=70, decimales=1)
        layout.addWidget(self.soportes, 10, 3, 1, 1)

        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Construction labor")), 11, 1, 1, 2)
        self.construccion = Entrada_con_unidades(float, width=70, decimales=1)
        layout.addWidget(self.construccion, 11, 3, 1, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Buildings")), 12, 1, 1, 2)
        self.edificios = Entrada_con_unidades(float, width=70, decimales=1)
        layout.addWidget(self.edificios, 12, 3, 1, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Engineering & supervision")), 13, 1, 1, 2)
        self.ingenieria = Entrada_con_unidades(float, width=70, decimales=1)
        layout.addWidget(self.ingenieria, 13, 3, 1, 1)
        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox, 14, 1, 1, 3)

        # Fill all widget data
        self.indices = []
        with open(os.environ["pychemqt"] + "dat/costindex.dat") as archivo:
            texto = archivo.readlines()
            for txt in texto:
                dato = txt.split()
                self.fecha.addItem(dato[0])
                self.indices.append(dato[1:])

        fecha = indiceActual[0]
        self.index.setValue(indiceActual[1])
        self.equipos.setValue(indiceActual[2])
        self.cambiadores_calor.setValue(indiceActual[3])
        self.maquinaria.setValue(indiceActual[4])
        self.tuberias.setValue(indiceActual[5])
        self.instrumentos.setValue(indiceActual[6])
        self.bombas.setValue(indiceActual[7])
        self.equipos_electricos.setValue(indiceActual[8])
        self.soportes.setValue(indiceActual[9])
        self.construccion.setValue(indiceActual[10])
        self.edificios.setValue(indiceActual[11])
        self.ingenieria.setValue(indiceActual[12])
        if fecha:
            self.fecha.setCurrentIndex(self.fecha.findText(fecha))
        else:
            self.fecha.setCurrentIndex(-1)

        self.fecha.currentIndexChanged.connect(self.loadData)

    def setCustom(self):
        """Set custom currentIndex"""
        self.fecha.setCurrentIndex(-1)
        self.custom = True

    def loadData(self, int):
        """Load costIndex data from file"""
        self.index.setValue(float(self.indices[int][0]))
        self.equipos.setValue(float(self.indices[int][1]))
        self.cambiadores_calor.setValue(float(self.indices[int][2]))
        self.maquinaria.setValue(float(self.indices[int][3]))
        self.tuberias.setValue(float(self.indices[int][4]))
        self.instrumentos.setValue(float(self.indices[int][5]))
        self.bombas.setValue(float(self.indices[int][6]))
        self.equipos_electricos.setValue(float(self.indices[int][7]))
        self.soportes.setValue(float(self.indices[int][8]))
        self.construccion.setValue(float(self.indices[int][9]))
        self.edificios.setValue(float(self.indices[int][10]))
        self.ingenieria.setValue(float(self.indices[int][11]))
        self.custom = False

    def closeEvent(self, event):
        """Override close event to ask data changes"""
        dialog = QtWidgets.QMessageBox.question(
            self,
            QtWidgets.QApplication.translate("pychemqt", "Unsaved changes"),
            QtWidgets.QApplication.translate("pychemqt",
                                             "Save unsaved changes?"),
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
            QtWidgets.QMessageBox.Yes)
        if dialog == QtWidgets.QMessageBox.Yes:
            self.accept()
        else:
            event.accept()

    def accept(self):
        """Overwrite accept signal to save changes"""
        with open(config.conf_dir+"CostIndex.dat", "w") as archivo:
            if self.custom:
                archivo.write("custom\n")
            else:
                archivo.write(self.fecha.currentText()+"\n")

            archivo.write(str(self.index.value)+"\n")
            archivo.write(str(self.equipos.value)+"\n")
            archivo.write(str(self.cambiadores_calor.value)+"\n")
            archivo.write(str(self.maquinaria.value)+"\n")
            archivo.write(str(self.tuberias.value)+"\n")
            archivo.write(str(self.instrumentos.value)+"\n")
            archivo.write(str(self.bombas.value)+"\n")
            archivo.write(str(self.equipos_electricos.value)+"\n")
            archivo.write(str(self.soportes.value)+"\n")
            archivo.write(str(self.construccion.value)+"\n")
            archivo.write(str(self.edificios.value)+"\n")
            archivo.write(str(self.ingenieria.value)+"\n")
            QtWidgets.QDialog.accept(self)


class CostData(QtWidgets.QWidget):
    """Common widget to equipment with cost section
    It have property to easy access to properties:
        factor: install factor
        base: base index (January 1982)
        actual: current index
        values: a tuple with all properties, (factor, base,actual)
    """
    valueChanged = QtCore.pyqtSignal(str, float)

    def __init__(self, equipment, parent=None):
        """constructor
        equipment: equipment class where the widget have to be put, define
        indiceCostos as a index in costIndex"""
        super(CostData, self).__init__(parent)
        self.indice = equipment.indiceCostos
        factor = equipment.kwargs["f_install"]
        gridLayout = QtWidgets.QGridLayout(self)
        gridLayout.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 1, 0, 1, 7)
        gridLayout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Instalation factor:")), 2, 0, 1, 1)
        self.factorInstalacion = Entrada_con_unidades(
            float, spinbox=True, decimales=1, step=0.1, width=50, value=factor)
        self.factorInstalacion.valueChanged.connect(partial(
            self.valueChanged.emit, "f_install"))
        gridLayout.addWidget(self.factorInstalacion, 2, 1, 1, 1)
        gridLayout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Base index:")), 2, 4, 1, 1)
        self.indiceBase = Entrada_con_unidades(
            float, readOnly=True, value=indiceBase[self.indice], decimales=1)
        gridLayout.addWidget(self.indiceBase, 2, 5, 1, 1)
        gridLayout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Current index:")), 3, 4, 1, 1)
        self.indiceActual = Entrada_con_unidades(
            float, readOnly=True, colorReadOnly="white",
            value=indiceActual[self.indice], decimales=1)
        gridLayout.addWidget(self.indiceActual, 3, 5, 1, 1)
        self.costIndex = QtWidgets.QToolButton()
        self.costIndex.setFixedSize(QtCore.QSize(24, 24))
        self.costIndex.clicked.connect(self.on_costIndex_clicked)
        self.costIndex.setText("...")
        self.costIndex.setVisible(False)
        gridLayout.addWidget(self.costIndex, 3, 5, 1, 1)
        gridLayout.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 4, 0, 1, 7)

    def on_costIndex_clicked(self):
        """Show costIndes dialog to show/change"""
        dialog = Ui_CostIndex()
        if dialog.exec_():
            with open(config.conf_dir+"CostIndex.dat", "r") as archivo:
                self.indiceActual.setValue(
                    float(archivo.readlines()[self.indice][:-1]))
            self.valueChanged.emit("Current_index", self.indiceActual.value)

    def enterEvent(self, event):
        self.costIndex.setVisible(True)

    def leaveEvent(self, event):
        self.costIndex.setVisible(False)

    @property
    def factor(self):
        return self.factorInstalacion.value

    @property
    def base(self):
        return self.indiceBase.value

    @property
    def actual(self):
        return self.indiceActual.value

    @property
    def values(self):
        return self.factorInstalacion.value, self.indiceBase.value,\
            self.indiceActual.value

    def setFactor(self, value):
        self.factorInstalacion.setValue(value)

    def setBase(self, value):
        self.indiceBase.setValue(value)

    def setActual(self, value):
        self.indiceActual.setValue(value)

    def setValues(self, factor, base, actual):
        self.factorInstalacion.setValue(factor)
        self.indiceBase.setValue(base)
        self.indiceActual.setValue(actual)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    CostIndex = Ui_CostIndex()
    CostIndex.show()
    sys.exit(app.exec_())
