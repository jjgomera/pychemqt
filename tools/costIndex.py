#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2023, Juan José Gómez Romera <jjgomera@gmail.com>

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

from tools.qt import QtCore, QtWidgets, tr

from UI.widgets import Entrada_con_unidades
from lib.config import conf_dir

indiceBase = ["Jan-1982", 313.95, 336.19, 326.01, 312.03, 383.18, 297.63,
              421.1, 235.42, 338.2, 263.92, 290.13, 303.26]

# Load data from config file
indiceActual = []
with open(os.path.join(conf_dir, "CostIndex.dat"), "r") as costFile:
    indiceActual.append(costFile.readline()[:-1])
    while True:
        data = costFile.readline().rstrip("\n")
        if data:
            indiceActual.append(float(data))
            if len(indiceActual) == 13:
                break


class Ui_CostIndex(QtWidgets.QDialog):
    """Dialog to show/configure costIndex"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(
            tr("pychemqt", "Cost Index"))
        self.custom = True
        layout = QtWidgets.QGridLayout(self)
        self.fecha = QtWidgets.QComboBox()
        layout.addWidget(self.fecha, 1, 1, 1, 3)
        layout.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "CE INDEX")), 2, 1, 1, 2)
        self.index = Entrada_con_unidades(float, width=70, decimales=1)
        layout.addWidget(self.index, 2, 3, 1, 1)
        layout.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Equipments")), 3, 1, 1, 2)
        self.equipos = Entrada_con_unidades(float, width=70, decimales=1)
        layout.addWidget(self.equipos, 3, 3, 1, 1)
        layout.addItem(QtWidgets.QSpacerItem(
            30, 0, QtWidgets.QSizePolicy.Policy.Fixed, QtWidgets.QSizePolicy.Policy.Fixed),
            4, 1, 1, 1)
        layout.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Heat exchangers & Tanks")), 4, 2, 1, 1)
        self.cambiadores_calor = Entrada_con_unidades(
            float, width=70, decimales=1)
        layout.addWidget(self.cambiadores_calor, 4, 3, 1, 1)
        layout.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Process machinery")), 5, 2, 1, 1)
        self.maquinaria = Entrada_con_unidades(float, width=70, decimales=1)
        layout.addWidget(self.maquinaria, 5, 3, 1, 1)
        layout.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Pipe, valves & fittings")), 6, 2, 1, 1)
        self.tuberias = Entrada_con_unidades(float, width=70, decimales=1)
        layout.addWidget(self.tuberias, 6, 3, 1, 1)
        layout.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Process instruments")), 7, 2, 1, 1)
        self.instrumentos = Entrada_con_unidades(float, width=70, decimales=1)
        layout.addWidget(self.instrumentos, 7, 3, 1, 1)
        layout.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Pumps & compressors")), 8, 2, 1, 1)
        self.bombas = Entrada_con_unidades(float, width=70, decimales=1)
        layout.addWidget(self.bombas, 8, 3, 1, 1)
        layout.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Electrical equipments")), 9, 2, 1, 1)
        self.equipos_electricos = Entrada_con_unidades(
            float, width=70, decimales=1)
        layout.addWidget(self.equipos_electricos, 9, 3, 1, 1)
        layout.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Structural supports & misc")), 10, 2, 1, 1)
        self.soportes = Entrada_con_unidades(float, width=70, decimales=1)
        layout.addWidget(self.soportes, 10, 3, 1, 1)

        layout.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Construction labor")), 11, 1, 1, 2)
        self.construccion = Entrada_con_unidades(float, width=70, decimales=1)
        layout.addWidget(self.construccion, 11, 3, 1, 1)
        layout.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Buildings")), 12, 1, 1, 2)
        self.edificios = Entrada_con_unidades(float, width=70, decimales=1)
        layout.addWidget(self.edificios, 12, 3, 1, 1)
        layout.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Engineering & supervision")), 13, 1, 1, 2)
        self.ingenieria = Entrada_con_unidades(float, width=70, decimales=1)
        layout.addWidget(self.ingenieria, 13, 3, 1, 1)
        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox, 14, 1, 1, 3)

        # Fill all widget data
        self.indices = []
        fdata = os.path.join(os.environ["pychemqt"], "dat", "costindex.dat")
        with open(fdata) as archivo:
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

    def loadData(self, i):
        """Load costIndex data from file"""
        self.index.setValue(float(self.indices[i][0]))
        self.equipos.setValue(float(self.indices[i][1]))
        self.cambiadores_calor.setValue(float(self.indices[i][2]))
        self.maquinaria.setValue(float(self.indices[i][3]))
        self.tuberias.setValue(float(self.indices[i][4]))
        self.instrumentos.setValue(float(self.indices[i][5]))
        self.bombas.setValue(float(self.indices[i][6]))
        self.equipos_electricos.setValue(float(self.indices[i][7]))
        self.soportes.setValue(float(self.indices[i][8]))
        self.construccion.setValue(float(self.indices[i][9]))
        self.edificios.setValue(float(self.indices[i][10]))
        self.ingenieria.setValue(float(self.indices[i][11]))
        self.custom = False

    def closeEvent(self, event):
        """Override close event to ask data changes"""
        dialog = QtWidgets.QMessageBox.question(
            self,
            tr("pychemqt", "Unsaved changes"),
            tr("pychemqt",
                                             "Save unsaved changes?"),
            QtWidgets.QMessageBox.StandardButton.Yes | QtWidgets.QMessageBox.StandardButton.No,
            QtWidgets.QMessageBox.StandardButton.Yes)
        if dialog == QtWidgets.QMessageBox.StandardButton.Yes:
            self.accept()
        else:
            event.accept()

    def accept(self):
        """Overwrite accept signal to save changes"""
        with open(os.path.join(conf_dir, "CostIndex.dat"), "w") as file:
            if self.custom:
                print()
                print("custom", file=file)
            else:
                print(self.fecha.currentText(), file=file)

            print(self.index.value, file=file)
            print(self.equipos.value, file=file)
            print(self.cambiadores_calor.value, file=file)
            print(self.maquinaria.value, file=file)
            print(self.tuberias.value, file=file)
            print(self.instrumentos.value, file=file)
            print(self.bombas.value, file=file)
            print(self.equipos_electricos.value, file=file)
            print(self.soportes.value, file=file)
            print(self.construccion.value, file=file)
            print(self.edificios.value, file=file)
            print(self.ingenieria.value, file=file)
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
        super().__init__(parent)
        self.indice = equipment.indiceCostos
        factor = equipment.kwargs["f_install"]
        gridLayout = QtWidgets.QGridLayout(self)
        gridLayout.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 1, 0, 1, 7)
        gridLayout.addWidget(QtWidgets.QLabel(
            tr(
                "pychemqt", "Instalation factor:")), 2, 0)
        self.factorInstalacion = Entrada_con_unidades(
            float, spinbox=True, decimales=1, step=0.1, width=50, value=factor)
        self.factorInstalacion.valueChanged.connect(partial(
            self.valueChanged.emit, "f_install"))
        gridLayout.addWidget(self.factorInstalacion, 2, 1)
        gridLayout.addWidget(QtWidgets.QLabel(
            tr("pychemqt", "Base index:")), 2, 4)
        self.indiceBase = Entrada_con_unidades(
            float, readOnly=True, value=indiceBase[self.indice], decimales=1)
        gridLayout.addWidget(self.indiceBase, 2, 5)
        gridLayout.addWidget(QtWidgets.QLabel(
            tr(
                "pychemqt", "Current index:")), 3, 4)
        self.indiceActual = Entrada_con_unidades(
            float, readOnly=True, colorReadOnly="white",
            value=indiceActual[self.indice], decimales=1)
        gridLayout.addWidget(self.indiceActual, 3, 5)
        self.costIndex = QtWidgets.QToolButton()
        self.costIndex.setFixedSize(QtCore.QSize(24, 24))
        self.costIndex.clicked.connect(self.on_costIndex_clicked)
        self.costIndex.setText("...")
        self.costIndex.setVisible(False)
        gridLayout.addWidget(self.costIndex, 3, 5)
        gridLayout.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 4, 0, 1, 7)

    def on_costIndex_clicked(self):
        """Show costIndes dialog to show/change"""
        dialog = Ui_CostIndex()
        if dialog.exec():
            with open(os.path.join(conf_dir, "CostIndex.dat"), "r") as file:
                self.indiceActual.setValue(
                    float(file.readlines()[self.indice][:-1]))
            self.valueChanged.emit("Current_index", self.indiceActual.value)

    def enterEvent(self, event=None):
        """Show costIndex button when mouse enter in widget"""
        self.costIndex.setVisible(True)

    def leaveEvent(self, event=None):
        """Hide costIndex button when mouse leave widget"""
        self.costIndex.setVisible(False)

    @property
    def factor(self):
        """Instalation factor value property"""
        return self.factorInstalacion.value

    def setFactor(self, value):
        """Set instalation factor value property"""
        self.factorInstalacion.setValue(value)

    @property
    def base(self):
        """Base index value"""
        return self.indiceBase.value

    def setBase(self, value):
        """Set base index value"""
        self.indiceBase.setValue(value)

    @property
    def actual(self):
        """Actual index value"""
        return self.indiceActual.value

    def setActual(self, value):
        """Set actual index value"""
        self.indiceActual.setValue(value)

    @property
    def values(self):
        """Tuple with the instalation factor, the base and actual index"""
        return self.factorInstalacion.value, self.indiceBase.value,\
            self.indiceActual.value

    def setValues(self, factor, base, actual):
        """Set values property"""
        self.factorInstalacion.setValue(factor)
        self.indiceBase.setValue(base)
        self.indiceActual.setValue(actual)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    CostIndex = Ui_CostIndex()
    CostIndex.show()
    sys.exit(app.exec())
