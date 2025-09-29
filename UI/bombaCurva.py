#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2025, Juan José Gómez Romera <jjgomera@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


:class:`UI_bombaCurva`: Dialog to manipulate pump working curve
'''


import pickle

from numpy import transpose
from tools.qt import QtCore, QtWidgets

from lib.plot import PlotWidget
from lib.unidades import Length, VolFlow, Power, Dimensionless
from UI.inputTable import InputTable
from UI.widgets import Entrada_con_unidades


class Ui_bombaCurva(QtWidgets.QDialog):
    """Dialog to manipulate pump working curve"""
    def __init__(self, curva=None, parent=None):
        """curva: Optional parameter to set characteristic pump curve"""
        super().__init__(parent)
        self.setWindowTitle(self.tr("Pump curves dialog"))

        lyt = QtWidgets.QGridLayout(self)
        self.buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Apply
            | QtWidgets.QDialogButtonBox.StandardButton.Open
            | QtWidgets.QDialogButtonBox.StandardButton.Save
            | QtWidgets.QDialogButtonBox.StandardButton.Reset,
            QtCore.Qt.Orientation.Vertical)
        self.buttons.clicked.connect(self.buttonsClicked)
        lyt.addWidget(self.buttons, 1, 1, 5, 1)

        lyt.addWidget(QtWidgets.QLabel(self.tr("Curves")), 1, 3)
        self.lst = QtWidgets.QComboBox()
        self.lst.currentIndexChanged.connect(self.fillCurve)
        lyt.addWidget(self.lst, 1, 4)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Diameter")), 2, 3)
        self.diameter = Entrada_con_unidades(int, width=60, textounidad='"')
        lyt.addWidget(self.diameter, 2, 4)
        lyt.addWidget(QtWidgets.QLabel(self.tr("RPM")), 3, 3)
        self.rpm = Entrada_con_unidades(int, width=60, textounidad="rpm")
        lyt.addWidget(self.rpm, 3, 4)
        self.grid = QtWidgets.QCheckBox(self.tr("Grid"))
        self.grid.toggled.connect(self.rejilla_toggled)
        lyt.addWidget(self.grid, 4, 3, 1, 2)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Minimum), 5, 4)

        header = [self.tr("Flowrate"), self.tr("Head"),
                  self.tr("Power"), self.tr("NPSH")]
        self.dataTable = InputTable(
            4, horizontalHeader=header, stretch=False, verticalHeader=False,
            unit=[VolFlow, Length, Power, Dimensionless])
        self.dataTable.setColumnWidth(0, 120)
        self.dataTable.setColumnWidth(1, 100)
        self.dataTable.setColumnWidth(2, 100)
        self.dataTable.setColumnWidth(3, 100)
        self.dataTable.setFixedWidth(440)
        self.dataTable.setConnected()
        self.dataTable.setVerticalScrollBarPolicy(
            QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        lyt.addWidget(self.dataTable, 6, 1, 1, 4)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 1, 5)

        self.Plot = PlotWidget(width=5, height=1, dpi=100, dim=0, parent=self)
        lyt.addWidget(self.Plot, 1, 6, 6, 5)
        lyt.addItem(QtWidgets.QSpacerItem(
            10000, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 1, 6, 6, 5)

        buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok
            | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)
        lyt.addWidget(buttonBox, 7, 9, 1, 2)

        if curva:
            self.curvas = curva
            for i in curva:
                self.lst.addItem(str(i[0])+'", '+str(i[1])+" rpm")
            self.lst.setCurrentIndex(self.lst.count()-1)
            # self.actualizarPlot()
        else:
            self.curva = []

    def buttonsClicked(self, boton):
        """Do actions when clicked same of buttons"""
        if boton == self.buttons.button(
                QtWidgets.QDialogButtonBox.StandardButton.Reset):
            # Reset button clicked
            self.lst.setCurrentIndex(-1)
            self.diameter.clear()
            self.rpm.clear()
            self.dataTable.clear()
            self.checkHead.setChecked(False)
            self.checkPower.setChecked(False)
            self.checkNPSH.setChecked(False)
            self.grid.setChecked(False)

        elif boton == self.buttons.button(
                QtWidgets.QDialogButtonBox.StandardButton.Open):
            # Open button clicked
            fname = QtWidgets.QFileDialog.getOpenFileName(
                self, self.tr("Open curve file"), "./", "cpickle file (*.pkl)")
            if fname[0]:
                with open(fname[0], "r") as archivo:
                    curvas = pickle.load(archivo)
                self.curvas = curvas
                self.curva = curvas[-1]
                self.lst.clear()
                self.lst.blockSignals(True)
                for i in curvas:
                    self.lst.addItem(str(i[0])+'", '+str(i[1])+" rpm")
                self.lst.blockSignals(False)
                self.lst.setCurrentIndex(self.lst.count()-1)
                self.fillCurve(self.lst.count()-1)
                self.actualizarPlot()

        elif boton == self.buttons.button(
                QtWidgets.QDialogButtonBox.StandardButton.Apply):
            # Apply button clicked
            txt = f'{self.diameter.value}", {self.rpm.value} rpm'
            index = self.lst.findText(txt)
            print(txt)
            if index == -1:
                self.lst.addItem(txt)
                self.curva = self.getCurve()
                self.curvas.append(self.curva)
                self.lst.setCurrentIndex(self.lst.count()-1)
                self.actualizarPlot()
            else:
                self.curva = self.getCurve()
                self.curvas[index] = self.curva
                self.lst.setCurrentIndex(index)
                self.actualizarPlot()

        elif boton == self.buttons.button(
                QtWidgets.QDialogButtonBox.StandardButton.Save):
            # Save button clicked
            fname = QtWidgets.QFileDialog.getSaveFileName(
                self, self.tr("Save curve to file"), "./", "cpickle file (*.pkl)")
            if fname[0]:
                if fname[0].split(".")[-1] != "pkl":
                    fname[0]+=".pkl"
                with open(fname[0], "w") as file:
                    pickle.dump(self.curvas, file)

    def getCurve(self):
        """Get curve from widget"""
        flow = self.dataTable.getUnifiedColumn(0)
        head = self.dataTable.getUnifiedColumn(1)
        power = self.dataTable.getUnifiedColumn(2)
        npsh = self.dataTable.getUnifiedColumn(3)
        return [self.diameter.value, self.rpm.value, flow, head, power, npsh]

    def fillCurve(self, ind):
        """Show selected curve data in widget"""
        self.curva = self.curvas[ind]
        self.diameter.setValue(self.curva[0])
        self.rpm.setValue(self.curva[1])

        self.dataTable.setData(transpose(self.curva[2:]))
        self.actualizarPlot()

    def actualizarPlot(self):
        """Update plot"""
        self.Plot.fig.clear()
        ax = []
        ax.append(self.Plot.fig.add_subplot(3, 1, 1))
        ax.append(self.Plot.fig.add_subplot(3, 1, 2, sharex=ax[-1]))
        ax.append(self.Plot.fig.add_subplot(3, 1, 3, sharex=ax[-1]))
        self.Plot.ax = ax

        for curva in self.curvas:
            q = curva[2]
            h = curva[3]
            self.Plot.ax[0].plot(
                q, h, label=str(curva[0])+'", '+str(curva[1])+" rpm")
        self.Plot.ax[0].grid(self.grid.isChecked())
        self.Plot.ax[0].set_title("Curva H/Q", size='14')
        self.Plot.ax[0].set_ylabel("H, m", size='12')
        self.Plot.ax[0].legend()

        for curva in self.curvas:
            q = curva[2]
            w = curva[4]
            self.Plot.ax[1].plot(
                q, w, label=str(curva[0])+'", '+str(curva[1])+" rpm")
        self.Plot.ax[1].grid(self.grid.isChecked())
        self.Plot.ax[1].set_ylabel("P, kW", size='12')
        self.Plot.ax[1].legend(loc='upper left')

        for curva in self.curvas:
            q = curva[2]
            npsh = curva[5]
            self.Plot.ax[2].plot(
                q, npsh, label=str(curva[0])+'", '+str(curva[1])+" rpm")
        self.Plot.ax[2].grid(self.grid.isChecked())
        self.Plot.ax[2].set_ylabel("NPSH, m", size='12')
        self.Plot.ax[2].legend(loc='upper left')
        self.Plot.ax[2].set_xlabel(
            "Q, $m^3/s$", horizontalalignment='left', size='12')

        self.Plot.draw()

    def rejilla_toggled(self, boolean):
        """Show/hide grid in plot"""
        for ax in self.Plot.ax:
            ax.grid(boolean)
        self.Plot.draw()


if __name__ == "__main__":
    import sys
    from numpy import zeros, arange
    app = QtWidgets.QApplication(sys.argv)
    Q = arange(0, 21, 1.)
    H = [15.5, 15.468115, 15.40372, 15.319705, 15.23896, 15.154375, 15.06884,
         14.945245, 14.84648, 14.725435, 14.605, 14.418065, 14.18752,
         13.906255, 13.56716, 13.163125, 12.68704, 12.131795, 11.49028,
         10.755385, 9.92]
    Pot = [0.5, 0.51, 0.53, 0.55, 0.575, 0.595, 0.605, 0.625, 0.66, 0.685,
           0.705, 0.73, 0.765, 0.795, 0.82, 0.855, 0.895, 0.92, 0.975, 1.005,
           1.06]
    NHPS = zeros(21)
    curv = [[212, 1400, Q, H, [1e3*p for p in Pot], NHPS], ]

    bombaCurva = Ui_bombaCurva(curv)
    bombaCurva.show()
    sys.exit(app.exec())
