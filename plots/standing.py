#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2017, Juan José Gómez Romera <jjgomera@gmail.com>

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
# Standing Katz plot
###############################################################################


from PyQt5 import QtWidgets
from scipy import arange

from UI.widgets import Entrada_con_unidades
from lib.petro import Z_list
from lib.plot import mpl


class Standing_Katz(QtWidgets.QDialog):
    title = QtWidgets.QApplication.translate("pychemqt", "Standing Katz chart")

    def __init__(self, parent=None):
        super(Standing_Katz, self).__init__(parent)
        self.setWindowTitle(self.title)
        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Method:")),1,1)
        self.metodos = QtWidgets.QComboBox()
        self.metodos.addItem("Hall Yarborough")
        self.metodos.addItem("Dranchuk Abu-Kassem")
        self.metodos.addItem("Dranchuk Purvis Robinson")
        self.metodos.addItem("Beggs Brill")
        self.metodos.addItem("Sarem")
        self.metodos.addItem("Gopal")
        self.metodos.addItem("Papay")
        self.metodos.currentIndexChanged.connect(self.plot_Z)
        layout.addWidget(self.metodos,1,2)
        layout.addItem(QtWidgets.QSpacerItem(10,10,QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Fixed),1, 3)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Pr<sub>min</sub>")),1,4)
        self.Prmin = Entrada_con_unidades(float, spinbox=True, value=0.0, width=60, decimales=1, step=0.1)
        layout.addWidget(self.Prmin,1,5)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Pr<sub>max</sub>")),1,6)
        self.Prmax = Entrada_con_unidades(float, spinbox=True, value=8.0, width=60, decimales=1, step=0.1)
        layout.addWidget(self.Prmax,1,7)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "Tr")),1, 8)
        self.Tr = QtWidgets.QLineEdit(", ".join(str(i) for i in [1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.2, 2.4, 2.6, 2.8, 3.]))
        self.Tr.setMinimumWidth(400)
        layout.addWidget(self.Tr, 1, 9)
        self.diagrama = mpl(self, dpi=90)
        layout.addWidget(self.diagrama, 2, 1, 1, 9)

        self.buttonBox = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Close)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox, 5, 1, 1, 6)

        self.plot_Z(0)

    def plot_Z(self, indice):
        Prmin = self.Prmin.value
        Prmax = self.Prmax.value
        self.diagrama.config(self.Prmin.value, self.Prmax.value)

        try:
            Tr = self.Tr.text().split()
        except:
            Tr = [1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.6,
                  1.7, 1.8, 1.9, 2., 2.2, 2.4, 2.6, 2.8, 3.]
        Z = Z_list[indice]

        P = arange(Prmin, Prmax, 0.1)
        for Tr in Tr:
            self.diagrama.plot(P, [Z(Tr, Pr) for Pr in P], "k")
        title = QtWidgets.QApplication.translate("pychemqt", "Standing and Katz compressivitity factors chart for natural gas")
        self.diagrama.axes2D.set_title(title, size='12')
        self.diagrama.draw()

