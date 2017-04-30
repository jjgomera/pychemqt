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
# Module for UI binary interaction parameter editor/viewer
###############################################################################

import os

from PyQt5 import QtWidgets

from lib.sql import getElement
from UI.widgets import Entrada_con_unidades


class Ui_BIP(QtWidgets.QDialog):
    """Dialog to view/edit the BIP for EoS available"""
    def __init__(self, i, j, parent=None):
        """Constructor
        i, j are the index for component to BIP show
        """
        super(Ui_BIP, self).__init__(parent)
        self.setWindowTitle(QtWidgets.QApplication.translate(
            "pychemqt", "BIP (Binary interaction parameters)"))

        lyt = QtWidgets.QGridLayout(self)
        title = QtWidgets.QApplication.translate("pychemqt", "Component")
        lyt.addWidget(QtWidgets.QLabel(title + " 1:"), 0, 0, 1, 1)
        lyt.addWidget(QtWidgets.QLabel(title + " 2:"), 1, 0, 1, 1)
        self.id1 = QtWidgets.QLabel()
        lyt.addWidget(self.id1, 0, 1, 1, 1)
        self.id2 = QtWidgets.QLabel()
        lyt.addWidget(self.id2, 1, 1, 1, 1)
        tabWidget = QtWidgets.QTabWidget()
        lyt.addWidget(tabWidget, 2, 0, 1, 3)
        button = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Cancel |
                                            QtWidgets.QDialogButtonBox.Ok)
        button.accepted.connect(self.accept)
        button.rejected.connect(self.reject)
        lyt.addWidget(button, 3, 2, 1, 1)

        self.nrtl = QtWidgets.QWidget()
        lyt_nrtl = QtWidgets.QGridLayout(self.nrtl)
        self.nrtl_widget = [
            Entrada_con_unidades(float), Entrada_con_unidades(float),
            Entrada_con_unidades(float), Entrada_con_unidades(float),
            Entrada_con_unidades(float), Entrada_con_unidades(float)]
        lyt_nrtl.addWidget(QtWidgets.QLabel("Bij"), 2, 1, 1, 1)
        lyt_nrtl.addWidget(self.nrtl_widget[0], 2, 2, 1, 1)
        lyt_nrtl.addWidget(QtWidgets.QLabel("Bji"), 2, 4, 1, 1)
        lyt_nrtl.addWidget(self.nrtl_widget[1], 2, 5, 1, 1)
        lyt_nrtl.addWidget(QtWidgets.QLabel(u"αij"), 3, 1, 1, 1)
        lyt_nrtl.addWidget(self.nrtl_widget[2], 3, 2, 1, 1)
        lyt_nrtl.addWidget(QtWidgets.QLabel(u"αji"), 3, 4, 1, 1)
        lyt_nrtl.addWidget(self.nrtl_widget[3], 3, 5, 1, 1)
        lyt_nrtl.addWidget(QtWidgets.QLabel("Aij"), 4, 1, 1, 1)
        lyt_nrtl.addWidget(self.nrtl_widget[4], 4, 2, 1, 1)
        lyt_nrtl.addWidget(QtWidgets.QLabel("Aji"), 4, 4, 1, 1)
        lyt_nrtl.addWidget(self.nrtl_widget[5], 4, 5, 1, 1)
        tabWidget.addTab(self.nrtl, "NRTL")

        self.wils = QtWidgets.QWidget()
        lyt_wils = QtWidgets.QGridLayout(self.wils)
        self.wils_widget = [
            Entrada_con_unidades(float), Entrada_con_unidades(float),
            Entrada_con_unidades(float), Entrada_con_unidades(float),
            Entrada_con_unidades(float), Entrada_con_unidades(float)]
        lyt_wils.addWidget(QtWidgets.QLabel("Aij"), 2, 1, 1, 1)
        lyt_wils.addWidget(self.wils_widget[0], 2, 2, 1, 1)
        lyt_wils.addWidget(QtWidgets.QLabel("Aji"), 2, 4, 1, 1)
        lyt_wils.addWidget(self.wils_widget[1], 2, 5, 1, 1)
        lyt_wils.addWidget(QtWidgets.QLabel("Bij"), 3, 1, 1, 1)
        lyt_wils.addWidget(self.wils_widget[2], 3, 2, 1, 1)
        lyt_wils.addWidget(QtWidgets.QLabel("Bji"), 3, 4, 1, 1)
        lyt_wils.addWidget(self.wils_widget[3], 3, 5, 1, 1)
        lyt_wils.addWidget(QtWidgets.QLabel("Cij"), 4, 1, 1, 1)
        lyt_wils.addWidget(self.wils_widget[4], 4, 2, 1, 1)
        lyt_wils.addWidget(QtWidgets.QLabel("Cji"), 4, 4, 1, 1)
        lyt_wils.addWidget(self.wils_widget[5], 4, 5, 1, 1)
        tabWidget.addTab(self.wils, "Wilson")

        self.uniq = QtWidgets.QWidget()
        lyt_uniq = QtWidgets.QGridLayout(self.uniq)
        self.uniq_widget = [
            Entrada_con_unidades(float), Entrada_con_unidades(float),
            Entrada_con_unidades(float), Entrada_con_unidades(float),
            Entrada_con_unidades(float), Entrada_con_unidades(float)]
        lyt_uniq.addWidget(QtWidgets.QLabel("Uij-Ujj"), 2, 1, 1, 1)
        lyt_uniq.addWidget(self.uniq_widget[0], 2, 2, 1, 1)
        lyt_uniq.addWidget(QtWidgets.QLabel("Uji-Uii"), 2, 4, 1, 1)
        lyt_uniq.addWidget(self.uniq_widget[1], 2, 5, 1, 1)
        lyt_uniq.addWidget(QtWidgets.QLabel("aij"), 3, 1, 1, 1)
        lyt_uniq.addWidget(self.uniq_widget[2], 3, 2, 1, 1)
        lyt_uniq.addWidget(QtWidgets.QLabel("aji"), 3, 4, 1, 1)
        lyt_uniq.addWidget(self.uniq_widget[3], 3, 5, 1, 1)
        lyt_uniq.addWidget(QtWidgets.QLabel("Cij"), 4, 1, 1, 1)
        lyt_uniq.addWidget(self.uniq_widget[4], 4, 2, 1, 1)
        lyt_uniq.addWidget(QtWidgets.QLabel("Cji"), 4, 4, 1, 1)
        lyt_uniq.addWidget(self.uniq_widget[5], 4, 5, 1, 1)
        tabWidget.addTab(self.uniq, "UNIQUAC")

        self.srk = QtWidgets.QWidget()
        lyt_srk = QtWidgets.QGridLayout(self.srk)
        self.srk_widget = [
            Entrada_con_unidades(float), Entrada_con_unidades(float),
            Entrada_con_unidades(float)]
        lyt_srk.addWidget(QtWidgets.QLabel("Aij"), 2, 1, 1, 1)
        lyt_srk.addWidget(self.srk_widget[0], 2, 2, 1, 1)
        lyt_srk.addWidget(QtWidgets.QLabel("Bij"), 3, 1, 1, 1)
        lyt_srk.addWidget(self.srk_widget[1], 3, 2, 1, 1)
        lyt_srk.addWidget(QtWidgets.QLabel("Cij"), 4, 1, 1, 1)
        lyt_srk.addWidget(self.srk_widget[2], 4, 2, 1, 1)
        tabWidget.addTab(self.srk, "SRK")

        self.pr = QtWidgets.QWidget()
        lyt_pr = QtWidgets.QGridLayout(self.pr)
        self.pr_widget = [
            Entrada_con_unidades(float), Entrada_con_unidades(float),
            Entrada_con_unidades(float)]
        lyt_pr.addWidget(QtWidgets.QLabel("Aij"), 2, 1, 1, 1)
        lyt_pr.addWidget(self.pr_widget[0], 2, 2, 1, 1)
        lyt_pr.addWidget(QtWidgets.QLabel("Bij"), 3, 1, 1, 1)
        lyt_pr.addWidget(self.pr_widget[1], 3, 2, 1, 1)
        lyt_pr.addWidget(QtWidgets.QLabel("Cij"), 4, 1, 1, 1)
        lyt_pr.addWidget(self.pr_widget[2], 4, 2, 1, 1)
        tabWidget.addTab(self.pr, "PR")

        self.bwrs = QtWidgets.QWidget()
        lyt_bwrs = QtWidgets.QGridLayout(self.bwrs)
        lyt_bwrs.addWidget(QtWidgets.QLabel("Aij"), 2, 1, 1, 1)
        self.bwrs_widget = [
            Entrada_con_unidades(float), Entrada_con_unidades(float),
            Entrada_con_unidades(float)]
        lyt_bwrs.addWidget(self.bwrs_widget[0], 2, 2, 1, 1)
        lyt_bwrs.addWidget(QtWidgets.QLabel("Bij"), 3, 1, 1, 1)
        lyt_bwrs.addWidget(self.bwrs_widget[1], 3, 2, 1, 1)
        lyt_bwrs.addWidget(QtWidgets.QLabel("Cij"), 4, 1, 1, 1)
        lyt_bwrs.addWidget(self.bwrs_widget[2], 4, 2, 1, 1)
        tabWidget.addTab(self.bwrs, "BWRS")

        self.id1.setText("%3i - %s" % (i, getElement(i)[2]))
        self.id2.setText("%3i - %s" % (j, getElement(j)[2]))

        for coef in ["nrtl", "wils", "uniq", "srk", "pr", "bwrs"]:
            dat = os.environ["pychemqt"] + "dat/bip/%s.dat" % coef
            widget = self.__getattribute__("%s_widget" % coef)
            load(dat, i, j, widget)


def load(file, i, j, widgets):
    with open(file, "r") as archivo:
        interacciones = archivo.readlines()
    for interaccion in interacciones:
        numeros = interaccion.split()
        if int(numeros[0]) == i and int(numeros[1]) == j:
            if len(widgets) == 3:
                widgets[0].setValue(str(float(numeros[2])))
                widgets[1].setValue(str(float(numeros[4])))
                widgets[2].setValue(str(float(numeros[6])))
            else:
                for k in range(len(widgets)):
                    widgets[k].setValue(str(float(numeros[k+2])))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = Ui_BIP(2, 49)
    Dialog.show()
    sys.exit(app.exec_())
