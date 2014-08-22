#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Module for UI binary interaction parameter editor/viewer
###############################################################################

import os

from PyQt4 import QtGui

#from lib.databank import base_datos
from UI.widgets import Entrada_con_unidades


class Ui_BIP(QtGui.QDialog):
    """Dialog to view/edit the BIP for EoS available"""
    def __init__(self, i, j, parent=None):
        """Constructor
        i, j are the index for component to BIP show"""
        super(Ui_BIP, self).__init__(parent)
        self.setWindowTitle(QtGui.QApplication.translate(
            "pychemqt", "BIP (Binary interaction parameters)"))

        gridLayout = QtGui.QGridLayout(self)
        gridLayout.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Component")+" 1:"), 0, 0, 1, 1)
        gridLayout.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Component")+" 2:"), 1, 0, 1, 1)
        self.id1 = QtGui.QLabel()
        gridLayout.addWidget(self.id1, 0, 1, 1, 1)
        self.id2 = QtGui.QLabel()
        gridLayout.addWidget(self.id2, 1, 1, 1, 1)
        self.tabWidget = QtGui.QTabWidget()
        gridLayout.addWidget(self.tabWidget, 2, 0, 1, 3)
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Cancel |
                                                QtGui.QDialogButtonBox.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        gridLayout.addWidget(self.buttonBox, 3, 2, 1, 1)

        self.NRTL = QtGui.QWidget()
        gridLayout_NRTL = QtGui.QGridLayout(self.NRTL)
        self.NRTL_widget = [
            Entrada_con_unidades(float), Entrada_con_unidades(float),
            Entrada_con_unidades(float), Entrada_con_unidades(float),
            Entrada_con_unidades(float), Entrada_con_unidades(float)]
        gridLayout_NRTL.addWidget(QtGui.QLabel("Bij"), 2, 1, 1, 1)
        gridLayout_NRTL.addWidget(self.NRTL_widget[0], 2, 2, 1, 1)
        gridLayout_NRTL.addWidget(QtGui.QLabel("Bji"), 2, 4, 1, 1)
        gridLayout_NRTL.addWidget(self.NRTL_widget[1], 2, 5, 1, 1)
        gridLayout_NRTL.addWidget(QtGui.QLabel(u"αij"), 3, 1, 1, 1)
        gridLayout_NRTL.addWidget(self.NRTL_widget[2], 3, 2, 1, 1)
        gridLayout_NRTL.addWidget(QtGui.QLabel(u"αji"), 3, 4, 1, 1)
        gridLayout_NRTL.addWidget(self.NRTL_widget[3], 3, 5, 1, 1)
        gridLayout_NRTL.addWidget(QtGui.QLabel("Aij"), 4, 1, 1, 1)
        gridLayout_NRTL.addWidget(self.NRTL_widget[4], 4, 2, 1, 1)
        gridLayout_NRTL.addWidget(QtGui.QLabel("Aji"), 4, 4, 1, 1)
        gridLayout_NRTL.addWidget(self.NRTL_widget[5], 4, 5, 1, 1)
        self.tabWidget.addTab(self.NRTL, "NRTL")

        self.Wilson = QtGui.QWidget()
        gridLayout_Wilson = QtGui.QGridLayout(self.Wilson)
        self.wilson_widget = [
            Entrada_con_unidades(float), Entrada_con_unidades(float),
            Entrada_con_unidades(float), Entrada_con_unidades(float),
            Entrada_con_unidades(float), Entrada_con_unidades(float)]
        gridLayout_Wilson.addWidget(QtGui.QLabel("Aij"), 2, 1, 1, 1)
        gridLayout_Wilson.addWidget(self.wilson_widget[0], 2, 2, 1, 1)
        gridLayout_Wilson.addWidget(QtGui.QLabel("Aji"), 2, 4, 1, 1)
        gridLayout_Wilson.addWidget(self.wilson_widget[1], 2, 5, 1, 1)
        gridLayout_Wilson.addWidget(QtGui.QLabel("Bij"), 3, 1, 1, 1)
        gridLayout_Wilson.addWidget(self.wilson_widget[2], 3, 2, 1, 1)
        gridLayout_Wilson.addWidget(QtGui.QLabel("Bji"), 3, 4, 1, 1)
        gridLayout_Wilson.addWidget(self.wilson_widget[3], 3, 5, 1, 1)
        gridLayout_Wilson.addWidget(QtGui.QLabel("Cij"), 4, 1, 1, 1)
        gridLayout_Wilson.addWidget(self.wilson_widget[4], 4, 2, 1, 1)
        gridLayout_Wilson.addWidget(QtGui.QLabel("Cji"), 4, 4, 1, 1)
        gridLayout_Wilson.addWidget(self.wilson_widget[5], 4, 5, 1, 1)
        self.tabWidget.addTab(self.Wilson, "Wilson")

        self.UNIQUAC = QtGui.QWidget()
        gridLayout_UNIQUAC = QtGui.QGridLayout(self.UNIQUAC)
        self.UNIQUAC_widget = [
            Entrada_con_unidades(float), Entrada_con_unidades(float),
            Entrada_con_unidades(float), Entrada_con_unidades(float),
            Entrada_con_unidades(float), Entrada_con_unidades(float)]
        gridLayout_UNIQUAC.addWidget(QtGui.QLabel("Uij-Ujj"), 2, 1, 1, 1)
        gridLayout_UNIQUAC.addWidget(self.UNIQUAC_widget[0], 2, 2, 1, 1)
        gridLayout_UNIQUAC.addWidget(QtGui.QLabel("Uji-Uii"), 2, 4, 1, 1)
        gridLayout_UNIQUAC.addWidget(self.UNIQUAC_widget[1], 2, 5, 1, 1)
        gridLayout_UNIQUAC.addWidget(QtGui.QLabel("aij"), 3, 1, 1, 1)
        gridLayout_UNIQUAC.addWidget(self.UNIQUAC_widget[2], 3, 2, 1, 1)
        gridLayout_UNIQUAC.addWidget(QtGui.QLabel("aji"), 3, 4, 1, 1)
        gridLayout_UNIQUAC.addWidget(self.UNIQUAC_widget[3], 3, 5, 1, 1)
        gridLayout_UNIQUAC.addWidget(QtGui.QLabel("Cij"), 4, 1, 1, 1)
        gridLayout_UNIQUAC.addWidget(self.UNIQUAC_widget[4], 4, 2, 1, 1)
        gridLayout_UNIQUAC.addWidget(QtGui.QLabel("Cji"), 4, 4, 1, 1)
        gridLayout_UNIQUAC.addWidget(self.UNIQUAC_widget[5], 4, 5, 1, 1)
        self.tabWidget.addTab(self.UNIQUAC, "UNIQUAC")

        self.SRK = QtGui.QWidget()
        gridLayout_SRK = QtGui.QGridLayout(self.SRK)
        self.SRK_widget = [
            Entrada_con_unidades(float), Entrada_con_unidades(float),
            Entrada_con_unidades(float)]
        gridLayout_SRK.addWidget(QtGui.QLabel("Aij"), 2, 1, 1, 1)
        gridLayout_SRK.addWidget(self.SRK_widget[0], 2, 2, 1, 1)
        gridLayout_SRK.addWidget(QtGui.QLabel("Bij"), 3, 1, 1, 1)
        gridLayout_SRK.addWidget(self.SRK_widget[1], 3, 2, 1, 1)
        gridLayout_SRK.addWidget(QtGui.QLabel("Cij"), 4, 1, 1, 1)
        gridLayout_SRK.addWidget(self.SRK_widget[2], 4, 2, 1, 1)
        self.tabWidget.addTab(self.SRK, "SRK")

        self.PR = QtGui.QWidget()
        gridLayout_PR = QtGui.QGridLayout(self.PR)
        self.PR_widget = [
            Entrada_con_unidades(float), Entrada_con_unidades(float),
            Entrada_con_unidades(float)]
        gridLayout_PR.addWidget(QtGui.QLabel("Aij"), 2, 1, 1, 1)
        gridLayout_PR.addWidget(self.PR_widget[0], 2, 2, 1, 1)
        gridLayout_PR.addWidget(QtGui.QLabel("Bij"), 3, 1, 1, 1)
        gridLayout_PR.addWidget(self.PR_widget[1], 3, 2, 1, 1)
        gridLayout_PR.addWidget(QtGui.QLabel("Cij"), 4, 1, 1, 1)
        gridLayout_PR.addWidget(self.PR_widget[2], 4, 2, 1, 1)
        self.tabWidget.addTab(self.PR, "Peng-Robinson")

        self.BWRS = QtGui.QWidget()
        gridLayout_BWRS = QtGui.QGridLayout(self.BWRS)
        gridLayout_BWRS.addWidget(QtGui.QLabel("Aij"), 2, 1, 1, 1)
        self.BWRS_widget = [
            Entrada_con_unidades(float), Entrada_con_unidades(float),
            Entrada_con_unidades(float)]
        gridLayout_BWRS.addWidget(self.BWRS_widget[0], 2, 2, 1, 1)
        gridLayout_BWRS.addWidget(QtGui.QLabel("Bij"), 3, 1, 1, 1)
        gridLayout_BWRS.addWidget(self.BWRS_widget[1], 3, 2, 1, 1)
        gridLayout_BWRS.addWidget(QtGui.QLabel("Cij"), 4, 1, 1, 1)
        gridLayout_BWRS.addWidget(self.BWRS_widget[2], 4, 2, 1, 1)
        self.tabWidget.addTab(self.BWRS, "BWRS")

        #self.id1.setText("%3s - %s" %(str(i), str(base_datos[i][1])))
        #self.id2.setText("%3s - %s" %(str(j), str(base_datos[j][1])))

        #load(os.environ["pychemqt"] + "dat/bip/nrtl.dat", i, j, self.NRTL_widget)
        #load(os.environ["pychemqt"] + "dat/bip/wils.dat", i, j, self.wilson_widget)
        #load(os.environ["pychemqt"] + "dat/bip/uniq.dat", i, j, self.UNIQUAC_widget)
        #load(os.environ["pychemqt"] + "dat/bip/srk.dat", i, j, self.SRK_widget)
        #load(os.environ["pychemqt"] + "dat/bip/pr.dat", i, j, self.PR_widget)
        #load(os.environ["pychemqt"] + "dat/bip/bwrs.dat", i, j, self.BWRS_widget)


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
    app = QtGui.QApplication(sys.argv)
    Dialog = Ui_BIP(2, 49)
    Dialog.show()
    sys.exit(app.exec_())
