#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Tools with units converter
###############################################################################

from PyQt4 import QtGui

from lib import unidades
from UI.conversor_unidades import UI_conversorUnidades, moneda


class UI_unitConverter(QtGui.QDialog):
    def __init__(self, parent=None):
        super(UI_unitConverter, self).__init__(parent)
        self.setWindowTitle(
            QtGui.QApplication.translate("pychemqt", "Units converter"))

        self.verticalLayout = QtGui.QVBoxLayout(self)
        self.lista = QtGui.QListWidget()
        self.lista.itemDoubleClicked.connect(self.mostrar_ventana_hijo)
        self.verticalLayout.addWidget(self.lista)
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Close)
        self.buttonBox.rejected.connect(self.reject)
        self.verticalLayout.addWidget(self.buttonBox)
        for unidad in unidades._all:
            self.lista.addItem(unidad.__title__)
        self.lista.setCurrentRow(-1)

    def mostrar_ventana_hijo(self):
        """Show child window with selected unit converter"""
        indice = self.lista.currentRow()
        if unidades._all[indice].__name__ == "Currency":
            dialog = moneda()
        else:
            dialog = UI_conversorUnidades(unidades._all[indice])
        dialog.exec_()

if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    conversion_unidades = UI_unitConverter()
    conversion_unidades.show()
    sys.exit(app.exec_())
