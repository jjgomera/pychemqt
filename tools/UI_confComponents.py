#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Component list tools
###############################################################################

import os
from PyQt5 import QtCore, QtGui, QtWidgets


from tools.UI_databank import UI_databank_widget


class UI_confComponents_widget(QtWidgets.QWidget):
    """Component list widget to use in dialor and in wizard"""
    componentChanged = QtCore.pyqtSignal("bool")

    def __init__(self, config=None, parent=None):
        """Constructor, opcional config parameter with project config"""
        super(UI_confComponents_widget, self).__init__(parent)
        layout = QtWidgets.QGridLayout(self)
        self.databank = UI_databank_widget()
        self.databank.BaseDatos.itemSelectionChanged.connect(
            self.comprobarBotones)
        layout.addWidget(self.databank, 1, 1, 17, 1)

        layout.addItem(QtWidgets.QSpacerItem(30, 30, QtWidgets.QSizePolicy.Expanding,
                                         QtWidgets.QSizePolicy.Fixed), 1, 2, 1, 2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Components list")), 2, 3)
        self.DeleteComponente = QtWidgets.QToolButton()
        self.DeleteComponente.setIcon(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"] + "/images/button/arrow-left.png")))
        self.DeleteComponente.clicked.connect(self.Delete)
        layout.addWidget(self.DeleteComponente, 4, 2)
        self.AddComponente = QtWidgets.QToolButton()
        self.AddComponente.setIcon(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"] + "/images/button/arrow-right.png")))
        self.AddComponente.clicked.connect(self.Add)
        layout.addWidget(self.AddComponente, 5, 2)
        self.Arriba = QtWidgets.QToolButton()
        self.Arriba.setIcon(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"] + "/images/button/arrow-up.png")))
        self.Arriba.clicked.connect(self.Up)
        layout.addWidget(self.Arriba, 6, 2)
        self.Abajo = QtWidgets.QToolButton()
        self.Abajo.setIcon(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"] + "/images/button/arrow-down.png")))
        self.Abajo.clicked.connect(self.Down)
        layout.addWidget(self.Abajo, 7, 2)
        self.clearComp = QtWidgets.QToolButton()
        self.clearComp.setIcon(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"] + "/images/button/clear.png")))
        self.clearComp.clicked.connect(self.clear)
        layout.addWidget(self.clearComp, 8, 2)

        self.ListaComponentes = QtWidgets.QListWidget()
        self.ListaComponentes.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                            QtWidgets.QSizePolicy.Expanding)
        self.ListaComponentes.itemSelectionChanged.connect(self.comprobarBotones)
        layout.addWidget(self.ListaComponentes, 3, 3, 7, 1)

        self.DeleteSolido = QtWidgets.QToolButton()
        self.DeleteSolido.setIcon(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"] + "/images/button/arrow-left.png")))
        self.DeleteSolido.clicked.connect(self.DeleteSolid)
        layout.addWidget(self.DeleteSolido, 13, 2)
        self.AddSolido = QtWidgets.QToolButton()
        self.AddSolido.setIcon(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"] + "/images/button/arrow-right.png")))
        self.AddSolido.clicked.connect(self.AddSolid)
        layout.addWidget(self.AddSolido, 14, 2)
        self.clearSolido = QtWidgets.QToolButton()
        self.clearSolido.setIcon(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"] + "/images/button/clear.png")))
        self.clearSolido.clicked.connect(self.clearSolids)
        layout.addWidget(self.clearSolido, 15, 2)

        layout.addItem(QtWidgets.QSpacerItem(20, 20, QtWidgets.QSizePolicy.Fixed,
                                         QtWidgets.QSizePolicy.Expanding), 10, 4, 1, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Solids", None)), 11, 3)
        self.ListaSolidos = QtWidgets.QListWidget()
        self.ListaSolidos.setFixedHeight(100)
        self.ListaSolidos.itemSelectionChanged.connect(self.comprobarBotones)
        layout.addWidget(self.ListaSolidos, 12, 3, 5, 1)

        self.indices = []
        self.solidos = []

        if config and config.has_section("Components"):
            try:
                self.indices = eval(config.get("Components", "Components"))
                self.solidos = eval(config.get("Components", "Solids"))
            except:
                self.indices = config.get("Components", "Components")
                self.solidos = config.get("Components", "Solids")
            for indice in self.indices:
                self.ListaComponentes.addItem("%03i - %s" % (
                    indice, self.databank.BaseDatos.item(indice-1, 1).text()))
            for solido in self.solidos:
                self.ListaSolidos.addItem("%03i - %s" % (
                    solido, self.databank.BaseDatos.item(solido-1, 1).text()))

        self.comprobarBotones()

    def Add(self):
        """Add selected component to component list"""
        insertar = self.ListaComponentes.currentRow()+1
        indice = self.databank.currentRow()+1
        if indice not in self.indices + self.solidos:
            self.indices.insert(insertar, indice)
            self.ListaComponentes.insertItem(insertar, "%03i - %s" % (
                indice, self.databank.BaseDatos.item(indice-1, 1).text()))
            self.comprobarBotones()
            self.ListaComponentes.setCurrentRow(insertar)
            self.componentChanged.emit(bool(self.indices))

    def Delete(self):
        """Remove component to component list"""
        del self.indices[self.ListaComponentes.currentRow()]
        self.ListaComponentes.takeItem(self.ListaComponentes.currentRow())
        self.comprobarBotones()
        self.componentChanged.emit(bool(self.indices))

    def Up(self):
        """Up selected item in component list"""
        indice = self.ListaComponentes.currentRow()
        self.ListaComponentes.insertItem(
            indice-1, self.ListaComponentes.item(indice).text())
        self.ListaComponentes.takeItem(indice+1)
        self.ListaComponentes.setCurrentRow(indice-1)
        self.indices.insert(indice-1, self.indices.pop(indice))

    def Down(self):
        """Down selected item in component list"""
        indice = self.ListaComponentes.currentRow()
        self.ListaComponentes.insertItem(
            indice+2, self.ListaComponentes.item(indice).text())
        self.ListaComponentes.takeItem(indice)
        self.ListaComponentes.setCurrentRow(indice+1)
        self.indices.insert(indice+1, self.indices.pop(indice))

    def AddSolid(self):
        """Add selected component to solid component list"""
        insertar = self.ListaSolidos.currentRow()+1
        indice = self.databank.currentRow()+1
        if indice not in self.indices+self.solidos:
            self.solidos.insert(insertar, indice)
            self.ListaSolidos.insertItem(insertar, "%03i - %s" % (
                indice, self.databank.BaseDatos.item(indice-1, 1).text()))
            self.comprobarBotones()
            self.ListaSolidos.setCurrentRow(insertar)

    def DeleteSolid(self):
        """Remove component to solid component list"""
        del self.solidos[self.ListaSolidos.currentRow()]
        self.ListaSolidos.takeItem(self.ListaSolidos.currentRow())
        self.comprobarBotones()

    def clear(self):
        """Clear component list"""
        self.ListaComponentes.clear()
        self.indices = []
        self.comprobarBotones()
        self.componentChanged.emit(bool(self.indices))

    def clearSolids(self):
        """Clear solid component list"""
        self.ListaSolidos.clear()
        self.solidos = []
        self.comprobarBotones()

    def comprobarBotones(self):
        """Check current status of button availibility"""
        if self.ListaComponentes.currentRow() == -1:
            self.Arriba.setEnabled(False)
            self.Abajo.setEnabled(False)
            self.DeleteComponente.setEnabled(False)
        else:
            self.DeleteComponente.setEnabled(True)
            if self.ListaComponentes.count() > 1:
                if self.ListaComponentes.currentRow() == 0:
                    self.Arriba.setEnabled(False)
                    self.Abajo.setEnabled(True)
                elif self.ListaComponentes.currentRow() == \
                        self.ListaComponentes.count()-1:
                    self.Abajo.setEnabled(False)
                    self.Arriba.setEnabled(True)
                else:
                    self.Arriba.setEnabled(True)
                    self.Abajo.setEnabled(True)
            else:
                self.Arriba.setEnabled(False)
                self.Abajo.setEnabled(False)

        self.clearComp.setEnabled(bool(self.indices))
        self.clearSolido.setEnabled(bool(self.solidos))

        if self.databank.currentRow() == -1 or \
                self.databank.currentRow() + 1 in self.indices + self.solidos:
            self.AddComponente.setEnabled(False)
            self.AddSolido.setEnabled(False)
        else:
            self.AddComponente.setEnabled(True)
            self.AddSolido.setEnabled(True)

        if self.ListaSolidos.currentRow() == -1:
            self.DeleteSolido.setEnabled(False)
        else:
            self.DeleteSolido.setEnabled(True)

    def value(self, config):
        """Function to result wizard"""
        if not config.has_section("Components"):
            config.add_section("Components")
        config.set("Components", "Components", self.indices)
        config.set("Components", "Solids", self.solidos)
        return config

    @classmethod
    def default(cls, config):
        config.add_section("Components")
        config.set("Components", "Components", "[]")
        config.set("Components", "Solids", "[]")
        return config


class Dialog(QtWidgets.QDialog):
    """Component list config dialog"""
    def __init__(self, config=None, parent=None):
        super(Dialog, self).__init__(parent)
        self.setWindowTitle(QtWidgets.QApplication.translate(
            "pychemqt", "Define project Components"))
        layout = QtWidgets.QVBoxLayout(self)
        self.datos = UI_confComponents_widget(config)
        layout.addWidget(self.datos)
        self.buttonBox = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Cancel |
                                                QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox)

    def value(self, config):
        """Function to result wizard"""
        config = self.datos.value(config)
        return config


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = Dialog()
    Dialog.show()
    sys.exit(app.exec_())
