#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Module to show optional dependences availability
###############################################################################

import os

from PyQt5 import QtGui, QtWidgets

optional_modules = (
    ("freesteam", QtWidgets.QApplication.translate(
        "pychemqt", "freesteam thermal option disabled")),
    ("oasa", QtWidgets.QApplication.translate(
        "pychemqt", "graphic formula disabled")),
    ("Elemental", QtWidgets.QApplication.translate(
        "pychemqt", "periodic table don´t available")),
    ("CoolProp", QtWidgets.QApplication.translate(
        "pychemqt", "coolprop thermal option disabled")),
    ("refprop", QtWidgets.QApplication.translate(
        "pychemqt", "refprop thermal option disabled")),
    ("ezodf", QtWidgets.QApplication.translate(
        "pychemqt", "openoffice/libreoffice interaction disabled")),
    ("openpyxl", QtWidgets.QApplication.translate(
        "pychemqt", "Microsoft Excel 2007/2010 interaction disabled")),
    ("xlwt", QtWidgets.QApplication.translate(
        "pychemqt", "Microsoft Excel 97/2000/XP/2003 interaction disabled")),
    )


class ShowDependences(QtWidgets.QDialog):
    """Dialog to show optional dependences availability"""
    def __init__(self, parent=None):
        super(ShowDependences, self).__init__(parent)
        self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"]+"/images/button/showPrograms.png")))
        self.setWindowTitle(
            QtWidgets.QApplication.translate("pychemqt", "External program"))
        layout = QtWidgets.QVBoxLayout(self)
        self.tree = QtWidgets.QTreeWidget()
        header = QtWidgets.QTreeWidgetItem(
            [QtWidgets.QApplication.translate("pychemqt", "Module"),
             QtWidgets.QApplication.translate("pychemqt", "Status")])
        self.tree.setHeaderItem(header)

        for module, txt in optional_modules:
            if os.environ[module] == "True":
                mod = __import__(module)
                estado = mod.__file__
            else:
                estado = QtWidgets.QApplication.translate("pychemqt", "not found")
            item = QtWidgets.QTreeWidgetItem([module, estado])
            self.tree.addTopLevelItem(item)

        layout.addWidget(self.tree)
        self.buttonBox = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Close)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox)


if __name__ == "__main__":
    import sys
    os.environ["pychemqt"] = "/home/jjgomera/pychemqt/"
    os.environ["freesteam"] = "False"
    os.environ["oasa"] = "False"
    os.environ["Elemental"] = "False"
    os.environ["CoolProp"] = "False"
    os.environ["refprop"] = "False"
    os.environ["ezodf"] = "False"
    os.environ["openpyxl"] = "False"
    os.environ["xlwt"] = "False"
    app = QtWidgets.QApplication(sys.argv)
    dialog = ShowDependences()
    dialog.show()
    app.exec_()
