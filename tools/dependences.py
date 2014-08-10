#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Module to show optional dependences availability
###############################################################################

import os

from PyQt4 import QtGui

optional_modules = (
    ("freesteam", QtGui.QApplication.translate("pychemqt", "freesteam thermal option disabled")),
    ("oasa", QtGui.QApplication.translate("pychemqt", "graphic formula disabled")),
    ("Elemental", QtGui.QApplication.translate("pychemqt", "periodic table don´t available")),
    ("CoolProp", QtGui.QApplication.translate("pychemqt", "coolprop thermal option disabled")),
    ("refprop", QtGui.QApplication.translate("pychemqt", "refprop thermal option disabled")),
    ("ezodf", QtGui.QApplication.translate("pychemqt", "openoffice/libreoffice interaction disabled")),
    ("openpyxl", QtGui.QApplication.translate("pychemqt", "Microsoft Excel 2007/2010 interaction disabled")),
    ("xlwt", QtGui.QApplication.translate("pychemqt", "Microsoft Excel 97/2000/XP/2003 interaction disabled")),
    )


class ShowDependences(QtGui.QDialog):
    """Dialog to show optional dependences availability"""
    def __init__(self, parent=None):
        super(ShowDependences, self).__init__(parent)
        self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"]+"/images/button/showPrograms.png")))
        self.setWindowTitle(
            QtGui.QApplication.translate("pychemqt", "External program"))
        layout = QtGui.QVBoxLayout(self)
        self.tree = QtGui.QTreeWidget()
        header = QtGui.QTreeWidgetItem(
            [QtGui.QApplication.translate("pychemqt", "Module"),
             QtGui.QApplication.translate("pychemqt", "Status")])
        self.tree.setHeaderItem(header)

        for module, txt in optional_modules:
            if os.environ[module] == "True":
                mod = __import__(module)
                estado = mod.__file__
            else:
                estado = QtGui.QApplication.translate("pychemqt", "not found")
            item = QtGui.QTreeWidgetItem([module, estado])
            self.tree.addTopLevelItem(item)

        layout.addWidget(self.tree)
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Close)
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
    app = QtGui.QApplication(sys.argv)
    dialog = ShowDependences()
    dialog.show()
    app.exec_()
