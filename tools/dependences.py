#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2022, Juan José Gómez Romera <jjgomera@gmail.com>

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
# Module to show optional dependences availability
###############################################################################

import os

from qt import QtGui, QtWidgets

optional_modules = (
    ("freesteam", QtWidgets.QApplication.translate(
        "pychemqt", "freesteam thermal option disabled")),
    ("CoolProp", QtWidgets.QApplication.translate(
        "pychemqt", "coolprop thermal option disabled")),
    ("refprop", QtWidgets.QApplication.translate(
        "pychemqt", "refprop thermal option disabled")),
    ("openbabel", QtWidgets.QApplication.translate(
        "pychemqt", "graphic formula disabled")),
    ("ezodf", QtWidgets.QApplication.translate(
        "pychemqt", "openoffice/libreoffice interaction disabled")),
    ("openpyxl", QtWidgets.QApplication.translate(
        "pychemqt", "Microsoft Excel 2007/2010 interaction disabled")),
    ("xlwt", QtWidgets.QApplication.translate(
        "pychemqt", "Microsoft Excel 97/2000/XP/2003 interaction disabled")),
    ("icu", QtWidgets.QApplication.translate(
        "pychemqt",
        "Unicode collation algorithm for improved string sorting disabled")),
    ("reportlab", QtWidgets.QApplication.translate(
        "pychemqt", "Pdf report exporting disabled")),
    ("PyQt5.Qsci", QtWidgets.QApplication.translate(
        "pychemqt", "Qscintilla custom module editor disabled")),
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
                st = mod.__file__
            else:
                st = QtWidgets.QApplication.translate("pychemqt", "not found")
            item = QtWidgets.QTreeWidgetItem([module, st])
            self.tree.addTopLevelItem(item)

        layout.addWidget(self.tree)
        button = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Close)
        button.rejected.connect(self.reject)
        layout.addWidget(button)


if __name__ == "__main__":
    import sys
    os.environ["pychemqt"] = "/home/jjgomera/pychemqt/"
    os.environ["freesteam"] = "False"
    os.environ["openbabel"] = "False"
    os.environ["CoolProp"] = "False"
    os.environ["refprop"] = "False"
    os.environ["ezodf"] = "False"
    os.environ["openpyxl"] = "False"
    os.environ["xlwt"] = "False"
    os.environ["icu"] = "False"
    os.environ["reportlab"] = "False"
    os.environ["PyQt5.Qsci"] = "True"
    app = QtWidgets.QApplication(sys.argv)
    dialog = ShowDependences()
    dialog.show()
    app.exec()
