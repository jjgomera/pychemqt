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
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


###############################################################################
# Module to show optional dependences availability
###############################################################################

import os

from tools.qt import QtGui, QtWidgets, __qt__, translate


optional_modules = (
    ("freesteam", translate(
        "dependences", "freesteam thermal option disabled")),
    ("CoolProp", translate("dependences", "coolprop thermal option disabled")),
    ("refprop", translate("dependences", "refprop thermal option disabled")),
    ("openbabel", translate("dependences", "graphic formula disabled")),
    ("ezodf", translate(
        "dependences", "openoffice/libreoffice interaction disabled")),
    ("openpyxl", translate(
        "dependences", "Microsoft Excel 2007/2010 interaction disabled")),
    ("xlwt", translate(
        "dependences", "Microsoft Excel 97/2000/XP/2003 interaction disabled")),
    ("icu", translate(
        "dependences", "Unicode collation algorithm for improved string "
        "sorting disabled")),
    ("reportlab", translate("dependences", "Pdf report exporting disabled")),
    ("Qsci", translate(
        "dependences", "Qscintilla custom module editor disabled")))


class ShowDependences(QtWidgets.QDialog):
    """Dialog to show optional dependences availability"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(os.path.join(
            os.environ["pychemqt"], "images", "button", "showPrograms.png"))))
        self.setWindowTitle(self.tr("External program"))
        layout = QtWidgets.QVBoxLayout(self)
        self.tree = QtWidgets.QTreeWidget()
        header = QtWidgets.QTreeWidgetItem(
            [self.tr("Module"), self.tr("Status")])
        self.tree.setHeaderItem(header)

        for module, txt in optional_modules:
            if os.environ[module] == "True":
                if module == "Qsci":
                    # Special case for Qsci, a optional module from qt
                    if __qt__ == 6:
                        mod = __import__("PyQt6.Qsci")
                    else:
                        mod = __import__("PyQt5.Qsci")
                else:
                    mod = __import__(module)
                st = mod.__file__
                icon = QtGui.QIcon(QtGui.QPixmap(os.path.join(
                    os.environ["pychemqt"], "images", "button", "ok.png")))
            else:
                st = self.tr("Module not found")
                st += ", " + txt
                icon = QtGui.QIcon(QtGui.QPixmap(os.path.join(
                    os.environ["pychemqt"], "images", "button",
                    "fileClose.png")))
            item = QtWidgets.QTreeWidgetItem([module, st])
            item.setIcon(0, icon)
            self.tree.addTopLevelItem(item)

        self.tree.resizeColumnToContents(0)
        self.tree.resizeColumnToContents(1)
        self.resize(800, 300)
        layout.addWidget(self.tree)
        button = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Close)
        button.rejected.connect(self.reject)
        layout.addWidget(button)
