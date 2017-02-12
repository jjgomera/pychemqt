#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2016, Juan José Gómez Romera <jjgomera@gmail.com>

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


import os
import inspect

from PyQt5 import QtCore, QtGui, QtWidgets


from lib import meos, mEoS
from equipment import equipments

objects = inspect.getmembers(mEoS)

for nombre, objeto in objects:
    if inspect.isclass(objeto):
        functions = inspect.getmembers(objeto)
#        print objeto, functions
        for name, function in functions:
            inspect.getdoc(function)
#            print inspect.getdoc(function)

# print os.path.curdir
# print glob.glob(os.path.curdir+"/*/*.py")


class ShowReference(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super(ShowReference, self).__init__(parent)
        self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"]+"/images/button/help.png")))
        self.setWindowTitle(QtWidgets.QApplication.translate(
            "pychemqt", "Reference Paper Show Dialog"))
        layout = QtWidgets.QVBoxLayout(self)

        self.tree = QtWidgets.QTreeWidget()
        header = QtWidgets.QTreeWidgetItem(
            [QtWidgets.QApplication.translate("pychemqt", "File"),
             QtWidgets.QApplication.translate("pychemqt", "Description")])
        self.tree.setHeaderItem(header)
        layout.addWidget(self.tree)

        bttBox = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Close)
        bttBox.rejected.connect(self.reject)
        layout.addWidget(bttBox)

        self.fill()
        self.tree.itemDoubleClicked.connect(self.open)

    def fill(self):
        """Fill tree with documentation entries"""
        # Equipment
        itemEquipment = QtWidgets.QTreeWidgetItem(
            [QtWidgets.QApplication.translate("pychemqt", "Equipments")])
        self.tree.addTopLevelItem(itemEquipment)
        for equip in equipments:
            itemequip = QtWidgets.QTreeWidgetItem([equip.__name__])
            itemEquipment.addChild(itemequip)
            for link in equip.__doi__:
                item = QtWidgets.QTreeWidgetItem([link["doi"], "%s: %s. %s" % (
                    link["autor"], link["title"], link["ref"])])
                itemequip.addChild(item)

        # MEoS
        link = []
        objects = [meos.MEoS]+mEoS.__all__

        for objeto in objects:
            item = QtWidgets.QTreeWidgetItem([objeto.__name__, ""])
            self.tree.addTopLevelItem(item)
            functions = inspect.getmembers(objeto)
            for name, function in functions:
                doc = inspect.getdoc(function)
                if doc:
                    lastline = doc.split("\n")[-1]
                    if lastline[:17] == "http://dx.doi.org" and \
                            lastline not in link:
                        link.append(lastline)
                        child = QtWidgets.QTreeWidgetItem([name, lastline])
                        item.addChild(child)

            listas = ["eq", "_viscosity", "_thermal"]
            for lista in listas:
                if lista in objeto.__dict__ and objeto.__dict__[lista]:
                    for eq in objeto.__dict__[lista]:
                        if eq and "__doi__" in eq:
                            item.addChild(QtWidgets.QTreeWidgetItem(
                                [eq["__doi__"]["doi"], "%s: %s, %s" % (
                                    eq["__doi__"]["autor"],
                                    eq["__doi__"]["title"],
                                    eq["__doi__"]["ref"])]))

    def open(self, item, int):
        """Open file if exist in doc/doi folder or open a browser with link"""
        if item.parent():
            text = item.text(0)
            code = str(text).replace("/", "_")
            file = "doc/doi/"+code+".pdf"
            if os.path.isfile(file):
                os.system('evince "'+file+'"')
            else:
                url = QtCore.QUrl("http://dx.doi.org/%s" % text)
                QtGui.QDesktopServices.openUrl(url)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    dialog = ShowReference()
    dialog.show()
    sys.exit(app.exec_())
