#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import inspect

from PyQt5 import QtCore, QtGui, QtWidgets


from lib import mEoS
from equipment import equipments
from tools import HelpView

objects = inspect.getmembers(mEoS)

for nombre, objeto in objects:
    if inspect.isclass(objeto):
        functions = inspect.getmembers(objeto)
#        print objeto, functions
        for name, function in functions:
            inspect.getdoc(function)
#            print inspect.getdoc(function)

#print os.path.curdir
#print glob.glob(os.path.curdir+"/*/*.py")


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

        self.buttonBox = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Close)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox)

        self.fill()
        self.tree.itemDoubleClicked.connect(self.open)

    def fill(self):
        """Fill tree with documentation entries"""
        # Equipment
        itemEquipment = QtWidgets.QTreeWidgetItem([QtWidgets.QApplication.translate(
            "pychemqt", "Equipments")])
        self.tree.addTopLevelItem(itemEquipment)
        for equip in equipments:
            itemequip = QtWidgets.QTreeWidgetItem([equip.__name__])
            itemEquipment.addChild(itemequip)
            for link in equip.__doi__:
                item = QtWidgets.QTreeWidgetItem([link["doi"], "%s: %s. %s" %(
                    link["autor"], link["title"], link["ref"])])
                itemequip.addChild(item)

        # MEoS
        link = []
        objects = [mEoS.MEoS]+mEoS.__all__

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
                        item.addChild(QtWidgets.QTreeWidgetItem([name, lastline]))

            listas = ["eq", "_viscosity", "_thermal"]
            for lista in listas:
                if lista in objeto.__dict__:
                    for eq in objeto.__dict__[lista]:
                        if eq and "__doi__" in eq:
                            item.addChild(QtWidgets.QTreeWidgetItem(
                                [eq["__name__"], eq["__doi__"]]))

    def open(self, item, int):
        """Open file if exist in doc/doi folder or open a browser in paper link"""
        if item.parent():
            text = item.text(0)
            code = str(text).replace("/", "_")
            file = "doc/doi/"+code+".pdf"
            if os.path.isfile(file):
                os.system('evince "'+file+'"')
            else:
                explorer = HelpView.HelpView(text, QtCore.QUrl(
                    "http://dx.doi.org/%s" %text))
                explorer.exec_()


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    dialog = ShowReference()
    dialog.show()
    sys.exit(app.exec_())
