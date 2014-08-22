#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import inspect

from PyQt4 import QtCore, QtGui

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


class ShowReference(QtGui.QDialog):
    def __init__(self, parent=None):
        super(ShowReference, self).__init__(parent)
        self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"]+"/images/button/help.png")))
        self.setWindowTitle(QtGui.QApplication.translate(
            "pychemqt", "Reference Paper Show Dialog"))
        layout = QtGui.QVBoxLayout(self)

        self.tree = QtGui.QTreeWidget()
        header = QtGui.QTreeWidgetItem(
            [QtGui.QApplication.translate("pychemqt", "File"),
             QtGui.QApplication.translate("pychemqt", "Description")])
        self.tree.setHeaderItem(header)
        layout.addWidget(self.tree)

        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Close)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox)

        self.fill()
        self.tree.itemDoubleClicked.connect(self.open)

    def fill(self):
        """Fill tree with documentation entries"""
        # Equipment
        itemEquipment = QtGui.QTreeWidgetItem([QtGui.QApplication.translate(
            "pychemqt", "Equipments")])
        self.tree.addTopLevelItem(itemEquipment)
        for equip in equipments:
            itemequip = QtGui.QTreeWidgetItem([equip.__name__])
            itemEquipment.addChild(itemequip)
            for link in equip.__doi__:
                item = QtGui.QTreeWidgetItem([link["doi"], "%s: %s. %s" %(
                    link["autor"], link["title"], link["ref"])])
                itemequip.addChild(item)

        # MEoS
        link = []
        objects = [mEoS.MEoS]+mEoS.__all__

        for objeto in objects:
            item = QtGui.QTreeWidgetItem([objeto.__name__, ""])
            self.tree.addTopLevelItem(item)
            functions = inspect.getmembers(objeto)
            for name, function in functions:
                doc = inspect.getdoc(function)
                if doc:
                    lastline = doc.split("\n")[-1]
                    if lastline[:17] == "http://dx.doi.org" and \
                            lastline not in link:
                        link.append(lastline)
                        item.addChild(QtGui.QTreeWidgetItem([name, lastline]))

            listas = ["eq", "_viscosity", "_thermal"]
            for lista in listas:
                if lista in objeto.__dict__:
                    for eq in objeto.__dict__[lista]:
                        if eq and "__doi__" in eq:
                            item.addChild(QtGui.QTreeWidgetItem(
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
    app = QtGui.QApplication(sys.argv)
    dialog = ShowReference()
    dialog.show()
    sys.exit(app.exec_())
