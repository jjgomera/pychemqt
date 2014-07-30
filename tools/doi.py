#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import inspect

from PyQt4 import QtCore, QtGui

from lib import mEoS
from tools import HelpView

objects=inspect.getmembers(mEoS)

for nombre, objeto in objects:
    if inspect.isclass(objeto):
        functions=inspect.getmembers(objeto)
#        print objeto, functions
        for name, function in functions:
            inspect.getdoc(function)
#            print inspect.getdoc(function)

#print os.path.curdir
#print glob.glob(os.path.curdir+"/*/*.py")

class ShowReference(QtGui.QDialog):
    def __init__(self, parent=None):
        super(ShowReference, self).__init__(parent)
        self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/help.png")))
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Reference Paper Show Dialog"))
        layout = QtGui.QVBoxLayout(self)

        self.tree=QtGui.QTreeWidget()
        header=QtGui.QTreeWidgetItem([QtGui.QApplication.translate("pychemqt", "File"), QtGui.QApplication.translate("pychemqt", "Link")])
        self.tree.setHeaderItem(header)

        link=[]
        objects=[mEoS.MEoS]+mEoS.__all__
#        objects.insert(0,("MEoS", mEoS.MEoS))

        for objeto in objects:
            item=QtGui.QTreeWidgetItem([objeto.__name__, ""])
            self.tree.addTopLevelItem(item)
            functions=inspect.getmembers(objeto)
            for name, function in functions:
                doc=inspect.getdoc(function)
                if doc:
                    lastline=doc.split("\n")[-1]
                    if lastline[:17]=="http://dx.doi.org" and lastline not in link:
                        link.append(lastline)
                        item.addChild(QtGui.QTreeWidgetItem([name, lastline]))

            listas=["eq", "_viscosity", "_thermal"]
            for lista in listas:
                if lista in objeto.__dict__:
                    for eq in objeto.__dict__[lista]:
                        if eq and "__doi__" in eq:
                            item.addChild(QtGui.QTreeWidgetItem([eq["__name__"], eq["__doi__"]]))

        layout.addWidget(self.tree)

        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Close)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox)

        self.tree.itemDoubleClicked.connect(self.openpdf)

    def openpdf(self, item, int):
        if item.parent():
            text=item.text(self.tree.columnCount()-1)
            code=str(text).replace("http://","").replace("/","_")
            file="doc/doi/"+code+".pdf"
            if os.path.isfile(file):
                os.system('evince "'+file+'"')
            else:
                explorer=HelpView.HelpView("aa", QtCore.QUrl(text))
                explorer.exec_()





if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    dialog = ShowReference()
    dialog.show()
    sys.exit(app.exec_())

