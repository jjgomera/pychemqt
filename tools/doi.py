#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2017, Juan José Gómez Romera <jjgomera@gmail.com>

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
import subprocess

from PyQt5 import QtCore, QtGui, QtWidgets

import lib
from lib.config import IMAGE_PATH
# from equipment import equipments


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
            ["id",
             QtWidgets.QApplication.translate("pychemqt", "Autor"),
             QtWidgets.QApplication.translate("pychemqt", "Title"),
             QtWidgets.QApplication.translate("pychemqt", "Reference"),
             QtWidgets.QApplication.translate("pychemqt", "doi")])
        self.tree.setHeaderItem(header)
        layout.addWidget(self.tree)

        bttBox = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Close)
        bttBox.rejected.connect(self.reject)
        layout.addWidget(bttBox)

        self.fill()
        self.tree.itemDoubleClicked.connect(self.open)

    def fill(self):
        """Fill tree with documentation entries"""
        # General population Library
        for library in lib.__all__:
            __import__("lib.%s" % library)
            module = lib.__getattribute__(library)
            if hasattr(module, "__doi__") and module.__doi__:
                itemModule = QtWidgets.QTreeWidgetItem(["lib/"+library])
                self.tree.addTopLevelItem(itemModule)
                for key in sorted(module.__doi__.keys()):

                    # Special case for submodules
                    if "EoS" in library:
                        itemSubModule = QtWidgets.QTreeWidgetItem([key])
                        itemModule.addChildren([itemSubModule])
                        for key2 in sorted(module.__doi__[key].keys()):
                            link = module.__doi__[key][key2]
                            if library == "EoS":
                                title = ""
                            else:
                                title = key2.replace("_", "")
                            item = QtWidgets.QTreeWidgetItem([
                                title, link["autor"],
                                link["title"], link["ref"], link["doi"]])
                            code = link["doi"].replace("/", "_")
                            file = os.path.join("doc", "doi", code) + ".pdf"
                            file2 = os.path.join(
                                    "doc", "doi", link["title"]) + ".pdf"
                            if os.path.isfile(file) or os.path.isfile(file2):
                                icon = QtGui.QIcon(QtGui.QPixmap(os.path.join(
                                    IMAGE_PATH, "button", "ok.png")))
                                item.setIcon(0, icon)
                            itemSubModule.addChild(item)

                    else:
                        link = module.__doi__[key]
                        item = QtWidgets.QTreeWidgetItem([
                            "", link["autor"], link["title"], link["ref"],
                            link["doi"]])

                        code = link["doi"].replace("/", "_")
                        file = os.path.join("doc", "doi", code) + ".pdf"
                        file2 = os.path.join(
                                "doc", "doi", link["title"]) + ".pdf"
                        if os.path.isfile(file) or os.path.isfile(file2):
                            icon = QtGui.QIcon(QtGui.QPixmap(os.path.join(
                                IMAGE_PATH, "button", "ok.png")))
                            item.setIcon(0, icon)
                        itemModule.addChild(item)

        # Equipment
        # itemEquipment = QtWidgets.QTreeWidgetItem(
            # [QtWidgets.QApplication.translate("pychemqt", "Equipments")])
        # self.tree.addTopLevelItem(itemEquipment)
        # for equip in equipments:
            # itemequip = QtWidgets.QTreeWidgetItem([equip.__name__])
            # itemEquipment.addChild(itemequip)
            # for link in equip.__doi__:
                # item = QtWidgets.QTreeWidgetItem([link["doi"], "%s: %s. %s" % (
                    # link["autor"], link["title"], link["ref"])])
                # itemequip.addChild(item)

    def open(self, item, int):
        """Open file if exist in doc/doi folder or open a browser with link"""
        if item.parent() and not item.icon(0).isNull():
            title = item.text(2)
            text = item.text(4)
            code = str(text).replace("/", "_")
            file = os.path.join("doc", "doi", code) + ".pdf"
            file2 = os.path.join("doc", "doi", title) + ".pdf"
            if os.path.isfile(file):
                subprocess.Popen(['atril', file])
            elif os.path.isfile(file2):
                os.system(['atril', file2])
        elif item.parent():
            url = QtCore.QUrl("http://dx.doi.org/%s" % item.text(4))
            QtGui.QDesktopServices.openUrl(url)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    dialog = ShowReference()
    dialog.show()
    sys.exit(app.exec_())
