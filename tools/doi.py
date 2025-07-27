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


import os

from tools.qt import QtCore, QtGui, QtWidgets
from tools.pdf import openPDF

import lib
from lib.config import IMAGE_PATH
from equipment import equipments


class QLineEditClickable(QtWidgets.QLineEdit):
    """Custom QLineEdit to catch Enter key and set focus to list and avoid
    close dialog"""
    def keyPressEvent(self, event):
        """Rewrite keyPressEvent to setFocus when Enter key is pressed"""
        if event.key() == QtCore.Qt.Key.Key_Return:
            self.parent().parent().tree.setFocus()
        else:
            super().keyPressEvent(event)


class QTreeWidgetClickable(QtWidgets.QTreeWidget):
    """Custom QTreeWidget to catch Enter key and open file with it"""
    def keyPressEvent(self, event):
        """Rewrite keyPressEvent to expand or open item"""
        if event.key() == QtCore.Qt.Key.Key_Return:
            if self.currentItem().childCount():
                self.currentItem().setExpanded(True)
            else:
                self.parent().open(self.currentItem())
        else:
            super().keyPressEvent(event)


class ShowReference(QtWidgets.QDialog):
    """Dialog to show the references used in the program"""

    searchIndex = -1
    searchResults = []

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"]+"/images/button/help.png")))
        self.setWindowTitle(self.tr("Reference Paper Show Dialog"))
        layout = QtWidgets.QGridLayout(self)

        self.tree = QTreeWidgetClickable()
        header = QtWidgets.QTreeWidgetItem(
            ["id", self.tr("Autor"), self.tr("Title"),
             self.tr("Reference"), self.tr("doi")])
        self.tree.setHeaderItem(header)
        layout.addWidget(self.tree, 1, 1, 2, 2)

        self.searchWidget = QtWidgets.QWidget()
        self.searchWidget.setMaximumSize(200, 25)
        searchlayout = QtWidgets.QHBoxLayout(self.searchWidget)
        searchlayout.setSpacing(0)
        searchlayout.setContentsMargins(0, 0, 0, 0)
        self.searchTxt = QLineEditClickable()
        self.searchTxt.textChanged.connect(self.search)
        searchlayout.addWidget(self.searchTxt)
        self.btnPrevious = QtWidgets.QToolButton()
        self.btnPrevious.setIcon(QtGui.QIcon(QtGui.QPixmap(
            os.path.join(IMAGE_PATH, "button", "arrow-left.png"))))
        self.btnPrevious.setEnabled(False)
        self.btnPrevious.clicked.connect(self.searchPrevious)
        searchlayout.addWidget(self.btnPrevious)
        self.btnNext = QtWidgets.QToolButton()
        self.btnNext.setIcon(QtGui.QIcon(QtGui.QPixmap(
            os.path.join(IMAGE_PATH, "button", "arrow-right.png"))))
        self.btnNext.setEnabled(False)
        self.btnNext.clicked.connect(self.searchNext)
        searchlayout.addWidget(self.btnNext)
        self.searchWidget.hide()
        layout.addWidget(self.searchWidget, 3, 1)

        bttBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Close)
        bttBox.rejected.connect(self.reject)
        layout.addWidget(bttBox, 3, 2)

        shSearch = QtGui.QShortcut(QtGui.QKeySequence.StandardKey.Find, self)
        shSearch.activated.connect(self.enableSearch)
        shHide = QtGui.QShortcut(QtGui.QKeySequence.StandardKey.Cancel, self)
        shHide.activated.connect(self.disableSearch)
        self.shNext = QtGui.QShortcut(
            QtGui.QKeySequence.StandardKey.FindNext, self)
        self.shNext.activated.connect(self.searchNext)
        self.shNext.setEnabled(False)
        self.shPrevious = QtGui.QShortcut(
            QtGui.QKeySequence.StandardKey.FindPrevious, self)
        self.shPrevious.activated.connect(self.searchPrevious)
        self.shPrevious.setEnabled(False)

        self.fill()
        self.tree.sortItems(0, QtCore.Qt.SortOrder.AscendingOrder)
        self.tree.itemDoubleClicked.connect(self.open)
        self.tree.setColumnWidth(0, 200)
        self.tree.setColumnWidth(1, 200)
        self.tree.setColumnWidth(2, 200)
        self.tree.setColumnWidth(3, 200)

    # Search functionality
    def enableSearch(self):
        """Show search widgets"""
        self.searchWidget.show()
        self.searchTxt.setFocus(QtCore.Qt.FocusReason.ShortcutFocusReason)

    def disableSearch(self):
        """Hide search widgets"""
        self.searchWidget.hide()
        self.searchIndex = -1
        self.searchResults = []
        self.shNext.setEnabled(False)
        self.shPrevious.setEnabled(False)

    def search(self, txt):
        """Search txt in tree widget contents"""
        self.searchIndex = -1
        self.searchResults = []
        for col in range(4):
            flags = (QtCore.Qt.MatchFlag.MatchContains
                     | QtCore.Qt.MatchFlag.MatchRecursive)
            self.searchResults += self.tree.findItems(txt, flags, col)

        # Enable navitation in search results if search if successful
        if self.searchResults:
            self.searchNext()

        if len(self.searchResults) > 1:
            self.shNext.setEnabled(True)
            self.shPrevious.setEnabled(True)
            self.searchNext()
            self.btnPrevious.setEnabled(True)
            self.btnNext.setEnabled(True)
        else:
            self.shNext.setEnabled(False)
            self.shPrevious.setEnabled(False)
            self.btnPrevious.setEnabled(False)
            self.btnNext.setEnabled(False)

    def searchPrevious(self):
        """Set previous result in search"""
        self.searchIndex -= 1
        if self.searchIndex < 0:
            self.searchIndex = len(self.searchResults)-1
        self.tree.setCurrentItem(self.searchResults[self.searchIndex])

    def searchNext(self):
        """Set next result in search"""
        self.searchIndex += 1
        if self.searchIndex >= len(self.searchResults):
            self.searchIndex = 0
        self.tree.setCurrentItem(self.searchResults[self.searchIndex])
        self.tree.scrollToItem(
            self.searchResults[self.searchIndex],
            QtWidgets.QAbstractItemView.ScrollHint.PositionAtCenter)

    def fill(self):
        """Fill tree with documentation entries"""
        # General population Library
        for library in lib.__all__:
            __import__(f"lib.{library}")
            module = lib.__getattribute__(library)
            if hasattr(module, "__doi__") and module.__doi__:
                itemModule = QtWidgets.QTreeWidgetItem(["lib/"+library])
                self.tree.addTopLevelItem(itemModule)
                for key in sorted(module.__doi__.keys()):

                    if key == "lib.EoS.Cubic":
                        itemSubModule = QtWidgets.QTreeWidgetItem([key])
                        itemModule.addChildren([itemSubModule])
                        for key2 in sorted(module.__doi__[key].keys()):
                            itemSubModule2 = QtWidgets.QTreeWidgetItem([key2])
                            itemSubModule.addChildren([itemSubModule2])
                            for link in module.__doi__[key][key2]:
                                item = QtWidgets.QTreeWidgetItem([
                                    "", link["autor"], link["title"],
                                    link["ref"], link["doi"]])
                                if findFile(link):
                                    icon = QtGui.QIcon(QtGui.QPixmap(
                                        os.path.join(
                                            IMAGE_PATH, "button", "ok.png")))
                                    item.setIcon(0, icon)
                                itemSubModule2.addChild(item)

                    # Special case for submodules
                    elif library in ["EoS", "mEoS", "newComponent"]:
                        itemSubModule = QtWidgets.QTreeWidgetItem([key])
                        itemModule.addChildren([itemSubModule])
                        for key2 in sorted(module.__doi__[key].keys()):
                            link = module.__doi__[key][key2]
                            if key == "lib.EoS.Cubic":
                                title = key2
                            elif library == "EoS":
                                title = ""
                            elif library == "newComponent":
                                title = ""
                            else:
                                title = key2.replace("_", "")
                            item = QtWidgets.QTreeWidgetItem([
                                title, link["autor"],
                                link["title"], link["ref"], link["doi"]])

                            if findFile(link):
                                icon = QtGui.QIcon(QtGui.QPixmap(os.path.join(
                                    IMAGE_PATH, "button", "ok.png")))
                                item.setIcon(0, icon)
                            itemSubModule.addChild(item)

                    else:
                        link = module.__doi__[key]
                        if isinstance(key, int):
                            header = ""
                        else:
                            header = key
                        item = QtWidgets.QTreeWidgetItem([
                            header, link["autor"], link["title"], link["ref"],
                            link["doi"]])

                        if findFile(link):
                            icon = QtGui.QIcon(QtGui.QPixmap(os.path.join(
                                IMAGE_PATH, "button", "ok.png")))
                            item.setIcon(0, icon)
                        itemModule.addChild(item)

                # Show item sorted by autor in libraries
                if "EoS" not in library:
                    itemModule.sortChildren(
                        1, QtCore.Qt.SortOrder.AscendingOrder)

        # Equipment
        itemEquipment = QtWidgets.QTreeWidgetItem([self.tr("Equipments")])
        self.tree.addTopLevelItem(itemEquipment)
        for equip in equipments:
            itemequip = QtWidgets.QTreeWidgetItem([equip.__name__])
            itemEquipment.addChild(itemequip)
            for link in equip.__doi__:
                item = QtWidgets.QTreeWidgetItem([
                    "", link["autor"], link["title"],
                    link["ref"], link["doi"]])

                if findFile(link):
                    icon = QtGui.QIcon(QtGui.QPixmap(os.path.join(
                        IMAGE_PATH, "button", "ok.png")))
                    item.setIcon(0, icon)

                itemequip.addChild(item)

    def open(self, item):
        """Open file if exist in doc/doi folder or open a browser with link"""
        if item.parent() and not item.icon(0).isNull():
            title = item.text(2)
            text = item.text(4)
            code = str(text).replace("/", "_")
            file = os.path.join("doc", code) + ".pdf"
            file2 = os.path.join("doc", title) + ".pdf"
            if os.path.isfile(file):
                openPDF(file, title)
            elif os.path.isfile(file2):
                openPDF(file2, title)
        elif item.parent():
            url = QtCore.QUrl(f"http://dx.doi.org/{item.text(4)}")
            QtGui.QDesktopServices.openUrl(url)


def findFile(ref):
    """Search reference paper path in documentation forder and return boolean
    if it's available"""
    code = ref["doi"].replace("/", "_")
    file = os.path.join("doc", code) + ".pdf"
    file2 = os.path.join("doc", ref["title"]) + ".pdf"
    return os.path.isfile(file) or os.path.isfile(file2)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    dialog = ShowReference()
    dialog.show()
    sys.exit(app.exec())
