#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Basic web browser to show help files
###############################################################################

import os
from PyQt5 import QtCore, QtGui, QtWidgets, QtWebKitWidgets



class HelpView(QtWidgets.QDialog):
    """HTML viewer to show help files"""
    def __init__(self, titulo, archivo, parent=None):
        super(HelpView, self).__init__(parent)
        self.setWindowState(QtCore.Qt.WindowMaximized)
        self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"]+"/images/button/help.png")))
        self.setWindowTitle(QtWidgets.QApplication.translate("pychemqt", "Help") +
                            " - " + titulo)
        layout = QtWidgets.QVBoxLayout(self)

        backAction = QtWidgets.QAction(QtGui.QIcon(
            os.environ["pychemqt"]+"/images/button/back.png"), "&Back", self)
        backAction.setShortcut(QtGui.QKeySequence.Back)
        homeAction = QtWidgets.QAction(QtGui.QIcon(
            os.environ["pychemqt"]+"/images/button/home.png"), "&Home", self)
        homeAction.setShortcut("Home")
        self.pageLabel = QtWidgets.QLabel()

        toolBar = QtWidgets.QToolBar()
        toolBar.addAction(backAction)
        toolBar.addAction(homeAction)
        toolBar.addWidget(self.pageLabel)
        layout.addWidget(toolBar)

        self.Viewer = QtWebKitWidgets.QWebView()
        self.Viewer.setUrl(archivo)
        layout.addWidget(self.Viewer)
        self.buttonBox = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Close)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox)

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    dialog = HelpView("Ciclón", QtCore.QUrl("help/index.html"))
    dialog.show()
    sys.exit(app.exec_())
