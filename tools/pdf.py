#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2023, Juan José Gómez Romera <jjgomera@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Library to pdf viewer functionality

    * :class:`PDFWidget`: Internal viewer using a QWebEngineView
    * :func:`openPDF`: Show a viewer with a pdf file
'''


import os
import subprocess

from lib.config import Preferences
from tools.qt import QtCore, QtWidgets, QtWebEngineWidgets


class PDFWidget(QtWidgets.QDialog):
    """Internal viewer using a QtWebEngineView"""
    def __init__(self, pdffile, title=None):
        super().__init__()

        if title:
            self.setWindowTitle(title)
        lyt = QtWidgets.QVBoxLayout(self)

        self.webView = QtWebEngineWidgets.QWebEngineView()
        lyt.addWidget(self.webView)
        self.webView.settings().setAttribute(
            self.webView.settings().WebAttribute.PluginsEnabled, True)
        self.webView.settings().setAttribute(
            self.webView.settings().WebAttribute.PdfViewerEnabled, True)
        wd = os.environ["pychemqt"]
        self.webView.setUrl(QtCore.QUrl(f"file://{wd}/{pdffile}"))
        self.showMaximized()


def openPDF(file, title=None):
    """Show a viewer with a pdf file

    Parameters
    ----------
    file : str
        Path for file to view
    title : str, optional
        Optional title for the window used only when use the internal viewer
    """

    if Preferences.getboolean("Applications", "PDF"):
        app = Preferences.get("Applications", 'PDFExternal')
        subprocess.Popen([app, file])
    else:
        winPDF = PDFWidget(file, title)
        winPDF.exec()
