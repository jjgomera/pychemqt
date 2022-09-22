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
# Moody plot
###############################################################################


from configparser import ConfigParser
import os

from PyQt5 import QtGui, QtWidgets
from matplotlib import image

from lib.config import conf_dir, IMAGE_PATH
from lib.plot import mpl


class Chart(QtWidgets.QDialog):
    """Generic chart dialog"""

    def __init__(self, parent=None):
        super(Chart, self).__init__(parent)
        self.setWindowTitle(self.title)
        layout = QtWidgets.QGridLayout(self)
        layout.setColumnStretch(3, 1)
        self.plt = mpl(self)
        self.plt.fig.canvas.mpl_connect('button_press_event', self.click)
        layout.addWidget(self.plt, 2, 1, 1, 4)

        btBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Close, self)
        butonPNG = QtWidgets.QPushButton(QtGui.QIcon(
            os.path.join(IMAGE_PATH, "button", "image.png")),
            QtWidgets.QApplication.translate("pychemqt", "Save as PNG"))
        butonPNG.clicked.connect(self.plt.savePNG)
        butonConfig = QtWidgets.QPushButton(QtGui.QIcon(
            os.path.join(IMAGE_PATH, "button", "configure.png")),
            QtWidgets.QApplication.translate("pychemqt", "Configure"))
        butonConfig.clicked.connect(self.configure)
        butonCalculate = QtWidgets.QPushButton(QtGui.QIcon(
            os.path.join(IMAGE_PATH, "button", "calculator.png")),
            QtWidgets.QApplication.translate("pychemqt", "Calculate point"))
        butonCalculate.clicked.connect(self.calculate)
        btBox.rejected.connect(self.reject)
        btBox.layout().insertWidget(0, butonPNG)
        btBox.layout().insertWidget(0, butonCalculate)
        btBox.layout().insertWidget(0, butonConfig)
        layout.addWidget(btBox, 3, 1, 1, 4)

        self.Preferences = ConfigParser()
        self.Preferences.read(conf_dir+"pychemqtrc")
        self.config()
        self.plot()

        logo = image.imread(os.path.join(IMAGE_PATH, "pychemqt_98.png"))
        newax = self.plt.fig.add_axes(self.locLogo, anchor='SW')
        newax.imshow(logo, alpha=0.5)
        newax.axis('off')

    def configure(self):
        dlg = self.configDialog(self.Preferences)
        if dlg.exec_():
            self.Preferences = dlg.value(self.Preferences)
            self.Preferences.write(open(conf_dir+"pychemqtrc", "w"))
            self.plot()

    def config(self):
        """Initialization action in plot don't neccesary to update in plot"""
        pass

    def click(self, event):
        """Define functionality when the mouse click in plot"""
        pass

    def plot(self):
        """Plot procedure"""
        pass

    def calculate(self):
        """Define the functionality when click the calculate point button"""
        pass
