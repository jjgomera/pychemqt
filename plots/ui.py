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
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


###############################################################################
# Common graphycal functionality for plots

#     * :class:`Chart`: Generic chart dialog
#     * :class:`ConfigDialog`: Dialog tool for standalone use
###############################################################################


import os

from matplotlib import image

from lib.config import conf_dir, IMAGE_PATH, Preferences
from lib.plot import PlotWidget
from tools.qt import QtGui, QtWidgets, tr


class Chart(QtWidgets.QDialog):
    """Generic chart dialog"""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(self.title)
        layout = QtWidgets.QGridLayout(self)
        layout.setColumnStretch(3, 1)
        self.plotWidget = QtWidgets.QWidget(self)
        lyt = QtWidgets.QGridLayout(self.plotWidget)
        self.plt = PlotWidget(parent=self)
        self.plt.fig.canvas.mpl_connect('button_press_event', self.click)
        lyt.addWidget(self.plt, 1, 1)
        layout.addWidget(self.plotWidget, 2, 1, 1, 4)

        btBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Close, self)
        self.butonPNG = QtWidgets.QPushButton(QtGui.QIcon(
            os.path.join(IMAGE_PATH, "button", "image.png")),
            tr("pychemqt", "Save as PNG"))
        self.butonPNG.clicked.connect(self.plt.savePNG)
        self.butonConf = QtWidgets.QPushButton(QtGui.QIcon(
            os.path.join(IMAGE_PATH, "button", "configure.png")),
            tr("pychemqt", "Configure"))
        self.butonConf.clicked.connect(self.configure)
        self.butonCalc = QtWidgets.QPushButton(QtGui.QIcon(
            os.path.join(IMAGE_PATH, "button", "calculator.png")),
            tr("pychemqt", "Calculate point"))
        self.butonCalc.clicked.connect(self.calculate)
        btBox.rejected.connect(self.reject)
        btBox.addButton(self.butonConf,
                        QtWidgets.QDialogButtonBox.ButtonRole.ResetRole)
        btBox.addButton(self.butonCalc,
                        QtWidgets.QDialogButtonBox.ButtonRole.ResetRole)
        btBox.addButton(self.butonPNG,
                        QtWidgets.QDialogButtonBox.ButtonRole.ResetRole)
        layout.addWidget(btBox, 3, 1, 1, 4)

        self.customUI()

        self.config()
        self.plot()
        self.set_logo(self.plt)

    def set_logo(self, plot):
        """Draw pychemqt logo as watermark in plot"""
        logo = image.imread(os.path.join(IMAGE_PATH, "pychemqt_98.png"))
        newax = plot.fig.add_axes(self.locLogo, anchor='SW')
        newax.imshow(logo, alpha=0.5)
        newax.axis('off')

    def configure(self):
        """Show configure dialog and save changes"""
        dlg = ConfigDialog(self.widgetConfig, Preferences)
        if dlg.exec():
            Pref = dlg.value(Preferences)
            with open(conf_dir+"pychemqtrc", "w", encoding="utf-8") as file:
                Pref.write(file)
            self.plot()

    def customUI(self):
        """Define custom UI element"""

    def config(self):
        """Initialization action in plot don't neccesary to update in plot"""

    def click(self, event):
        """Define functionality when the mouse click in plot"""

    def plot(self):
        """Plot procedure"""

    def calculate(self):
        """Define the functionality when click the calculate point button"""


class ConfigDialog(QtWidgets.QDialog):
    """Dialog to configure a chart"""
    def __init__(self, widget, config=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle(widget.TITLECONFIG)
        layout = QtWidgets.QVBoxLayout(self)
        self.widget = widget(config)
        layout.addWidget(self.widget)
        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Cancel
            | QtWidgets.QDialogButtonBox.StandardButton.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox)

    def value(self, config):
        """Function result for wizard"""
        config = self.widget.value(config)
        return config
