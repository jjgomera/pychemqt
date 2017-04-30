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


###############################################################################
# Library to configure the moody chart
#
#   - Widget: Moody chart configuration
#   - Dialog: Dialog tool for standalone use
###############################################################################


from PyQt5 import QtWidgets

from lib.friction import f_list
from UI.widgets import LineConfig


class Widget(QtWidgets.QWidget):
    """Moody chart configuration"""
    def __init__(self, config=None, parent=None):
        super(Widget, self).__init__(parent)
        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Method:")), 1, 1)
        self.metodos = QtWidgets.QComboBox()
        for f in f_list:
            line = f.__doc__.split("\n")[1]
            year = line.split(" ")[-1]
            name = line.split(" ")[-3]
            doc = name + " " + year
            self.metodos.addItem(doc)
        layout.addWidget(self.metodos, 1, 2)
        self.fanning = QtWidgets.QCheckBox(QtWidgets.QApplication.translate(
            "pychemqt", "Calculate fanning friction factor"))
        layout.addWidget(self.fanning, 2, 1, 1, 2)

        layout.addWidget(QtWidgets.QLabel("ε/d:"), 3, 1)
        self.ed = QtWidgets.QLineEdit()
        layout.addWidget(self.ed, 3, 2)
        self.lineconfig = LineConfig(
            "line", QtWidgets.QApplication.translate(
                "pychemqt", "Relative roughtness style line"))
        layout.addWidget(self.lineconfig, 4, 1, 1, 2)
        self.cruxconfig = LineConfig(
            "crux", QtWidgets.QApplication.translate(
                "pychemqt", "Crux style line"))
        layout.addWidget(self.cruxconfig, 5, 1, 1, 2)

        layout.addItem(QtWidgets.QSpacerItem(
            10, 0, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 10, 1, 1, 3)

        if config and config.has_section("Moody"):
            self.metodos.setCurrentIndex(config.getint("Moody", 'method'))
            self.fanning.setChecked(config.getboolean("Moody", 'fanning'))
            self.ed.setText(config.get("Moody", "ed"))
            self.lineconfig.setConfig(config, "Moody")
            self.cruxconfig.setConfig(config, "Moody")

    def value(self, config):
        if not config.has_section("Moody"):
            config.add_section("Moody")
        config.set("Moody", "method", str(self.metodos.currentIndex()))
        config.set("Moody", "fanning", str(self.fanning.isChecked()))
        config.set("Moody", "ed", self.ed.text())
        config = self.lineconfig.value(config, "Moody")
        config = self.cruxconfig.value(config, "Moody")
        return config


class Dialog(QtWidgets.QDialog):
    """Dialog to config thermal method calculations"""
    def __init__(self, config=None, parent=None):
        super(Dialog, self).__init__(parent)
        self.setWindowTitle(QtWidgets.QApplication.translate(
            "pychemqt", "Moody diagram configuration"))
        layout = QtWidgets.QVBoxLayout(self)
        self.widget = Widget(config)
        layout.addWidget(self.widget)
        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Cancel | QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox)

    def value(self, config):
        """Function result for wizard"""
        config = self.widget.value(config)
        return config


if __name__ == "__main__":
    import os
    import sys
    from configparser import ConfigParser
    app = QtWidgets.QApplication(sys.argv)

    conf_dir = os.path.expanduser('~') + "/.pychemqt/"
    config = ConfigParser()
    config.read(conf_dir+"pychemqtrc")

    Dialog = Dialog(config)
    Dialog.show()
    sys.exit(app.exec_())
