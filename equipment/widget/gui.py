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


from UI.widgets import SimpleStatus
from tools.qt import QtCore, QtWidgets


class ToolGui(QtWidgets.QWidget):
    """Parent widget with common functionality for equipment tools"""

    # Signal to propagate changes to parent widgets
    toggled = QtCore.pyqtSignal(bool)
    valueChanged = QtCore.pyqtSignal(object)

    def __init__(self, parent=None):
        super().__init__(parent)
        lyt = QtWidgets.QGridLayout(self)
        self.check = QtWidgets.QCheckBox(self.title)
        self.check.toggled.connect(self.setEnabled)
        lyt.addWidget(self.check, 0, 1, 1, 3)
        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 1, 1)

        self.loadUI()

        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 9, 4)

        self.msg = SimpleStatus()
        lyt.addWidget(self.msg, 10, 1, 1, 3)

        self.setEnabled(False)

    def setEnabled(self, boolean):
        """Toggled enable/disable state for all children widget except
        checkbox used to change this"""
        self.toggled.emit(boolean)
        for wdg in self.children():
            if wdg is not self.check and wdg is not self.layout():
                wdg.setEnabled(boolean)
        self.populate(self.Entity)

    def isChecked(self):
        """Level up check state as global value for widget"""
        return self.check.isChecked()

    def setChecked(self, boolean):
        self.check.setChecked(boolean)

    def loadUI(self):
        """Rewrite in child class to add widget"""

    def changeParams(self, key, value):
        """Change any kwargs value"""
        self.Entity(**{key: value})

    def populate(self, entity):
        if self.isChecked():
            self.msg.setState(entity)
        else:
            self.msg.clear()


class CallableEntity(QtCore.QObject):
    """Class with callable capability to add each input in a call"""
    def __init__(self, **kwargs):
        """Class constructor, copy kwargs for child class, it can be customize
        for child class to add functionality"""
        super().__init__()
        self.kw = self.__class__.kw.copy()

        self.__call__(**kwargs)

    def __call__(self, **kw):
        """Add callable functionality, so it can be possible add kwargs,
        advanced functionality can be added in subclass"""
        self._oldkw = self.kw.copy()

        self.kw.update(kw)

        if self.isCalculable and self._oldkw != self.kw:
            self.calculo()
        self.inputChanged.emit(self)
