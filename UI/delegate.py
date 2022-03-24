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
# Module to implement delegate special editiing in tables
#   -CellEditor
#   -SpinEditor
#   -CheckEditor
#   -SpinEditor
###############################################################################


from PyQt5 import QtCore, QtGui, QtWidgets


class CellEditor(QtWidgets.QItemDelegate):
    """Numeric editor of tableitem, with numeric validator"""
    def __init__(self, parent=None):
        super(CellEditor, self).__init__(parent)

    def createEditor(self, parent, option, index):
        widget = QtWidgets.QLineEdit(parent)
        widget.setAlignment(QtCore.Qt.AlignRight)
        validator = QtGui.QDoubleValidator(self)
        locale = QtCore.QLocale("en")
        validator.setLocale(locale)
        widget.setValidator(validator)
        return widget


class SpinEditor(QtWidgets.QItemDelegate):
    """Spinbox editor for tableitem"""
    def __init__(self, parent=None):
        super(SpinEditor, self).__init__(parent)

    def createEditor(self, parent, option, index):
        widget = QtWidgets.QSpinBox(parent)
        widget.setAlignment(QtCore.Qt.AlignRight)
        widget.setMinimum(1)
        return widget


class CheckEditor(QtWidgets.QItemDelegate):
    """Checkbox editor for tableitem"""
    def __init__(self, parent=None):
        super(CheckEditor, self).__init__(parent)

    def createEditor(self, parent, option, index):
        widget = QtWidgets.QCheckBox(parent)
        return widget

    def setEditorData(self, editor, index):
        value = bool(index.data(QtCore.Qt.DisplayRole))
        editor.setChecked(value)

    def setModelData(self, editor, model, index):
        value = editor.isChecked()
        model.setData(index, QtCore.QVariant(value), QtCore.Qt.DisplayRole)


class ComboEditor(QtWidgets.QItemDelegate):
    """Combobox Editor for tableitem"""
    def __init__(self, owner, items=None):
        super(ComboEditor, self).__init__(owner)
        self.setItems(items)

    def setItems(self, items):
        self.items = items

    def createEditor(self, parent, option, index):
        self.editor = QtWidgets.QComboBox(parent)
        self.editor.addItems(self.items)
        return self.editor

    def setEditorData(self, editor, index):
        value = str(index.data(QtCore.Qt.DisplayRole).toString())
        try:
            num = self.items.index(value)
        except ValueError:
            num = -1
        editor.setCurrentIndex(num)

    def setModelData(self, editor, model, index):
        value = editor.currentText()
        model.setData(index, QtCore.QVariant(value), QtCore.Qt.DisplayRole)
