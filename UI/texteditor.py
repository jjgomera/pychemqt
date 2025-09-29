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
along with this program.  If not, see <http://www.gnu.org/licenses/>.


Module to define TextEditor widget, define common text functionality: \
fontcolor, fontstyle, fontweigh, alignment, bold, italic...
This widget is useful to show/edit note properties of stream, equipment... or \
whatever a rich text could be used
'''


import os
from tools.qt import QtCore, QtGui, QtWidgets
from UI.widgets import createAction


class TextEditor(QtWidgets.QWidget):
    """Text editor widget"""
    textChanged = QtCore.pyqtSignal(str, str)

    def __init__(self, text="", parent=None):
        """Constructor, opcional parameter text to set initial value"""
        super().__init__(parent)
        self.text = text
        self.setWindowTitle(self.tr("Notes"))
        gridLayout = QtWidgets.QVBoxLayout(self)

        toolbar = QtWidgets.QToolBar()
        toolbar.setIconSize(QtCore.QSize(16, 16))
        gridLayout.addWidget(toolbar)
        self.fontComboBox = QtWidgets.QFontComboBox()
        self.fontComboBox.setToolTip(self.tr("Font name"))
        self.fontComboBox.textActivated.connect(self.font)
        toolbar.addWidget(self.fontComboBox)
        self.fontColor = QtWidgets.QPushButton()
        self.fontColor.setFixedSize(22, 22)
        self.fontColor.setPalette(QtGui.QPalette(QtGui.QColor("black")))
        self.fontColor.setToolTip(self.tr("Font color"))
        self.fontColor.clicked.connect(self.colordialog)
        toolbar.addWidget(self.fontColor)
        self.fontSize = QtWidgets.QComboBox()
        for i in QtGui.QFontDatabase.standardSizes():
            self.fontSize.addItem(str(i))
        self.fontSize.setToolTip(self.tr("Font size"))
        self.fontSize.textActivated.connect(self.pointSize)
        toolbar.addWidget(self.fontSize)

        self.actionBold = createAction(
            icon=os.path.join("button", "format-text-bold.png"),
            text=self.tr("Bold"),
            slot=self.bold, checkable=True)
        toolbar.addAction(self.actionBold)
        self.actionItalic = createAction(
            icon=os.path.join("button", "format-text-italic.png"),
            text=self.tr("Italic"),
            slot=self.italic, checkable=True)
        toolbar.addAction(self.actionItalic)
        self.actionUnderline = createAction(
            icon=os.path.join("button", "format-text-underline.png"),
            text=self.tr("Underline"),
            slot=self.underline, checkable=True)
        toolbar.addAction(self.actionUnderline)
        self.actionStrikeout = createAction(
            icon=os.path.join("button", "format-text-strikethrough.png"),
            text=self.tr("Strike through"),
            slot=self.strikethrough, checkable=True)
        toolbar.addAction(self.actionStrikeout)
        self.actionSuperScript = createAction(
            icon=os.path.join("button", "font-superscript.png"),
            text=self.tr("Superscript"),
            slot=self.superscript, checkable=True)
        toolbar.addAction(self.actionSuperScript)
        self.actionSubScript = createAction(
            icon=os.path.join("button", "font-subscript.png"),
            text=self.tr("Subscript"),
            slot=self.subscript, checkable=True)
        toolbar.addAction(self.actionSubScript)
        toolbar.addSeparator()
        self.actionAlignLeft = createAction(
            icon=os.path.join("button", "format-justify-left.png"),
            text=self.tr("Align left"),
            slot=self.left, checkable=True)
        toolbar.addAction(self.actionAlignLeft)
        self.actionCenter = createAction(
            icon=os.path.join("button", "format-justify-center.png"),
            text=self.tr("Center"),
            slot=self.center, checkable=True)
        toolbar.addAction(self.actionCenter)
        self.actionJustify = createAction(
            icon=os.path.join("button", "format-justify-fill.png"),
            text=self.tr("Justify"),
            slot=self.justify, checkable=True)
        toolbar.addAction(self.actionJustify)
        self.actionAlignRight = createAction(
            icon=os.path.join("button", "format-justify-right.png"),
            text=self.tr("Align right"),
            slot=self.right, checkable=True)
        toolbar.addAction(self.actionAlignRight)

        self.notes = QtWidgets.QTextEdit(text)
        self.notes.setMinimumWidth(450)
        self.notes.textChanged.connect(self.setText)
        self.notes.cursorPositionChanged.connect(self.updateUI)
        gridLayout.addWidget(self.notes)

        group = QtGui.QActionGroup(self)
        group.addAction(self.actionAlignLeft)
        group.addAction(self.actionCenter)
        group.addAction(self.actionJustify)
        group.addAction(self.actionAlignRight)

        self.notes.setFocus()
        self.pointSize("10")
        self.fontSize.setCurrentIndex(self.fontSize.findText("10"))
        self.font(self.fontComboBox.currentFont().family())

    def setText(self, text=None):
        """Set text value"""
        if text:
            self.notes.setHtml(text)
        else:
            text = self.notes.toHtml()
            self.textChanged.emit(text, self.notes.toPlainText())
        self.text = text

    def MergeFormat(self, fmt):
        """Merge format to current text"""
        cursor = self.notes.textCursor()
        if not cursor.hasSelection():
            cursor.select(QtGui.QTextCursor.SelectionType.Document)
        cursor.mergeCharFormat(fmt)
        self.notes.mergeCurrentCharFormat(fmt)

    def font(self, family):
        """Change font of text"""
        fmt = QtGui.QTextCharFormat()
        fmt.setFontFamily(family)
        self.MergeFormat(fmt)

    def left(self):
        """Left align of text"""
        self.notes.setAlignment(QtCore.Qt.AlignmentFlag.AlignLeft
                                | QtCore.Qt.AlignmentFlag.AlignAbsolute)

    def center(self):
        """Center align of text"""
        self.notes.setAlignment(QtCore.Qt.AlignmentFlag.AlignHCenter)

    def right(self):
        """Right align of text"""
        self.notes.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight
                                | QtCore.Qt.AlignmentFlag.AlignAbsolute)

    def justify(self):
        """Justify text to end of line"""
        self.notes.setAlignment(QtCore.Qt.AlignmentFlag.AlignJustify)

    def pointSize(self, size):
        """Change point size of font"""
        points = int(size)
        if points> 0:
            fmt = QtGui.QTextCharFormat()
            fmt.setFontPointSize(points)
            self.MergeFormat(fmt)

    def bold(self):
        """Set bold property of font"""
        fmt = QtGui.QTextCharFormat()
        if self.actionBold.isChecked():
            fmt.setFontWeight(QtGui.QFont.Weight.Bold)
        else:
            fmt.setFontWeight(QtGui.QFont.Weight.Normal)
        self.MergeFormat(fmt)

    def italic(self):
        """Set italic property of font"""
        fmt = QtGui.QTextCharFormat()
        fmt.setFontItalic(self.actionItalic.isChecked())
        self.MergeFormat(fmt)

    def underline(self):
        """Set underline property of font"""
        fmt = QtGui.QTextCharFormat()
        fmt.setFontUnderline(self.actionUnderline.isChecked())
        self.MergeFormat(fmt)

    def strikethrough(self):
        """Set strikethrough property of font"""
        fmt = QtGui.QTextCharFormat()
        fmt.setFontStrikeOut(self.actionStrikeout.isChecked())
        self.MergeFormat(fmt)

    def superscript(self):
        """Set superscript property of font"""
        if self.actionSubScript.isChecked():
            self.actionSubScript.blockSignals(True)
            self.actionSubScript.setChecked(False)
            self.actionSubScript.blockSignals(False)
        fmt = QtGui.QTextCharFormat()
        if self.actionSuperScript.isChecked():
            fmt.setVerticalAlignment(
                QtGui.QTextCharFormat.VerticalAlignment.AlignSuperScript)
        else:
            fmt.setVerticalAlignment(
                QtGui.QTextCharFormat.VerticalAlignment.AlignNormal)
        self.MergeFormat(fmt)

    def subscript(self):
        """Set subscript property of font"""
        if self.actionSuperScript.isChecked():
            self.actionSuperScript.blockSignals(True)
            self.actionSuperScript.setChecked(False)
            self.actionSuperScript.blockSignals(False)
        fmt = QtGui.QTextCharFormat()
        if self.actionSubScript.isChecked():
            fmt.setVerticalAlignment(
                QtGui.QTextCharFormat.VerticalAlignment.AlignSubScript)
        else:
            fmt.setVerticalAlignment(
                QtGui.QTextCharFormat.VerticalAlignment.AlignNormal)
        self.MergeFormat(fmt)

    def updateUI(self):
        """Update button status when cursor move"""
        self.fontComboBox.setCurrentIndex(self.fontComboBox.findText(
            self.notes.fontFamily()))
        self.fontColor.setPalette(QtGui.QPalette(self.notes.textColor()))
        self.fontSize.setCurrentIndex(self.fontSize.findText(
            str(int(self.notes.fontPointSize()))))
        self.actionBold.setChecked(
            self.notes.fontWeight() == QtGui.QFont.Weight.Bold)
        self.actionItalic.setChecked(self.notes.fontItalic())
        self.actionUnderline.setChecked(self.notes.fontUnderline())
        fmt = self.notes.currentCharFormat()
        self.actionStrikeout.setChecked(fmt.fontStrikeOut())
        self.actionSuperScript.setChecked(False)
        self.actionSubScript.setChecked(False)
        if fmt.verticalAlignment() == \
                QtGui.QTextCharFormat.VerticalAlignment.AlignSuperScript:
            self.actionSuperScript.setChecked(True)
        elif fmt.verticalAlignment() == \
                QtGui.QTextCharFormat.VerticalAlignment.AlignSubScript:
            self.actionSubScript.setChecked(True)
        self.actionAlignLeft.setChecked(
            self.notes.alignment() == QtCore.Qt.AlignmentFlag.AlignLeft)
        self.actionCenter.setChecked(
            self.notes.alignment() == QtCore.Qt.AlignmentFlag.AlignHCenter)
        self.actionJustify.setChecked(
            self.notes.alignment() == QtCore.Qt.AlignmentFlag.AlignJustify)
        self.actionAlignRight.setChecked(
            self.notes.alignment() == QtCore.Qt.AlignmentFlag.AlignRight)

    def colordialog(self):
        """Show dialog to choose font color"""
        dialog = QtWidgets.QColorDialog(self.notes.textColor(), self)
        if dialog.exec():
            self.fontColor.setPalette(QtGui.QPalette(dialog.currentColor()))
            self.notes.setTextColor(dialog.currentColor())
            fmt = QtGui.QTextCharFormat()
            fmt.setForeground(dialog.currentColor())
            self.MergeFormat(fmt)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Form = TextEditor("Ejemplo")
    Form.show()
    sys.exit(app.exec())
