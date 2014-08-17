#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Module to define TextEditor widget, define common text functionality:
# fontcolor, fontstyle, fontweigh, alignment, bold, italic...
# This widget is useful to show/edit note properties of stream, equipment...
###############################################################################

import os

from PyQt4 import QtCore, QtGui

from widgets import createAction


class TextEditor(QtGui.QWidget):
    """Text editor widget"""
    textChanged = QtCore.pyqtSignal(str, str)

    def __init__(self, texto="", parent=None):
        """Constructor, opcional parameter texto to set initial value"""
        super(TextEditor, self).__init__(parent)
        self.texto = texto
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Notes"))
        gridLayout = QtGui.QVBoxLayout(self)

        toolbar = QtGui.QToolBar()
        toolbar.setIconSize(QtCore.QSize(16, 16))
        gridLayout.addWidget(toolbar)
        self.fontComboBox = QtGui.QFontComboBox()
        self.fontComboBox.setToolTip(
            QtGui.QApplication.translate("pychemqt", "Font name"))
        self.fontComboBox.activated[str].connect(self.font)
        toolbar.addWidget(self.fontComboBox)
        self.FontColor = QtGui.QPushButton()
        self.FontColor.setFixedSize(22, 22)
        self.FontColor.setPalette(QtGui.QPalette(QtGui.QColor("black")))
        self.FontColor.setToolTip(
            QtGui.QApplication.translate("pychemqt", "Font color"))
        self.FontColor.clicked.connect(self.colordialog)
        toolbar.addWidget(self.FontColor)
        self.FontSize = QtGui.QComboBox()
        for i in QtGui.QFontDatabase.standardSizes():
            self.FontSize.addItem(str(i))
        self.FontSize.setToolTip(
            QtGui.QApplication.translate("pychemqt", "Font size"))
        self.FontSize.activated[str].connect(self.PointSize)
        toolbar.addWidget(self.FontSize)

        self.actionNegrita = createAction(
            icon=os.environ["pychemqt"]+"/images/button/format-text-bold.png",
            text=QtGui.QApplication.translate("pychemqt", "Bold"),
            slot=self.Negrita, checkable=True)
        toolbar.addAction(self.actionNegrita)
        self.actionCursiva = createAction(
            icon=os.environ["pychemqt"]+"/images/button/format-text-italic.png",
            text=QtGui.QApplication.translate("pychemqt", "Italic"),
            slot=self.Cursiva, checkable=True)
        toolbar.addAction(self.actionCursiva)
        self.actionSubrayado = createAction(
            icon=os.environ["pychemqt"]+"/images/button/format-text-underline.png",
            text=QtGui.QApplication.translate("pychemqt", "Underline"),
            slot=self.Subrayado, checkable=True)
        toolbar.addAction(self.actionSubrayado)
        self.actionTachado = createAction(
            icon=os.environ["pychemqt"]+"/images/button/format-text-strikethrough.png",
            text=QtGui.QApplication.translate("pychemqt", "Strike through"),
            slot=self.Tachado, checkable=True)
        toolbar.addAction(self.actionTachado)
        self.actionSuperScript = createAction(
            icon=os.environ["pychemqt"]+"/images/button/font-superscript.png",
            text=QtGui.QApplication.translate("pychemqt", "Superscript"),
            slot=self.Superindice, checkable=True)
        toolbar.addAction(self.actionSuperScript)
        self.actionSubScript = createAction(
            icon=os.environ["pychemqt"]+"/images/button/font-subscript.png",
            text=QtGui.QApplication.translate("pychemqt", "Subscript"),
            slot=self.Subindice, checkable=True)
        toolbar.addAction(self.actionSubScript)
        toolbar.addSeparator()
        self.actionAlinearIzquierda = createAction(
            icon=os.environ["pychemqt"]+"/images/button/format-justify-left.png",
            text=QtGui.QApplication.translate("pychemqt", "Align left"),
            slot=self.izquierda, checkable=True)
        toolbar.addAction(self.actionAlinearIzquierda)
        self.actionCentrar = createAction(
            icon=os.environ["pychemqt"]+"/images/button/format-justify-center.png",
            text=QtGui.QApplication.translate("pychemqt", "Center"),
            slot=self.centrar, checkable=True)
        toolbar.addAction(self.actionCentrar)
        self.actionJustificar = createAction(
            icon=os.environ["pychemqt"]+"/images/button/format-justify-fill.png",
            text=QtGui.QApplication.translate("pychemqt", "Justify"),
            slot=self.justificar, checkable=True)
        toolbar.addAction(self.actionJustificar)
        self.actionAlinearDerecha = createAction(
            icon=os.environ["pychemqt"]+"/images/button/format-justify-right.png",
            text=QtGui.QApplication.translate("pychemqt", "Align right"),
            slot=self.derecha, checkable=True)
        toolbar.addAction(self.actionAlinearDerecha)

        self.notas = QtGui.QTextEdit(texto)
        self.notas.setMinimumWidth(450)
        self.notas.textChanged.connect(self.setText)
        self.notas.cursorPositionChanged.connect(self.updateUI)
        gridLayout.addWidget(self.notas)

        group = QtGui.QActionGroup(self)
        group.addAction(self.actionAlinearIzquierda)
        group.addAction(self.actionCentrar)
        group.addAction(self.actionJustificar)
        group.addAction(self.actionAlinearDerecha)

        self.notas.setFocus()
        self.PointSize("10")
        self.FontSize.setCurrentIndex(self.FontSize.findText("10"))
        self.font(self.fontComboBox.currentFont().family())

    def setText(self, texto=None):
        """Set text value"""
        if texto:
            self.notas.setHtml(texto)
        else:
            texto = self.notas.toHtml()
            self.textChanged.emit(texto, self.notas.toPlainText())
        self.texto = texto

    def MergeFormat(self, format):
        """Merge format to current text"""
        cursor = self.notas.textCursor()
        if not cursor.hasSelection():
            cursor.select(QtGui.QTextCursor.Document)
        cursor.mergeCharFormat(format)
        self.notas.mergeCurrentCharFormat(format)

    def font(self, family):
        format = QtGui.QTextCharFormat()
        format.setFontFamily(family)
        self.MergeFormat(format)

    def izquierda(self):
        self.notas.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignAbsolute)

    def centrar(self):
        self.notas.setAlignment(QtCore.Qt.AlignHCenter)

    def derecha(self):
        self.notas.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignAbsolute)

    def justificar(self):
        self.notas.setAlignment(QtCore.Qt.AlignJustify)

    def PointSize(self, size):
        puntos = int(size)
        if puntos > 0:
            format = QtGui.QTextCharFormat()
            format.setFontPointSize(puntos)
            self.MergeFormat(format)

    def Negrita(self):
        format = QtGui.QTextCharFormat()
        if self.actionNegrita.isChecked():
            format.setFontWeight(QtGui.QFont.Bold)
        else:
            format.setFontWeight(QtGui.QFont.Normal)
        self.MergeFormat(format)

    def Cursiva(self):
        format = QtGui.QTextCharFormat()
        format.setFontItalic(self.actionCursiva.isChecked())
        self.MergeFormat(format)

    def Subrayado(self):
        format = QtGui.QTextCharFormat()
        format.setFontUnderline(self.actionSubrayado.isChecked())
        self.MergeFormat(format)

    def Tachado(self):
        format = QtGui.QTextCharFormat()
        format.setFontStrikeOut(self.actionTachado.isChecked())
        self.MergeFormat(format)

    def Superindice(self):
        if self.actionSubScript.isChecked():
            self.actionSubScript.blockSignals(True)
            self.actionSubScript.setChecked(False)
            self.actionSubScript.blockSignals(False)
        format = QtGui.QTextCharFormat()
        if self.actionSuperScript.isChecked():
            format.setVerticalAlignment(QtGui.QTextCharFormat.AlignSuperScript)
        else:
            format.setVerticalAlignment(QtGui.QTextCharFormat.AlignNormal)
        self.MergeFormat(format)

    def Subindice(self):
        if self.actionSuperScript.isChecked():
            self.actionSuperScript.blockSignals(True)
            self.actionSuperScript.setChecked(False)
            self.actionSuperScript.blockSignals(False)
        format = QtGui.QTextCharFormat()
        if self.actionSubScript.isChecked():
            format.setVerticalAlignment(QtGui.QTextCharFormat.AlignSubScript)
        else:
            format.setVerticalAlignment(QtGui.QTextCharFormat.AlignNormal)
        self.MergeFormat(format)

    def updateUI(self):
        """Update button status when cursor move"""
        self.fontComboBox.setCurrentIndex(self.fontComboBox.findText(
            self.notas.fontFamily()))
        self.FontColor.setPalette(QtGui.QPalette(self.notas.textColor()))
        self.FontSize.setCurrentIndex(self.FontSize.findText(
            QtCore.QString.number(self.notas.fontPointSize())))
        self.actionNegrita.setChecked(
            self.notas.fontWeight() == QtGui.QFont.Bold)
        self.actionCursiva.setChecked(self.notas.fontItalic())
        self.actionSubrayado.setChecked(self.notas.fontUnderline())
        format = self.notas.currentCharFormat()
        self.actionTachado.setChecked(format.fontStrikeOut())
        self.actionSuperScript.setChecked(False)
        self.actionSubScript.setChecked(False)
        if format.verticalAlignment() == QtGui.QTextCharFormat.AlignSuperScript:
            self.actionSuperScript.setChecked(True)
        elif format.verticalAlignment() == QtGui.QTextCharFormat.AlignSubScript:
            self.actionSubScript.setChecked(True)
        self.actionAlinearIzquierda.setChecked(
            self.notas.alignment() == QtCore.Qt.AlignLeft)
        self.actionCentrar.setChecked(
            self.notas.alignment() == QtCore.Qt.AlignHCenter)
        self.actionJustificar.setChecked(
            self.notas.alignment() == QtCore.Qt.AlignJustify)
        self.actionAlinearDerecha.setChecked(
            self.notas.alignment() == QtCore.Qt.AlignRight)

    def colordialog(self):
        """Show dialog to choose font color"""
        dialog = QtGui.QColorDialog(self.notas.textColor(), self)
        if dialog.exec_():
            self.FontColor.setPalette(QtGui.QPalette(dialog.currentColor()))
            self.notas.setTextColor(dialog.currentColor())
            format = QtGui.QTextCharFormat()
            format.setForeground(dialog.currentColor())
            self.MergeFormat(format)

if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    Form = TextEditor("Ejemplo")
    Form.show()
    sys.exit(app.exec_())
