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


###############################################################################
# qsci_simple_pythoneditor.pyw
#
# QScintilla sample with PyQt
#
# Eli Bendersky (eliben@gmail.com)
# This code is in the public domain
#
# For now only integrated in pychemqt in meos tools for view custom method
# Other possible use to let user define custom functionality, define custom
# equipment, units, mEoS...
###############################################################################


import os
import sys

from tools.qt import QtGui, QtWidgets, Qsci


if os.environ["Qsci"] == "True":
    # With scintilla available use as python editor
    class SimplePythonEditor(Qsci.QsciScintilla):
        """Code editor for python code using Qscintilla"""
        ARROW_MARKER_NUM = 8

        def __init__(self, parent=None):
            super().__init__(parent)

            # For now set read-only property for all its uses,
            # possible modify when add same programatic functionality and use
            # this at code editor
            self.setReadOnly(True)

            # Set the default font
            font = QtGui.QFont()
            font.setFamily('Courier')
            font.setFixedPitch(True)
            font.setPointSize(10)
            self.setFont(font)
            self.setMarginsFont(font)

            # Margin 0 is used for line numbers
            fontmetrics = QtGui.QFontMetrics(font)
            self.setMarginsFont(font)
            self.setMarginWidth(0, fontmetrics.horizontalAdvance("000") + 6)
            self.setMarginLineNumbers(0, True)
            self.setMarginsBackgroundColor(QtGui.QColor("#cccccc"))

            # Clickable margin 1 for showing markers
            self.setMarginSensitivity(1, True)
            self.marginClicked.connect(self.on_margin_clicked)
            self.markerDefine(Qsci.QsciScintilla.MarkerSymbol.RightArrow, self.ARROW_MARKER_NUM)
            self.setMarkerBackgroundColor(
                QtGui.QColor("#ee1111"), self.ARROW_MARKER_NUM)

            # Brace matching: enable for a brace immediately before or after
            # the current position
            self.setBraceMatching(Qsci.QsciScintilla.BraceMatch.SloppyBraceMatch)

            # Current line visible with special background color
            self.setCaretLineVisible(True)
            self.setCaretLineBackgroundColor(QtGui.QColor("#ffe4e4"))

            # Set Python lexer
            # Set style for Python comments (style number 1) to a fixed-width
            # courier.
            lexer = Qsci.QsciLexerPython()
            lexer.setDefaultFont(font)
            self.setLexer(lexer)
            # self.SendScintilla(Qsci.QsciScintilla.SCI_STYLESETFONT, 1, 'Courier')
            # self.SendScintilla(Qsci.QsciScintilla.SCI_STYLESETFONT, 1)

            # Don't want to see the horizontal scrollbar at all
            # Use raw message to Scintilla here (all messages are documented
            # here: http://www.scintilla.org/ScintillaDoc.html)
            self.SendScintilla(Qsci.QsciScintilla.SCI_SETHSCROLLBAR, 0)

            # not too small
            # self.setMinimumSize(700, 450)

        def on_margin_clicked(self, nmargin, nline, modifiers):
            """Toggle marker for the line the margin was clicked on"""
            if self.markersAtLine(nline) != 0:
                self.markerDelete(nline, self.ARROW_MARKER_NUM)
            else:
                self.markerAdd(nline, self.ARROW_MARKER_NUM)

else:
    # Use a normal Qt widget
    class SimplePythonEditor(QtWidgets.QPlainTextEdit):
        """Simple text editor without scintilla dependence"""
        def __init__(self, *args):
            super().__init__(*args)
            self.setReadOnly(True)

        def setText(self, txt):
            """Define the same functionality as for scintilla api"""
            self.setPlainText(txt)


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    editor = SimplePythonEditor()
    editor.show()
    editor.setText(open(sys.argv[0]).read())
    app.exec()
