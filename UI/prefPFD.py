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


Library to configure process flow diagram (PFD)

    * :class:`Widget`: PFD configuration
    * :class:`Dialog`: Dialog tool for standalone use
    * :class:`ConfLine`: Composite widget with line format configuration tools
    * :class:`ConfLineDialog`: ConfLine dialog for standalone use
'''


import os

from UI.widgets import ColorSelector, Entrada_con_unidades
from UI.widgets import PFDLineCombo
from tools import UI_confResolution
from tools.qt import QtCore, QtGui, QtWidgets


class Widget(QtWidgets.QDialog):
    """Flow Diagram configuration"""
    def __init__(self, config=None, parent=None):
        super().__init__(parent)

        lyt = QtWidgets.QGridLayout(self)
        lyt.setContentsMargins(0, 0, 0, 0)
        scroll = QtWidgets.QScrollArea()
        scroll.setFrameStyle(QtWidgets.QFrame.Shape.NoFrame)
        lyt.addWidget(scroll)
        dlg = QtWidgets.QWidget()
        layout = QtWidgets.QGridLayout(dlg)

        layout.addWidget(QtWidgets.QLabel(self.tr("Brush color")), 1, 1)
        self.ColorButtonBrush = ColorSelector()
        layout.addWidget(self.ColorButtonBrush, 1, 2)
        layout.addWidget(QtWidgets.QLabel(self.tr("Brush style")), 2, 1)
        self.brush = BrushCombo()
        layout.addWidget(self.brush, 2, 2)
        layout.addWidget(QtWidgets.QLabel(self.tr("Input color")), 4, 1)
        self.ColorButtonEntrada = ColorSelector()
        layout.addWidget(self.ColorButtonEntrada, 4, 2)
        layout.addWidget(QtWidgets.QLabel(self.tr("Output color:")), 5, 1)
        self.ColorButtonSalida = ColorSelector()
        layout.addWidget(self.ColorButtonSalida, 5, 2)

        group = QtWidgets.QGroupBox(self.tr("Line format"))
        layout.addWidget(group, 6, 1, 1, 3)
        lyt = QtWidgets.QHBoxLayout(group)
        self.lineFormat = ConfLine()
        lyt.addWidget(self.lineFormat)

        group = QtWidgets.QGroupBox(self.tr("PFD resolution"))
        layout.addWidget(group, 7, 1, 1, 3)
        lyt = QtWidgets.QHBoxLayout(group)
        self.resolution = UI_confResolution.UI_confResolution_widget(config)
        lyt.addWidget(self.resolution)

        layout.addItem(QtWidgets.QSpacerItem(
            10, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 14, 1, 1, 4)
        scroll.setWidget(dlg)

        if config and config.has_section("PFD"):
            self.ColorButtonBrush.setColor(
                config.get("PFD", 'brushColor'))
            self.brush.setCurrentIndex(config.getint("PFD", "brush"))
            self.ColorButtonEntrada.setColor(
                config.get("PFD", 'Color_Entrada'))
            self.ColorButtonSalida.setColor(config.get("PFD", 'Color_Salida'))
            self.lineFormat.ColorButtonLine.setColor(
                config.get("PFD", 'Color_Stream'))
            self.lineFormat.groupJoint.button(
                (config.getint("PFD", 'Union')+2)*-1).setChecked(True)
            self.lineFormat.mitterLimit.setValue(
                config.getfloat("PFD", 'Miter_limit'))
            self.lineFormat.groupCap.button(
                (config.getint("PFD", 'Punta')+2)*-1).setChecked(True)
            self.lineFormat.guion.setCurrentIndex(
                config.getint("PFD", 'Guion'))
            self.lineFormat.dashOffset.setValue(
                config.getfloat("PFD", 'Dash_offset'))
            self.lineFormat.width.setValue(config.getfloat("PFD", 'Width'))

    def value(self, config):
        """Update ConfigParser instance with the config"""
        if not config.has_section("PFD"):
            config.add_section("PFD")
        config = self.resolution.value(config)
        config.set("PFD", "brushColor",
                   self.ColorButtonBrush.color.name())
        config.set("PFD", "brush", str(self.brush.currentIndex()))
        config.set("PFD", "Color_Entrada",
                   self.ColorButtonEntrada.color.name())
        config.set("PFD", "Color_Salida", self.ColorButtonSalida.color.name())
        config.set("PFD", "Color_Stream",
                   self.lineFormat.ColorButtonLine.color.name())
        config.set("PFD", "Width", str(self.lineFormat.width.value))
        config.set("PFD", "Union",
                   str(abs(self.lineFormat.groupJoint.checkedId())-2))
        config.set("PFD", "Miter_limit",
                   str(self.lineFormat.mitterLimit.value))
        config.set("PFD", "Punta",
                   str(abs(self.lineFormat.groupCap.checkedId())-2))
        config.set("PFD", "Guion", str(self.lineFormat.guion.currentIndex()))
        config.set("PFD", "Dash_offset", str(self.lineFormat.dashOffset.value))
        return config


class Dialog(QtWidgets.QDialog):
    """Dialog to config thermal method calculations"""
    def __init__(self, config=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle(self.tr("PFD configuration"))
        layout = QtWidgets.QVBoxLayout(self)
        self.widget = Widget(config)
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


class BrushCombo(QtWidgets.QComboBox):
    """Custom QComboBox to choose a QBrush style shown the selected style"""

    BRUSH = [
        QtCore.Qt.BrushStyle.NoBrush,
        QtCore.Qt.BrushStyle.Dense7Pattern,
        QtCore.Qt.BrushStyle.Dense6Pattern,
        QtCore.Qt.BrushStyle.Dense5Pattern,
        QtCore.Qt.BrushStyle.Dense4Pattern,
        QtCore.Qt.BrushStyle.Dense3Pattern,
        QtCore.Qt.BrushStyle.Dense2Pattern,
        QtCore.Qt.BrushStyle.Dense1Pattern,
        QtCore.Qt.BrushStyle.SolidPattern,
        QtCore.Qt.BrushStyle.HorPattern,
        QtCore.Qt.BrushStyle.VerPattern,
        QtCore.Qt.BrushStyle.CrossPattern,
        QtCore.Qt.BrushStyle.BDiagPattern,
        QtCore.Qt.BrushStyle.FDiagPattern,
        QtCore.Qt.BrushStyle.DiagCrossPattern]

    def __init__(self, *args, **kwargs):
        """Autofill of widget with accepted values"""
        super().__init__(*args, **kwargs)

        for style in self.BRUSH:
            pix = QtGui.QPixmap(50, 50)
            pix.fill(QtGui.QColorConstants.White)

            painter = QtGui.QPainter()
            painter.begin(pix)
            brush = QtGui.QBrush(style)
            painter.setBrush(brush)
            painter.drawRect(0, 0, 50, 50)
            icon = QtGui.QIcon(pix)
            painter.end()
            self.addItem(icon, str(style).split(".")[-1])

    def paintEvent(self, event):
        """Paint the widget with the selected QBrush"""
        # Paint default element
        painter = QtWidgets.QStylePainter(self)
        opt = QtWidgets.QStyleOptionComboBox()
        self.initStyleOption(opt)
        painter.drawComplexControl(
            QtWidgets.QStyle.ComplexControl.CC_ComboBox, opt)

        # Apply selected brush
        brush = self.BRUSH[self.currentIndex()]
        painter.setBrush(brush)

        # Don't paint the arrow region in QComboBox
        rect = event.rect()
        rect.setRight(rect.right()-24)
        painter.drawRect(rect)


class ConfLine(QtWidgets.QWidget):
    """Composite widget with line format configuration tools"""
    join = [QtCore.Qt.PenJoinStyle.MiterJoin,
            QtCore.Qt.PenJoinStyle.BevelJoin,
            QtCore.Qt.PenJoinStyle.RoundJoin]
    cap = [QtCore.Qt.PenCapStyle.FlatCap,
           QtCore.Qt.PenCapStyle.RoundCap,
           QtCore.Qt.PenCapStyle.SquareCap]
    line = [QtCore.Qt.PenStyle.SolidLine,
            QtCore.Qt.PenStyle.DashLine,
            QtCore.Qt.PenStyle.DotLine,
            QtCore.Qt.PenStyle.DashDotLine,
            QtCore.Qt.PenStyle.DashDotDotLine]

    def __init__(self, pen=None, parent=None):
        super().__init__(parent)
        lyt = QtWidgets.QVBoxLayout(self)

        lyt1 = QtWidgets.QHBoxLayout()
        lyt1.addWidget(QtWidgets.QLabel(
            self.tr("Line")))
        self.ColorButtonLine = ColorSelector()
        self.ColorButtonLine.setToolTip(
            self.tr("Default line color"))
        lyt1.addWidget(self.ColorButtonLine)
        self.width = Entrada_con_unidades(
            float, width=50, decimales=1, spinbox=True, step=0.1,
            textounidad="px")
        self.width.entrada.setToolTip(self.tr("Line Width"))
        lyt1.addWidget(self.width)
        lyt.addLayout(lyt1)

        lyt2 = QtWidgets.QHBoxLayout()
        lyt2.addWidget(QtWidgets.QLabel(self.tr("Join")))
        self.mitterLimit = Entrada_con_unidades(
            float, width=50, decimales=1, spinbox=True, step=0.1)
        self.mitterLimit.entrada.setToolTip(self.tr("Mitter Limit"))
        lyt2.addWidget(self.mitterLimit)
        toolJoinMitter = QtWidgets.QToolButton()
        toolJoinMitter.setIcon(QtGui.QIcon(QtGui.QPixmap(os.path.join(
            os.environ["pychemqt"], "images", "button",
            "stroke-join-miter.png"))))
        toolJoinMitter.setIconSize(QtCore.QSize(24, 24))
        toolJoinMitter.setCheckable(True)
        toolJoinMitter.setToolTip(self.tr(
            "Join mitter: The triangular notch between the two lines is not "
            "filled"))
        lyt2.addWidget(toolJoinMitter)
        toolJoinBevel = QtWidgets.QToolButton()
        toolJoinBevel.setIcon(QtGui.QIcon(QtGui.QPixmap(os.path.join(
            os.environ["pychemqt"], "images", "button",
            "stroke-join-bevel.png"))))
        toolJoinBevel.setIconSize(QtCore.QSize(24, 24))
        toolJoinBevel.setCheckable(True)
        toolJoinBevel.setToolTip(self.tr(
            "Join bevel: The triangular notch between the two lines is "
            "filled"))
        lyt2.addWidget(toolJoinBevel)
        toolJoinRound = QtWidgets.QToolButton()
        toolJoinRound.setIcon(QtGui.QIcon(QtGui.QPixmap(os.path.join(
            os.environ["pychemqt"], "images", "button",
            "stroke-join-round.png"))))
        toolJoinRound.setIconSize(QtCore.QSize(24, 24))
        toolJoinRound.setCheckable(True)
        toolJoinRound.setToolTip(self.tr(
            "Join round: A circular arc between the two lines is filled"))
        lyt2.addWidget(toolJoinRound)

        self.groupJoint = QtWidgets.QButtonGroup()
        self.groupJoint.addButton(toolJoinMitter)
        self.groupJoint.addButton(toolJoinBevel)
        self.groupJoint.addButton(toolJoinRound)
        self.groupJoint.idClicked.connect(self.mitterlimitEnabled)
        lyt.addLayout(lyt2)

        lyt3 = QtWidgets.QHBoxLayout()
        lyt3.addWidget(QtWidgets.QLabel(self.tr("Cap")))
        toolCapFlat = QtWidgets.QToolButton()
        toolCapFlat.setIcon(QtGui.QIcon(QtGui.QPixmap(os.path.join(
            os.environ["pychemqt"], "images", "button",
            "stroke-cap-butt.png"))))
        toolCapFlat.setIconSize(QtCore.QSize(24, 24))
        toolCapFlat.setCheckable(True)
        toolCapFlat.setToolTip(self.tr(
            "Flat Cap: A square line end that does not cover the end point of "
            "the line"))
        lyt3.addWidget(toolCapFlat)
        toolCapRound = QtWidgets.QToolButton()
        toolCapRound.setIcon(QtGui.QIcon(QtGui.QPixmap(os.path.join(
            os.environ["pychemqt"], "images", "button",
            "stroke-cap-round.png"))))
        toolCapRound.setIconSize(QtCore.QSize(24, 24))
        toolCapRound.setCheckable(True)
        toolCapRound.setToolTip(self.tr("Round Cap: A rounded line end"))
        lyt3.addWidget(toolCapRound)
        toolCapSquare = QtWidgets.QToolButton()
        toolCapSquare.setIcon(QtGui.QIcon(QtGui.QPixmap(os.path.join(
            os.environ["pychemqt"], "images", "button",
            "stroke-cap-square.png"))))
        toolCapSquare.setIconSize(QtCore.QSize(24, 24))
        toolCapSquare.setCheckable(True)
        toolCapSquare.setToolTip(self.tr(
            "Square Cap: A square line end that covers the end point and "
            "extends beyond it by half the line width"))
        lyt3.addWidget(toolCapSquare)

        self.groupCap = QtWidgets.QButtonGroup()
        self.groupCap.addButton(toolCapFlat)
        self.groupCap.addButton(toolCapRound)
        self.groupCap.addButton(toolCapSquare)
        lyt.addLayout(lyt3)

        lyt4 = QtWidgets.QHBoxLayout()
        lyt4.addWidget(QtWidgets.QLabel(self.tr("Dash")))
        self.guion = PFDLineCombo()
        lyt4.addWidget(self.guion)
        self.dashOffset = Entrada_con_unidades(
            float, width=50, decimales=1, spinbox=True, step=0.1)
        self.dashOffset.entrada.setToolTip(self.tr("Dash offset"))
        lyt4.addWidget(self.dashOffset)
        lyt.addLayout(lyt4)
        lyt.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding))

        if pen:
            self.ColorButtonLine.setColor(pen.color().name())
            self.groupJoint.button((self.join.index(
                pen.joinStyle())+2)*-1).setChecked(True)
            self.mitterLimit.setValue(pen.miterLimit())
            self.groupCap.button((self.cap.index(
                pen.capStyle())+2)*-1).setChecked(True)
            self.guion.setCurrentIndex(self.line.index(pen.style()))
            self.dashOffset.setValue(pen.dashOffset())
            self.width.setValue(pen.widthF())

    def mitterlimitEnabled(self, index):
        """Enable mitterLimit widget only if joinmitter is enabled"""
        self.mitterLimit.setEnabled(index == -2)

    def pen(self):
        """Return a QPen with the live configuration"""
        pen = QtGui.QPen(QtGui.QColor(self.ColorButtonLine.color.name()))
        pen.setWidthF(self.width.value)
        pen.setJoinStyle(self.join[abs(self.groupJoint.checkedId())-2])
        pen.setMiterLimit(self.mitterLimit.value)
        pen.setCapStyle(self.cap[abs(self.groupCap.checkedId())-2])
        pen.setStyle(self.line[self.guion.currentIndex()])
        pen.setDashOffset(self.dashOffset.value)
        return pen


class ConfLineDialog(QtWidgets.QDialog, ConfLine):
    """Dialog derived ConfLine for standalone use"""

    def __init__(self, pen=None, parent=None):
        super().__init__(pen, parent)
        self.setWindowTitle(self.tr("Edit format line"))
        buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok
            | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)
        self.layout().addWidget(buttonBox)
