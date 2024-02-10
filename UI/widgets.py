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
# Module to define common graphics widget
#   - Status: Label with status (for equipment, stream)
#   - Entrada_con_unidades: Composite widget for unit values for input/view
#   - Tabla: Custom tablewidget tablewidget with added functionality
#   - ClickableLabel: Label with custom clicked signal
#   - ColorSelector: Composite widget for colour definition
#   - DragButton: Button with drag & drop support
#   - PathConfig: Custom widget for file path show and configure
#   - LineConfig: Custom QGroupbox with all matplotlib Line configuration
#   - CustomCombo: General custom QComboBox
#   - LineStyleCombo: Custom QComboBox for select matplotlib line styles
#   - MarkerCombo: Custom QComboBox for select matplotlib line marker
#   - NumericFactor: Numeric format configuration dialog
#   - InputFont: Custom widget to edit a text input with font and color support
#   - Table_Graphics: Custom widget to implement a popup in PFD
#   - QLabelMath: Customized QLabel to show a pixmap from a latex code

#   - createAction
#   - okToContinue: Function to ask user if any unsaved change
#   - mathTex2QPixmap: Convert a latex text to a QPixmap
###############################################################################


from configparser import ConfigParser
from math import pi
import os
import sys

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

from lib.config import conf_dir, IMAGE_PATH
from lib.corriente import Corriente
from lib.utilities import representacion
from tools.qt import QtCore, QtGui, QtWidgets, tr
from tools.UI_unitConverter import UI_conversorUnidades, moneda
from UI.delegate import CellEditor


class Status(QtWidgets.QLabel):
    """Widget with status of dialog, equipment, stream, project, ..."""
    status = (
        (0, tr("pychemqt", "Underspecified"), "yellow"),
        (1, tr("pychemqt", "Solved"), "green"),
        (2, tr("pychemqt", "Ignored"), "Light gray"),
        (3, tr("pychemqt", "Warning"), "green"),
        (4, tr("pychemqt", "Calculating..."), "Cyan"),
        (5, tr("pychemqt", "Error"), "red"))

    def __init__(self, state=0, text="", parent=None):
        """
        state:
            0   -   Not solved
            1   -   OK
            2   -   Ignore
            3   -   Warning (Recommend: Use text parameter to explain)
            4   -   Calculating
            5   -   Error
        """
        super().__init__(parent)
        self.setState(state, text)
        self.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.setFrameShadow(QtWidgets.QFrame.Shadow.Sunken)
        self.setFrameShape(QtWidgets.QFrame.Shape.Panel)
        self.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding,
                           QtWidgets.QSizePolicy.Policy.Preferred)
        self.state = 0
        self.oldState = 0
        self.oldText = ""

    def setState(self, state, text=""):
        """Change the state"""

        if state == 2:
            self.oldState = self.state
            oldtext = self.text().split(": ")
            if len(oldtext) == 1:
                self.oldText = ""
            else:
                self.oldText = oldtext[1:].join(": ")
        if text:
            self.setText(self.status[state][1]+": "+text)
        else:
            self.setText(self.status[state][1])

        color = self.status[state][2]
        self.setStyleSheet(f"QLabel {{background-color: {color}}}")
        QtWidgets.QApplication.processEvents()
        self.state = state

    def restaurar(self):
        """Restore old stade"""
        self.setState(self.oldState, self.oldText)


class Entrada_con_unidades(QtWidgets.QWidget):
    """Customized widget with unit functionality"""

    valueChanged = QtCore.pyqtSignal(float)

    def __init__(self, unidad, UIconfig=None, retornar=True, readOnly=False,
                 boton=True, texto=True, textounidad="", title="", value=None,
                 start=0, max=float("inf"), min=0, decimales=4, tolerancia=4,
                 parent=None, width=85, resaltado=False, spinbox=False,
                 suffix="", step=0.01, colorReadOnly=None, colorResaltado=None,
                 frame=True, showNull=False):
        """
        Units:
            unidad: The unit (lib/unidades class) to use, mandatory
            UIconfig: Magnitud necessary if the main unit have several meaning
            title: Update unit title property
            retornar: Boolean to let or avoid the conversion window update
                the value of widget
            value: Inicial value of widget
            max: Maximum value for widget
            min: Minimum value for widget
            decimales: Decimal number count to show of value
            tolerancia: Value of exponent over than to use exponential notation
        UI:
            readOnly: Boolean, set widget readOnly property
            frame: Boolean, show the frame of widget or not
            width: Width of value widget
            boton: Boolean, show or not the button for unit conversion dialog
            texto: Boolean, show the unit text at right of value
            textounidad: Alternate text to show as unit text
            suffix: Text added to value in value representation
            showNull: Boolean, show value if it's 0
            resaltado: Boolean, use base color in widget
            colorResaltado: Color to use as base color if value
            colorReadOnly: Color to use is the widget is readOnly
        Spinbox functionality:
            spinbox: boolean to specified a QSpinbox use, with mouse response
            start: initial value for spinbox mouse interaction
            step: value of step at mouse spingox interaction
        """
        super().__init__(parent)
        self.resize(self.minimumSize())
        self.unidad = unidad

        if title:
            self.unidad.__title__ = title
        if unidad in (float, int):
            self.magnitud = None
        else:
            self.magnitud = unidad.__name__
        if unidad == int and spinbox and step == 0.01:
            step = 1
        self.decimales = decimales
        self.tolerancia = tolerancia
        self.step = step
        self.spinbox = spinbox
        self.max = max
        self.suffix = suffix
        self.min = min
        self.start = start
        self.textounidad = textounidad
        self.boton = boton
        self.resaltado = resaltado
        self.showNull = showNull

        Config = ConfigParser()
        Config.read(conf_dir+"pychemqtrc")
        if colorReadOnly:
            self.colorReadOnly = colorReadOnly
        else:
            self.colorReadOnly = QtGui.QColor(
                Config.get("General", 'Color_ReadOnly'))
        if colorResaltado:
            self.colorResaltado = colorResaltado
        else:
            self.colorResaltado = QtGui.QColor(
                Config.get("General", 'Color_Resaltado'))

        if UIconfig:
            self.UIconfig = UIconfig
        else:
            self.UIconfig = self.magnitud
        self.retornar = retornar
        layout = QtWidgets.QGridLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.entrada = QtWidgets.QLineEdit()
        self.entrada.setFixedSize(width, 24)
        self.entrada.editingFinished.connect(self.entrada_editingFinished)
        self.entrada.setAlignment(
            QtCore.Qt.AlignmentFlag.AlignRight
            | QtCore.Qt.AlignmentFlag.AlignVCenter)
        if unidad == int:
            if max == float("inf"):
                max = 1000000000
            validator = QtGui.QIntValidator(min, max, self)
        else:
            validator = QtGui.QDoubleValidator(min, max, decimales, self)
            locale = QtCore.QLocale("en")
            validator.setLocale(locale)
        self.entrada.setValidator(validator)
        self.setReadOnly(readOnly)
        self.setRetornar(self.retornar)
        self.setFrame(frame)
        layout.addWidget(self.entrada, 0, 1, 1, 3)

        if value is None:
            self.value = self.unidad(0)
        else:
            self.setValue(value)
        if self.magnitud:
            if boton:
                self.unidades = QtWidgets.QPushButton(".")
                self.unidades.setFixedSize(12, 24)
                self.unidades.setVisible(False)
                self.unidades.clicked.connect(self.unidades_clicked)
                layout.addWidget(self.unidades, 0, 1)

        if boton:
            self.botonClear = QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(
                os.path.join(IMAGE_PATH, "button", "editDelete.png"))), "")
            self.botonClear.setFixedSize(12, 24)
            self.botonClear.setVisible(False)
            self.botonClear.clicked.connect(self.clear)
            layout.addWidget(self.botonClear, 0, 3)

        if texto:
            self.texto = QtWidgets.QLabel()
            self.texto.setAlignment(QtCore.Qt.AlignmentFlag.AlignVCenter)
            self.texto.setIndent(5)
            txt = ""
            if self.UIconfig:
                txt += self.value.text(self.UIconfig)
            if textounidad:
                txt += textounidad
            self.texto.setText(txt)
            layout.addWidget(self.texto, 0, 4)

        layout.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed), 0, 5)
        self.setResaltado(resaltado)

    def unidades_clicked(self):
        """Show the unit converter dialog"""
        if self.magnitud == "Currency":
            dialog = moneda(self.value)
        else:
            dialog = UI_conversorUnidades(self.unidad, self.value)

        if dialog.exec() and self.retornar:
            # Change the value if change and retornar if active
            self.entrada.setText(
                representacion(dialog.value.config(self.UIconfig))+self.suffix)
            oldvalue = self.value
            self.value = dialog.value
            if oldvalue != self.value:
                self.valueChanged.emit(self.value)

    def entrada_editingFinished(self):
        """Change the value at finish of edit"""
        if not self.readOnly:
            # Filter suffix and fix bad numeric , interpretation
            if self.suffix:
                txt = self.entrada.text().split(self.suffix).replace(',', '.')
            else:
                txt = self.entrada.text().replace(',', '.')
            if self.unidad != int:
                self.entrada.setText(
                    representacion(float(txt), decimales=self.decimales,
                                   tol=self.tolerancia)+self.suffix)
            oldvalue = self.value
            if self.magnitud:
                self.value = self.unidad(
                    float(txt), "conf", magnitud=self.UIconfig)
            else:
                self.value = self.unidad(txt)
            if self.value != oldvalue or self.value == 0:
                self.valueChanged.emit(self.value)
                self.setToolTip()

    def clear(self):
        """Clear value"""
        self.entrada.setText("")
        self.value = None

    def setResaltado(self, boolean):
        """Set resalted background color for input widget"""
        self.resaltado = boolean
        paleta = QtGui.QPalette()
        if boolean:
            paleta.setColor(QtGui.QPalette.ColorRole.Base,
                            QtGui.QColor(self.colorResaltado))
        elif self.readOnly:
            paleta.setColor(QtGui.QPalette.ColorRole.Base,
                            QtGui.QColor(self.colorReadOnly))
        else:
            paleta.setColor(QtGui.QPalette.ColorRole.Base,
                            QtGui.QColor("white"))
        self.entrada.setPalette(paleta)

    def setReadOnly(self, readOnly):
        """Set read only state for widget"""
        self.entrada.setReadOnly(readOnly)
        self.readOnly = readOnly
        self.setResaltado(self.resaltado)

    def setNotReadOnly(self, editable):
        self.setReadOnly(not editable)

    def setRetornar(self, retornar):
        self.retornar = retornar

    def setValue(self, value):
        if value is None:
            value = 0
        self.value = self.unidad(value)
        if value or self.showNull:
            if self.magnitud:
                self.entrada.setText(
                    self.value.format(magnitud=self.UIconfig)+self.suffix)
            elif self.unidad == float:
                self.entrada.setText(
                    representacion(self.value, decimales=self.decimales,
                                   tol=self.tolerancia)+self.suffix)
            else:
                self.entrada.setText(str(self.value)+self.suffix)
            self.setToolTip()

    def setFrame(self, frame):
        self.entrada.setFrame(frame)
        self.frame = frame

    def setToolTip(self):
        """Define the tooltip with the values in confguration"""
        Preferences = ConfigParser()
        Preferences.read(conf_dir+"pychemqtrc")
        if Preferences.getboolean("Tooltip", "Show"):
            Config = ConfigParser()
            Config.read(conf_dir+"pychemqtrc")
            try:
                lst = map(int, Config.get("Tooltip", self.magnitud).split(","))
            except AttributeError:
                lst = []
            if lst:
                valores = []
                for i in lst:
                    valores.append(representacion(
                        self.value.__getattribute__(self.value.__units__[i]),
                        self.decimales, self.tolerancia) + " " +
                        self.value.__text__[i])
                self.entrada.setToolTip(os.linesep.join(valores))

    def keyPressEvent(self, e):
        """Manage the key press to emulate a QSpinbox"""
        if not self.readOnly:
            if e.key() in [QtCore.Qt.Key.Key_Insert, QtCore.Qt.Key.Key_Backspace]:
                self.clear()
            if self.spinbox:
                if not self.value:
                    self.value = self.start
                if e.key() == QtCore.Qt.Key.Key_Up:
                    valor = self.value+self.step
                    if valor > self.max:
                        self.setValue(self.max)
                    else:
                        self.setValue(valor)
                elif e.key() == QtCore.Qt.Key.Key_Down:
                    valor = self.value-self.step
                    if valor < self.min:
                        self.setValue(self.min)
                    else:
                        self.setValue(valor)
                self.valueChanged.emit(self.value)

    def enterEvent(self, event):
        """When mouse enter in widget show the unidades and clear button, and
        add margin to let space to clear button"""
        if self.magnitud and self.boton:
            self.unidades.setVisible(True)
        if self.value and self.boton and not self.readOnly:
            self.botonClear.setVisible(True)
            self.entrada.setTextMargins(0, 0, 10, 0)

    def leaveEvent(self, event):
        """When mouse leave the widget undo the enterEvent actions"""
        if self.magnitud and self.boton:
            self.unidades.setVisible(False)
        if self.value and self.boton and not self.readOnly:
            self.botonClear.setVisible(False)
            self.entrada.setTextMargins(0, 0, 0, 0)


class Tabla(QtWidgets.QTableWidget):
    """QTableWidget with custom functionality"""
    editingFinished = QtCore.pyqtSignal()
    rowFinished = QtCore.pyqtSignal(list)

    def __init__(self, columnas=0, filas=0, stretch=True, dinamica=False,
                 readOnly=False, columnReadOnly=None,
                 horizontalHeader=None, verticalHeader=True,
                 verticalHeaderLabels=None, verticalHeaderModel="",
                 verticalOffset=0, orientacion=QtCore.Qt.AlignmentFlag.AlignRight,
                 delegate=CellEditor, delegateforRow=None,
                 parent=None):
        """
        columnas: Column count of widget
        filas: Row count, this value is initial and can be changed
        stretch: Boolean, stretch the last column to fill all space available
        dinamica: Boolean, let user fill the data adding new row when last
            are filled
        readOnly: Boolean, set the readOnly state of widget
        columnReadOnly: Array with boolean for column readOnly state, used when
            the readOnly state is different for each column

        horizontalHeader: Array with text for top header,
            Null don't show the top header
        verticalHeader: Boolean, to show or hide the right header
        verticalHeaderLabels: Array with text for right header
        verticalHeaderModel: If the table is dinamica and the right header is
            show, this are the text model to generate new vertical header label

        verticalOffset: Index, row don't included in fill data, used for show
            info, widget...
        orientation: Horizontal text alignment general for cell in tables,
            default is right as normal for numbers

        delegate: QItemDelegate subclass to configure cell editor
            default CellEditor with float appropiate functionality
        delegateforRow: delegate with specific values for row
        """
        super().__init__(parent)

        # Dimensions
        self.columnas = columnas
        self.filas = filas+verticalOffset
        self.verticalOffset = verticalOffset
        self.setColumnCount(self.columnas)

        # Header labels
        if not verticalHeader:
            self.verticalHeader().hide()
        if not horizontalHeader:
            self.horizontalHeader().hide()
        else:
            self.setHorizontalHeaderLabels(horizontalHeader)

        self.horizontalHeaderLabel = horizontalHeader
        self.verticalHeaderLabel = verticalHeaderLabels
        self.verticalHeaderBool = verticalHeader
        self.verticalHeaderModel = verticalHeaderModel
        self.horizontalHeader().setStretchLastSection(stretch)

        # readOnly state
        self.readOnly = readOnly
        if readOnly:
            self.setEditTriggers(
                QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        else:
            self.setEditTriggers(
                QtWidgets.QAbstractItemView.EditTrigger.AllEditTriggers)
        if columnReadOnly is None:
            self.columnReadOnly = [self.readOnly]*self.columnas
        else:
            self.columnReadOnly = columnReadOnly
            for i in range(self.columnCount()):
                self.setColumnReadOnly(i, columnReadOnly[i])

        # Delegate functionality
        if delegate:
            self.setItemDelegate(delegate(self))
        self.delegateforRow = delegateforRow

        if dinamica:
            self.cellChanged.connect(self.tabla_cellChanged)
        self.dinamica = dinamica

        # self.setAlternatingRowColors(True)
        self.setGridStyle(QtCore.Qt.PenStyle.DotLine)
        self.orientacion = orientacion
        for i in range(filas):
            self.addRow()

    def setConnected(self):
        """The dynamic state can be defined at start or call this procedure"""
        self.cellChanged.connect(self.tabla_cellChanged)
        self.dinamica = True
        if self.rowCount() == 0:
            self.addRow()

    def addRow(self, data=None, index=None):
        """Add row to widget
        data: Array with data to fill new row
        index: Index to add row, default add row at endo of table"""
        if not data:
            data = [""]*self.columnas
        else:
            data = [representacion(i) for i in data]
        if index is not None:
            row = index
        else:
            row = self.rowCount()
        self.insertRow(row)
        self.setRowHeight(row, 22)

        if self.delegateforRow:
            delegate = self.delegateforRow(self.parent())
            self.setItemDelegateForRow(row, delegate)

        Config = ConfigParser()
        Config.read(conf_dir+"pychemqtrc")
        inactivo = QtGui.QColor(Config.get("General", 'Color_ReadOnly'))
        for j in range(self.columnCount()):
            self.setItem(row, j, QtWidgets.QTableWidgetItem(data[j]))
            self.item(row, j).setTextAlignment(
                self.orientacion | QtCore.Qt.AlignmentFlag.AlignVCenter)

            if self.columnReadOnly[j]:
                flags = QtCore.Qt.ItemFlag.ItemIsEnabled \
                    | QtCore.Qt.ItemFlag.ItemIsSelectable
                self.item(row, j).setBackground(inactivo)
            else:
                flags = QtCore.Qt.ItemFlag.ItemIsEditable \
                    | QtCore.Qt.ItemFlag.ItemIsEnabled \
                    | QtCore.Qt.ItemFlag.ItemIsSelectable
            self.item(row, j).setFlags(flags)

        self.setVHeader(row)

        # Set focus to first editable cell in new row
        if self.dinamica and self.rowCount() > 1:
            columna = self.columnReadOnly.index(False)
            self.setCurrentCell(row, columna)

    def setVHeader(self, row):
        """Set vertical header text"""
        if self.verticalHeaderBool:
            if self.verticalHeaderLabel:
                txt = self.verticalHeaderLabel[row]
            elif self.verticalHeaderModel:
                txt = self.verticalHeaderModel+str(row)
            else:
                txt = str(row+1)
            self.setVerticalHeaderItem(row, QtWidgets.QTableWidgetItem(txt))

    def tabla_cellChanged(self, i, j):
        """When edit a cell, check status to add new row"""
        new_line = True
        col = 0
        while col < self.columnas:
            if self.item(i, col).text() != "" or self.columnReadOnly[col]:
                col += 1
            else:
                new_line = False
                break
        if new_line and i == self.rowCount()-1:
            self.addRow()
            fila = self.getRow(i)
            self.rowFinished.emit(fila)

    def getValue(self, row, column):
        """Get value from cell in row and column"""
        txt = self.item(row, column).text()
        try:
            value = float(txt)
        except ValueError:
            value = txt
        return value

    def setValue(self, row, column, value, **fmt):
        """Set value for cell in row and column with text formating"""
        if isinstance(value, (float, int)):
            value = representacion(value, **fmt)
        self.item(row, column).setText(value)
        self.item(row, column).setTextAlignment(
            self.orientacion | QtCore.Qt.AlignmentFlag.AlignVCenter)

    def setColumn(self, column, data, **fmt):
        """Set data for a complete column"""
        while len(data) > self.rowCount()-self.verticalOffset:
            self.addRow()
        for row, dato in enumerate(data):
            self.setValue(row, column, dato, **fmt)

    def getColumn(self, column):
        """Get column data as array"""
        lst = []
        for row in range(self.verticalOffset, self.rowCount()):
            value = self.getValue(row, column)
            if isinstance(value, float):
                lst.append(value)
        return lst

    def getRow(self, row):
        """Get row data as array"""
        lst = []
        for column in range(self.columnCount()):
            lst.append(self.getValue(row, column))
        return lst

    def getData(self):
        """Get all table data as array"""
        matriz = []
        for row in range(self.verticalOffset, self.rowCount()-1):
            lst = self.getRow(row)
            matriz.append(lst)
        return matriz

    def setData(self, matriz):
        """Set table data"""
        for i in range(self.rowCount(), len(matriz)+self.verticalOffset):
            self.addRow()
        for fila in range(self.rowCount()-self.verticalOffset):
            for columna, dato in enumerate(matriz[fila]):
                # self.setVerticalHeaderItem(fila, QtWidgets.QTableWidgetItem(
                    # self.verticalHeaderModel+str(fila)))
                self.item(fila, columna).setText(str(dato))
                self.item(fila+self.verticalOffset, columna).setTextAlignment(
                    self.orientacion | QtCore.Qt.AlignmentFlag.AlignVCenter)
        for i in range(self.verticalOffset, self.rowCount()):
            self.setRowHeight(i+self.verticalOffset, 20)

    def clear(self, size=True):
        """Clear table, remove all data and optionally remove all row"""
        if size:
            self.setRowCount(1+self.verticalOffset)
        for fila in range(self.rowCount()):
            for columna in range(self.columnas):
                self.item(fila+self.verticalOffset, columna).setText("")

    def setColumnReadOnly(self, column, boolean):
        """Set readonly estate per column"""
        if boolean:
            flags = QtCore.Qt.ItemFlag.ItemIsEnabled \
                | QtCore.Qt.ItemFlag.ItemIsSelectable
        else:
            flags = QtCore.Qt.ItemFlag.ItemIsEditable \
                | QtCore.Qt.ItemFlag.ItemIsEnabled \
                | QtCore.Qt.ItemFlag.ItemIsSelectable

        for row in range(self.rowCount()):
            self.item(row, column).setFlags(flags)

    def leaveEvent(self, event):
        if self.isEnabled():
            self.editingFinished.emit()


class ClickableLabel(QtWidgets.QLabel):
    """Custom QLabel with clicked functionality"""
    clicked = QtCore.pyqtSignal()

    def mousePressEvent(self, event):
        self.clicked.emit()


class ColorSelector(QtWidgets.QWidget):
    """Color selector widget"""
    valueChanged = QtCore.pyqtSignal('QString')

    def __init__(self, color="#ffffff", alpha=255, isAlpha=False, parent=None):
        super().__init__(parent)

        lyt = QtWidgets.QHBoxLayout(self)
        lyt.setContentsMargins(0, 0, 0, 0)
        lyt.setSpacing(0)

        self.RGB = QtWidgets.QLineEdit()
        self.RGB.editingFinished.connect(self.rgbChanged)
        self.RGB.setFixedSize(80, 24)
        lyt.addWidget(self.RGB)
        self.button = QtWidgets.QToolButton()
        self.button.setFixedSize(24, 24)
        self.button.clicked.connect(self.ColorButtonClicked)
        lyt.addWidget(self.button)
        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed))

        self.isAlpha = isAlpha
        if isAlpha:
            self.name = QtGui.QColor.NameFormat.HexArgb
        else:
            self.name = QtGui.QColor.NameFormat.HexRgb

        r = int(color[1:3], 16)
        g = int(color[3:5], 16)
        b = int(color[5:7], 16)
        color = QtGui.QColor(r, g, b, alpha)
        self.setColor(color)

    def setColor(self, color, alpha=255):
        """Set new color value and update text and button color"""
        # Accept color args as a #rgb string too
        if isinstance(color, str):
            color = QtGui.QColor(color)
            color.setAlpha(alpha)
        self.color = color
        self.button.setStyleSheet(f"background: {color.name(self.name)};")
        self.RGB.setText(color.name(self.name))

    def ColorButtonClicked(self):
        """Show the QColorDialog to let user choose new color"""
        dlg = QtWidgets.QColorDialog(self.color, self)
        if self.isAlpha:
            dlg.setOption(
                QtWidgets.QColorDialog.ColorDialogOption.ShowAlphaChannel)
        if dlg.exec():
            self.setColor(dlg.currentColor())
            self.valueChanged.emit(dlg.currentColor().name())

    def rgbChanged(self):
        """Let user define the color manually"""
        txt = self.RGB.text()

        # Avoid the editing finished with no changes
        if txt == self.color.name(self.name):
            return

        # Define the new color from text
        # From Qt 6.4 qcolor can read a color by name including the SVG color
        # keyword names
        # https://www.w3.org/TR/SVG11/types.html#ColorKeywords
        # Anyway it's retain this manual reading code for qt5 compatibility
        if self.isAlpha:
            alpha = int(txt[1:3], 16)
            r = int(txt[3:5], 16)
            g = int(txt[5:7], 16)
            b = int(txt[7:9], 16)
            color = QtGui.QColor(r, g, b, alpha)
        else:
            color = QtGui.QColor(txt)

        # Only accept new value if it's valid
        if color.isValid():
            self.setColor(color)
            self.valueChanged.emit(color.name(self.name))


class DragButton(QtWidgets.QToolButton):
    """Clase que define un botón especial que permite arrastrar"""

    def __init__(self, parent=None):
        super().__init__(parent)

    # def mouseMoveEvent(self, event):
        # self.startDrag()
        # QtWidgets.QToolButton.mouseMoveEvent(self, event)

    # def startDrag(self):
        # if self.icon().isNull():
        #     return
        # data = QtCore.QByteArray()
        # stream = QtCore.QDataStream(data, QtCore.QIODevice.WriteOnly)
        # stream << self.icon()
        # mimeData = QtCore.QMimeData()
        # mimeData.setData("application/x-equipment", data)
        # drag = QtGui.QDrag(self)
        # drag.setMimeData(mimeData)
        # pixmap = self.icon().pixmap(24, 24)
        # drag.setHotSpot(QtCore.QPoint(12, 12))
        # drag.setPixmap(pixmap)
        # drag.exec(QtCore.Qt.DropAction.CopyAction)


class PathConfig(QtWidgets.QWidget):
    """Custom widget for a file path show and configure functionality"""
    valueChanged = QtCore.pyqtSignal('QString')

    def __init__(self, title="", path="", patron="", msg="", folder=False,
                 save=False, parent=None):
        """
        title: Optional text an right of widget
        path: Inicial value for file path
        patron: File format to filter in file search dialog
        msg: Title of dialog file
        """
        super().__init__(parent)

        self.folder = folder
        self.save = save

        layout = QtWidgets.QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        if title:
            layout.addWidget(QtWidgets.QLabel(title))
            layout.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
                QtWidgets.QSizePolicy.Policy.Fixed))

        self.path = QtWidgets.QLineEdit()
        self.path.setFixedHeight(24)
        self.path.textEdited.connect(self.pathEdited)
        layout.addWidget(self.path)
        self.boton = QtWidgets.QPushButton(
            tr("pychemqt", "Browse"))
        self.boton.setFixedHeight(24)
        self.boton.clicked.connect(self.select)
        layout.addWidget(self.boton)

        # Define default values for parameters don't defined
        if not patron:
            patron = tr(
                "pychemqt", "All files") + "(*)"
        elif patron == "exe":
            if sys.platform == "win32":
                patron = tr(
                    "pychemqt", "Executable files") + "(*.exe *.bat)"
            else:
                patron = tr(
                    "pychemqt", "All files") + "(*)"
        self.patron = patron

        if not msg:
            msg = tr(
                "pychemqt", "Select path of file")
        self.msg = msg
        self.setText(path)

    def text(self):
        return self.path.text()

    def setText(self, text):
        self.path.setText(text)

    def select(self):
        """Open the QFileDialog to select the file path"""
        folder = os.path.dirname(str(self.path.text()))
        if self.save:
            ruta = QtWidgets.QFileDialog.getSaveFileName(
                self, self.msg, folder, self.patron)[0]
        elif self.folder:
            ruta = QtWidgets.QFileDialog.getExistingDirectory(
                self, self.msg, folder)
        else:
            ruta = QtWidgets.QFileDialog.getOpenFileName(
                self, self.msg, folder, self.patron)[0]
        if ruta:
            self.path.setText(ruta)
            self.valueChanged.emit(ruta)

    def pathEdited(self, path):
        if os.path.isfile(path):
            self.valueChanged.emit(path)


class LineConfig(QtWidgets.QGroupBox):
    """Custom QGroupbox with all matplotlib Line configuration"""

    def __init__(self, confSection, title, parent=None):
        """
        confSection: Name key to identify the line
        title: Title to use in QGroupbox
        """
        super().__init__(title, parent)
        self.conf = confSection

        layout = QtWidgets.QVBoxLayout(self)
        lyt1 = QtWidgets.QHBoxLayout()
        self.Width = QtWidgets.QDoubleSpinBox()
        self.Width.setFixedWidth(60)
        self.Width.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight)
        self.Width.setRange(0.1, 5)
        self.Width.setDecimals(1)
        self.Width.setSingleStep(0.1)
        self.Width.setToolTip(tr("pychemqt", "Line width"))
        lyt1.addWidget(self.Width)
        self.Style = LineStyleCombo()
        self.Style.setToolTip(tr("pychemqt", "Line style"))
        lyt1.addWidget(self.Style)
        self.Color = ColorSelector(isAlpha=True)
        self.Color.setToolTip(tr("pychemqt", "Line color"))
        lyt1.addWidget(self.Color)
        self.Marker = MarkerCombo()
        self.Marker.setToolTip(tr("pychemqt", "Line marker"))
        self.Marker.currentIndexChanged.connect(self.changeMarker)
        lyt1.addWidget(self.Marker)
        lyt1.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed))
        layout.addLayout(lyt1)

        lyt2 = QtWidgets.QHBoxLayout()
        self.MarkerSize = QtWidgets.QDoubleSpinBox()
        self.MarkerSize.setFixedWidth(60)
        self.MarkerSize.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight)
        self.MarkerSize.setRange(0.1, 5)
        self.MarkerSize.setDecimals(1)
        self.MarkerSize.setSingleStep(0.1)
        self.MarkerSize.setToolTip(tr("pychemqt", "Marker size"))
        lyt2.addWidget(self.MarkerSize)
        self.MarkerColor = ColorSelector()
        self.MarkerColor.setToolTip(tr("pychemqt", "Marker face color"))
        lyt2.addWidget(self.MarkerColor)
        self.EdgeSize = QtWidgets.QDoubleSpinBox()
        self.EdgeSize.setFixedWidth(60)
        self.EdgeSize.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight)
        self.EdgeSize.setRange(0.1, 5)
        self.EdgeSize.setDecimals(1)
        self.EdgeSize.setSingleStep(0.1)
        self.EdgeSize.setToolTip(tr("pychemqt", "Marker edge width"))
        lyt2.addWidget(self.EdgeSize)
        self.EdgeColor = ColorSelector()
        self.EdgeColor.setToolTip(tr("pychemqt", "Marker edge color"))
        lyt2.addWidget(self.EdgeColor)
        layout.addLayout(lyt2)

        self.changeMarker(0)

    def changeMarker(self, index):
        self.MarkerSize.setVisible(index)
        self.MarkerColor.setVisible(index)
        self.EdgeSize.setVisible(index)
        self.EdgeColor.setVisible(index)

    def setConfig(self, config, section="MEOS"):
        alfa = config.getint(section, self.conf+"alpha")
        self.Color.setColor(config.get(section, self.conf+'Color'), alfa)
        self.Width.setValue(config.getfloat(section, self.conf+'lineWidth'))
        self.Style.setCurrentValue(config.get(section, self.conf+'lineStyle'))
        if self.Marker.isEnabled():
            self.Marker.setCurrentValue(config.get(section, self.conf+'marker'))
            self.MarkerSize.setValue(
                config.getfloat(section, self.conf+'markersize'))
            self.MarkerColor.setColor(
                config.get(section, self.conf+'markerfacecolor'), alfa)
            self.EdgeSize.setValue(
                config.getfloat(section, self.conf+'markeredgewidth'))
            self.EdgeColor.setColor(
                config.get(section, self.conf+'markeredgecolor'), alfa)

    def value(self, config, section="MEOS"):
        config.set(section, self.conf+"Color", self.Color.color.name())
        config.set(section, self.conf+"alpha", str(self.Color.color.alpha()))
        config.set(section, self.conf+"lineWidth", str(self.Width.value()))
        config.set(section, self.conf+"lineStyle", self.Style.currentValue())
        if self.Marker.isEnabled():
            config.set(section, self.conf+"marker", self.Marker.currentValue())
            config.set(section, self.conf+"markersize",
                       str(self.MarkerSize.value()))
            config.set(section, self.conf+"markerfacecolor",
                       self.MarkerColor.color.name())
            config.set(section, self.conf+"markeredgewidth",
                       str(self.EdgeSize.value()))
            config.set(section, self.conf+"markeredgecolor",
                       self.EdgeColor.color.name())

        return config


class GridConfig(LineConfig):
    """Extended version of LineConfig for grid specific line"""

    WHICH = ["both", "major", "minor"]
    AXIS = ["both", "x", "y"]

    def __init__(self, confSection, title, parent=None):
        """
        confSection: Name key to identify the line
        title: Title to use in QGroupbox
        """
        super().__init__(confSection, title, parent)

        # Disable marker functionality for grid line
        # Raise a error by duplicate marker kwarg in plot, anyway marker in
        # grid line are useless
        self.Marker.setEnabled(False)

        self.grid = QtWidgets.QCheckBox(tr(
            "pychemqt", "Show grid"))
        self.layout().insertWidget(0, self.grid)

        lyt = QtWidgets.QHBoxLayout()
        lyt.addWidget(QtWidgets.QLabel(
            tr("pychemqt", "Which") + ":"))
        self.gridWhich = QtWidgets.QComboBox()
        for name in self.WHICH:
            self.gridWhich.addItem(name)
        lyt.addWidget(self.gridWhich)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed))
        self.layout().insertLayout(1, lyt)

        lyt = QtWidgets.QHBoxLayout()
        lyt.addWidget(QtWidgets.QLabel(
            tr("pychemqt", "Axis") + ":"))
        self.gridAxis = QtWidgets.QComboBox()
        for name in self.AXIS:
            self.gridAxis.addItem(name)
        lyt.addWidget(self.gridAxis)
        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed))
        self.layout().insertLayout(2, lyt)

    def setConfig(self, config, section):
        LineConfig.setConfig(self, config, section)
        self.grid.setChecked(config.getboolean(section, 'grid'))
        txt = config.get(section, 'gridwhich')
        self.gridWhich.setCurrentIndex(self.WHICH.index(txt))
        txt = config.get(section, 'gridaxis')
        self.gridAxis.setCurrentIndex(self.AXIS.index(txt))

    def value(self, config, section):
        """Update ConfigParser instance with the config"""
        LineConfig.value(self, config, section)
        config.set(section, "grid", str(self.grid.isChecked()))
        config.set(section, "gridwhich", self.gridWhich.currentText())
        config.set(section, "gridaxis", self.gridAxis.currentText())
        return config


class CustomCombo(QtWidgets.QComboBox):
    """General custom QComboBox"""
    valueChanged = QtCore.pyqtSignal("QString")

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setIconSize(QtCore.QSize(35, 18))
        self.currentIndexChanged.connect(self.emit)
        self._populate()

    def setCurrentValue(self, value):
        """Overload QComboBox setCurrentValue to let text input from keymap"""
        ind = self.key.index(value)
        self.setCurrentIndex(ind)

    def currentValue(self):
        """Convert index in value from keymap"""
        return self.key[self.currentIndex()]

    def emit(self, ind):
        """Emit signal change integer for keymap value"""
        self.valueChanged.emit(self.key[ind])


class LineStyleCombo(CustomCombo):
    """Custom QComboBox for select matplotlib line styles"""
    key = ["None", "-", "--", ":", "-."]
    image = {
        "None": "",
        "-": os.path.join("images", "button", "solid_line.png"),
        "--": os.path.join("images", "button", "dash_line.png"),
        ":": os.path.join("images", "button", "dot_line.png"),
        "-.": os.path.join("images", "button", "dash_dot_line.png")}

    def _populate(self):
        for key in self.key:
            self.addItem(QtGui.QIcon(QtGui.QPixmap(
                os.environ["pychemqt"] + self.image[key])), "")


class PFDLineCombo(LineStyleCombo):
    """Custom QComboBox for select PFD line styles for stream"""
    key = [QtCore.Qt.PenStyle.SolidLine, QtCore.Qt.PenStyle.DashLine,
           QtCore.Qt.PenStyle.DotLine, QtCore.Qt.PenStyle.DashDotLine,
           QtCore.Qt.PenStyle.DashDotDotLine]
    image = {
        key[0]: os.path.join("images", "button", "solid_line.png"),
        key[1]: os.path.join("images", "button", "dash_line.png"),
        key[2]: os.path.join("images", "button", "dot_line.png"),
        key[3]: os.path.join("images", "button", "dash_dot_line.png"),
        key[4]: os.path.join("images", "button", "dash_dot_dot_line.png")}

    valueChanged = QtCore.pyqtSignal(int)


class MarkerCombo(CustomCombo):
    """Custom QComboBox for select matplotlib line marker"""
    key = ["None", ".", ",", "o", "v", "^", "<", ">", "1", "2", "3", "4", "8",
           "s", "p", "*", "h", "H", "+", "x", "D", "d", "|", "_"]
    text = {"None": "", ".": "point", ",": "pixel", "o": "circle",
            "v": "triangle_down", "^": "triangle_up", "<": "triangle_left",
            ">": "triangle_right", "1": "tri_down", "2": "tri_up",
            "3": "tri_left", "4": "tri_right", "8": "octagon", "s": "square",
            "p": "pentagon", "*": "star", "h": "hexagon1", "H": "hexagon2",
            "+": "plus", "x": "x", "D": "diamond", "d": "thin_diamond",
            "|": "vline", "_": "hline"}

    def _populate(self):
        for key in self.key:
            txt = self.text[key]
            if txt:
                image = os.environ["pychemqt"] + \
                    os.path.join("images", "marker", f"{txt}.png")
                self.addItem(QtGui.QIcon(QtGui.QPixmap(image)), self.text[key])
            else:
                self.addItem(self.text[key])


class NumericFactor(QtWidgets.QDialog):
    """Numeric format configuration dialog"""
    def __init__(self, config, unit=None, order=0, parent=None):
        super().__init__(parent)
        self.setWindowTitle(
            tr("pychemqt", "Format"))
        layout = QtWidgets.QGridLayout(self)
        self.checkFixed = QtWidgets.QRadioButton(
            tr(
                "pychemqt", "Fixed decimal point"))
        layout.addWidget(self.checkFixed, 1, 1, 1, 3)
        layout.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Total digits")), 2, 2)
        self.TotalDigits = Entrada_con_unidades(
            int, width=45, value=0, boton=False, spinbox=True, min=0, max=12,
            showNull=True)
        layout.addWidget(self.TotalDigits, 2, 3)
        layout.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Decimal digits")), 3, 2)
        self.DecimalDigits = Entrada_con_unidades(
            int, width=45, value=4, boton=False, spinbox=True, min=1, max=12)
        layout.addWidget(self.DecimalDigits, 3, 3)
        self.checkSignificant = QtWidgets.QRadioButton(
            tr(
                "pychemqt", "Significant figures"))
        layout.addWidget(self.checkSignificant, 4, 1, 1, 3)
        layout.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Figures")), 5, 2)
        self.FiguresSignificatives = Entrada_con_unidades(
            int, width=45, value=5, boton=False, spinbox=True, min=1, max=12)
        layout.addWidget(self.FiguresSignificatives, 5, 3)
        self.checkExp = QtWidgets.QRadioButton(
            tr(
                "pychemqt", "Exponential preferred"))
        layout.addWidget(self.checkExp, 6, 1, 1, 3)
        layout.addWidget(QtWidgets.QLabel(tr(
            "pychemqt", "Figures")), 7, 2)
        self.FiguresExponential = Entrada_con_unidades(
            int, width=45, value=5, boton=False, spinbox=True, min=1, max=12)
        layout.addWidget(self.FiguresExponential, 7, 3)
        layout.addItem(QtWidgets.QSpacerItem(
            30, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 8, 1)
        self.checkExpVariable = QtWidgets.QCheckBox(
            tr(
                "pychemqt", "Exponential for big/small values"))
        layout.addWidget(self.checkExpVariable, 9, 1, 1, 3)
        self.labelTolerancia = QtWidgets.QLabel(
            tr("pychemqt", "Tolerance"))
        layout.addWidget(self.labelTolerancia, 10, 2)
        self.Tolerance = Entrada_con_unidades(
            int, width=45, value=4, boton=False, spinbox=True, min=0, max=12)
        layout.addWidget(self.Tolerance, 10, 3)
        self.checkSign = QtWidgets.QCheckBox(tr(
            "pychemqt", "Show sign in positive values"))
        layout.addWidget(self.checkSign, 11, 1, 1, 3)
        self.checkThousand = QtWidgets.QCheckBox(
            tr(
                "pychemqt", "Show thousand separator"))
        layout.addWidget(self.checkThousand, 12, 1, 1, 3)

        self.checkFixed.toggled.connect(self.TotalDigits.setNotReadOnly)
        self.checkFixed.toggled.connect(self.DecimalDigits.setNotReadOnly)
        self.checkSignificant.toggled.connect(
            self.FiguresSignificatives.setNotReadOnly)
        self.checkExp.toggled.connect(self.ExpToggled)
        self.checkExp.toggled.connect(self.FiguresExponential.setNotReadOnly)
        self.checkExpVariable.toggled.connect(self.Tolerance.setNotReadOnly)
        self.checkExpVariable.toggled.connect(self.labelTolerancia.setEnabled)

        layout.addItem(QtWidgets.QSpacerItem(
            20, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 13, 1)
        self.muestra = QtWidgets.QLabel()
        layout.addWidget(self.muestra, 14, 1, 1, 3)

        buttonBox = QtWidgets.QDialogButtonBox()
        buttonBox.setStandardButtons(
            QtWidgets.QDialogButtonBox.StandardButton.Cancel
            | QtWidgets.QDialogButtonBox.StandardButton.Ok)
        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)
        layout.addWidget(buttonBox, 20, 1, 1, 3)

        self.checkFixed.setChecked(config["fmt"] == 0)
        self.TotalDigits.setNotReadOnly(config["fmt"] == 0)
        self.DecimalDigits.setNotReadOnly(config["fmt"] == 0)
        self.checkSignificant.setChecked(config["fmt"] == 1)
        self.FiguresSignificatives.setNotReadOnly(config["fmt"] == 1)
        self.checkExp.setChecked(config["fmt"] == 2)
        self.FiguresExponential.setNotReadOnly(config["fmt"] == 2)
        if config["fmt"] == 0:
            self.DecimalDigits.setValue(config["decimales"])
        elif config["fmt"] == 1:
            self.FiguresSignificatives.setValue(config["decimales"])
        else:
            self.FiguresExponential.setValue(config["decimales"])
        if "total" in config:
            self.TotalDigits.setValue(config["total"])
        if "exp" in config:
            self.checkExpVariable.setChecked(config["exp"])
        if "tol" in config:
            self.Tolerance.setValue(config["tol"])
        self.Tolerance.setNotReadOnly(config.get("exp", False))
        if "signo" in config:
            self.checkSign.setChecked(config["signo"])
        if "thousand" in config:
            self.checkThousand.setChecked(config["thousand"])

        self.updateMuestra()
        self.checkFixed.toggled.connect(self.updateMuestra)
        self.checkSignificant.toggled.connect(self.updateMuestra)
        self.checkExp.toggled.connect(self.updateMuestra)
        self.checkExpVariable.toggled.connect(self.updateMuestra)
        self.TotalDigits.valueChanged.connect(self.updateMuestra)
        self.DecimalDigits.valueChanged.connect(self.updateMuestra)
        self.FiguresSignificatives.valueChanged.connect(self.updateMuestra)
        self.FiguresExponential.valueChanged.connect(self.updateMuestra)
        self.Tolerance.valueChanged.connect(self.updateMuestra)
        self.checkSign.toggled.connect(self.updateMuestra)
        self.checkThousand.toggled.connect(self.updateMuestra)

        if unit and unit.__text__:
            layout.addItem(QtWidgets.QSpacerItem(
                20, 10, QtWidgets.QSizePolicy.Policy.Fixed,
                QtWidgets.QSizePolicy.Policy.Fixed), 15, 1, 1, 3)
            self.muestra = QtWidgets.QLabel()
            layout.addWidget(QtWidgets.QLabel(tr(
                "pychemqt", "Convert units")), 16, 1)
            self.unit = QtWidgets.QComboBox()
            for txt in unit.__text__:
                self.unit.addItem(txt)
            self.unit.setCurrentIndex(order)
            layout.addWidget(self.unit, 16, 2, 1, 2)

    def ExpToggled(self, boolean):
        """Exponential checkbox toggled"""
        self.FiguresExponential.setNotReadOnly(boolean)
        self.checkExpVariable.setDisabled(boolean)
        if self.checkExpVariable.isChecked():
            self.labelTolerancia.setDisabled(False)
            self.Tolerance.setReadOnly(True)

    def updateMuestra(self):
        """Update numeric sample"""
        kwargs = self.args()
        txt = tr("pychemqt", "Sample") + ": " + \
            representacion(pi*1e4, **kwargs)
        self.muestra.setText(txt)

    def args(self):
        """Return a dict with the selected numeric representation parameters"""
        kwarg = {}
        if self.checkFixed.isChecked():
            kwarg["fmt"] = 0
            kwarg["total"] = self.TotalDigits.value
            kwarg["decimales"] = self.DecimalDigits.value
        elif self.checkSignificant.isChecked():
            kwarg["fmt"] = 1
            kwarg["decimales"] = self.FiguresSignificatives.value
        else:
            kwarg["fmt"] = 2
            kwarg["decimales"] = self.FiguresExponential.value
        kwarg["eng"] = self.checkExpVariable.isEnabled() and \
            self.checkExpVariable.isChecked()
        kwarg["tol"] = self.Tolerance.value
        kwarg["signo"] = self.checkSign.isChecked()
        kwarg["thousand"] = self.checkThousand.isChecked()
        return kwarg


class InputFont(QtWidgets.QWidget):
    """Custom widget to edit a text input with font and color support"""
    textChanged = QtCore.pyqtSignal("QString")
    fontChanged = QtCore.pyqtSignal("QFont")
    colorChanged = QtCore.pyqtSignal("QString")

    def __init__(self, text="", font=None, color="#000000", parent=None):
        """
        text: Initial txt for widget
        font: QFont instance to initialize widget
        color: Inicial color to widget, in code #rrggbb
        """
        super().__init__(parent)

        layout = QtWidgets.QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        self.lineEdit = QtWidgets.QLineEdit()
        self.lineEdit.setFixedHeight(24)
        layout.addWidget(self.lineEdit)
        self.fontButton = QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(
            os.path.join(IMAGE_PATH, "button", "font.png"))), "")
        self.fontButton.setAutoDefault(False)
        self.fontButton.setDefault(False)
        self.fontButton.setFixedSize(24, 24)
        self.fontButton.setIconSize(QtCore.QSize(24, 24))
        self.fontButton.clicked.connect(self.fontButtonClicked)
        layout.addWidget(self.fontButton)
        self.colorButton = QtWidgets.QToolButton()
        self.colorButton.setFixedSize(24, 24)
        self.colorButton.clicked.connect(self.colorButtonClicked)
        layout.addWidget(self.colorButton)

        self.font = font
        self.color = color
        self.setRGB(color)
        self.setText(text)
        self.lineEdit.textChanged.connect(self.textChanged.emit)

    def setText(self, txt):
        """Set text for QLineEdit"""
        self.txt = txt
        self.lineEdit.setText(txt)

    def setRGB(self, rgb):
        """Wrap method to set color with a #rrggbb code"""
        color = QtGui.QColor(rgb)
        if color.isValid():
            self.setColor(color)

    def setColor(self, color):
        """Set color for widget"""
        self.colorButton.setPalette(QtGui.QPalette(color))
        paleta = QtGui.QPalette()
        paleta.setColor(QtGui.QPalette.ColorRole.Text, color)
        self.lineEdit.setPalette(paleta)
        self.colorChanged.emit(color.name())
        self.color = color

    def setFont(self, font):
        """Set font for widget"""
        self.font = font
        self.lineEdit.setFont(font)
        self.fontChanged.emit(font)

    def colorButtonClicked(self):
        """Show QColorDialog to change the color"""
        dlg = QtWidgets.QColorDialog(self.color, self)
        if dlg.exec():
            self.setColor(dlg.currentColor())

    def fontButtonClicked(self):
        """Show QFontDialog to choose the font"""
        dlg = QtWidgets.QFontDialog(self.lineEdit.font())
        if dlg.exec():
            self.setFont(dlg.currentFont())


class Table_Graphics(QtWidgets.QWidget):
    """Custom widget to implement as popup in PFD when mouse over stream and
    equipment graphic item, to show the status of entity and the properties
    desired if availables"""
    def __init__(self, entity, idx, preferences, parent=None):
        super().__init__(parent)
        self.setWindowFlags(QtCore.Qt.WindowType.Popup)
        layout = QtWidgets.QVBoxLayout(self)
        if isinstance(entity, Corriente):
            title = f"Stream {idx}"
        else:
            title = f"Equipment {idx}"
        label = QtWidgets.QLabel(title)
        label.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(label)
        line = QtWidgets.QFrame()
        line.setFrameShape(QtWidgets.QFrame.Shape.HLine)
        line.setFrameShadow(QtWidgets.QFrame.Shadow.Sunken)
        layout.addWidget(line)
        if entity:
            if entity.status:
                textos = entity.popup(preferences)
                for txt, tooltip, j in textos:
                    label = QtWidgets.QLabel(txt)
                    label.setToolTip(tooltip)
                    if j:
                        label.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight)
                    layout.addWidget(label)
            else:
                layout.addWidget(QtWidgets.QLabel(entity.msg))
        else:
            layout.addWidget(QtWidgets.QLabel(
                tr("pychemqt", "Undefined")))


def createAction(text, **kw):
    """Create a QAction and a QToolButton if its requested

    Parameters
    ----------
    text : str
        Text to define action
    button : boolean
        Return too a QToolButton with same properties
    icon : str
        Local path of icon in pychemqt images folder
    slot : obj
        Function to use as slot when action is triggered
    shortcut : QtGui.QKeySequence
        Shortcut keyboard to associate with action
    tip : str
        Text to use as tooltip
    checkable : boolean
        Set the created action checkable
    parent : QtCore.QObject
        Define object parent of created action and button
    """
    slot = kw.get("slot", None)
    shortcut = kw.get("shortcut", None)
    icon = kw.get("icon", None)
    checkable = kw.get("checkable", False)
    button = kw.get("button", False)
    parent = kw.get("parent", None)

    # Use text as default tip if no avialable
    tip = kw.get("tip", text)

    # Define action and its properties as requested
    action = QtGui.QAction(text, parent)
    if icon:
        action.setIcon(QtGui.QIcon(IMAGE_PATH + icon))
    if shortcut:
        action.setShortcut(shortcut)
    action.setToolTip(tip)
    action.setStatusTip(tip)
    if slot:
        action.triggered.connect(slot)
    if checkable:
        action.setCheckable(True)

    # Create too the button if its requested
    if button:
        boton = DragButton(parent)

        boton.setIcon(QtGui.QIcon(IMAGE_PATH + icon))
        boton.setToolTip(tip)
        boton.setStatusTip(tip)
        if slot:
            boton.clicked.connect(slot)
        boton.setCheckable(checkable)
        boton.setIconSize(QtCore.QSize(36, 36))
        boton.setFixedSize(QtCore.QSize(36, 36))
        return action, boton

    return action


def okToContinue(parent, dirty, func, parameters):
    """Function to ask user if any unsaved change

    Parameters
    ----------
    parent : QtWidgets.QWidget
        Parent widget to dialog
    dirty: boolean
        Indicates if there are unsaved changes
    func: object
        function to run if user want to save changes
    parameters: dict
        parameter of func
    """
    if not dirty:
        return True

    dialog = QtWidgets.QMessageBox.question(
        parent,
        tr("pychemqt", "Unsaved changes"),
        tr("pychemqt", "Save unsaved changes?"),
        QtWidgets.QMessageBox.StandardButton.Yes
        | QtWidgets.QMessageBox.StandardButton.No
        | QtWidgets.QMessageBox.StandardButton.Cancel,
        QtWidgets.QMessageBox.StandardButton.Yes)

    if dialog == QtWidgets.QMessageBox.StandardButton.Cancel:
        return False

    if dialog == QtWidgets.QMessageBox.StandardButton.No:
        return True

    if dialog == QtWidgets.QMessageBox.StandardButton.Yes:
        func(*parameters)
        return True

    return None


def mathTex2QPixmap(mathTex, fs):
    """Convert a latex text to a QPixmap to display in any qt widget. Code
    snippet from https://stackoverflow.com/questions/32035251

    Parameters
    ----------
    mathTex : str
        Latex text of expression to show
    fs : int
        Font size of text to show

    Returns
    -------
    qpixmap : :class:`QtGui.QPixmap`
        QPixmap ready to show in any other qt widget
    """

    # set up a mpl figure instance
    fig = Figure()
    fig.patch.set_facecolor('none')
    fig.set_canvas(FigureCanvasAgg(fig))
    renderer = fig.canvas.get_renderer()

    # plot the mathTex expression
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')
    ax.patch.set_facecolor('none')
    t = ax.text(0, 0, mathTex, ha='left', va='bottom', fontsize=fs)

    # fit figure size to text artist
    fwidth, fheight = fig.get_size_inches()
    fig_bbox = fig.get_window_extent(renderer)

    text_bbox = t.get_window_extent(renderer)

    tight_fwidth = text_bbox.width * fwidth / fig_bbox.width
    tight_fheight = text_bbox.height * fheight / fig_bbox.height

    fig.set_size_inches(tight_fwidth, tight_fheight)

    # convert mpl figure to QPixmap
    buf, size = fig.canvas.print_to_buffer()
    qimage = QtGui.QImage(
        buf, size[0], size[1], QtGui.QImage.Format.Format_ARGB32)
    qpixmap = QtGui.QPixmap(qimage)

    return qpixmap


class QLabelMath(QtWidgets.QLabel):
    """Customized QLabel to show a pixmap with a mathematical formulae"""
    def __init__(self, tex, *args, fs=12, **kw):
        """
        Parameters
        ----------
        tex : str
            Latex code text of mathematical expresion to show
        fs : int
            Font size used in image
        """
        super().__init__(*args, **kw)
        self.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        self.setFrameStyle(QtWidgets.QFrame.Shadow.Plain)
        self.fs = fs
        if tex:
            pixmap = mathTex2QPixmap(tex, fs)
            self.setPixmap(pixmap)

    def setTex(self, tex):
        """Set the image to show from a latex text source"""
        pixmap = mathTex2QPixmap(tex, self.fs)
        self.setPixmap(pixmap)


def demo():

    from lib import unidades

    app = QtWidgets.QApplication(sys.argv)

    ui = QtWidgets.QDialog()
    layout = QtWidgets.QVBoxLayout(ui)

    w = Entrada_con_unidades(unidades.Pressure)
    layout.addWidget(w)
    w2 = ColorSelector(isAlpha=True)
    layout.addWidget(w2)
    w3 = PathConfig()
    layout.addWidget(w3)
    w4 = LineConfig("saturation", "Line Style")
    layout.addWidget(w4)
    w5 = InputFont(text="foo bar", color="#0000ff")
    layout.addWidget(w5)
    w6 = Tabla(columnas=1, filas=1, verticalHeaderModel="C", dinamica=True)
    layout.addWidget(w6)
    w7 = QLabelMath(r"$\frac{\pi r^2}{2\epsilon \mu}$", fs=15)
    layout.addWidget(w7)
    w8 = Status()
    layout.addWidget(w8)

    ui.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    demo()
