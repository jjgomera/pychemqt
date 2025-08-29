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


Module to common data entry

* :class:`InputTableWidget`: Widget for table data
* :class:`InputTableDialog`: Dialog for table data
* :class:`eqDIPPR`: Widget for select DIPPR equation

'''


from functools import partial
import os

from numpy import loadtxt
from tools.qt import QtCore, QtGui, QtWidgets

from lib import unidades
from UI.widgets import Entrada_con_unidades, Tabla


class eqDIPPR(QtWidgets.QWidget):
    """Custom widget to define DIPPR equation input"""
    def __init__(self, value, parent=None):
        super().__init__(parent)
        layout = QtWidgets.QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(QtWidgets.QLabel(self.tr("Eq DIPPR") + " "))
        self.eqDIPPR = QtWidgets.QSpinBox()
        self.eqDIPPR.setValue(value)
        self.eqDIPPR.setRange(1, 9)
        self.eqDIPPR.setAlignment(
            QtCore.Qt.AlignmentFlag.AlignRight
            | QtCore.Qt.AlignmentFlag.AlignVCenter)
        self.eqDIPPR.setFixedWidth(50)
        txt = self.tr("Equation") + "\n"
        txt += "\t1:\tY = A+B*T+C*T²+D*T³+E*T⁴\n"
        txt += "\t2:\tY = exp(A+B*T+C*ln(T)+D*T^E)\n"
        txt += "\t3:\tY = A*T^(B/(1+C*T+D*T^2))\n"
        txt += "\t4:\tY = A+B*exp(-C/T^D)\n"
        txt += "\t5:\tY = A + B/T + C/T³ + D/T⁸ + E/T⁹\n"
        txt += "\t6:\tY = A/(B^(1+(1-T/C)^D)\n"
        txt += "\t7:\tY = A*(1-Tr)^(B+C*Tr+D*Tr²+E*Tr³)\n"
        txt += "\t8:\tY = A+ B*((C/T)/sinh(C/T))² + D*((E/T)/cosh(E/T))²\n"
        txt += "\t9:\tY = A²/Tr+B-2ACTr-ADTr²-C²Tr³/3-CDTr⁴/2-D²Tr⁵/5\n"
        txt += self.tr("where") + ":\n"
        txt += "\t" + self.tr("Y Property to fit") + "\n"
        txt += "\t" + self.tr("T temperature in Kelvin") + "\n"
        txt += "\t" + self.tr("Tr: reduced temperature T/Tc") + "\n"
        txt += "\t" + self.tr("A,B,C,D,E parameters")
        self.eqDIPPR.setToolTip(txt)
        layout.addWidget(self.eqDIPPR)
        layout.addStretch()

    @property
    def value(self):
        """Return value of widget"""
        return self.eqDIPPR.value

    def setValue(self, value):
        """Set value property of widget"""
        self.eqDIPPR.setValue(value)

    def clear(self):
        """Clear widget value"""
        self.eqDIPPR.clear()


class InputTableWidget(QtWidgets.QWidget):
    """Table data input dialog"""
    def __init__(self, columns, data=None, t=None, prop=None,
                 horizontalHeader=None, title="", DIPPR=False, hasTc=0,
                 Tc=None, eq=1, unit=None, parent=None):
        """
        data: mrray with original data
        t: values for x column, generally temperature
        prop: values for 2...n columns
        horizontalHeader: List with column title
        DIPPR: boolean to show DIPPR widget
        hasTc: boolean to show critical temperature (some DIPPR eq need it)
        Tc: value for critical temperature
        eq: Value for DIPPR equation
        unit: List of unidades classes for column definition
        """
        super().__init__(parent)
        self.columns = columns
        self.title = title
        self.unit = unit
        gridLayout = QtWidgets.QGridLayout(self)
        gridLayout.setContentsMargins(0, 0, 0, 0)
        openButton = QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"]+"/images/button/fileOpen.png")), "")
        openButton.setToolTip(self.tr("Load data from a file"))
        openButton.clicked.connect(self.open)
        gridLayout.addWidget(openButton, 1, 1)
        saveButton = QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"]+"/images/button/fileSave.png")), "")
        saveButton.setToolTip(self.tr("Save data to a file"))
        saveButton.clicked.connect(self.save)
        gridLayout.addWidget(saveButton, 1, 2)
        clearButton = QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"]+"/images/button/clear.png")), "")
        clearButton.setToolTip(self.tr("Clear data"))
        clearButton.clicked.connect(self.delete)
        gridLayout.addWidget(clearButton, 1, 3)
        gridLayout.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed), 1, 4)

        self.tabla = Tabla(self.columns, horizontalHeader=horizontalHeader,
                           verticalHeader=False, stretch=False)
        self.tabla.setConnected()
        if unit:
            hHeader = []
            for un, head in zip(self.unit, horizontalHeader):
                hHeader.append(f"{head}, {un.text()}")
            self.tabla.setHorizontalHeaderLabels(hHeader)
            self.tabla.horizontalHeader().sectionClicked.connect(self.editUnit)

        if data:
            self.tabla.setData(data)
            self.tabla.addRow()
        elif t and prop:
            self.tabla.setColumn(0, t)
            self.tabla.setColumn(1, prop)
        gridLayout.addWidget(self.tabla, 2, 1, 1, 4)

        if DIPPR:
            self.eqDIPPR = eqDIPPR(eq)
            gridLayout.addWidget(self.eqDIPPR, 3, 1, 1, 4)
            self.eqDIPPR.eqDIPPR.valueChanged.connect(self.showTc)

            self.labelTc = QtWidgets.QLabel("Tc: ", self)
            gridLayout.addWidget(self.labelTc, 4, 1)
            self.tc = Entrada_con_unidades(unidades.Temperature, value=Tc)
            gridLayout.addWidget(self.tc, 4, 2, 1, 3)
            self.showTc(1)

    def showTc(self, value):
        """Show/hide Tc widget"""
        self.labelTc.setVisible(value in (7, 9))
        self.tc.setVisible(value in (7, 9))

    def open(self):
        """Load data from a test file"""
        fname, ext = QtWidgets.QFileDialog.getOpenFileName(
            self, self.tr("Open text file"), "./")
        if fname:
            try:
                # Numpy raise error if use the fname directly and find a
                # non-latin1 character, inclusive in comment lines
                with open(fname, "rb") as file:
                    data = loadtxt(file)
                self.delete()
                self.tabla.setData(data)
            except ValueError as er:
                # Raise a error msg if the file can load by loadtxt, the user
                # can select any type of file and the input error is possible
                title = self.tr("Failed to load file")
                msg = fname + "\n" + er.args[0]
                QtWidgets.QMessageBox.critical(self, title, msg)

    def save(self):
        """Save currend data of table to a file"""
        fname, fileFilter = QtWidgets.QFileDialog.getSaveFileName(
            self, self.tr("Save data to file"), "./")
        if fname:
            with open(fname, 'w') as file:
                file.write("#"+self.title+"\n")
                file.write("#")
                for i in range(self.tabla.columnCount()):
                    item = self.tabla.horizontalHeaderItem(i)
                    file.write(item.text()+"\t")
                file.write("\n")
                data = self.data
                for fila in data:
                    for dat in fila:
                        file.write(str(dat)+"\t")
                    file.write("\n")

    def delete(self):
        """Clear table"""
        self.tabla.setRowCount(0)
        self.tabla.clearContents()
        self.tabla.addRow()

    @property
    def data(self):
        """Table data"""
        return self.tabla.getData()

    def column(self, column, magnitud=None, unit="conf"):
        """
        column: column to get
        magnitud: magnitud to get the values
        unit: unit of the values in table"""
        data = self.tabla.getColumn(column)
        if self.unit:
            magnitud = self.unit[column]
            tx = self.tabla.horizontalHeaderItem(column).text().split(", ")[-1]
            unit = magnitud.__units__[magnitud.__text__.index(tx)]

        if magnitud is not None:
            data = [magnitud(x, unit) for x in data]
        return data

    def editUnit(self, col):
        """Show dialog to config input unit"""
        unit = self.unit[col]
        widget = QtWidgets.QComboBox(self.tabla)
        for txt in unit.__text__:
            widget.addItem(txt)
        txt = self.tabla.horizontalHeaderItem(col).text().split(", ")[-1]
        widget.setCurrentText(txt)

        # Define Unit combobox geometry
        size = self.tabla.horizontalHeader().sectionSize(col)
        pos = self.tabla.horizontalHeader().sectionPosition(col)
        h = self.tabla.horizontalHeader().height()
        geometry = QtCore.QRect(pos, 0, size, h)
        widget.setGeometry(geometry)
        widget.currentIndexChanged.connect(partial(self.updateHeader, col))
        widget.show()
        widget.showPopup()

    def updateHeader(self, col):
        """Change the text in header"""
        widget = self.sender()
        txt = self.tabla.horizontalHeaderItem(col).text()
        newtxt = f"{txt.split(",")[0]}, {widget.currentText()}"
        self.tabla.setHorizontalHeaderItem(
                col, QtWidgets.QTableWidgetItem(newtxt))
        widget.close()


class InputTableDialog(QtWidgets.QDialog):
    """Dialog to config thermal method calculations"""
    def __init__(self, columns=2, showHelp=False, helpFile="", **kwargs):
        """
        title: window title
        showHelp: boolean to show help button
        helpFile: Path for help file, file or url
        """
        parent = kwargs.get("parent", None)
        super().__init__(parent)
        title = kwargs.get("title", "")
        self.setWindowTitle(title)
        self.helpFile = helpFile
        layout = QtWidgets.QVBoxLayout(self)
        self.widget = InputTableWidget(columns, **kwargs)
        layout.addWidget(self.widget)
        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Cancel
            | QtWidgets.QDialogButtonBox.StandardButton.Ok)
        if showHelp:
            self.buttonBox.addButton(
                QtWidgets.QDialogButtonBox.StandardButton.Help)
            self.buttonBox.helpRequested.connect(self.ayuda)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox)

    def ayuda(self):
        """Show help file"""
        url = QtCore.QUrl(self.helpFile)
        QtGui.QDesktopServices.openUrl(url)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    header = ["T", "μ"]
    ui = InputTableDialog(
        title="titulo", horizontalHeader=header,
        unit=[unidades.Temperature, unidades.Viscosity], showHelp=True)
    ui.show()
    sys.exit(app.exec())
