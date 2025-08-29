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


Module for UI binary interaction parameter viewer

* :class:`Ui_BIP`: Dialog to view the BIP for selected component and EoS

'''


from tools.qt import QtWidgets

from lib.bip import Kij, EoSBIP
from lib.sql import databank
from UI.widgets import Tabla


class Ui_BIP(QtWidgets.QDialog):
    """Dialog to view the BIP between selected component for EoS available

    Parameters
    ----------
    ids : list
        Index of component
    EoSIndex: integer
        Index of EoS with BIP to show
    """
    def __init__(self, ids, EoSIndex=0, parent=None):
        """Constructor"""
        super().__init__(parent)
        self.setWindowTitle(self.tr("BIP (Binary interaction parameters)"))

        lyt = QtWidgets.QGridLayout(self)
        lyt.addWidget(QtWidgets.QLabel("EoS"), 1, 1)
        self.eos = QtWidgets.QComboBox()
        lyt.addWidget(self.eos, 1, 2)
        lyt.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed), 1, 3)
        self.stacked = QtWidgets.QStackedWidget()
        lyt.addWidget(self.stacked, 2, 1, 1, 3)

        # Get component names to show in table header
        names = []
        for cmp in ids:
            databank.execute(
                f"SELECT id, name FROM compuestos WHERE id=={cmp}")
            names.append("%4i - %s" % databank.fetchone())

        kw = {"stretch": False, "readOnly": True, "horizontalHeader": names,
              "verticalHeaderLabels": names}

        title = {"WILSON": "Aij", "UNIQUAC": "ΔUij", "NRTL": "Gij"}

        # Iterate over the EoS available
        for EoS in EoSBIP:
            self.eos.addItem(EoS)
            k = Kij(ids, EoS)

            widget = QtWidgets.QWidget()
            lyt2 = QtWidgets.QVBoxLayout(widget)
            lyt2.addWidget(QtWidgets.QLabel(title.get(EoS, "Kij")))
            table1 = Tabla(len(ids), len(ids), **kw)
            lyt2.addWidget(table1)

            # Special case for NRTL with two interaction parameters matrix
            if EoS == "NRTL":
                lyt2.addItem(QtWidgets.QSpacerItem(
                    20, 20, QtWidgets.QSizePolicy.Policy.Fixed,
                    QtWidgets.QSizePolicy.Policy.Fixed))
                lyt2.addWidget(QtWidgets.QLabel("α"))
                table2 = Tabla(len(ids), len(ids), **kw)
                lyt2.addWidget(table2)
                table1.setData(k[0])
                table2.setData(k[1])
                table2.resizeColumnsToContents()

            else:
                table1.setData(k)

            table1.resizeColumnsToContents()

            lyt2.addItem(QtWidgets.QSpacerItem(
                0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
                QtWidgets.QSizePolicy.Policy.Expanding))

            width = table1.verticalHeader().sizeHint().width() + 2
            for i in range(table1.columnCount()):
                width += table1.columnWidth(i)
            table1.setFixedWidth(width)

            height = table1.horizontalHeader().sizeHint().height() + 2
            for i in range(table1.rowCount()):
                height += table1.rowHeight(i)
            table1.setFixedHeight(height)

            if EoS == "NRTL":
                table2.setFixedWidth(width)
                table2.setFixedHeight(height)

            self.stacked.addWidget(widget)

        self.eos.currentIndexChanged.connect(self.stacked.setCurrentIndex)
        button = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Cancel
            | QtWidgets.QDialogButtonBox.StandardButton.Ok)
        button.accepted.connect(self.accept)
        button.rejected.connect(self.reject)
        lyt.addWidget(button, 4, 1, 1, 4)

        self.eos.setCurrentIndex(EoSIndex)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = Ui_BIP([38, 44, 45], 4)
    Dialog.show()
    sys.exit(app.exec())
