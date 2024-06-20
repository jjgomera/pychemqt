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
# Module for UI binary interaction parameter viewer
###############################################################################


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
        self.eleccion = QtWidgets.QComboBox()
        lyt.addWidget(self.eleccion, 1, 2)
        lyt.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed), 1, 3)
        self.stacked = QtWidgets.QStackedWidget()
        lyt.addWidget(self.stacked, 2, 1, 1, 3)

        # Get component names to show in table header
        names = []
        for cmp in ids:
            databank.execute(
                "SELECT id, name FROM compuestos WHERE id==%i" % cmp)
            names.append("%4i - %s" % databank.fetchone())

        args = (len(ids), len(ids))
        kw = {"stretch": False, "readOnly": True, "horizontalHeader": names,
              "verticalHeaderLabels": names}

        title = {"WILSON": "Aij", "UNIQUAC": "ΔUij", "NRTL": "Gij"}

        # Iterate over the EoS available
        for EoS in EoSBIP:
            self.eleccion.addItem(EoS)
            k = Kij(ids, EoS)

            widget = QtWidgets.QWidget()
            lyt2 = QtWidgets.QVBoxLayout(widget)
            lyt2.addWidget(QtWidgets.QLabel(title.get(EoS, "Kij")))
            table1 = Tabla(*args, **kw)
            lyt2.addWidget(table1)

            # Special case for NRTL with two interaction parameters matrix
            if EoS == "NRTL":
                lyt2.addWidget(QtWidgets.QLabel("α"))
                table2 = Tabla(*args, **kw)
                lyt2.addWidget(table2)
                kij, aij = k
                table1.setData(kij)
                table2.setData(aij)
                table1.resizeColumnsToContents()
                table2.resizeColumnsToContents()

            else:
                table1.setData(k)
                table1.resizeColumnsToContents()

            self.stacked.addWidget(widget)

        self.eleccion.currentIndexChanged.connect(self.stacked.setCurrentIndex)
        button = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Cancel
            | QtWidgets.QDialogButtonBox.StandardButton.Ok)
        button.accepted.connect(self.accept)
        button.rejected.connect(self.reject)
        lyt.addWidget(button, 3, 1, 1, 3)

        self.eleccion.setCurrentIndex(EoSIndex)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = Ui_BIP([38, 44, 45], 4)
    Dialog.show()
    sys.exit(app.exec())
