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


from functools import partial
from math import pi, log, cos, atan

from tools.qt import QtCore, QtWidgets, translate

from equipment.widget.gui import ToolGui, CallableEntity
from lib.unidades import Length
from lib.utilities import refDoc
from UI.widgets import Entrada_con_unidades


__doi__ = {
    1:
        {"autor": "García, A., Vicente, P.G., Viedma, A.",
         "title": "Experimental study of heat transfer enhancement with "
                  "wire coil inserts in laminar-transition-turbulent regimes "
                  "at different Prandtl numbers",
         "ref": "Int. J. Heat Mass Transfer 48(21-22) (2005) 4640-4651",
         "doi": "10.1016/j.ijheatmasstransfer.2005.04.024"},
    # 2:
    #     {"autor": "",
    #      "title": "",
    #      "ref": "",
    #      "doi": ""}
        }


@refDoc(__doi__, [1])
def f_wire_Garcia(Re, p, e):
    """Calculate friction factor for a pipe with a wire coil using the Garcia
    et al. correlation (2005).

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    p : float
        helical pitch for twist of π radians (180º), [m]
    e : float
        Wire diameter, [m]

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    # Eq 4
    f = 9.35*(p/e)**-1.16*Re**-0.217

    return f


@refDoc(__doi__, [1])
def Nu_wire_Garcia(Re, Pr, p, e):
    """Calculate Nusselt number for a pipe with a wire coil using the Garcia
    et al. correlation (2005).

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    Pr : float
        Prandtl number, [-]
    p : float
        helical pitch for twist of π radians (180º), [m]
    e : float
        Wire diameter, [m]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    # Eq 10
    Nu = 0.132*(p/e)**-0.372*Re**0.72*Pr**0.37

    return Nu


class WireCoil(CallableEntity):
    """Wire coil insert for pipe to improve heat transfer

    Parameters
    ----------
    p : float
        helical pitch for twist of π radians (180º), [m]
    e : float
        Wire diameter, [m]
    """
    status = 0
    msg = ""
    kw = {
        "p": 0,
        "e": 0}

    valueChanged = QtCore.pyqtSignal(object)
    inputChanged = QtCore.pyqtSignal(object)

    @property
    def isCalculable(self):
        """Check if all input are defined"""
        if not self.kw["p"]:
            self.msg = translate("equipment", "undefined wire coil pitch")
            self.status = 0
            return False
        if not self.kw["e"]:
            self.msg = translate("equipment", "undefined wire coil diameter")
            self.status = 0
            return False

        self.msg = ""
        self.status = 1
        return True

    def calculo(self):
        """Definition of twisted tape inserts for annuli sections"""
        self.e = self.kw["e"]
        self.p = self.kw["p"]

        self.valueChanged.emit(self)

    def Nu(self, Re, Pr):
        """Calculate nusselt number"""
        Nu = Nu_wire_Garcia(Re, Pr, self.p, self.e)
        return Nu

    def f(self, Re):
        """Calculate friction factor"""
        f= f_wire_Garcia(Re, self.p, self.e)
        return f


class UI_WireCoil(ToolGui):
    """Wire-coil insert dialog"""

    title = translate("equipment", "Use wire coil insert")

    def loadUI(self):
        """Add widget"""
        self.Entity = WireCoil()

        lyt = self.layout()

        label = QtWidgets.QLabel(self.tr("Wire pitch"))
        label.setToolTip(self.tr("Wire pitch for twist of π radians (180º)"))
        lyt.addWidget(label, 2, 1)
        self.p = Entrada_con_unidades(Length)
        self.p.valueChanged.connect(partial(self.changeParams, "p"))
        lyt.addWidget(self.p, 2, 2)
        label = QtWidgets.QLabel("Wire diameter")
        lyt.addWidget(label, 3, 1)
        self.e = Entrada_con_unidades(Length, "Thickness")
        self.e.valueChanged.connect(partial(self.changeParams, "e"))
        lyt.addWidget(self.e, 3, 2)

        self.Entity.valueChanged.connect(self.valueChanged.emit)
        self.Entity.inputChanged.connect(self.populate)


class Dialog(QtWidgets.QDialog):
    """Component list config dialog"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(self.tr("Twisted-tape insert"))
        layout = QtWidgets.QVBoxLayout(self)
        self.datos = UI_WireCoil()
        layout.addWidget(self.datos)
        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Cancel
            | QtWidgets.QDialogButtonBox.StandardButton.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = Dialog()
    Dialog.show()
    sys.exit(app.exec())
