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
    2:
        {"autor": "Uttarwar, S.B., Raja Rao, M.",
         "title": "Augmentation of Laminar Flow Heat Transfer in Tubes by "
                  "Means of Wire Coil Inserts",
         "ref": "J. Heat Transfer 107(4) (1985) 930-935",
         "doi": "10.1115/1.3247523"},
    3:
        {"autor": "Inaba, H., Ozaki, K., Kanaoka, S.",
         "title": "A Fundamental Study of Heat-Transfer Enhancement and "
                  "Flow-Drag Reduction in Tubes by Means of Wire Coil Insert",
         "ref": "Trans. Jpn. Soc. Mech. Eng. 60 (1994) 240-247",
         "doi": "10.1299/kikaib.60.240"},
    4:
        {"autor": "Naphon, P.",
         "title": "Effect of coil-wire insert on heat trasnfer enhancement "
                  "pressure drop of the horizontal concentric tubes",
         "ref": "Int. Comm. Heat Mass Transfer 33(6) (2006) 753-763",
         "doi": "10.1016/j.icheatmasstransfer.2006.01.020"},
    # 5:
    #     {"autor": "",
    #      "title": "",
    #      "ref": "",
    #      "doi": ""},
        }


# Friction factor correlations
@refDoc(__doi__, [1])
def f_wire_Garcia(Re, P, e):
    """Calculate friction factor for a pipe with a wire coil using the Garcia
    et al. correlation (2005).

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    P : float
        helical pitch for twist of 2π radians (360º), [m]
    e : float
        Wire diameter, [m]

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    # Eq 4
    f = 9.35*(P/e)**-1.16*Re**-0.217

    return f


@refDoc(__doi__, [3])
def f_wire_Inaba(Re, P, e):
    """Calculate friction factor for a pipe with a wire coil using the Inaba
    correlation (1994).

    Valid for P/e > 10, for bellow ratio the pipe is like a rough pipe

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    P : float
        helical pitch for twist of 2π radians (360º), [m]
    e : float
        Wire diameter, [m]

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    # Eq 7
    f = 11.5*Re**-0.39*(P/e)**-0.87

    return f


@refDoc(__doi__, [4])
def f_wire_Naphon(Re, P, D):
    """Calculate friction factor for a pipe with a wire coil using the Naphon
    correlation (2006).

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    P : float
        helical pitch for twist of 2π radians (360º), [m]
    D : float
        Internal diameter of tube, [m]

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    # Eq 7
    f = 322.92*log(Re)**-1.849*(P/D)**0.061

    return f


# Heat Transfer coefficient correlations
@refDoc(__doi__, [1])
def Nu_wire_Garcia(Re, Pr, P, e):
    """Calculate Nusselt number for a pipe with a wire coil using the Garcia
    et al. correlation (2005).

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    Pr : float
        Prandtl number, [-]
    P : float
        helical pitch for twist of 2π radians (360º), [m]
    e : float
        Wire diameter, [m]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    # Eq 10
    Nu = 0.132*(P/e)**-0.372*Re**0.72*Pr**0.37

    return Nu


@refDoc(__doi__, [2])
def Nu_wire_Uttarwar(Re, Pr, P, D, mu=None, muW=None):
    """Calculate Nusselt number for a pipe with a wire coil using the Uttarwar-
    Raja Rao correlation (1985).

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    Pr : float
        Prandtl number, [-]
    P : float
        helical pitch for twist of 2π radians (360º), [m]
    e : float
        Wire diameter, [m]
    D : float
        Internal diameter of tube, [m]
    mu : float
        Bulk flow temperature viscosity, [Pa·s]
    muW : float
        Wall flow temperature viscosity, [Pa·s]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    # Helix angle
    tan_alpha = P/pi/D

    # Eq 8
    Nu = 1.65*tan_alpha * Re**(0.25*tan_alpha**-0.38)*Pr**0.35

    if mu and muW:
        Nu *= (mu/muW)**0.14

    return Nu


@refDoc(__doi__, [3])
def Nu_wire_Inaba(Re, Pr, P, e):
    """Calculate Nusselt number for a pipe with a wire coil using the Inaba
    correlation (1994).

    Valid for P/e > 10, for bellow ratio the pipe is like a rough pipe

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    Pr : float
        Prandtl number, [-]
    P : float
        helical pitch for twist of 2π radians (360º), [m]
    e : float
        Wire diameter, [m]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    if Re < 2000:
        # Eq 11
        Nu = 0.225*Re**0.8*Pr**(1/3)*(P/e)**-0.48
    else:
        # Eq 10
        Nu = 0.803*Re**0.63*Pr**(1/3)*(P/e)**-0.48

    return Nu


@refDoc(__doi__, [4])
def Nu_wire_Naphon(Re, Pr, P, D):
    """Calculate Nusselt number for a pipe with a wire coil using the Naphon
    correlation (2006).

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    Pr : float
        Prandtl number, [-]
    P : float
        helical pitch for twist of 2π radians (360º), [m]
    D : float
        Internal diameter of tube, [m]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    if Re < 5000:
        raise NotImplementedError("Input out of bound")

    # Eq 6
    Nu = 0.156*Re**0.512*Pr**(1/3)*(P/D)**0.253

    return Nu


class WireCoil(CallableEntity):
    """Wire coil insert for pipe to improve heat transfer

    Parameters
    ----------
    P : float
        helical pitch for twist of 2π radians (360º), [m]
    e : float
        Wire diameter, [m]
    """
    TEXT_FRICTION = (
        "García (2005)",
        "Inaba (1994)",
        "Naphon (2006)"
        )

    TEXT_HEAT = (
        "García (2005)",
        "Uttarwar-Raja Rao (1985)",
        "Inaba (1994)",
        "Naphon (2006)"
        )

    status = 0
    msg = ""
    kw = {
        "methodFriction": 0,
        "methodHeat": 0,
        "P": 0,
        "e": 0}

    valueChanged = QtCore.pyqtSignal(object)
    inputChanged = QtCore.pyqtSignal(object)

    @property
    def isCalculable(self):
        """Check if all input are defined"""
        if not self.kw["P"]:
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
        self.P = self.kw["P"]

        self.valueChanged.emit(self)

    def Nu(self, Re, Pr, D, mu, muW):
        """Calculate nusselt number"""
        if self.kw["methodHeat"] == 1:
            # Uttarwar-Raja Rao (1985)
            Nu = Nu_wire_Uttarwar(Re, Pr, self.P, D, mu, muW)

        elif self.kw["methodHeat"] == 2:
            # Inaba (1994)
            Nu = Nu_wire_Inaba(Re, Pr, self.P, self.e)

        elif self.kw["methodHeat"] == 3:
            # Naphon (2006)
            Nu = Nu_wire_Naphon(Re, Pr, self.P, D)

        else:
            # García (2005)
            Nu = Nu_wire_Garcia(Re, Pr, self.P, self.e)

        return Nu

    def f(self, Re, D):
        """Calculate friction factor"""
        if self.kw["methodFriction"] == 1:
            # Inaba (1994)
            f = f_wire_Inaba(Re, self.P, self.e)

        elif self.kw["methodHeat"] == 2:
            # Naphon (2006)
            f = f_wire_Naphon(Re, self.P, D)

        else:
            # García (2005)
            f = f_wire_Garcia(Re, self.P, self.e)

        return f


class UI_WireCoil(ToolGui):
    """Wire-coil insert dialog"""

    title = translate("equipment", "Use wire coil insert")

    def loadUI(self):
        """Add widget"""
        self.Entity = WireCoil()

        lyt = self.layout()

        lytH = QtWidgets.QHBoxLayout()
        lytH.addWidget(QtWidgets.QLabel(
            self.tr("Friction factor calculation method")))
        self.methodFriction = QtWidgets.QComboBox()
        for method in WireCoil.TEXT_FRICTION:
            self.methodFriction.addItem(method)
        self.methodFriction.currentIndexChanged.connect(
            partial(self.changeParams, "methodFriction"))
        lytH.addWidget(self.methodFriction)
        lyt.addLayout(lytH, 2, 1, 1, 2)

        lytH = QtWidgets.QHBoxLayout()
        lytH.addWidget(QtWidgets.QLabel(
            self.tr("Heat transfer calculation method")))
        self.methodHeat = QtWidgets.QComboBox()
        for method in WireCoil.TEXT_HEAT:
            self.methodHeat.addItem(method)
        self.methodHeat.currentIndexChanged.connect(
            partial(self.changeParams, "methodHeat"))
        lytH.addWidget(self.methodHeat)
        lyt.addLayout(lytH, 3, 1, 1, 2)


        label = QtWidgets.QLabel(self.tr("Wire pitch"))
        label.setToolTip(self.tr("Wire pitch for twist of 2π radians (360º)"))
        lyt.addWidget(label, 4, 1)
        self.P = Entrada_con_unidades(Length)
        self.P.valueChanged.connect(partial(self.changeParams, "P"))
        lyt.addWidget(self.P, 4, 2)
        label = QtWidgets.QLabel("Wire diameter")
        lyt.addWidget(label, 5, 1)
        self.e = Entrada_con_unidades(Length, "Thickness")
        self.e.valueChanged.connect(partial(self.changeParams, "e"))
        lyt.addWidget(self.e, 5, 2)

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
