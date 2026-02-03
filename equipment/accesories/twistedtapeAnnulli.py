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


from math import pi, log, cos, atan

from tools.qt import QtWidgets, translate

from lib.unidades import Dimensionless, Area, Length
from lib.utilities import refDoc
from UI.widgets import Entrada_con_unidades
from equipment.accesories.gui import ToolGui


__doi__ = {
    1:
        {"autor": "Coetzee, H., Liebenberg, L., Meyer, J.P.",
         "title": "Heat Transfer and Pressure Drop Characteristics of Angled "
                  "Spiraling Tape Inserts in a Heat Exchanger Annulus",
         "ref": "Heat Transfer Engineering 24(6) (2003) 29-39",
         "doi": "10.1080/714044412"},
    2:
        {"autor": "Gupte, N.S., Date, A.W.",
         "title": "Friction and Heat Transfer Characteristics of Helical "
                  "Turbulent Air Flow in Annuli",
         "ref": "J. Heat Transfer 111(2) (1989) 337-344",
         "doi": "10.1115/1.3250682"},
    # 3:
    #     {"autor": "",
    #      "title": "",
    #      "ref": "",
    #      "doi": ""},
        }


# Friction factor correlations
@refDoc(__doi__, [1])
def f_twistedAnnulli_Coetzee(Re, D, H, orientation=0):
    """Calculate friction factor for a annulus section for a double pipe with
    a twisted-tape insert using the Coetzee correlation (2003).

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    D : float
        Width of tape, [m]
    H : float
        Tape pitch for twist of π radians (180º), [m]
    orientation : boolean
        Set flow orientation of flow in the annulus with the curvature of tape
          0(along flow), 1(against flow)

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    y = H/D

    if orientation:
        g1 = 0.3618*y**2 - 1.047*y + 0.9186                             # Eq 13
        g2 = -0.0669*y**2 + 0.1656*y - 0.2282                           # Eq 14

    else:
        g1 = 0.6283*y**2 - 1.6519*y + 1.1939                            # Eq 15
        g2 = -0.0797*y**2 + 0.1814*y - 0.2392                           # Eq 16

    # Eq 12
    f = g1*Re**g2
    return f


@refDoc(__doi__, [2])
def f_twistedAnnulli_Gupte(Re, Do, Di, H):
    """Calculate friction factor for a annulus section for a double pipe with
    a twisted-tape insert using the Gupte correlation (1989).

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    Do : float
        Outer diameter of annulii, [m]
    Di : float
        Inner diameter of annulii, [m]
    H : float
        Tape pitch for twist of π radians (180º), [m]

    Returns
    -------
    f : float
        Friction factor, [-]

    Examples
    --------
    Selected point from Table 2 in [1]_

    >>> st = f_twistedAnnulli_Gupte(1e4, 1, 0.61, 3)
    """
    phi0 = 2*pi
    ri = Di/2
    ro = Do/2
    y = H/(ro-ri)

    r = ri/ro                                                           # Eq 9
    tau = (r**0.686-r**2)/(r*(1-r**0.686))                              # Eq 11
    finf = 0.046/Re**0.2                                                # Eq 13
    Dh_ = 2*phi0*(1-r**2)/(phi0*(1+r)+2*(1-r))                          # Eq 15

    K1 = ((1+r)*phi0 + 2*(1-r)) / ((1+r*tau)*phi0 + (1+tau)*(1-r))      # Eq 16
    K2 = phi0*(1+r*tau) + (1-r)*(1+tau)                                 # Eq 17

    # Eq 18
    K3 = (pi**2/(4*y**2*(1-r)**2))/(phi0*(1+r**2*tau)+(1-r**3)*(1+tau)/3)

    K5 = (K1**0.5/(2*y*(1-r)))*(2/finf)**0.5/Re**2                      # Eq 20
    K6 = (5.5+2.5*log(Re*(finf/2)**0.5/2))/K1**0.5                      # Eq 21
    K7 = 1/Dh_*(finf/2)**0.5*Re/K1**0.5                                 # Eq 22

    # Eq 19
    K4 = K5*(K6*(890*(1/tau-1)+30*K7*(r/tau**0.5+1)) + 5180*(1-1/tau**0.5)
             + 281.4*K7*(r-1))

    # Eq 23
    f = Dh_ * ((0.5*K1*(K2+K3)*finf + pi*Dh_**2*K4)/(phi0*(1-r**2)))
    return f


# Heat Transfer coefficient correlations
@refDoc(__doi__, [1])
def Nu_twistedAnnulli_Coetzee(Re, Pr, D, H, mu=1, muW=1, orientation=0):
    """Calculate Nusselt number for a annulus section for a double pipe with
    a twisted-tape insert using the Coetzee correlation (2003).

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    Pr : float
        Prandtl number, [-]
    D : float
        Width of tape, [m]
    H : float
        Tape pitch for twist of π radians (180º), [m]
    mu : float, optional
        Bulk flow temperature viscosity, [Pa·s]
    muW : float, optional
        Wall flow temperature viscosity, [Pa·s]
    orientation: boolean
        Set flow orientation of flow in the annulus with the curvature of tape
          0(along flow), 1(against flow)

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    y = H/D

    if orientation:
        f1 = 2.256e-9*y**2 - 10.989e-9*y + 16.03e-9                     # Eq 6
        f2 = -21.04e-6*y**2 + 125.27e-6*y - 211.33e-6                   # Eq 7
        f3 = 0.449*y**2 - 2.329*y + 4.503                               # Eq 8

    else:
        f1 = -0.3969e-9*y**2 + 1.233e-9*y + 3.369e-9                    # Eq 9
        f2 = 4.866e-6*y**2 - 26.27e-6*y - 54.17e-6                      # Eq 10
        f3 = 0.449*y**2 - 1.969*y + 4.021                               # Eq 11

    # Eq 5
    Nu = 0.0726*Re**0.8*Pr**0.333*(mu/muW)**0.14 * (f1*Re**2 + f2*Re + f3)
    return Nu


@refDoc(__doi__, [2])
def Nu_twistedAnnulli_Gupte(Re, Pr, Do, Di, H, boundary=0):
    """Calculate nusselt number for a annulus section for a double pipe with
    a twisted-tape insert using the Gupte correlation (1989).

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    Pr : float
        Prandtl number, [-]
    Do : float
        Outer diameter of annulii, [m]
    Di : float
        Inner diameter of annulii, [m]
    H : float
        Tape pitch for twist of π radians (180º), [m]
    boundary : integer
        0 - inner
        1 - outer

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    phi0 = 2*pi
    ri = Di/2
    ro = Do/2
    r = ri/ro                                                           # Eq 9
    tau = (r**0.686-r**2)/(r*(1-r**0.686))                              # Eq 11
    y = H/(ro-ri)
    finf = 0.046/Re**0.2                                                # Eq 13
    Dh_ = 2*phi0*(1-r**2)/(phi0*(1+r)+2*(1-r))                          # Eq 15
    K1 = ((1+r)*phi0 + 2*(1-r)) / ((1+r*tau)*phi0 + (1+tau)*(1-r))      # Eq 16
    K8 = pi*Dh_**2/(y*(1-r**2)*phi0*K1*finf*Re**2)                      # Eq 39
    ro_ = Re/Dh_*(finf*K1/2)**0.5                                       # Eq 44
    ri_ = Re*r*K1/Dh_*(finf*tau/2)**0.5                                 # Eq 45
    alfai = atan(pi*ri/H)                                               # Eq 48
    alfao = atan(pi*ro/H)                                               # Eq 49
    PF = 9.24*((Pr/0.9)**0.75-1)                                        # Eq 46
    A = 25/Pr*(((1+5*Pr)/Pr) * (log(1+5*Pr)-1) + 1)                     # Eq 33

    if boundary == 0:
        # Only outer wall heated
        K12 = 5250*Pr**0.731 - (137.5*Pr+A)*ro_                         # Eq 43

        # Eq 38a
        Stro = (K1*finf/2/cos(alfao))**0.5/0.9/(PF+(2*cos(alfao)/K1/finf)**0.5)

        # Eq 52
        Nu = (Stro/(1+r)-K8*(450-30*ro_)) / (1/(1+r)-K8*K12*(2/K1*finf)**0.5)
    else:
        # Only inner wall heated
        K11 = ((137.5*Pr+A)*ri_ + 5250*Pr**0.731)/tau**0.5              # Eq 42

        # Eq 38b
        Stri = (K1*tau**finf/2/cos(alfai))**0.5/0.9 / \
            (PF+(2*cos(alfai)/tau/K1/finf)**0.5)

        # Eq 53
        Nu = (Stri*r/(1+r) + 0.5*K8*K1*finf*(450+30*ri_)) / \
            (r/(1+r) + K8*K11*(2/finf/K1)**0.5)

    return Nu*Re*Pr



class TwistedTapeAnnuli():
    """Twisted-tape insert specially for annulli section of double-pipe
    section."""

    TEXT_FRICTION = (
        "Coetzee (2003)",
        "Gupte (1989)")

    TEXT_HEAT = (
        "Coetzee (2003)",
        "Gupte (1989)")

    TEXT_ORIENTACION = (
        translate("equipment", "Along flow"),
        translate("equipment", "Against flow"))

    def __init__(self, H, D, delta):
        """
        Definition of twisted tape inserts

        Parameters
        ----------
        H : float
            Tape pitch for twist of π radians (180º), [m]
        D : float
            Internal diameter of tube, [m]
        delta : float
            Tape thickness, [m]
        """
        # Area tube without tape
        self.A = pi*D**2/4

        # Tape twist parameter
        self.y = Dimensionless(H/D)


class UI_TwistedTapeAnnuli(ToolGui):
    """Twisted-tape insert dialog"""

    title = translate("equipment", "Use twisted tape insert in annuli section")

    def loadUI(self):
        """Add widget"""
        lyt = self.layout()

        label = QtWidgets.QLabel(self.tr("Tape pitch"))
        label.setToolTip(self.tr("Tape pitch for twist of π radians (180º)"))
        lyt.addWidget(label, 2, 1)
        self.H = Entrada_con_unidades(Length)
        lyt.addWidget(self.H, 2, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Tape diameter")), 3, 1)
        self.Dt = Entrada_con_unidades(Length)
        lyt.addWidget(self.Dt, 3, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Tape thickness")), 4, 1)
        self.delta = Entrada_con_unidades(Length, "Thickness")
        lyt.addWidget(self.delta, 4, 2)

        self.angled = QtWidgets.QCheckBox(self.tr("Angled twisted-tape"))
        self.angled.toggled.connect(self.setEnableOrientation)
        lyt.addWidget(self.angled, 5, 1, 1, 2)

        lytH = QtWidgets.QHBoxLayout()
        self.labelOrientation = QtWidgets.QLabel(self.tr(
            "Direction of flow relative to the tape curvature"))
        lytH.addWidget(self.labelOrientation)
        self.orientation = QtWidgets.QComboBox()
        for method in TwistedTapeAnnuli.TEXT_ORIENTACION:
            self.orientation.addItem(method)
        lytH.addWidget(self.orientation)
        lyt.addLayout(lytH, 7, 1, 1, 2)

    def setEnabled(self, boolean):
        """Add logic to parent setEnabled for orientation option"""
        ToolGui.setEnabled(self, boolean)
        self.setEnableOrientation(boolean and self.angled.isChecked())

    def setEnableOrientation(self, boolean):
        """Change Enable/Disable state for orientation of twisted tape"""
        self.labelOrientation.setEnabled(boolean)
        self.orientation.setEnabled(boolean)


class Dialog(QtWidgets.QDialog):
    """Component list config dialog"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(self.tr("Twisted-tape insert"))
        layout = QtWidgets.QVBoxLayout(self)
        self.datos = UI_TwistedTapeAnnuli()
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
