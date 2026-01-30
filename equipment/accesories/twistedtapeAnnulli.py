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

from tools.qt import QtWidgets

from lib.unidades import Dimensionless, Area, Length
from lib.utilities import refDoc
from UI.widgets import Entrada_con_unidades


__doi__ = {
    1:
        {"autor": "Coetzee, H., Liebenberg, L., Meyer, J.P.",
         "title": "Heat Transfer and Pressure Drop Characteristics of Angled "
                  "Spiraling Tape Inserts in a Heat Exchanger Annulus",
         "ref": "Heat Transfer Engineering 24(6) (2003) 29-39",
         "doi": "10.1080/714044412"},
        }


# Friction factor correlations
@refDoc(__doi__, [1])
def f_twistedAnnulli_Coetzee(Re, D, H, opposite=0):
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
    opposite : boolean
        Set flow orientation of flow in the annulus with the curvature of tape
          0(along flow), 1(against flow)

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    y = H/D

    if opposite:
        g1 = 0.3618*y**2 - 1.047*y + 0.9186                             # Eq 13
        g2 = -0.0669*y**2 + 0.1656*y - 0.2282                           # Eq 14

    else:
        g1 = 0.6283*y**2 - 1.6519*y + 1.1939                            # Eq 15
        g2 = -0.0797*y**2 + 0.1814*y - 0.2392                           # Eq 16

    # Eq 12
    f = g1*Re**g2
    return f


# Heat Transfer coefficient correlations
@refDoc(__doi__, [1])
def Nu_twistedAnnulli_Coetzee(Re, Pr, D, H, mu=1, muW=1, opposite=0):
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
    opposite : boolean
        Set flow orientation of flow in the annulus with the curvature of tape
          0(along flow), 1(against flow)

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    y = H/D

    if opposite:
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

# class TwistedTape():
#     """Twisted-tape insert used in heat exchanger to improve efficiency.
#     This tape, generally a thin metal strip, is twisted about its longitudinal
#     axis"""

#     TEXT_FRICTION = (
#         "Manglik-Bergles (1993)",
#         "Plassis-Kröger (1984)",
#         "Lopina-Bergles (1969)",
#         "Shah-London (1978)")

#     TEXT_HEAT = (
#         "HTRI",
#         "Lopina-Bergles (1969)",
#         "Manglik-Bergles (1993)",
#         "Plessis-Kröger (1987)",
#         "Hong-Bergles (1976)")

#     def __init__(self, H, D, delta):
#         """
#         Definition of twisted tape accesory

#         Parameters
#         ----------
#         H : float
#             Tape pitch for twist of π radians (180º), [m]
#         D : float
#             Internal diameter of tube, [m]
#         delta : float
#             Tape thickness, [m]
#         """
#         # Geometrical definition of parameters in [1]_

#         # Helical factor, Eq 2
#         self.G = Dimensionless((1+pi**2/D**2/4/H**2)**0.5)

#         # Effective cross-sectional flow area, Eq 3
#         self.Ae = Area(2*H**2/pi*(self.G-1) - D*delta)

#         # Effective wetted perimeter, Eq 4
#         self.Pe = Length(2 * (D - delta + pi*D/2/self.G))

#         # Effective hydraulic diameter, Eq 7
#         self.De = Length(4*self.Ae/self.Pe)

#         # Area tube without tape
#         self.A = pi*D**2/4

#         # Tape twist parameter
#         self.y = Dimensionless(H/D)


# class UI_TwistedTape(QtWidgets.QWidget):

#     """Custom widget to define DIPPR equation input"""
#     def __init__(self, parent=None):
#         super().__init__(parent)
#         lyt = QtWidgets.QGridLayout(self)
#         self.check = QtWidgets.QCheckBox(self.tr("Use twisted tape insert"))
#         self.check.toggled.connect(self.setEnabled)
#         lyt.addWidget(self.check, 1, 1)
#         label = QtWidgets.QLabel(self.tr("Tape pitch"))
#         label.setToolTip(self.tr("Tape pitch for twist of π radians (180º)"))
#         lyt.addWidget(label, 2, 1)
#         self.H = Entrada_con_unidades(Length)
#         lyt.addWidget(self.H, 2, 2)
#         lyt.addWidget(QtWidgets.QLabel(self.tr("Tape diameter")), 3, 1)
#         self.Dt = Entrada_con_unidades(Length)
#         lyt.addWidget(self.Dt, 3, 2)
#         lyt.addWidget(QtWidgets.QLabel(self.tr("Tape thickness")), 4, 1)
#         self.delta = Entrada_con_unidades(Length, "Thickness")
#         lyt.addWidget(self.delta, 4, 2)

#         lytH = QtWidgets.QHBoxLayout()
#         lytH.addWidget(QtWidgets.QLabel(self.tr("Friction factor calculation method")))
#         self.delta = Entrada_con_unidades(Length, "Thickness")
#         self.methodFriction = QtWidgets.QComboBox()
#         for method in TwistedTape.TEXT_FRICTION:
#             self.methodFriction.addItem(method)
#         lytH.addWidget(self.methodFriction)
#         lyt.addLayout(lytH, 5, 1, 1, 2)

#         lytH = QtWidgets.QHBoxLayout()
#         lytH.addWidget(QtWidgets.QLabel(self.tr("Heat transfer calculation method")))
#         self.methodHeat = QtWidgets.QComboBox()
#         for method in TwistedTape.TEXT_HEAT:
#             self.methodHeat.addItem(method)
#         lytH.addWidget(self.methodHeat)
#         lyt.addLayout(lytH, 6, 1, 1, 2)

#         lyt.addItem(QtWidgets.QSpacerItem(
#             10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
#             QtWidgets.QSizePolicy.Policy.Expanding), 10, 3)

#         # self.fill()

#     def setEnabled(self, boolean):
#         """Toggled enable/disable state for all children widget except
#         checkbox used to change this"""
#         for wdg in self.children():
#             if wdg is not self.check:
#                 wdg.setEnabled(boolean)

#     # def fill(self):
#     #     self.check.setChecked(True)
#         # self.check.setChecked(False)


# class Dialog(QtWidgets.QDialog):
#     """Component list config dialog"""
#     def __init__(self, parent=None):
#         super().__init__(parent)
#         self.setWindowTitle(self.tr("Twisted-tape insert"))
#         layout = QtWidgets.QVBoxLayout(self)
#         self.datos = UI_TwistedTape()
#         layout.addWidget(self.datos)
#         self.buttonBox = QtWidgets.QDialogButtonBox(
#             QtWidgets.QDialogButtonBox.StandardButton.Cancel
#             | QtWidgets.QDialogButtonBox.StandardButton.Ok)
#         self.buttonBox.accepted.connect(self.accept)
#         self.buttonBox.rejected.connect(self.reject)
#         layout.addWidget(self.buttonBox)

#     # def value(self, config):
#     #     """Function to result wizard"""
#     #     config = self.datos.value(config)
#     #     return config


if __name__ == "__main__":
    import sys
    # app = QtWidgets.QApplication(sys.argv)
    # Dialog = Dialog()
    # Dialog.show()
    # sys.exit(app.exec())
    print(f_twistedAnnulli_Gupte(4e4, 1, 0.61, pi*2))

