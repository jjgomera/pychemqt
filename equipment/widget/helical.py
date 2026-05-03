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
from lib.adimensional import Dean
from lib.friction import f_friccion
from lib.unidades import Length
from lib.utilities import refDoc
from UI.widgets import Entrada_con_unidades


__doi__ = {
    1:
        {"autor": "El-Genk, M.S., Timothy, M.S.",
         "title": "A Review and Correlations for Convection Heat Transfer and "
                  "Pressure Losses in Toroidal and Helically Coiled Tubes",
         "ref": "Heat Transfer Eng. 38(5) (2017) 447-474",
         "doi": "10.1080/01457632.2016.1194693"},
    2:
        {"autor": "",
         "title": "Perry's Chemical Engineers' Handbook 9th Edition",
         "ref": "McGraw-Hill (2019)",
         "doi": ""},
    3:
        {"autor": "Schmidt, E.F.",
         "title": "Wärmeübergand und Druckverlust in Rohrschlangen",
         "ref": "Chemie Ingenieur Technik 39(13) (1967) 781-789",
         "doi": "10.1002/cite.330391302"},
    4:
        {"autor": "Ito, H.",
         "title": "Friction Factors for Turbulent Flow in Curved Pipes",
         "ref": "J. Basic Eng. 81 (1959) 123-134",
         "doi": "10.1115/1.4008390"},
    5:
        {"autor": "Kubair, V., Kuloor, N.R.",
         "title": "Heat Transfer to Newtonian Fluids in Coiled Pipes in "
                  "Laminar Flow",
         "ref": "Int. J. Heat Mass Transfer 9 (1966) 63-75",
         "doi": "10.1016/0017-9310(66)90057-3"},
    6:
        {"autor": "Srinivasan, P.S., Nandapurkar, S.S., Holland, F.A.",
         "title": "Pressure Drop and Heat Transfer in Coils",
         "ref": "Chemical Engineer, vol. 218, CE131–119, 1968.",
         "doi": ""}
    7:
        {"autor": "White, C.M.",
         "title": "Streamline Flow through Curved Pipes",
         "ref": "Proc. R .Soc. London A 123 (1929) 645-63",
         "doi": "10.1098/rspa.1929.0089"}

    # 5:
        # {"autor": "",
         # "title": "",
         # "ref": "",
         # "doi": ""}

    20:
        {"autor": "Xin, R.C., Ebadian, M.A.  ",
         "title": "The Effects of Prandtl Numbers on Local and Average "
                  "Convective Heat Transfer Characteristics in Helical Pipes",
         "ref": "J. Heat Transfer 119(3) (1997) 467-73",
         "doi": "10.1115/1.2824120."}
}


# Critical Reynolds number correlations
@refDoc(__doi__, [2])
def helical_transition_Re_Srinivasan(Di, Dc):
    r'''Calculates the transition Reynolds number for flow inside a curved or
    helical coil between laminar and turbulent flow, using the method of [1]_,
    also shown in [2]_ and [3]_. Correlation recommended in [3]_.

    .. math::
        Re_{crit} = 2100\left[1 + 12\left(\frac{D_i}{D_c}\right)^{0.5}\right]

    Parameters
    ----------
    Di : float
        Inner diameter of the coil, [m]
    Dc : float
        Diameter of the helix/coil measured from the center of the tube on one
        side to the center of the tube on the other side, [m]

    Returns
    -------
    Re_crit : float
        Transition Reynolds number between laminar and turbulent [-]

    Notes
    -----
    At very low curvatures, converges to Re = 2100.
    Recommended for :math:`0.004 < d_i/D_c < 0.1`.

    Examples
    --------
    >>> helical_transition_Re_Srinivasan(1, 7.)
    11624.704719832524

    References
    ----------
    .. [1] Srinivasan, P. S., Nandapurkar, S. S., and Holland, F. A., "Pressure
       Drop and Heat Transfer in Coils", Chemical Engineering, 218, CE131-119,
       (1968).
    .. [2] El-Genk, Mohamed S., and Timothy M. Schriener. "A Review and
       Correlations for Convection Heat Transfer and Pressure Losses in
       Toroidal and Helically Coiled Tubes." Heat Transfer Engineering 0, no. 0
       (June 7, 2016): 1-28. doi:10.1080/01457632.2016.1194693.
    .. [3] Rohsenow, Warren and James Hartnett and Young Cho. Handbook of Heat
       Transfer, 3E. New York: McGraw-Hill, 1998.
    '''
    return 2100.*(1. + 12.*sqrt(Di/Dc))



@refDoc(__doi__, [3])
def Rec_Schmidt(di, Dc):
    r"""Calculates critical Reynolds to define transition between laminar and
    turbulent flow using using the correlation of Schmidt (1967)

    .. math::
        Re_c = 2300 \left(1+8.6\left(\frac{di}{Dc}\right)^{0.45}\right)

    Parameters
    ----------
    di : float
        Inner diameter of the pipe, [m]
    Dc : float
        Diameter of the helix, [m]

    Returns
    -------
    Rec : float
        Critical reynolds number, [-]
    """
    # Eq 14
    Rec = 2300*(1+8.6*(di/Dc)**0.45)
    return Rec


@refDoc(__doi__, [4])
def Rec_Ito(di, Dc):
    r"""Calculates critical Reynolds to define transition between laminar and
    turbulent flow using using the correlation of Ito (1959)

    .. math::
        Re_c = 2x10^4 \left(\frac{di}{Dc}\right)^{0.32}

    Parameters
    ----------
    di : float
        Inner diameter of the pipe, [m]
    Dc : float
        Diameter of the helix, [m]

    Returns
    -------
    Rec : float
        Critical reynolds number, [-]
    """
    # Eq 11
    Rec = 2e4*(di/Dc)**0.32
    return Rec


@refDoc(__doi__, [5])
def Rec_Kubair(di, Dc):
    r"""Calculates critical Reynolds to define transition between laminar and
    turbulent flow using using the correlation of Kubair-Kuloor (1966)

    .. math::
        Re_c = 1.273x10^4 \left(\frac{di}{Dc}\right)^{0.2}

    Parameters
    ----------
    di : float
        Inner diameter of the pipe, [m]
    Dc : float
        Diameter of the helix, [m]

    Returns
    -------
    Rec : float
        Critical reynolds number, [-]
    """
    Rec = 1.273e4*(di/Dc)**0.2
    return Rec


@refDoc(__doi__, [5, 2])
def Rec_Srinivasan(di, Dc):
    r"""Calculates critical Reynolds to define transition between laminar and
    turbulent flow using using the correlation of Srinivasan (1968). Recomended
    method by [2]_.

    .. math::
        Re_c = 2100 \left(1 + 12\sqrt{\frac{d_i}{D_c}}\right)

    Parameters
    ----------
    di : float
        Inner diameter of the pipe, [m]
    Dc : float
        Diameter of the helix, [m]

    Returns
    -------
    Rec : float
        Critical reynolds number, [-]
    """
    Rec = 2100 * (1 + 12*(di/Dc)**0.5)
    return Rec


# Friction factor correlations
@refDoc(__doi__, [3])
def f_Schmidt(Re, di, Dc):
    """Calculate friction factor for internal flow of a helical coil using
    the correlation of Schmidt (1967)

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    di : float
        Inner diameter of the pipe, [m]
    Dc : float
        Diameter of the helix, [m]

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    Rec = Rec_Schmidt(di, Dc)

    if Re < Rec:
        # Laminar flow, Eq 15
        f = 64/Re * (1+0.14*(di/Dc)**0.97*Re**(1-0.644*(di/Dc)**0.312))
    elif Re < 2.2e4:
        # Eq 16
        f = 0.3164/Re**0.25 * (1+2.88e4/Re*(di/Dc)**0.62)
    else:
        # Eq 17
        f = 0.3164/Re**0.25 * (1+0.0823*(1+di/Dc)*(di/Dc)**0.53*Re**0.25)

    return f


@refDoc(__doi__, [7])
def f_laminar_White(Re, di, Dc):
    r"""Calculates friction factor for internal flow of a helical coil in
    laminar flow using the method of White (1929).

    .. math::
        f_c = \frac{f_{s,L}} {1 - \left(1-\left(\frac{11.6}{De}\right)^{0.45}
        \right)^{\frac{1}{0.45}}

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    di : float
        Inner diameter of the pipe, [m]
    Dc : float
        Diameter of the helix, [m]

    Returns
    -------
    f : float
        Friction factor, [-]

    """
    De = Dean(Re=Re, Di=Di, D=Dc)
    fd = f_friccion(Re)
    if De > 11.6:
        C = 1 - (1 - (11.6/De)**0.45)**(1/0.45)
    else:
        C = 1
    return fd/C


# Heat Transfer coefficient correlations
@refDoc(__doi__, [3])
def Nu_Schmidt(Re, Pr, di, Dc):
    r"""Calculates Nusselt number for internal flow of a helical coil using the
    correlation of Schmidt (1967)

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    Pr : float
        Prandtl number, [-]
    di : float
        Inner diameter of the pipe, [m]
    Dc : float
        Diameter of the helix, [m]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    Rec = Rec_Schmidt(di, Dc)

    if Re < Rec:
        # Laminar flow, Eq 18
        Nu = 3.65 + Pr**0.8 * 0.08*(1+0.8*(di/Dc)**0.9) \
            * Re**(0.5+0.2903*(di/Dc)**0.194)
    elif Re < 2.2e4:
        # Eq 21
        Nu = 0.023 * Pr**(1/3) * (1+14.8*(1+di/Dc)*(di/Dc)**(1/3)) \
            * Re**(0.8-0.22*(di/Dc)**0.1)
    else:
        # Eq 22
        Nu = 0.023 * (1+3.6*(1-di/Dc)*(di/Dc)**0.8) * Re**0.8 * Pr**(1/3)

    return Nu


@refDoc(__doi__, [20, 1])
def Nu_XinEbadian(Re, Pr, di, Dc):
    r"""Calculates Nusselt number for internal flow of a helical coil using the
    correlation of Xin-Ebadian (1997)

    For laminar flow:

    .. math::
        Nu = \left(2.153 + 0.318 \left(Re \frac{d_i}{D_c}\right)^{0.643}\right)
        Pr^{0.177}

    For turbulent flow:

    .. math::
        Nu = 0.00619 Re^{0.92} Pr^{0.4} \left(1 + 3.455 \frac{d_i}{D_c}\right)

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    Pr : float
        Prandtl number, [-]
    di : float
        Inner diameter of the pipe, [m]
    Dc : float
        Diameter of the helix, [m]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    # Laminar flow
    if Re < 5000:
        De = Re*(di/Dc)**0.5

        # Eq 5
        Nu = (2.153+0.318*De**0.643) * Pr**0.177
    else:
        # Eq 6
        Nu = 0.00619 * Re**0.92 * Pr**0.4 * (1+3.455*di/Dc)

    return Nu


class Helical(CallableEntity):
    """Helical coil tube used as anhancing heat transfer equipment.

    Parameters
    ----------
    H : float
        Tape pitch for twist of π radians (180º), [m]
    di : float
        Inner diameter of the pipe, [m]
    Dc : float
        Diameter of the helix, [m]
    """

    TEXT_REYNOLDS_CRITICAL = (
        "Schmidt (1967)"
        "Ito (1959)",
        "Kubair-Kuloor (1966)",
        "Srinivasan (1968)")

    status = 0
    msg = ""
    kw = {
        "methodReCritic": 0,
        "methodFrictionLaminar": 0,
        "methodFrictionTurbulent": 0,
        "methodHeatLaminar": 0,
        "methodHeatTurbulent": 0,

        "H": 0,
        "di": 0,
        "Dc": 0}

    valueChanged = QtCore.pyqtSignal(object)
    inputChanged = QtCore.pyqtSignal(object)

    @property
    def isCalculable(self):
        """Check if all input are defined"""
        if not self.kw["di"]:
            self.msg = translate("equipment", "undefined internal diameter")
            self.status = 0
            return False
        if not self.kw["Dc"]:
            self.msg = translate("equipment", "undefined diameter of helix")
            self.status = 0
            return False

        self.msg = ""
        self.status = 1
        return True

    def calculo(self):
        """Definition of twisted tape inserts for annuli sections"""
        self.di = self.kw["di"]
        self.Dc = self.kw["Dc"]

        self.valueChanged.emit(self)

    @property
    def ReCritical(self):
        """Calculate critical Reynolds number to define transition of regimen
        flow from laminar to turbulent"""
        if self.kw["methodReCritic"] == 1:
            # Ito (1959)
            Rec = Rec_Ito(self.di, self.Dc)

        elif self.kw["methodReCritic"] == 2:
            # Kubair-Kuloor (1966)
            Rec = Rec_Kubair(self.di, self.Dc)

        elif self.kw["methodReCritic"] == 2:
            # Srinivasan (1968)
            Rec = Rec_Srinivasan(self.di, self.Dc)

        else:
            # Schmidt (1967)
            Rec = Rec_Schmidt(self.di, self.Dc)

        return Rec

    def Nu(self, Re, Pr):
        """Calculate nusselt number"""
        Rec = self.ReCritical

        if Re < Rec:
            # Laminar flow
            pass
        else:
            Nu = Nu_XinEbadian(Re, Pr, self.di, self.Dc)

        return Nu

    def f(self, Re):
        """Calculate friction factor"""
        f = f_Schmidt(Re, self.di, self.Dc)

        return f


class UI_Helical(ToolGui):
    """Helical coil dialog"""

    title = translate("equipment", "Use helical coil")

    def loadUI(self):
        """Add widget"""
        self.Entity = Helical()

        lyt = self.wdg.layout()

        # label = QtWidgets.QLabel(self.tr("Tape pitch, H"))
        # label.setToolTip(self.tr("Tape pitch for twist of π radians (180º)"))
        # lyt.addWidget(label, 2, 1)
        # self.H = Entrada_con_unidades(Length)
        # self.H.valueChanged.connect(partial(self.changeParams, "H"))
        # lyt.addWidget(self.H, 2, 2)
        # label = QtWidgets.QLabel("di")
        # label.setToolTip(self.tr("Internal diameter of annuli section"))
        # lyt.addWidget(label, 3, 1)
        # self.di = Entrada_con_unidades(Length, "PipeDiameter")
        # self.di.valueChanged.connect(partial(self.changeParams, "di"))
        # lyt.addWidget(self.di, 3, 2)
        # label = QtWidgets.QLabel("Do")
        # label.setToolTip(self.tr("External diameter of annuli section"))
        # lyt.addWidget(label, 4, 1)
        # self.Do = Entrada_con_unidades(Length, "PipeDiameter")
        # self.Do.valueChanged.connect(partial(self.changeParams, "Do"))
        # lyt.addWidget(self.Do, 4, 2)

        # self.angled = QtWidgets.QCheckBox(self.tr("Angled twisted-tape"))
        # self.angled.toggled.connect(self.setEnableOrientation)
        # lyt.addWidget(self.angled, 5, 1, 1, 2)

        # lytH = QtWidgets.QHBoxLayout()
        # self.labelOrientation = QtWidgets.QLabel(self.tr(
            # "Direction of flow relative to the tape curvature"))
        # lytH.addWidget(self.labelOrientation)
        # self.orientation = QtWidgets.QComboBox()
        # for method in TwistedTapeAnnuli.TEXT_ORIENTACION:
            # self.orientation.addItem(method)
        # self.orientation.currentIndexChanged.connect(
            # partial(self.changeParams, "orientation"))
        # lytH.addWidget(self.orientation)
        # lyt.addLayout(lytH, 7, 1, 1, 2)
        self.Entity.valueChanged.connect(self.valueChanged.emit)
        self.Entity.inputChanged.connect(self.populate)

    def setEnableOrientation(self, boolean):
        """Change Enable/Disable state for orientation of twisted tape"""
        self.labelOrientation.setEnabled(boolean)
        self.orientation.setEnabled(boolean)
        self.changeParams("angled", boolean)


class Dialog(QtWidgets.QDialog):
    """Component list config dialog"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(self.tr("Twisted-tape insert"))
        layout = QtWidgets.QVBoxLayout(self)
        self.datos = UI_Helical()
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
