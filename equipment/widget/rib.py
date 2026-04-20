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
from math import atan, log, pi

from scipy.optimize import newton
from tools.qt import QtCore, QtWidgets, translate

from equipment.widget.gui import ToolGui, CallableEntity
from lib.unidades import Length
from lib.utilities import refDoc
from UI.widgets import Entrada_con_unidades


__doi__ = {
    1:
        {"autor": "Naphon, P., Nuchjapo, M., Kurujareon, J.",
         "title": "Tube side heat transfer coefficient and friction factor "
                  "characteristics of horizontal tubes with helical rib",
         "ref": "Energy Conv. Management 47(18-19) (2006) 3031-3044",
         "doi": "10.1016/j.enconman.2006.03.023"},
    2:
        {"autor": "Vicente, P.G., García, A., Viedma, A.",
         "title": "Experimental investigation on heat transfer and frictional "
                  "characteristics of spirally corrugated tubes in turbulent "
                  "flow at different Prandtl numbers",
         "ref": "Int. J. Heat Mass Transfer 47(4) (2004) 671-681",
         "doi": "10.1016/j.ijheatmasstransfer.2003.08.005"},
    3:
        {"autor": "Vicente, P.G., García, A., Viedma, A.",
         "title": "Mixed convection heat transfer and isothermal pressure "
                  "drop in corruageted tubes for laminar and transition flow",
         "ref": "Int. Comm. Heat Mass Transfer 31(5) (2004) 651-662",
         "doi": "10.1016/S0735-1933(04)00052-1"},
    4:
        {"autor": "Sethumadhavan, R., Raja Rao, M.",
         "title": "Turbulent Flow Friction and Heat Transfer Characteristics "
                  "of Single and Multistart Spirally Enghanced Tubes",
         "ref": "J. Heat Transfer 108(1) (1986) 55-61",
         "doi": "10.1115/1.3246905"},
    5:
        {"autor": "Dong, Y., Huixiong, L., Tingkuan, C.",
         "title": "Pressure drop, heat transfer and performance of "
                  "single-phase turbulent flow in spirally corrugated tubes",
         "ref": "Exp. Thermal Fluid Sci. 24(3-4) (2001) 131-138",
         "doi": "10.1016/S0894-1777(01)00047-4"},
    6:
        {"autor": "Srinivasan, V., Christensen, R.N.",
         "title": "Experimental Investigation of Heat Transfer and Pressure "
                  "Drop Characteristics of Flow Through Spirally Fluted Tubes",
         "ref": "Exp. Thermal Fluid Sci. 5(6) (1992) 820-827",
         "doi": "10.1016/0894-1777(92)90126-P"},
    # 7:
    #     {"autor": "",
    #      "title": "",
    #      "ref": "",
    #      "doi": ""},
        }


# Friction factor correlations
@refDoc(__doi__, [2, 3])
def f_corrugated_Vicente(Re, Di, p, h):
    """Calculate friction factor for a corrugated tube using the Vicente et al.
    correlation (2004).

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    Di : float
        Internal diameter of tube, [m]
    p : float
        Helical pitch for twist of 2π radians (360º), [m]
    h : float
        Roughness height, [m]

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    # Tube severity factor
    phi = h**2/p/Di

    if Re > 2000:
        # Turbulent flow,  Eq 2
        f = 1.53 * phi**0.46 / Re**0.16
    else:
        # Laminar flow, Eq 3 in [3]_
        f = 29.8 * phi**0.11 / Re**0.97

    return f


@refDoc(__doi__, [4])
def f_corrugated_Sethumadhavan(Re, Di, P, h):
    """Calculate friction factor for a corrugated pipe using the
    Sethumadhavan-Raja Rao correlation (1986).

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    D : float
        Internal diameter of tube, [m]
    P : float
        helical pitch for twist of 2π radians (360º), [m]
    h : float
        Roughness height, [m]

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    if Re < 5000:
        raise NotImplementedError("Input out of bound")

    Deq = Di-h

    # Eq 11
    def f_res(f):
        """Iterative solution of intrinsic equation"""

        R = (2/f)**0.5 + 2.5*log(2*h/Deq) + 3.75
        h_ = h/Deq * Re * (f/2)**0.5
        return R*h**2/(P*Deq)**0.33 - 0.40*h_**0.164

    fo = f_corrugated_Vicente(Re, Di, P, h)
    f = newton(f_res, fo)

    if isinstance(f, complex):
        raise ValueError("Solution don't converge")

    return f


@refDoc(__doi__, [5])
def f_corrugated_Dong(Re, Di, P, h):
    """Calculate friction factor for a corrugated pipe using the Dong et al.
    correlation (2001).

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    D : float
        Internal diameter of tube, [m]
    P : float
        helical pitch for twist of 2π radians (360º), [m]
    h : float
        Roughness height, [m]

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    if Re < 6000:
        raise NotImplementedError("Input out of bound")

    # Helix angle
    alpha = atan(P/pi/Di)*180/pi

    # Tube severity factor
    phi = h**2/P/Di

    # Eq 13
    def f_res(f):
        """Iterative solution of intrinsic equation"""

        R = (2/f)**0.5 + 2.5*log(2*h/Di) + 3.75
        h_ = h/Di * Re * (f/2)**0.5
        return R * phi**0.317 * (alpha/50)**0.16 - 0.455*h_**0.169

    fo = f_corrugated_Vicente(Re, Di, P, h)
    f = newton(f_res, fo)

    if isinstance(f, complex):
        raise ValueError("Solution don't converge")

    return f


# Heat Transfer coefficient correlations
@refDoc(__doi__, [2, 3])
def Nu_corrugated_Vicente(Re, Pr, Di, p, h):
    """Calculate friction factor for a corrugated tube using the Vicente et al.
    correlation (2004).

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    Pr : float
        Prandtl number, [-]
    Di : float
        Internal diameter of tube, [m]
    p : float
        Helical pitch for twist of 2π radians (360º), [m]
    h : float
        Roughness height, [m]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    # Tube severity factor
    phi = h**2/p/Di

    if Re > 2000:
        # Turbulent flow, Eq 9
        Nu = 0.374 * phi**0.25 * (Re-1500)**0.74 * Pr**0.44
    else:
        # Laminar flow, for simplicity use the constant infinite value to avoid
        # Rayleigh input parameters
        Nu = 4.36

    return Nu


@refDoc(__doi__, [4])
def Nu_corrugated_Sethumadhavan(Re, Pr, Di, P, h):
    """Calculate Nusselt number for a corrugated pipe using the
    Sethumadhavan-Raja Rao correlation (1986).

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    Pr : float
        Prandtl number, [-]
    Di : float
        Internal diameter of tube, [m]
    p : float
        Helical pitch for twist of 2π radians (360º), [m]
    h : float
        Roughness height, [m]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    if Re < 5000:
        raise NotImplementedError("Input out of bound")

    Deq = Di-h

    f = f_corrugated_Vicente(Re, Di, P, h)
    h_ = h/Deq * Re * (f/2)**0.5
    R = (2/f)**0.5 + 2.5*log(2*h/Deq) + 3.75

    # Eq 15
    G = 8.6 * h_**0.13 * Pr**0.55

    # Eq 7
    St = 1/((((G-R)*(f/2)**0.5)+1)*2/f)

    return St*Re*Pr


@refDoc(__doi__, [5])
def Nu_corrugated_Dong(Re, Pr, Di, P, h):
    """Calculate nusselt number for a corrugated pipe using the Dong et al.
    correlation (2001).

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    Pr : float
        Prandtl number, [-]
    D : float
        Internal diameter of tube, [m]
    P : float
        helical pitch for twist of 2π radians (360º), [m]
    h : float
        Roughness height, [m]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    if Re < 6000:
        raise NotImplementedError("Input out of bound")

    # Helix angle
    alpha = atan(P/pi/Di)*180/pi
    Deq = Di-h

    f = f_corrugated_Vicente(Re, Di, P, h)
    h_ = h/Deq * Re * (f/2)**0.5
    R = (2/f)**0.5 + 2.5*log(2*h/Deq) + 3.75

    # Eq 16
    G = 7.33*h_**0.175 / (alpha/50)**0.16 * Pr**0.548

    # Eq 7
    St = 1/((((G-R)*(f/2)**0.5)+1)*2/f)

    return St*Re*Pr


# Fluted tubes
@refDoc(__doi__, [1])
def f_fluted_Srinivasan(Re, Di, p, h, N):
    """Calculate friction factor for a spirally flutud tube using the
    Srinivasan correlation (1992).

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    Di : float
        Internal diameter of tube, [m]
    p : float
        Helical rib pitch for twist of 2π radians (360º), [m]
    h : float
        Flute depth, [m]
    N : integer
        Number of flute starts, [-]

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    # Helix angle
    alpha = atan(p/pi/Di/N)*180/pi
    tita = alpha/90

    # Eq 2
    f = 12.745 * Re**(-0.474-0.209*p/Di+0.685*tita) * \
        (h/Di)**(1.292+0.031*p/Di) * p/Di**(9.908+0.331e-5*Re-12.074*tita)

    return f


@refDoc(__doi__, [6])
def Nu_fluted_Srinivasan(Re, Pr, Di, p, h, N):
    """Calculate Nusselt number for a tube with helical rib using the Naphon
    correlation (1992).

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    Pr : float
        Prandtl number, [-]
    Di : float
        Internal diameter of tube, [m]
    p : float
        Helical rib pitch for twist of 2π radians (360º), [m]
    h : float
        Helical rib depth, [m]
    N : integer
        Number of flute starts, [-]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    # Helix angle
    alpha = atan(p/pi/Di/N)*180/pi
    tita = alpha/90

    if Re < 5000:
        # Eq 4, laminar flow
        Nu = 0.014*Re**0.847*Pr**0.4/(h/Di)**0.067/(p/Di)**0.293/tita**0.705
    else:
        # Eq 5, turbulent flow
        Nu = 0.064*Re**0.773*Pr**0.4/(h/Di)**0.242/(p/Di)**0.108/tita**0.599

    return Nu


# Rib correlation
@refDoc(__doi__, [1])
def f_rib_Naphon(Re, Di, p, h):
    """Calculate friction factor for a tube with helical rib using the Naphon
    correlation (2006).

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    Di : float
        Internal diameter of tube, [m]
    p : float
        Helical rib pitch for twist of 2π radians (360º), [m]
    h : float
        Helical rib depth, [m]

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    # Eq 11
    f = 7.85 / Re**0.21 * (h/Di)**1.68 / (p/Di)**0.54

    return f


@refDoc(__doi__, [1])
def Nu_rib_Naphon(Re, Pr, Di, p, h):
    """Calculate Nusselt number for a tube with helical rib using the Naphon
    correlation (2006).

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    Pr : float
        Prandtl number, [-]
    Di : float
        Internal diameter of tube, [m]
    p : float
        Helical rib pitch for twist of 2π radians (360º), [m]
    h : float
        Helical rib depth, [m]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    # Eq 8
    Nu = 44.26 * (Re-1500)**0.27 / Pr**0.26 * (h/Di)**0.89 / (p/Di)**0.96

    return Nu


class Rib(CallableEntity):
    """Helical rip tube

    Parameters
    ----------
    p : float
        Helical rib pitch for twist of 2π radians (360º), [m]
    h : float
        Helical rib depth, [m]
    """

    TEXT_FRICTION = (
        "Vicente (2004)",
        "Sethumadhavan-Raja Rao (1986)",
        "Dong (2001)")
    TEXT_HEAT = (
        "Vicente (2004)",
        "Sethumadhavan-Raja Rao (1986)",
        "Dong (2001)")

    status = 0
    msg = ""
    kw = {
        "p": 0,
        "h": 0,

        "methodCorrugatedFriction": 0,
        "methodCorrugatedHeat": 0,

        "isRib": False,
        "isFluted": False,
        "N": 1}

    valueChanged = QtCore.pyqtSignal(object)
    inputChanged = QtCore.pyqtSignal(object)

    @property
    def isCalculable(self):
        """Check if all input are defined"""
        if not self.kw["p"]:
            self.msg = translate("equipment", "undefined rib pitch")
            self.status = 0
            return False
        if not self.kw["h"]:
            self.msg = translate("equipment", "undefined rib depth")
            self.status = 0
            return False

        self.msg = ""
        self.status = 1
        return True

    def calculo(self):
        """Definition of twisted tape inserts for annuli sections"""
        self.p = self.kw["p"]
        self.h = self.kw["h"]

        self.valueChanged.emit(self)

    def Nu(self, Re, Pr, Di):
        """Calculate nusselt number"""
        msg = ""
        if self.kw["isRib"]:
            Nu = Nu_rib_Naphon(Re, Pr, Di, self.p, self.h)
            return Nu
        if self.kw["isFluted"]:
            Nu = Nu_fluted_Srinivasan(Re, Pr, Di, self.p, self.h, self.kw["N"])
            return Nu

        if self.kw["methodCorrugatedHeat"] == 1:
            # "Sethumadhavan-Raja Rao (1986)",
            try:
                Nu = Nu_corrugated_Sethumadhavan(Re, Pr, Di, self.p, self.h)
            except RuntimeError:
                Nu = Nu_corrugated_Vicente(Re, Pr, Di, self.p, self.h)
                msg = "Sethumadhavan method don't converge, using Vicente method instead"

        elif self.kw["methodCorrugatedHeat"] == 2:
            # "Dong (2001)")
            Nu = Nu_corrugated_Dong(Re, Pr, Di, self.p, self.h)

        else:
            # "Vicente (2004)",
            Nu = Nu_corrugated_Vicente(Re, Pr, Di, self.p, self.h)

        if msg:
            self.status = 3
            self.msg = translate("equipment", msg)
            self.inputChanged.emit(self)

        return Nu

    def f(self, Re, Di):
        """Calculate friction factor"""
        msg = ""

        if self.kw["isRib"]:
            f = f_rib_Naphon(Re, Di, self.p, self.h)
            return f
        if self.kw["isFluted"]:
            f = f_fluted_Srinivasan(Re, Di, self.p, self.h, self.kw["N"])
            return f

        if self.kw["methodCorrugatedFriction"] == 1:
            # "Sethumadhavan-Raja Rao (1986)",
            try:
                f = f_corrugated_Sethumadhavan(Re, Di, self.p, self.h)
            except RuntimeError:
                f = f_corrugated_Vicente(Re, Di, self.p, self.h)
                msg = "Sethumadhavan method don't converge, using Vicente method instead"

        elif self.kw["methodCorrugatedFriction"] == 2:
            # "Dong (2001)")
            f = f_corrugated_Dong(Re, Di, self.p, self.h)

        else:
            # "Vicente (2004)",
            f = f_corrugated_Vicente(Re, Di, self.p, self.h)

        if msg:
            self.status = 3
            self.msg = translate("equipment", msg)
            self.inputChanged.emit(self)

        return f


class UI_Rib(ToolGui):
    """Helical rib dialog"""

    title = translate("equipment", "Use helical rib in tube")

    def loadUI(self):
        """Add widget"""
        self.Entity = Rib()

        lyt = self.wdg.layout()

        self.corrugated = QtWidgets.QRadioButton(self.tr("Corrugated tube"))
        self.corrugated.setChecked(True)
        lyt.addWidget(self.corrugated, 1, 1, 1, 3)

        groupMethods = QtWidgets.QWidget()
        lytM = QtWidgets.QGridLayout(groupMethods)
        lytM.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 1, 0, 2, 1)
        lytM.addWidget(QtWidgets.QLabel(
            self.tr("Friction factor method")), 1, 1)
        self.methodCorrugatedFriction = QtWidgets.QComboBox()
        for method in Rib.TEXT_FRICTION:
            self.methodCorrugatedFriction.addItem(method)
        self.methodCorrugatedFriction.currentIndexChanged.connect(
            partial(self.changeParams, "methodCorrugatedFriction"))
        lytM.addWidget(self.methodCorrugatedFriction, 1, 2)
        lytM.addWidget(QtWidgets.QLabel(
            self.tr("Heat transfer method")), 2, 1)
        self.methodCorrugatedHeat = QtWidgets.QComboBox()
        for method in Rib.TEXT_HEAT:
            self.methodCorrugatedHeat.addItem(method)
        self.methodCorrugatedHeat.currentIndexChanged.connect(
            partial(self.changeParams, "methodCorrugatedHeat"))
        lytM.addWidget(self.methodCorrugatedHeat, 2, 2)
        lyt.addWidget(groupMethods, 2, 1, 1, 2)
        self.corrugated.toggled.connect(groupMethods.setEnabled)

        self.rib = QtWidgets.QRadioButton(self.tr("Rib tube"))
        self.rib.toggled.connect(partial(self.changeParams, "isRib"))
        lyt.addWidget(self.rib, 4, 1, 1, 3)

        self.fluted = QtWidgets.QRadioButton(self.tr("Fluted tube"))
        self.fluted.toggled.connect(self.setEnableFluted)
        lyt.addWidget(self.fluted, 5, 1, 1, 3)

        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 6, 1)
        label = QtWidgets.QLabel(self.tr("Rib pitch"))
        label.setToolTip(self.tr("Pitch for twist of 2π radians (360º)"))
        lyt.addWidget(label, 8, 1)
        self.p = Entrada_con_unidades(Length)
        self.p.valueChanged.connect(partial(self.changeParams, "p"))
        lyt.addWidget(self.p, 8, 2)
        lyt.addWidget(QtWidgets.QLabel("Rib depth"), 9, 1)
        self.h = Entrada_con_unidades(Length, "Thickness")
        self.h.valueChanged.connect(partial(self.changeParams, "h"))
        lyt.addWidget(self.h, 9, 2)
        self.lblN = QtWidgets.QLabel(self.tr("Rib count"))
        lyt.addWidget(self.lblN, 10, 1)
        self.N = Entrada_con_unidades(int, spinbox=True)
        self.N.valueChanged.connect(partial(self.changeParams, "N"))
        lyt.addWidget(self.N, 10, 2)

        self.Entity.valueChanged.connect(self.valueChanged.emit)
        self.Entity.inputChanged.connect(self.populate)
        self.lblN.setVisible(False)
        self.N.setVisible(False)

    def setEnableFluted(self, boolean):
        self.lblN.setVisible(boolean)
        self.N.setVisible(boolean)
        self.changeParams("isFluted", boolean)

class Dialog(QtWidgets.QDialog):
    """Component list config dialog"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(self.tr("Twisted-tape insert"))
        layout = QtWidgets.QVBoxLayout(self)
        self.datos = UI_Rib()
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
