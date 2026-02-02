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


from math import pi

from tools.qt import QtWidgets

from lib.unidades import Dimensionless, Area, Length
from lib.utilities import refDoc
from UI.widgets import Entrada_con_unidades


__doi__ = {
    1:
        {"autor": "du Plessis, J.P., Kröger, D.G.",
         "title": "Friction factor prediction for fully developed laminar "
                  "twisted-tape flow",
         "ref": "Int. J. Heat Mass Transfer 27(11) (1984) 2095-2100",
         "doi": "10.1016/0017-9310(84)90196-0"},
    2:
        {"autor": "du Plessis, J.P., Kröger, D.G.",
         "title": "Heat transfer correlation for thermally developing laminar "
                  "flow in a smooth tube with a twisted-tape insert",
         "ref": "Int. J. Heat Mass Transfer 30(3) (1987) 509-515",
         "doi": "10.1016/0017-9310(87)90265-1"},
    3:
        {"autor": "Shah, R.K., London, A.L.",
         "title": "Laminar Flow Forced Convection in Ducts: A Source Book for "
                  "Compact Heat Exchanger Analytical Data",
         "ref": "Academic Press 1978",
         "doi": ""},
    4:
        {"autor": "Manglik, R.M., Bergles, A.E.",
         "title": "Heat Transfer and Pressure Drop Correlations for "
                  "Twisted-Tape Inserts in Isothermal Tubes: Part I - Laminar "
                  "Flows",
         "ref": "J. Heat Transfer 115(4) (1993) 881-889",
         "doi": "10.1115/1.2911383"},
    5:
        {"autor": "Manglik, R.M., Bergles, A.E.",
         "title": "Heat Transfer and Pressure Drop Correlations for "
                  "Twisted-Tape Inserts in Isothermal Tubes: Part II - "
                  "Transition and Turbulent Flows",
         "ref": "J. Heat Transfer 115(4) (1993) 890-896",
         "doi": "10.1115/1.2911384"},
    6:
        {"autor": "Hong, S.W., Bergles, A.E.",
         "title": "Augmenttion of Laminar Flow Heat Transfer in Tubes by "
                  "Means of Twisted-Tape Inserts",
         "ref": "J. Heat Transfer 98(2) (1976) 251-256",
         "doi": "10.1115/1.3450527"},
    7:
        {"autor": "Lopina, R.F., Bergles, A.E.",
         "title": "Heat Transfer and Pressure Drop in Tape-Generaged Swirl "
                  "Flow of Single-Phase Water",
         "ref": "ASME J. Heat Transfer 91(3) (1969) 434-442",
         "doi": "10.1115/1.3580212"},
    8:
        {"autor": "",
         "title": "HTRI Design Manual",
         "ref": "",
         "doi": ""},
    9:
        {"autor": "Naphon, P.",
         "title": "Heat transfer and pressure drop in the horizontal double "
                  "pipes with and without twisted tape insert",
         "ref": "Int. Comm. Heat Mass Transfer 33 (2006) 166-175",
         "doi": "10.1016/j.icheatmasstransfer.2005.09.007"},
    10:
        {"autor": "Agarwal, S.K., Raja Rao, M.",
         "title": "Heat transfer augmentation for the flow of a viscous "
                  "liquid in circular tubes using twisted tape inserts",
         "ref": "Int. J. Heat Mass Transfer 39(17) (1996) 3547-3557",
         "doi": ""},
    # 12:
        # {"autor": "",
         # "title": "",
         # "ref": "",
         # "doi": ""},
        }


# Friction factor correlations
@refDoc(__doi__, [1])
def f_twisted_Plessis(Re, D, H, delta, Ae, De):
    """Calculate friction factor for a pipe with a twisted-tape insert using
    the Plessis and Kröger correlation (1984)

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    D : float
        Internal diameter of tube, [m]
    delta : float
        Tape thickness, [m]
    H : float
        Tape pitch for twist of π radians (180º), [m]
    Ae : float
        Effective flow area, [m²]
    De : float
        Effective hydraulic diameter, [m]

    Returns
    -------
    f : float
        Friction factor, [-]

    Examples
    --------
    Selected point from Table 2 in [1]_

    >>> st = TwistedTape(10, 1, 0)
    >>> print("%0.3f" % f_twisted_Plessis(50, 1, 10, 0, st.Ae, st.De))
    0.849

    >>> print("%0.4f" % f_twisted_Plessis(2000, 1, 10, 0, st.Ae, st.De))
    0.0296
    """
    A = pi*D**2/4
    y = H/D

    # Eq 15
    Deltae = A*D**2/Ae/De**2

    # Eq 16
    fe = Deltae/Re*(15.767-0.14706*delta/D)

    # Eq 17
    f = fe * (1 + (Re/(70*y**1.3))**1.5)**(1/3)
    return f


@refDoc(__doi__, [3])
def f_twisted_Shah(Re, D, H, delta):
    """Calculate friction factor for a pipe with a twisted-tape insert using
    the Shah and London correlation (1978)

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    D : float
        Internal diameter of tube, [m]
    H : float
        Tape pitch for twist of π radians (180º), [m]
    delta : float
        Tape thickness, [m]

    Returns
    -------
    f : float
        Friction factor, [-]

    """
    # Chapter XVI: Longitudinal Fins and Twisted Tapes within Ducts
    # F: Circular Duct with a Twisted Tape, Pag 379 and on next

    Xl = H/D

    # Eq 559
    C = 8.8201*Xl - 2.1193*Xl**2 + 0.2108*Xl**3 - 0.0069*Xl**4

    # Eq 563
    Xi = (pi/(pi+2))**2 * ((pi+2-2*delta/D)/(pi-4*delta/D))**2 \
        * (pi/(pi-4*delta/D))

    if Re/Xl < 6.7:
        # Eq 560
        fRe = 42.23*Xi
    elif Re/Xl > 100:
        # Eq 562
        fRe = C*(Re/Xl)**0.3*Xi
    else:
        # Eq 561
        fRe = 38.4*(Re/Xl)**0.05*Xi

    return fRe/Re


@refDoc(__doi__, [4, 5])
def f_twisted_Manglik(Re, D, H, delta, Dh):
    """Calculate friction factor for a pipe with a twisted-tape insert using
    the Manglik and Bergles correlation (1993)

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    D : float
        Internal pipe diameter, [m]
    H : float
        Tape pitch for twist of π radians (180º), [m]
    delta : float
        Tape thickness, [m]
    Dh : float
        Hydraulic diameter, [m]
    """

    if Re <= 2000:
        # Laminar flow

        Resw = Re*(pi/(pi-4*delta/D)) * ((pi*D/H)**2)**0.5
        Sw = Resw/(H/2/D)**0.5

        fsw = 15.767/Resw * ((pi+2-2*delta/D)/(pi-4*delta/D))**2 \
            * (1+1e-6*Sw**2.55)**(1/6)

        f = fsw * Dh/D*(1+(pi*D/H)**2)**1.5

    elif Re > 1e4:
        # Turbulent flow
        f = 0.0791/Re**0.25 * (pi/(pi-4*delta/D))**1.75 \
            * ((pi+2-2*delta/D)/(pi-4*delta/D))**1.25 * (1+2.752/(H/D)**1.29)

    else:
        # Transition flow
        fl = f_twisted_Manglik(2000, D, H, delta, Dh)
        ft = f_twisted_Manglik(1e4, D, H, delta, Dh)

        # Eq 11 in 5
        f = (fl**10 + ft**10)**0.1

    return f


@refDoc(__doi__, [7, 8])
def f_twisted_Lopina(Re, D, H, Dh):
    """Calculate friction factor for a pipe with a twisted-tape insert using
    the Lopina and Bergles correlation (1969). Only valid for turbulent flow
    with Re > 5000

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    D : float
        Internal diameter of tube, [m]
    H : float
        Tape pitch for twist of π radians (180º), [m]
    Dh : float
        Hydraulic diameter, [m]

    Returns
    -------
    f : float
        Friction factor, [-]

    Examples
    --------
    B2.2.1.2.3 in [8]_
    >>> print("%0.2f" % f_twisted_Lopina(24491, 0.61, 6.1, 1.07))
    0.01
    """
    # Using the modified parameter in HTRI Design Manual, section B2.1.2.1

    Reh = Re*Dh/D
    y = H/D

    f = 3.8/y**0.406*(0.046/Reh**0.2)
    return f


@refDoc(__doi__, [9])
def f_twisted_Naphon(Re, D, H):
    """Calculate friction factor for a pipe with a twisted-tape insert using
    the Naphon correlation (2006). Only valid for turbulent flow with Re > 7000

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    D : float
        Internal diameter of tube, [m]
    H : float
        Tape pitch for twist of π radians (180º), [m]

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    # Eq 9
    f = 3.517/Re**0.414*(1+D/H)**1.045
    return f


@refDoc(__doi__, [10])
def f_twisted_Agarwal(Re, D, H):
    """Calculate friction factor for a pipe with a twisted-tape insert using
    the Agarwal and Raja Rao correlation (1996).

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    D : float
        Internal diameter of tube, [m]
    H : float
        Tape pitch for twist of π radians (180º), [m]

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    y = H/D

    if Re/y < 9 or Re/y > 1000:
        raise NotImplementedError("Input out of bound")

    # 29
    f = 1/Re/y**0.28*((75.74*(Re/y)**0.0216)**10 + (19.48*(Re/y)**0.3481)**10)

    return f

# Heat Transfer coefficient correlations
@refDoc(__doi__, [2])
def Nu_twisted_Plessis(Re, Pr, D, H, delta, Ae, De, x=None):
    """Calculate Nusselt number for a pipe with a twisted-tape insert using
    the Plessis and Kröger correlation (1987).

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    D : float
        Internal diameter of tube, [m]
    H : float
        Tape pitch for twist of π radians (180º), [m]
    delta : float
        Tape thickness, [m]
    Ae : float
        Effective flow area, [m²]
    De : float
        Effective hydraulic diameter, [m]
    x : float, optional
        Length in axial flow, [m]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    y = H/D
    A = pi*D**2/4
    Ac = A-D*delta

    # Eq 2
    Psye = (D/De)**2*Ae/A

    Ree = Re*A/Ae*De/D
    Omge = Ree/y

    if x is not None:
        # Correction for flow don't fully developed
        xe = x*Ac/Ae
        xe_ = xe/(Ree*Pr*De)
        xinf = xe_ * (1+0.04*(Omge*Pr)**3)**(1/3)
        Psye *= (1+0.153*xinf**-1.05)**(1/3)

    # Eq 14
    Nu = 1.58*Psye*(1+6.4e-5*(Omge*Pr)**3)**0.117*(1+0.002*Omge**1.4)**(1/7)

    return Nu


@refDoc(__doi__, [6])
def Nu_twisted_Hong(Re, Pr, D, H):
    """Calculate Nusselt number for a pipe with a twisted-tape insert using
    the Hong and Bergles correlation (1976). Valid only for laminar region
    Re < 2500

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    Pr : float
        Prandtl number, [-]
    D : float
        Internal diameter of tube, [m]
    H : float
        Tape pitch for twist of π radians (180º), [m]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    y = H/D

    # Eq 3
    Nu = 5.172*(1+5.484e-3*Pr**0.7*(Re/y)**1.25)**0.5
    return Nu


@refDoc(__doi__, [4, 5])
def Nu_twisted_Manglik(Re, Pr, Gr, Gz, D, H, delta, Dh, mu, muW):
    """Calculate Nusselt number for a pipe with a twisted-tape insert using
    the Manglik and Bergles correlation (1993)

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    Pr : float
        Prandtl number, [-]
    Gr : float
        Grashof number, [-]
    Gz : float
        Graetz number, [-]
    D : float
        Internal pipe diameter, [m]
    H : float
        Tape pitch for twist of π radians (180º), [m]
    delta : float
        Tape thickness, [m]
    Dh : float
        Hydraulic diameter, [m]
    mu : float
        Bulk flow temperature viscosity, [Pa·s]
    muW : float
        Wall flow temperature viscosity, [Pa·s]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """

    if Re <= 2000:
        # Laminar flow

        y = H/D
        Ra = Gr*Pr

        A = pi*D**2/4
        Ac = A-D*delta

        Resw = Re*(pi/(pi-4*delta/D)) * ((pi*D/H)**2)**0.5
        Sw = Resw/(H/2/D)**0.5
        Reax = Resw/Vs*Va*(1+(pi/2/y)**2)**0.5

        # Eq 17
        Nu = 4.612*(mu/muW)**0.14*(
            ((1+0.0951*Gz**0.894)**2.5 + 6.413e-9*(Sw*Pr**0.391)**3.835)**2
            + 2.132e-14*(Reax*Ra)**2.23)**0.1

    elif Re >= 5000:
        # Transition flow
        # Eq 8-9
        Nu = (1+0.769/y)*0.023*Re**0.8*Pr**0.4*(pi/(pi-4*delta/D))**0.8 \
            * ((pi+2-2*delta/D)/(pi-4*delta/D))**0.2

        if mu < muW:
            # Cooling
            n = 0.3
        else:
            n = 0.18
        Nu *= (mu/muW)**n

    else:
        NuL = Nu_twisted_Manglik(2000, Pr, Gr, Gz, D, H, delta, Dh, mu, muW)
        NuH = Nu_twisted_Manglik(2000, Pr, Gr, Gz, D, H, delta, Dh, mu, muW)
        Nu = NuL*(5000-Re)/3000 + NuH*(Re-2000)/3000
    return Nu


@refDoc(__doi__, [7, 8])
def Nu_twisted_Lopina(Re, Pr, D, H, Dh, mu, muW, beta, DT, HTRI=False):
    """Calculate Nusselt number for a pipe with a twisted-tape insert using
    the Lopina and Bergles correlation (1969). Only valid for turbulent flow
    with Re > 5000

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    Pr : float
        Prandtl number, [-]
    D : float
        Internal diameter of tube, [m]
    H : float
        Tape pitch for twist of π radians (180º), [m]
    Dh : float
        Hydraulic diameter, [m]
    mu : float
        Bulk flow temperature viscosity, [Pa·s]
    muW : float
        Wall flow temperature viscosity, [Pa·s]
    beta : float, optional
        Volumetric expansion coefficient, [1/K]
    dT : float, optional
        Temperature difference between bulk and wall, [K]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    Reh = Re*Dh/D
    y = H/D
    alpha = (1+(pi/y)**2)**0.5

    if HTRI:
        # Use parameter used in HTRI manual
        F = 1.03
        x = 0.79
    else:
        F = 1.137
        x = 0.8

    if mu > muW:
        # In liquids viscosity decrease with temperature, so cooling processes
        Nc = 0
    else:
        # Heating processes
        Nc = 0.193*((2*Reh/y)**2*Dh/D*beta*DT*Pr)**(1/3)

    Nu = F*(0.023*(alpha*Reh)**x * Pr**0.4 + Nc)
    return Nu


@refDoc(__doi__, [8])
def Nu_twisted_HTRI(Re, Pr, D, H, Dh, mu, muW, beta=None, dT=None, L=None):
    """Calculate Nusselt number for a pipe with a twisted-tape insert using
    the correlation used in HTRI® software.

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    Pr : float
        Prandtl number, [-]
    Gz : float
        Graetz number, [-]
    D : float
        Internal pipe diameter, [m]
    H : float
        Tape pitch for twist of π radians (180º), [m]
    Dh : float
        Hydraulic diameter, [m]
    mu : float
        Bulk flow temperature viscosity, [Pa·s]
    muW : float
        Wall flow temperature viscosity, [Pa·s]
    L : float, optional
        Length of heated pipe, [m]
    beta : float, optional
        Volumetric expansion coefficient, [1/K]
    dT : float, optional
        Temperature difference between bulk and wall, [K]

    Returns
    -------
    Nu : float
        Nusselt number, [-]

    Examples
    --------
    B3.1.2.4 turbulent flow
    >>> beta = -2/(527+609)*(527-609)/(106-36.7)
    >>> args =(23390, 4.41, 0.0158, 0.158, 0.0096, 0.000208, 1e-3, beta, 69.3)
    >>> Nu = Nu_twisted_HTRI(*args)
    >>> print("%0.0f" % (Nu*0.11/0.0096))
    1306

    laminar flow
    >>> args =(1656, 8.69, 0.0211, 0.211, None, 0.000343, None, None, 12.2)
    >>> Nu = Nu_twisted_HTRI(*args)
    >>> print("%0.1f" % (Nu*0.107/0.0211))
    157.8
    """
    y = H/D

    if Re <= 2000:
        # Laminar flow

        if Pr <= 200:
            # Hong-Bergles modified correlation
            Nu = 5.172*(1+5.484e-3*Pr**0.7*(2*Re/y)**1.25)**0.5
        elif Pr >= 2000:
            # Marner-Bergles correlation
            Nu = 1.322*(pi/4*Re*Pr*D/L)**0.458*(mu/muW)**0.14
        else:
            NuL = Nu_twisted_HTRI(Re, 200, D, H, Dh, mu, muW)
            NuH = Nu_twisted_HTRI(Re, 2000, D, H, Dh, mu, muW)
            # Eq B3.1-42
            Nu = Pr*(NuH-NuL)/1800 + NuL - (NuH-NuL)/9

    elif Re >= 5000:
        Nu = Nu_twisted_Lopina(Re, Pr, D, H, Dh, mu, muW, beta, dT, HTRI=True)

    else:
        NuL = Nu_twisted_HTRI(2000, Pr, D, H, Dh, mu, muW)
        NuH = Nu_twisted_HTRI(5000, Pr, D, H, Dh, mu, muW)
        Nu = NuL*(5000-Re)/3000 + NuH*(Re-2000)/3000

    return Nu


@refDoc(__doi__, [9])
def Nu_twisted_Naphon(Re, Pr, D, H):
    """Calculate Nusselt number for a pipe with a twisted-tape insert using
    the Naphon correlation (2006). Only valid for turbulent flow with Re > 7000

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    Pr : float
        Prandtl number, [-]
    D : float
        Internal diameter of tube, [m]
    H : float
        Tape pitch for twist of π radians (180º), [m]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    # Eq 8
    Nu = 0.648*Re**0.36*Pr**(1/3)*(1+D/H)**2.475
    return Nu


@refDoc(__doi__, [10])
def Nu_twisted_Agarwal(Re, Pr, D, H, mu, muW):
    """Calculate Nusselt number for a pipe with a twisted-tape insert using
    the Agarwal and Raja Rao correlation (1996).

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    Pr : float
        Prandtl number, [-]
    D : float
        Internal diameter of tube, [m]
    H : float
        Tape pitch for twist of π radians (180º), [m]
    mu : float
        Bulk flow temperature viscosity, [Pa·s]
    muW : float
        Wall flow temperature viscosity, [Pa·s]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    y = H/D

    if Re/y < 9 or Re/y > 1000:
        raise NotImplementedError("Input out of bound")

    if mu > muW:
        # In liquids viscosity decrease with temperature, so cooling processes
        Nu = 1.365*Re**0.517*y**-1.05*Pr**(1/3)*(mu/muW)**0.14          # Eq 31
    else:
        # Heating processes
        Nu = 0.725*Re**0.568*y**-0.788*Pr**(1/3)*(mu/muW)**0.14         # Eq 30

    return Nu



class TwistedTape():
    """Twisted-tape insert used in heat exchanger to improve efficiency.
    This tape, generally a thin metal strip, is twisted about its longitudinal
    axis"""

    TEXT_FRICTION = (
        "Manglik-Bergles (1993)",
        "Plassis-Kröger (1984)",
        "Lopina-Bergles (1969)",
        "Shah-London (1978)",
        "Naphon (2006)",
        "Agarwal-Rao (1996)")

    TEXT_HEAT = (
        "HTRI",
        "Lopina-Bergles (1969)",
        "Manglik-Bergles (1993)",
        "Plessis-Kröger (1987)",
        "Hong-Bergles (1976)",
        "Naphon (2006)",
        "Agarwal-Rao (1996)",t
        "Kidd (1969)")

    def __init__(self, H, D, delta):
        """
        Definition of twisted tape accesory

        Parameters
        ----------
        H : float
            Tape pitch for twist of π radians (180º), [m]
        D : float
            Internal diameter of tube, [m]
        delta : float
            Tape thickness, [m]
        """
        # Geometrical definition of parameters in [1]_

        # Helical factor, Eq 2
        self.G = Dimensionless((1+pi**2/D**2/4/H**2)**0.5)

        # Effective cross-sectional flow area, Eq 3
        self.Ae = Area(2*H**2/pi*(self.G-1) - D*delta)

        # Effective wetted perimeter, Eq 4
        self.Pe = Length(2 * (D - delta + pi*D/2/self.G))

        # Effective hydraulic diameter, Eq 7
        self.De = Length(4*self.Ae/self.Pe)

        # Area tube without tape
        self.A = pi*D**2/4

        # Tape twist parameter
        self.y = Dimensionless(H/D)


class UI_TwistedTape(QtWidgets.QWidget):

    """Custom widget to define DIPPR equation input"""
    def __init__(self, parent=None):
        super().__init__(parent)
        lyt = QtWidgets.QGridLayout(self)
        self.check = QtWidgets.QCheckBox(self.tr("Use twisted tape insert"))
        self.check.toggled.connect(self.setEnabled)
        lyt.addWidget(self.check, 1, 1)
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

        lytH = QtWidgets.QHBoxLayout()
        lytH.addWidget(QtWidgets.QLabel(self.tr("Friction factor calculation method")))
        self.delta = Entrada_con_unidades(Length, "Thickness")
        self.methodFriction = QtWidgets.QComboBox()
        for method in TwistedTape.TEXT_FRICTION:
            self.methodFriction.addItem(method)
        lytH.addWidget(self.methodFriction)
        lyt.addLayout(lytH, 5, 1, 1, 2)

        lytH = QtWidgets.QHBoxLayout()
        lytH.addWidget(QtWidgets.QLabel(self.tr("Heat transfer calculation method")))
        self.methodHeat = QtWidgets.QComboBox()
        for method in TwistedTape.TEXT_HEAT:
            self.methodHeat.addItem(method)
        lytH.addWidget(self.methodHeat)
        lyt.addLayout(lytH, 6, 1, 1, 2)

        lyt.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 10, 3)

        # self.fill()

    def setEnabled(self, boolean):
        """Toggled enable/disable state for all children widget except
        checkbox used to change this"""
        for wdg in self.children():
            if wdg is not self.check:
                wdg.setEnabled(boolean)

    # def fill(self):
    #     self.check.setChecked(True)
        # self.check.setChecked(False)


class Dialog(QtWidgets.QDialog):
    """Component list config dialog"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(self.tr("Twisted-tape insert"))
        layout = QtWidgets.QVBoxLayout(self)
        self.datos = UI_TwistedTape()
        layout.addWidget(self.datos)
        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Cancel
            | QtWidgets.QDialogButtonBox.StandardButton.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox)

    # def value(self, config):
    #     """Function to result wizard"""
    #     config = self.datos.value(config)
    #     return config


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = Dialog()
    Dialog.show()
    sys.exit(app.exec())
