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
from math import pi, log10

from tools.qt import QtCore, QtWidgets, translate

from lib.unidades import Dimensionless, Area, Length
from lib.utilities import refDoc
from UI.widgets import Entrada_con_unidades
from equipment.widget.gui import CallableEntity, ToolGui


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
         "doi": "10.1016/0017-9310(96)00039-7"},
    11:
        {"autor": "Kidd, G.J. Jr.",
         "title": "Heat Transfer and Pressure Drop for Nitrogen Flowing in "
                  "Tubes Containing Twisted Tapes",
         "ref": "AIChE J. 15(4) (1969) 581-585.",
         "doi": "10.1002/aic.690150420"},
    12:
        {"autor": "Smithberg, E., Landis, F.",
         "title": "Friction and Forced Convection Heat-Transfer "
                  "Characteristics in Tubes With Twisted Tape Swirl Generators",
         "ref": "J. Heat Transfer. 86(1) (1964) 39-48",
         "doi": "10.1115/1.3687060"},
    13:
        {"autor": "Sivashanmugam, P., Suresh, S.",
         "title": "Experimental studies on heat transfer and friction factor "
                  "characteristics of laminar flow through a circular tube "
                  "fitted with helical screw-tape inserts",
         "ref": "App. Thermal Eng. 26(16) (2006) 1990-1997",
         "doi": "10.1016/j.applthermaleng.2006.01.008"},
    14:
        {"autor": "Sivashanmugam, P., Suresh, S.",
         "title": "Experimental studies on heat transfer and friction factor "
                  "characteristics of turbulent flow through a circular tube "
                  "fitted with helical screw-tape inserts",
         "ref": "Chem. Eng. Processing 46(12) (2007) 1292-1298",
         "doi": "10.1016/j.cep.2006.10.009"},
    15:
        {"autor": "Sarma, P.K., Kishore, P.S., Rao, V.D., Subrahnamyam, T.",
         "title": "A Conbined approach to predict friction coefficients and "
                  "convective heat transfer characteristics in A tube with "
                  "twisted tape inserts for a wide range of Re and Pr",
         "ref": "Int. J. Therm. Sciences 44(4) (2005) 393-398",
         "doi": "10.1016/j.ijthermalsci.2004.12.001"},
    16:
        {"autor": "Murugesan, P., Mayilsamy, K., Suresh, S.",
         "title": "Heat Transfer and Friction Factor Studies in a Circular "
                  "Tube Fitted with Twisted Tape Consisting of Wire-nails",
         "ref": "Chin. J. Chem. Eng. 18(6) (2010) 1038-1042",
         "doi": "10.1016/S1004-9541(09)60166-X"},
    17:
        {"autor": "Murugesan, P., Mayilsamy, K., Suresh, S.",
         "title": "Turbulent Heat Transfer and Pressure Drop in Tube Fitted "
                  "with Square-cut Twisted Tape",
         "ref": "Chin. J. Chem. Eng. 18(4) (2010) 609-617",
         "doi": "10.1016/s1004-9541(10)60264-9"},
    18:
        {"autor": "Murugesan, P., Mayilsamy, K., Suresh, S., Srinivasan, P.S.S",
         "title": "Heat transfer and pressure drop characteristics in a "
                  "circular tube fitted with and without V-cut twisted tape"
                  "insert",
         "ref": "Int. Comm. Heat Mass Transfer 38(3) (2011) 329-334",
         "doi": "10.1016/j.icheatmasstransfer.2010.11.010"},
    19:
        {"autor": "Murugesan, P., Mayilsamy, K., Suresh, S.",
         "title": "Heat Transfer in Tubes Fitted with Trapezoidal-Cut and "
                  "Plain Twisted Tape Inserts",
         "ref": "Chem. Eng. Communications 198(7) (2011) 886-904",
         "doi": "10.1080/00986445.2011.545294"},
    20:
        {"autor": "Murugesan, P., Mayilsamy, K., Suresh, S.",
         "title": "Heat Transfer in a Tube Fitted with Vertical and "
                  "Horizontal Wing-cut Twisted Tapes",
         "ref": "Exp. Heat Transfer 25(1) (2012) 30-47",
         "doi": "10.1080/08916152.2011.559567"},
    # 21:
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


@refDoc(__doi__, [15])
def f_twisted_Sarma(Re, D, H):
    """Calculate friction factor for a pipe with a twisted-tape insert using
    the Sarma et al. correlation (2005)

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
    # Eq 4
    rhs = 0.474 - 0.3*log10(Re) + 0.065*log10(Re)**2 - 4.66e-3*log10(Re)**3
    f = (1+D/H)**3.378 * rhs
    return f


@refDoc(__doi__, [12])
def f_twisted_Smithberg(Re, D, H):
    """Calculate friction factor for a pipe with a twisted-tape insert using
    the Smithberg-Landis correlation (1964). Valid for turbulent flow

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
    # Eq 17
    n = 0.2*(1+1.7*(H/D)**-0.5)
    f = (0.046 + 2.1/(H/D-0.5)**1.2) / Re**n
    return f


@refDoc(__doi__, [16, 17, 18, 19, 20])
def f_twisted_Murugesan(Re, D, H, mod="", de=None, w=None):
    """Calculate friction factor for a pipe with a twisted-tape insert using
    the Murugesan-Mayilsamy-Suresh correlation (2010). Valid in turbulent flow

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    D : float
        Twisted tape diameter, [m]
    H : float
        Tape pitch for twist of π radians (180º), [m]
    mod : string
        Name of modification code of twisted tape
        Nails | Square cut | V cut | Trapezoidal cut
    de : float, optional
        Depth of V cut, [m]
    w : float
        Width of V cut, [m]

    Returns
    -------
    f : float
        Friction factor, [-]
    """

    if Re < 2000:
        raise NotImplementedError("Input out of bound")

    if mod == "Nails":
        # Eq 6
        f = 28.91*Re**-0.731*(H/D)**-0.255
    elif mod == "Square cut":
        # Eq 18 in 17_
        f = 6.936*Re**-0.579*(H/D)**-0.259
    elif mod == "V cut":
        # Eq 9 in 18_
        f = 8.632*Re**-0.615*(H/D)**-0.269*(1+de/D)**2.477/(1+w/D)**1.914
    elif mod == "Trapezoidal cut":
        # Eq 22 in 19_
        f = 7.401*Re**-0.587*(H/D)**-0.278
    elif mod == "Vertical wings":
        # Eq 19 in 20_
        f = 13.769*Re**-0.65*(H/D)**-0.257
    elif mod == "Horizontal wings":
        # Eq 21 in 20_
        f = 31.477*Re**-0.74*(H/D)**-0.24
    else:
        # Eq 4
        f = 2.642*Re**-0.474*(H/D)**-0.302

    return f


@refDoc(__doi__, [13, 14])
def f_helical_Sivashanmugam(Re, D, H):
    """Calculate friction factor for a pipe with a twisted-tape insert using
    the Sivashanmugam-Suresh correlation (2006)

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
    if Re > 2000:
        # Turbulent flow
        # Eq 6 in 14_
        f = 32.415*Re**-0.598*(H/D)**-0.7986

    else:
        # Laminar flow
        # Eq 6 in 13_
        f = 10.7564*Re**0.387*(H/D)**-1.054

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
        # Turbulent flow
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

    if mu >= muW:
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


@refDoc(__doi__, [11])
def Nu_twisted_Kidd(Re, Pr, D, H, L, T, Tw):
    """Calculate Nusselt number for a pipe with a twisted-tape insert using
    the Kidd correlation (1969).

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
    L : float, optional
        Length of heated pipe, [m]
    T : float
        Bulk flow temperature, [K]
    Tw : float
        Wall flow temperature, [K]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    if Re < 2e4:
        raise NotImplementedError("Input out of bound")

    y = H/D

    # Eq 3
    Nu = 0.024*Re**0.8*Pr**0.4*(T/Tw)**0.7*(1+(L/D)**-0.55)*(y/(y-1))**1.1

    return Nu


@refDoc(__doi__, [15])
def Nu_twisted_Sarma(Re, Pr, D, H):
    """Calculate Nusselt number for a pipe with a twisted-tape insert using
    the Sarma et al. correlation (2005)

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
    rhs = 0.974 - 0.783*log10(Re) + 0.35*log10(Re)**2 - 0.0273*log10(Re)**3
    Nu = Pr**(1/3) * (1+D/H)**2 * 10**rhs
    return Nu


@refDoc(__doi__, [12])
def Nu_twisted_Smithberg(Re, Pr, D, H, Dh):
    """Calculate Nusselt number for a pipe with a twisted-tape insert using
    the Smithberg-Landis correlation (1964)

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

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    f = f_twisted_Smithberg(Re, D, H)

    # Eq 38
    P = (1+0.0219/(H/D)**2/f)**0.5
    A = 50.9*D/H/Re/f**0.5 + 0.023*D/Dh/Re**0.2/Pr**(2/3)*P
    B = 1 + 700/Re/f * D/H * Dh/D * Pr**0.731
    Nu = Re*Pr*A/B

    return Nu


@refDoc(__doi__, [16, 17, 18, 19, 20])
def Nu_twisted_Murugesan(Re, Pr, D, H, mod="", de=None, w=None):
    """Calculate Nusselt number for a pipe with a twisted-tape insert using
    the Murugesan-Mayilsamy-Suresh correlation (2010). Valid in turbulent flow

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    Pr : float
        Prandtl number, [-]
    D : float
        Tape diameter, [m]
    H : float
        Tape pitch for twist of π radians (180º), [m]
    mod : string
        Name of modification code of twisted tape
        Nails|Square cut|V cut
    de : float, optional
        Depth of V cut, [m]
    w : float
        Width of V cut, [m]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """

    if Re < 2000:
        raise NotImplementedError("Input out of bound")

    if mod == "Nails":
        # Eq 5
        Nu = 0.063*Re**0.789*Pr**0.33*(H/D)**-0.257
    elif mod == "Square cut":
        # Eq 18 in 17_
        Nu = 0.041*Re**0.862*Pr**0.33*(H/D)**-0.228
    elif mod == "V cut":
        # Eq 8 in 18_
        Nu = 0.0296*Re**0.853*Pr**0.33*(H/D)**-0.222 * (1+de/D)**1.148 \
            * (1+w/D)**-0.751
    elif mod == "Trapezoidal cut":
        # Eq 21 in 19_
        Nu = 0.034*Re**0.841*Pr**0.33*(H/D)**-0.226
    elif mod == "Vertical wings":
        # Eq 18 in 20_
        Nu = 0.0484*Re**0.817*Pr**0.33*(H/D)**-0.263
    elif mod == "Horizontal wings":
        # Eq 20 in 20_
        Nu = 0.071*Re**0.775*Pr**0.33*(H/D)**-0.236
    else:
        # Eq 4
        Nu = 0.027*Re**0.862*Pr**0.33*(H/D)**-0.215

    return Nu


@refDoc(__doi__, [13, 14])
def Nu_helical_Sivashanmugam(Re, Pr, D, H):
    """Calculate Nusselt number for a pipe with a twisted-tape insert using
    the Sivashanmugam-Suresh correlation (2006).

    Correlation for helical screw-tape insert valid for all flow regimen, using
    turbulent correlation for transition regimen.

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
    # There is a extension study to regularly spaced helical screw-tape insert
    # Don't implement because laminar paper don´t give correlation

    # Sivashanmugam, P., Suresh, S.
    # Experimental studies on heat transfer and friction factor characteristics
    # of turbulent flow through a circular tube fitted with regularly spaced
    # helical screw-tape inserts
    # App. Thermal Eng. 27(8-9) (2007) 1311-1319
    # doi: 10.1016/j.applthermaleng.2006.10.035

    # Experimental studies on heat transfer and friction factor characteristics
    # of laminar flow through a circular tube fitted with regularly spaced
    # helical screw-tape inserts
    # Exp. Thermal Fluid Sci. 31(4) (2007) 301-308
    # doi: 10.1016/j.expthermflusci.2006.05.005

    if Re > 2000:
        # Turbulent flow
        # Eq 5 in 14_
        Nu = 0.4675*Re**0.4774*Pr*(H/D)**-0.2138

    else:
        # Laminar flow
        # Eq 6 in 13_
        Nu = 0.017*Re**0.996*Pr*(H/D)**-0.5437

    return Nu


class TwistedTape(CallableEntity):
    """Twisted-tape insert used in heat exchanger to improve efficiency.
    This tape, generally a thin metal strip, is twisted about its longitudinal
    axis

    Parameters
    ----------
    methodFriction : integer
        Index of method used for friction factor calculation
    methodHeat : integer
        Index of method used for heat transfer calculation
    H : float
        Tape pitch for twist of π radians (180º), [m]
    Dt : float
        Internal diameter of tube, [m]
    delta : float
        Tape thickness, [m]

    Lopina: Turbulent flow
    Naphon: Turbulent flow
    Kidd: Turbulent flow, developed with gas data
    Sivashanmugam: Laminar flow Re < 2000
    """

    TEXT_FRICTION = (
        "Manglik-Bergles (1993)",
        "Plessis-Kröger (1984)",
        "Lopina-Bergles (1969)",
        "Shah-London (1978)",
        "Naphon (2006)",
        "Agarwal-Rao (1996)",
        "Sarma (2005)",
        "Smithberg-Landis (1964)",
        "Murugesan (2010)",
        )

    TEXT_HEAT = (
        "HTRI",
        "Manglik-Bergles (1993)",
        "Plessis-Kröger (1987)",
        "Lopina-Bergles (1969)",
        "Hong-Bergles (1976)",
        "Naphon (2006)",
        "Agarwal-Rao (1996)",
        "Kidd (1969)",
        "Sarma (2005)",
        "Smithberg-Landis (1964)",
        "Murugesan (2010)",
        )

    TEXT_MURUGESAN = (
        "",
        "Nails",
        "Square cut",
        "V cut",
        "Trapezoidal cut",
        "Vertical wings",
        "Horizontal wings")

    # Helical method
    # "Sivashanmugam-Suresh (2006)"

    status = 0
    msg = ""
    kw = {
        "methodFriction": 0,
        "methodHeat": 0,
        "H": 0,
        "Dt": 0,
        "delta": 0,
        "isHelical": False,
        "modMurugesan": "",
        "Vcut_w": 0,
        "Vcut_De": 0
        }

    valueChanged = QtCore.pyqtSignal(object)
    inputChanged = QtCore.pyqtSignal(object)

    @property
    def isCalculable(self):
        """Check if all input are defined"""
        if not self.kw["H"]:
            self.msg = translate("equipment", "undefined tape pitch")
            self.status = 0
            return False
        if not self.kw["Dt"]:
            self.msg = translate("equipment", "undefined tape diameter")
            self.status = 0
            return False
        if not self.kw["delta"]:
            self.msg = translate("equipment", "undefined tape thickness")
            self.status = 0
            return False

        self.msg = ""
        self.status = 1
        return True

    def calculo(self):
        """Definition of twisted tape inserts"""
        # Geometrical definition of parameters in [1]_

        self.H = self.kw["H"]
        self.Dt = self.kw["Dt"]
        self.delta = self.kw["delta"]

        # Helical factor, Eq 2
        self.G = Dimensionless((1+pi**2/self.Dt**2/4/self.H**2)**0.5)

        # Effective cross-sectional flow area, Eq 3
        self.Ae = Area(2*self.H**2/pi*(self.G-1) - self.Dt*self.delta)

        # Effective wetted perimeter, Eq 4
        self.Pe = Length(2 * (self.Dt - self.delta + pi*self.Dt/2/self.G))

        # Effective hydraulic diameter, Eq 7
        self.De = Length(4*self.Ae/self.Pe)

        # Area tube without tape
        self.A = pi*self.Dt**2/4

        # Tape twist parameter
        self.y = Dimensionless(self.H/self.Dt)

        self.valueChanged.emit(self)

    def Nu(self, Re, Pr, mu, muW, beta, dT, L, method=None):
        """Calculate nusselt number"""

        msg = ""
        if method is None:
            method = self.kw["methodHeat"]

        print(method)
        if method == 0:
            # HTRI
            Nu = Nu_twisted_HTRI(
                Re, Pr, self.Dt, self.H, self.De, mu, muW, beta, dT, L)

        elif method == 1:
            # Manglik-Bergles (1993)
            Nu = Nu_twisted_Manglik(
                Re, Pr, Gr, Gz, self.Dt, self.H, self.delta, self.De, mu, muW)

        elif method == 2:
            # Plessis-Kröger (1987)
            Nu = Nu_twisted_Plessis(
                Re, Pr, self.Dt, self.H, self.delta, self.Ae, self.De)

        elif method == 3:
            # Lopina-Bergles (1969)
            if Re < 5000:
                Nu = self.Nu(Re, Pr, mu, muW, beta, dT, L, method=0)
                msg = "Lopina correlation only valid in turbulent flow, using "
                msg += "HTRI instead"
            else:
                Nu = Nu_twisted_Lopina(
                    Re, Pr, self.Dt, self.H, self.De, mu, muW, beta, dT)

        elif method == 4:
            # Hong-Bergles (1976)
            if Re > 2500:
                Nu = self.Nu(Re, Pr, mu, muW, beta, dT, L, method=0)
                msg = "Hong correlation only valid in laminar flow, using "
                msg += "HTRI instead"
            else:
                Nu = Nu_twisted_Hong(Re, Pr, self.Dt, self.H)

        elif method == 5:
            # Naphon (2006)
            if Re < 5000:
                Nu = self.Nu(Re, Pr, mu, muW, beta, dT, L, method=0)
                msg = "Naphon correlation only valid in laminar flow, using "
                msg += "HTRI instead"
            else:
                Nu = Nu_twisted_Naphon(Re, Pr, self.Dt, self.H)

        elif method == 6:
            # Agarwal-Rao (1996)
            try:
                Nu = Nu_twisted_Agarwal(Re, Pr, self.Dt, self.H, mu, muW)
            except NotImplementedError:
                Nu = self.Nu(Re, Pr, mu, muW, beta, dT, L, method=0)
                msg = "Agarwal correlation out of range, using HTRI instead"

        elif method == 7:
            # Kidd (1969)
            if Re < 2e4:
                Nu = self.Nu(Re, Pr, mu, muW, beta, dT, L, method=0)
                msg = "Kidd correlation only valid in turbulent flow, using "
                msg += "HTRI instead"
            else:
                Nu = Nu_twisted_Kidd(Re, Pr, self.Dt, self.H, L, 1, 1)

        elif method == 8:
            # Sarma (2005)
            Nu = Nu_twisted_Sarma(Re, Pr, self.Dt, self.H)

        elif method == 9:
            # Smithberg-Landis (1964)
            if Re < 5000:
                Nu = self.Nu(Re, Pr, mu, muW, beta, dT, L, method=0)
                msg = "Smithberg correlation only valid in turbulent flow, "
                msg += "using HTRI instead"
            else:
                Nu = Nu_twisted_Smithberg(Re, Pr, self.Dt, self.H, self.De)

        elif method == 10:
            # Murugesan (2010)
            if Re < 2000:
                Nu = self.Nu(Re, Pr, mu, muW, beta, dT, L, method=0)
                msg = "Murugesan correlation only valid in turbulent flow, "
                msg += "using HTRI instead"
            elif self.kw["modMurugesan"] == "V cut" and self.kw["Vcut_De"] \
                    and self.kw["Vcut_w"]:
                Nu = Nu_twisted_Murugesan(
                    Re, Pr, self.Dt, self.H, self.kw["modMurugesan"],
                    self.kw["Vcut_De"], self.kw["Vcut_w"])
            elif self.kw["modMurugesan"] == "V cut":
                Nu = Nu_twisted_Murugesan(Re, Pr, self.Dt, self.H)
                msg = "V cut twisted tape geometry don't defined, using plain "
                msg += "twisted tape instead"
            else:
                Nu = Nu_twisted_Murugesan(
                    Re, Pr, self.Dt, self.H, self.kw["modMurugesan"])

        return Nu

    def f(self, Re, method=None):
        """Calculate friction factor"""
        msg = ""
        if method is None:
            method = self.kw["methodFriction"]

        if self.kw["isHelical"]:
            f = f_helical_Sivashanmugam(Re, self.Dt, self.H)
            return f

        if method == 0:
            # Manglik-Bergles (1993)
            f = f_twisted_Manglik(Re, self.Dt, self.H, self.delta, self.De)

        elif method == 1:
            # Plessis-Kröger (1984)
            f = f_twisted_Plessis(
                Re, self.Dt, self.H, self.delta, self.Ae, self.De)

        elif method == 2:
            # Lopina-Bergles (1969)
            if Re < 5000:
                f = self.f(Re, 0)
                msg = "Lopina correlation only valid in turbulent flow, using "
                msg += "Manglik instead"
            else:
                f = f_twisted_Lopina(Re, self.Dt, self.H, self.De)

        elif method == 3:
            # Shah-London (1978)
            f = f_twisted_Shah(Re, self.Dt, self.H, self.delta)

        elif method == 4:
            # Naphon (2006)
            if Re < 5000:
                f = self.f(Re, 0)
                msg = "Naphon correlation only valid in turbulent flow, using "
                msg += "Manglik instead"
            else:
                f = f_twisted_Naphon(Re, self.Dt, self.H)

        elif method == 5:
            # Agarwal-Rao (1996)
            try:
                f = f_twisted_Agarwal(Re, self.Dt, self.H)
            except NotImplementedError:
                f = self.f(Re, 0)
                msg = "Agarwal correlation out of range, using Manglik instead"

        elif method == 6:
            # Sarma (2005)
            f = f_twisted_Sarma(Re, self.Dt, self.H)

        elif method == 7:
            # Smithberg-Landis (1964)
            if Re < 3000:
                f = self.f(Re, 0)
                msg = "Smithberg correlation only valid in turbulent flow, "
                msg += "using Manglik instead"
            else:
                f = f_twisted_Smithberg(Re, self.Dt, self.H)

        elif method == 8:
            # Murugesan (2010)
            if Re < 2000:
                f = self.f(Re, 0)
                msg = "Murugesan correlation only valid in turbulent flow, "
                msg += "using Manglik instead"
            elif self.kw["modMurugesan"] == "V cut" and self.kw["Vcut_De"] \
                    and self.kw["Vcut_w"]:
                f = f_twisted_Murugesan(
                    Re, self.Dt, self.H, self.kw["modMurugesan"],
                    self.kw["Vcut_w"], self.kw["Vcut_De"])
            elif self.kw["modMurugesan"] == "V cut":
                f = f_twisted_Murugesan(Re, self.Dt, self.H)
                msg = "V cut twisted tape geometry don't defined, using plain "
                msg += "twisted tape instead"
            else:
                f = f_twisted_Murugesan(
                    Re, self.Dt, self.H, self.kw["modMurugesan"])

        if msg:
            self.status = 3
            self.msg = translate("equipment", msg)
            self.inputChanged.emit(self)

        return f


class UI_TwistedTape(ToolGui):
    """Twisted-tape insert dialog"""

    title = translate("equipment", "Use twisted tape insert")

    def loadUI(self):
        """Add widget"""
        self.Entity = TwistedTape()

        lyt = self.layout()

        lytH = QtWidgets.QHBoxLayout()
        lytH.addWidget(QtWidgets.QLabel(
            self.tr("Friction factor calculation method")))
        self.methodFriction = QtWidgets.QComboBox()
        for method in TwistedTape.TEXT_FRICTION:
            self.methodFriction.addItem(method)
        self.methodFriction.currentIndexChanged.connect(
            partial(self.changeParams, "methodFriction"))
        self.methodFriction.currentTextChanged.connect(self.setVisibleMod)
        lytH.addWidget(self.methodFriction)
        lyt.addLayout(lytH, 2, 1, 1, 2)

        lytH = QtWidgets.QHBoxLayout()
        lytH.addWidget(QtWidgets.QLabel(
            self.tr("Heat transfer calculation method")))
        self.methodHeat = QtWidgets.QComboBox()
        for method in TwistedTape.TEXT_HEAT:
            self.methodHeat.addItem(method)
        self.methodHeat.currentIndexChanged.connect(
            partial(self.changeParams, "methodHeat"))
        self.methodHeat.currentTextChanged.connect(self.setVisibleMod)
        lytH.addWidget(self.methodHeat)
        lyt.addLayout(lytH, 3, 1, 1, 2)

        label = QtWidgets.QLabel(self.tr("Tape pitch"))
        label.setToolTip(self.tr("Tape pitch for twist of π radians (180º)"))
        lyt.addWidget(label, 4, 1)
        self.H = Entrada_con_unidades(Length)
        self.H.valueChanged.connect(partial(self.changeParams, "H"))
        lyt.addWidget(self.H, 4, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Tape diameter")), 5, 1)
        self.Dt = Entrada_con_unidades(Length, "PipeDiameter")
        self.Dt.valueChanged.connect(partial(self.changeParams, "Dt"))
        lyt.addWidget(self.Dt, 5, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Tape thickness")), 6, 1)
        self.delta = Entrada_con_unidades(Length, "Thickness")
        self.delta.valueChanged.connect(partial(self.changeParams, "delta"))
        lyt.addWidget(self.delta, 6, 2)

        self.groupMurugesan = QtWidgets.QWidget()
        lytMuru = QtWidgets.QGridLayout(self.groupMurugesan)
        lytMuru.addWidget(QtWidgets.QLabel(self.tr(
            "Murugesan correlation modification")), 1, 1, 1, 2)
        self.modMurugesan = QtWidgets.QComboBox()
        for method in TwistedTape.TEXT_MURUGESAN:
            self.modMurugesan.addItem(method)
        self.modMurugesan.currentTextChanged.connect(
            partial(self.changeParams, "modMurugesan"))
        self.modMurugesan.currentTextChanged.connect(self.setEnable_Murugesan)
        lytMuru.addWidget(self.modMurugesan, 1, 3)
        self.lblDe = QtWidgets.QLabel(self.tr("Depth of V cut"))
        lytMuru.addWidget(self.lblDe, 2, 2)
        self.De = Entrada_con_unidades(Length, "thickness")
        self.De.valueChanged.connect(partial(self.changeParams, "Vcut_De"))
        lytMuru.addWidget(self.De, 2, 3)
        self.lblw = QtWidgets.QLabel(self.tr("Widgth of V cut"))
        lytMuru.addWidget(self.lblw, 3, 2)
        self.w = Entrada_con_unidades(Length, "thickness")
        self.w.valueChanged.connect(partial(self.changeParams, "Vcut_w"))
        lytMuru.addWidget(self.w, 3, 3)

        lyt.addWidget(self.groupMurugesan, 8, 1, 1, 2)

        self.helical = QtWidgets.QCheckBox(self.tr(
            "Helical screw-tape inserts"))
        self.helical.toggled.connect(
            partial(self.changeParams, "isHelical"))
        lyt.addWidget(self.helical, 9, 1, 1, 3)

        self.Entity.valueChanged.connect(self.valueChanged.emit)
        self.Entity.inputChanged.connect(self.populate)
        self.setVisibleMod()

    def setVisibleMod(self):
        if self.methodHeat.currentText() == "Murugesan (2010)" or \
                self.methodFriction.currentText() == "Murugesan (2010)":
            self.groupMurugesan.setVisible(True)
        else:
            self.groupMurugesan.setVisible(False)
        self.setEnable_Murugesan(self.modMurugesan.currentText())

    def setEnable_Murugesan(self, mod):
        """Change Enable/Disable state for Murugesan aditional parameters"""
        self.lblDe.setEnabled(mod == "V cut")
        self.De.setEnabled(mod == "V cut")
        self.lblw.setEnabled(mod == "V cut")
        self.w.setEnabled(mod == "V cut")


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


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = Dialog()
    Dialog.show()
    sys.exit(app.exec())
