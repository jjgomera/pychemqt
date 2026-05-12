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
         "doi": ""},
    7:
        {"autor": "Kutateladze, S.S., Borishanskii, V.M. ",
         "title": "A Concise Encyclopedia of Heat Transfer",
         "ref": "Pergamon Press (1966)",
         "doi": ""},
    8:
        {"autor": "White, C.M.",
         "title": "Streamline Flow through Curved Pipes",
         "ref": "Proc. R .Soc. London A 123 (1929) 645-63",
         "doi": "10.1098/rspa.1929.0089"},
    9:
        {"autor": "Mori, Y., Nakayama, W.",
         "title": "Study on Forced Convective Heat Transfer in Curved Pipes "
                  "(1st Report, Laminar Region)",
         "ref": "Int. J. Heat Mass Transfer 8(1) (1965) 67-82",
         "doi": "10.1016/0017-9310(65)90098-0"},
    10:
        {"autor": "Mori, Y., Nakayama, W.",
         "title": "Study on Forced Convective Heat Transfer in Curved Pipes "
                  "(2nd Report, Turbulent Region)",
         "ref": "Int. J. Heat Mass Transfer 10(1) (1967) 37-59",
         "doi": "10.1016/0017-9310(67)90182-2"},
    11:
        {"autor": "Hart, J., Ellenberger, J., Hamersma, P.J.",
         "title": "Single- and Two-Phase Flow Through Helically Coiled Tubes",
         "ref": "Chem. Eng. Sci. 43(4) (1988) 775-783",
         "doi": "10.1016/0009-2509(88)80072-1"},
    12:
        {"autor": "Ju, H., Huang, Z., Xu, Y., Duan, B, Yu, Y.",
         "title": "Hydraulic Performance of Small Bending Radius Helical "
                  "Coil-Pipe",
         "ref": "J. Nuclear Sci. Eng. 38(10) (2001) 826-831",
         "doi": "10.1080/18811248.2001.9715102"},
    13:
        {"autor": "Mishra, P., Gupta, S.N.",
         "title": "Momentum Transfer in Curved Pipes. 1. Newtonian Fluids",
         "ref": "Ind. Eng. Chem. Process Des. Dev. 18(1) (1979) 130-137",
         "doi": "10.1021_i260069a017"},
    14:
        {"autor": "Czop, V., Barbier, D., Dong, S.",
         "title": "Pressure drop, void fraction and shear stress measurements "
                  "in an adiabatic two-phase flow in a coiled tube",
         "ref": "Nuclear Eng. Design 149 (1994) 323-333",
         "doi": "10.1016/0029-5493(94)90298-4"},
    15:
        {"autor": "Xin, R.C., Ebadian, M.A.  ",
         "title": "The Effects of Prandtl Numbers on Local and Average "
                  "Convective Heat Transfer Characteristics in Helical Pipes",
         "ref": "J. Heat Transfer 119(3) (1997) 467-73",
         "doi": "10.1115/1.2824120."},
    16:
        {"autor": "Seban R.A., McLaughlin, E.F.",
         "title": "Heat Transfer in Tube Coils with Laminar and Turbulent Flow",
         "ref": "Int. J. Heat Mass Transfer 6() (1963) 387-395",
         "doi": "10.1016/0017-9310(63)90100-5"},
    17:
        {"autor": "Manlapaz, R.L., Churchill, S.W.",
         "title": "Fully Developed Laminar Flow in a Helically Coiled Tube of "
                  "Finite Pitch",
         "ref": "Chem. Eng. Communications 7 (1980) 57-78",
         "doi": "10.1080/00986448008912549"},

    # 18:
        # {"autor": "",
         # "title": "",
         # "ref": "",
         # "doi": ""},

}


# Critical Reynolds number correlations
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


@refDoc(__doi__, [7])
def Rec_Kutateladze(di, Dc):
    r"""Calculates critical Reynolds to define transition between laminar and
    turbulent flow using using the correlation of Kutateladze (1966).

    .. math::
        Re_c = 2300 + 10500 \left(\frac{d_i}{D_c}\right)^{0.3}

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
    # Eq 7.26
    Rec = 2300 + 10500*(di/Dc)**0.3
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


@refDoc(__doi__, [9, 10])
def f_MoriNakayama(Re, di, Dc):
    r"""Calculates friction factor for internal flow of a helical coil in
    laminar flow using the method of Mori-Nakayama (1965).

    .. math::
        \frac{f_c}{f_{s}=\left(\frac{0.108De^{0.5}}{1-3.253 De^{-0.5}}\right)

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

    Rec = Rec_Ito(di, Dc)

    # Limit between both turbulent phases, Eq. 43
    Re_ = 6.5e5 * (di/Dc)**0.5

    if Re < Rec:
        # Laminar flow
        De = Dean(Re, di, Dc)
        fd = f_friccion(Re)

        fI = 0.108*De**0.5                                            # Eq 1.33
        fII = fI / (1-3.253/De**0.5)                                  # Eq 1.34
        f = fII * fd
    elif Re < Re_:
        # Low turbulent region, # Eq 40 from [10]_
        f = 0.3/(Re*(di/Dc)**2)**0.2*(1+0.112/(Re*(di/Dc)**2)**0.2)/(di/Dc)**0.5
    else:
        # Hith turbulent region, # Eq 41 from [10]_
        f = 0.192/(Re*(di/Dc)**2.5)**(1/6) * \
            (1+0.068/(Re*(di/Dc)**2.5)**(1/6)) / (di/Dc)**0.5

    return f


@refDoc(__doi__, [12])
def f_Ju(Re, di, Dc):
    r"""Calculates friction factor for internal flow of a helical coil using
    the method of Ju et al. (2001).

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
    Rec = Rec_Srinivasan(di, Dc)

    De = Dean(Re, di, Dc)
    fd = f_friccion(Re)

    if De < 11.6:
        # Laminar flow, Eq 12, straight tube correlation
        f = fd
    elif Re < Rec:
        # Laminar with big vortex, Eq. 13
        f = fd * (1 + 0.015*Re**0.75*(di/Dc)**0.4)
    else:
        # Turbulent flow, Eq 14
        f = fd * (1 + 0.011*Re**0.23*(di/Dc)**0.14)

    return f


@refDoc(__doi__, [13])
def f_MishraGupta(Re, di, Dc):
    r"""Calculates friction factor for internal flow of a helical coil using
    the method of Mishra-Gupta (1979).

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
    Rec = Rec_Ito(di, Dc)

    fd = f_friccion(Re)

    if Re < Rec:
        # Laminar flow, Eq 5.
        De = Dean(Re, di, Dc)
        f = fd * (1 + 0.033*log(De)**4)
    else:
        # Turbulent flow, Eq 10.
        f = fd + 0.03*(di/Dc)**0.5

    return f

@refDoc(__doi__, [17])
def f_laminar_ManlapazChurchill(Re, di, Dc, p):
    r"""Calculates friction factor in laminar regimen for internal flow of a
    helical coil using the method of Manlapaz-Churchill (1980).

    .. math::
        \frac{f_c}{f_{s,L}} = \left[\left(1 -
        \frac{0.18}{\left(1+\left(\frac{35}{De}\right)^2\right)^{1/2}}\right)^m
        + \left(1+\frac{d_i}{3 D_c}\right)^2 \frac{De}{88.33}\right]^{1/2}

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    di : float
        Inner diameter of the pipe, [m]
    Dc : float
        Diameter of the helix, [m]
    p : float
        Pitch for twist of 2π radians (360º), [m]

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    fd = f_friccion(Re)
    De = Dean(Re, di, Dc)

    if De <= 20:
        m = 2
    elif De <= 40:
        m = 1
    else:
        m = 0

    if p:
        # Eq 25
        X = De*(1/(1+(p/Dc/2/pi)**2))**0.5
    else:
        X = De

    # Eq 29
    f = fd * ((1-0.18/(1+(35/X)**2)**0.5)**m + (1+di/Dc/3)**2*X/88.33)**0.5

    return f



@refDoc(__doi__, [8])
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
    De = Dean(Re, di, Dc)
    fd = f_friccion(Re)
    if De > 11.6:
        C = 1 - (1 - (11.6/De)**0.45)**(1/0.45)
    else:
        C = 1
    return fd/C


@refDoc(__doi__, [11, 2])
def f_laminar_Hart(Re, di, Dc):
    r"""Calculates friction factor for internal flow of a helical coil in
    laminar flow using the method of Hart (1988).

    .. math::
        f_c = \frac{f_{s,L}} {1 + 0.09 \frac{De^{1.5}}{70+De}

    Recomended method in [2]_ for friction factor in laminar flow.

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
    De = Dean(Re, di, Dc)
    fd = f_friccion(Re)

    f = 1 + 0.09 * De**1.5 / (70+De)
    return f * fd


@refDoc(__doi__, [14])
def f_turbulent_Czop(Re, di, Dc):
    r"""Calculates friction factor for internal flow of a helical coil in
    turbulent flow using the method of Czop (1994).

    .. math::
        f_c =  \frac{0.096}{De^{-1517}}

    The paper give this correlation for single phase flow. Give too
    correlations for two phase flow.

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
    De = Dean(Re, di, Dc)

    # Eq 9
    f = 0.096 / De**0.1517
    return f



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


@refDoc(__doi__, [9, 10])
def Nu_MoriNakayama(Re, Pr, di, Dc):
    r"""Calculates Nusselt number for internal flow of a helical coil in
    laminar flow using the method of Mori-Nakayama (1965).

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
    Rec = Rec_Ito(di, Dc)

    if Re < Rec:
        # Laminar flow
        De = Dean(Re, di, Dc)

        if Pr >= 1:
            # Eq 2.15
            Z = 2/11*(1+(1+77/4/Pr**2)**0.5)
        else:
            # Eq 2.18
            Z = (2+(10/Pr**2-1)**0.5)/5

        # Eq 2.23
        NuI = 0.1979*De**0.5/Z

        if Pr >= 1:
            # Eq 2.24
            f = 1 + 37.05/Z * (1/40 - 17/120*Z + (1/10/Z + 13/30)/10/Pr)*De**-0.5
        else:
            # Eq 2.25
            f = 1 - 37.05/Z * (Z**2/12 + 1/24 - 1/120/Z
                               - (4/3*Z - 1/3/Z + 1/15/Z**2)/20/Pr)*De**-0.5
        Nu = 48/11 * NuI/f

    else:
        # Turbulent flow
        # Eq 91 in [10]_
        Nu = Pr/(26.2*(Pr**(2/3)-0.074)) * Re**0.8 * (di/Dc)**0.1 * \
            (1+0.098/(Re*(di/Dc)**2)**0.2)

    return Nu


@refDoc(__doi__, [15])
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


@refDoc(__doi__, [16])
def Nu_SebanMcLaughlin(Re, Pr, di, Dc):
    r"""Calculates Nusselt number for internal flow of a helical coil using the
    correlation of Seban-McLaughlin (1963)

    For laminar flow:

    .. math::
        Nu = 1.04 \left(\frac{Re}{1-\left(1-\left(1-\frac{11.6}{De}\right)
        ^{0.45}\right)^{1/0.45}}\right)^{1/3} Pr^{1/3}

    For turbulent flow:

    .. math::
        Nu = 0.023 Re^{0.85} Pr^{0.4} \left(\frac{d_i}{D_c}\right)^{0.1}

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
    Rec = Rec_Ito(di, Dc)

    if Re < Rec:
        # Laminar flow
        # Use Whie correlation for friction factor
        f = f_laminar_White(Re, di, Dc)

        # Eq 3
        Nu = 0.13*(f/8*Re**2)**(1/3)*Pr**(1/3)

    else:
        # Use friction factor for a straight tube given in paper
        fs = 0.023 / Re**0.2

        # Eq 4, Friction factor
        f = fs * (Re*(di/Dc)**2)**0.05

        # Eq 6
        Nu = f * Re * Pr**0.4

    return Nu


class Helical(CallableEntity):
    """Helical coil tube used as anhancing heat transfer equipment.

    Parameters
    ----------
    di : float
        Inner diameter of the pipe, [m]
    Dc : float
        Diameter of the helix, [m]
    p : float, optional
        Pitch for twist of 2π radians (360º), [m]
    """

    TEXT_REYNOLDS_CRITICAL = (
        "Schmidt (1967)",
        "Ito (1959)",
        "Kubair-Kuloor (1966)",
        "Srinivasan (1968)",
        "Kutateladze (1966)")

    TEXT_LAMINAR_FRICTION = (
        "Schmidt (1967)",
        "White (1929)",
        "Mori-Nakayama (1965)",
        "Hart (1988)",
        "Ju (2001)",
        "Mishra-Gupta (1979)",
        "Manlapaz-Churchill (1980)",
    )

    TEXT_TURBULENT_FRICTION = (
        "Schmidt (1967)",
        "Mori-Nakayama (1965)",
        "Ju (2001)",
        "Mishra-Gupta (1979)",
        "Czop (1994)",
    )

    TEXT_LAMINAR_HEAT = (
        "Schmidt (1967)",
        "Xin-Ebadian (1997)",
        "Mori-Nakayama (1965)",
        "Seban-McLaughlin (1963)",
    )

    TEXT_TURBULENT_HEAT = (
        "Schmidt (1967)",
        "Xin-Ebadian (1997)",
        "Mori-Nakayama (1965)",
        "Seban-McLaughlin (1963)",
    )

    status = 0
    msg = ""
    kw = {
        "methodReCritic": 0,
        "methodFrictionLaminar": 0,
        "methodFrictionTurbulent": 0,
        "methodHeatLaminar": 0,
        "methodHeatTurbulent": 0,

        "di": 0,
        "Dc": 0,
        "p": 0
    }

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

        elif self.kw["methodReCritic"] == 3:
            # Kutateladze (1966)
            Rec = Rec_Kutateladze(self.di, self.Dc)

        else:
            # Schmidt (1967)
            Rec = Rec_Schmidt(self.di, self.Dc)

        return Rec

    def Nu(self, Re, Pr):
        """Calculate nusselt number"""
        Rec = self.ReCritical

        if Re < Rec:
            # Laminar flow
            if self.kw["methodHeatLaminar"] == 1:
                # Xin-Ebadian (1997)
                Nu = Nu_XinEbadian(Re, Pr, self.di, self.Dc)

            elif self.kw["methodHeatLaminar"] == 2:
                # Mori-Nakayama (1965)
                Nu = Nu_MoriNakayama(Re, Pr, self.di, self.Dc)

            elif self.kw["methodHeatLaminar"] == 3:
                # Seban-McLaughlin (1963)
                Nu = Nu_SebanMcLaughlin(Re, Pr, self.di, self.Dc)

            else:
                # Schmidt (1967)
                Nu = Nu_Schmidt(Re, Pr, self.di, self.Dc)

        else:
            # Turbulent flow
            if self.kw["methodHeatTurbulent"] == 1:
                # Xin-Ebadian (1997)
                Nu = Nu_XinEbadian(Re, Pr, self.di, self.Dc)

            elif self.kw["methodHeatTurbulent"] == 2:
                # Mori-Nakayama (1965)
                Nu = Nu_MoriNakayama(Re, Pr, self.di, self.Dc)

            elif self.kw["methodHeatTurbulent"] == 3:
                # Seban-McLaughlin (1963)
                Nu = Nu_SebanMcLaughlin(Re, Pr, self.di, self.Dc)

            else:
                # Schmidt (1967)
                Nu = Nu_Schmidt(Re, Pr, self.di, self.Dc)

        return Nu

    def f(self, Re):
        """Calculate friction factor"""
        msg = ""
        Rec = self.ReCritical

        if Re < Rec:
            # Laminar flow
            if self.kw["methodFrictionLaminar"] == 1:
                # White (1929)
                f = f_laminar_White(Re, self.di, self.Dc)

            elif self.kw["methodFrictionLaminar"] == 2:
                # Mori-Nakayama (1965)
                f = f_MoriNakayama(Re, self.di, self.Dc)

            elif self.kw["methodFrictionLaminar"] == 3:
                # Hart (1988)
                f = f_laminar_Hart(Re, self.di, self.Dc)

            elif self.kw["methodFrictionLaminar"] == 4:
                # Ju (2001)
                f = f_Ju(Re, self.di, self.Dc)

            elif self.kw["methodFrictionLaminar"] == 5:
                # Mishra-Gupta (1979)
                f = f_MishraGupta(Re, self.di, self.Dc)

            elif self.kw["methodFrictionLaminar"] == 6:
                # Manlapaz-Churchill (1980)
                f = f_laminar_ManlapazChurchill(
                    Re, self.di, self.Dc, self.kw["p"])
                if not self.kw["p"]:
                    msg = "Helical pitch undefined, using Manlapaz correlation"
                    msg += "with Dean number"

            else:
                # Schmidt (1967)
                f = f_Schmidt(Re, self.di, self.Dc)

        else:
            # Turbulent flow
            if self.kw["methodFrictionTurbulent"] == 1:
                # Mori-Nakayama (1965)
                f = f_MoriNakayama(Re, self.di, self.Dc)

            elif self.kw["methodFrictionTurbulent"] == 2:
                # Ju (2001)
                f = f_Ju(Re, self.di, self.Dc)

            elif self.kw["methodFrictionTurbulent"] == 3:
                # Mishra-Gupta (1979)
                f = f_MishraGupta(Re, self.di, self.Dc)

            elif self.kw["methodFrictionTurbulent"] == 4:
                # Czop (1994)
                f = f_turbulent_Czop(Re, self.di, self.Dc)

            else:
                # Schmidt (1967)
                f = f_Schmidt(Re, self.di, self.Dc)

        if msg:
            self.status = 3
            self.msg = translate("equipment", msg)
            self.inputChanged.emit(self)

        return f


class UI_Helical(ToolGui):
    """Helical coil dialog"""

    title = translate("equipment", "Use helical coil")

    def loadUI(self):
        """Add widget"""
        self.Entity = Helical()

        lyt = self.wdg.layout()

        groupMethods = QtWidgets.QWidget()
        lytM = QtWidgets.QGridLayout(groupMethods)
        lytM.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 1, 0)
        lbl = QtWidgets.QLabel(self.tr("Laminar Flow"))
        lbl.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter
                         | QtCore.Qt.AlignmentFlag.AlignVCenter)
        lytM.addWidget(lbl, 1, 2)
        lbl = QtWidgets.QLabel(self.tr("Turbulent Flow"))
        lbl.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter
                         | QtCore.Qt.AlignmentFlag.AlignVCenter)
        lytM.addWidget(lbl, 1, 3)
        lytM.addWidget(QtWidgets.QLabel(
            self.tr("Friction factor method")), 2, 1)
        self.methodFrictionLaminar = QtWidgets.QComboBox()
        for method in Helical.TEXT_LAMINAR_FRICTION:
            self.methodFrictionLaminar.addItem(method)
        self.methodFrictionLaminar.currentIndexChanged.connect(
            partial(self.changeParams, "methodFrictionLaminar"))
        lytM.addWidget(self.methodFrictionLaminar, 2, 2)
        self.methodFrictionTurbulent = QtWidgets.QComboBox()
        for method in Helical.TEXT_TURBULENT_FRICTION:
            self.methodFrictionTurbulent.addItem(method)
        self.methodFrictionTurbulent.currentIndexChanged.connect(
            partial(self.changeParams, "methodFrictionTurbulent"))
        lytM.addWidget(self.methodFrictionTurbulent, 2, 3)
        lytM.addWidget(QtWidgets.QLabel(
            self.tr("Heat transfer method")), 3, 1)
        self.methodHeatLaminar = QtWidgets.QComboBox()
        for method in Helical.TEXT_LAMINAR_HEAT:
            self.methodHeatLaminar.addItem(method)
        self.methodHeatLaminar.currentIndexChanged.connect(
            partial(self.changeParams, "methodHeatLaminar"))
        lytM.addWidget(self.methodHeatLaminar, 3, 2)
        self.methodHeatTurbulent = QtWidgets.QComboBox()
        for method in Helical.TEXT_TURBULENT_HEAT:
            self.methodHeatTurbulent.addItem(method)
        self.methodHeatTurbulent.currentIndexChanged.connect(
            partial(self.changeParams, "methodHeatTurbulent"))
        lytM.addWidget(self.methodHeatTurbulent, 3, 3)
        lytM.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 4, 1)
        lyt.addWidget(groupMethods, 1, 1, 1, 2)

        lytH = QtWidgets.QHBoxLayout()
        lytH.addWidget(QtWidgets.QLabel(
            self.tr("Critical Reynolds correlation")))
        self.methodReCritic = QtWidgets.QComboBox()
        for method in Helical.TEXT_REYNOLDS_CRITICAL:
            self.methodReCritic.addItem(method)
        self.methodReCritic.currentIndexChanged.connect(
            partial(self.changeParams, "methodReCritic"))
        lytH.addWidget(self.methodReCritic)
        lytH.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed))
        lyt.addLayout(lytH, 2, 1, 1, 2)

        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 3, 1)

        lyt.addWidget(QtWidgets.QLabel(self.tr("Pipe internal diameter")), 4, 1)
        self.di = Entrada_con_unidades(Length, "PipeDiameter")
        self.di.valueChanged.connect(partial(self.changeParams, "di"))
        lyt.addWidget(self.di, 4, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Helical coil diameter")), 5, 1)
        self.Dc = Entrada_con_unidades(Length)
        self.Dc.valueChanged.connect(partial(self.changeParams, "Dc"))
        lyt.addWidget(self.Dc, 5, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Helical pitch")), 6, 1)
        self.p = Entrada_con_unidades(Length)
        self.p.valueChanged.connect(partial(self.changeParams, "p"))
        lyt.addWidget(self.p, 6, 2)

        self.Entity.valueChanged.connect(self.valueChanged.emit)
        self.Entity.inputChanged.connect(self.populate)


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
