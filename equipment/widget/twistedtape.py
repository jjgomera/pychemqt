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
from math import atan, exp, log10, pi, tan

from tools.qt import QtCore, QtWidgets, translate

from lib.unidades import Dimensionless, Area, Length, Angle
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
    21:
        {"autor": "Jaisankar, S., Radhakrishnan, T.;., Sheeba, K.N.",
         "title": "Experimental studies on heat transfer and friction factor "
                  "characteristics of forced circulation solar water heater "
                  "system fitted with helical twisted tapes",
         "ref": "Solar Energy 83(11) (2009) 1943-1952",
         "doi": "10.1016/j.solener.2009.07.006"},
    22:
        {"autor": "Sivashanmugam, P., Suresh, S.",
         "title": "Experimental studies on heat transfer and friction factor "
                  "characteristics of turbulent flow through a circular tube "
                  "fitted with regularly spaced helical screw-tape inserts",
         "ref": "App. Thermal Eng. 27(8-9) (2007) 1311-1319",
         "doi": "10.1016/j.applthermaleng.2006.10.035"},
    23:
        {"autor": "Ibrahim, E.Z.",
         "title": "Augmentation of laminar flow and heat transfer in flat "
                  "tubes by means of helical screw-tape inserts",
         "ref": "Energy Conv. Management 52(1) (2011) 250-257",
         "doi": "10.1016/j.enconman.2010.06.065"},
    24:
        {"autor": "Saha, S.K., Gaitonde, U.N., Date, A.W.",
         "title": "Heat Transfer and Pressure Drop Characteristics of Laminar "
                  "Flow in a Circular Tube Fitted with Regularly Spaced "
                  "Twisted-Tape Elements",
         "ref": "Exp. Thermal Fluid Sci. 2(3) (1989) 310-322",
         "doi": "10.1016/0894-1777(89)90020-4"},
    25:
        {"autor": "Date, A.W., Gaitonde, U.N.",
         "title": "Development of Correlations for Predicting Characteristics "
                  "of Laminar Flow in a Tube Fitted with Regularly Spaced "
                  "Twisted-Tape Elements",
         "ref": "Exp. Thermal Fluid Sci. 3(4) (1990) 373-382",
         "doi": "10.1016/0894-1777(90)90035-6"},
    26:
        {"autor": "Klaczak, A.",
         "title": "Heat transfer by laminar flow in a vertical pipe with "
                  "twisted-tape inserts",
         "ref": "Heat Mass Transfer 36 (2000) 195-199",
         "doi": "10.1007/s002310050384"},
    27:
        {"autor": "Chang, S.W., Guo, M.H.",
         "title": "Thermal perfomances of enhanced smooth and spiky twisted "
                  "tapes for laminar and turbulent tubular flows",
         "ref": "Int. J. Heat Mass Transfer 55(25-26) (2012) 7651-7667",
         "doi": "10.1016/j.ijheatmasstransfer.2012.07.077"},
    28:
        {"autor": "Chang, S.W., Jan, Y.J., Liou, J.S.",
         "title": "Turbulent heat transfer and pressure drop in tube fitted"
                  "with serrated twisted tape",
         "ref": "Int. J. Thermal Sci. 46(5) (2007) 506-518",
         "doi": "10.1016/j.ijthermalsci.2006.07.009"},
    29:
        {"autor": "Chang, S.W., Yang, T.L., Liou, J.S.",
         "title": "Heat transfer and pressure drop in tube with broken "
                  "twisted tape insert",
         "ref": "Exp. Thermal Fluid Sci. 32(2) (2007) 489-501",
         "doi": "10.1016/j.expthermflusci.2007.06.002"},
    30:
        {"autor": "Eiamsa-ard, S., Thianpong, C., Eiamsa-ard, P.",
         "title": "Turbulent heat transfer enhancement by counter/co-swirling "
                  "flow in a tube fitted with twin twisted tapes",
         "ref": "Exp. Thermal Fluid Sci. 34(1) (2010) 53-62",
         "doi": "10.1016/j.expthermflusci.2009.09.002"},
    31:
        {"autor": "Eiamsa-ard, S., Wongcharee, K., Eiamsa-ard, P., Thianpong, C.",
         "title": "Heat transfer enhancement in a tube using delta-winglet "
                  "twisted tape inserts",
         "ref": "Applied Thermal Eng. 30(4) (2010) 310-318",
         "doi": "10.1016/j.applthermaleng.2009.09.006"},
    32:
        {"autor": "Eiamsa-ard, S., Seemawute, P., Wongcharee, K.",
         "title": "Influences of peripherally-cut twisted tape insert on "
                  "heat transfer and thermal performance characteristics in "
                  "laminar and turbulent tube flows",
         "ref": "Exp. Thermal Fluid Sci. 34(6) (2010) 711-719",
         "doi": "10.1016/j.expthermflusci.2009.12.013"},
    33:
        {"autor": "Eiamsa-ard, P., Piriyarungrod, N., Thianpong, C., "
                  "Eiamsa-ard, S.",
         "title": "A case study on thermal performance assessment of a heat "
                  "exchanger tube equipped with regularly-spaced twisted "
                  "tapes as swirl generators",
         "ref": "Case Studies Thermal Eng. 3 (2014) 86-102",
         "doi": "10.1016/j.csite.2014.04.002"},
    34:
        {"autor": "Eiamsa-ard, S., Wongcharee, K., Eiamsa-ard, P., "
                  "Thianpong, C.",
         "title": "Thermohydraulic investigation of turbulent flow through a "
                  "round tube equipped with twisted tapes consisting of "
                  "centre wings and alternate-axes",
         "ref": "Exp. Thermal Fluid Sci. 34(8) (2010) 1151-1161",
         "doi": "10.1016/j.expthermflusci.2010.04.004"},
    35:
        {"autor": "Eiamsa-ard, S., Thianpong, C., Eiamsa-ard, P., Promvonge, P.",
         "title": "Thermal characteristics in a heat exchanger tube fitted "
                  "with dual twisted tape elements in tandem",
         "ref": "Int. Comm. Heat Mass Transfer 37(1) (2010) 39-46",
         "doi": "10.1016/j.icheatmasstransfer.2009.08.010"},
    36:
        {"autor": "Ponnada, S., Subrahmanyam, T., Naidu, S.V.",
         "title": "A comparative study on the thermal performance of water in "
                  "a circular tube with twisted tapes, perforated twisted "
                  "tapes and perforated twisted tapes with alternate axis",
         "ref": "Int. J. Thermal Sci. 136 (2019) 530-538",
         "doi": "10.1016/j.ijthermalsci.2018.11.008"},
    37:
        {"autor": "He, Y., Liu, L., Li, P., Ma, L.",
         "title": "Experimental study on Heat transfer enhancement "
                  "characteristics of tube with cross hollow twisted tape ",
         "ref": "Applied Thermal Eng. 131 (2018) 743-749",
         "doi": "10.1016/j.applthermaleng.2017.12.029"},
    38:
        {"autor": "Piriyarungrod, N., Eiamsa-ard, S., Thianpong, C., Pimsarn, "
                  "M., Nanan, K.",
         "title": "Heat transfer enhancement by tapered twisted tape inserts",
         "ref": "Chem. Eng. Process. 96 (2015) 62-71",
         "doi": "10.1016/j.cep.2015.08.002"},
    39:
        {"autor": "Eiamsa-ard, S., Somkleang, P., Nuntadusit, C., Thianpong, C.",
         "title": "Heat transfer enhancement in tube by inserting "
                  "uniform/non-uniform twisted-tapes with alternate axes: "
                  "Effect of rotated-axis length",
         "ref": "Applied Thermal Eng. 54 (2013) 289-309",
         "doi": "10.1016/j.applthermaleng.2013.01.041"},
    40:
        {"autor": "Thianpong, C., Eiamsa-ard, S., Somkleang, P.",
         "title": "Heat transfer and thermal performance characteristics of "
                  "heat exchanger tube fitted with perforated twisted-tapes",
         "ref": "Heat Mass Transfer 48(6) (2012) 881-892",
         "doi": "10.1007/s00231-011-0943-0"},

    # 41:
    #     {"autor": "",
    #      "title": "",
    #      "ref": "",
    #      "doi": ""},

        }


# Friction factor correlations
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


@refDoc(__doi__, [27, 28, 29])
def f_twisted_Chang(Re, D, H, mod="", bf=False):
    """Calculate friction factor a pipe with a twisted-tape insert using
    the Chang et al. correlation (2012).

    The twisted-tape have geometrical modifications:

        * PT: Perforated twisted tape
        * PJT: Perforated twisted tape with jaggedness
        * PST: Perforated spiky twisted tape
        * PJST: Perforated spicy twisted tape with jaggedness
        * VST: V-notched spicy twisted tape
        * SR: Serrated roughened
        * BT: Broken twisted tape

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    D : float
        Internal diameter of tube, [m]
    H : float
        Tape pitch for twist of π radians (180º), [m]
    mod : string
        Name of modification code of twisted tape
        PT | PJT | PST | PJSJ | VST
    bf : boolean
        In jaggedness mod flow orientation is relevant, set backward flow state

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    y = H/D

    # Coefficient from Table 5
    if mod == "PT":
        ci = ((0.0174, 0.13, 0.339), (0.02, 4.41, 0.725),
              (3.4e-4, 1.05e-3, 0.512))
    elif mod == "PJT":
        if bf:
            ci = ((0.0174, 0.35, 0.567), (0.02, 0.648, 0.209),
                  (3.4e-4, 3.14e-4, 0.719))
        else:
            ci = ((0.0174, 0.311, 0.493), (0.02, 2.21, 0.565),
                  (3.4e-4, 1.22e-3, 0.699))
    elif mod == "PST":
        ci = ((0.0174, 0.249, 0.34), (0.02, 2.77, 0.527),
              (3.4e-4, 1.06e-3, 0.702))
    elif mod == "PJST":
        if bf:
            ci = ((0.0174, 3.01, 1.21), (0.02, 1.14, 0.295),
                  (3.4e-4, 0.177, 3.04))
        else:
            ci = ((0.0174, 1.89, 0.96), (0.02, 0.364, 0.134),
                  (3.4e-4, 1.09e-3, 0.873))
    elif mod == "VST":
        ci = ((0.0174, 0.561, 0.471), (0.02, 0.706, 0.453),
              (3.4e-4, 2.23, 1.13), (0, 0.00878, 4.64))

    elif mod == "BT":
        # From [29]_, Eq 9
        ci = ((0.0174, 0.21, 0.332), (0.02, 0.161, 1.27),
              (3.4e-4, -3.21e-4, 0.548))

    # From [28]_
    elif mod == "SR":
        # Eq 9 in [28]_
        ci = ((0, 0, 0), (0.033, 0.756, 0.765), (0.166, -0.235, 0.524))
    else:
        # Smooth twisted tape
        # Eq 8 in [28]_
        ci = ((0, 0, 0), (0.07, 9.87, 1.81), (-0.08, -0.94, 1.23))

    c0 = ci[0][0]+ci[0][1]*exp(-ci[0][2]*y)
    c1 = ci[1][0]+ci[1][1]*exp(-ci[1][2]*y)
    c2 = ci[2][0]+ci[2][1]*exp(-ci[2][2]*y)

    # Eq 6
    if not mod or mod == "SR":
        f = c1*Re**c2
    else:
        f = c0 + c1*exp(-c2*Re)

    # Aditional term for VST
    if mod == "VST":
        c3 = ci[3][0]+ci[3][1]*exp(-ci[3][2]*y)
        f += c3*Re

    return f


@refDoc(__doi__, [1])
def f_twisted_laminar_Plessis(Re, D, H, delta, Ae, De):
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
    >>> print("%0.3f" % f_twisted_laminar_Plessis(50, 1, 10, 0, st.Ae, st.De))
    0.849

    >>> print("%0.4f" % f_twisted_laminar_Plessis(2000, 1, 10, 0, st.Ae, st.De))
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
def f_twisted_laminar_Shah(Re, D, H, delta):
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


@refDoc(__doi__, [10])
def f_twisted_laminar_Agarwal(Re, D, H):
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


@refDoc(__doi__, [24])
def f_twisted_laminar_Saha(Re, D, H, delta, S):
    """Calculate friction factor a pipe with a twisted-tape insert using
    the Saha-Gaitonde-Date correlation (1989).

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    D : float
        Internal diameter of tube, [m]
    H : float
        Tape pitch for twist of π radians (180º), [m]
    S : float
        Spacer length without twisted section, [m]

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    y = H/D
    s = S/D

    if Re > 2300 or y > 10 or y < 3 or s < 2.5 or s > 10:
        raise NotImplementedError("Input out of bound")

    Dh0 = (pi*D**2*y+pi*(D**2-delta**2)*s)/((pi+2)*y*D+pi*(D+delta)*s)  # Eq 21
    Dh1 = ((pi*D**2-4*delta*D)*y+pi*(D**2-delta**2)*s) / \
        ((pi+2-2*delta/D)*y*D+pi*(D+delta)*s)                           # Eq 22
    Ac0 = (pi*(D**2*y+(D**2-delta**2)*s))/(4*(y+s))                     # Eq 23
    Ac1 = ((pi*D**2-4*delta*D)*y+pi*(D**2-delta**2)*s)/(4*(y+s))        # Eq 24

    xi = Dh0**2*Ac0/Dh1**2*Ac1                                          # Eq 20

    if s <= 2.5:
        if 7.5 <= y <= 10:
            C = 0.0678*exp(-0.0631*y)*s-0.9936*exp(0.0069*y)+1          # Eq 25
        else:
            C = 0.1998*exp(-0.0631*y)*s + 0.011*y - 0.3175              # Eq 26
    elif s <= 5:
        if 7.5 <= y <= 10:
            C = -0.0031*exp(0.1649*y)*s + 0.02812*exp(0.092*y)          # Eq 27
        else:
            C = -3.97e-3*y*s + 0.01*s + 0.018*y - 7.15e-3               # Eq 30
    elif s <= 7.5:
        if 7.5 <= y <= 10:
            C = -2.45e-5*exp(0.1649*y)*s - 3.51e-3*exp(0.092*y)         # Eq 28
        else:
            C = -4.05e-3*y*s + 0.01*s + 0.018*y - 7.15e-4               # Eq 31
    else:
        if 7.5 <= y <= 10:
            C = -6.39e-4*exp(0.1649*y)*s + 8.7e-3*exp(0.092*y)          # Eq 29
        else:
            C = -2.96e-3*y*s + 0.01*s + 0.018*y - 0.0516                # Eq 32

    C1 = 8.8201*y - 2.1193*y**2 + 0.2108*y**3 - 0.0069*y**4             # Eq 33

    if Re/y <= 100:
        # Eq 17
        f = 38.4 * xi * Re**-0.95 * y**-0.05 * (1+C*s)
    elif Re/y <= 155:
        # Eq 18
        f = 0.5*(38.4*xi*Re**-0.95*y**-0.05*C1*xi*Re**-0.07*y**-0.3)*(1+C*s)
    else:
        # Eq 19
        f = C1*xi*Re**-0.07*y**-0.3*(1+C*s)

    return f


@refDoc(__doi__, [25])
def f_twisted_laminar_Date(Re, D, H):
    """Calculate friction factor a pipe with a twisted-tape insert using
    the Date-Gaitonde correlation (1990).

    The paper include a complex correalation for regularly spaced twisted-tape

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

    # Eq 52
    f = 10.666*(9.818-17.96/y**2)*(y-0.7)/(y-1)*(6e-4/y+0.44*Re)

    return f


@refDoc(__doi__, [7, 8])
def f_twisted_turbulent_Lopina(Re, D, H, Dh):
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
    >>> print("%0.2f" % f_twisted_turbulent_Lopina(24491, 0.61, 6.1, 1.07))
    0.01
    """
    # Using the modified parameter in HTRI Design Manual, section B2.1.2.1

    Reh = Re*Dh/D
    y = H/D

    f = 3.8/y**0.406*(0.046/Reh**0.2)
    return f


@refDoc(__doi__, [9])
def f_twisted_turbulent_Naphon(Re, D, H):
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


@refDoc(__doi__, [12])
def f_twisted_turbulent_Smithberg(Re, D, H):
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
    return 4*f


@refDoc(__doi__, [16, 17, 18, 19, 20])
def f_twisted_turbulent_Murugesan(Re, D, H, mod="", de=None, w=None):
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


@refDoc(__doi__, [21])
def f_twisted_turbulent_Jaisankar(Re, D, H):
    """Calculate friction factor a pipe with a twisted-tape insert using
    the Jaisankar et al. correlation (2009).

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
    if Re < 3000:
        raise NotImplementedError("Input out of bound")

    # Eq 10
    f = 271.1*Re**-0.947*(H/D)**-0.584

    return f


@refDoc(__doi__, [30, 31, 32, 33, 34, 35, 38, 39, 40])
def f_twisted_turbulent_Eiamsaard(Re, D, H, mod="", **kw):
    """Calculate friction factor a pipe with a twisted-tape insert using
    the Eiamsa-ard et al. correlation (2010).

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    D : float
        Internal diameter of tube, [m]
    H : float
        Tape pitch for twist of π radians (180º), [m]
    mod : string
        Name of modification code of twisted tape
        CT | CoT | oDWT | sDWT | PCT | ST | WT | AWT | DST | TT | AT | PT
    dW : float
        depth of wing cut, [m]
    w : float
        Peripherally-cut width
    S : float
        Spacer length without twisted section, [m]
    beta : float
        Attack angle, [ºdeg]
    l : float
        Length of alternate axis, [m]
    sP : float
        Spaced-pitch length of perforated, [m]
    dP : float
        Diameter of perforated, [m]

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    # Be careful with nomenclature in papers, use y as tape width
    y = H/D

    if mod == "CT":
        # Eq 14 from [30]_
        f = 72.29 / Re**0.53 / y**1.01
    elif mod == "CoT":
        # Eq 16 from [30]_
        f = 41.7 / Re**0.52 / y**0.84
    elif mod == "oDWT":
        # Eq 16 from [31]_
        dW = kw.get("dW", 0)
        f = 24.8 / Re**0.51 / y**0.566 * (1+dW/D)**1.87
    elif mod == "sDWT":
        # Eq 19 from [31]_
        dW = kw.get("dW", 0)
        f = 21.7 / Re**0.45 / y**0.564 * (1+dW/D)**1.41
    elif mod == "PCT":
        # Eq 18 from [32]_
        w = kw.get("w", 0)
        f = 39.46 / Re**0.591 * y**0.195 / (w/D)**0.201
    elif mod == "ST":
        # Eq 20 from [33]_
        S = kw.get("S", 0)
        f = 3.044 / y**0.556 / (S/D+1)**0.34
    elif mod == "WT":
        # Table 2 from [34]_
        beta = kw.get("beta", 0)
        f = 14.039 / Re**0.505 * (1+tan(beta))**0.406
    elif mod == "AWT":
        # Table 2 from [34]_
        beta = kw.get("beta", 0)
        f = 20.445 / Re**0.504 * (1+tan(beta))**0.283
    elif mod == "DST":
        # Eq 24 from [35]_
        S = kw.get("S", 0)
        f = 30.5 / Re**0.56 / y**0.54 / (1.5*S/D+1)**0.2
    elif mod == "TT":
        # Eq 20 from [38]_
        teta = kw.get("teta", 0)
        f = 16.559 / Re**0.49 / y**0.51 / (1+teta)**0.53
    elif mod == "AT":
        # Eq 18 from [39]_
        l = kw.get("l", 0)
        f = 75.18 / Re**0.612 / (1+l/H)**0.623
    elif mod == "PT":
        # Eq 13 from [40]_
        sP = kw.get("sP", 0)
        dP = kw.get("dP", 0)
        f = 9.03 / Re**0.272 / y**0.631 / (sP/D)**0.204 * (dP/D)**0.428
    else:
        # Eq 12 from [30]_
        f = 65.4 / Re**0.52 / y**1.31

    return f


@refDoc(__doi__, [36])
def f_twisted_turbulent_Ponnada(Re, D, H, mod=""):
    """Calculate friction factor a pipe with a twisted-tape insert using
    the Ponnada et al. correlation (2019).

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    D : float
        Internal diameter of tube, [m]
    H : float
        Tape pitch for twist of π radians (180º), [m]
    mod : string
        Name of modification code of twisted tape
        PTT | PATT

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    # Be careful with nomenclature in papers, use y as tape width
    y = H/D

    if mod == "PTT":
        # Eq 23
        f = 0.881 / Re**0.359 / y**0.046
    elif mod == "PATT":
        # Eq 25
        f = 1.295 / Re**0.399 / y**0.049
    else:
        # Eq 21
        f = 1.093 / Re**0.39 * y**0.004

    return f


# Heat Transfer coefficient correlations
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
        Nu = Nu_twisted_turbulent_Lopina(Re, Pr, D, H, Dh, mu, muW, beta, dT, HTRI=True)

    else:
        NuL = Nu_twisted_HTRI(2000, Pr, D, H, Dh, mu, muW)
        NuH = Nu_twisted_HTRI(5000, Pr, D, H, Dh, mu, muW)
        Nu = NuL*(5000-Re)/3000 + NuH*(Re-2000)/3000

    return Nu


@refDoc(__doi__, [4, 5])
def Nu_twisted_Manglik(Re, Pr, D, H, delta, Dh, mu, muW):
    r"""Calculate Nusselt number for a pipe with a twisted-tape insert using
    the Manglik and Bergles correlation (1993)

    In laminar flow for simplicity use only fully developed swirl flow
    correlation and without thermal entrance effects and combined forced and
    free convection.

    .. math::
        Nu = 4.612 \left(1+6.413e-9\left(Sw Pr^{0.391}\right)^{3.835}\right)^{0.2}

    For turbulent flow:

    .. math::
        \frac{Nu}{Nu_{y=\infty}} = 1+\frac{0.769}{y}

    .. math::
        Nu_{y=\infty} = 0.023 Re^{0.8} Pr^{0.4}
        \left(\frac{\pi}{pi-4\delta/d}\right)^{0.8}
        \left(\frac{\pi+2-2\delta/d}{\pi-4\delta/d}\right)^{0.2}


    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    Pr : float
        Prandtl number, [-]
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

        Resw = Re*(pi/(pi-4*delta/D)) * ((pi*D/H)**2)**0.5
        Sw = Resw/(H/2/D)**0.5

        # Eq 14
        Nu = 4.612*(mu/muW)**0.14*((1 + 6.413e-9*(Sw*Pr**0.391)**3.835)**0.2)

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
        NuL = Nu_twisted_Manglik(2000, Pr, D, H, delta, Dh, mu, muW)
        NuH = Nu_twisted_Manglik(2000, Pr, D, H, delta, Dh, mu, muW)
        Nu = NuL*(5000-Re)/3000 + NuH*(Re-2000)/3000
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


@refDoc(__doi__, [27, 28, 29])
def Nu_twisted_Chang(Re, Pr, D, H, mod="", bf=False):
    """Calculate friction factor a pipe with a twisted-tape insert using
    the Chang et al. correlation (2012).

    The twisted-tape have geometrical modifications:

        * PT: Perforated twisted tape
        * PJT: Perforated twisted tape with jaggedness
        * PST: Perforated spiky twisted tape
        * PJST: Perforated spicy twisted tape with jaggedness
        * VST: V-notched spicy twisted tape
        * SR: Serrated roughened twisted tape
        * BT: Broken twisted tape

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
    mod : string
        Name of modification code of twisted tape
        PT | PJT | PST | PJSJ | VST
    bf : boolean
        In jaggedness mod flow orientation is relevant, set backward flow state

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    y = H/D

    # Coef for developed flow from Table 3(b)
    if mod == "PT":
        ai = (0.0364, 0.771, 0.363)
        bi = (0.782, 0.239, 0.145)
    elif mod == "PJT":
        if bf:
            ai = (0.0364, 0.581, 0.282)
            bi = (0.782, 0.199, 0.103)
        else:
            ai = (0.0364, 0.732, 0.367)
            bi = (0.782, 0.223, 0.161)
    elif mod == "PST":
        ai = (0.0452, 0.56, 0.286)
        bi = (0.77, 0.184, 0.147)
    elif mod == "PJST":
        if bf:
            ai = (0.0452, 0.416, 0.323)
            bi = (0.77, 0.144, 0.098)
        else:
            ai = (0.0452, 0.333, 0.323)
            bi = (0.77, 0.126, 0.226)
    elif mod == "VST":
        ai = (0.0452, 0.337, 0.204)
        bi = (0.77, 0.122, 0.109)

    elif mod == "BT":
        # From [29]_, Eq 4
        # Using only the parameter for developed flow
        ai = (0.0452, 0.3, 0.141)
        bi = (0.77, 0.157, 0.13)

    # From [28]_
    elif mod == "SR":
        # Eq 7 in [28]_
        ai = (0.118, 5.84, 1.83)
        bi = (0.73, 0.695, 1.26)
    else:
        # Smooth twisted tape
        # Eq 6 in [28]_
        ai = (0.0364, 3.66, 1.11)
        bi = (0.8, 0.375, 0.31)

    A = ai[0]+ai[1]*exp(-ai[2]*y)
    B = bi[0]-bi[1]*exp(-bi[2]*y)

    Nu = A*Re**B*Pr**(1/3)

    return Nu


@refDoc(__doi__, [2])
def Nu_twisted_laminar_Plessis(Re, Pr, D, H, delta, Ae, De, x=None):
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
def Nu_twisted_laminar_Hong(Re, Pr, D, H):
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


@refDoc(__doi__, [10])
def Nu_twisted_laminar_Agarwal(Re, Pr, D, H, mu, muW):
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


@refDoc(__doi__, [24])
def Nu_twisted_laminar_Saha(Re, Pr, D, H, delta, S):
    """Calculate Nusselt number for a pipe with a twisted-tape insert using
    the Saha-Gaitonde-Date correlation (1989).

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
    delta : float
        Tape thickness, [m]
    S : float
        Spacer length without twisted section, [m]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    y = H/D
    s = S/D

    if Re > 2300 or y > 10 or y < 3 or s < 2.5 or s > 10:
        raise NotImplementedError("Input out of bound")

    K1 = pi*D**2*(y+s)/((pi*D**2-4*delta*D)*y+pi*(D**2-delta**2)*s)     # Eq 35

    if y < 7.5:
        C = (0.057*y*s+0.3622)*exp((-0.0296*y-0.305)*s)                 # Eq 38
    elif s <= 5:
        C = 0.0112*y*s - 0.1233*s - 0.0629*y + 0.6948                   # Eq 36
    else:
        C = 0.00015*y*s - 0.00377*s - 0.0056*y + 0.0751                 # Eq 37

    if Re < 700:
        X = 1-4.0422e-2*s                                               # Eq 39
    else:
        X = 1

    Nu = 5.172*(1+6.7482e-3*Pr**0.7*(K1*Re/y)**1.25)**0.5*(1+C*s)*X     # Eq 34

    return Nu


@refDoc(__doi__, [24])
def Nu_twisted_laminar_Klaczak(Re, Pr, D, H, delta, mu=None, muW=None):
    """Calculate Nusselt number for a pipe with a twisted-tape insert using
    the Klaczak correlation (2000).

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
    delta : float
        Tape thickness, [m]
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

    Sw = Re/y**0.5*(pi/(pi-4*delta/D)) * (1+(pi/2/y)**2)**0.5           # Eq 10

    if Sw < 58 or Sw > 2300 or y > 5.29 or y < 1.62 or Pr < 2.06 or Pr > 2.73:
        raise NotImplementedError("Input out of bound")

    Nu = 0.858 * Pr**0.3 * Sw**0.3                                      # Eq 11

    if mu and muW:
        Nu *= (mu/muW)**0.14

    return Nu


@refDoc(__doi__, [7, 8])
def Nu_twisted_turbulent_Lopina(Re, Pr, D, H, Dh, mu, muW, beta, DT, HTRI=False):
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


@refDoc(__doi__, [9])
def Nu_twisted_turbulent_Naphon(Re, Pr, D, H):
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


@refDoc(__doi__, [11])
def Nu_twisted_turbulent_Kidd(Re, Pr, D, H, L, T, Tw):
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


@refDoc(__doi__, [12])
def Nu_twisted_turbulent_Smithberg(Re, Pr, D, H, Dh):
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
    f = f_twisted_turbulent_Smithberg(Re, D, H)

    # Eq 38
    P = (1+0.0219/(H/D)**2/f)**0.5
    A = 50.9*D/H/Re/f**0.5 + 0.023*D/Dh/Re**0.2/Pr**(2/3)*P
    B = 1 + 700/Re/f * D/H * Dh/D * Pr**0.731
    Nu = Re*Pr*A/B

    return Nu


@refDoc(__doi__, [16, 17, 18, 19, 20])
def Nu_twisted_turbulent_Murugesan(Re, Pr, D, H, mod="", de=None, w=None):
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


@refDoc(__doi__, [21])
def Nu_twisted_turbulent_Jaisankar(Re, Pr, D, H):
    """Calculate Nusselt number for a pipe with a twisted-tape insert using
    the Jaisankar et al. correlation (2009).

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
    if Re < 3000:
        raise NotImplementedError("Input out of bound")

    # Eq 9
    Nu = 0.000115*Re**1.169*Pr**2.424*(H/D)**-0.511

    return Nu


@refDoc(__doi__, [30, 31, 32, 33, 34, 35, 38, 39, 40])
def Nu_twisted_turbulent_Eiamsaard(Re, Pr, D, H, mod="", **kw):
    """Calculate nusselt number for a pipe with a twisted-tape insert using
    the Eiamsa-ard et al. correlation (2010).

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
    mod : string
        Name of modification code of twisted tape
        CT | CoT | oDWT | sDWT | PCT | ST | WT | AWT | DST | TT | AT | PT
    dW : float
        depth of wing cut, [m]
    w : float
        Peripherally-cut width
    S : float
        Spacer length without twisted section, [m]
    beta : float
        Attack angle, [ºdeg]
    teta : float
        Taper angle, [ºdeg]
    l : float
        Length of alternate axis, [m]
    sP : float
        Spaced-pitch length of perforated, [m]
    dP : float
        Diameter of perforated, [m]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    # Be careful with nomenclature in papers, use y as tape width
    y = H/D

    if mod == "CT":
        # Eq 13 from [30]_
        Nu = 0.473 * Re**0.66 * Pr**0.4 / y**0.9
    elif mod == "CoT":
        # Eq 15 from [30]_
        Nu = 0.264 * Re**0.66 * Pr**0.4 / y**0.61
    elif mod == "oDWT":
        # Eq 15 from [31]_
        dW = kw.get("dW", 0)
        Nu = 0.18 * Re**0.67 * Pr**0.4 / y**0.423 * (1+dW/D)**0.982
    elif mod == "sDWT":
        # Eq 18 from [31]_
        dW = kw.get("dW", 0)
        Nu = 0.184 * Re**0.675 * Pr**0.4 / y**0.465 * (1+dW/D)**0.76
    elif mod == "PCT":
        # Eq 17 from [32]_
        w = kw.get("w", 0)
        Nu = 0.244 * Re**0.625 * Pr**0.4 * y**0.168 / (w/D)**0.112
    elif mod == "ST":
        # Eq 19 from [33]_
        S = kw.get("S", 0)
        Nu = 0.144 * Re**0.697 * Pr**0.4 / y**0.228 / (S/D+1)**0.179
    elif mod == "WT":
        # Table 2 from [34]_
        beta = kw.get("beta", 0)
        Nu = 0.232 * Re**0.595 * Pr**0.4 * (1+tan(beta))**0.202
    elif mod == "AWT":
        # Table 2 from [34]_
        beta = kw.get("beta", 0)
        Nu = 0.385 * Re**0.568 * Pr**0.4 * (1+tan(beta))**0.129
    elif mod == "DST":
        # Eq 23 from [35]_
        S = kw.get("S", 0)
        Nu = 0.069 * Re**0.74 * Pr**0.4 / y**0.26 / (1.5*S/D+1)**0.1
    elif mod == "TT":
        # Eq 19 from [38]_
        teta = kw.get("teta", 0)
        Nu = 0.076 * Re**0.75 * Pr**0.4 / y**0.39 / (1+teta)**0.1
    elif mod == "AT":
        # Eq 17 from [39]_
        l = kw.get("l", 0)
        Nu = 1.364 * Re**0.472 * Pr**0.4 / (1+l/H)**0.437
    elif mod == "PT":
        # Eq 12 from [40]_
        sP = kw.get("sP", 0)
        dP = kw.get("dP", 0)
        Nu = 0.09 * Re**0.768 * Pr**0.4 / y**0.325 / (sP/D)**0.133 * (dP/D)**0.114
    else:
        # Eq 11 from [30]_
        Nu = 0.224 * Re**0.66 * Pr**0.4 / y**0.6

    return Nu


@refDoc(__doi__, [36])
def Nu_twisted_turbulent_Ponnada(Re, Pr, D, H, mod=""):
    """Calculate nusselt number for a pipe with a twisted-tape insert using
    the Ponnada et al. correlation (2019).

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
    mod : string
        Name of modification code of twisted tape
        PTT | PATT

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    # Be careful with nomenclature in papers, use y as tape width
    y = H/D

    if mod == "PTT":
        # Eq 22
        Nu = 1.199 * Re**0.507 / Pr**0.243 / y**0.132
    elif mod == "PATT":
        # Eq 24
        Nu = 1.161 * Re**0.509 / Pr**0.214 / y**0.141
    else:
        # Eq 20
        Nu = 0.388 * Re**0.541 * Pr**0.212 / y**0.145

    return Nu


# Helical screw-tape
@refDoc(__doi__, [13, 14, 22, 23])
def f_helical_Sivashanmugam(Re, D, H, S=None):
    """Calculate friction factor for a pipe with a twisted-tape insert using
    the Sivashanmugam-Suresh correlation (2006) with lamanar flow correlation
    for spacer from Ibrahim correlation (2011).

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    D : float
        Internal diameter of tube, [m]
    H : float
        Tape pitch for twist of π radians (180º), [m]
    S : float
        Spacer length without twisted section, [m]

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    if Re > 2000:
        # Turbulent flow
        if S:
            # Eq 6 in 22_
            f = Re**-0.384*(H/D)**-0.852*(1+S/D)**-0.047
        else:
            # Eq 6 in 14_
            f = 32.415*Re**-0.598*(H/D)**-0.7986

    else:
        # Laminar flow
        if S:
            # Eq 12 in 23_
            Nu = 54.41*Re**-0.87*(1+S/D)**-0.045*(H/D)**-0.146
        else:
            # Eq 6 in 13_
            f = 10.7564*Re**0.387*(H/D)**-1.054

    return f


@refDoc(__doi__, [13, 14, 22, 23])
def Nu_helical_Sivashanmugam(Re, Pr, D, H, S=None):
    """Calculate Nusselt number for a pipe with a twisted-tape insert using
    the Sivashanmugam-Suresh correlation (2006) with lamanar flow correlation
    for spacer from Ibrahim correlation (2011).

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
    S : float
        Spacer length without twisted section, [m]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    if Re > 2000:
        # Turbulent flow
        if S:
            # Eq 5 in 22_
            Nu = 0.258*Re**0.554*Pr*(H/D)**-0.242*(1+S/D)**-0.042
        else:
            # Eq 5 in 14_
            Nu = 0.4675*Re**0.4774*Pr*(H/D)**-0.2138

    else:
        # Laminar flow
        if S:
            # Eq 11 in 23_
            Nu = 6.11*Re**0.199*(1+S/D)**-0.064*(H/D)**-0.318
        else:
            # Eq 6 in 13_
            Nu = 0.017*Re**0.996*Pr*(H/D)**-0.5437

    return Nu


# Hollow twisted-tape
@refDoc(__doi__, [37])
def f_hollow_He(Re, D, C):
    """Calculate friction factor for a pipe with a hollow twisted-tape insert
    using the He et al. correlation (2018)

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    D : float
        Internal diameter of tube, [m]
    C : float
        Hollow width of the cross hollow twisted tape, [m]

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    c = C/D

    # Eq 17
    f = 9.348 / Re**0.3959 * (5.53*c**3 + 2.578*c**2 - 7.307*c + 3.499)

    return f


@refDoc(__doi__, [37])
def Nu_hollow_He(Re, Pr, D, C):
    """Calculate nusselt number for a pipe with a hollow twisted-tape insert
    using the He et al. correlation (2018)

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    Pr : float
        Prandtl number, [-]
    D : float
        Internal diameter of tube, [m]
    C : float
        Hollow width of the cross hollow twisted tape, [m]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    c = C/D

    # Eq 16
    Nu = 0.3415 * Re**0.5911 * Pr**0.32 * \
        (0.9058*c**3 + 0.5439*c**2 - 1.345*c + 1.271)

    return Nu



class TwistedTape(CallableEntity):
    """Twisted-tape insert used in heat exchanger to improve efficiency.
    This tape, generally a thin metal strip, is twisted about its longitudinal
    axis

    Parameters
    ----------
    methodFrictionLaminar : integer
        Index of method used for friction factor calculation in laminar flow
    methodFTurbulent: integer
        Index of method used for friction factor calculation in turbulent flow
    methodHeatLaminar : integer
        Index of method used for heat transfer calculation in laminar flow
    methodHeatTurbulent : integer
        Index of method used for heat transfer calculation in turbulent flow
    H : float
        Tape pitch for twist of π radians (180º), [m]
    Dt : float
        Internal diameter of tube, [m]
    delta : float
        Tape thickness, [m]

    """

    TEXT_LAMINAR_FRICTION = (
        "Manglik-Bergles (1993)",
        "Plessis-Kröger (1984)",
        "Shah-London (1978)",
        "Agarwal-Rao (1996)",
        "Sarma (2005)",
        "Saha-Gaitonde-Date (1989)",
        "Date-Gaitonde (1990)",
        "Chang (2012)"
        )

    TEXT_TURBULENT_FRICTION = (
        "Manglik-Bergles (1993)",
        "Lopina-Bergles (1969)",
        "Naphon (2006)",
        "Sarma (2005)",
        "Smithberg-Landis (1964)",
        "Murugesan (2010)",
        "Jaisankar (2009)",
        "Chang (2012)",
        "Eiamsa-ard (2010)",
        "Ponnada (2019)"
        )

    TEXT_LAMINAR_HEAT = (
        "HTRI",
        "Manglik-Bergles (1993)",
        "Plessis-Kröger (1984)",
        "Hong-Bergles (1976)",
        "Agarwal-Rao (1996)",
        "Sarma (2005)",
        "Saha-Gaitonde-Date (1989)",
        "Klaczak (2000)",
        "Chang (2012)"
        )

    TEXT_TURBULENT_HEAT = (
        "HTRI",
        "Manglik-Bergles (1993)",
        "Lopina-Bergles (1969)",
        "Naphon (2006)",
        "Kidd (1969)",
        "Sarma (2005)",
        "Smithberg-Landis (1964)",
        "Murugesan (2010)",
        "Jaisankar (2009)",
        "Chang (2012)",
        "Eiamsa-ard (2010)",
        "Ponnada (2019)"
        )

    TEXT_MURUGESAN = (
        "",
        "Nails",
        "Square cut",
        "V cut",
        "Trapezoidal cut",
        "Vertical wings",
        "Horizontal wings")

    TEXT_CHANG = ("", "PT", "PJT", "PST", "PJST", "VST", "SR", "BT")
    TEXT_CHANG_TOOLTIP = (
        "",
        translate("twistedtape", "Perforated twisted tape"),
        translate("twistedtape", "Perforated twisted tape with jaggedness"),
        translate("twistedtape", "Perforated spiky twisted tape"),
        translate("twistedtape", "Perforated spicy twisted tape with jaggedness"),
        translate("twistedtape", "V-notched spicy twisted tape"),
        translate("twistedtape", "Serrated roughened twisted tape"),
        translate("twistedtape", "Broken twisted tape"))

    TEXT_EIAMSA = ("", "CT", "CoT", "oDWT", "sDWT", "PCT","WT", "AWT", "ST",
                   "DST", "TT", "AT", "PT")
    TEXT_EIAMSA_TOOLTIP = (
        "",
        translate("twistedtape", "Twin counter twisted tape"),
        translate("twistedtape", "Twin co-twisted tape"),
        translate("twistedtape", "Oblique delta-winglet twisted tape"),
        translate("twistedtape", "Straight delta-winglet twisted tape"),
        translate("twistedtape", "Peripherally-cut twisted tape"),
        translate("twistedtape", "Twisted tape with centre wings"),
        translate("twistedtape", "Twisted tape with centre wings and alternate axes"),
        translate("twistedtape", "Regularly-spaced twisted tape"),
        translate("twistedtape", "Regularly-spaced dual twisted tape"),
        translate("twistedtape", "Tapered twisted tape"),
        translate("twistedtape", "Alternate axes twisted tape"),
        translate("twistedtape", "Perforated twisted tape"))

    TEXT_PONNADA = ("", "PTT", "PATT")
    TEXT_PONNADA_TOOLTIP = (
        "",
        translate("twistedtape", "Perforated twisted tape"),
        translate("twistedtape", "Perforated twisted tape with alternate axis"))

    # Helical method
    # "Sivashanmugam-Suresh-Ibrahim (2006)"

    status = 0
    msg = ""
    kw = {
        "methodFrictionLaminar": 0,
        "methodFTurbulent": 0,
        "methodHeatLaminar": 0,
        "methodHeatTurbulent": 0,

        "H": 0,
        "Dt": 0,
        "delta": 0,
        "S": 0,

        "isHelical": False,
        "isHollow": False,
        "C": 0,

        "modMurugesan": "",
        "Vcut_w": 0,
        "Vcut_De": 0,

        "modChang": "",
        "bf": False,

        "modEiamsa": "",
        "dW": 0,
        "w": 0,
        "beta": 0,
        "teta": 0,
        "l": 0,
        "sP": 0,
        "dP": 0,

        "modPonnada": ""}

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

        if self.kw["isHollow"] and not self.kw["C"]:
            self.msg = translate("equipment", "undefined hollow width")
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

        # Helix angle
        self.alpha = atan(pi/2/self.y)

        self.valueChanged.emit(self)

    def Nu(self, Re, Pr, mu, muW, beta, dT, L, method=None):
        """Calculate nusselt number"""
        msg = ""

        if self.kw["isHelical"]:
            Nu = Nu_helical_Sivashanmugam(Re, Pr, self.Dt, self.H, self.kw["S"])
            return Nu

        if self.kw["isHollow"]:
            Nu = Nu_hollow_He(Re, Pr, self.Dt, self.kw["C"])
            return Nu

        if Re < 2000:
            # Laminar methods

            if method is None:
                method = self.kw["methodHeatLaminar"]

            if method == 1:
                # Manglik-Bergles (1993)
                Nu = Nu_twisted_Manglik(Re, Pr, self.Dt, self.H,
                                        self.delta, self.De, mu, muW)

            elif method == 2:
                # Plessis-Kröger (1987)
                Nu = Nu_twisted_laminar_Plessis(
                    Re, Pr, self.Dt, self.H, self.delta, self.Ae, self.De)

            elif method == 3:
                # Hong-Bergles (1976)
                Nu = Nu_twisted_laminar_Hong(Re, Pr, self.Dt, self.H)

            elif method == 4:
                # Agarwal-Rao (1996)
                try:
                    Nu = Nu_twisted_laminar_Agarwal(
                        Re, Pr, self.Dt, self.H, mu, muW)
                except NotImplementedError:
                    Nu = self.Nu(Re, Pr, mu, muW, beta, dT, L, method=0)
                    msg = "Agarwal correlation out of range, using HTRI instead"

            elif method == 5:
                # Sarma (2005)
                Nu = Nu_twisted_Sarma(Re, Pr, self.Dt, self.H)

            elif method == 6:
                # Saha-Gaitonde-Date (1989)
                try:
                    Nu = Nu_twisted_laminar_Saha(
                        Re, Pr, self.Dt, self.H, self.delta, self.kw["S"])
                except NotImplementedError:
                    Nu = self.Nu(Re, Pr, mu, muW, beta, dT, L, method=0)
                    msg = "Saha-Gaitonde-Date correlation out of range, "
                    msg += "using HTRI instead"

            elif method == 7:
                # Klaczak (2000)
                try:
                    Nu = Nu_twisted_laminar_Klaczak(
                        Re, Pr, self.Dt, self.H, self.delta, mu, muW)
                except NotImplementedError:
                    Nu = self.Nu(Re, Pr, mu, muW, beta, dT, L, method=0)
                    msg = "Klaczak correlation out of range, using HTRI instead"

            elif method == 8:
                # Chang (2012)
                Nu = Nu_twisted_Chang(Re, Pr, self.Dt, self.H,
                                      self.kw["modChang"], self.kw["bf"])

            else:
                # HTRI
                Nu = Nu_twisted_HTRI(
                    Re, Pr, self.Dt, self.H, self.De, mu, muW, beta, dT, L)

        else:
            # Turbulent methods
            if method is None:
                method = self.kw["methodHeatTurbulent"]

            if method == 1:
                # Manglik-Bergles (1993)
                Nu = Nu_twisted_Manglik(Re, Pr, self.Dt, self.H,
                                        self.delta, self.De, mu, muW)

            elif method == 2:
                # Lopina-Bergles (1969)
                Nu = Nu_twisted_turbulent_Lopina(
                    Re, Pr, self.Dt, self.H, self.De, mu, muW, beta, dT)

            elif method == 3:
                # Naphon (2006)
                Nu = Nu_twisted_turbulent_Naphon(Re, Pr, self.Dt, self.H)

            elif method == 4:
                # Kidd (1969)
                Nu = Nu_twisted_turbulent_Kidd(
                    Re, Pr, self.Dt, self.H, L, 1, 1)

            elif method == 5:
                # Sarma (2005)
                Nu = Nu_twisted_Sarma(Re, Pr, self.Dt, self.H)

            elif method == 6:
                # Smithberg-Landis (1964)
                Nu = Nu_twisted_turbulent_Smithberg(
                    Re, Pr, self.Dt, self.H, self.De)

            elif method == 7:
                # Murugesan (2010)
                if self.kw["modMurugesan"] == "V cut" and self.kw["Vcut_De"] \
                        and self.kw["Vcut_w"]:
                    Nu = Nu_twisted_turbuletn_Murugesan(
                        Re, Pr, self.Dt, self.H, self.kw["modMurugesan"],
                        self.kw["Vcut_De"], self.kw["Vcut_w"])
                elif self.kw["modMurugesan"] == "V cut":
                    Nu = Nu_twisted_turbuletn_Murugesan(
                        Re, Pr, self.Dt, self.H)
                    msg = "V cut twisted tape geometry don't defined, using "
                    msg += "plain twisted tape instead"
                else:
                    Nu = Nu_twisted_turbulent_Murugesan(
                        Re, Pr, self.Dt, self.H, self.kw["modMurugesan"])

            elif method == 8:
                # Jaisankar
                Nu = Nu_twisted_turbulent_Jaisankar(Re, Pr, self.Dt, self.H)

            elif method == 9:
                # Chang (2012)
                Nu = Nu_twisted_Chang(Re, Pr, self.Dt, self.H,
                                      self.kw["modChang"], self.kw["bf"])

            elif method == 10:
                # Eiamsa-ard (2010)
                kw = {"mod": self.kw["modEiamsa"],
                      "dW": self.kw["dW"],
                      "w": self.kw["w"],
                      "S": self.kw["S"],
                      "beta": self.kw["beta"],
                      "l": self.kw["l"],
                      "teta": self.kw["teta"],
                      "sP": self.kw["sP"],
                      "dP": self.kw["dP"]}

                if "DWT" in self.kw["modEiamsa"] and not self.kw["dW"]:
                    kw["mod"] = ""
                    msg = "Depth of wing cut don't defined, using plain "
                    msg += "twisted tape instead"
                elif self.kw["modEiamsa"] == "PCT" and not self.kw["w"]:
                    kw["mod"] = ""
                    msg = "Peripherally-cut width don't defined, using plain "
                    msg += "twisted tape instead"
                elif self.kw["modEiamsa"] == "PT" and \
                        (not self.kw["sP"] or not self.kw["dP"]):
                    kw["mod"] = ""
                    msg = "Perforated geometry don't defined, using plain "
                    msg += "twisted tape instead"

                Nu = Nu_twisted_turbulent_Eiamsaard(
                    Re, Pr, self.Dt, self.H, **kw)

            elif method == 11:
                # Ponnada (2019)
                Nu = Nu_twisted_turbulent_Ponnada(
                    Re, Pr, self.Dt, self.H, self.kw["modPonnada"])

            else:
                # HTRI
                Nu = Nu_twisted_HTRI(
                    Re, Pr, self.Dt, self.H, self.De, mu, muW, beta, dT, L)

        if msg:
            self.status = 3
            self.msg = translate("equipment", msg)
            self.inputChanged.emit(self)

        return Nu

    def f(self, Re, method=None):
        """Calculate friction factor"""
        msg = ""

        if self.kw["isHelical"]:
            f = f_helical_Sivashanmugam(Re, self.Dt, self.H, self.kw["S"])
            return f

        if self.kw["isHollow"]:
            f = f_hollow_He(Re, self.Dt, self.kw["C"])
            return f

        if Re < 2000:
            # Laminar methods

            if method is None:
                method = self.kw["methodFrictionLaminar"]

            if method == 1:
                # Plessis-Kröger (1984)
                f = f_twisted_laminar_Plessis(
                    Re, self.Dt, self.H, self.delta, self.Ae, self.De)

            elif method == 2:
                # Shah-London (1978)
                f = f_twisted_laminar_Shah(Re, self.Dt, self.H, self.delta)

            elif method == 3:
                # Agarwal-Rao (1996)
                try:
                    f = f_twisted_laminar_Agarwal(Re, self.Dt, self.H)
                except NotImplementedError:
                    f = self.f(Re, 0)
                    msg = "Agarwal correlation out of range, "
                    msg += "using Manglik instead"

            elif method == 4:
                # Sarma (2005)
                f = f_twisted_Sarma(Re, self.Dt, self.H)

            elif method == 5:
                # Saha-Gaitonde-Date (1989)
                try:
                    f = f_twisted_laminar_Saha(
                        Re, self.Dt, self.H, self.delta, self.kw["S"])
                except NotImplementedError:
                    f = self.f(Re, 0)
                    msg = "Saha-Gaitonde-Date correlation out of range, "
                    msg += "using Manglik instead"

            elif method == 6:
                # Date-Gaitonde (1990)
                f = f_twisted_laminar_Date(Re, self.Dt, self.H)

            elif method == 7:
                # Chang (2012)
                f = f_twisted_Chang(
                    Re, self.Dt, self.H, self.kw["modChang"], self.kw["bf"])

            else:
                # Manglik-Bergles (1993)
                f = f_twisted_Manglik(Re, self.Dt, self.H, self.delta, self.De)

        else:
            # Turbulent methods
            if method is None:
                method = self.kw["methodFTurbulent"]

            if method == 1:
                # Lopina-Bergles (1969)
                f = f_twisted_turbulent_Lopina(Re, self.Dt, self.H, self.De)

            elif method == 2:
                # Naphon (2006)
                f = f_twisted_turbulent_Naphon(Re, self.Dt, self.H)

            elif method == 3:
                # Sarma (2005)
                f = f_twisted_Sarma(Re, self.Dt, self.H)

            elif method == 4:
                # Smithberg-Landis (1964)
                f = f_twisted_turbulent_Smithberg(Re, self.Dt, self.H)

            elif method == 5:
                # Murugesan (2010)
                if self.kw["modMurugesan"] == "V cut" and self.kw["Vcut_De"] \
                        and self.kw["Vcut_w"]:
                    f = f_twisted_turbulent_Murugesan(
                        Re, self.Dt, self.H, self.kw["modMurugesan"],
                        self.kw["Vcut_w"], self.kw["Vcut_De"])
                elif self.kw["modMurugesan"] == "V cut":
                    f = f_twisted_turbulent_Murugesan(Re, self.Dt, self.H)
                    msg = "V cut twisted tape geometry don't defined, using "
                    msg += "plain twisted tape instead"
                else:
                    f = f_twisted_turbulent_Murugesan(
                        Re, self.Dt, self.H, self.kw["modMurugesan"])

            elif method == 6:
                # Jaisankar (2009)
                f = f_twisted_turbulent_Jaisankar(Re, self.Dt, self.H)

            elif method == 7:
                # Chang (2012)
                f = f_twisted_Chang(
                    Re, self.Dt, self.H, self.kw["modChang"], self.kw["bf"])

            elif method == 8:
                # Eiamsa-ard (2010)
                kw = {"mod": self.kw["modEiamsa"],
                      "dW": self.kw["dW"],
                      "w": self.kw["w"],
                      "S": self.kw["S"],
                      "beta": self.kw["beta"],
                      "l": self.kw["l"],
                      "teta": self.kw["teta"],
                      "sP": self.kw["sP"],
                      "dP": self.kw["dP"]}

                if "DWT" in self.kw["modEiamsa"] and not self.kw["dW"]:
                    kw["mod"] = ""
                    msg = "Depth of wing cut don't defined, using plain "
                    msg += "twisted tape instead"
                elif self.kw["modEiamsa"] == "PCT" and not self.kw["w"]:
                    kw["mod"] = ""
                    msg = "Peripherally-cut width don't defined, using plain "
                    msg += "twisted tape instead"
                elif self.kw["modEiamsa"] == "PT" and \
                        (not self.kw["sP"] or not self.kw["dP"]):
                    kw["mod"] = ""
                    msg = "Perforated geometry don't defined, using plain "
                    msg += "twisted tape instead"

                f = f_twisted_turbulent_Eiamsaard(Re, self.Dt, self.H, **kw)

            elif method == 9:
                # Ponnada (2019)
                f = f_twisted_turbulent_Ponnada(
                    Re, self.Dt, self.H, self.kw["modPonnada"])

            else:
                # Manglik-Bergles (1993)
                f = f_twisted_Manglik(Re, self.Dt, self.H, self.delta, self.De)

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

        lyt = self.wdg.layout()


        self.twisted = QtWidgets.QRadioButton(self.tr(
            "Tipical twisted-tape inserts"))
        self.twisted.setChecked(True)
        self.twisted.toggled.connect(self.setVisibleMod)
        lyt.addWidget(self.twisted, 1, 1, 1, 2)

        self.groupMethods = QtWidgets.QWidget()
        lytM = QtWidgets.QGridLayout(self.groupMethods)
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
        for method in TwistedTape.TEXT_LAMINAR_FRICTION:
            self.methodFrictionLaminar.addItem(method)
        self.methodFrictionLaminar.currentIndexChanged.connect(
            partial(self.changeParams, "methodFrictionLaminar"))
        self.methodFrictionLaminar.currentTextChanged.connect(
            self.setVisibleMod)
        lytM.addWidget(self.methodFrictionLaminar, 2, 2)
        self.methodFTurbulent = QtWidgets.QComboBox()
        for method in TwistedTape.TEXT_TURBULENT_FRICTION:
            self.methodFTurbulent.addItem(method)
        self.methodFTurbulent.currentIndexChanged.connect(
            partial(self.changeParams, "methodFTurbulent"))
        self.methodFTurbulent.currentTextChanged.connect(
            self.setVisibleMod)
        lytM.addWidget(self.methodFTurbulent, 2, 3)

        lytM.addWidget(QtWidgets.QLabel(
            self.tr("Heat transfer method")), 3, 1)
        self.methodHeatLaminar = QtWidgets.QComboBox()
        for method in TwistedTape.TEXT_LAMINAR_HEAT:
            self.methodHeatLaminar.addItem(method)
        self.methodHeatLaminar.currentIndexChanged.connect(
            partial(self.changeParams, "methodHeatLaminar"))
        self.methodHeatLaminar.currentTextChanged.connect(self.setVisibleMod)
        lytM.addWidget(self.methodHeatLaminar, 3, 2)
        self.methodHeatTurbulent = QtWidgets.QComboBox()
        for method in TwistedTape.TEXT_TURBULENT_HEAT:
            self.methodHeatTurbulent.addItem(method)
        self.methodHeatTurbulent.currentIndexChanged.connect(
            partial(self.changeParams, "methodHeatTurbulent"))
        self.methodHeatTurbulent.currentTextChanged.connect(self.setVisibleMod)
        lytM.addWidget(self.methodHeatTurbulent, 3, 3)
        lytM.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 4, 1)
        lyt.addWidget(self.groupMethods, 2, 1, 1, 2)
        self.twisted.toggled.connect(self.groupMethods.setEnabled)

        self.helical = QtWidgets.QRadioButton(self.tr(
            "Helical screw-tape inserts"))
        self.helical.toggled.connect(partial(self.changeParams, "isHelical"))
        self.helical.toggled.connect(self.setEnableSpacer)
        lyt.addWidget(self.helical, 3, 1, 1, 3)

        self.hollow = QtWidgets.QRadioButton(self.tr(
            "Hollow twisted-tape inserts"))
        self.hollow.toggled.connect(self.setEnableHollow)
        lyt.addWidget(self.hollow, 4, 1, 1, 3)
        lyt.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 5, 1, 1, 2)

        label = QtWidgets.QLabel(self.tr("Tape pitch"))
        label.setToolTip(self.tr("Tape pitch for twist of π radians (180º)"))
        lyt.addWidget(label, 6, 1)
        self.H = Entrada_con_unidades(Length)
        self.H.valueChanged.connect(partial(self.changeParams, "H"))
        lyt.addWidget(self.H, 6, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Tape diameter")), 7, 1)
        self.Dt = Entrada_con_unidades(Length, "PipeDiameter")
        self.Dt.valueChanged.connect(partial(self.changeParams, "Dt"))
        lyt.addWidget(self.Dt, 7, 2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Tape thickness")), 8, 1)
        self.delta = Entrada_con_unidades(Length, "Thickness")
        self.delta.valueChanged.connect(partial(self.changeParams, "delta"))
        lyt.addWidget(self.delta, 8, 2)

        self.lblS = QtWidgets.QLabel(self.tr("Spacer length"))
        lyt.addWidget(self.lblS, 9, 1)
        self.S = Entrada_con_unidades(Length)
        self.S.valueChanged.connect(partial(self.changeParams, "S"))
        lyt.addWidget(self.S, 9, 2)

        self.lblC = QtWidgets.QLabel(self.tr("Hollow width"))
        lyt.addWidget(self.lblC, 10, 1)
        self.C = Entrada_con_unidades(Length)
        self.C.valueChanged.connect(partial(self.changeParams, "C"))
        lyt.addWidget(self.C, 10, 2)

        # Murugesan additional parameters
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
        self.lblVcut_w = QtWidgets.QLabel(self.tr("Widgth of V cut"))
        lytMuru.addWidget(self.lblVcut_w, 3, 2)
        self.Vcut_w = Entrada_con_unidades(Length, "thickness")
        self.Vcut_w.valueChanged.connect(partial(self.changeParams, "Vcut_w"))
        lytMuru.addWidget(self.Vcut_w, 3, 3)
        lytMuru.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed), 4, 4)
        lyt.addWidget(self.groupMurugesan, 11, 1, 1, 2)

        # Chang-Guo additional parameters
        self.groupChang = QtWidgets.QWidget()
        lytChang = QtWidgets.QGridLayout(self.groupChang)
        lytChang.addWidget(QtWidgets.QLabel(self.tr(
            "Chan-Guo correlation modification")), 1, 1, 1, 2)
        self.modChang = QtWidgets.QComboBox()
        for method, txt in zip(TwistedTape.TEXT_CHANG, TwistedTape.TEXT_CHANG_TOOLTIP):
            if method and txt:
                self.modChang.addItem(f"{method} - {txt}")
            else:
                self.modChang.addItem("")
        self.modChang.currentTextChanged.connect(self.changeModChang)
        lytChang.addWidget(self.modChang, 1, 3)
        lytChang.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 2, 1)
        self.checkBF = QtWidgets.QCheckBox(self.tr("Backward flow"))
        self.checkBF.toggled.connect(partial(self.changeParams, "bf"))
        lytChang.addWidget(self.checkBF, 2, 2)
        lytChang.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed), 3, 4)
        lyt.addWidget(self.groupChang, 12, 1, 1, 2)

        # Eiamsa-ard additional parameters
        self.groupEiamsa = QtWidgets.QWidget()
        lytEiamsa = QtWidgets.QGridLayout(self.groupEiamsa)
        lytEiamsa.addWidget(QtWidgets.QLabel(self.tr(
            "Eiamsa-ard correlation modification")), 1, 1, 1, 2)
        self.modEiamsa = QtWidgets.QComboBox()
        for method, txt in zip(TwistedTape.TEXT_EIAMSA, TwistedTape.TEXT_EIAMSA_TOOLTIP):
            if method and txt:
                self.modEiamsa.addItem(f"{method} - {txt}")
            else:
                self.modEiamsa.addItem("")
        self.modEiamsa.currentTextChanged.connect(self.changeModEiamsa)
        lytEiamsa.addWidget(self.modEiamsa, 1, 3)

        self.lbldW = QtWidgets.QLabel(self.tr("Depth of wing cut"))
        lytEiamsa.addWidget(self.lbldW, 2, 2)
        self.dW = Entrada_con_unidades(Length, "thickness")
        self.dW.valueChanged.connect(partial(self.changeParams, "dW"))
        lytEiamsa.addWidget(self.dW, 2, 3)
        self.lblw = QtWidgets.QLabel(self.tr("Peripherally-cut width"))
        lytEiamsa.addWidget(self.lblw, 3, 2)
        self.w = Entrada_con_unidades(Length, "thickness")
        self.w.valueChanged.connect(partial(self.changeParams, "w"))
        lytEiamsa.addWidget(self.w, 3, 3)
        self.lblbeta = QtWidgets.QLabel(self.tr("Attack angle"))
        lytEiamsa.addWidget(self.lblbeta, 4, 2)
        self.beta = Entrada_con_unidades(Angle)
        self.beta.valueChanged.connect(partial(self.changeParams, "beta"))
        lytEiamsa.addWidget(self.beta, 4, 3)
        self.lblteta = QtWidgets.QLabel(self.tr("Taper angle"))
        lytEiamsa.addWidget(self.lblteta, 5, 2)
        self.teta = Entrada_con_unidades(Angle)
        self.teta.valueChanged.connect(partial(self.changeParams, "teta"))
        lytEiamsa.addWidget(self.teta, 5, 3)
        self.lbll = QtWidgets.QLabel(self.tr("Alternate axes length"))
        lytEiamsa.addWidget(self.lbll, 6, 2)
        self.l = Entrada_con_unidades(Length)
        self.l.valueChanged.connect(partial(self.changeParams, "l"))
        lytEiamsa.addWidget(self.l, 6, 3)
        self.lbldP = QtWidgets.QLabel(self.tr("Diameter of perforated"))
        lytEiamsa.addWidget(self.lbldP, 7, 2)
        self.dP = Entrada_con_unidades(Length, "thickness")
        self.dP.valueChanged.connect(partial(self.changeParams, "dP"))
        lytEiamsa.addWidget(self.dP, 7, 3)
        self.lblsP = QtWidgets.QLabel(self.tr("Spaced-pitch length of perforated"))
        lytEiamsa.addWidget(self.lblsP, 8, 2)
        self.sP = Entrada_con_unidades(Length, "thickness")
        self.sP.valueChanged.connect(partial(self.changeParams, "sP"))
        lytEiamsa.addWidget(self.sP, 8, 3)


        lyt.addWidget(self.groupEiamsa, 13, 1, 1, 2)

        # Ponnada additional parameters
        self.groupPonnada = QtWidgets.QWidget()
        lytPonnada = QtWidgets.QGridLayout(self.groupPonnada)
        lytPonnada.addWidget(QtWidgets.QLabel(self.tr(
            "Ponnada correlation modification")), 1, 1, 1, 2)
        self.modPonnada = QtWidgets.QComboBox()
        for method, txt in zip(TwistedTape.TEXT_PONNADA, TwistedTape.TEXT_PONNADA_TOOLTIP):
            if method and txt:
                self.modPonnada.addItem(f"{method} - {txt}")
            else:
                self.modPonnada.addItem("")
        self.modPonnada.currentTextChanged.connect(self.changeModPonnada)
        lytPonnada.addWidget(self.modPonnada, 1, 3)
        lytPonnada.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed), 3, 4)
        lyt.addWidget(self.groupPonnada, 14, 1, 1, 2)

        self.Entity.valueChanged.connect(self.valueChanged.emit)
        self.Entity.inputChanged.connect(self.populate)
        self.setVisibleMod()
        self.lblC.setVisible(False)
        self.C.setVisible(False)

    def changeModChang(self, txt):
        """Extract code from txt"""
        if txt:
            txt = txt.split(" - ")[0]
        self.changeParams("modChang", txt)
        self.checkBF.setEnabled("J" in txt)

    def changeModEiamsa(self, txt):
        """Extract code from txt"""
        self.setEnable_Eiamsa(txt)
        if txt:
            txt = txt.split(" - ")[0]
        self.changeParams("modEiamsa", txt)

    def changeModPonnada(self, txt):
        """Extract code from txt"""
        if txt:
            txt = txt.split(" - ")[0]
        self.changeParams("modPonnada", txt)

    def setVisibleMod(self):
        """Enable widget with special parameters for selected method"""
        twisted = self.twisted.isChecked()
        # Murugesan method
        if twisted and (self.methodHeatTurbulent.currentText() == "Murugesan (2010)" or
                        self.methodFTurbulent.currentText() == "Murugesan (2010)"):
            self.groupMurugesan.setVisible(True)
        else:
            self.groupMurugesan.setVisible(False)
        self.setEnable_Murugesan(twisted and self.modMurugesan.currentText())

        # Saha use Spacer special parameter
        self.setEnableSpacer()

        # Chang method
        if twisted and (self.methodHeatTurbulent.currentText() == "Chang (2012)" or
                        self.methodFTurbulent.currentText() == "Chang (2012)" or
                        self.methodHeatLaminar.currentText() == "Chang (2012)" or
                        self.methodFrictionLaminar.currentText() == "Chang (2012)"):
            self.groupChang.setVisible(True)
        else:
            self.groupChang.setVisible(False)
        self.checkBF.setEnabled(twisted and "J" in self.modChang.currentText())

        # Eiamsa-ard method
        if twisted and (self.methodHeatTurbulent.currentText() == "Eiamsa-ard (2010)" or
                        self.methodFTurbulent.currentText() == "Eiamsa-ard (2010)"):
            self.groupEiamsa.setVisible(True)
        else:
            self.groupEiamsa.setVisible(False)
        self.setEnable_Eiamsa(self.modEiamsa.currentText())

        # Ponnada method
        if twisted and (self.methodHeatTurbulent.currentText() == "Ponnada (2019)" or
                        self.methodFTurbulent.currentText() == "Ponnada (2019)"):
            self.groupPonnada.setVisible(True)
        else:
            self.groupPonnada.setVisible(False)

    def setEnable_Murugesan(self, mod):
        """Change Enable/Disable state for Murugesan aditional parameters"""
        self.lblDe.setVisible(mod == "V cut")
        self.De.setVisible(mod == "V cut")
        self.lblVcut_w.setVisible(mod == "V cut")
        self.Vcut_w.setVisible(mod == "V cut")

    def setEnable_Eiamsa(self, mod):
        """Change Enable/Disable state for Eiamsa-ard aditional parameters"""
        self.lbldW.setVisible("DWT" in mod)
        self.dW.setVisible("DWT" in mod)
        self.lblw.setVisible("PCT" in mod)
        self.w.setVisible("PCT" in mod)
        self.lblbeta.setVisible(mod[:3] in ("AWT", "WT "))
        self.beta.setVisible(mod[:3] in ("AWT", "WT "))
        self.lblteta.setVisible("TT" in mod)
        self.teta.setVisible("TT" in mod)
        self.lbll.setVisible("AT" in mod)
        self.l.setVisible("AT" in mod)
        self.lblsP.setVisible("PT" in mod)
        self.sP.setVisible("PT" in mod)
        self.lbldP.setVisible("PT" in mod)
        self.dP.setVisible("PT" in mod)
        self.setEnableSpacer()

    def setEnabled(self, boolean):
        """Add logic to parent setEnabled for orientation option"""
        ToolGui.setEnabled(self, boolean)
        self.setEnableSpacer()

    def setEnableSpacer(self):
        method = self.methodHeatLaminar.currentText() == "Saha-Gaitonde-Date (1989)" or \
            self.methodFrictionLaminar.currentText() == "Saha-Gaitonde-Date (1989)"
        eiamsa = self.groupEiamsa.isVisible() and "ST" in self.modEiamsa.currentText()
        boolean = method or eiamsa or self.helical.isChecked()
        self.lblS.setEnabled(boolean)
        self.S.setEnabled(boolean)
        # self.setVisibleMod()

    def setEnableHollow(self, boolean):
        self.lblC.setVisible(boolean)
        self.C.setVisible(boolean)
        self.changeParams("isHollow", boolean)
        self.setVisibleMod()


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
