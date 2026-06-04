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
         "ref": "Chem. Eng. 218 (1968) 113-119",
         "doi": ""},
    7:
        {"autor": "Kutateladze, S.S., Borishanskii, V.M. ",
         "title": "A Concise Encyclopedia of Heat Transfer",
         "ref": "Pergamon Press (1966)",
         "doi": ""},
    8:
        {"autor": "White, C.M.",
         "title": "Streamline Flow through Curved Pipes",
         "ref": "Proc. R .Soc. London A 123 (1929) 645-663",
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
         "doi": "10.1115/1.2824120"},
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
    18:
        {"autor": "Prasad, B.V.S.S.S., Das, D.H., Prabhaker, A.K.",
         "title": "Pressure Drop, Heat Transfer and Performance of a "
                  "Helically Coiled Tubular Exchanger",
         "ref": "Heat Recovery Systems & CHP 9(3) (1989) 249-256",
         "doi": "10.1016/0890-4332(89)90008-2"},
    19:
        {"autor": "Liu, S., Masliyah, J.H.",
         "title": "Axially invariant laminar flow in helical pipes with a "
                  "finite pitch",
         "ref": "J. Fluid Mech. 251 (1993) 315-353",
         "doi": "10.1017/S002211209300343X"},
    20:
        {"autor": "Seth, K.K., Stahel, E.P.",
         "title": "Heat Transfer from Helical Coils Immersed in Agitated "
                  "Vessels",
         "ref": "Ind. Eng. Chem. 61(6) (1969) 39-49",
         "doi": "10.1021/ie50714a007"},
    21:
        {"autor": "Ali, S.",
         "title": "Pressure drop correlations for flow through regular "
                  "helical coil tubes",
         "ref": "Fluid Dyn. Research 28 (2001) 295-310",
         "doi": "10.1016/s0169-5983(00)00034-4"},
    22:
        {"autor": "Guo, L., Feng, Z, Chen, X.",
         "title": "An experimental investigation of the frictional pressure "
                  "drop of steam–water two-phase flow in helical coils",
         "ref": "Int. J. Heat Mass Transfer 44(14) (2001): 2601-2610",
         "doi": "10.1016/S0017-9310(00)00312-4"},
    23:
        {"autor": "Ito, H.",
         "title": "Laminar Flow in Curved Pipes",
         "ref": "Z. Angew. Math. Mech. 11 (1969) 653-663",
         "doi": "10.1002/zamm.19690491104"},
    24:
        {"autor": "Tarbell, J.M., Samuels, M.R.",
         "title": "Momentum and Heat Transfer in Helical Coils",
         "ref": "Chem. Eng. J. 5(2) (1973) 117-127",
         "doi": "10.1016/0300-9467(73)80002-4"},
    25:
        {"autor": "Mandal, M. M., Nigam, K.D.P.",
         "title": "Experimental Study on Pressure Drop and Heat Transfer of "
                  "Turbulent Flow in Tube in Tube Helical Heat Exchanger",
         "ref": "Ind. Eng. Chem. Res. 48(20) (2009) 9318-9324",
         "doi": "10.1021/ie9002393"},
    26:
        {"autor": "Ciencolini, A., Santini, L.",
         "title": "An experimental investigation regarding the laminar to "
                  "turbulent flow transition in helically coiled pipes",
         "ref": "Exp. Thermal Fluid Sci. 30 (2006) 367-380",
         "doi": "10.1016/j.expthermflusci.2005.08.005"},
    27:
        {"autor": "Adler, M.",
         "title": "Strömung in gekrümmten Rohren",
         "ref": "Z. Angew. Math. Mech. 14(5) 257-275",
         "doi": "10.1002/zamm.19340140502"},
    28:
        {"autor": "Barua, S.N.",
         "title": "On Secondary Flow in Stationary Curved Pipes",
         "ref": "Quart. J. Mech. Appl. Math. 16(1) (1963) 61-77",
         "doi": "10.1093/qjmam_16.1.61"},
    29:
        {"autor": "Pimenta, T.A., Campos, J.B.L.M.",
         "title": "Friction losses of Newtonian and non-Newtonian fluids "
                  "flowing in laminar regime in a helical coil",
         "ref": "Exp. Thermal Fluid Sci. 36 (2012) 194-204",
         "doi": "10.1016/j.expthermflusci.2011.09.013"},
    30:
        {"autor": "Yanase, S., Goto, N., Yamamoto, K.",
         "title": "Dual solutions of the flow through a curved tube",
         "ref": "Fluid Dyn. Research 5 (1989) 191-201",
         "doi": "10.1016/0169-5983(89)90021-x"},
    31:
        {"autor": "Dennis, S.C.R.",
         "title": "Calculation of the steady flow through a curved tube "
                  "using a new finite-difference method",
         "ref": "J. Fluid Mech. 99(3) (1980) 449-467",
         "doi": "10.1017/S0022112080000705"},
    32:
        {"autor": "Van Dyke, M.",
         "title": "Extended Stokes series: laminar flow through a loosely "
                  "coiled pipe",
         "ref": "J. Fluid Mech. 86(1) 129-145",
         "doi": "10.1017/S0022112078001032"},
    33:
        {"autor": "Collins, W.M., Dennis, S.C.R.",
         "title": "The Steady Motion of a Viscous Fluid in a Curved Tube",
         "ref": "Q. J. Mech. Appl. Math. 28(2) (1975) 133-156",
         "doi": "10.1093/qjmam_28.2.133"},
    34:
        {"autor": "Dean, W.R.",
         "title": "The stream-line motion of fluid in a curved pipe "
                  "(Second paper)",
         "ref": "London Edinburgh Dublin Phil. Mag. J. Sci. Serie 7 5(30) "
                "(1928) 673-695",
         "doi": "10.1080/14786440408564513"},
    35:
        {"autor": "Abushammala, O., Hreiz, R., Lamaître, C., Favre, E.",
         "title": "Laminar flow friction factor in highly curved helical "
                  "pipes: Numerical investigation, predictive correlation "
                  "and experimental validation using a 3D-printed model",
         "ref": "Chem. Eng. Sci. 207(7) (2019) 1030-1039",
         "doi": "10.1016/j.ces.2019.07.018"},
    36:
        {"autor": "Kalb, C.E., Seader, J.D.",
         "title": "Heat and Mass Transfer Phenomena for Viscous Flow in "
                  "Curved Circular Tubes",
         "ref": "Int. J. Heat Mass TRansfer 15() (1972) 801-817",
         "doi": "10.1016/0017-9310(72)90122-6"},

    37:
        {"autor": "Dravid, A.N., Smith, K.A., Merrill, E.W., Brian, P.L.T.",
         "title": "Effect of Secondary Fluid Motion on Laminar Flow Heat "
                  "Transfer in Helically Coiled Tubes",
         "ref": "AIChE J. 17(5) (1971) 1114-1122",
         "doi": "10.1002/aic.690170517"},
    38:
        {"autor": "Janssen, L.A.M., Hoogendoorn, C.J.",
         "title": "Laminar Convective Heat Transfer in Helical Coiled Tubes",
         "ref": "Int. J. Heat Mass Transfer 21(9) (1978) 1197-1206",
         "doi": "10.1016/0017-9310(78)90138-2"},
    39:
        {"autor": "Manlapaz, R.L., Churchill, S.W.",
         "title": "Fully Developed Laminar Convection From a Helical Coil",
         "ref": "Chem. Eng. Commun. 9 (1981) 185-200",
         "doi": "10.1080/00986448108911023"},
    40:
        {"autor": "Salimpour, M.R.",
         "title": "Heat transfer coefficients of shell and coiled tube heat "
                  "exchangers",
         "ref": "Exp. Thermal Fluid Sci. 33(2) (2009) 203-207",
         "doi": "10.1016/j.expthermflusci.2008.07.015"},
    41:
        {"autor": "Pimenta, T.A., Campos, J.B.L.M.",
         "title": "Heat transfer coefficients from Newtonian and non-Newtonian"
                  " fluids flowing in laminar regime in a helical coil",
         "ref": "Int. J. Heat Mass Transfer 58 (2013) 676-690",
         "doi": "10.1016/j.ijheatmasstransfer.2012.10.078"},
    42:
        {"autor": "Pawar, S.S., Sunnapwar, V.K.",
         "title": "Studies on convective heat transfer through helical coils",
         "ref": "Heat Mass Transfer 49(12) (2013) 1741-1754",
         "doi": "10.1007/s00231-013-1210-3"},
    43:
        {"autor": "Hardik, B.K., Baburajan, P.K., Prabhu, S.V.",
         "title": "Local heat transfer coefficient in helical coils with "
                  "single phase flow",
         "ref": "Int. J. Heat Mass Transf. 89 (2015) 522-538",
         "doi": "10.1016/j.ijheatmasstransfer.2015.05.069"},
    44:
        {"autor": "Rogers, G.F.C., Mayhew, Y.R.",
         "title": "Heat Transfer and Pressure Loss in Helically Coiled Tubes "
                  "with Turbulent Flow",
         "ref": "Int. J. Heat Mass Transfer 7(11) (1964) 1207-1216",
         "doi": "10.1016/0017-9310(64)90062-6"},
    45:
        {"autor": "Pawar, S.S., Sunnapwar, V.K.",
         "title": "Experimental studies on heat transfer to Newtonian and "
                  "non-Newtonian fluids in helical coils with laminar and "
                  "turbulent flow",
         "ref": "Exp. Thermal Fluid Sci. 44 (2013) 792-804",
         "doi": "10.1016/j.expthermflusci.2012.09.024"},
    46:
        {"autor": "Mori, Y., Nakayama, W.",
         "title": "Study on Forced Convective Heat Transfer in Curved Pipes "
                  "(3rd Report, Theoretical Analysis under the Condition of "
                  "Uniform Wall Temperature and Practical Formulae)",
         "ref": "Int. J. Heat Mass Transfer 10(5) (1967) 681-695",
         "doi": "10.1016_0017-9310(67)90113-5"},
    47:
        {"autor": "Shchukin, V.K.",
         "title": "Correlation of Experimental Data on Heat Transfer in "
                  "Curved Pipes",
         "ref": "Teploenergetika 16(2) (1969) 72-76",
         "doi": ""},
    48:
        {"autor": "Guo, L., Chen, X., Feng, Z., Bai, B.",
         "title": "Transie:nt convective heat transfer in a helical coiled "
                  "tube with pulsatile fully developed turbulent flow",
         "ref": "Int. J. Heat Mass Transfer 41() (1998) 2867-2875",
         "doi": "10.1016/s0017-9310(98)80003-3"},
    49:
        {"autor": "Ghobadi, M., Muzychka, Y.S.",
         "title": "A Review of Heat Transfer and Pressure Drop Correlations "
                  "for Laminar Flow in Curved Circular Ducts",
         "ref": "Heat Transfer Eng. 37(10) (2016) 815-839",
         "doi": "10.1080/01457632.2015.1089735"},
    50:
        {"autor": "Jayakumar, J.S., Mahajani, S.M., Mandal, J.C., Vijayan, "
                  "P.K., Bhoi, R.",
         "title": "Experimental and CFD estimation of heat transfer in "
                  "helically coiled heat exchangers",
         "ref": "Chem. Eng. Res. Design 86(3) (2008) 221-232",
         "doi": "10.1016/j.cherd.2007.10.021"},
    51:
        {"autor": "Yildiz, C., Biçer, Y., Pehlivan, D.",
         "title": "Heat Transfer and Pressure Drop in a Heat Exchanger with "
                  "a Helical Pipe Containing Inside Springs",
         "ref": "Energy Convers. Management 38(6) (1997) 619-624",
         "doi": "10.1016/S0196-8904(96)00040-4"},
    52:
        {"autor": "Wu, Z., Li, K., Zhang, K., Tian, W.",
         "title": "Single-phase flow heat transfer characteristics in helical "
                  "coils with large coil diameters",
         "ref": "Appl. Thermal Eng. 266 (2025) 125776",
         "doi": "10.1016/j.applthermaleng.2025.125776"},
    53:
        {"autor": "Acharya, N., Sen, M., Chang, H.-C.",
         "title": "Analysis of heat transfer enhancement in coiled-tube heat "
                  "exchangers",
         "ref": "Int. J. Heat Mass Transfer 44(17) (2001) 3189-3199",
         "doi": "10.1016/S0017-9310(01)00002-3"},
    54:
        {"autor": "Akiyama, M., Chen g, K.C.",
         "title": "Boundary Vorticity Method for Lamniar Forced Convection "
                  "Heat Transfer in Curved Pipes",
         "ref": "Int. J. Heat Mass Transfer 14(10) (1971) 1659-1675",
         "doi": "10.1016/0017-9310(71)90075-5"},
    55:
        {"autor": "Moawed, M.",
         "title": "Experimental study of forced convection from helical "
                  "coiled tubes with different parameters",
         "ref": "Energy Conv. Management 52(2) (2011) 1150-1156",
         "doi": "10.1016/j.enconman.2010.09.009"},

    # 56:
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


@refDoc(__doi__, [6, 1, 2])
def Rec_Srinivasan(di, Dc):
    r"""Calculates critical Reynolds to define transition between laminar and
    turbulent flow using using the correlation of Srinivasan (1968) as shown in
    [1]_. Recomended method by [2]_.

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


@refDoc(__doi__, [20])
def Rec_SethStahel(di, Dc):
    r"""Calculates critical Reynolds to define transition between laminar and
    turbulent flow using using the correlation of Seth-Stahel (1969).

    .. math::
        Re_c = 1900 \left(1 + 8 \sqrt{\frac{d_i}{D_c}}\right)

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

    # Eq 32
    Rec = 1900 * (1 + 8*(di/Dc)**0.5)
    return Rec


@refDoc(__doi__, [26])
def Rec_Cioncolini(di, Dc):
    r"""Calculates critical Reynolds to define transition between laminar and
    turbulent flow using using the correlation of Cioncolini-Santini (2006).

    .. math::
        Re_c = 30,000 \left(\frac{d_i}{D_c}\right)^{0.47}

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

    # Eq 13
    Rec = 30000 * (di/Dc)**0.47
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
        f = 16/Re * (1+0.14*(di/Dc)**0.97*Re**(1-0.644*(di/Dc)**0.312))
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
        \frac{f_c}{f_s}=\left(\frac{0.108De^{0.5}}{1-3.253 De^{-0.5}}\right)

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
        \frac{f_c}{f_s} = 1 - \left(1-\left(\frac{11.6}{De}\right)^{0.45}
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
        f = fd * (1 + 0.033*log10(De)**4)
    else:
        # Turbulent flow, Eq 10.
        f = fd + 0.03*(di/Dc)**0.5

    return f


@refDoc(__doi__, [18])
def f_Prasad(Re, di, Dc):
    r"""Calculates friction factor for internal flow of a helical coil using
    the method of Prasad et al. (1989).

    For laminar flow use a modified Write correlation:

    .. math::
        \frac{f}{f_s} = \frac{1}{1-\left(1-\left(\frac{B}{De}\right)^{0.45}
        \right)^{\frac{1}{0.45}}}

    For turbulent flow use a modified Ito correlation:

    .. math::
        \frac{f}{f_s} = 1 +
        0.18\left(Re \left(\frac{d_i}{D_c}\right)^2\right)^{0.25}

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
        # Laminar flow
        De = Dean(Re, di, Dc)
        if De < 500:
            B = 11.6
        else:
            B = 6

        f = fd / (1-(1-(B/De)**0.45)**(1/0.45))

    else:
        # Turbulent flow
        f = fd * (1+0.18*(Re*(di/Dc)**2)**0.25)

    return f


@refDoc(__doi__, [21])
def f_Ali(Re, di, Dc, p):
    r"""Calculates friction factor for internal flow of a helical coil using
    the method of Ali (2001).

    Define four flow regime fitted with the correlation:

    .. math::
        f = Eu\frac{d}{L}=\alpha\left(\frac{d_i}{D_{eq}}\right)^{0.15}Re^{\beta}

    with the equivalent diameter of coil:

    .. math::
        D_{eq} = \sqrt{\frac{p^2+\left(\pi D\right)^2}{pi}}

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
    # Equivalent diameter of coil
    Deq = ((p**2 + (pi*Dc)**2)/pi)**0.5

    if Re < 500:
        # Low Laminar flow, Eq 15
        f = 21.88 * (di/Deq)**0.15 / Re**0.9

    elif Re < 6300:
        # Laminar flow, Eq 16
        f = 5.25 * (di/Deq)**0.15 / Re**(2/3)

    elif Re < 10000:
        # Mixed flow, Eq 17
        f = 0.56 * (di/Deq)**0.15 / Re**(2/5)

    else:
        # Turbulent flow, Eq 18
        f = 0.09 * (di/Deq)**0.15 / Re**(1/5)

    return f


@refDoc(__doi__, [1])
def f_ElGenkSchriener(Re, di, Dc, p):
    r"""Calculates friction factor for internal flow of a helical coil using
    the method of ElGenk-Schriener (2017).

    .. math::
        \frac{f_c}{f_s} = 1 + 0.00325 De_m

    where modified Dean number is defined as:

    .. math::
        De_m = De^{0.86} \delta^{0.09} \left(\frac{d_i}{D_c}\right)^{-0.38}

    δ is the curvature defined as:

    .. math::
        \delta = \frac{d_i/D_c}{1+4\pi^2 \tan^2 \alpha}

    α is the helix angle:

    .. math::
        \alpha = \tan^{-1}{\frac{p}{\pi D}}

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
    # Helix angle
    alpha = atan(p/pi/Dc)

    # Curvature
    delta = (di/Dc)/(1+4*pi**2*tan(alpha)**2)

    De = Dean(Re, di, Dc)
    fd = f_friccion(Re)

    # Modified Dean number
    Dem = De**0.86 * delta**0.09 / (di/Dc)**0.38

    # Eq 50
    f = fd * (1+0.00325*Dem)

    return f


@refDoc(__doi__, [6, 49])
def f_Srinivasan(Re, di, Dc):
    r"""Calculates friction factor for internal flow of a helical coil using
    the method of Srinivasan (1968) as explain in [49]_.

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

    if De < 30:
        f = fd

    elif De < 300:
        f = fd * 0.419 * De**0.275

    else:
        f = fd * 0.1125 * De**0.5

    return f


@refDoc(__doi__, [23, 4])
def f_Ito(Re, di, Dc):
    r"""Calculates friction factor for internal flow of a helical coil using
    the method of Ito (1969).

    For laminar flow:

    .. math::
        \frac{f_c}{f_s} = 0.1033 De^{0.5} \left(\left(1+\frac{1.729}{De}\right)
        ^{0.5} - \frac{1.315}{De^{0.5}}\right)^{-3}

    For turbulent flow:

    .. math::
        f_c = 4 \left(0.029 sqrt{\frac{d_i}{D_c}} + 0.304 Re^{-0.25}\right)

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

    if Re < Rec:
        # Laminar flow, Eq 57 from [23]_
        De = Dean(Re, di, Dc)
        fd = f_friccion(Re)

        f = fd * 0.1033 * De**0.5 / ((1+1.729/De)**0.5 - 1.315/De**0.5)**3
    else:
        # Turbulent flow, Eq 2 from [4]_
        f = 0.029*(di/Dc)**0.5 + 0.304/Re**0.25

    return f


@refDoc(__doi__, [8])
def f_laminar_White(Re, di, Dc):
    r"""Calculates friction factor for internal flow of a helical coil in
    laminar flow using the method of White (1929).

    .. math::
        \frac{f_c}{f_s} = 1 - \left(1-\left(\frac{11.6}{De}\right)^{0.45}
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
        \frac{f_c}{f_s} = 1 + 0.09 \frac{De^{1.5}}{70+De}

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


@refDoc(__doi__, [19])
def f_laminar_LiuMasliyah(Re, di, Dc, p):
    r"""Calculates friction factor for internal flow of a helical coil in
    laminar flow using the method of Liu-Masliyah (1993).

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
    Rc = Dc/di
    De = Dean(Re, di, Dc)

    # Curvature ratio, Eq 10
    l = Rc/(Rc**2+(p/2/pi)**2)

    # Torsion, Eq 11
    nu = (p/2/pi)/(Rc**2+(p/2/pi)**2)

    # Eq 43
    f = (16 + (0.378*De*l**0.25 + 12.1)*De**0.5*l**0.5*nu**2) / Re * \
        (1+((0.0908+0.0233*l**0.5)*De**0.5-0.132*l**0.5+0.37*l-0.2)/(1+49/De))

    return f


@refDoc(__doi__, [24])
def f_laminar_TarbellSamuels(Re, di, Dc):
    r"""Calculates friction factor for internal flow of a helical coil in
    laminar flow using the method of Tarbell-Samuels (1973).

    .. math::
        \frac{f_c}{f_s} = 1 + \left(
        8.279e^{-4} + \frac{7.964e{-3}}{d_i/D_c}\right) Re - 2.096e-7 Re^2

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
    fd = f_friccion(Re)

    # Eq 27
    f = fd * (1 + (8.279e-4 + 7.964e-3/(Dc/di))*Re - 2.096e-7*Re**2)
    return f


@refDoc(__doi__, [27])
def f_laminar_Adler(Re, di, Dc):
    r"""Calculates friction factor for internal flow of a helical coil in
    laminar flow using the method of Adler (1934).

    .. math::
        \frac{f_c}{f_s} = 0.1064 \sqrt{De}

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

    # Eq 29
    f = fd * 0.1064*De**0.5
    return f


@refDoc(__doi__, [28])
def f_laminar_Barua(Re, di, Dc):
    r"""Calculates friction factor for internal flow of a helical coil in
    laminar flow using the method of Barua (1963).

    .. math::
        \frac{f_c}{f_s} = 0.509 + 0.0918 \sqrt{De}

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

    # Eq 29
    f = fd * (0.509 + 0.0918*De**0.5)
    return f


@refDoc(__doi__, [29])
def f_laminar_PimentaCampos(Re, di, Dc):
    r"""Calculates friction factor for internal flow of a helical coil in
    laminar flow using the method of Pimenta-Campos (2012)

    .. math::
        \frac{f_c}{f_s} = 1 + \frac{0.028 De^{1.68}}{70+De}

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

    # Eq 10
    f = fd * (1 + 0.028*De**1.68/(70+De))
    return f


@refDoc(__doi__, [30])
def f_laminar_Yanase(Re, di, Dc):
    r"""Calculates friction factor for internal flow of a helical coil in
    laminar flow using the method of Yanase et al. (1989)

    .. math::
        \frac{f_c}{f_s} = 0.557 + 0.0938 \sqrt{De}

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

    # Eq 18
    f = fd * (0.557 + 0.0938*De**0.5)
    return f


@refDoc(__doi__, [31])
def f_laminar_Dennis(Re, di, Dc):
    r"""Calculates friction factor for internal flow of a helical coil in
    laminar flow using the method of Dennis (1980)

    .. math::
        \frac{f_c}{f_s} = 0.388 + 0.1015 \sqrt{De}

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

    # Eq 18
    f = fd * (0.388 + 0.1015*De**0.5)
    return f


@refDoc(__doi__, [32])
def f_laminar_vanDyke(Re, di, Dc):
    r"""Calculates friction factor for internal flow of a helical coil in
    laminar flow using the method of Van Dyke (1978)

    .. math::
        \frac{f_c}{f_s} = 0.47136 De^{0.25}

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

    # Eq 7.6
    f = fd * 0.47136 * De**0.25
    return f


@refDoc(__doi__, [33])
def f_laminar_CollinsDennis(Re, di, Dc):
    r"""Calculates friction factor for internal flow of a helical coil in
    laminar flow using the method of Collins-Dennis (1975)

    .. math::
        \frac{f_c}{f_s} = 0.38036 + 0.1028 \sqrt{De}

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

    # Eq 34
    f = fd * (0.38036 + 0.1028*De**0.5)
    return f


@refDoc(__doi__, [34])
def f_laminar_Dean(Re, di, Dc):
    r"""Calculates friction factor for internal flow of a helical coil in
    laminar flow using the method of Dean (1928)

    .. math::
        \frac{f_c}{f_s} = 0.38036 + 0.1028 \sqrt{De}

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

    Notes
    -----
    Correlation only valid for De < 20

    """
    De = Dean(Re, di, Dc)

    if De > 20:
        raise ValueError("Input out of range")

    fd = f_friccion(Re)

    # Eq 29
    # The K parameter in original paper really is 2*De²
    f = fd * (1 - 0.03058*(De**2/288)**2 + 0.01195*(De**2/288)**4)
    return f


@refDoc(__doi__, [35])
def f_laminar_Abushammala(Re, di, Dc, p):
    r"""Calculates friction factor for internal flow of a helical coil in
    laminar flow using the method of Abushammala et al. (2019)

    .. math::
        f_c = f_s + A B e^{-C}

    with:

    .. math::
        A = p_1 D \left(\frac{D}{Re}\right)^{p_2}

    .. math::
        B = \left(\frac{D_c}{2 d_i} + \frac{2 d_i}{D_c}\right)^{p_3}

    .. math::
        C = p_4 D \frac{p}{d_i} \left(\frac{D_c}{2 d_i}\right)^{-p_5}

    .. math::
        D = \left(\left(\frac{D_c}{2 d_i}\right)^{-p6} \left(1 + \left(
        \frac{p/d_i}{2 \pi D_c/2/d_i}\right)^2\right)\right)^{-p_7}


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

    Notes
    -----
    Correlation only valid for De < 20

    """
    Rh = Dc/2/di
    ph = p/di

    fd = f_friccion(Re)

    # Table 3, parameters
    if Re < 400:
        p = (1.98, 4.07e-1, 8.49e-1, 8.71e-2, 8.91e-1, 2.31, 3.67e-1)
    else:
        p = (2.88, 3.82e-1, 9.16e-3, 2.48e-3, 2.62, 1.1, 3.23e-1)

    # Eq 7
    D = (Rh**-p[5]*(1+(ph/2/pi/Rh)**2))**-p[6]
    C = p[3]*D*ph*Rh**-p[4]
    B = (Rh+1/Rh)**p[2]
    A = p[0]*D*(D/Re)**p[1]
    f = fd + A*B*exp(-C)/4
    return f


@refDoc(__doi__, [14])
def f_turbulent_Czop(Re, di, Dc):
    r"""Calculates friction factor for internal flow of a helical coil in
    turbulent flow using the method of Czop (1994).

    .. math::
        f_c =  \frac{0.096}{De^{-0.1517}}

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


@refDoc(__doi__, [22])
def f_turbulent_Guo(Re, di, Dc):
    r"""Calculates friction factor for internal flow of a helical coil in
    turbulent flow using the method of Guo et al. (2001).

    .. math::
        f_c = 2.552 Re^{-0.15} \left(\frac{d_i}{D_c}\right)^{0.51}

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
    # Eq 20
    # Convert fanning friction factor
    f = 2.552/4 / Re**0.15 * (di/Dc)**0.51
    return f


@refDoc(__doi__, [25])
def f_turbulent_MandalNigam(Re, di, Dc):
    r"""Calculates friction factor for internal flow of a helical coil in
    turbulent flow using the method of Mandal-Nigam (2009).

    .. math::
        \frac{f_c}{f_s} = 1 + 0.03{De}^{0.27}

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

    # Eq 11
    f = fd * (1 + 0.03*De**0.27)
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


@refDoc(__doi__, [9, 10, 46])
def Nu_MoriNakayama(Re, Pr, di, Dc, simple=False):
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

        if simple:
            # Simplified formulae from [46]_

            # Eq 55
            Nu = 0.864/Z * De**0.5 * (1+2.35/De**0.5)

        else:
            # Eq 2.23
            NuI = 0.1979*De**0.5/Z

            if Pr >= 1:
                # Eq 2.24
                f = 1 + 37.05/Z * (1/40 - 17/120*Z + (1/10/Z+13/30)/10/Pr)*De**-0.5
            else:
                # Eq 2.25
                f = 1 - 37.05/Z * (Z**2/12 + 1/24 - 1/120/Z
                                   - (4/3*Z - 1/3/Z + 1/15/Z**2)/20/Pr)*De**-0.5
            Nu = 48/11 * NuI/f


    else:
        # Turbulent flow
        if Pr < 10:
            # Eq 91 in [10]_, for gases
            Nu = Pr/(26.2*(Pr**(2/3)-0.074)) * Re**0.8 * (di/Dc)**0.1 * \
                (1+0.098/(Re*(di/Dc)**2)**0.2)
        else:
            # Eq 94 in [10], for liquids
            Nu = Re**(5/6)/41*(di/Dc)**(1/12)*(1+0.061/(Re*(di/Dc)**2.5)**(1/6))


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


@refDoc(__doi__, [18])
def Nu_Prasad(Re, Pr, di, Dc):
    r"""Calculates nusselt number for internal flow of a helical coil using
    the method of Prasad et al. (1989).

    For laminar flow:

    .. math::
        Nu = 0.25 \left(\frac{f}{8} Re^2\right)^{1/3} Pr^{1/3}

    For turbulent flow:

    .. math::
        Nu = \frac{f}{8} Re


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

    f = f_Prasad(Re, di, Dc)

    if Re < Rec:
        # Laminar flow
        Nu = 0.25 * (f/8*Re**2)**(1/3) * Pr

    else:
        # Turbulent flow
        Nu = f/8 * Re

    return Nu


@refDoc(__doi__, [42, 45])
def Nu_PawarSunnapwar(Re, Pr, di, Dc):
    r"""Calculates nusselt number for internal flow at constant heat flux
    boundary condition of a helical coil using the method of Pawar-Sunnapwar
    (2013).

    For laminar flow:

    .. math::
        Nu = 0.02198 Re^{0.9314} Pr^{0.4} \left(\frac{d_i}{D_c}\right)^{0.391}

    For turbulent flow:

    .. math::
        Nu = 0.0472 De^{0.8346} Pr^{0.4}

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
        # Laminar flow, Eq 23 in [42]_
        Nu = 0.02198 * Re**0.9314 * Pr**0.4 * (di/Dc)**0.391

    else:
        # Turbulent flow, Eq 16 in [45]_
        De = Dean(Re, di, Dc)
        Nu = 0.0472 * De**0.8346 * Pr**0.4

    return Nu


@refDoc(__doi__, [1])
def Nu_ElGenkSchriener(Re, Pr, di, Dc, p):
    r"""Calculates nusselt number for internal flow of a helical coil using
    the method of ElGenk-Schriener (2017).

    For fluids with Pr < 15:

    .. math::
        Nu_c = 3.66 + 0.014 Re_m^{0.86} Pr^{0.4}

    For fluids with Pr > 15:

    .. math::
        Nu_c = 3.66 + 0.02 Re_m^{0.7} Pr^{0.4}

    using a modified Reynolds number:

    .. math::

        Re_m = Re \left(1+3.4 \delta\right)

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
    p : float
        Pitch for twist of 2π radians (360º), [m]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    # Helix angle
    alpha = atan(p/pi/Dc)

    # Curvature
    delta = (di/Dc)/(1+4*pi**2*tan(alpha)**2)

    Rem = Re * (1+3.4*delta)

    if Pr < 15:
        # Eq 51
        Nu = 3.66 + 0.014*Re**0.86*Pr**0.4
    else:
        # Eq 52
        Nu = 3.66 + 0.02*Re**0.7*Pr**0.4

    return Nu


@refDoc(__doi__, [36])
def Nu_laminar_KalbSeader(Re, Pr, di, Dc):
    r"""Calculates nusselt number for internal flow at constant heat flux
    boundary condition of a helical coil in laminar flow using the method of
    Kalb-Seader (1972).

    For Pr < 0.05:

    .. math::
        Nu = 3.31 De^{0.115} Pr^{0.0108}

    For Pr > 0.7

    .. math::
        Nu = 0.913 De^{0.476} Pr^{0.2}

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
    De = Dean(Re, di, Dc)

    if Pr < 0.5:
        # Eq 22
        Nu = 3.31 * De**0.115 * Pr**0.0108
    else:
        # Eq 23
        Nu = 0.913 * De**0.476 * Pr**0.2

    return Nu


@refDoc(__doi__, [37])
def Nu_laminar_Dravid(Re, Pr, di, Dc):
    r"""Calculates nusselt number for internal flow at constant heat flux
    boundary condition of a helical coil in laminar flow using the method of
    Dravid et al. (1971).

    .. math::
        Nu = \left(0.76 + 0.65 De^{0.5}\right) Pr^{0.175}

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
    De = Dean(Re, di, Dc)

    # Eq 23
    Nu = (0.76 + 0.65*De**0.5) * Pr**0.175

    return Nu


@refDoc(__doi__, [38])
def Nu_laminar_JanssenHoogendoorn(Re, Pr, di, Dc, f):
    r"""Calculates nusselt number for internal flow at constant heat flux
    boundary condition of a helical coil in laminar flow using the method of
    Janssen-Hoogendoorn (1978).

    .. math::
        Nu = 0.6166 \left(f Re^2\right)^{0.26} Pr^{1/6}

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
    f : float
        Friction factor, [-]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """

    De = Dean(Re, di, Dc)

    if De > 20:
        # Eq 19
        # The original equation use the Darcy-Weisbach friction factor, so
        # convert to Fanning:
        # 0.43*4**0.26
        Nu = 0.6166 * (f*Re**2)**0.26 * Pr**(1/6)
    else:
        # Eq 22
        Nu = 1.7 * (De**2*Pr)**(1/6)

    return Nu


@refDoc(__doi__, [39])
def Nu_laminar_ManlapazChurchill(Re, Pr, di, Dc, p):
    r"""Calculates nusselt number for internal flow at constant heat flux
    boundary condition of a helical coil in laminar flow using the method of
    Manlapaz-Churchill (1981).

    .. math::
        Nu = \left(\left(3.657 + \frac{4.343}{\left(1+\frac{957}{Pr He^2}
        \right)^2}\right)^3 + 1.158 \left(\frac{He}{1+\frac{0.477}{Pr}}\right)
        ^{1.5}\right)^{1/3}

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
    p : float, optional
        Pitch for twist of 2π radians (360º), [m]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """

    De = Dean(Re, di, Dc)
    He = De/(1+(p/2/pi/di)**2)**0.5

    # Eq 39, Uniform wall temperature
    Nu = ((3.657 + 4.343/(1+957/Pr/He**2)**2)**3
          + 1.158*(He/(1+0.477/Pr))**1.5)**(1/3)

    # Paper give too a correlation for uniform heat flux

    return Nu


@refDoc(__doi__, [40])
def Nu_laminar_Salimpour(Re, Pr, di, Dc, p):
    r"""Calculates nusselt number for internal flow at constant heat flux
    boundary condition of a helical coil in laminar flow using the method of
    Salimpour (2009).

    .. math::
        Nu = 0.152De^{0.431}Pr^{1.06}\left(\frac{b}{2 \pi D_c}\right)^{-0.277}

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
    p : float, optional
        Pitch for twist of 2π radians (360º), [m]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """

    De = Dean(Re, di, Dc)

    # Eq 5
    Nu = 0.152 * De**0.431 * Pr**1.06 * (2*pi*di/p)**0.277

    return Nu


@refDoc(__doi__, [41])
def Nu_laminar_PimentaCampos(Re, Pr, di, Dc):
    r"""Calculates nusselt number for internal flow at constant heat flux
    boundary condition of a helical coil in laminar flow using the method of
    Pimenta-Campos (2013).

    .. math::
        Nu = \left(0.5 De^{0.481} - 0.465\right) Pr^{0.367}

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

    De = Dean(Re, di, Dc)

    # Eq 38
    Nu = (0.5*De**0.481 - 0.465) * Pr**0.367

    return Nu


@refDoc(__doi__, [43])
def Nu_laminar_Hardik(Re, Pr, di, Dc):
    r"""Calculates nusselt number for internal flow at constant heat flux
    boundary condition of a helical coil in laminar flow using the method of
    Hardik et al. (2015).

    .. math::
        Nu = 0.0456 Re^{0.8} Pr^{0.4} \left(\frac{d_i}{D_c}\right)^{0.16}

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
    # Eq 22
    Nu = 0.0456 * (di/Dc)**0.16 * Re**0.8 * Pr**0.4

    return Nu


@refDoc(__doi__, [53])
def Nu_laminar_Acharya(Re, Pr, di, Dc, AA=False):
    r"""Calculates friction factor for internal flow of a helical coil in
    laminar flow using the method of Acharya et al. (2001)

    For Pr > 1:

    .. math::
        Nu = 0.67 Re^{0.5} Pr^{0.21} \left(\frac{d_i}{D_c}\right)^{0.13}

    For Pr ≤ 1:

    .. math::
        Nu = 0.69 Re^{0.5} Pr^{0.43} \left(\frac{d_i}{D_c}\right)^{0.13}


    Include too correlations for alternate axis coil geometric configuration

    For Pr > 1:

    .. math::
        Nu = 0.7 Re^{0.5} Pr^{0.3} \left(\frac{d_i}{D_c}\right)^{0.18}

    For Pr ≤ 1:

    .. math::
        Nu = 0.7 Re^{0.5} Pr^{0.375} \left(\frac{d_i}{D_c}\right)^{0.18}


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
    AA : boolean
        Set alternalte axis configuraiton for helical coil

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    if AA:
        if Pr <= 1:
            # Eq 15
            Nu = 0.7 * Re**0.5 * Pr**0.375 * (di/Dc)**0.18
        else:
            # Eq 16
            Nu = 0.7 * Re**0.5 * Pr**0.3 * (di/Dc)**0.18

    else:
        if Pr <= 1:
            # Eq 17
            Nu = 0.69 * Re**0.5 * Pr**0.43 * (di/Dc)**0.13
        else:
            # Eq 18
            Nu = 0.67 * Re**0.5 * Pr**0.21 * (di/Dc)**0.13

    return Nu


@refDoc(__doi__, [54])
def Nu_laminar_AkiyamaCheng(Re, Pr, di, Dc):
    r"""Calculates friction factor for internal flow of a helical coil in
    laminar flow using the method of Akiyama-Cheng (1971)

    .. math::
        \frac{Nu}{Nu_o} = 0.181 Q
        \left(1 - 0.839Q^{-1} + 35.4Q^{-2} - 207Q^{-3} + 419Q^{-4}\right)

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
    De = Dean(Re, di, Dc)

    Q = (De**2*Pr)**0.25

    # Eq 20
    Nur = 0.181*Q*(1 - 0.839/Q + 35.4/Q**2 -207/Q**3 + 419/Q**4)

    return Nur * 48/11


@refDoc(__doi__, [55])
def Nu_laminar_Moawed(Re, do, Dc, p):
    r"""Calculates friction factor for internal flow of a helical coil in
    laminar flow using the method of Moawed (2011)

    .. math::
        Nu = 0.0345 Re^{0.48} \left(\frac{D_c}{d_o}\right)^{0.914}
        \left(\frac{p}{d_o}\right)^{0.281}

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    do : float
        Outer diameter of the pipe, [m]
    Dc : float
        Diameter of the helix, [m]
    p : float, optional
        Pitch for twist of 2π radians (360º), [m]

    This correlation use external diameter of pipe, without Prandtl dependence

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    # Eq 8
    Nu = 0.0345 * Re**0.48 * (Dc/do)**0.914 * (p/do)**0.281

    return Nu


@refDoc(__doi__, [25])
def Nu_turbulent_MandalNigam(Re, Pr, di, Dc):
    r"""Calculates nusselt number for internal flow of a helical coil in
    turbulent flow using the method of Mandal-Nigam (2009).

    .. math::
        Nu = 0.55 De^{0.637} Pr^{0.4}

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
    De = Dean(Re, di, Dc)

    # Eq 13
    Nu = 0.55 * De**0.637 * Pr**0.4
    return Nu


@refDoc(__doi__, [44])
def Nu_turbulent_RogersMayhew(Re, Pr, di, Dc):
    r"""Calculates nusselt number for internal flow of a helical coil in
    turbulent flow using the method of Rogers-Mayhew (1964).

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
    # Eq 12
    Nu = 0.023 * Re**0.85 * Pr**0.4 * (di/Dc)**0.1
    return Nu


@refDoc(__doi__, [44, 1])
def Nu_turbulent_Shchukin(Re, Pr, di, Dc):
    r"""Calculates nusselt number for internal flow of a helical coil in
    turbulent flow using the method of Shchukin (1969) as show in [1]_

    For math:`Re (d_i/D_c)^2 < 20`

    .. math::
        Nu = 0.0316 Re^{0.8} Pr^{0.4} \left(\frac{d_i}{D_c}\right)^{0.05}

    For math:`Re (d_i/D_c)^2 > 20`

    .. math::
        Nu = 0.0266 Re^{0.85} Pr^{0.4} \left(\frac{d_i}{D_c}\right)^{0.15}

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
    if Re*(di/Dc)**2 < 20:
        Nu = 0.0316 * Re**0.8 * Pr**0.4 * (di/Dc)**0.05
    else:
        Nu = 0.0266 * Re**0.85 * Pr**0.4 * (di/Dc)**0.15

    return Nu


@refDoc(__doi__, [48])
def Nu_turbulent_Guo(Re, Pr):
    r"""Calculates nusselt number for internal flow of a helical coil in
    turbulent flow using the method of Guo (1998).

    .. math::
        Nu = 0.023 Re^{0.58} Pr^{0.4}

    This correlation don't include any helical coil geometrical parameters
    dependence.

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    Pr : float
        Prandtl number, [-]

    Returns
    -------
    Nu : float
        Nusselt number, [-]
    """
    # Eq 12
    Nu = 0.328 * Re**0.58 * Pr**0.4
    return Nu


@refDoc(__doi__, [50])
def Nu_turbulent_Jayakumar(Re, Pr, di, Dc):
    r"""Calculates nusselt number for internal flow of a helical coil in
    turbulent flow using the method of Jayakumar et al. (2008)

    .. math::
        Nu = 0.025 De^{0.9112} Pr^{0.4}

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
    De = Dean(Re, di, Dc)

    # Eq 8
    Nu = 0.025 * De**0.9112 * Pr**0.4

    return Nu


@refDoc(__doi__, [51])
def Nu_turbulent_Yildiz(Re, Pr, di, Dc):
    r"""Calculates nusselt number for internal flow of a helical coil in
    turbulent flow using the method of Yildiz et al. (1997)

    .. math::
        Nu = 0.0551 De^{0.864} Pr^{0.4}

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
    De = Dean(Re, di, Dc)

    # Eq 5
    Nu = 0.0551 * De**0.864 * Pr**0.4

    return Nu


@refDoc(__doi__, [52])
def Nu_turbulent_Wu(Re, Pr, di, Dc):
    r"""Calculates nusselt number for internal flow of a helical coil in
    turbulent flow using the method of Wu et al. (2025)

    .. math::
        Nu = 0.023 Re^{0.759} Pr^{0.4} \left(\frac{d_i}{D_c}\right)^{-0.079}

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

    # Eq 26
    Nu = 0.023 * Re**0.759 * Pr**0.4 / (di/Dc)**0.079

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
    MoriSimple : boolean, optional
        Use Simple correlation for Mori-Nakayama nusselt number correlation
    AA : boolean, optional
        Use alternate axis configuration for Acharya nusselt number correlation
    """

    TEXT_REYNOLDS_CRITICAL = (
        "Ito (1959)",
        "Schmidt (1967)",
        "Kubair-Kuloor (1966)",
        "Srinivasan (1968)",
        "Kutateladze (1966)",
        "Seth-Stahel (1969)",
        "Cioncolini-Santini (2006)")

    TEXT_LAMINAR_FRICTION = (
        "Schmidt (1967)",
        "White (1929)",
        "Mori-Nakayama (1965)",
        "Hart (1988)",
        "Ju (2001)",
        "Mishra-Gupta (1979)",
        "Manlapaz-Churchill (1980)",
        "Prasad (1989)",
        "Liu-Masliyah (1993)",
        "Ali (2001)",
        "Ito (1969)",
        "Tarbell-Samuels (1973)",
        "Pimenta-Campos (2012)",
        "Adler (1934)",
        "Barua (1963)",
        "Yanase (1989)",
        "Dennis (1980)",
        "van Dyke (1978)",
        "Collins-Dennis (1975)",
        "Dean (1928)",
        "Abushammala (2019)",
        "ElGenk-Schriener (2017)",
        "Srinivasan (1968)",
    )

    TEXT_TURBULENT_FRICTION = (
        "Schmidt (1967)",
        "Mori-Nakayama (1965)",
        "Ju (2001)",
        "Mishra-Gupta (1979)",
        "Czop (1994)",
        "Prasad (1989)",
        "Ali (2001)",
        "Guo (2001)",
        "Mandal-Nigam (2009)",
        "ElGenk-Schriener (2017)",
        "Srinivasan (1968)",
        "Ito (1959)",
    )

    TEXT_LAMINAR_HEAT = (
        "Schmidt (1967)",
        "Xin-Ebadian (1997)",
        "Mori-Nakayama (1965)",
        "Seban-McLaughlin (1963)",
        "Prasad (1989)",
        "Kalb-Seader (1972)",
        "Dravid (1971)",
        "Janssen-Hoogendoorn (1978)",
        "Manlapaz-Churchill (1981)",
        "Salimpour (2009)",
        "Pimenta-Campos (2013)",
        "Pawar-Sunnapwar (2013)",
        "Hardik (2015)",
        "ElGenk-Schriener (2017)",
        "Acharya (2001)",
        "Akiyama-Cheng (1971)",
        "Moawed (2011)",
    )

    TEXT_TURBULENT_HEAT = (
        "Schmidt (1967)",
        "Xin-Ebadian (1997)",
        "Mori-Nakayama (1965)",
        "Seban-McLaughlin (1963)",
        "Prasad (1989)",
        "Mandal-Nigam (2009)",
        "Rogers-Mayhew (1964)",
        "Pawar-Sunnapwar (2013)",
        "ElGenk-Schriener (2017)",
        "Shchukin (1969)",
        "Guo (1998)",
        "Jayakumar (2008)",
        "Yildiz (1997)",
        "Wu (2025)",
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
        "p": 0,

        "MoriSimple": False,
        "AA": False
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
            # Schmidt (1967)
            Rec = Rec_Schmidt(self.di, self.Dc)

        elif self.kw["methodReCritic"] == 2:
            # Kubair-Kuloor (1966)
            Rec = Rec_Kubair(self.di, self.Dc)

        elif self.kw["methodReCritic"] == 2:
            # Srinivasan (1968)
            Rec = Rec_Srinivasan(self.di, self.Dc)

        elif self.kw["methodReCritic"] == 3:
            # Kutateladze (1966)
            Rec = Rec_Kutateladze(self.di, self.Dc)

        elif self.kw["methodReCritic"] == 4:
            #  Seth-Stahel (1969)
            Rec = Rec_SethStahel(self.di, self.Dc)

        elif self.kw["methodReCritic"] == 5:
            # Cioncolini-Santini (2006)
            Rec = Rec_Cioncolini(self.di, self.Dc)

        else:
            # Ito (1959)
            Rec = Rec_Ito(self.di, self.Dc)

        return Rec

    def Nu(self, Re, Pr):
        """Calculate nusselt number"""
        msg = ""
        Rec = self.ReCritical

        if Re < Rec:
            # Laminar flow
            if self.kw["methodHeatLaminar"] == 1:
                # Xin-Ebadian (1997)
                Nu = Nu_XinEbadian(Re, Pr, self.di, self.Dc)

            elif self.kw["methodHeatLaminar"] == 2:
                # Mori-Nakayama (1965)
                Nu = Nu_MoriNakayama(
                    Re, Pr, self.di, self.Dc, self.kw["MoriSimple"])

            elif self.kw["methodHeatLaminar"] == 3:
                # Seban-McLaughlin (1963)
                Nu = Nu_SebanMcLaughlin(Re, Pr, self.di, self.Dc)

            elif self.kw["methodHeatLaminar"] == 4:
                # Prasad (1989)
                Nu = Nu_Prasad(Re, Pr, self.di, self.Dc)

            elif self.kw["methodHeatLaminar"] == 5:
                # Kalb-Seader (1972)
                Nu = Nu_laminar_KalbSeader(Re, Pr, self.di, self.Dc)

            elif self.kw["methodHeatLaminar"] == 6:
                # Dravid (1971)
                Nu = Nu_laminar_Dravid(Re, Pr, self.di, self.Dc)

            elif self.kw["methodHeatLaminar"] == 7:
                # Janssen-Hoogendoorn (1978)
                f = self.f(Re)
                Nu = Nu_laminar_JanssenHoogendoorn(Re, Pr, self.di, self.Dc, f)

            elif self.kw["methodHeatLaminar"] == 8:
                # Manlapaz-Churchill (1981)
                Nu = Nu_laminar_ManlapazChurchill(
                    Re, Pr, self.di, self.Dc, self.kw["p"])

            elif self.kw["methodHeatLaminar"] == 9:
                # Salimpour (2009)
                if self.kw["p"]:
                    Nu = Nu_laminar_Salimpour(
                        Re, Pr, self.di, self.Dc, self.kw["p"])
                else:
                    Nu = Nu_Schmidt(Re, Pr, self.di, self.Dc)
                    msg = "Helical pitch undefined, using Schmidt correlation"
                    msg += "instead."

            elif self.kw["methodHeatLaminar"] == 10:
                # Pimenta-Campos (2013)
                Nu = Nu_laminar_PimentaCampos(Re, Pr, self.di, self.Dc)

            elif self.kw["methodHeatLaminar"] == 11:
                # Pawar-Sunnapwar (2013)
                Nu = Nu_PawarSunnapwar(Re, Pr, self.di, self.Dc)

            elif self.kw["methodHeatLaminar"] == 12:
                # Hardik (2015)
                Nu = Nu_laminar_Hardik(Re, Pr, self.di, self.Dc)

            elif self.kw["methodHeatLaminar"] == 13:
                # ElGenk-Schriener (2017)
                Nu = Nu_ElGenkSchriener(Re, Pr, self.di, self.Dc, self.kw["p"])

            elif self.kw["methodHeatLaminar"] == 14:
                # Acharya (2001)
                Nu = Nu_laminar_Acharya(
                    Re, Pr, self.di, self.Dc, self.kw["AA"])

            elif self.kw["methodHeatLaminar"] == 15:
                # Akiyama-Cheng (1971)
                Nu = Nu_laminar_AkiyamaCheng(Re, Pr, self.di, self.Dc)

            elif self.kw["methodHeatLaminar"] == 16:
                # Moawed (2011)
                if self.kw["p"]:
                    Nu = Nu_laminar_Moawed(Re, self.di, self.Dc, self.kw["p"])
                else:
                    Nu = Nu_Schmidt(Re, Pr, self.di, self.Dc)
                    msg = "Helical pitch undefined, using Schmidt correlation"
                    msg += "instead."

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
                Nu = Nu_MoriNakayama(
                    Re, Pr, self.di, self.Dc, self.kw["MoriSimple"])

            elif self.kw["methodHeatTurbulent"] == 3:
                # Seban-McLaughlin (1963)
                Nu = Nu_SebanMcLaughlin(Re, Pr, self.di, self.Dc)

            elif self.kw["methodHeatTurbulent"] == 4:
                # Prasad (1989)
                Nu = Nu_Prasad(Re, Pr, self.di, self.Dc)

            elif self.kw["methodHeatTurbulent"] == 5:
                # Mandal-Nigam (2009)
                Nu = Nu_turbulent_MandalNigam(Re, Pr, self.di, self.Dc)

            elif self.kw["methodHeatTurbulent"] == 6:
                # Rogers-Mayhew (1964)
                Nu = Nu_turbulent_RogersMayhew(Re, Pr, self.di, self.Dc)

            elif self.kw["methodHeatTurbulent"] == 7:
                # Pawar-Sunnapwar (2013)
                Nu = Nu_PawarSunnapwar(Re, Pr, self.di, self.Dc)

            elif self.kw["methodHeatTurbulent"] == 8:
                # ElGenk-Schriener (2017)
                Nu = Nu_ElGenkSchriener(Re, Pr, self.di, self.Dc, self.kw["p"])

            elif self.kw["methodHeatTurbulent"] == 9:
                # Shchukin (1969)
                Nu = Nu_turbulent_Shchukin(Re, Pr, self.di, self.Dc)

            elif self.kw["methodHeatTurbulent"] == 10:
                # Guo (1998)
                Nu = Nu_turbulent_Guo(Re, Pr)

            elif self.kw["methodHeatTurbulent"] == 11:
                # Jayakumar (2008)
                Nu = Nu_turbulent_Jayakumar(Re, Pr, self.di, self.Dc)

            elif self.kw["methodHeatTurbulent"] == 12:
                # Yildiz (1997)
                Nu = Nu_turbulent_Yildiz(Re, Pr, self.di, self.Dc)

            elif self.kw["methodHeatTurbulent"] == 13:
                # Wu (2025)
                Nu = Nu_turbulent_Wu(Re, Pr, self.di, self.Dc)

            else:
                # Schmidt (1967)
                Nu = Nu_Schmidt(Re, Pr, self.di, self.Dc)

        if msg:
            self.status = 3
            self.msg = translate("equipment", msg)
            self.inputChanged.emit(self)

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

            elif self.kw["methodFrictionLaminar"] == 7:
                # Prasad (1989)
                f = f_Prasad(Re, self.di, self.Dc)

            elif self.kw["methodFrictionLaminar"] == 8:
                # Liu-Masliyah (1993)
                f = f_laminar_LiuMasliyah(Re, self.di, self.Dc, self.kw["p"])

            elif self.kw["methodFrictionLaminar"] == 9:
                # Ali (2001)
                f = f_Ali(Re, self.di, self.Dc, self.kw["p"])

            elif self.kw["methodFrictionLaminar"] == 10:
                # Ito (1969)
                f = f_Ito(Re, self.di, self.Dc)

            elif self.kw["methodFrictionLaminar"] == 11:
                # Tarbell-Samuels (1973)
                f = f_laminar_TarbellSamuels(Re, self.di, self.Dc)

            elif self.kw["methodFrictionLaminar"] == 12:
                # Pimenta-Campos (2012)
                f = f_laminar_PimentaCampos(Re, self.di, self.Dc)

            elif self.kw["methodFrictionLaminar"] == 13:
                # Adler (1934)
                f = f_laminar_Adler(Re, self.di, self.Dc)

            elif self.kw["methodFrictionLaminar"] == 14:
                # Barua (1963)
                f = f_laminar_Barua(Re, self.di, self.Dc)

            elif self.kw["methodFrictionLaminar"] == 15:
                # Yanase (1989)
                f = f_laminar_Yanase(Re, self.di, self.Dc)

            elif self.kw["methodFrictionLaminar"] == 16:
                # Dennis (1980)
                f = f_laminar_Dennis(Re, self.di, self.Dc)

            elif self.kw["methodFrictionLaminar"] == 17:
                # Van Dyke (1978)
                f = f_laminar_vanDyke(Re, self.di, self.Dc)

            elif self.kw["methodFrictionLaminar"] == 18:
                # Collins-Dennis (1975)
                f = f_laminar_CollinsDennis(Re, self.di, self.Dc)

            elif self.kw["methodFrictionLaminar"] == 19:
                # Dean (1928)
                try:
                    f = f_laminar_Dean(Re, self.di, self.Dc)
                except ValueError:
                    f = f_Schmidt(Re, self.di, self.Dc)
                    msg = "Dean correlation out of range, using Schmidt instead"

            elif self.kw["methodFrictionLaminar"] == 20:
                # Abushammala (2019)
                f = f_laminar_Abushammala(Re, self.di, self.Dc, self.kw["p"])

            elif self.kw["methodFrictionLaminar"] == 21:
                # ElGenk-Schriener (2017)
                f = f_ElGenkSchriener(Re, self.di, self.Dc, self.kw["p"])

            elif self.kw["methodFrictionLaminar"] == 22:
                # Srinivasan (1968)
                f = f_Srinivasan(Re, self.di, self.Dc)

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

            elif self.kw["methodFrictionTurbulent"] == 5:
                # Prasad (1989)
                f = f_Prasad(Re, self.di, self.Dc)

            elif self.kw["methodFrictionTurbulent"] == 6:
                # Ali (2001)
                f = f_Ali(Re, self.di, self.Dc, self.kw["p"])

            elif self.kw["methodFrictionTurbulent"] == 7:
                # Guo (2001)
                f = f_turbulent_Guo(Re, self.di, self.Dc)

            elif self.kw["methodFrictionTurbulent"] == 8:
                # Mandal-Nigam (2009)
                f = f_turbulent_MandalNigam(Re, self.di, self.Dc)

            elif self.kw["methodFrictionTurbulent"] == 9:
                # ElGenk-Schriener (2017)
                f = f_ElGenkSchriener(Re, self.di, self.Dc, self.kw["p"])

            elif self.kw["methodFrictionTurbulent"] == 10:
                # Srinivasan (1968)
                f = f_Srinivasan(Re, self.di, self.Dc)

            elif self.kw["methodFrictionTurbulent"] == 11:
                # Ito (1959)
                f = f_Ito(Re, self.di, self.Dc)

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
        self.methodHeatLaminar.currentTextChanged.connect(self.setVisibleMod)
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

        # Mori-Nakayama additional parameters
        self.MoriSimple = QtWidgets.QCheckBox(self.tr(
            "Use simple correlation for laminar nusselt number"))
        self.MoriSimple.toggled.connect(
            partial(self.changeParams, "MoriSimple"))
        lyt.addWidget(self.MoriSimple, 7, 1, 1, 2)

        # Acharya additional parameters
        self.AA = QtWidgets.QCheckBox(self.tr(
            "Use alternate axis geometric configuration"))
        self.AA.toggled.connect(
            partial(self.changeParams, "AA"))
        lyt.addWidget(self.AA, 8, 1, 1, 2)

        self.Entity.valueChanged.connect(self.valueChanged.emit)
        self.Entity.inputChanged.connect(self.populate)
        self.setVisibleMod()

    def setVisibleMod(self):
        """Enable widget with special parameters for selected method"""
        # Mori-Nakayama
        self.MoriSimple.setVisible(
            self.methodHeatLaminar.currentText() == "Mori-Nakayama (1965)")

        # Acharya
        self.AA.setVisible(
            self.methodHeatLaminar.currentText() == "Acharya (2001)")


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
