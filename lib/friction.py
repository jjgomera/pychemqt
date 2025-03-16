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
along with this program.  If not, see <http://www.gnu.org/licenses/>.


Module for implement friction factor related functionality

**Friction factor**

Global function with all functionality, :func:`f_friccion`

Friction factor for rough pipes, Colebrook-White intrinsic equation is solved
by iteration so so many method have been implement to get direct equations:

    * :func:`f_colebrook`
    * :func:`f_chen`
    * :func:`f_Vatankhah`
    * :func:`f_buzzelli`
    * :func:`f_romeo`
    * :func:`f_serghides`
    * :func:`f_zigrang`
    * :func:`f_Samadianfard`
    * :func:`f_brkic`
    * :func:`f_fang`
    * :func:`f_ghanbari`
    * :func:`f_haaland`
    * :func:`f_round`
    * :func:`f_swamee`
    * :func:`f_jain`
    * :func:`f_barr`
    * :func:`f_shacham`
    * :func:`f_tsal`
    * :func:`f_manadilli`
    * :func:`f_goudar`
    * :func:`f_goudar2007`
    * :func:`f_avci`
    * :func:`f_papaevangelou`
    * :func:`f_churchill`
    * :func:`f_chen1979`
    * :func:`f_moody`
    * :func:`f_wood`
    * :func:`f_eck`
    * :func:`f_altshul`

'''


from math import log, log10

from scipy.optimize import fsolve

from lib.unidades import Dimensionless
from lib.utilities import refDoc


__doi__ = {
    1:
        {"autor": "Colebrook, C.F., White, C.M.",
         "title": "Experiments with Fluid Friction in Roughened Pipes",
         "ref": "Proc. R. Soc. Lond. A 161 (1937) 367-381.",
         "doi": "10.1098/rspa.1937.0150"},
    2:
        {"autor": "Chen, H.J.",
         "title": "An Explicit Equation for Friction Factor in Pipe",
         "ref": "Ind. Eng. Chem. Fundam. 18(3) (1979) 296-297",
         "doi": "10.1021/i160071a019"},
    3:
        {"autor": "Chen, H.J.",
         "title": "An Exact Solution to the Colebrook Equation",
         "ref": "Chem. Eng. 94(2) (1987) 196-198",
         "doi": ""},
    4:
        {"autor": "Moody, L. F.",
         "title": "An approximate formula for pipe friction factors",
         "ref": "Trans. ASME, 69(12) (1947) 1005-1006.",
         "doi": ""},
    5:
        {"autor": "Churchill, S.W.",
         "title": "Friction-factor equation spans all fluid-flow regimes",
         "ref": "Chem. Eng. 84 (1977) 94-95",
         "doi": ""},
    6:
        {"autor": "Wood D.J.",
         "title": "An explicit friction factor relationship",
         "ref": "Civil Eng. ASCE 60, 1966",
         "doi": ""},
    7:
        {"autor": "Haaland, S.E.",
         "title": "Simple and explicit formulas for the friction factor in"
                  "turbulent flow",
         "ref": "J. Fluids Eng., 105(1) (1983) 89-90.",
         "doi": "10.1115/1.3240948"},
    8:
        {"autor": "Serghides, T.K.",
         "title": "Estimate friction factor accurately",
         "ref": "Chem. Eng., 91(5) (1984) 63-64.",
         "doi": ""},
    9:
        {"autor": "Round, G.F.",
         "title": "An Explicit Approximation for the Friction Factor-Reynolds"
                  "Number Relation for Rough and Smooth Pipes",
         "ref": "Can. J. Chem. Eng. 58 (1980) 122-123",
         "doi": "10.1002/cjce.5450580119"},
    10:
        {"autor": "Swamee, P.K.; Jain, A.K.",
         "title": "Explicit equations for pipe-flow problems",
         "ref": "J. Hydraulics Division (ASCE) 102(5) (1976) 657-664.",
         "doi": ""},
    11:
        {"autor": "Jain, A.K.",
         "title": "Accurate Explicit Equation for Friction Factor",
         "ref": "J. Hydraulics Division 102(5) (1976) 674-77",
         "doi": ""},
    12:
        {"autor": "Barr, D.I.H.",
         "title": "Solutions of the Colebrook-White functions for resistance "
                  "to uniform turbulent flows.",
         "ref": "Proc Inst Civil Eng 71, 1981, 529-536.",
         "doi": "10.1680/iicep.1981.1895"},
    13:
        {"autor": "Zigrang, D.J., Sylvester, N.D.",
         "title": "Explicit approximations to the solution of Colebrook's"
                  "friction factor equation",
         "ref": "AICHE J 28, 514-515.",
         "doi": "10.1002/aic.690280323"},
    14:
        {"autor": "Tsal, R.J.",
         "title": "Altshul-Tsal friction factor equation",
         "ref": "Heating Piping Air Conditioning 8 (1989), 30-45.",
         "doi": ""},
    15:
        {"autor": "Eck, B.",
         "title": "Technische Stromungslehre",
         "ref": "Springer, New York, 1973.",
         "doi": ""},
    16:
        {"autor": "Shacham. M.",
         "title": "An explicit equation for friction factor in pipe",
         "ref": "Ind. Eng. Chem. Fund. 19 (1981) 228-229.",
         "doi": ""},
    17:
        {"autor": "Manadilli, G.",
         "title": "Replace implicit equations with signomial functions.",
         "ref": "Chem. Eng. 104 (1997) 129-132.",
         "doi": ""},
    18:
        {"autor": "Romeo, E., Royo, C., Monzon, A.",
         "title": "Improved explicit equation for estimation of the friction"
                  "factor in rough and smooth pipes.",
         "ref": "Chem. Eng. J. 86(3) (2002) 369-374",
         "doi": "10.1016/S1385-8947(01)00254-6"},
    19:
        {"autor": "Sonnad, J.R., Goudar, C.T.",
         "title": "Explicit Reformulation of the Colebrook-White Equation for "
                  "Turbulent Flow Frcition Factor Calculation",
         "ref": "Ind. Eng. Chem. Res. 46(8) (2007) 2593-2600",
         "doi": "10.1021/ie0340241"},
    20:
        {"autor": "Buzzelli, D.",
         "title": "Calculating friction in one step",
         "ref": "Machine Design, 80 (2008), 54–55.",
         "doi": ""},
    21:
        {"autor": "Vatankhah, A.R., Kouchakzadeh, S",
         "title": "Full-range pipe-flow equations",
         "ref": "Journal of Hydraulic Research 46(4) (2008) 559",
         "doi": ""},
    22:
        {"autor": "Avci, A., Karagoz, I.",
         "title": "A Novel Explicit Equation for Friction Factor in Smooth and"
                  "Rough Pipes",
         "ref": "J. Fluids Eng 131(6) (2009) 061203",
         "doi": "10.1115/1.3129132"},
    23:
        {"autor": "Papaevangelou, G., Evangelides, C., Tzimopoulos, C.,",
         "title": "A new explicit relation for friction coefficient in the "
                  "Darcy-Weisbach equation",
         "ref": "Proceedings of the Tenth Conference on Protection and "
                "Restoration of the Environment 166,1-7pp, PRE10 July 6-09 "
                "2010 Corfu, Greece.",
         "doi": ""},
    24:
        {"autor": "Brkić, D.",
         "title": "An Explicit Approximation of Colebrook’s equation for fluid"
                  "flow friction factor",
         "ref": "Petroleum Science and Technology 29 (15): 1596–1602. ",
         "doi": "10.1080/10916461003620453"},
    25:
        {"autor": "Fang, X,, Xu, Y., Zhou, Z.",
         "title": "New correlations of single-phase friction factor for"
                  "turbulent pipe flow and evaluation of existing single-phase"
                  "friction factor correlations.",
         "ref": "Nucl. Eng. Des. 241 (2011) 897-902",
         "doi": "10.1016/j.nucengdes.2010.12.019"},
    26:
        {"autor": "Ghanbari, A., Farshad, F., Rieke, H.H.",
         "title": "Newly developed friction factor correlation for pipe flow "
                  "and flow assurance",
         "ref": "J Chem Eng Mat Sci 2 (2011), 83-86.",
         "doi": ""},
    27:
        {"autor": "Goudar, C.T., Sonnad J.R.",
         "title": "Comparison of the iterative approximations of the "
                  "Colebrook-White equation",
         "ref": "Hydrocarb. Process. 87 (2008) 79-83",
         "doi": ""},
    28:
        {"autor": "Samadianfard, S.",
         "title": "Gene expression programming analysis of implicit "
                  "Colebrook-White equation in turbulent flow friction factor "
                  "calculation",
         "ref": "J. Pet. Sci. Eng. 92-93 (2012) 48-55",
         "doi": "10.1016/j.petrol.2012.06.005"},
    29:
        {"autor": "Darby, R., Chhabra, R.P.",
         "title": "Chemical Engineering Fluid Mechanics, 3rd Edition",
         "ref": "CRC Press, 2017",
         "doi": ""},

    # 30:
    #     {"autor": "",
    #      "title": "",
    #      "ref": "",
    #      "doi": ""},
}


# Friction factor for pipes
# All this function return the darcy friction factor, fanning friction factor
# can be obtained divided darcy factor by 4.
@refDoc(__doi__, [1])
def f_colebrook(Re, eD):
    r"""Calculates friction factor `f` with Colebrook-White correlation (1939)

    .. math::
        \frac{1}{\sqrt{f}}=-2\log\left(\frac{\epsilon/D}{3.7}+
        \frac{2.51}{Re\sqrt{f}}\right)

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]

    Notes
    -----
    This is the original, implicit expression, slowlest to solve
    """
    fo = f_chen(Re, eD)
    if eD:
        f = fsolve(lambda x: 1/x**0.5+2.0*log10(eD/3.7+2.51/Re/x**0.5), fo)
    else:
        f = fsolve(lambda x: 1/x**0.5-2.0*log10(Re*x**0.5)+0.8, fo)
    return Dimensionless(f[0])


@refDoc(__doi__, [2])
def f_chen1979(Re, eD):
    r"""Calculates friction factor `f` with Chen correlation (1979)

    .. math::
        \frac{1}{\sqrt{f}}=-2\log\left(\frac{\epsilon/D}{3.7065}
        -\frac{5.0452}{Re}\log\left(\frac{\left(\epsilon/D\right)
        ^{1.1098}}{2.8257}+\left(\frac{7.149}{Re}\right)^{0.8981}\right)\right)

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]

    Notes
    -----
    Range of validity:

        * 4e3 <= Re <= 4e8
        * 1e-7 <= eD <= 0.05
    """
    # Eq 7
    A = eD**1.1098/2.8257+5.8506/Re**0.8981
    f = 1/(-2*log10(eD/3.7065-5.0452/Re*log10(A)))**2
    return Dimensionless(f)


@refDoc(__doi__, [3])
def f_chen(Re, eD):
    r"""Calculates friction factor `f` with Chen correlation (1987)

    .. math::
        \frac{1}{\sqrt{f}}=-4\log\left(\frac{\epsilon/D}{3.7}-
        \frac{5.02}{Re}\log\left(\frac{\epsilon/D}{3.7}+
        \left(\frac{6.7}{Re}\right)^{0.9}\right)\right)

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]

    Notes
    -----
    The most satisfactory explicit friction factor correlation by [2].
    """
    A = eD/3.7+(6.7/Re)**0.9
    f = 1/(-2*log10(eD/3.7-5.02/Re*log10(A)))**2
    return Dimensionless(f)


@refDoc(__doi__, [4])
def f_moody(Re, eD):
    r"""Calculates friction factor `f` with Moody correlation (1947)

    .. math::
        f = 5.5 10^{-3}\left[1+\left(2 10^4\frac{\epsilon}{D} +
        \frac{10^6}{Re}\right)^{1/3}\right]

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]

    Notes
    -----
    Range of validity:

        * 4e3 <= Re <= 1e8
        * 0<= eD < 0.01.
    """
    f = 5.5e-3*(1+(2e4*eD+1e6/Re)**(1./3))
    return Dimensionless(f)


@refDoc(__doi__, [5])
def f_churchill(Re, eD):
    r"""Calculates friction factor `f` with Churchill correlation (1977)

    .. math::
        f = 2\left[(\frac{8}{Re})^{12} + (A + B)^{-1.5}\right]^{1/12}

    .. math::
        A = \left\{2.457\ln\left[(\frac{7}{Re})^{0.9}
        + 0.27\frac{\epsilon}{D}\right]\right\}^{16}

    .. math::
        B = \left( \frac{37530}{Re}\right)^{16}

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]

    Notes
    -----
    Represent fanning friction factor over the entire range of Reynolds numbers
    including intermediate region between laminar and turbulent flow.
    """
    # Eq 18
    A = (2.457*log(1/(0.27*eD+(7./Re)**0.9)))**16
    B = (37530./Re)**16
    f = 8.*((8./Re)**12+(A+B)**-1.5)**(1./12)
    return Dimensionless(f)


@refDoc(__doi__, [6])
def f_wood(Re, eD):
    r"""Calculates friction factor `f` with Wood correlation (1966)

    .. math::
        f_d = 0.094\left(\epsilon/D\right)^{0.225} + 0.53\left(
        \epsilon/D\right) + 88\left(\epsilon/D\right)^{0.4}Re^{-A_1}

    .. math::
        A_1 = 1.62\left(\epsilon/D\right)^{0.134}

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]

    Notes
    -----
    Range of validity:

        * Re > 10000
        * 1e-5 < eD < 0.04
    """
    a = 0.094*eD**0.225+0.53*eD
    b = 88*eD**0.44
    c = -1.62*eD**0.134
    f = a + b*Re**c
    return Dimensionless(f)


@refDoc(__doi__, [7])
def f_haaland(Re, eD):
    r"""Calculates friction factor `f` with Haaland correlation (1983)

    .. math::
        f = \left(-1.8\log_{10}\left[\left(\frac{\epsilon/D}{3.7}
        \right)^{1.11} + \frac{6.9}{Re}\right]\right)^{-2}

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]

    Notes
    -----
    Range of validity:

        * 4e3 <= Re <= 1e8
        * 1e-6 <= eD <= 0.05
    """
    # Eq 8
    f = 1/(-1.8*log10((eD/3.75)**1.11+6.9/Re))**2
    return Dimensionless(f)


@refDoc(__doi__, [8])
def f_serghides(Re, eD):
    r"""Calculates friction factor `f` with Serguides correlation (1984)

    .. math::
        f=\left[A-\frac{(B-A)^2}{C-2B+A}\right]^{-2}

    .. math::
        A=-2\log_{10}\left[\frac{\epsilon/D}{3.7}+\frac{12}{Re}\right]

    .. math::
        B=-2\log_{10}\left[\frac{\epsilon/D}{3.7}+\frac{2.51A}{Re}\right]

    .. math::
        C=-2\log_{10}\left[\frac{\epsilon/D}{3.7}+\frac{2.51B}{Re}\right]

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    A = -2*log10(eD/3.7+12/Re)
    B = -2*log10(eD/3.7+2.51*A/Re)
    C = -2*log10(eD/3.7+2.51*B/Re)
    f = (A-(B-A)**2/(C-2*B+A))**-2
    return Dimensionless(f)


@refDoc(__doi__, [9])
def f_round(Re, eD):
    r"""Calculates friction factor `f` with Round correlation (1980)

    .. math::
        \frac{1}{\sqrt{f}} = 1.8\log\left[\frac{Re}{0.135Re
        \epsilon/D+6.5}\right]

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]

    Notes
    -----
    Range of validity:

        * 4e3 <= Re <= 4e8
        * eD <= 0.05
    """
    # Eq 8
    f = 1/(1.8*log10(Re/(0.135*Re*eD+6.5)))**2
    return Dimensionless(f)


@refDoc(__doi__, [10])
def f_swamee(Re, eD):
    r"""Calculates friction factor `f` with Swamee-Jain correlation (1976)

    .. math::
        \frac{1}{\sqrt{f}} = -2\log\left[\left(\frac{6.97}{Re}\right)^{0.9}
        + (\frac{\epsilon}{3.7D})\right]

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]

    Notes
    -----
    Range of validity:

        * 5e3 <= Re <= 1e8
        * 1e-6 <= eD <= 5e-2
    """
    f = 1/(-2*log10(eD/3.7+(6.97/Re)**0.9))**2
    return Dimensionless(f)


@refDoc(__doi__, [11])
def f_jain(Re, eD):
    r"""Calculates friction factor `f` with Jain correlation (1976)

    .. math::
        \frac{1}{\sqrt{f}} = 1.14 - 2\log\left[\epsilon/D +
        \left(\frac{29.843}{Re}\right)^{0.9}\right]

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]

    Notes
    -----
    Range of validity:

        * 5e3 <= Re <= 1e7
        * 4e-5 <= eD <= 0.05
    """
    f = 1/(1.14-2*log10(eD+(29.843/Re)**0.9))**2
    return Dimensionless(f)


@refDoc(__doi__, [12])
def f_barr(Re, eD):
    r"""Calculates friction factor `f` with Barr correlation (1981)

    .. math::
        \frac{1}{\sqrt{f}} = -2\log\left\{\frac{\epsilon}{3.7D} +
        \frac{4.518\log(\frac{Re}{7})}{Re\left[1+\frac{Re^{0.52}}{29}
        \left(\frac{\epsilon}{D}\right)^{0.7}\right]}\right\}

    Parameters
    ----------
        Re : float
            Reynolds number, [-]
        eD : float
            Relative roughness of a pipe, [-]

    Returns
    -------
        f : float
            Friction factor, [-]
    """
    # Eq 6
    f = 1/(2*log10(eD/3.7+4.518*log10(Re/7)/Re/(1+Re**0.52/29*eD**0.7)))**2
    return Dimensionless(f)


@refDoc(__doi__, [13])
def f_zigrang(Re, eD):
    r"""Calculates friction factor `f` with Zigrang-Sylvester correlation (1982)

    .. math::
        \frac{1}{\sqrt{f}} = -2\log\left[\frac{\epsilon}{3.7D}
        - \frac{5.02}{Re}\log A\right]

    .. math::
        A = \frac{\epsilon}{3.7D} + \frac{13}{Re}

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]

    Notes
    -----
    Range of validity:

        * 4e3 <= Re <= 1e8
        * 4e-5 <= eD <= 5e-2
    """
    # Eq 12
    A = log10(eD/3.7-5.02/Re*log10(eD/3.7+13./Re))
    f = 1/(-2*log10(eD/3.7-5.02*A/Re))**2
    return Dimensionless(f)


@refDoc(__doi__, [14])
def f_altshul(Re, eD):
    r"""Calculates friction factor `f` with Altshul correlation (1975)

    .. math::
        f = 0.11\left( \frac{68}{Re} + \frac{\epsilon}{D}\right)^{0.25}

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    f = 0.11*(eD+68/Re)**0.25
    return Dimensionless(f)


@refDoc(__doi__, [14])
def f_tsal(Re, eD):
    r"""Calculates friction factor `f` with Tsal correlation (1989)

    .. math::
        A = 0.11(\frac{68}{Re} + \frac{\epsilon}{D})^{0.25}

    if A >= 0.018 then f = A

    if A < 0.018 then f = 0.0028 + 0.85 A

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]

    Notes
    -----
    Range of validity:

        * 4e3 <= Re <= 1e8
        * eD <= 0.05
    """
    f = 0.11*(68/Re+eD)**0.25
    if f < 0.018:
        f = 0.0028+0.85*f
    return f


@refDoc(__doi__, [15])
def f_eck(Re, eD):
    r"""Calculates friction factor `f` with Eck correlation (1973)

    .. math::
        \frac{1}{\sqrt{f_d}} = -2\log\left[\frac{\epsilon}{3.715D}
        + \frac{15}{Re}\right]

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    f = 1/(-2*log10(eD/3.71+15/Re))**2
    return Dimensionless(f)


@refDoc(__doi__, [16])
def f_shacham(Re, eD):
    r"""Calculates friction factor `f` with Shacham correlation (1980)

    .. math::
        \frac{1}{\sqrt{f}} = -2\log\left[\frac{\epsilon}{3.7D} -
        \frac{5.02}{Re} \log\left(\frac{\epsilon}{3.7D}
        + \frac{14.5}{Re}\right)\right]

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]

    Notes
    -----
    Range of validity:

        * 4e3 <= Re <= 4e8
    """
    f = 1/(-2*log10(eD/3.7-5.02/Re*log10(eD/3.7+14.5/Re)))**2
    return Dimensionless(f)


@refDoc(__doi__, [17])
def f_manadilli(Re, eD):
    r"""Calculates friction factor `f` with Manadilli correlation (1997)

    .. math::
        \frac{1}{\sqrt{f}} = -2\log\left[\frac{\epsilon}{3.7D} +
        \frac{95}{Re^{0.983}} - \frac{96.82}{Re}\right]

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]

    Notes
    -----
    Range of validity:

        * 5.245e3 <= Re <= 1e8
        * eD <= 0.05
    """
    f = 1/(-2*log10(eD/3.7+95./Re**0.983-96.82/Re))**2
    return Dimensionless(f)


@refDoc(__doi__, [18])
def f_romeo(Re, eD):
    r"""Calculates friction factor `f` with Monzon-Romeo-Royo correlation (2002)

    .. math::
        \frac{1}{\sqrt{f_d}} = -2\log\left\{\frac{\epsilon}{3.7065D}\times
        \frac{5.0272}{Re}\times\log\left[\frac{\epsilon}{3.827D} -
        \frac{4.567}{Re}\times\log\left(\frac{\epsilon}{7.7918D}^{0.9924} +
        \left(\frac{5.3326}{208.815+Re}\right)^{0.9345}\right)\right]\right\}

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]

    Notes
    -----
    Range of validity:

        * 3e3 <= Re <= 1.5e8
        * eD <= 0.05
    """
    A = log10((eD/7.7918)**0.9924+(5.3326/(208.815+Re))**0.9345)
    B = log10(eD/3.827-4.567/Re*A)
    f = 1/(-2*log10(eD/3.7065-5.0272*B/Re))**2
    return Dimensionless(f)


@refDoc(__doi__, [19])
def f_goudar2007(Re, eD):
    r"""Calculates friction factor `f` with Goudar-Sonnad correlation (2007)

    .. math::
        \frac{1}{\sqrt{f}} = 0.8686\ln\left(\frac{0.4587Re}{S^{S/(S+1)}}\right)

    .. math::
        S = 0.1240\times\frac{\epsilon}{D}\times Re + \ln(0.4587Re)

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]

    Notes
    -----
    Range of validity:

        * 4e3 <= Re <= 1e8
        * 1e-6 <= eD <= 0.05
    """
    C = 0.124*Re*eD+log(0.4587*Re)
    f = 1/(0.8686*log(0.4587*Re/(C-0.31)**(C/(C+1))))**2
    return Dimensionless(f)


@refDoc(__doi__, [27])
def f_goudar(Re, eD):
    r"""Calculates friction factor `f` with Goudar-Sonnad correlation (2008)

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    a = 2/log(10)
    b = eD/3.7
    d = log(10)*Re/5.2
    s = b*d+log(d)
    q = s**(s/(s+1))
    g = b*d+log(d/q)
    z = log(q/g)
    Dla = g/(g+1)*z
    Dcfa = Dla*(1+z/2/((g+1)**2+z/3*(2*g-1)))

    f = (a*(log(d/q)+Dcfa))**-2
    return Dimensionless(f)


@refDoc(__doi__, [20])
def f_buzzelli(Re, eD):
    r"""Calculates friction factor `f` with Buzzelli correlation (2008)

    .. math::
        \frac{1}{\sqrt{f}} = A - \left[\frac{A +2\log(\frac{B}{Re})}
        {1 + \frac{2.18}{B}}\right]

    .. math::
        A = \frac{0.774\ln(Re)-1.41}{1+1.32\sqrt{\frac{\epsilon}{D}}}

    .. math::
        B = \frac{\epsilon}{3.7D}Re+2.51\times B_1

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    # Eq 2
    A = (0.744*log(Re)-1.41)/(1+1.32*eD**0.5)
    B = eD/3.7*Re+2.51*A
    f = 1/(A-((A+2*log10(B/Re))/(1+2.18/B)))**2
    return Dimensionless(f)


@refDoc(__doi__, [21])
def f_Vatankhah(Re, eD):
    r"""Calculates friction factor `f` with Vatankhah-Kouchakzadeh corr (2008)

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    S = 0.124*Re*eD+log(0.4587*Re)
    f = 1/(0.8686*log(0.4587*Re/(S-0.31)**(S/(S+0.9633))))**2
    return Dimensionless(f)


@refDoc(__doi__, [22])
def f_avci(Re, eD):
    r"""Calculates friction factor `f` with Avci-Karagoz correlation (2009)

    .. math::
        f = \frac{6.4} {\left\{\ln(Re) - \ln\left[
        1 + 0.01Re\frac{\epsilon}{D}\left(1 + 10(\frac{\epsilon}{D})^{0.5}
        \right)\right]\right\}^{2.4}}

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    # Eq 17
    f = 6.4/(log(Re)-log(1+0.01*Re*eD*(1+10*eD**0.5)))**2.4
    return Dimensionless(f)


@refDoc(__doi__, [23])
def f_papaevangelou(Re, eD):
    r"""Calculates friction factor `f` with Papaevangelou correlation (2009)

    .. math::
        f = \frac{0.2479 - 0.0000947(7-\log Re)^4}{\left[\log\left
        (\frac{\epsilon}{3.615D} + \frac{7.366}{Re^{0.9142}}\right)\right]^2}

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]

    Notes
    -----
    Range of validity:

        * 1e4 <= Re <= 1e7
        * 1e-5 <= eD <= 1e-3
    """
    # Eq 12
    f = (0.2479-9.47e-5*(7-log10(Re))**4)/log10(eD/3.615+7.366/Re**0.9142)**2
    return Dimensionless(f)


@refDoc(__doi__, [24])
def f_brkic(Re, eD, alternate=False):
    r"""Calculates friction factor `f` with Brkić correlation (2010)

    .. math::
        f = [-2\log(10^{-0.4343\beta} + \frac{\epsilon}{3.71D})]^{-2}

    .. math::
        \beta = \ln \frac{Re}{1.816\ln\left(\frac{1.1Re}{\ln(1+1.1Re)}\right)}

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]
    alternate : boolean
        Choose the alternate correlation from the paper

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    S = log(Re/1.816/log(1.1*Re/log(1+1.1*Re)))                          # Eq 8

    # Eq 7
    if alternate:
        f = 1/(-2*log10(2.18*S/Re+eD/3.71))**2
    else:
        f = 1/(-2*log10(10**(-0.4343*S)+eD/3.71))**2

    return Dimensionless(f)


@refDoc(__doi__, [25])
def f_fang(Re, eD):
    r"""Calculates friction factor `f` with Fang-Xua-Zhou correlation (2011)

    .. math::
        f = 1.613\left\{\ln\left[0.234\frac{\epsilon}{D}^{1.1007}
        - \frac{60.525}{Re^{1.1105}}
        + \frac{56.291}{Re^{1.0712}}\right]\right\}^{-2}

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]

    Notes
    -----
    Range of validity:

        * 3e3 <= Re <= 1e8
        * eD <= 0.05
    """
    # Eq 13
    f = 1.613/(log(0.234*eD**1.1007-60.525/Re**1.1105+56.291/Re**1.0712))**2
    return Dimensionless(f)


@refDoc(__doi__, [26])
def f_ghanbari(Re, eD):
    r"""Calculates friction factor `f` with Ghanbari correlation (2011)
    Second correlation

    .. math::
        f = \left(-1.52\log\left(\left(\frac{\epsilon/D}{7.21}\right)^{1.042}
        + \left(\frac{2.731}{Re}\right)^{0.9152}\right)\right)^{-2.169}

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    # Eq 10
    return (-1.52*log10((eD/7.21)**1.042+(2.731/Re)**0.9152))**-2.169


@refDoc(__doi__, [28])
def f_Samadianfard(Re, eD):
    r"""Calculates friction factor `f` with Samadianfard correlation (2012)

    .. math::
        f = \frac{Re^{\epsilon/D}-0.6315093}{Re^{1/3}+Re·\epsilon/D}
        + 0.0275308\left(\frac{6.929841}{Re}+\epsilon/D\right)^{1/9}
        + \frac{10^{\epsilon/D}}{\epsilon/D+4.781616} \left(\epsilon/D
        +\frac{9.99701}{Re}\right)

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]
    """
    # Eq 29
    f = (Re**eD-0.6315093)/(Re**(1/3)+Re*eD) + \
        0.0275308*(6.929841/Re+eD)**(1/9) + \
        10**eD/(eD+4.481616)*(eD**0.5+9.99701/Re)
    return f


f_list = (f_colebrook, f_chen, f_Vatankhah, f_buzzelli, f_romeo, f_serghides,
          f_zigrang, f_Samadianfard, f_brkic, f_fang, f_ghanbari, f_haaland,
          f_round, f_swamee, f_jain, f_barr, f_shacham, f_tsal, f_manadilli,
          f_goudar, f_goudar2007, f_avci, f_papaevangelou, f_churchill,
          f_chen1979, f_moody, f_wood, f_eck, f_altshul)


def f_blasius(Re):
    """
    Friction factor by Blasius, for bare tubes

    Input parameters:
    Re: Reynolds number
    """
    return 0.079/Re**0.25


def f_Gnielinsky(Re):
    """
    Friction factor by Gnielinsky, for annulli section in bare tubes

    Input parameters:
    Re: Reynolds number, based in hydraulic number
    """
    return (1.8*log(Re)-1.5)**-2.


@refDoc(__doi__, [29])
def f_friccion(Re, eD=0, method=0, geometry=0, *args):
    """
    Generalized method for calculate friction factor for laminar or turbulent
    flux in several geometries

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]
    method: int
        Index of method to use (default 0 for use Colebrook original function):

            * 0 - Colebrook (default)
            * 1 - Chen (1987)
            * 2 - Vatankhah-Kouchakzadeh (2008)
            * 3 - Buzzelli (2008)
            * 4 - Romeo (2002)
            * 5 - Serghides (1984)
            * 6 - Zigrang-Sylvester (1982)
            * 7 - Samadianfard (2012)
            * 8 - Brkić (2011)
            * 9 - Fang-Xu-Zhou (2011)
            * 10 - Ghanbari (2011)
            * 11 - Haaland (1983)
            * 12 - Round (1980)
            * 13 - Swamee-Jain (1976)
            * 14 - Jain (1976)
            * 15 - Barr (1981)
            * 16 - Shacham (1980)
            * 17 - Tsal (1989)`
            * 18 - Manadilli (1997)
            * 19 - Goudar (2008)
            * 20 - Goudar (2007)
            * 21 - Avci (2009)
            * 22 - Papaevangelou (2010)
            * 23 - Churchill (1977)
            * 24 - Chen (1979)
            * 25 - Moody (1947)
            * 26 - Wood (1966)
            * 27 - Eck (1973)
            * 28 - Altshul (1975)

    geometry: int
        Index for duct geometry (Table 7.1, Pag 184 in [29]_)

            * 0 - Circular section (default)
            * 1 - Square
            * 2 - Isosceles triangle
            * 3 - Rectangle
            * 4 - Ellipse
            * 5 - Right triangle
            * 6 - Anulli
    args: float
        Other parameter necessary for noncircular geometries

            * Isosceles triangle: Angle, [º]
            * Rectangle: Ratio width/height, [-]
            * Ellipse: Both diameters of ellipse, [m]
            * Right triangle: Angle, [º]
            * Anulli: Internal and external diameter, [m]
    """
    if Re < 2100:
        if geometry == 0:
            # Circle
            f_friccion = 16./Re
        elif geometry == 1:
            # Square
            f_friccion = 14.2/Re
        elif geometry == 2:
            # Isosceles triangle
            pass
        elif geometry == 3:
            # Rectangle
            D, d = args[1], args[0]
            f_friccion = 16/(2/3+11/24*d/D*(2-d/D))/Re
        elif geometry == 4:
            # Ellipse
            D, d = args[1], args[0]
            c = (D-d)/(D+d)
            Dh = 4*d*D*(64-16*c**2)/((d+D)*(64-3*c**4))
            f_friccion = 2*Dh**2*(D**2+d**2)/D**2/d**2/Re
        elif geometry == 5:
            pass
        elif geometry == 6:
            # Annulli
            # Eq 7.9, 7.10, pag 183
            Di, Do = args
            alpha = (Do-Di)**2/(Do**2+Di**2-(Do**2-Di**2)/log(Do/Di))
            f_friccion = 16*alpha/Re

    else:
        if geometry == 6:
            f_friccion = f_Gnielinsky(Re)
        else:
            f_friccion = f_list[method](Re, eD)

    return Dimensionless(f_friccion)


@refDoc(__doi__, [1])
def eD(Re, f):
    """Calculates relative roughness

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    f : float
        Friction factor, [-]

    Returns
    -------
    eD : float
        Relative roughness of a pipe, [-]
    """
    eD = (10**(-0.5/f**0.5)-2.51/Re/f**0.5)*3.7
    return eD


if __name__ == "__main__":
    for f in f_list:
        line = f.__doc__.split("\n")[0]
        year = line.split(" ")[-1]
        name = line.split(" ")[-3]
        doc = name + " " + year
        print(f(1e7, 0.0002), doc)
