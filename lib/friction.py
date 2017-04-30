#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2017, Juan José Gómez Romera <jjgomera@gmail.com>

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


###############################################################################
# Module for implement friction factor
#   -Friction factor for rough pipes
#       · f_colebrook
#       · f_chen (1991)
#       · f_chen (1979)
#       · f_moody (1947)
#       · f_churchill (1977)
#       · f_wood (1966)
#       · f_eck (1973)
#       · f_altshul (1975)
#       · f_haaland (1983)
#       · f_serghides (1984)
#       · f_round (1980)
#       · f_swamee (1976)
#       · f_jain (1976)
#       · f_barr (1981)
#       · f_zigrang (1982)
#       · f_shacham (1980)
#       · f_tsal (1989)
#       · f_manadilli (1997)
#       · f_romeo (2002)
#       · f_goudar (2008)
#       · f_goudar2006 (2006)
#       · f_buzzelli (2008)
#       · f_Vatankhah (2008)
#       · f_avci (2009)
#       · f_papaevangelou (2009)
#       · f_brkic (2010)
#       · f_brkic2 (2010)
#       · f_fang (2011)
#       · f_ghanbari (2011)
#
#
#
#   -Fitting K
###############################################################################

from math import exp, log, log10, sqrt, sin, pi

from scipy.optimize import fsolve

from lib.unidades import Dimensionless


__doi__ = [
    {"autor": "Colebrook, C. F. and White, C. M.",
     "title": "Experiments with Fluid Friction in Roughened Pipes",
     "ref": "Proceedings of the Royal Society of London. Series A,"
            "Mathematical and Physical Sciences 161 (906): 367–381.",
     "doi": "10.1098/rspa.1937.0150"},
    {"autor": "Chen, H.J.",
     "title": "An Explicit Equation for Friction Factor in Pipe",
     "ref": "Ind. Eng. Chem. Fundamen., 1979, 18 (3), pp 296–297",
     "doi": "10.1021/i160071a019"},
    {"autor": "Haaland, S. E.",
     "title": "Simple and explicit formulas for the friction factor in"
              "turbulent flow",
     "ref": "J. Fluids Eng., 105(1), 89–90.",
     "doi": "10.1115/1.3240948"},
    {"autor": "Round, G. F.",
     "title": "An explicit approximation for the friction factor-reynolds"
              "number relation for rough and smooth pipes",
     "ref": "Can. J. Chem. Eng., 58: 122-123",
     "doi": "10.1002/cjce.5450580119"},
    {"autor": "Zigrang, D.J., Sylvester, N.D.",
     "title": "Explicit approximations to the solution of Colebrook's"
              "friction factor equation",
     "ref": "AICHE J 28, 514-515.",
     "doi": "10.1002/aic.690280323"},
    {"autor": "Avci, A.; Karagoz, I.",
     "title": "A Novel Explicit Equation for Friction Factor in Smooth and"
              "Rough Pipes",
     "ref": "J. Fluids Eng 131(6), 061203 (May 13, 2009)",
     "doi": "10.1115/1.3129132"},
    {"autor": "Brkić, D.",
     "title": "An Explicit Approximation of Colebrook’s equation for fluid"
              "flow friction factor",
     "ref": "Petroleum Science and Technology 29 (15): 1596–1602. ",
     "doi": "10.1080/10916461003620453"},
    {"autor": "Romeo, E., Royo, C., Monzon, A.",
     "title": "Improved explicit equation for estimation of the friction"
              "factor in rough and smooth pipes.",
     "ref": "Chem. Eng. J. 86 (3), 369–374. (2002)",
     "doi": "10.1016/S1385-8947(01)00254-6"},
    {"autor": "Fang X, Xua Y, Zhou Z",
     "title": "New correlations of single-phase friction factor for"
              "turbulent pipe flow and evaluation of existing single-phase"
              "friction factor correlations.",
     "ref": "Nucl Eng Des 241, 897-902. (2011)",
     "doi": "10.1016/j.nucengdes.2010.12.019"}]


# Friction factor for pipes
# All this function return the darcy friction factor, fanning friction factor
# can be obtained divided darcy factor by 4.
def f_colebrook(Re, eD):
    """
    Calculates friction factor `f` with Colebrook-White correlation (1939)
    This is the original impicit correlation

    .. math::
        $\frac{1}{\sqrt{f}}=-2\log\left(\frac{\nicefrac{\epsilon}{D}}{3.7}+
        \frac{2.51}{Re\sqrt{f}}\right)$

    Parameters
    ------------
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

    References
    ----------
    [1] .. Colebrook, C. F. and White, C. M. (1937). "Experiments with Fluid
    Friction in Roughened Pipes". Proceedings of the Royal Society of London.
    Series A, Mathematical and Physical Sciences 161 (906): 367–381.
    """
    fo = f_chen(Re, eD)
    if eD:
        f = fsolve(lambda x: 1/x**0.5+2.0*log10(eD/3.7+2.51/Re/x**0.5), fo)
    else:
        f = fsolve(lambda x: 1/x**0.5-2.0*log10(Re*x**0.5)+0.8, fo)
    return f[0]


def f_chen(Re, eD):
    """
    Calculates friction factor `f` with Chen correlation (1991)

    .. math::
        $\frac{1}{\sqrt{f}}=-4\log\left(\frac{\nicefrac{\epsilon}{D}}{3.7}-\frac{5.02}{Re}\log\left(\frac{\nicefrac{\epsilon}{D}}{3.7}+\left(\frac{6.7}{Re}\right)^{0.9}\right)\right)$

    Parameters
    ------------
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

    References
    ----------
    [1] .. Chen, H-J, “An Exact Solution to the Colebrook Equation”, Chem.
    Eng., Vol. 94, No. 2, 1987, p. 196..
    [2] .. Gregory, G.A. and M. Fogarasi, “Alternate to Standard Friction
    Factor Equation”, Oil Gas J., Apr 1, 1985, p. 120.
    """
    A = eD/3.7+(6.7/Re)**0.9
    f = 1/(-2*log10(eD/3.7-5.02/Re*log10(A)))**2
    return f


def f_chen1979(Re, eD):
    """
    Calculates friction factor `f` with Chen correlation (1979)

    .. math::
        $\frac{1}{\sqrt{f}}=-2\log\left(\frac{\nicefrac{\epsilon}{D}}{3.7065}-\frac{5.0452}{Re}\log\left(\frac{\left(\nicefrac{\epsilon}{D}\right)^{1.1098}}{2.8257}+\left(\frac{7.149}{Re}\right)^{0.8981}\right)\right)$

    Parameters
    ------------
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
    Range of validity
        4e3 <= Re <= 4e8
        1e-7 <= eD <= 0.05

    References
    ----------
    [1] .. Chen, H-J, “An Explicit Equation for Friction Factor in Pipe”, Ind.
    Eng. Chem. Fundam., Vol. 18, No. 3, 1979, p. 297..
    """
    A = eD**1.1098/2.8257+5.8506/Re**0.8981
    f = 1/(-2*log10(eD/3.7065-5.0452/Re*log10(A)))**2
    return f


def f_moody(Re, eD):
    """
    Calculates friction factor `f` with Moody correlation (1947)

    .. math::
        f = 5.5 10^{-3}\left[1+\left(2 10^4\frac{\epsilon}{D} +
        \frac{10^6}{Re}\right)^{1/3}\right]

    Parameters
    ------------
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
        4e3 <= Re <= 1e8
        0<= eD < 0.01.

    References
    ----------
    [1] .. Moody, L. F. (1947). “An approximate formula for pipe friction
    factors.” Trans. ASME, 69(12), 1005–1006.
    """
    return 5.5e-3*(1+(2e4*eD+1e6/Re)**(1./3))


def f_churchill(Re, eD):
    """
    Calculates friction factor `f` with Churchill correlation (1977)

    .. math::
        f = 2\left[(\frac{8}{Re})^{12} + (A_2 + A_3)^{-1.5}\right]^{1/12}

        A = \left\{2.457\ln\left[(\frac{7}{Re})^{0.9}
        + 0.27\frac{\epsilon}{D}\right]\right\}^{16}

        B = \left( \frac{37530}{Re}\right)^{16}

    Parameters
    ------------
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

    References
    ----------
    [1] .. Churchill, S. W. (1977). “Friction factor equation that spans all
    fluid-flow regimes.” Chem. Eng., 84, 91–92.
    """
    A = (2.457*log(1/(0.27*eD+(7./Re)**0.9)))**16
    B = (37530./Re)**16
    return 8.*((8./Re)**12+(A+B)**-1.5)**(1./12)


def f_wood(Re, eD):
    """
    Calculates friction factor `f` with Wood correlation (1966)

    .. math::
        f_d = 0.094(\frac{\epsilon}{D})^{0.225} + 0.53(\frac{\epsilon}{D})
        + 88(\frac{\epsilon}{D})^{0.4}Re^{-A_1}

        A_1 = 1.62(\frac{\epsilon}{D})^{0.134}

    Parameters
    ------------
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
    Range of validity
        Re > 10000
        1e-5 < eD < 0.04

    References
    ----------
    [1] .. Wood DJ (1966). An explicit friction factor relationship. Civil Eng.
    ASCE 60
    """
    a = -1.62*eD**0.134
    return 0.094*eD**0.225+0.53*eD+88*eD**0.44*Re**a


def f_haaland(Re, eD):
    """
    Calculates friction factor `f` with Haaland correlation (1983)

    .. math::
        f = \left(-1.8\log_{10}\left[\left(\frac{\epsilon/D}{3.7}
        \right)^{1.11} + \frac{6.9}{Re}\right]\right)^{-2}

    Parameters
    ------------
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
    Range of validity
        4e3 <= Re <= 1e8
        1e-6 <= eD <= 0.05

    References
    ----------
    [1] .. Haaland, S. E. (1983). “Simple and explicit formulas for the
    friction factor in turbulent flow.” J. Fluids Eng., 105(1), 89–90.
    """
    return 1/(-1.8*log10((eD/3.7)**1.11+6.9/Re))**2


def f_serghides(Re, eD):
    """
    Calculates friction factor `f` with Serguides correlation (1984)

    .. math::
        f=\left[A-\frac{(B-A)^2}{C-2B+A}\right]^{-2}

        A=-2\log_{10}\left[\frac{\epsilon/D}{3.7}+\frac{12}{Re}\right]

        B=-2\log_{10}\left[\frac{\epsilon/D}{3.7}+\frac{2.51A}{Re}\right]

        C=-2\log_{10}\left[\frac{\epsilon/D}{3.7}+\frac{2.51B}{Re}\right]

    Parameters
    ------------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]

    References
    ----------
    [1] .. Serghides, T. K. (1984). “Estimate friction factor accurately.”
    Chem. Eng., 91(5), 63–64.
    """
    A = -2*log10(eD/3.7+12/Re)
    B = -2*log10(eD/3.7+2.51*A/Re)
    C = -2*log10(eD/3.7+2.51*B/Re)
    return (A-(B-A)**2/(C-2*B+A))**-2


def f_round(Re, eD):
    """
    Calculates friction factor `f` with Round correlation (1980)

    .. math::
        \frac{1}{\sqrt{f}} = 1.8\log\left[\frac{Re}{0.135Re
        \frac{\epsilon}{D}+6.5}\right]

    Parameters
    ------------
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
    Range of validity
        4e3 <= Re <= 4e8
        eD <= 0.05

    References
    ----------
    [1] .. Round, G. F. (1980), An explicit approximation for the friction
    factor-reynolds number relation for rough and smooth pipes. Can. J. Chem.
    Eng., 58: 122–123.
    """
    return 1/(1.8*log10(Re/(0.135*Re*eD+6.5)))**2


def f_swamee(Re, eD):
    """
    Calculates friction factor `f` with Swamee-Jain correlation (1976)

    .. math::
        \frac{1}{\sqrt{f_f}} = -2\log\left[\left(\frac{6.97}{Re}\right)^{0.9}
        + (\frac{\epsilon}{3.7D})\right]

    Parameters
    ------------
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
    Range of validity
        5e3 <= Re <= 1e8
        1e-6 <= eD <= 5e-2

    References
    ----------
    [1] ..Swamee, P.K.; Jain, A.K. (1976). "Explicit equations for pipe-flow
    problems". Journal of the Hydraulics Division (ASCE) 102 (5): 657–664.
    """
    return 1/(-2*log10(eD/3.7+(6.97/Re)**0.9))**2


def f_jain(Re, eD):
    """
    Calculates friction factor `f` with Jain correlation (1976)

    .. math::
        \frac{1}{\sqrt{f_f}} = 1.14 - 2\log\left[ \frac{\epsilon}{D} +
        \left(\frac{29.843}{Re}\right)^{0.9}\right]

    Parameters
    ------------
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
        5e3 <= Re <= 1e7
        4e-5 <= eD <= 0.05

    References
    ----------
    [1] .. Jain, Akalank K."Accurate Explicit Equation for Friction Factor."
       Journal of the Hydraulics Division 102, no. 5 (May 1976): 674-77.
    """
    return 1/(1.14-2*log10(eD+(29.843/Re)**0.9))**2


def f_barr(Re, eD):
    """
    Calculates friction factor `f` with Barr correlation (1981)

    .. math::
        \frac{1}{\sqrt{f}} = -2\log\left\{\frac{\epsilon}{3.7D} +
        \frac{4.518\log(\frac{Re}{7})}{Re\left[1+\frac{Re^{0.52}}{29}
        \left(\frac{\epsilon}{D}\right)^{0.7}\right]}\right\}

    Parameters
    ------------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]

    References
    ----------
    [1] .. Barr DIH (1981). Solutions of the Colebrook-White functions for
    resistance to uniform turbulent flows.  Proc Inst Civil Eng 71, 529-536.
    """
    return 1/(2*log10(eD/3.7+4.518*log10(Re/7)/Re/(1+Re**0.52/29*eD**0.7)))**2


def f_zigrang(Re, eD):
    """
    Calculates friction factor `f` with Zigrang-Sylvester correlation (1982)

    .. math::
        \frac{1}{\sqrt{f}} = -2\log\left[\frac{\epsilon}{3.7D}
        - \frac{5.02}{Re}\log A\right]

        A = \frac{\epsilon}{3.7D} + \frac{13}{Re}

    Parameters
    ------------
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
    Range of validity
        4e3 <= Re <= 1e8
        4e-5 <= eD <= 5e-2

    References
    ----------
    [1] .. Zigrang DJ, Sylvester ND (1982). Explicit approximations to the
    Colebrook’s friction factor. AICHE J 28, 514-515.
    """
    A = log10(eD/3.7-5.02/Re*log10(eD/3.7+13./Re))
    return 1/(-2*log10(eD/3.7-5.02*A/Re))**2


def f_altshul(Re, eD):
    """
    Calculates friction factor `f` with Altshul correlation (1975)

    .. math::
        f = 0.11\left( \frac{68}{Re} + \frac{\epsilon}{D}\right)^{0.25}

    Parameters
    ------------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]

    References
    ----------
    [1] .. Tsal RJ (1989). Altshul-Tsal friction factor equation. Heating
    Piping Air Conditioning 8, 30-45.
    """
    return 0.11*(eD+68/Re)**0.25


def f_tsal(Re, eD):
    """
    Calculates friction factor `f` with Tsal correlation (1989)

    .. math::
        A = 0.11(\frac{68}{Re} + \frac{\epsilon}{D})^{0.25}

    if A >= 0.018 then f = A
    if A < 0.018 then f = 0.0028 + 0.85 A

    Parameters
    ------------
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
    Range of validity
        4e3 <= Re <= 1e8
        eD <= 0.05

    References
    ----------
    [1] .. Tsal RJ (1989). Altshul-Tsal friction factor equation. Heating
    Piping Air Conditioning 8, 30-45.
    """
    f = 0.11*(68/Re+eD)**0.25
    if f < 0.018:
        f = 0.0028+0.85*f
    return f


def f_eck(Re, eD):
    """
    Calculates friction factor `f` with Eck correlation (1973)

    .. math::
        \frac{1}{\sqrt{f_d}} = -2\log\left[\frac{\epsilon}{3.715D}
        + \frac{15}{Re}\right]

    Parameters
    ------------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]

    References
    ----------
    [1] .. Eck B (1973). Technische Stromungslehre. Springer, New York.
    """
    return 1/(-2*log10(eD/3.71+15/Re))**2


def f_shacham(Re, eD):
    """
    Calculates friction factor `f` with Shacham correlation (1980)

    .. math::
        \frac{1}{\sqrt{f}} = -2\log\left[\frac{\epsilon}{3.7D} -
        \frac{5.02}{Re} \log\left(\frac{\epsilon}{3.7D}
        + \frac{14.5}{Re}\right)\right]

    Parameters
    ------------
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
    Range of validity
        4e3 <= Re <= 4e8

    References
    ----------
    [1] .. Shacham M (1980). An explicit equation for friction factor in pipe.
    Ind Eng Chem Fund 19, 228 - 229.
    """
    return 1/(-2*log10(eD/3.7-5.02/Re*log10(eD/3.7+14.5/Re)))**2


def f_manadilli(Re, eD):
    """
    Calculates friction factor `f` with Manadilli correlation (1997)

    .. math::
        \frac{1}{\sqrt{f}} = -2\log\left[\frac{\epsilon}{3.7D} +
        \frac{95}{Re^{0.983}} - \frac{96.82}{Re}\right]

    Parameters
    ------------
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
    Range of validity
        5.245e3 <= Re <= 1e8
        eD <= 0.05

    References
    ----------
    [1] .. Manadilli G (1997). Replace implicit equations with signomial
    functions. Chem Eng 104, 129-132.
    """
    return 1/(-2*log10(eD/3.7+95./Re**0.983-96.82/Re))**2


def f_romeo(Re, eD):
    """
    Calculates friction factor `f` with Monzon-Romeo-Royo correlation (2002)

    .. math::
        \frac{1}{\sqrt{f_d}} = -2\log\left\{\frac{\epsilon}{3.7065D}\times
        \frac{5.0272}{Re}\times\log\left[\frac{\epsilon}{3.827D} -
        \frac{4.567}{Re}\times\log\left(\frac{\epsilon}{7.7918D}^{0.9924} +
        \left(\frac{5.3326}{208.815+Re}\right)^{0.9345}\right)\right]\right\}

    Parameters
    ------------
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
    Range of validity
        3e3 <= Re <= 1.5e8
        eD <= 0.05

    References
    ----------
    [1] .. Romeo, E., Royo, C., Monzon, A., 2002. Improved explicit equation
    for estimation of the friction factor in rough and smooth pipes. Chem. Eng.
    J. 86 (3), 369–374.
    """
    A = log10((eD/7.7918)**0.9924+(5.3326/(208.815+Re))**0.9345)
    B = log10(eD/3.827-4.567/Re*A)
    return 1/(-2*log10(eD/3.7065-5.0272*B/Re))**2


def f_goudar2006(Re, eD):
    """
    Calculates friction factor `f` with Goudar-Sonnad correlation (2006)

    .. math::
        \frac{1}{\sqrt{f}} = 0.8686\ln\left(\frac{0.4587Re}{S^{S/(S+1)}}\right)

        S = 0.1240\times\frac{\epsilon}{D}\times Re + \ln(0.4587Re)

    Parameters
    ------------
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
    Range of validity
        4e3 <= Re <= 1e8
        1e-6 <= eD <= 0.05

    References
    ----------
    [1] .. Goudar CT, Sonnad JR (2008). Comparison of the iterative
    approximations of the Colebrook-White equation. Hydrocarb Process 87, 79-83
    """
    C = 0.124*Re*eD+log(0.4587*Re)
    return 1/(0.8686*log(0.4587*Re/(C-0.31)**(C/(C+1))))**2


def f_goudar(Re, eD):
    """
    Calculates friction factor `f` with Goudar-Sonnad correlation (2008)

    Parameters
    ------------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]

    References
    ----------
    [1] .. Goudar CT, Sonnad JR (2008). Comparison of the iterative
    approximations of the Colebrook-White equation. Hydrocarb Process 87, 79-83
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

    return (a*(log(d/q)+Dcfa))**-2


def f_buzzelli(Re, eD):
    """
    Calculates friction factor `f` with Buzzelli correlation (2008)

    .. math::
        \frac{1}{\sqrt{f}} = A - \left[\frac{A +2\log(\frac{B}{Re})}
        {1 + \frac{2.18}{B}}\right]

        A = \frac{0.774\ln(Re)-1.41}{1+1.32\sqrt{\frac{\epsilon}{D}}}

        B = \frac{\epsilon}{3.7D}Re+2.51\times B_1

    Parameters
    ------------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]

    References
    ----------
    [1] .. Buzzelli D (2008). Calculating friction in one step. Mach Des, 80,
    54–55.
    """
    A = (0.744*log(Re)-1.41)/(1+1.32*eD**0.5)
    B = eD/3.7*Re+2.51*A
    return 1/(A-((A+2*log10(B/Re))/(1+2.18/B)))**2


def f_Vatankhah(Re, eD):
    """
    Calculates friction factor `f` with Vatankhah-Kouchakzadeh corr (2008)

    Parameters
    ------------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]

    References
    ----------
    [1] .. Vatankhah, A. R., Kouchakzadeh, S; Full-range pipe-flow equations";
    Journal of Hydraulic Research Vol. 46, Iss. 4, 2008
    """
    S = 0.124*Re*eD+log(0.4587*Re)
    return 1/(0.8686*log(0.4587*Re/(S-0.31)**(S/(S+0.9633))))**2


def f_avci(Re, eD):
    """
    Calculates friction factor `f` with Avci-Karagoz correlation (2009)

    .. math::
        f = \frac{6.4} {\left\{\ln(Re) - \ln\left[
        1 + 0.01Re\frac{\epsilon}{D}\left(1 + 10(\frac{\epsilon}{D})^{0.5}
        \right)\right]\right\}^{2.4}}

    Parameters
    ------------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]

    References
    ----------
    [1] .. Avci, A.; Karagoz, I.; ""A Novel Explicit Equation for Friction
    Factor in Smooth and Rough Pipes; J. Fluids Eng 131(6), 061203 (2009)
    """
    return 6.4/(log(Re)-log(1+0.01*Re*eD*(1+10*eD**0.5)))**2.4


def f_papaevangelou(Re, eD):
    """
    Calculates friction factor `f` with Papaevangelou correlation (2009)

    .. math::
        f = \frac{0.2479 - 0.0000947(7-\log Re)^4}{\left[\log\left
        (\frac{\epsilon}{3.615D} + \frac{7.366}{Re^{0.9142}}\right)\right]^2}

    Parameters
    ------------
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
    Range of validity
        1e4 <= Re <= 1e7
        1e-5 <= eD <= 1e-3

    References
    ----------
    [1] .. Papaevangelou G, Evangelides C, Tzimopoulos C (2010). A new explicit
    relation for friction coefficient in the Darcy-Weisbach equation.
    Proceedings of the Tenth Conference on Protection and Restoration of the
    Environment 166,1-7pp, PRE10 July 6-09 2010 Corfu, Greece.
    """
    return (0.2479-9.47e-5*(7-log10(Re))**4)/log10(eD/3.615+7.366/Re**0.9142)**2


def f_brkic(Re, eD):
    """
    Calculates friction factor `f` with Brkić correlation (2010a)

    .. math::
        f = [-2\log(10^{-0.4343\beta} + \frac{\epsilon}{3.71D})]^{-2}

        \beta = \ln \frac{Re}{1.816\ln\left(\frac{1.1Re}{\ln(1+1.1Re)}\right)}

    Parameters
    ------------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]

    References
    ----------
    [1] .. Brkić, Dejan (2011). "An Explicit Approximation of Colebrook’s
    equation for fluid flow friction factor". Petroleum Science and Technology
    29 (15): 1596–1602.
    """
    S = log(Re/1.816*log(1.1*Re/log(1+1.1*Re)))
    return 1/(-2*log10(10**(-0.4343*S)+eD/3.71))**2


def f_brkic2(Re, eD):
    """
    Calculates friction factor `f` with Brkić correlation (2010b)
    Second correlation

    .. math::
        f = [-2\log(\frac{2.18\beta}{Re}+ \frac{\epsilon}{3.71D})]^{-2}

        \beta = \ln \frac{Re}{1.816\ln\left(\frac{1.1Re}{\ln(1+1.1Re)}\right)}

    Parameters
    ------------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]

    References
    ----------
    [1] .. Brkić, Dejan (2011). "An Explicit Approximation of Colebrook’s
    equation for fluid flow friction factor". Petroleum Science and Technology
    29 (15): 1596–1602.
    """
    S = log(Re/1.816*log(1.1*Re/log(1+1.1*Re)))
    return 1/(-2*log10(2.18*S/Re+eD/3.71))**2


def f_fang(Re, eD):
    """
    Calculates friction factor `f` with Fang-Xua-Zhou correlation (2011)

    .. math::
        f = 1.613\left\{\ln\left[0.234\frac{\epsilon}{D}^{1.1007}
        - \frac{60.525}{Re^{1.1105}}
        + \frac{56.291}{Re^{1.0712}}\right]\right\}^{-2}

    Parameters
    ------------
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
    Range of validity
        3e3 <= Re <= 1e8
        eD <= 0.05

    References
    ----------
    [1] .. Fang X, Xua Y, Zhou Z (2011). New correlations of single-phase
    friction factor for turbulent pipe flow and evaluation of existing
    single-phase friction factor correlations. Nucl Eng Des 241, 897-902.
    """
    return 1.613/(log(0.234*eD**1.1007-60.525/Re**1.1105+56.291/Re**1.0712))**2


def f_ghanbari(Re, eD):
    """
    Calculates friction factor `f` with Ghanbari correlation (2011)
    Second correlation

    Parameters
    ------------
    Re : float
        Reynolds number, [-]
    eD : float
        Relative roughness of a pipe, [-]

    Returns
    -------
    f : float
        Friction factor, [-]

    References
    ----------
    [1] .. Ghanbari A, Farshad F, Rieke HH (2011). Newly developed friction
    factor correlation for pipe flow and flow assurance. J Chem Eng Mat Sci 2,
    83-86.
    """
    return (-1.52*log10((eD/7.21)**1.042+(2.731/Re)**0.9152))**-2.169


f_list = (f_colebrook, f_chen, f_chen1979, f_moody, f_wood, f_eck, f_altshul,
          f_churchill, f_haaland, f_serghides, f_round, f_swamee, f_jain,
          f_barr, f_zigrang, f_shacham, f_tsal, f_manadilli, f_romeo,
          f_goudar, f_goudar2006, f_buzzelli, f_Vatankhah, f_avci,
          f_papaevangelou, f_brkic, f_brkic2, f_fang, f_ghanbari)


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


def f_friccion(Re, eD=0, metodo=0, geometria=0, adicional=0):
    """
    Generalized method for calculate friction factor for laminar or turbulent
    flux in several geometries

    Input parameters:
    Re: Reynolds number
    eD: Ratio between porosity and internal diameter of pipe, e/D
    metodo
        0   -   Colebrook (default)
        1   -   Chen (1979)
        2   -   Romeo (2002)
        3   -   Goudar-Sonnad
        4   -   Manadilli (1997)
        5   -   Serghides
        6   -   Churchill (1977)
        7   -   Zigrang-Sylvester (1982)
        8   -   Swamee-Jain (1976)
    geometria: Index for duct geometry Ref, Darby pag. 215
        0   -   Circular section (default)
        1   -   Square
        2   -   Isosceles triangle
        3   -   rectangular
        4   -   ellipse
        5   -   scalen triangle
        6   -   Anulli
    adicional: other parameter necessary for geometries not circular
        Triangulo isosceles: angulo del vértice superior
        Rectangulo: realación de longitudes de los lados
        Elipse: Array con ambos diametros
        Triángulo rectangular: ángulo del vértice inferior
    """
    if Re < 2100:
        if geometria == 0:
            f_friccion = 16./Re
        elif geometria == 1:
            f_friccion = 14.2/Re
        elif geometria == 2:
            pass
        elif geometria == 3:
            pass
        elif geometria == 4:
            D, d = adicional[1], adicional[0]
            c = (D-d)/(D+d)
            Dh = 4*d*D*(64-16*c**2)/((d+D)*(64-3*c**4))
            f_friccion = 2*Dh**2*(D**2+d**2)/D**2/d**2/Re
        elif geometria == 5:
            pass
        elif geometria == 6:
            f_friccion = 64./Re

    else:
        if geometria == 6:
            f_friccion = f_Gnielinsky(Re)
        else:
            f_friccion = f_list[metodo](Re, eD)

    return Dimensionless(f_friccion)


def eD(Re, f):
    """Calculates relative roughness

    Parameters
    ------------
    Re : float
        Reynolds number, [-]
    f : float
        Friction factor, [-]

    Returns
    -------
    eD : float
        Relative roughness of a pipe, [-]

    References
    ----------
    [1] .. Colebrook, C. F. and White, C. M. (1937). "Experiments with Fluid
    Friction in Roughened Pipes". Proceedings of the Royal Society of London.
    Series A, Mathematical and Physical Sciences 161 (906): 367–381.
    """
    eD = (10**(-0.5/f**0.5)-2.51/Re/f**0.5)*3.7
    return eD


# Fitting K
# Crane, Flow-of-Fluids-Through-Valve Pag 107
def K_contraction(tita, beta):
    """Tita: angulo de la contracción en grados
    beta: razón entre los diametros antes y despues de la contracción"""
    if tita < 45.:
        K = 0.8*sin(tita*pi/360)*(1-beta**2)/beta**4
    else:
        K = 0.5*sqrt(sin(tita*pi/360))*(1.-beta**2)/beta**4
    return K


def K_enlargement(tita, beta):
    """Tita: angulo del ensanchamiento en grados
    beta: razón entre los diametros antes y despues del ensanchamiento"""
    if tita < 45.:
        K = 2.6*sin(tita*pi/360)*(1.-beta**2)**2/beta**4
    else:
        K = (1.-beta**2)**2/beta**4
    return K


def K_flush(rd):
    """Usando el ajuste exponencial de la tabla disponible"""
    if rd <= 1e-4:
        K = 0.5
    elif rd >= 0.15:
        K = 0.04
    else:
        K = 0.038756579558111+0.45581466480399*exp(-rd/0.041195038092995)
    return K


def K_MitreBend(tita):
    """Usando el ajuste gausiano de la tabla disponible"""
    return -0.31591884532927+29830.477527796*sqrt(2/pi)/122.94894071438*exp(-2*((tita-183.88854928482)/122.94894071438)**2)


def K_longBend(rD):
    if rD <= 1.2:
        K = 20
    elif rD <= 1.7:
        K = 14
    elif rD <= 3.5:
        K = 12
    elif rD <= 5:
        K = 14
    elif rD <= 7:
        K = 17
    elif rD <= 9:
        K = 24
    elif rD <= 11:
        K = 30
    elif rD <= 13:
        K = 34
    elif rD <= 15:
        K = 38
    elif rD <= 18:
        K = 42
    else:
        K = 50
    return K


def Ft(D):
    """Función dependiente del diametro que sirve para definir las K de
    valvulas y accesorios de tuberías D debe ser introducido en mm"""
    if D <= 15:
        ft = 0.027
    elif D <= 20:
        ft = 0.025
    elif D <= 25:
        ft = 0.023
    elif D <= 32:
        ft = 0.022
    elif D <= 40:
        ft = 0.021
    elif D <= 50:
        ft = 0.019
    elif D <= 80:
        ft = 0.018
    elif D <= 100:
        ft = 0.017
    elif D <= 125:
        ft = 0.016
    elif D <= 150:
        ft = 0.015
    elif D <= 250:
        ft = 0.014
    elif D <= 400:
        ft = 0.013
    else:
        ft = 0.012
    return ft


if __name__ == "__main__":
    for f in f_list:
        line = f.__doc__.split("\n")[1]
        year = line.split(" ")[-1]
        name = line.split(" ")[-3]
        doc = name + " " + year
        print(f(1e7, 0.0002), doc)

#    print K_MitreBend(60)
