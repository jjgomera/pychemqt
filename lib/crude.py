#!/usr/bin/python3
# -*- coding: utf-8 -*-


'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2023, Juan José Gómez Romera <jjgomera@gmail.com>

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
# Pseudocomponent definition from crude database
###############################################################################


from math import pi, exp, log, sin
import warnings

from numpy.lib.scimath import log10
from scipy.optimize import fsolve

from lib import unidades
from lib.compuestos import Componente
from lib.petro import Petroleo
from lib.physics import R_atml
from lib.sql import databank
from lib.utilities import refDoc
from tools.qt import tr


__doi__ = {
    1:
        {"autor": "",
         "title": "",
         "ref": "",
         "doi": ""},

    9:
        {"autor": "Ahmed, T.",
         "title": "Equations of State and PVT Analysis: Applications for"
                  "Improved Reservoir Modeling, 2nd Edition",
         "ref": "Gulf Professional Publishing, 2016, ISBN 9780128015704,",
         "doi": "10.1016/B978-0-12-801570-4.00002-7"},
    30:
        {"autor": "Papay, J.A.,",
         "title": "Termelestechnologiai Parameterek Valtozasa a Gazlelepk"
                  "Muvelese",
         "ref": "Soran. OGIL MUSZ, Tud, Kuzl. [Budapest], 1985. pp. 267–273.",
         "doi": ""},
    31:
        {"autor": "Hall, K.R., Yarborough, L.",
         "title": "A New Equation of State for Z-factor Calculations",
         "ref": "Oil and Gas Journal (June 18, 1973): 82–92.",
         "doi": ""},
    32:
        {"autor": "Dranchuk, P.M., Abu-Kassem, J.H.",
         "title": "Calculate of Z factors for Natural Gases Using Equations "
                  "of State",
         "ref": "Journal of Canadian Petroleum Technology (July–September "
                "1975): 34-36",
         "doi": "10.2118/75-03-03"},
    33:
        {"autor": "Dranchuk, P.M., Purvis, R.A., Robinson, D.B.",
         "title": "Computer Calculations of Natural Gas Compressibility "
                  "Factors Using the Standing and Katz Correlation",
         "ref": "Technical Series, no. IP 74–008. Institute of Petroleum, "
                "Alberta, Canada, 1974.",
         "doi": "10.2118/73-112"},
    34:
        {"autor": "Brill, J.P., Beggs, H.D.",
         "title": "Two-Phase Flow in Pipes",
         "ref": "University of Tulsa, INTERCOMP Course, The Hague, 1974",
         "doi": ""},
    35:
        {"autor": "Kumar, N.",
         "title": "Compressibility factors for natural and sour reservoir "
                  "gases by correlations and cubic equations of state",
         "ref": "Thesis of master of science in Petroleum Engineering, 2004, "
                "Texas Tech University.",
         "doi": ""},
    36:
        {"autor": "Gopal, V.N.",
         "title": "Gas Z-Factor Equations Developed for Computer",
         "ref": "Oil and Gas J. (Aug. 8, 1977) 58-60",
         "doi": ""},
    37:
        {"autor": "Sarem, A.M.",
         "title": "Z-Factor Equation Developed for Use in Digital Computers.",
         "ref": "Oil and Gas J. (Sept. 18, 1961) 118",
         "doi": ""},
    38:
        {"autor": "Dranchuk, P.M., Quon, D.",
         "title": "A General Solution of the Equations Describing Steady State"
                  "Turbulent Compressible Flow in Circular Conduits",
         "ref": "Journal of Canadian Petroleum Technology 3(2):60-65, 1964",
         "doi": "10.2118/64-02-04"},
    39:
        {"autor": "Burnett, R.R.",
         "title": "Calculator gives compressibility factors",
         "ref": "Oil & Gas Journal, June 11, 1979, pp. 70-74.",
         "doi": ""},
    40:
        {"autor": "Takacs., G.",
         "title": "Comparing Methods for Calculating Z-factor",
         "ref": "Oil & Gas Journal, May 15, 1989, pp. 43-46.",
         "doi": ""},
    41:
        {"autor": "Sanjari E, Lay E.N.",
         "title": "An accurate empirical correlation for predicting natural "
                  "gas compressibility factors.",
         "ref": "Journal of Natural Gas Chemistry 21(2012):184-188.",
         "doi": "10.1016/s1003-9953(11)60352-6"},
    42:
        {"autor": "Heidaryan E, Moghadasi J, Rahimi M.",
         "title": "New correlations to predict natural gas viscosity and "
                  "compressibility factor.",
         "ref": "Journal of Petroleum Science and Engineering 73 (2010):67-72",
         "doi": "10.1016/j.petrol.2010.05.008"},
    43:
        {"autor": "Heidaryan, E., Salarabadi, A., Moghadasi, J.",
         "title": "A novel correlation approach for prediction of natural gas "
                  "compressibility factor.",
         "ref": "J. Nat. Gas Chem. 19 (2) 2010, 189–192.",
         "doi": "10.1016/s1003-9953(09)60050-5"},
    44:
        {"autor": "Azizi N, Behbahani R, Isazadeh M A.",
         "title": "An efficient correlation for calculating compressibility "
                  "factor of natural gases",
         "ref": "Journal of Natural Gas Chemistry 19 (2010) 642-645",
         "doi": "10.1016/s1003-9953(09)60081-5"},
    45:
        {"autor": "Hall, K.R., Iglesias-Silva, G.A.",
         "title": "Improved equations for the StandingeKatz tables",
         "ref": "Hydrocarb. Process 86 (4), 2007. 107-110",
         "doi": ""},
    46:
        {"autor": "Shokir, Eissa M.El-M., El-Awad, Musaed N., Al-Quraishi, "
                  "Adulhrahman A., Al-Mahdy, Osama A.",
         "title": "Compressibility factor model of sweet, sour, and condensate"
                  " gases using genetic programming",
         "ref": "Chem. Eng. Res. Des. 90 (2012), 785-792.",
         "doi": "10.1016/j.cherd.2011.10.006"},
    47:
        {"autor": "Bahadori, A., Mokhatab, S., Towler, B.F.",
         "title": "Rapidly estimating natural gas compressibility factor",
         "ref": "J. Nat. Gas Chem. 16 (4) 2007, 349-353.",
         "doi": "10.1016/s1003-9953(08)60003-1"},
    48:
        {"autor": "Londono, F.E., Archer, R.A., Blasingame, T.A.",
         "title": "Correlations for hydrocarbon-gas viscosity and gas "
                  "density-validation and correlation of behavior using a "
                  "large-scale database",
         "ref": "SPE Reserv. Evalu. Eng. 8 (6) 2005, 561–572.",
         "doi": "10.2118/75721-PA"},
    49:
        {"autor": "Chankalani, A., Mae'soumi, A., Sameni, A.",
         "title": "An Intelligent Approach for Optimal Prediction of Gas "
                  "Deviation Factor Using Particle Swarm Optimization and "
                  "Genetic Algorithm",
         "ref": "Journal of Natural Gas Science and Engineering 14(2013) "
                "132-143",
         "doi": "10.1016/j.jngse.2013.06.002"},
    50:
        {"autor": "Elsharkawy, A.M.",
         "title": "Efficient methods for calculations of compressibility, "
                  "density and viscosity of natural gases",
         "ref": "Fluid Phase Equilibria 218:1 (2004) 1-13",
         "doi": "10.1016/j.fluid.2003.02.003"}}


# Gas compresibility factor
@refDoc(__doi__, [30, 9])
def Z_Papay(Tr, Pr):
    """Calculate gas compressibility factor using the correlation of Papay
    (1985)

    Parameters
    ----------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]

    Examples
    --------
    >>> "%0.4f" % Z_Papay(2, 3)
    '0.9422'
    """
    Z = 1 - 3.53*Pr/10**(0.9813*Tr) + 0.274*Pr**2/10**(0.8157*Tr)
    return unidades.Dimensionless(Z)


@refDoc(__doi__, [31, 9])
def Z_Hall_Yarborough(Tr, Pr):
    """Calculate gas compressibility factor using the correlation of Hall &
    Yarborough (1973)

    Parameters
    ----------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]
    """
    # Check input in range of validity
    if Tr <= 1:
        warnings.warn("Using extrapolated values")

    X1 = -0.06125*Pr/Tr*exp(-1.2*(1-1/Tr)**2)
    X2 = 14.76/Tr - 9.76/Tr**2 + 4.58/Tr**3
    X3 = 90.7/Tr - 242.2/Tr**2 + 42.4/Tr**3
    X4 = 2.18+2.82/Tr

    def f(Y):
        return X1+(Y+Y**2+Y**3-Y**4)/(1-Y)**3-X2*Y**2+X3*Y**X4

    Yo = 0.0125*Pr/Tr*exp(-1.2*(1-1/Tr)**2)
    Y = fsolve(f, Yo, full_output=True)
    if Y[2] == 1 and abs(Y[1]["fvec"][0]) < 1e-5:
        Z = 0.06125*Pr/Tr/Y[0][0]*exp(-1.2*(1-1/Tr)**2)
    else:
        raise ValueError("Iteration not converge")

    return unidades.Dimensionless(Z)


@refDoc(__doi__, [32])
def Z_Dranchuk_Abu_Kassem(Tr, Pr):
    """Calculate gas compressibility factor using the correlation of Dranchuk-
    Abu-Kassem (1975)

    Parameters
    ----------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]

    Notes
    -----
    Raise :class:`NotImplementedError` if input pair isn't in limit:

        * 1 ≤ Tr ≤ 3
        * 0.2 ≤ Pr ≤ 30
    """
    # Check input in range of validity
    if Tr < 1 or Tr > 3 or Pr < 0.2 or Pr > 30:
        raise NotImplementedError("Incoming out of bound")

    C1 = 0.3265 - 1.07/Tr - 0.5339/Tr**3 + 0.01569/Tr**4 - 0.05165/Tr**5
    C2 = 0.5475 - 0.7361/Tr + 0.1844/Tr**2
    C3 = 0.1056*(-0.7361/Tr + 0.1844/Tr**2)

    # Eq 2
    def f(rho):
        C4 = 0.6134*(1+0.721*rho**2)*rho**2/Tr**3
        Z = 0.27*Pr/Tr/rho
        return 1 + C1*rho + C2*rho**2 - C3*rho**5 + C4*exp(-0.721*rho**2) - Z

    rho0 = 0.27*Pr/Tr
    rho = fsolve(f, rho0, full_output=True)
    if rho[2] == 1:
        Z = 0.27*Pr/Tr/rho[0]
    else:
        raise ValueError("Iteration not converge")

    return unidades.Dimensionless(Z)


@refDoc(__doi__, [33])
def Z_Dranchuk_Purvis_Robinson(Tr, Pr):
    """Calculate gas compressibility factor using the correlation of Dranchuk-
    Purvis-Robinson (1974)

    Parameters
    ----------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]

    Notes
    -----
    Raise :class:`NotImplementedError` if input pair isn't in limit:

        * 1.05 ≤ Tr ≤ 3
        * 0.2 ≤ Pr ≤ 30
    """
    # Check input in range of validity
    if Tr < 1.05 or Tr > 3 or Pr < 0.2 or Pr > 30:
        raise NotImplementedError("Incoming out of bound")

    C1 = 0.31506237 - 1.0467099/Tr - 0.5783272/Tr**3
    C2 = 0.53530771 - 0.61232032/Tr
    C3 = 0.61232032*0.10488813/Tr
    A8 = 0.68446549

    # Eq 3
    def f(rho):
        C4 = 0.68157001*rho**2/Tr**3*(1+A8*rho**2)
        Z = 0.27*Pr/Tr/rho
        return 1 + C1*rho + C2*rho**2 + C3*rho**5 + C4*exp(-A8*rho**2) - Z

    rho0 = 0.27*Pr/Tr
    rho = fsolve(f, rho0, full_output=True)
    if rho[2] == 1:
        Z = 0.27*Pr/Tr/rho[0]
    else:
        raise ValueError("Iteration not converge")

    return unidades.Dimensionless(Z)


@refDoc(__doi__, [34, 35])
def Z_Brill_Beggs(Tr, Pr):
    """Calculate gas compressibility factor using the correlation of Brill-
    Beggs (1974)

    Parameters
    ----------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]

    Notes
    -----
    Raise :class:`NotImplementedError` if input pair isn't in limit:

        * 1.15 ≤ Tr ≤ 2.4
        * 0.2 ≤ Pr ≤ 15
    """
    # Check input in range of validity
    if Tr < 1.15 or Tr > 2.4 or Pr < 0.2 or Pr > 15:
        raise NotImplementedError("Incoming out of bound")

    A = 1.39*(Tr-0.92)**0.5 - 0.36*Tr - 0.101
    B = (0.62-0.23*Tr)*Pr + (0.066/(Tr-0.86)-0.037)*Pr**2 + \
        0.32/10**(9*(Tr-1))*Pr**6
    C = 0.132 - 0.32*log10(Tr)
    D = 10**(0.3016-0.49*Tr+0.1824*Tr**2)
    Z = A + (1-A)/exp(B) + C*Pr**D
    return unidades.Dimensionless(Z)


@refDoc(__doi__, [36, 35])
def Z_Gopal(Tr, Pr):
    """Calculate gas compressibility factor using the correlation of Gopal
    (1974)

    Parameters
    ----------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]

    Notes
    -----
    Raise :class:`NotImplementedError` if input pair isn't in limit:

        * 1.05 ≤ Tr ≤ 3
        * 0.2 ≤ Pr ≤ 15
    """
    if Pr <= 1.2 and Tr <= 3:
        if Tr <= 1.2:
            A, B, C, D = 1.6643, -2.2114, -0.3647, 1.4385
        elif 1.2 <= Tr < 1.4:
            A, B, C, D = 0.0522, -0.8511, -0.0364, 1.0490
        elif 1.4 <= Tr < 2.0:
            A, B, C, D = 0.1391, -0.2988, 0.0007, 0.9969
        elif 2.0 <= Tr <= 3.0:
            A, B, C, D = 0.0295, -0.0825, 0.0009, 0.9967
    elif 1.2 <= Pr < 2.8 and Tr <= 3:
        if Tr <= 1.2:
            A, B, C, D = -1.3570, 1.4942, 4.6315, -4.7009
        elif 1.2 <= Tr < 1.4:
            A, B, C, D = 0.1717, -0.3232, 0.5869, 0.1229
        elif 1.4 <= Tr < 2.0:
            A, B, C, D = 0.0984, -0.2053, 0.0621, 0.8580
        elif 2.0 <= Tr <= 3.0:
            A, B, C, D = 0.0211, -0.0527, 0.0127, 0.9549
    elif 2.8 <= Pr <= 5.4 and Tr <= 3:
        if Tr <= 1.2:
            A, B, C, D = -0.3278, 0.4752, 1.8223, -1.9036
        elif 1.2 <= Tr < 1.4:
            A, B, C, D = -0.2521, 0.3871, 1.6087, -1.6635
        elif 1.4 <= Tr < 2.0:
            A, B, C, D = -0.0284, 0.0625, 0.4714, -0.0011
        elif 2.0 <= Tr <= 3.0:
            A, B, C, D = 0.0041, 0.0039, 0.0607, 0.7927
    elif 5.4 <= Pr < 15:
        Z = Pr*(0.711 + 3.66*Tr)**-1.4667 - 1.637/(0.319*Tr+0.522) + 2.071
        return unidades.Dimensionless(Z)
    else:
        # Input not in range of validity
        raise NotImplementedError("Incoming out of bound")

    Z = Pr*(A*Tr+B) + C*Tr + D
    return unidades.Dimensionless(Z)


@refDoc(__doi__, [35])
def Z_ShellOil(Tr, Pr):
    """Calculate gas compressibility factor using the Shell Oil Company
    correlation (2004)

    Parameters
    ----------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]
    """
    A = -0.101 - 0.36*Tr + 1.3868*(Tr-0.919)**0.5
    B = 0.021 + 0.04275/(Tr-0.65)
    C = 0.6222 - 0.224*Tr
    D = 0.0657/(Tr-0.86) - 0.037
    E = 0.32*exp(-19.53*(Tr-1))
    F = 0.122*exp(-11.3*(Tr-1))
    G = Pr*(C + D*Pr + E*Pr**4)
    Z = A + B*Pr + (1-A)*exp(-G) - F*(Pr/10)**4
    return unidades.Dimensionless(Z)


@refDoc(__doi__, [37])
def Z_Sarem(Tr, Pr):
    """Calculate gas compressibility factor using the Sarem (1969) correlation

    Parameters
    ----------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]

    Notes
    -----
    Raise :class:`NotImplementedError` if input pair isn't in limit:

        * 1.05 ≤ Tr ≤ 2.95
        * 0.1 ≤ Pr ≤ 14.9
    """
    # Check input in range of validity
    if Tr < 1.05 or Tr > 2.95 or Pr < 0.1 or Pr > 14.9:
        raise NotImplementedError("Incoming out of bound")

    x = (2.*Pr-15)/14.8
    y = (2.*Tr-4)/1.9
    Aij = [
        [2.1433504, .0831762, -.0214670, -.0008714, .0042846, -.0016595],
        [.3312352, -.1340361, .0668810, -.0271743, .0088512, -.002152],
        [.1057287, -.0503937, .0050925, .0105513, -.0073182, .0026960],
        [.0521840, .0443121, -.0193294, .0058973, .0015367, -.0028327],
        [.0197040, -.0263834, .019262, -.0115354, .0042910, -.0081303],
        [.0053096, .0089178, -.0108948, .0095594, -.0060114, .0031175]]

    P = [lambda a: 0.7071068,
         lambda a: 1.224745*a,
         lambda a: 0.7905695*(3*a**2-1),
         lambda a: 0.9354145*(5*a**3-3*a),
         lambda a: 0.265165*(35*a**4-30*a**2+3),
         lambda a: 0.293151*(63*a**5-70*a**3+15*a)]

    Z = 0
    for i in range(6):
        for j in range(6):
            Z += Aij[i][j]*P[i](x)*P[j](y)
    return unidades.Dimensionless(Z)


@refDoc(__doi__, [38])
def Z_Leung(Tr, Pr):
    """Calculate gas compressibility factor using the Leung (1964) correlation

    Parameters
    ----------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]

    Notes
    -----
    The correlation is in cited reference, the parameters are least square
    fitting by Leung.

    Raise :class:`NotImplementedError` if input pair isn't in limit:

        * 1.1 ≤ Tr ≤ 2.6
        * 0.5 ≤ Pr ≤ 11
    """
    # Check input in range of validity
    if Tr < 1.1 or Tr > 2.6 or Pr < 0.5 or Pr > 11:
        raise NotImplementedError("Incoming out of bound")

    Bij = [
        [1.877, -4.936, 8.987, -5.215],
        [0.6562, 3.692, -6.477, 3.077],
        [0.1015, -0.5242, 0.8359, -0.3192],
        [-0.00422, 0.0205, -0.0288, 0.00742]]

    Z = 0
    for i in range(4):
        for j in range(4):
            Z += Bij[i][j] * Pr**(i-1) * Tr**(1-j)
    return unidades.Dimensionless(Z)


@refDoc(__doi__, [39, 40])
def Z_Burnett(Tr, Pr):
    """Calculate gas compressibility factor using the Burnett (1979)
    correlation

    Parameters
    ----------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]

    Notes
    -----
    The correlation is in cited reference, the parameters are least square
    fitting by Leung.

    Raise :class:`NotImplementedError` if input pair isn't in limit:

        * 1.3 ≤ Tr ≤ 3
        * 0.2 ≤ Pr ≤ 4
    """
    # FIXME: Don't work
    # Check input in range of validity
    if Tr < 1.1 or Tr > 2.6 or Pr < 0.5 or Pr > 11:
        raise NotImplementedError("Incoming out of bound")

    Zo = 0.3379*log(log(Tr)) + 1.091
    Po = 21.46*Zo - 11.9*Zo**2 - 5.9
    N = (1.1 + 0.26*Tr + (1.04-1.42*Tr)*Pr/Po)*exp(Pr/Po)/Tr
    Z = 1 + (Zo-1) * sin(pi/2*Pr/Po)**N
    return unidades.Dimensionless(Z)


@refDoc(__doi__, [41])
def Z_Sanjari_Lay(Tr, Pr):
    """Calculate gas compressibility factor using the Sanjari-Lay (2012)
    correlation

    Parameters
    ----------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]

    Notes
    -----
    Raise :class:`NotImplementedError` if input pair isn't in limit:

        * 1.01 ≤ Tr ≤ 3
        * 0.01 ≤ Pr ≤ 15
    """
    # Check input in range of validity
    if Tr < 1.01 or Tr > 3 or Pr < 0.01 or Pr > 15:
        raise NotImplementedError("Incoming out of bound")

    # Table 1
    if Pr < 3:
        A = [0, 0.007698, 0.003839, -0.467212, 1.018801, 3.805723, -0.087361,
             7.138305, 0.083440]
    else:
        A = [0, 0.015642, 0.000701, 2.341511, -0.657903, 8.902112, -1.136000,
             3.543614, 0.134041]

    # Eq 16
    Z = 1 + A[1]*Pr + A[2]*Pr**2 + A[3]*Pr**A[4]/Tr**A[5] + \
        A[6]*Pr**(A[4]+1)/Tr**A[7] + A[8]*Pr**(A[4]+2)/Tr**(A[7]+1)
    return unidades.Dimensionless(Z)


@refDoc(__doi__, [43])
def Z_Heidaryan_Salarabadi(Tr, Pr):
    """Calculate gas compressibility factor using the Heidaryan-Salarabadi-
    Moghadasi (2010) correlation

    Parameters
    ----------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]

    Notes
    -----
    Raise :class:`NotImplementedError` if input pair isn't in limit:

        * 1.2 ≤ Tr ≤ 3
        * 0.2 ≤ Pr ≤ 15
    """
    # Check input in range of validity
    if Tr < 1.2 or Tr > 3 or Pr < 0.2 or Pr > 15:
        raise NotImplementedError("Incoming out of bound")

    # Table 1
    A = [0, 1.11532372699824, -.0790395208876, .01588138045027, .0088613449601,
         -2.16190792611599, 1.1575311867207, -0.05367780720737,
         0.01465569989618, -1.80997374923296, 0.95486038773032]

    # Eq 5
    num = A[1] + A[2]*log(Pr) + A[3]*log(Pr)**2 + A[4]*log(Pr)**3 + \
        A[5]/Tr + A[6]/Tr**2
    dem = 1 + A[7]*log(Pr) + A[8]*log(Pr)**2 + A[9]/Tr + A[10]/Tr**2
    Z = log(num/dem)
    return unidades.Dimensionless(Z)


@refDoc(__doi__, [42])
def Z_Heidaryan_Moghadasi(Tr, Pr):
    """Calculate gas compressibility factor using the Heidaryan-Moghadasi-
    Rahimi (2010) correlation

    Parameters
    ----------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]

    Notes
    -----
    Raise :class:`NotImplementedError` if input pair isn't in limit:

        * 1.2 ≤ Tr ≤ 3
        * 0.2 ≤ Pr ≤ 15
    """
    # Check input in range of validity
    if Tr < 1.2 or Tr > 3 or Pr < 0.2 or Pr > 15:
        raise NotImplementedError("Incoming out of bound")

    # Table 1
    if Pr < 3:
        A = [0, 2.827793, -4.688191e-1, -1.262288, -1.536524, -4.535045,
             6.895104e-2, 1.903869e-1, 6.200089e-1, 1.838479, 4.052367e-1,
             1.073574]
    else:
        A = [0, 3.252838, -1.306424e-1, -6.449194e-1, -1.518028, -5.391019,
             -1.379588e-2, 6.600633e-2, 6.120783e-1, 2.317431, 1.632223e-1,
             5.660595e-1]

    # Eq 8
    num = A[1] + A[3]*log(Pr) + A[5]/Tr + A[7]*log(Pr)**2 + A[9]/Tr**2 + \
        A[11]/Tr*log(Pr)
    dem = 1 + A[2]*log(Pr) + A[4]/Tr + A[6]*log(Pr)**2 + A[8]/Tr**2 + \
        A[10]/Tr*log(Pr)
    Z = log(num/dem)
    return unidades.Dimensionless(Z)


@refDoc(__doi__, [44])
def Z_Azizi(Tr, Pr):
    """Calculate gas compressibility factor using the Azizi-Behbahani-Isazadeh
    (2010) correlation

    Parameters
    ----------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]

    Notes
    -----
    Raise :class:`NotImplementedError` if input pair isn't in limit:

        * 1.1 ≤ Tr ≤ 2
        * 0.2 ≤ Pr ≤ 11
    """
    # Check input in range of validity
    if Tr < 1.1 or Tr > 2 or Pr < 0.2 or Pr > 11:
        raise NotImplementedError("Incoming out of bound")

    # Table 1
    a = 0.0373142485385592
    b = -0.0140807151485369
    c = 0.0163263245387186
    d = -0.0307776478819813
    e = 13843575480.943800
    f = -16799138540.763700
    g = 1624178942.6497600
    h = 13702270281.086900
    i = -41645509.896474600
    j = 237249967625.01300
    k = -24449114791.1531
    l = 19357955749.3274
    m = -126354717916.607
    n = 623705678.385784
    o = 17997651104.3330
    p = 151211393445.064
    q = 139474437997.172
    r = -24233012984.0950
    s = 18938047327.5205
    t = -141401620722.689

    A = a*Tr**2.16 + b*Pr**1.028 + c*Pr**1.58/Tr**2.1 + d*log(Tr)**0.5   # Eq 2
    B = e + f*Tr**2.4 + g*Pr**1.56 + h*Pr**0.124*Tr**3.033               # Eq 3
    C = i/log(Tr)**1.28 + j*log(Tr)**1.37 + k*log(Pr) + l*log(Pr)**2 + \
        m*log(Pr)*log(Tr)                                                # Eq 4
    D = 1 + n*Tr**5.55 + o*Pr**0.68*Tr**0.33                             # Eq 5
    E = p*log(Tr)**1.18 + q*log(Tr)**2.1 + r*log(Pr) + s*log(Pr)**2 + \
        t*log(Pr)*log(Tr)                                                # Eq 6
    Z = A + (B + C) / (D + E)                                            # Eq 1
    return unidades.Dimensionless(Z)


@refDoc(__doi__, [46])
def Z_Shokir(Tr, Pr):
    """Calculate gas compressibility factor using the Shokir-Awad-Quraishi
    (2012) correlation

    Parameters
    ----------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]
    """
    # Eq 8
    A = 2.679562*(2*Tr-Pr-1)/((Pr**2+Tr**3)/Pr)
    B = -7.686825*((Pr*Tr+Pr**2)/(Tr*Pr+2*Tr**2+Tr**3))
    C = -0.000624*(Tr**2*Pr-Tr*Pr**2+Tr*Pr**3+2*Tr*Pr-2*Pr**2+2*Pr**3)
    D = 3.067747*(Tr-Pr)/(Pr**2+Tr+Pr)
    E = 0.068059/Tr/Pr + 0.139489*Tr**2 + 0.081873*Pr**2 - 0.041098*Tr/Pr + \
        8.152325*Pr/Tr - 1.63028*Pr + 0.24287*Tr - 2.64988
    Z = A + B + C + D + E
    return unidades.Dimensionless(Z)


@refDoc(__doi__, [47])
def Z_Bahadori(Tr, Pr):
    """Calculate gas compressibility factor using the Bahadori-Mokhatab-Towler
    (2007) correlation

    Parameters
    ----------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]

    Examples
    --------
    Case study in paper
    >>> "%0.4f" % Z_Bahadori(297/197.98, 13860/4287.73)
    '0.7689'

    Notes
    -----
    Raise :class:`NotImplementedError` if input pair isn't in limit:

        * 1.05 ≤ Tr ≤ 2.4
        * 0.2 ≤ Pr ≤ 16
    """
    # Check input in range of validity
    if Tr < 1.05 or Tr > 2.4 or Pr < 0.2 or Pr > 16:
        raise NotImplementedError("Incoming out of bound")

    a = 0.969469 - 1.349238*Tr + 1.443959*Tr**2 - 0.36860*Tr**3         # Eq 8
    b = -0.107783 - 0.127013*Tr + 0.100828*Tr**2 - 0.012319*Tr**3       # Eq 9
    c = 0.0184810 + 0.0523405*Tr - 0.050688*Tr**2 + 0.010870*Tr**3      # Eq 10
    d = -0.000584 - 0.002146*Tr + 0.0020961*Tr**2 - 0.000459*Tr**3      # Eq 11
    Z = a + b*Pr + c*Pr**2 + d*Pr**3                                    # Eq 7
    return unidades.Dimensionless(Z)


@refDoc(__doi__, [48])
def Z_Londono_DAK(Tr, Pr, pure=False):
    """Calculate gas compressibility factor using the Londono-Archer-Blasingame
    (2007) correlation
    This method implement the Dranckuk-Abu-Kassem optimized version from paper

    Parameters
    ----------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]
    pure : boolean
        Aditional parameter to use the parameters of combined dabase

    Returns
    -------
    Z : float
        Gas compressibility factor [-]
    """
    if pure:
        # Eq 34b
        C1 = 0.2965749 - 1.032952/Tr - 0.05394955/Tr**3 - 0.7694/Tr**4 + \
            0.2183666/Tr**5
        C2 = 0.6226256 - 1.006653/Tr + 0.3116857/Tr**2
        C3 = 0.09506539*(-1.006653/Tr + 0.3116857/Tr**2)
        A10 = 0.7544825
        A11 = 0.788

    else:
        # Eq 34a
        C1 = 0.3024696 - 1.046964/Tr - 0.1078916/Tr**3 - 0.7694186/Tr**4 + \
            0.1965439/Tr**5
        C2 = 0.6527819 - 1.118884/Tr + 0.3951957/Tr**2
        C3 = 0.09313593*(-1.118884/Tr + 0.3951957/Tr**2)
        A10 = 0.8483081
        A11 = 0.7880011

    # Eq 2
    def f(rho):
        C4 = A10*(1+A11*rho**2)*rho**2/Tr**3
        Z = 0.27*Pr/Tr/rho
        return 1 + C1*rho + C2*rho**2 - C3*rho**5 + C4*exp(-A11*rho**2) - Z

    rho0 = 0.27*Pr/Tr
    rho = fsolve(f, rho0, full_output=True)
    if rho[2] == 1:
        Z = 0.27*Pr/Tr/rho[0]
    else:
        raise ValueError("Iteration not converge")

    return unidades.Dimensionless(Z)


@refDoc(__doi__, [48])
def Z_Londono_NS(Tr, Pr, pure=False):
    """Calculate gas compressibility factor using the Londono-Archer-Blasingame
    (2007) correlation
    This method implement the Nishiumi-Saito optimized version from paper

    Parameters
    ----------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]
    pure : boolean
        Aditional parameter to use the parameters of combined dabase

    Returns
    -------
    Z : float
        Gas compressibility factor [-]
    """
    if pure:
        # Eq 37b
        A1 = 4.645095e-1
        A2 = 1.627089
        A3 = -9.830729e-1
        A4 = 5.954591e-1
        A5 = 6.183499e-1
        A6 = 4.109793e-1
        A7 = 8.148481e-2
        A8 = 3.541591e-1
        A9 = -1.941089e-2
        A10 = -4.314707e-3
        A11 = 2.789035e-1
        A12 = 7.277907e-1
        A13 = -3.207280e-1
        A14 = 1.756311e-1
        A15 = 7.905733e-1

    else:
        # Eq 37a
        A1 = 2.669857e-1
        A2 = 1.048341
        A3 = -1.516869
        A4 = 4.435926
        A5 = -2.407212
        A6 = 6.089671e-1
        A7 = 5.174665e-1
        A8 = 1.296739
        A9 = -2.892824e-2
        A10 = -1.684037e-2
        A11 = 2.120655
        A12 = -5.046405e-1
        A13 = 1.802678e-1
        A14 = 8.563869e-2
        A15 = 4.956134e-1

    C1 = A1 - A2/Tr - A3/Tr**3 - A4/Tr**4 - A5/Tr**5
    C2 = A6 - A7/Tr - A8/Tr**2 - A9/Tr**5 - A10/Tr**24
    C3 = A11*(A7/Tr + A8/Tr**2 + A9/Tr**5 + A10/Tr**24)

    def f(rho):
        C4 = (A12/Tr**3+A13/Tr**9+A14/Tr**18) * rho**2 * (1+A15*rho**2)
        Z = 0.27*Pr/Tr/rho
        return 1 + C1*rho + C2*rho**2 - C3*rho**5 + C4*exp(-A15*rho**2) - Z

    rho0 = 0.27*Pr/Tr
    rho = fsolve(f, rho0, full_output=True)
    if rho[2] == 1:
        Z = 0.27*Pr/Tr/rho[0]
    else:
        raise ValueError("Iteration not converge")

    return unidades.Dimensionless(Z)


@refDoc(__doi__, [45, 49])
def Z_Hall_Iglesias(Tr, Pr):
    """Calculate gas compressibility factor using the correlation of Hall-
    Iglesias (2007). This is a extension of Hall-Yarborough correlation to
    reduced values of temperature.

    Parameters
    ----------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]
    """
    X1 = 14.54/Tr-8.23/Tr**2+3.39*Tr**3.5
    X2 = 90.7/Tr-242.2/Tr**2+42.4/Tr**3
    X3 = 1.18+2.82/Tr
    k1 = 1/(1.87+0.001*Tr**63)
    k2 = 171.8*Tr**13
    k3 = -3525*(1-exp(-219/Tr**63))

    def f(Z):
        Y = (0.06125*Pr*exp(-1.2*(1-1/Tr)**2))/Tr/Z
        return Z - (1+Y+Y**2-Y**3)/(1-Y)**3 - X1*Y + X2*Y**X3 + \
            k1*Y*exp(-k2*(Y-0.421)**2) + k3*Y**10*exp(-69279*(Y-0.374)**4)

    Z = fsolve(f, 0.5, full_output=True)
    if Z[2] != 1:
        raise ValueError("Iteration not converge")

    return unidades.Dimensionless(Z[0])


@refDoc(__doi__, [50])
def Z_Elsharkawy(Tr, Pr):
    """Calculate gas compressibility factor using the Elsharkawy (2003)
    correlation

    Parameters
    ----------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]
    """
    C1 = 0.3265 - 1.07/Tr - 0.5339/Tr**3 + 0.01569/Tr**4 - 0.05165/Tr**5
    C2 = 0.5475 - 0.7361/Tr + 0.1844/Tr**2
    C3 = 0.1056*(0.7361/Tr + 0.1844/Tr**2)

    # Eq 24
    def f(rho):
        C4 = 0.6134*(1+0.721*rho**2)*rho**2/Tr**3*3
        Z = 0.27*Pr/Tr/rho
        return 1 + C1*rho + C2*rho**2 - C3*rho**5 + C4*exp(-0.721*rho**2) - Z

    rho0 = 0.27*Pr/Tr
    rho = fsolve(f, rho0, full_output=True)
    if rho[2] == 1:
        Z = 0.27*Pr/Tr/rho[0]
    else:
        raise ValueError("Iteration not converge")

    return unidades.Dimensionless(Z)


Z_list = (Z_Hall_Yarborough, Z_Papay, Z_Dranchuk_Abu_Kassem,
          Z_Dranchuk_Purvis_Robinson, Z_ShellOil, Z_Brill_Beggs, Z_Sarem,
          Z_Gopal, Z_Leung, Z_Burnett, Z_Sanjari_Lay, Z_Heidaryan_Salarabadi,
          Z_Heidaryan_Moghadasi, Z_Azizi, Z_Shokir, Z_Bahadori, Z_Londono_DAK,
          Z_Londono_NS, Z_Elsharkawy, Z_Hall_Iglesias)


class Crudo(Petroleo):
    """Class to model a hypotetical component from a defined crude oil
    Clase que define una fracción de petroleo a partir de la base de datos

    Parameters
    ----------
    index : index
        index of crude in crude databank
    Cplus : index
        Number of fraction to split the crude oil for heavier component than C6

    Notes
    -----
    The resultant instancecan be used as hypotetical component as Petroleo.
    """
    kwargs = Petroleo.kwargs.copy()
    kwarg = {"index": 0,
             "Cplus": 0,

             "Rgo": 0.0,
             "gas": None,
             "water": None}
    kwargs.update(kwarg)

    status = 0
    _bool = False
    msg = ""
    hasCurve = False
    hasSG = True
    hasRefraction = False

    def isCalculable(self):
        if self.kwargs["index"]:
            self.status = 1
            self.msg = ""
            return True
        else:
            self.status = 0
            self.msg = tr("pychemqt", "Undefined petrol")

    def calculo(self):
        id = self.kwargs["index"]
        databank.execute("SELECT * FROM CrudeOil WHERE id=='%i'" % id)
        prop = databank.fetchone()

        API = prop[4]
        SG = 141.5/(API+131.5)
        PP = unidades.Temperature(prop[8], "F")
        v100 = prop[10]
        if v100:
            Tb = unidades.Temperature(
                (-753 - 136*(1-exp(-0.15*v100)) + 572*SG - 0.0512*v100 + PP.R)
                / 0.139, "R")
        else:
            Tb = None

        self.definicion = 1
        self.kwargs["name"] = ", ".join(prop[1:3])
        self.kwargs["SG"] = SG
        self.kwargs["Tb"] = Tb
        self.kwargs["S"] = prop[5]
        self.kwargs["N"] = prop[6]
        self.kwargs["v100"] = prop[10]
        Petroleo.calculo(self)

        self.vanadium = prop[11]
        self.nickel = prop[12]
        self.carbonResid = prop[13]
        self.asphaltene = prop[14]
        self.S2 = prop[6]
        self.H2S = prop[18]
        self.nNeutralization = prop[19]
        self.ash = prop[21]
        self.salt = prop[22]
        self.water = prop[20]
        self.NPentane = prop[15]
        if prop[16]:
            self.reidVP = unidades.Pressure(prop[16], "psi")
        else:
            self.reidVP = None
        if prop[17]:
            self.FlashP = unidades.Temperature(prop[17], "F")
        self.PourP = PP

        self.C1 = prop[23]/100.
        self.C2 = prop[24]/100.
        self.C3 = prop[25]/100.
        self.iC4 = prop[26]/100.
        self.nC4 = prop[27]/100.
        self.iC5 = prop[28]/100.
        self.nC5 = prop[29]/100.

        # SGo = 0.7
        # SG_ = (SG-SGo)/SGo
        # B = 3.
        # A = SG_**3/0.619**3
        # Cplus = int(self.kwargs["Cplus"])
        # Tbi=[unidades.Temperature(1090-exp(6.9955-0.11193*Nc**(2./3))) for Nc in range(6, Cplus)]
        # SGi=[1.07-exp(3.65097-3.8864*Nc**0.1) for Nc in range(6, Cplus)]
        # x=[1-1/exp(A/B*(SG-SGo)**B/SGo**B) for SG in SGi]
        # Mi=[prop_Riazi_Alsahhaf(1, g, reverse=True) for g in SGi]
        # APIi=[141.5/SG-131.5 for SG in SGi]
        # Kwi=[Tbi[i].R**(1./3)/SGi[i] for i in range(len(Mi))]
        # di=[unidades.Density(prop_Riazi_Alsahhaf(2, M), "gcc") for M in Mi]
        # Ii=[prop_Riazi_Alsahhaf(3, M) for M in Mi]
        # Tci=[unidades.Temperature(Tbi[i]/prop_Riazi_Alsahhaf(4, M)) for i, M in enumerate(Mi)]
        # Pci=[unidades.Pressure(prop_Riazi_Alsahhaf(5, M), "bar") for M in Mi]
        # Vci=[unidades.SpecificVolume(1/prop_Riazi_Alsahhaf(6, M), "ccg") for M in Mi]
        # Wi=[prop_Riazi_Alsahhaf(7, M) for M in Mi]
        # Tensioni=[unidades.Tension(prop_Riazi_Alsahhaf(8, M), "dyncm") for M in Mi]
        # ParSoli=[unidades.SolubilityParameter(prop_Riazi_Alsahhaf(9, M), "calcc") for M in Mi]


    def pb_Standing(self, T):
        """Standing, M.B.: Volumetric and Phase Behavior of Oil Field Hydrocarbon Systems, SPE, Dallas (1977)"""
        t=unidades.Temperature(T)
        F=(self.Rgo.ft3bbl/self.gas.SG)**0.83*10**(0.00091*t.F-0.0125*self.API)
        return unidades.Pressure(18.2*(F-1.4), "psi")

    def pb_Lasater(self, T):
        """Lasater, J.A: "Bubble Point Pressure Correlation," Trans., AIME (1958) 213, 379-381"""
        t=unidades.Temperature(T)
        if self.API<=40:
            M=630-10.*self.API
        else:
            M=73110.*self.API**-1.562
        yg=self.Rgo.ft3bbl/379.3/(self.Rgo.ft3bbl/379.3+350*self.SG/M)
        if yg<=0.6:
            pb=0.679*exp(2.786*yg)-0.323
        else:
            pb=8.26*yg**3.56+1.95
        return unidades.Pressure(pb*t.R/self.gas.SG, "psi")

    def pb_Vazquez_Beggs(self, T, ts=350, ps=10):
        """Vázquez, M.E. and Beggs, H.D.: "Correlations for Fluid Physical Property Prediction," J.Pet. Tech. (June 1980), 968-970"""
        t=unidades.Temperature(T)
        ts=unidades.Temperature(ts)
        ps=unidades.Pressure(ps, "atm")
        if self.API<=30:
            C1, C2, C3=0.0362, 1.0937, 25.724
        else:
            C1, C2, C3=0.0178, 1.187, 23.931

        gravity_corr=self.gas.SG*(1.+5.912e-5*self.API*ts.F*log10(ps.psi/114.7))
        return unidades.Pressure((self.Rgo.ft3bbl/C1/gravity_corr/exp(C3*self.API/T.R))**(1./C2), "psi")

    def pb_Glaso(self, T):
        """Glaso, O.: "Generalized Pressure-Volume-Temperature Correlations," J. Pet. Tech. (May 1980), 785-795"""
        t=unidades.Temperature(T)
        F=(self.Rgo.ft3bbl/self.gas.SG)**0.816*t.F**0.172/self.API**0.989
        return unidades.Pressure(10**(1.7669+1.7447*log10(F)-0.30218*log10(F)**2), "psi")

    def pb_Total(self, T):
        """TOTAL Compagnie Francaise Des Petroles: "Proyectos de Inyección de Fluidos - Correlaciones PVT para Crudos del Oriente de Venezuela," S.A. MENEVEN, Sept. 1983"""
        t=unidades.Temperature(T)
        if self.API<=10:
            C1, C2, C3, C4=12.847, 0.9636, 0.000993, 0.03417
        elif self.API<=35:
            C1, C2, C3, C4=25.2755, 0.7617, 0.000835, 0.011292
        else:
            C1, C2, C3, C4=216.4711, 0.6922, -0.000427, 0.02314
        return unidades.Pressure(C1*(self.Rgo.ft3bbl/self.gas.SG)**C2*10**(C3*t.F-C4*self.API), "psi")

    def pb_Al_Marhoun(self, T):
        """Al-Marhoun, M.A.: "PVT Correlation for Middle East Crude Oils," J. Pet. Tech (May 1988), 650-666"""
        t=unidades.Temperature(T)
        return unidades.Pressure(5.38088e-3*self.Rgo.ft3bbl**0.715082*self.gas.SG**-1.87784*self.SG**3.1437*t.R**1.32657, "psi")

    def pb_Dokla_Osman(self, T):
        """Dokla, M.E. and Osman, M.E.: "Correlation of PVT properties for UAE Crudes," Trans., AIME (1992) 293, 41-46"""
        t=unidades.Temperature(T)
        return unidades.Pressure(0.836386e4*self.Rgo.ft3bbl**0.724047*self.gas.SG**-1.01049*self.SG**0.107991*t.R**-0.952584, "psi")

    def pb_Petrosky_Farshad(self, T):
        """Petrosky, G.E., Jr. and Farshad, F.F.: "Pressure-Volume-Temperature Correlations for Gulf of Mexico Crude Oils," paper SPE 26644 presented at the 68th Annual Technical Conference and Exhibition, Houston, Texas, Oct. 3-6,1993."""
        t=unidades.Temperature(T)
        F=self.Rgo.ft3bbl**0.5774/self.gas.SG**0.8439*10**(4.561e-5*t.F**1.3911-7.916e-4*self.API**1.541)
        return unidades.Pressure(112.727*(F-12.34), "psi")

    def pb_Kartoatmodjo_Schmidt(self, T, ts=350, ps=10):
        """Kartoatmodjo, T. and Schmidt, Z.: "Large Data Bank Improve Crude Physical Property Correlations," Oil and Gas J. (July 4, 1994) 51-55"""
        t=unidades.Temperature(T)
        ts=unidades.Temperature(ts)
        ps=unidades.Pressure(ps, "atm")
        if self.API<=30:
            C1, C2, C3, C4=0.05958, 0.7972, 13.1405, 0.9986
        else:
            C1, C2, C3, C4=0.0315, 0.7587, 11.2895, 0.9143

        gravity_corr=self.gas.SG*(1.+0.1595*self.API**0.4078*ts.F**-0.2506*log10(ps.psi/114.7))
        return unidades.Pressure((self.Rgo.ft3bbl/C1/gravity_corr**C2/10**(C3*self.API/t.R))**C4, "psi")

    def pb(self, T):
        """Presión de burbujeo del crudo"""
        t=unidades.Temperature(T)
        metodo_pb=[self.pb_Standing, self.pb_Lasater, self.pb_Vazquez_Beggs, self.pb_Glaso, self.pb_Total, self.pb_Al_Marhoun, self.pb_Dokla_Osman, self.pb_Petrosky_Farshad, self.pb_Kartoatmodjo_Schmidt][Preferences.getint("petro", "pb")]
        pb=metodo_pb(t)
        CN2=1.+((-2.65e-4*self.API+5.5e-3)*t.F+(0.0931*self.API-0.895))*self.gas.N2
        CCO2=1.-693.8*self.gas.CO2*t.F**-1.553
        CH2S=1.-(0.9035+0.0015*self.API)*self.gas.H2S+0.019*(45-self.API)*self.gas.H2S**2
        return unidades.Pressure(CN2*CH2S*CCO2*pb)


    def B_Standing(self, T):
        """Standing, M.B.: Volumetric and Phase Behavior of Oil Field Hydrocarbon Systems, SPE, Dallas (1977)"""
        t=unidades.Temperature(T)
        F=self.Rgo.ft3bbl*(self.gas.SG/self.SG)**0.5+1.25*t.F
        return 0.9759+12e-5*F**1.2

    def B_Vazquez_Beggs(self, T, ts=350, ps=10):
        """Vázquez, M.E. and Beggs, H.D.: "Correlations for Fluid Physical Property Prediction," J.Pet. Tech. (June 1980), 968-970"""
        t=unidades.Temperature(T)
        ts=unidades.Temperature(ts)
        ps=unidades.Pressure(ps, "atm")
        if self.API<=30:
            C1, C2, C3=4.667e-4, 1.751e-5, -1.8106e-6
        else:
            C1, C2, C3=4.67e-4, 1.1e-5, 1.337e-9

        gravity_corr=self.gas.SG*(1.+5.912e-5*self.API*ts.F*log10(ps.psi/114.7))
        return 1.+C1*self.Rgo.ft3bbl+C2*(t.F-60)*self.API/gravity_corr+C3*self.Rgo.ft3bbl*(t.F-60)*self.API/gravity_corr

    def B_Glaso(self, T):
        """Glaso, O.: "Generalized Pressure-Volume-Temperature Correlations," J. Pet. Tech. (May 1980), 785-795"""
        t=unidades.Temperature(T)
        F=self.Rgo.ft3bbl(self.gas.SG/self.SG)**0.526*+0.968*t.F
        return 1.+10**(-6.58511+2.91329*log10(F)-0.27683*log10(F)**2)

    def B_Total(self, T):
        """TOTAL Compagnie Francaise Des Petroles: "Proyectos de Inyección de Fluidos - Correlaciones PVT para Crudos del Oriente de Venezuela," S.A. MENEVEN, Sept. 1983"""
        t=unidades.Temperature(T)
        return 1.022+4.857e-4*self.Rgo.ft3bbl-2.009e-6*(t.F-60)*self.API/self.gas.SG+17.569e-9*self.Rgo.ft3bbl*(t.F-60)*self.API/self.gas.SG

    def B_Al_Marhoun(self, T):
        """Al-Marhoun, M.A.: "PVT Correlation for Middle East Crude Oils," J. Pet. Tech (May 1988), 650-666"""
        t=unidades.Temperature(T)
        F=self.Rgo.ft3bbl**0.74239*self.gas.SG**0.323294*self.SG**-1.20204
        return 0.497069+0.862963e-3*t.R+0.182594e-2*F+0.318099e-5*F**2

    def B_Dokla_Osman(self, T):
        """Dokla, M.E. and Osman, M.E.: "Correlation of PVT properties for UAE Crudes," Trans., AIME (1992) 293, 41-46"""
        t=unidades.Temperature(T)
        F=self.Rgo.ft3bbl**0.773572*self.gas.SG**0.40402*self.SG**-0.882605
        return 0.431935e-1+0.156667e-2*t.R+0.139775e-2*F+0.380525e-5*F**2

    def B_Petrosky_Farshad(self, T):
        """Petrosky, G.E., Jr. and Farshad, F.F.: "Pressure-Volume-Temperature Correlations for Gulf of Mexico Crude Oils," paper SPE 26644 presented at the 68th Annual Technical Conference and Exhibition, Houston, Texas, Oct. 3-6,1993."""
        t=unidades.Temperature(T)
        F=self.Rgo.ft3bbl**0.3738*self.gas.SG**0.2914/self.SG**0.6265+0.24626*t.F**0.5371
        return 1.0113+7.2046e-5*F**3.0936

    def B_Kartoatmodjo_Schmidt(self, T, ts=350, ps=10):
        """Kartoatmodjo, T. and Schmidt, Z.: "Large Data Bank Improve Crude Physical Property Correlations," Oil and Gas J. (July 4, 1994) 51-55"""
        t=unidades.Temperature(T)
        ts=unidades.Temperature(ts)
        ps=unidades.Pressure(ps, "atm")

        gravity_corr=self.gas.SG*(1.+0.1595*self.API**0.4078*ts.F**-0.2506*log10(ps.psi/114.7))
        F=self.Rgo.ft3bbl**0.755*gravity_corr**0.25*self.SG**-1.5+0.45*t.F
        return 0.98496+1e-4*F**1.5

    def B(self, T, P):
        """Factor volumétrico, relación entre el volumen a las condiciones del yacimiento y las condiciones normales a la presión de burbujeo"""
        t=unidades.Temperature(T)
        p=unidades.Pressure(P, "atm")
        metodo_B=[self.pb_Standing, self.pb_Vazquez_Beggs, self.pb_Glaso, self.pb_Total, self.pb_Al_Marhoun, self.pb_Dokla_Osman, self.pb_Petrosky_Farshad, self.pb_Kartoatmodjo_Schmidt][Preferences.getint("petro", "Vol_factor")]
        Bpb=metodo_B(t)
        pb=self.pb(T)
        if p>pb:
            c=self.compressibilidad(T)
            return Bpb*exp(c*(pb-p))
        else:
            return Bpb


    def Mu_Gas(self, T):
        """Método de cáluclo de la viscosidad de vapores, API procedure 11B3.1, pag 1105"""
        #FIXME: no sale
        t=unidades.Temperature(T)
        return unidades.Viscosity(-0.0092696+t.R**0.5*(0.0010310+4.4507e-5*self.M**0.5)+1.1249e-5*self.M, "cP")
#        return unidades.Viscosity(-0.0092696+T**0.5*(0.001383-5.9712e-5*self.M**0.5)+1.1249e-5*self.M, "cP")


    def Mu_Beal(self, T):
        """Beal, C.: "The Viscosity of Air, Water, Natural Gas, crude Oil and its Associated Gases at Oil-Field Temperatures and Pressures," Trans., AIME (1946) 165, 94-115"""
        t=unidades.Temperature(T)
        a=10**(0.43+8.33/self.API)
        mu=(0.32+1.8e7/self.API**4.53)*(360/(t.F+200))**a
        return unidades.Viscosity(mu, "cP")

    def Mu_Beggs_Robinson(self, T):
        """Beggs, H.D. and Robinson, J.R.: "Estimating the Viscosity of Crude Oil Systems," J. Pet. Tech. Forum (Sept. 1975), 1140-1141"""
        t=unidades.Temperature(T)
        z=3.0324-0.02023*self.API
        y=10**z
        x=y*t.F**-1.163
        return unidades.Viscosity(10**x-1, "cP")

    def Mu_Glaso(self, T):
        """Glaso, O.: "Generalized Pressure-Volume-Temperature Correlations," J. Pet. Tech. (May 1980), 785-795"""
        t=unidades.Temperature(T)
        mu=3.141e10*t.F**-3.444*log10(self.API)**(10.313*log10(t.F)-36.447)
        return unidades.Viscosity(mu, "cP")

    def Mu_Egbogah(self, T):
        """Egbogah, E.O.: "An Improved Temperature-Viscosity Correlation for Crude Oil Systems," paper 83-34-32 presented at the 1983 Annual Technical Meeting of the Petroleum Society of CIM, Banff, Alberta, May 10-13, 1983"""
        t=unidades.Temperature(T)
        mu=1.8653-0.025086*self.API-0.5644*log10(t.F)
        return unidades.Viscosity(10**(10**mu)-1, "cP")

    def Mu_Kartoatmodjo_Schmidt(self, T):
        """Kartoatmodjo, T. and Schmidt, Z.: "Large Data Bank Improve Crude Physical Property Correlations," Oil and Gas J. (July 4, 1994) 51-55"""
        t=unidades.Temperature(T)
        mu=16e8*t.F**-2.8177*log10(self.API)**(5.7526*log10(t.F)-26.9718)
        return unidades.Viscosity(mu, "cP")

    def Mu_Muerto(self, T):
        """Viscosidad de petroleos muertos (sin gas disuelto)"""
        metodos=[self.Mu_Beal, self.Mu_Beggs_Robinson, self.Mu_Glaso, self.Mu_Egbogah, self.Mu_Kartoatmodjo_Schmidt][Preferences.getint("petro", "mu_dead")]
        return metodos(T)

    def Mu_Chew_Connally(self, T, R):
        """Chew, J.N. and Connally, C.A. Jr.: "A Viscosity Correlation for Gas-Saturated Crude Oils," Trans, AIME (1959) 216, 23-25"""
        t=unidades.Temperature(T)
        r=unidades.V2V(R)
        muo=self.Mu_Muerto(T)
        A=10**(r.ft3bbl*(2.2e-7*r.ft3bbl-7.4e-4))
        b=0.68/10**(8.62e-5*r.ft3bbl)+0.25/10**(1.1e-3*r.ft3bbl)+0.062/10**(3.74e-3*r.ft3bbl)
        return unidades.Viscosity(A*muo.cP**b, "cP")

    def Mu_Beggs_Robinson_vivo(self, T, R):
        """Beggs, H.D. and Robinson, J.R.: "Estimating the Viscosity of Crude Oil Systems," J. Pet. Tech. Forum (Sept. 1975), 1140-1141"""
        t=unidades.Temperature(T)
        r=unidades.V2V(R)
        muo=self.Mu_Muerto(T)
        A=10.715*(r.ft3bbl+100)**-0.515
        b=5.44*(r.ft3bbl+150)**-0.338
        return unidades.Viscosity(A*muo.cP**b, "cP")

    def Mu_Kartoatmodjo_Schmidt_vivo(self, T, R):
        """Kartoatmodjo, T. and Schmidt, Z.: "Large Data Bank Improve Crude Physical Property Correlations," Oil and Gas J. (July 4, 1994) 51-55"""
        t=unidades.Temperature(T)
        r=unidades.V2V(R)
        muo=self.Mu_Muerto(T)
        b=10**(-0.00081*r.ft3bbl)
        A=(0.2001+0.8428*10**(-0.000845*r.ft3bbl))*muo.cP**(0.43+0.5165*b)
        return unidades.Viscosity(-0.06821+0.9824*A+40.34e-5*A**2, "cP")

    def Mu_Vivo(self, T, R):
        """Viscosidad de petroleos vivos (con gas disuelto)"""
        metodos=[self.Mu_Chew_Connally, self.Mu_Beggs_Robinson_vivo, self.Mu_Kartoatmodjo_Schmidt_vivo][Preferences.getint("petro", "mu_live")]
        return metodos(T, R)


    def Mu_Beal_presion(self, T, P, R):
        """Beal, C.: "The Viscosity of Air, Water, Natural Gas, crude Oil and its Associated Gases at Oil-Field Temperatures and Pressures," Trans., AIME (1946) 165, 94-115"""
        p=unidades.Pressure(P, "atm")
        pb=self.pb(T)
        muo=self.Mu_Vivo(T, R)
        mu=(0.024*muo.cP**1.6+0.038*muo.cP**0.56)*0.001*(p.psi-pb.psi)+muo.cP
        return unidades.Viscosity(mu, "cP")

    def Mu_Vazquez_Beggs(self, T, P, R):
        """Vázquez, M.E. and Beggs, H.D.: "Correlations for Fluid Physical Property Prediction," J.Pet. Tech. (June 1980), 968-970"""
        p=unidades.Pressure(P, "atm")
        pb=self.pb(T)
        muo=self.Mu_Vivo(T, R)
        m=2.6*p.psi**1.187*exp(-11.513-8.98e-5*p.psi)
        return unidades.Viscosity(muo.cP*(p.psi/pb.psi)**m, "cP")

    def Mu_Kartoatmodjo_Schmidt_presion(self, T, P, ):
        """Kartoatmodjo, T. and Schmidt, Z.: "Large Data Bank Improve Crude Physical Property Correlations," Oil and Gas J. (July 4, 1994) 51-55"""
        p=unidades.Pressure(P, "atm")
        pb=self.pb(T)
        muo=self.Mu_Vivo(T, R)
        mu=1.00081*muo.cP+1.127e-3*(p.psi-pb.psi)*(-65.17e-4*muo.cP**1.8148+0.038*muo.cP**1.59)
        return unidades.Viscosity(mu, "cP")

    def Mu_Presion(self, T, P):
        """Viscosidad de petroleos vivos (con gas disuelto)"""
        metodos=[self.Mu_Beal_presion, self.Mu_Vazquez_Beggs, self.Mu_Kartoatmodjo_Schmidt_presion][Preferences.getint("petro", "mu_live")]
        return metodos(T, P)


class Water(Componente):
    """Clase que define el agua específica que acompaña al petroleo, con las propiedades específicas"""
    def __init__(self):
        Componente.__init__(self, 62)

    def Solubilidad_Culberson_McKetta(self, T, P, S=0):
        """Culberson, O.L. and McKetta, J.J., Jr.: "Phase Equilibria in Hydrocarbon-Water Systems III - The solubility of Methane in Water at Pressures to 10,000 psia," Trans., AIME (1951) 192, 223-226
        S: salinidad en % en peso"""
        t=unidades.Temperature(T)
        p=unidades.Pressure(P, "atm")
        A=8.15839-6.12265e-2*t.F+1.91663e-4*t.F**2-2.1654e-7*t.F**3
        B=1.01021e-2-7.44241e-5*t.F+3.05553e-7*t.F**2-2.94883e-10*t.F**3
        C=(-9.02505+0.130237*t.F-8.53425e-4*t.F**2+2.34122e-6*t.F**3-2.37049e-9*t.F**4)*1e-7
        R=A+B*p.psi+C*p.psi**2
        Rs=R*10**(-0.0840655*S*t.F**-0.285854)
        return unidades.V2V(Rs, "ft3bbl")

    def Solubilidad_McCoy(self, T, P, S=0):
        """McCoy, R.L.: Microcomputer Programs for Petroleum Engineers: Vol. 1, Reservoir Engineering and Formation Evaluation, Gulf Publishing Co., Houston (1983).
        S: salinidad en % en peso"""
        t=unidades.Temperature(T)
        p=unidades.Pressure(P, "atm")
        A=2.12+3.45e-3*t.F-3.59e-5*t.F**2
        B=0.0107-5.26e-5*t.F+1.48e-7*t.F**2
        C=-8.75e-7+3.9e-9*t.F-1.02e-11*t.F**2
        R=A+B*p.psi+C*p.psi**2
        Rs=R*(1-(0.0753-1.73e-4*t.F)*S)
        return unidades.V2V(Rs, "ft3bbl")


    def Factor_Volumetrico_McCain(self, T, P, S=0):
        """McCain, W.D., Jr: The Properties of Petroleum Fluids, 2nd ed. Tulsa, OK: PennWell Books, 1990.
        S: salinidad en % en peso"""
        t=unidades.Temperature(T)
        p=unidades.Pressure(P, "atm")
        DVp=-1.0001e-2+1.33391e-4*t.F+5.50654e-7*t.F**2
        DVt=-1.95301e-9*p.psi*t.F-1.72834e-13*p.psi**2*t.F-3.58922e-7*p.psi-2.25341e-10*p.psi**2
        B=(1+DVp)*(1+DVt)
        return B*(1+S*(5.1e-8*p.psi+(5.47e-6-1.95e-10*p.psi)*(t.F-60)-(3.23e-8-8.5e-13*p.psi)*(t.F-60)**2))

    def Factor_Volumetrico_McCoy(self, T, P, S=0, R=1):
        """McCoy, R.L.: Microcomputer Programs for Petroleum Engineers: Vol. 1, Reservoir Engineering and Formation Evaluation, Gulf Publishing Co., Houston (1983).
        R: razon de gas disuelto
        S: salinidad en % en peso"""
        t=unidades.Temperature(T)
        p=unidades.Pressure(P, "atm")
        if R==0:
            A=0.9947+5.8e-6*t.F+1.02e-6*t.F**2
            B=-4.228e-6+1.8276e-8*t.F-6.77e-11*t.F**2
            C=1.3e-10-1.3855e-12*t.F+4.285e-15*t.F**2
        else:
            A=0.9911+6.35e-5*t.F+8.5e-7*t.F**2
            B=-1.093e-6-3.497e-9*t.F+4.57e-12*t.F**2
            C=-5e-11+6.429e-13*t.F-1.43e-15*t.F**2
        B=A+B*p.psi+C*p.psi**2
        return B*(1+S*(5.1e-8*p.psi+(5.47e-6-1.95e-10*p.psi)*(t.F-60)-(3.23e-8-8.5e-13*p.psi)*(t.F-60)**2))

    def Rho(self, T, P, S=0):
        B=self.Factor_Volumetrico_McCain(T, P, S)
        s=S*1e7/58443
        g=1.+0.695e-6*s
        return unidades.Density(62.4*g/B, "lbft3")

    def Rho_McCain(self, T, P, S=0):
        """McCain, W.D., Jr: The Properties of Petroleum Fluids, 2nd ed. Tulsa, OK: PennWell Books, 1990."""
        B=self.Factor_Volumetrico_McCain(T, P, S)
        g=62.368+0.438603*S+1.60074e-3*S**2
        return unidades.Density(g/B, "lbft3")

    def Compresibilidad_Dodson_Standing(self, T, P, S=0, R=0):
        """Dodson, C.R. and Standing, M.B.: "Pressure-Volume-Temperature and Solubility RElations for Natural Gas-Water-Mixtures," Drill. and Prod. Prac., API (1944) 173-179"""
        t=unidades.Temperature(T)
        p=unidades.Pressure(P, "atm")
        R=unidades.V2V(R)
        a=3.8546-1.34e-4*p.psi
        b=-0.01052+4.77e-7*p.psi
        c=3.9267e-5-8.8e-10*p.psi
        cw=(a+b*t.F+c*t.F**2)/1e6
        cr=1.+8.9e-3*R.ft3bbl
        cs=1+S**0.7*(-5.2e-2+2.7e-4*t.F-1.14e-6*t.F**2+1.121e-9*t.F**3)
        return cw*cr*cs

    def Compresibilidad_Osif(self, T, P, S=0):
        """Osif, T.L.: "The Effects of Salt, Gas, Temperature and Pressure on the Compressibility of Water," SPE Res.Eng. (Feb. 1988) 3, No.1. 175-181"""
        t=unidades.Temperature(T)
        p=unidades.Pressure(P, "atm")
        S=S*1e4/58443
        return 1/(7.033*p.psi+541.5*S-537.*t.F+403300.)


    def Mu_Van_Wingen(self, T):
        """Van Wingen, N.: "Viscosity of Air, Water, Natural Gas, and Crude Oil at Varying Pressure and Temperatures," Secondary Recovery of Oil in the United States, API (1950) 127."""
        t=unidades.Temperature(T)
        return unidades.Viscosity(exp(1.003-1.479e-2*t.F+1.982e-5*t.F**2), "cP")

    def Mu_Mattews_Russel(self, T, P, S=0):
        """Mathews, C.S and Russel, D.G.: Pressure Buildup and Flow Text in Wells. Monograph Series. Society of Petroleum Engineers of AIME, Dallas (1967) 1, Appendix G."""
        t=unidades.Temperature(T)
        p=unidades.Pressure(P, "atm")
        A=-0.04518+0.009313*S-0.000383*S**2
        B=70.634+0.09576*S**2
        mu=A+B/t.F
        f=1.+3.5e-12*p.psi**2*(t.F-40.)
        return unidades.Viscosity(mu*f, "cP")

    def Mu_McCain(self, T, P, S=0):
        """McCain, W.D., Jr: The Properties of Petroleum Fluids, 2nd ed. Tulsa, OK: PennWell Books, 1990."""
        t=unidades.Temperature(T)
        p=unidades.Pressure(P, "atm")
        A=109.574-8.40564*S+0.313314*S**2+8.72213e-3*S**3
        B=-1.12166+2.63951e-2*S-6.79461e-4*S**2-5.47119e-5*S**3+1.55586e-6*S**4
        mu=A*t.F**B
        f=0.9994+4.0295e-5*p.psi+3.1062e-9*p.psi**2
        return unidades.Viscosity(mu*f, "cP")

    def Mu_McCoy(self, T, S=0):
        """McCoy, R.L.: Microcomputer Programs for Petroleum Engineers: Vol. 1, Reservoir Engineering and Formation Evaluation, Gulf Publishing Co., Houston (1983)."""
        t=unidades.Temperature(T)
        mu=0.02414*10**(247.8/(t-140))
        f=1.-1.87e-3*S**0.5+2.18e-4*S**2.5+(t.F**0.5-1.35e-2*t.F)*(2.76e-3*S-3.44e-4*S**1.5)
        return unidades.Viscosity(mu*f, "cP")

    def Tension_Jennings_Newman(self, T, P):
        """Jennings, H.Y., Jr. and Newman, G.I.L.:"The Effect of Temperature and Pressure on the Interfacial Tension of Water Against Mechane-Normal Decane Mixtures," Trans., AIME (1971) 251, 171-175."""
        t=unidades.Temperature(T)
        p=unidades.Pressure(P, "atm")
        A=79.1618-0.118978*t.F
        B=-5.28473e-3+9.87913e-6*t.F
        C=(2.33814-4.57194e-4*t.F-7.52678e-6*t.F**2)*1e-7
        return unidades.Tension(A+B*p.psi+C*p.psi**2, "dyncm")


class Natural_Gas(object):
    """Clase que define un gas natural como fracción desconocida de petroleo"""
    def __init__(self, SG=None, composicion=[], wet=False, CO2=0., H2S=0., N2=0.):
        """
        g: gravedad específica
        composicion: en la forma de array [[componentes],[fracciones molares]]
        wet: parámetro opcional que indica si se trata de un gas húmedo (con una pequeña fracción de fase líquida
        CO2: fracción molar de CO2 en el gas
        H2S: fracción molar de H2S en el gas
        """
        self.SG=SG
        self.wet=wet
        self.CO2=CO2
        self.H2S=H2S
        self.N2=N2
        if wet:
            self.tpc=unidades.Temperature(187.+330.*SG-72.5*SG**2, "R")
            self.ppc=unidades.Pressure(706-51.7*SG-11.1*SG**2, "psi")
        else:
            self.tpc=unidades.Temperature(168.+325.*SG-12.5*SG**2, "R")
            self.ppc=unidades.Pressure(677+15.*SG-37.5*SG**2, "psi")
        if N2!=0:
            self.tpc, self.ppc=self.Critical_Carr_Kobayashi_Burrows()
        if CO2+H2S!=0.:
            self.tpc, self.ppc=self.Critical_Wichert_Aziz()
        self.M=28.96*SG

    def Critical_Wichert_Aziz(self):
        """Wichert, E., and K. Aziz. “Calculation of Z’s for Sour Gases.” Hydrocarbon Processing 51, no. 5 (1972): 119–122."""
        """Wichert, E., Aziz, K., 1972. Calculation of Z’s for sour gases. Hydrocarb. Process. 51 (5), 119–122."""
        A=self.CO2+self.H2S
        e=120.*(A**0.9-A**1.6)+15*(self.H2S**0.5-self.H2S**4)
        tpc=self.tpc.R-e
        ppc=self.ppc.psi*tpc/(self.tpc.R+self.H2S*(1-self.H2S)*e)
        return unidades.Temperature(tpc, "R"), unidades.Pressure(ppc, "psi")

    def Critical_Carr_Kobayashi_Burrows(self):
        """Carr, N., Kobayashi, R., Burrows, D., 1954. Viscosity of hydrocarbon gases under pres-
        sure. Trans. AIME 201, 270–275."""
        tpc=self.tpc.R-80*self.CO2+130*self.H2S-250*self.N2
        ppc=self.ppc.psi+440*self.CO2+600*self.H2S-170*self.N2
        return unidades.Temperature(tpc, "R"), unidades.Pressure(ppc, "psi")

    def Critical_Whitson_Brule(self):
        """Whitson, C. H., and M. R. Brule. Phase Behavior. Richardson, TX: Society of Petroleum Engineers, 2000."""
        CO2=Componente(49)
        H2S=Componente(50)
        N2=Componente(46)
        g=(28.96*self.SG-(N2.M*self.N2+CO2.M*self.CO2+H2S.M*self.H2S))/28.96/(1-self.N2-self.CO2-self.H2S)
        tpcHC=168.+325.*g-12.5*g**2
        ppcHC=677+15.*g-37.5*g**2
        tpc=(1-self.N2-self.CO2-self.H2S)*tpcHC+N2.tc.R*self.N2+CO2.tc.R*self.CO2+H2S.tc.R*self.H2S
        ppc=(1-self.N2-self.CO2-self.H2S)*ppcHC+N2.pc.psi*self.N2+CO2.pc.psi*self.CO2+H2S.pc.psi*self.H2S
        return unidades.Temperature(tpc, "R"), unidades.Pressure(ppc, "psi")

    def tc_gas(self):
        """Cálculo de la temperatura crítica de un gas natural de composición conocida con metano como componente principal, API procedure 4C1.1 pag 329"""
        Tc1=exp(-5.624853)*exp(-0.0105852*self.MABP().R-1.4401126*self.SG+0.013200830*self.SG*self.MABP().R)
        Tc2=self.MABP().R**2.4289880
        Tc3=self.SG**-0.299808
        return unidades.Temperature(Tc1*Tc2*Tc3, "R")

    def Z_factor(self, T, P, Z_method=0):
        """Calculo del factor de compresibilidad del gas
        T: temperatura, en kelvin
        P: presión, en atm
        Z_method: método de cálculo del factor de compresibilidad:
            0   -   Hall_Yarborough
            1   -   Dranchuk_Abu_Kassem
            2   -   Dranchuk_Purvis_Robinson
            3   -   Shell-Oil Company
            4   -   Beggs_Brill
            5   -   Sarem
            6   -   Gopal
            7   -   Papay
        """
        Z=[Z_Hall_Yarborough, Z_Dranchuk_Abu_Kassem, Z_Dranchuk_Purvis_Robinson, Z_ShellOil, Z_Beggs_Brill, Z_Sarem, Z_Gopal, Z_Papay][Z_method]
        return Z(T/self.tpc, P/self.ppc.atm)

    def RhoG(self, T, P):
        Z=self.Z_factor(T, P)
        return unidades.Density(P*self.M/Z/R_atml/T, "gl")

    def Compressibility_Mattar_Brar_Aziz(self, T, P):
        """Mattar, L. G., S. Brar, and K. Aziz. “Compressibility of Natural Gases.” Journal of Canadian Petroleum Technology (October–November 1975): 77–80."""
        Tr=T/self.tpc
        Z=self.Z_factor(T, P)
        g=0.27*P/self.ppc.atm/Tr/Z
        T1=0.31506237-1.0467099/Tr-0.5783272/Tr**3
        T2=0.53530771-0.61232032/Tr
        T3=0.61232032*0.10488813/Tr
        T4=0.68157001/Tr**3
        T5=0.27*Pr/Tr
        dZ=T1+2*T2*g+5*T3*g**4+2*T4*g*(1+0.68446549*g**2-0.68446549**2*g**4)*exp(-0.68446549*g**2)-T5/g
        return self.ppc.atm/P-0.27/Z**2/Tr*dz/(1+g/Z*dz)

    def Gas_Formation_Volume_Factor(self, T, P):
        return Z*T/288.9/P

    def Viscosity_Carr_Kobayashi_Burrows(self, T):
        """Carr, N., R. Kobayashi, and D. Burrows. “Viscosity of Hydrocarbon Gases under Pressure.” Transactions of the AIME 201 (1954): 270–275."""
        muo=8.118e-3-6.15e-3*log10(self.SG)+(1.709e-5-2.062e-6*self.SG)*unidades.Temperature(T, "R").F
        muN2=self.N2*(8.49e-3*log10(self.SG)+9.59e-3)
        muCO2=self.CO2*(9.08e-3*log10(self.SG)+6.24e-3)
        muH2S=self.H2S*(8.49e-3*log10(self.SG)+3.73e-3)
        return unidades.Viscosity(muo+muN2+muCO2+muH2S, "cP")

    def Viscosity_Lee_Gonzalez_Eakin(self, T, P):
        """Lee, A. L., M. H. Gonzalez, and B. E. Eakin. “The Viscosity of Natural Gases.” Journal of Petroleum Technology (August 1966): 997–1000."""
        t=unidades.Temperature(T).R
        K=(9.4+0.02*self.M)*t**1.5/(209+19*self.M+t)
        X=3.5+986/t+0.01*self.M
        Y=2.4-0.2*X
        return unidades.Viscosity(1e-4*K*exp(X*(self.RhoG(T, P).lbft3/62.4)**Y), "cP")



if __name__ == '__main__':
    Crudo(index=1)
