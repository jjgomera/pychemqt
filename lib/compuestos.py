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


import math
import os
import re
import tempfile

from scipy import exp, cosh, sinh, log, log10, roots, absolute, sqrt
from scipy.optimize import fsolve
from scipy.constants import R, Avogadro

from lib.physics import R_atml, R_Btu, R_cal, factor_acentrico_octano
from lib import unidades, config, eos, sql


__doi__ = {
    1:
        {"autor": "Poling, B.E, Prausnitz, J.M, O'Connell, J.P",
         "title": "The Properties of Gases and Liquids 5th Edition",
         "ref": "McGraw-Hill, New York, 2001",
         "doi": ""},
    2:
        {"autor": "Tarek Ahmed",
         "title": "Equations of State and PVT Analysis: Applications for"
                  "Improved Reservoir Modeling, 2nd Edition",
         "ref": "Gulf Professional Publishing, 2016, ISBN 9780128015704,",
         "doi": "10.1016/B978-0-12-801570-4.00002-7"},
    3:
        {"autor": "Antoine, C.",
         "title": "Tensions des Vapeurs: Nouvelle Relation Entre les Tensions "
                  "et les Tempé",
         "ref": "Compt.Rend. 107:681-684 (1888)",
         "doi": ""},
    4:
        {"autor": "Lee, B. I. and Kesler, M. G.",
         "title": "A Generalized Thermodynamic Correlation Based on"
                  "Three-Parameter Corresponding States",
         "ref": "American Institute of Chemical Engineers Journal, 21, 1975",
         "doi": "10.1002/aic.690210313"},
    5:
        {"autor": "API",
         "title": "Technical Data book: Petroleum Refining 6th Edition",
         "ref": "",
         "doi": ""},
    6:
        {"autor": "Wagner, W.",
         "title": "New Vapour Pressure Measurements for Argon and Nitrogen "
                  "and a New Method for Establishing Rational Vapour Pressure "
                  "Equations",
         "ref": "Cryogenics 13, 8 (1973) 470-82",
         "doi": "10.1016/0011-2275(73)90003-9"},
    7:
        {"autor": "McGarry, J.",
         "title": "Correlation and Perediction of the Vapor Pressures of Pure"
                  "Liquids over Large Pressure Ranges",
         "ref": "Ind. Eng. Chem. Process. Des. Dev. 22 (1983) 313-322",
         "doi": "10.1021/i200021a023"},
    8:
        {"autor": "Ambrose, D., Walton, J.",
         "title": "Vapour Pressures up to Their Critical Temperatures of "
                  "Normal Alkanes and 1-Alkanols",
         "ref": "Pure & Appl. Chem. 61(8) 1395-1403 (1989)",
         "doi": "10.1351/pac198961081395"},
    9:
        {"autor": "Sanjari, E., Honarmand, M., Badihi, H., Ghaheri, A.",
         "title": "An Accurate Generalized Model for Predict Vapor Pressure "
                  "of Refrigerants",
         "ref": "International Journal of Refrigeration 36 (2013) 1327-1332",
         "doi": "10.1016/j.ijrefrig.2013.01.007"},
    10:
        {"autor": "Letsou, A., Stiel, L.I.",
         "title": "Viscosity of Saturated Nonpolar Liquids at Elevated "
                  "Pressures",
         "ref": "AIChE Journal 19(2) (1973) 409-411",
         "doi": "10.1002/aic.690190241"},
    11:
        {"autor": "Riazi, M.R., Faghri, A.",
         "title": "Thermal Conductivity of Liquid and Vapor Hydrocarbon "
                  "Systems: Pentanes and Heavier at Low Pressures",
         "ref": "Ind. Eng. Chem. Process Des. Dev. 24 (1985) 398-401",
         "doi": "10.1021/i200029a030"},
    12:
        {"autor": "Gharagheizi, F., Ilani-Kashkouli, P., Sattari, M., "
                  "Mohammadi, A.H., Ramjugernath, D., Richon, D.",
         "title": "Development of a General Model for Determination of "
                  "Thermal Conductivity of Liquid Chemical Compounds at "
                  "Atmospheric Pressure",
         "ref": "AIChE Journal 59 (2013) 1702-1708",
         "doi": "10.1002/aic.13938"},
    13:
        {"autor": "Lakshmi, D.S., Prasad, D.H.L.",
         "title": "A Rapid Estimation Method for Thermal Conductivity of Pure "
                  "Liquids",
         "ref": "The Chemical Engineering Journal 48 (1992) 211-14",
         "doi": "10.1016/0300-9467(92)80037-B"},
    14:
        {"autor": "Di Nicola, G., Ciarrocchi, E., Coccia, G., Pierantozzi, M.",
         "title": "Correlations of Thermal Conductivity for Liquid "
                  "Refrigerants at Atmospheric Pressure or near Saturation",
         "ref": "International Journal of Refrigeration, 2014",
         "doi": "10.1016/j.ijrefrig.2014.06.003"},
    15:
        {"autor": "Pachaiyappan, V., Ibrahim, S.H., Kuloor, N.R.",
         "title": "Thermal Conductivities of Organic Liquids: A New "
                  "Correlation",
         "ref": "J. Chem. Eng. Data, 11 (1966) 73-76",
         "doi": "10.1021/je60028a021"},
    16:
        {"autor": "Kanitkar, D., Thodos, G.",
         "title": "The Thermal Conductivity of Liquid Hydrocarbons",
         "ref": "Can. J. Chem. Eng. 47 (1969) 427-430",
         "doi": "10.1002/cjce.5450470502"},
    17:
        {"autor": "Riedel, L., Chem. Ingr. Tech., 26 (1954): 679.",
         "title": "",
         "ref": "",
         "doi": ""},
    18:
        {"autor": "Riedel, L.",
         "title": "Die Zustandsfunktion des realen Gases: Untersuchungen über "
                  "eine Erweiterung des Theorems der übereinstimmenden "
                  "Zustände",
         "ref": "Chem. Ings-Tech. 28 (1956) 557-562",
         "doi": "10.1002/cite.330280809"},
    19:
        {"autor": "Edmister, W.C.",
         "title": "Applied Hydrocarbon Thermodynamics, Part 4, Compressibility"
                  "Factors and Equations of State",
         "ref": "Petroleum Refiner. 37 (April, 1958), 173–179",
         "doi": ""},
    20:
        {"autor": "Hankinson, R.W., Thomson, G.H.",
         "title": "A New Correlation for Saturated Densities of Liquids and "
                  "Their Mixtures",
         "ref": "AIChE Journal 25(4) (1979) 653-663",
         "doi": "10.1002/aic.690250412"},
    21:
        {"autor": "Rackett, H.G.",
         "title": "Equation of State for Saturated Liquids",
         "ref": "J. Chem. Eng. Data 15(4) (1970) 514-517",
         "doi": "10.1021/je60047a012"},
    22:
        {"autor": "Thomson, G.H., Brobst, K.R., Hankinson, R.W.",
         "title": "An Improved Correlation for Densities of Compressed Liquids"
                  " and Liquid Mixtures",
         "ref": "AIChE Journal 28(4) (1982): 671-76",
         "doi": "10.1002/aic.690280420"},
    23:
        {"autor": "Yamada, T., Gunn. R.",
         "title": "Saturated Liquid Molar Volumes: The Rackett Equation",
         "ref": "Journal of Chemical Engineering Data 18(2) (1973): 234–236",
         "doi": "10.1021/je60057a006"},
    24:
        {"autor": "Stiel, L. I., Thodos, G.",
         "title": "The Viscosity of Nonpolar Gases at Normal Pressures",
         "ref": "AIChE J. 7(4) (1961) 611-615",
         "doi": "10.1002/aic.690070416"},




    25:
        {"autor": "",
         "title": "",
         "ref": "",
         "doi": ""},
}


def atomic_decomposition(cmp):
    """Procedure to decompose a molecular string representation in its atomic
    composition

    Parameters
    ------------
    cmp : string
        Compound test representation

    Returns
    -------
    kw : dict
        Dictionary with the atomic decomposition of compound

    Examples
    --------
    >>> kw = atomic_decomposition("CH3")
    >>> "%i %i" % (kw["C"], kw["H"])
    '1 3'
    >>> kw = atomic_decomposition("COO")
    >>> "%i %i" % (kw["C"], kw["O"])
    '1 2'
    >>> kw = atomic_decomposition("CH3COOCl")
    >>> "%i %i %i %i" % (kw["C"], kw["H"], kw["O"], kw["Cl"])
    '2 3 2 1'
    """
    kw = {}
    for element, c in re.findall("([A-Z][a-z]*)([0-9]*)", cmp):
        if element not in kw:
            kw[element] = 0

        if not c:
            c = 1
        else:
            c = int(c)
        kw[element] += c
    return kw


def DIPPR(prop, T, args, Tc=None, M=None):
    """Procedure to implement the DIPPR equations valid to calculate several
    physical properties of compounds. The parameters of the

    Parameters
    ----------
    prop : string
        Property to calculate, any of:
        rhoS, rhoL, rhoG, Hv, cpS, cpL, cpG, muL, muG, kL, kG, sigma
    T : float
        Temperature, [K]
    args : list
        Coefficients for DIPPR equation, [eq, A, B, C, D, E, Tmin, Tmax]
    Tc : float, optional
        Critical temperature, [K]
    M : float, optional
        Molecular weight, [g/mol]

    Notes
    -----
    The properties this method can calculate, and the units for the calculated
    properties are:

        -rhoS: Solid density, [kmol/m³]
        -rhoL: Liquid density, [kmol/m³]
        -rhoG: Vapor density, [kmol/m³]
        -Pv: Vapor pressure, [Pa]
        -Hv: Heat of vaporization, [J/kmol]
        -cpS: Solid heat capacity, [J/kmol·K]
        -cpL: Liquid heat capacity, [J/kmol·K]
        -cpG: Ideal gas heat capacity, [J/kmol·K]
        -muL: Liquid viscosity, [Pa·s]
        -muG: Vapor viscosity, [Pa·s]
        -kL: Liquid thermal conductivity, [W/m·K]
        -kG: Vapor thermal conductivity, [W/m·K]
        -sigma: Surface Tension, [N/m]

    The first element in args define the equation to use:

        Eq 1:   Y = A+B*T+C*T^2+D*T^3+E*T^4
        Eq 2:   Y = exp(A+B*T+C*ln(T)+D*T^E)
        Eq 3:   Y = A*T^B/(1+C*T+D*T^2)
        Eq 4:   Y = A+B*exp(-C/T^D)
        Eq 5:   Y = A + B*T + C*T^3 + D*T^8 + E*T^9
        Eq 6:   Y = A/(B^(1+(1-T/C)^D)
        Eq 7:   Y = A*(1-Tr)^(B+C*Tr+D*Tr^2+E*Tr^3)
        Eq 8:   Y = A+ B*((C/T)/sinh(C/T))^2 + D*((E/T)/cosh(E/T))^2
        Eq 9:   Y = A²/Tr + B - 2ACTr - ADTr² - C²Tr³/3 - CDTr⁴/2 - D²Tr⁵/5

        where: T: Temperature, [K]
            Tr: Reduced temperature T/Tc
            A,B,C,D,E: Parameters of equation

    This parameters are available in the pychemqt database for many compounds
    Some equation as 7 and 9 need aditional parameter Tc of compound
    """
    # Multiplier for return the properties in mass base
    mul = 1

    # Select unit of property to return
    if "rho" in prop:
        unit = unidades.Density
        mul = M
    elif prop == "Pv":
        unit = unidades.Pressure
    elif prop == "Hv":
        unit = unidades.MolarEnthalpy
        mul = M
    elif "cp" in prop:
        unit = unidades.SpecificHeat
        mul = M
    elif "mu" in prop:
        unit = unidades.Viscosity
    elif "k" in prop:
        unit = unidades.ThermalConductivity
    elif prop == "sigma":
        unit = unidades.Tension

    eq, A, B, C, D, E, Tmin, Tmax = args
    if eq == 1:
        value = A + B*T + C*T**2 + D*T**3 + E*T**4
    elif eq == 2:
        value = exp(A + B/T + C*log(T) + D*T**E)
    elif eq == 3:
        value = A*T**B/(1+C/T+D/T**2)
    elif eq == 4:
        value = A + B*exp(-C/T**D)
    elif eq == 5:
        value = A + B/T + C/T**3 + D/T**8 + E/T**9
    elif eq == 6:
        value = A/(B**(1+((1-T/C)**D)))
    elif eq == 7:
        Tr = T/Tc
        value = A*(1-Tr)**(B+C*Tr + D*Tr**2 + E*Tr**3)
    elif eq == 8:
        value = A + B*(C/T/sinh(C/T))**2 + D*(E/T/cosh(E/T))**2
    elif eq == 9:
        Tr = T/Tc
        value = A**2/Tr + B - 2*A*C*Tr - A*D*Tr**2 - C**2*Tr**3/3 - \
            C*D*Tr**4/2 - D**2*Tr**5/5

    return unit(value*mul)


# Liquid density correlations
def RhoL_Rackett(T, Tc, Pc, Zra, M):
    """Calculates saturated liquid densities of pure components using the
    modified Rackett equation, referenced too in API procedure 6A2.13 pag. 454

    .. math::
        \frac{1}{\rho_s} = \frac{RT_c}{P_c}Z_{RA}^{1+(1-{T/T_c})^{2/7}}

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tc : float
        Critical temperature, [K]
    Pc : float
        Critical pressure, [Pa]
    Zra : float
        Racket constant, [-]
    M : float
        Molecular weight, [kg/m³]

    Returns
    -------
    rho : float
        Saturated liquid density at T, [kg/m³]

    Examples
    --------
    Example from [5]_; propane at 30ºF

    >>> T = unidades.Temperature(30, "F")
    >>> Tc = unidades.Temperature(206.06, "F")
    >>> Pc = unidades.Pressure(616, "psi")
    >>> "%0.3f" % RhoL_Rackett(T, Tc, Pc, 0.2763, 44.1).kgl
    '0.531'

    References
    ----------
    .. [21] Rackett, H.G. Equation of State for Saturated Liquids. J. Chem.
        Eng. Data 15(4) (1970) 514-517
    .. [5] API. Technical Data book: Petroleum Refining 6th Edition
    """
    Pc_atm = Pc/101325
    Tr = T/Tc
    V = R_atml*Tc/Pc_atm*Zra**(1+(1-Tr)**(2/7))
    return unidades.Density(M/V)


def RhoL_Costald(T, Tc, w, Vc):
    """Calculates saturated liquid densities of pure components using the
    Corresponding STAtes Liquid Density (COSTALD) method, developed by
    Hankingon and Thomson, referenced too in API procedure 6A2.15 pag. 462

    .. math::
        \frac{V}{V^{o}}=V_{R}^{(0)}\left(1-\omega_{SRK}V_{R}^{(1)}\right)

        V_{R}^{(0)}=1+a\left(1-T_{R}\right)^{1/3}+b\left(1-T_{R}\right)^{2/3}
        +c\left(1-T_{R}\right)+d\left(1-T_{R}\right)^{4/3}

        V_{R}^{(1)}=\frac{e+fT_{R}+gT_{R}^{2}+hT_{R}^{3}}{Tr-1.00001}

    Parameters
    ----------
    T : float
        Temperature [K]
    Tc : float
        Critical temperature [K]
    w : float
        Acentric factor optimized to SRK, [-]
    Vc : float
        Characteristic volume, [m³/kg]

    Returns
    -------
    rho : float
        Saturated liquid density at T, [kg/m³]

    Examples
    --------
    Example 1 from [5]_; propane at 30ºF

    >>> T = unidades.Temperature(30, "F")
    >>> Tc = unidades.Temperature(206.01, "F")
    >>> Vc = unidades.SpecificVolume(3.205/44.097, "ft3lb")
    >>> "%0.3f" % RhoL_Costald(T, Tc, 0.1532, Vc).kgl
    '0.530'

    References
    ----------
    .. [20] Hankinson, R.W., Thomson, G.H. A New Correlation for Saturated
        Densities of Liquids and Their Mixtures. AIChE Journal 25(4) (1979)
        653-663
    .. [5] API. Technical Data book: Petroleum Refining 6th Edition
    """
    a = -1.52816
    b = 1.43907
    c = -0.81446
    d = 0.190454
    e = -0.296123
    f = 0.386914
    g = -0.0427258
    h = -0.0480645
    Tr = T/Tc

    # Eq 17
    Vr0 = 1 + a*(1-Tr)**(1/3) + b*(1-Tr)**(2/3) + c*(1-Tr) + d*(1-Tr)**(4/3)
    Vr1 = (e + f*Tr + g*Tr**2 + h*Tr**3)/(Tr-1.00001)                  # Eq 18

    # TODO: Add V* to the database
    V = Vc*Vr0*(1-w*Vr1)                                               # Eq 16
    return unidades.Density(1/V)


def RhoL_ThomsonBrobstHankinson(T, P, Tc, Pc, w, Ps, rhos):
    """Calculates compressed-liquid density, using the Thomson-Brobst-
    Hankinson correlation, also referenced in API procedure 6A2.23 pag. 477

    .. math::
        V = V_s\left(1-C\ln\frac{B + P}{B + P_s}\right)

        \frac{B}{P_c} = -1 + a\tau^{1/3} + b\tau^{2/3} + d\tau + e\tau^{4/3}

        e = \exp(f + g\omega_{SRK} + h \omega_{SRK}^2)

        C = j + k \omega_{SRK}

    Parameters
    ----------
    T : float
        Temperature, [K]
    P : float
        Pressure, [Pa]
    Tc : float
        Critical temperature, [K]
    Pc : float
        Critical pressure, [Pa]
    Ps : float
        Saturation pressure, [Pa]
    omega : float
        Acentric factor (SRK optimized), [-]
    rhos : float
        Saturation liquid volume, [kg/m^3]

    Returns
    -------
    rho : float
        High-pressure liquid density, [kg/m^3]

    Examples
    --------
    Example from [5]_; n-octane at 212ºF and 4410 psi

    >>> T = unidades.Temperature(212, "F")
    >>> P = unidades.Pressure(4410, "psi")
    >>> Tc = unidades.Temperature(564.22, "F")
    >>> Pc = unidades.Pressure(360.6, "psi")
    >>> Ps = unidades.Pressure(6.74, "psi")
    >>> rs = RhoL_Rackett(T, Tc, Pc, 0.2569, 114.232)
    >>> "%0.3f" % (1/rs.lbft3*114.232)
    '2.874'
    >>> "%0.3f" % RhoL_ThomsonBrobstHankinson(T, P, Tc, Pc, 0.3962, Ps, rs).kgl
    '0.676'

    References
    ----------
    .. [22] Thomson, G.H., Brobst, K.R., Hankinson, R.W. An Improved
        Correlation for Densities of Compressed Liquids and Liquid Mixtures.
        AIChE Journal 28(4) (1982): 671-76
    .. [5] API. Technical Data book: Petroleum Refining 6th Edition
    """
    Tr = T/Tc
    C = 0.0861488 + 0.0344483*w                                         # Eq 8
    e = exp(4.79594 + 0.250047*w + 1.14188*w**2)                        # Eq 7
    B = Pc*(-1 - 9.070217*(1-Tr)**(1/3) + 62.45326*(1-Tr)**(2/3) -
            135.1102*(1-Tr) + e*(1-Tr)**(4/3))

    # Eq 5
    rho = rhos/(1-C*log((B+P)/(B+Ps)))
    return unidades.Density(rho, "gl")


# Vapor pressure correlation
def Pv_Antoine(T, args, Tc=None, base=math.e, Punit="mmHg"):
    """Vapor Pressure calculation procedure using the Antoine equation

    .. math::
        \log_{\text{base}} P^{\text{sat}} = A - \frac{B}{T+C}

    The method implement too the extended Antoine Equation

    .. math::
        \log_{10} P^{sat} = A - \frac{B}{T + C} + 0.43429x^n + Ex^8 + Fx^{12}

        x = \max \left(\frac{T-t_o-273.15}{T_c}, 0 \right)

    Parameters
    ----------
    T : float
        Temperature of fluid, [K]
    args : list
        Coefficients for Antoine equation
    Tc : float, optional
        Critical temperature, [K]
    base : float, optional
        The base of logarithm in equation, default e
    Punit : string, optional
        Code of pressure unit calculated

    Returns
    -------
    Pv : float
        Vapor pressure, [Pa]

    Notes
    -----
    The length of args define the method to use, if args has three elements use
    the original version, if has seven element and define Tc use the advanced
    method.
    The coefficient of equation saved in database are for pressure in mmHg and
    with a exponential dependence. If it defines parameters for a new component
    it can configure this values, the saved equation will be converted to the
    appropiate format in database

    Examples
    --------
    Example 7-1 in [1]_, furan at 309.429 K
    >>> P = Pv_Antoine(309.429, (4.1199, 1070.2, -44.32), base=10, Punit="bar")
    >>> "%0.4f" % P.bar
    '1.2108'

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    .. [2] Antoine, C. 1888. Tensions des Vapeurs: Nouvelle Relation Entre les
       Tensions et les Tempé. Compt.Rend. 107:681-684.
    """
    if len(args) == 3:
        A, B, C = args
        Pv = base**(A-B/(T+C))
    elif len(args) == 7 and Tc is not None:
        A, B, C, n, E, F, to = args
        x = (T-to)/Tc
        if x <= 0:
            Pv = base**(A-B/(T+C))
        else:
            Pv = base**(A - B/(T+C) + 0.43429*x**n + E*x**8 + F*x**12)
    return unidades.Pressure(Pv, Punit)


def Pv_Lee_Kesler(T, Tc, Pc, w):
    """Calculates vapor pressure of a fluid using the Lee-Kesler correlation

    The vapor pressure is given by:

    .. math::
        \ln P_r = f^{(0)} + \omega f^{(1)}

        f^{(0)} = 5.92714-\frac{6.09648}{T_r}-1.28862\ln T_r + 0.169347T_r^6

        f^{(1)} = 15.2518-\frac{15.6875}{T_r} - 13.4721 \ln T_r + 0.43577T_r^6

    Parameters
    ----------
    T : float
        Temperature [K]
    Tc : float
        Critical temperature [K]
    Pc : float
        Critical pressure [Pa]
    w : float
        Acentric factor [-]

    Returns
    -------
    Pv : float
        Vapor pressure at T [Pa]

    Examples
    --------
    Example 1.2 from [4]_; propane at 80ºF

    >>> T = unidades.Temperature(80, "F")
    >>> Tc = unidades.Temperature(666.01, "R")
    >>> Pc = unidades.Pressure(616.3, "psi")
    >>> "%0.0f" % Pv_Lee_Kesler(T, Tc, Pc, 0.1522).psi
    '144'

    References
    ----------
    .. [4] Lee, B. I. and Kesler, M. G., A Generalized Thermodynamic
       Correlation Based on Three-Parameter Corresponding States. American
       Institute of Chemical Engineers Journal, Vot. 21, 1975
    .. [2] Tarek Ahmed. Equations of State and PVT Analysis: Applications for
       Improved Reservoir Modeling, 2nd Edition. Gulf Professional Publishing,
       2016, ISBN 9780128015704
    """
    # Eq 17, pag 525
    Tr = T/Tc
    f0 = 5.92714 - 6.09648/Tr - 1.28862*log(Tr) + 0.169347*Tr**6
    f1 = 15.2518 - 15.6875/Tr - 13.4721*log(Tr) + 0.43577*Tr**6
    return unidades.Pressure(exp(f0 + w*f1)*Pc)


def Pv_Wagner(T, args, Tc, Pc):
    """Calculates vapor pressure of a fluid using the Wagner correlation

    .. math::
        \ln P^{v}= \ln P_c + \frac{a\tau + b \tau^{1.5} + c\tau^{3}
        + d\tau^6} {T_r}

        \tau = 1 - \frac{T}{T_c}

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tc : float
        Critical temperature, [K]
    Pc : float
        Critical pressure, [Pa]
    args : list
        Coefficients for equation

    Returns
    -------
    Pv : float
        Vapor pressure, [Pa]

    Notes
    -----
    Same compound has the parameters of this equations saved in database. This
    method implement the origintal form of Wagner as in [6]_, with the
    parameters from McGarry. API use other same different form.

    References
    ----------
    .. [6] Wagner, W. New Vapour Pressure Measurements for Argon and Nitrogen
        and a New Method for Establishing Rational Vapour Pressure Equations.
        Cryogenics 13, 8 (1973) 470-82.
    .. [2] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
        New York: McGraw-Hill Professional, 2000.
    .. [5] API. Technical Data book: Petroleum Refining 6th Edition
    .. [7] McGarry, J. Correlation and Perediction of the Vapor Pressures of
        Pure Liquids over Large Pressure Ranges. Ind. Eng. Chem. Process. Des.
        Dev. 22 (1983) 313-322
    """
    a, b, c, d = args
    Tr = T/Tc
    tau = 1-Tr
    Pv = Pc/Tr*exp(a*tau + b*tau**1.5 + c*tau**3 + d*tau**6)            # Eq 14
    return unidades.Pressure(Pv)


def Pv_AmbroseWalton(T, Tc, Pc, w):
    """Calculates vapor pressure of a fluid using the Ambrose-Walton
    corresponding-states correlation

    .. math::
        \ln P_r=f^{(0)}+\omegaf^{(1)}+\omega^2f^{(2)}

        f^{(0)}=\frac{-5.97616\tau + 1.29874\tau^{1.5}- 0.60394\tau^{2.5}
        -1.06841\tau^5}{T_r}

        f^{(1)}=\frac{-5.03365\tau + 1.11505\tau^{1.5}- 5.41217\tau^{2.5}
        -7.46628\tau^5}{T_r}

        f^{(2)}=\frac{-0.64771\tau + 2.41539\tau^{1.5}- 4.26979\tau^{2.5}
        +3.25259\tau^5}{T_r}

        \tau = 1-T_{r}

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tc : float
        Critical temperature, [K]
    Pc : float
        Critical pressure, [Pa]
    w : float
        Acentric factor, [-]

    Returns
    -------
    Pv : float
        Vapor pressure, [Pa]

    Examples
    --------
    Example 7-3 from [2]_; ethylbenzene at 347.25 K.

    >>> "%0.4f" % Pv_AmbroseWalton(347.25, 617.15, 36.09E5, 0.304).bar
    '0.1328'
    >>> "%0.3f" % Pv_AmbroseWalton(460, 617.15, 36.09E5, 0.304).bar
    '3.325'

    References
    ----------
    .. [8] Ambrose, D., Walton, J. Vapour Pressures up to Their Critical
        Temperatures of Normal Alkanes and 1-Alkanols. Pure & Appl. Chem. 61(8)
        1395-1403 (1989)
    .. [2] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    """
    Tr = T/Tc
    t = 1 - T/Tc
    f0 = (-5.97616*t + 1.29874*t**1.5 - 0.60394*t**2.5 - 1.06841*t**5)/Tr
    f1 = (-5.03365*t + 1.11505*t**1.5 - 5.41217*t**2.5 - 7.46628*t**5)/Tr
    f2 = (-0.64771*t + 2.41539*t**1.5 - 4.26979*t**2.5 + 3.25259*t**5)/Tr

    # Eq 8
    Pv = Pc*exp(f0 + w*f1 + w**2*f2)
    return unidades.Pressure(Pv)


def Pv_Riedel(T, Tc, Pc, Tb):
    """Calculate vapor pressure of a fluid using the Rieel corresponding-states
    correlation

    .. math::
        \ln P_{\text{vp}} = A - \frac{B}{T} + C\ln T + DT^{6}

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tc : float
        Critical temperature, [K]
    Pc : float
        Critical pressure, [Pa]
    Tb : float
        Normal boiling temperature, [K]

    Returns
    -------
    Pv : float
        Vapor pressure, [Pa]

    Examples
    --------
    Example 7-4 from [2]_; ethylbenzene

    >>> "%0.3f" % Pv_Riedel(347.25, 617.15, 36.09E5, 409.36).bar
    '0.131'
    >>> "%0.2f" % Pv_Riedel(460, 617.15, 36.09E5, 409.36).bar
    '3.35'

    References
    ----------
    .. [2] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    """
    Tr = T/Tc
    Tbr = Tb/Tc

    # TODO: Define the K values by vetere dependences of chemical
    K = 0.0838

    fib = -35 + 36/Tbr + 42*log(Tbr) - Tbr**6
    alfa_c = (3.758*K*fib+log(Pc*1e-5/1.01325))/(K*fib-log(Tbr))
    Q = K*(3.758-alfa_c)
    Pv = Pc*exp(-35*Q + 36*Q/Tr + (42*Q+alfa_c)*log(Tr) - Q*Tr**6)
    return unidades.Pressure(Pv)


def Pv_MaxwellBonnel(T, Tb, Kw):
    """Calculates vapor pressure of a fluid using the Maxell-Bonnel correlation
    as explain in [5]_, procedure 5A1.18, Pag. 394

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tb : float
        Normal boiling temperature, [K]
    Kw : float
        Watson factor, [-]

    Returns
    -------
    Pv : float
        Vapor pressure, [Pa]

    Notes
    -----
    This method isn't recomended, only when none of other methods are
    applicable

    Examples
    --------
    Example in [5]_, tetralin at 302ºF
    >>> T = unidades.Temperature(302, "F")
    >>> Tb = unidades.Temperature(405.7, "F")
    >>> Pv = Pv_MaxwellBonnel(T, Tb, 9.78)
    >>> "%0.1f" % Pv.psi
    '3.1'

    References
    ----------
    .. [5] API. Technical Data book: Petroleum Refining 6th Edition
    """
    # Convert input Tb in Kelvin to Fahrenheit to use in the correlation
    Tb_F = unidades.K2F(Tb)
    Tb_R = unidades.K2R(Tb)
    T_R = unidades.K2R(T)

    if Tb_F > 400:
        f = 1.0
    elif Tb_F < 200:
        f = 0.0
    else:
        f = (Tb_R-659.7)/200

    def P(Tb):
        X = (Tb/T_R-0.0002867*Tb)/(748.1-0.2145*Tb)
        if X > 0.0022:
            p = 10**((3000.538*X-6.761560)/(43*X-0.987672))
        elif X < 0.0013:
            p = 10**((2770.085*X-6.412631)/(36*X-0.989679))
        else:
            p = 10**((2663.129*X-5.994296)/(95.76*X-0.972546))
        return p

    Tb = fsolve(lambda Tb: Tb-Tb_R+2.5*f*(Kw-12)*log10(P(Tb)/760), Tb)
    p = P(Tb)
    return unidades.Pressure(p, "mmHg")


def Pv_Sanjari(T, Tc, Pc, w):
    """Calculates vapor pressure of a fluid using the Sanjari correlation
    pressure, and acentric factor.

    The vapor pressure of a chemical at `T` is given by:

    .. math::
        P_{v} = P_c\exp(f^{(0)} + \omegaf^{(1)} + \omega^2f^{(2)})

        f^{(0)} = a_1 + \frac{a_2}{T_r} + a_3\ln T_r + a_4 T_r^{1.9}

        f^{(1)} = a_5 + \frac{a_6}{T_r} + a_7\ln T_r + a_8 T_r^{1.9}

        f^{(2)} = a_9 + \frac{a_{10}}{T_r} + a_{11}\ln T_r + a_{12} T_r^{1.9}

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tc : float
        Critical temperature, [K]
    Pc : float
        Critical pressure, [Pa]
    w : float
        Acentric factor, [-]

    Returns
    -------
    Pv : float
        Vapor pressure, [Pa]

    Notes
    -----
    This method have been developed fitting data of refrigerants, be careful
    when use with other type of compound.

    References
    ----------
    .. [9] Sanjari, E., Honarmand, M., Badihi, H., Ghaheri, A. An Accurate
        Generalized Model for Predict Vapor Pressure of Refrigerants
        International Journal of Refrigeration 36 (2013) 1327-1332
    """
    # Table 2
    a = [None, 6.83377, -5.76051, 0.90654, -1.16906, 5.32034, -28.1460,
         -58.0352, 23.57466, 18.19967, 16.33839, 65.6995, -35.9739]

    Tr = T/Tc
    f0 = a[1] + a[2]/Tr + a[3]*log(Tr) + a[4]*Tr**1.9
    f1 = a[5] + a[6]/Tr + a[7]*log(Tr) + a[8]*Tr**1.9
    f2 = a[9] + a[10]/Tr + a[11]*log(Tr) + a[12]*Tr**1.9
    Pv = Pc*exp(f0 + w*f1 + w**2*f2)
    return unidades.Pressure(Pv)


def Tension_Parametric(T, args, Tc):
    """Calculates surface tension of fluid using a paremtric equation

    .. math::
        $\sigma=A\left(1-T_{r}\right)^{B}$

        Tr = \frac{T}{T_c}

    Parameters
    ----------
    T : float
        Temperature, [K]
    args : list
        Coefficients for equation
    Tc : float
        Critical temperature, [K]

    Returns
    -------
    sigma : float
        Surface tension, [N/m]

    Notes
    -----
    The parameters for several compound are in database
    """
    A, B = args
    Tr = T/Tc
    sigma = A*(1-Tr)**B
    return unidades.Tension(sigma)


def MuL_Parametric(T, args):
    """Calculates liquid viscosity using a paremtric equation

    .. math::
        \log\mu = A\left(\frac{1}{T}-\frac{1}{B}\right)

    Parameters
    ----------
    T : float
        Temperature, [K]
    args : list
        Coefficients for equation

    Returns
    -------
    mu : float
        Liquid viscosity, [Pa·s]

    Notes
    -----
    The parameters for several compound are in database
    """
    A, B = args
    mu = 10**(A*(1/T-1/B))
    return unidades.Viscosity(mu, "cP")


def MuL_LetsouStiel(T, M, Tc, Pc, w):
    r"""Calculate the viscosity of a liquid using the Letsou-Stiel correlation

    .. math::
        \mu = (\xi^{(0)} + \omega \xi^{(1)})/\xi

    Parameters
    ----------
    T : float
        Temperature, [K]
    M : float
        Molecular weight, [g/mol]
    Tc : float
        Critical temperature, [K]
    Pc : float
        Critical pressure, [Pa]
    w : float
        Acentric factor, [-]

    Returns
    -------
    mu : float
        Viscosity, [Pa·s]

    References
    ----------
    .. [10] Letsou, A., Stiel, L.I. Viscosity of Saturated Nonpolar Liquids at
        Elevated Pressures. AIChE Journal 19(2) (1973) 409-411
    """
    Pc_atm = unidades.Pressure(Pc).atm
    Tr = T/Tc

    x = Tc**(1/6)/M**0.5/Pc_atm**(2/3)
    x0 = 0.015178 - 0.021351*Tr + 0.007503*Tr**2
    x1 = 0.042559 - 0.07675*Tr + 0.034007*Tr**2
    mu = (x0+w*x1)/x
    return unidades.Viscosity(mu, "cP")


def MuG_StielThodos(T, Tc, Pc, M):
    r"""Calculate the viscosity of a gas using the Stiel-Thodos correlation,
    also referenced in API procedure 11B1.3, pag 1099

    .. math::
        \mu=N/\xi

        \xi=\frac{T_{c}^{1/6}}{M^{1/2}P_{c}^{2/3}}

        N=3.4e^{-4}T_{r}^{0.94}   for Tr ≤ 1.5

        N=1.778e^{-4}\left(4.58T_{r}-1.67\right)^{0.625} for T_r > 1.5

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tc : float
        Critical temperature, [K]
    Pc : float
        Critical pressure, [Pa]
    M : float
        Molecular weight, [g/mol]

    Returns
    -------
    mu : float
        Viscosity of gas, [Pa·s]

    Examples
    --------
    Example A in [5]_, Propane at 176ºF
    >>> T = unidades.Temperature(176, "F")
    >>> Tc = unidades.Temperature(206, "F")
    >>> Pc = unidades.Pressure(616, "psi")
    >>> "%0.4f" % MuG_StielThodos(T, Tc, Pc, 44.1).cP
    '0.0100'

    Example B in [5]_, Methane at 543ºF
    >>> T = unidades.Temperature(543, "F")
    >>> Tc = unidades.Temperature(-116.67, "F")
    >>> Pc = unidades.Pressure(667, "psi")
    >>> "%0.4f" % MuG_StielThodos(T, Tc, Pc, 16.04).cP
    '0.0176'

    References
    ----------
    .. [24] Stiel, L. I., Thodos, G. The Viscosity of Nonpolar Gases at Normal
        Pressures. AIChE J. 7(4) (1961) 611-615
    .. [5] API. Technical Data book: Petroleum Refining 6th Edition
    """
    Pc_atm = Pc/101325
    Tr = T/Tc
    T_R = unidades.K2R(T)

    if M < 2:
        # Special case for hydrogen
        if Tr <= 1.5:
            mu = 3.7e-5*T_R**0.94
        else:
            mu = 9.071e-4*(7.639e-2*T_R-1.67)**0.625
    else:
        if Tr <= 1.5:
            N = 3.5e-4*Tr**0.94
        else:
            N = 1.778e-4*(4.58*Tr-1.67)**0.625
        x = Tc**(1/6)/M**0.5/Pc_atm**(2/3)
        mu = N/x
    return unidades.Viscosity(mu, "cP")


def facent_LeeKesler(Tb, Tc, Pc):
    """Calculates acentric factor of a fluid using the Lee-Kesler correlation

    Parameters
    ----------
    Tb : float
        Boiling temperature [K]
    Tc : float
        Critical temperature [K]
    Pc : float
        Critical pressure [Pa]

    Returns
    -------
    w : float
        Acentric factor [-]

    References
    ----------
    .. [4] Lee, B. I. and Kesler, M. G., A Generalized Thermodynamic
       Correlation Based on Three-Parameter Corresponding States. American
       Institute of Chemical Engineers Journal, Vot. 21, 1975
    """
    Tr = Tb/Tc
    Pr = 101325/Pc
    w = (log(Pr) - 5.92714 + 6.09648/Tr + 1.28862*log(Tr) - 0.169347*Tr**6)/(
        15.2518 - 15.6875/Tr - 13.4721*log(Tr) + 0.43577*Tr**6)

    return unidades.Dimensionless(w)


def prop_Edmister(**kwargs):
    """Calculate the missing parameters between Tc, Pc, Tb and acentric factor
    from the Edmister (1958) correlations

    Parameters
    ------------
    Tc : float
        Critic temperature, [ºR]
    Pc : float
        Critic pressure, [psi]
    Tb : float
        Boiling temperature, [ºR]
    w : float
        Acentric factor, [-]

    Returns
    -------
    prop : Dict with the input parameter and the missing parameter in input

    References
    ----------
    [19] .. Edmister, W. C. Applied Hydrocarbon Thermodynamics, Part 4,
        Compressibility Factors and Equations of State. Petroleum Refiner 37
        (April 1958): 173–179.
    """
    count_available = 0
    if "Tc" in kwargs and kwargs["Tc"]:
        count_available += 1
        Tc = unidades.Temperature(kwargs["Tc"])
    else:
        unknown = "Tc"

    if "Pc" in kwargs and kwargs["Pc"]:
        count_available += 1
        Pc = unidades.Pressure(kwargs["Pc"])
    else:
        unknown = "Pc"

    if "Tb" in kwargs and kwargs["Tb"]:
        count_available += 1
        Tb = unidades.Temperature(kwargs["Tb"])
    else:
        unknown = "Tb"
    if "w" in kwargs:
        count_available += 1
        w = unidades.Dimensionless(kwargs["w"])
    else:
        unknown = "w"

    if count_available != 3:
        raise ValueError("Bad incoming variables input")

    if unknown == "Tc":
        Tc = unidades.Temperature(Tb.R*(3*log10(Pc.psi)/7/(w+1)+1), "R")
    elif unknown == "Pc":
        Pc = unidades.Pressure(10**(7/3.*(w+1)*(Tc.R/Tb.R-1)), "atm")
    elif unknown == "Tb":
        Tb = unidades.Temperature(Tb.R/(3*log10(Pc.atm)/7/(w+1)+1), "R")
    elif unknown == "w":
        w = unidades.Dimensionless(3/7*log10(Pc.atm)/(Tc.R/Tb.R-1)-1)

    prop = {}
    prop["Tc"] = Tc
    prop["Pc"] = Pc
    prop["Tb"] = Tb
    prop["w"] = w
    return prop


def facent_AmbroseWalton(Pvr):
    """Calculates acentric factor of a fluid using the Ambrose-Walton
    corresponding-states correlation

    Parameters
    ----------
    Pvr : float
        Reduced vapor pressure of compound at 0.7Tc, [-]

    Returns
    -------
    w : float
        Acentric factor [-]

    References
    ----------
    .. [8] Ambrose, D., Walton, J. Vapour Pressures up to Their Critical
        Temperatures of Normal Alkanes and 1-Alkanols. Pure & Appl. Chem. 61(8)
        1395-1403 (1989)
    .. [2] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    """
    Tr = 0.7
    t = 1-Tr
    f0 = (-5.97616*t + 1.29874*t**1.5 - 0.60394*t**2.5 - 1.06841*t**5)/Tr
    f1 = (-5.03365*t + 1.11505*t**1.5 - 5.41217*t**2.5 - 7.46628*t**5)/Tr
    f2 = (-0.64771*t + 2.41539*t**1.5 - 4.26979*t**2.5 + 3.25259*t**5)/Tr
    coef = roots([f2, f1, f0-log(Pvr)])

    if absolute(coef[0]) < absolute(coef[1]):
        return coef[0]
    else:
        return coef[1]


def Vc_Riedel(Tc, Pc, w, M):
    """Calculates critical volume of a fluid using the Riedel correlation
    as explain in [5]_, procedure 4A3.1, Pag. 302

    Parameters
    ----------
    Tc : float
        Critical temperature, [K]
    Pc : float
        Critical pressure, [Pa]
    w : float
        Acentric factor, [-]
    M : float
        Molecular weight, [g/mol]

    Returns
    -------
    Vc : float
        Critical volume, [Pa]

    Examples
    --------
    Example in [5]_, n-nonane
    >>> Tc = unidades.Temperature(610.68, "F")
    >>> Pc = unidades.Pressure(331.8, "psi")
    >>> "%0.3f" % Vc_Riedel(Tc, Pc, 0.4368, 128.2551).ft3lb
    '0.068'

    References
    ----------
    .. [17]
    .. [18] Riedel, L. Die Zustandsfunktion des realen Gases: Untersuchungen
        über eine Erweiterung des Theorems der übereinstimmenden Zustände.
        Chem. Ings-Tech. 28 (1956) 557-562
    .. [5] API. Technical Data book: Petroleum Refining 6th Edition
    """
    # Eq 1 in [18]
    alfa = 5.811 + 4.919*w

    Vc = R*1000*Tc/Pc/(3.72+0.26*(alfa-7))/M
    return unidades.SpecificVolume(Vc, "lg")


def ThL_RiaziFaghri(T, Tb, SG):
    """Calculates thermal conductivity of liquid hydrocarbon at low pressure
    using the Riazi-Faghri correlation.

    .. math::
        \kappa = aT_{b}^{b}SG^{c}

        a = \exp\left(-4.5093-0.6844t-0.1305t^{2}\right)

        b = 0.3003+0.0918t+0.0195t^{2}

        c = 0.0129+0.0894t+0.0292t^{2}

    where t = T(F)/100

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tb : float
        Normal boiling temperature, [K]
    SG : float
        Specific gravity, [-]

    Returns
    -------
    k : float
        Thermal conductivity, [Btu/hftºF]

    Notes
    -----
    Range of validity:
        0ºF ≤ T ≤ 300ºF

    References
    ----------
    .. [11] Riazi, M.R., Faghri, A. Thermal Conductivity of Liquid and Vapor
        Hydrocarbon Systems: Pentanes and Heavier at Low Pressures. Ind. Eng.
        Chem. Process Des. Dev. 24 (1985) 398-401
    """
    # Convert input Tb in Kelvin to Fahrenheit to use in the correlation
    Tb_R = unidades.K2R(Tb)
    t = unidades.K2F(T)/100
    A = exp(-4.5093-0.6844*t-0.1305*t**2)                               # Eq 9
    B = 0.3003+0.0918*t+0.0195*t**2                                    # Eq 10
    C = 0.1029+0.0894*t+0.0292*t**2                                    # Eq 11
    k = 1.7307*A*Tb_R**B*SG**C                                          # Eq 7
    return unidades.ThermalConductivity(k, "BtuhftF")


def ThL_Gharagheizi(T, Pc, Tb, M, w):
    """Calculates the thermal conductivity of liquid using the Gharagheizi
    correlation.

    .. math::
        \kappa = 10^{-4}\left(10\omega+2P_c-2T+4+1.908\left(T_b+\frac{1.009B^2}
        {M^2}\right)+\frac{3.9287M^4}{B^4}+\frac{A}{B^8}\right)

        A = 3.8588M^8\left(1.0045B+6.5152M-8.9756\right)

        B = 16.0407M+2T_b-27.9074

    Parameters
    ----------
    T : float
        Temperature, [K]
    Pc : float
        Critical pressure, [Pa]
    Tb : float
        Normal boiling temperature, [K]
    M : float
        Molecular weight, [g/mol]
    w : float
        Acentric factor, [-]

    Returns
    -------
    k : float
        Thermal conductivity [W/m·k]

    References
    ----------
    .. [12] Gharagheizi, F., Ilani-Kashkouli, P., Sattari, M., Mohammadi, A.H.,
        Ramjugernath, D., Richon, D. Development of a General Model for
        Determination of Thermal Conductivity of Liquid Chemical Compounds at
        Atmospheric Pressure. AIChE Journal 59 (2013) 1702-1708
    """

    Pc_bar = Pc*1e-5
    B = 16.0407*M+2*Tb-27.9074                                          # Eq 6
    A = 3.8588*M**8*(1.0045*B + 6.5152*M - 8.9756)                      # Eq 5
    k = 1e-4*(10*w + 2*Pc_bar - 2*T + 4 + 1.908*(Tb+1.009*B**2/M**2) +
              3.9287*M**4/B**4 + A/B**8)                                # Eq 4
    return unidades.ThermalConductivity(k)


def ThL_LakshmiPrasad(T, M):
    """Calculates the thermal conductivity of liquid using the Lakshmi-Prasad
    correlation.

    .. math::
        \lambda = 0.0655-0.0005T + \frac{1.3855-0.00197T}{M^{0.5}}

    Parameters
    ----------
    T : float
        Temperature, [K]
    M : float
        Molecular weight, [g/mol]

    Returns
    -------
    k : float
        Thermal conductivity, [W/m/k]

    References
    ----------
    .. [13] Lakshmi, D.S., Prasad, D.H.L. A Rapid Estimation Method for Thermal
        Conductivity of Pure Liquids. The Chemical Engineering Journal 48
        (1992) 211-14
    """
    k = 0.0655 - 0.0005*T + (1.3855 - 0.00197*T)/M**0.5
    return unidades.ThermalConductivity(k)


def ThL_Nicola(T, M, Tc, Pc, w, mu=None):
    """Calculates the thermal conductivity of liquid using the Nicola
    correlation.

    .. math::
        \frac{\lambda}{\lambda_o} = aT_r + bPc + c\omega +
        \left(\frac{e}{M}\right)^{d}

        \frac{\lambda}{\lambda_o} = aT_r + bPc + c\omega +
        \left(\frac{e}{M}\right)^{d} + f\mu

    Parameters
    ----------
    T : float
        Temperature, [K]
    M : float
        Molecular weight, [g/mol]
    Tc : float
        Critical temperature, [K]
    Pc : float
        Critical pressure, [Pa]
    w : float
        Acentric factor, [-]
    mu : float
        Dipole moment, [Debye]

    Returns
    -------
    k : float
        Thermal conductivity [W/m·k]

    References
    ----------
    .. [14] Di Nicola, G., Ciarrocchi, E., Coccia, G., Pierantozzi, M.
        Correlations of Thermal Conductivity for Liquid Refrigerants at
        Atmospheric Pressure or near Saturation. International Journal of
        Refrigeration, 2014
    """
    Pc_bar = unidades.Pressure(Pc).bar
    if mu:
        # Eq 4 using dipole moment of compound
        k = 0.6542*(-0.2034*T/Tc+0.0013*Pc_bar+0.1714*w+(1/M)**0.3539-0.007*mu)
    else:
        # Eq 3
        k = 0.5147*(-0.2537*T/Tc+0.0017*Pc_bar+0.1501*w+(1/M)**-0.2999)
    return unidades.ThermalConductivity(k)


def ThL_Pachaiyappan(T, Tc, M, rho, branched=True):
    """Calculates the thermal conductivity of liquid using the Pachaiyappan
    correlation as explain in [5]_, procedure 12A1.2, pag 1141

    .. math::
        \kappa=\frac{CM^{n}}{Vm}\frac{3+20\left(1-Tr\right){}^{2/3}}
        {3+20\left(1-\frac{527.67}{Tc}\right)^{2/3}}

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tc : float
        Critical temperature, [K]
    M : float
        Molecular weight, [g/mol]
    rho : float
        Density of commpound at 68ºF, []
    branched : boolean, optional
        Linear or branched compound, default True

    Returns
    -------
    k : float
        Thermal conductivity [W/m·k]

    Examples
    --------
    Example in [5]_, n-butylbenzene at 140ºF
    >>> T = unidades.Temperature(140, "F")
    >>> Tc = unidades.Temperature(729.32, "F")
    >>> rho = unidades.Density(53.76, "lbft3")
    >>> k = ThL_Pachaiyappan(T, Tc, 134.22, rho)
    >>> "%0.4f" % k.BtuhftF
    '0.0673'

    References
    ----------
    .. [15] Pachaiyappan, V., Ibrahim, S.H., Kuloor, N.R. Thermal
        Conductivities of Organic Liquids: A New Correlation. J. Chem. Eng.
        Data, 11 (1966) 73-76
    .. [5] API. Technical Data book: Petroleum Refining 6th Edition
    """

    if branched:
        n = 0.7717
        C = 4.079e-3
    else:
        n = 1.001
        C = 1.676e-3

    Tr = T/Tc
    Tc_R = unidades.K2R(Tc)
    rhom = unidades.Density(rho/M).lbft3
    Vm = 1/rhom
    k = C*M**n/Vm*(3+20*(1-Tr)**(2./3))/(3+20*(1-527.67/Tc_R)**(2./3))
    return unidades.ThermalConductivity(k, "BtuhftF")


def ThL_KanitkarThodos(T, P, Tc, Pc, Vc, M, rho):
    """Calculates the thermal conductivity of liquid using the Kanitkar-Thodos
    correlation as explain in [5]_, procedure 12A1.3, pag 1143

    .. math::
        \kappa\lambda = -1.884e-6P_r^2 + 1.442e-3P_r +
        \alfa\exp\left(\beta\rho_r\right)

        \alfa = \frac{7.137e-3}{\beta^{3.322}}

        \beta = 0.4 + \frac{0.986}{\exp{0.58\lambda}}

        \lambda = \frac{Tc^{1/6}M^{1/2}}{Pc}^{2/3}

    Parameters
    ----------
    T : float
        Temperature, [K]
    P : float
        Pressure, [Pa]
    Tc : float
        Critical temperature, [K]
    Pc : float
        Critical pressure, [Pa]
    Vc : float
        Critical density, [m³/kg]
    M : float
        Molecular weight, [g/mol]
    rho : float
        Density of commpound at P and T, [kg/m³]

    Returns
    -------
    k : float
        Thermal conductivity [W/m·k]

    Notes
    -----
    This method let calculate the thermal conductivity of liquid hydrocarbons
    at any pressure

    Examples
    --------
    Example in [5]_, n-heptane at 320ºF and 197.4 atm
    >>> T = unidades.Temperature(320, "F")
    >>> P = unidades.Pressure(197.4, "atm")
    >>> Tc = unidades.Temperature(512.69, "F")
    >>> Pc = unidades.Pressure(397.41, "psi")
    >>> Vc = unidades.SpecificVolume(0.0684, "ft3lb")
    >>> rho = unidades.Density(37.93, "lbft3")
    >>> k = ThL_KanitkarThodos(T, P, Tc, Pc, Vc, 100.2, rho)
    >>> "%0.5f" % k.BtuhftF
    '0.06957'

    References
    ----------
    .. [16] Kanitkar, D., Thodos, G. The Thermal Conductivity of Liquid
        Hydrocarbons. Can. J. Chem. Eng. 47 (1969) 427-430
    .. [5] API. Technical Data book: Petroleum Refining 6th Edition
    """

    Pc_atm = Pc/101325
    Tc_R = unidades.K2R(Tc)
    Pr = P/Pc
    rhor = rho*Vc

    l = Tc_R**(1/6)*M**0.5/Pc_atm**(2/3)
    b = 0.4 + 0.986/exp(0.58*l)                                         # Eq 8
    alfa = 7.137e-3/b**3.322                                            # Eq 7
    k = (-1.884e-6*Pr**2+1.442e-3*Pr+alfa*exp(b*rhor))/l                # Eq 6
    return unidades.ThermalConductivity(k, "BtuhftF")


def ThG_RiaziFaghri(T, Tb, SG):
    """Calculates thermal conductivity of gas hydrocarbon at low pressure
    using the Riazi-Faghri correlation.

    .. math::
        \kappa = aT_{b}^{b}SG^{c}

        a = \exp\left(-4.5093-0.6844t-0.1305t^{2}\right)

        b = 0.3003+0.0918t+0.0195t^{2}

        c = 0.0129+0.0894t+0.0292t^{2}

    where t = T(F)/100

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tb : float
        Normal boiling temperature, [K]
    SG : float
        Specific gravity, [-]

    Returns
    -------
    k : float
        Thermal conductivity, [Btu/hftºF]

    Notes
    -----
    Range of validity:
        150ºF ≤ T ≤ 550ºF
        0.65 ≤ SG ≤ 0.9

    References
    ----------
    .. [11] Riazi, M.R., Faghri, A. Thermal Conductivity of Liquid and Vapor
        Hydrocarbon Systems: Pentanes and Heavier at Low Pressures. Ind. Eng.
        Chem. Process Des. Dev. 24 (1985) 398-401
    """
    # Convert input Tb in Kelvin to Fahrenheit to use in the correlation
    Tb_R = unidades.K2R(Tb)
    t = unidades.K2F(T)/100
    A = exp(21.78-8.07986*t+1.12981*t**2-0.05309*t**3)                 # Eq 16
    B = -4.13948+1.29924*t-0.17813*t**2+0.00833*t**3                   # Eq 17
    C = 0.19876-0.0313*t-0.00567*t**2                                  # Eq 18
    k = A*Tb_R**B*SG**C                                                # Eq 7
    return unidades.ThermalConductivity(k, "BtuhftF")


def Rackett(w):
    """Calculate the rackett constant using the Yamada-Gunn generalized
    correlation

    Parameters
    ----------
    w : float
        Acentric factor, [-]

    Returns
    -------
    Zra : float
        Rackett compressibility factor, [-]

    References
    ----------
    .. [23] Yamada, T., Gunn. R. Saturated Liquid Molar Volumes: The Rackett
        Equation. Journal of Chemical Engineering Data 18(2) (1973): 234–236
    """
    Zra = 0.29056 - 0.08775*w                                           # Eq 1
    return Zra



class Componente(object):
    """Class to define a chemical compound from the database"""
    _bool = False

    METHODS_RhoL = ["DIPPR", "Rackett", "Cavett", "COSTALD"]
    METHODS_RhoLP = ["Thomson-Brobst-Hankinson", "API"]
    METHODS_Pv = ["DIPPR", "Wagner", "Antoine", "Ambrose-Walton", "Lee-Kesler",
                  "Riedel", "Sanjari", "Maxwell-Bonnel"]
    METHODS_facent = ["Lee-Kesler", "Edmister", "Ambrose-Walton"]

    def __init__(self, id=None):
        """id: index of compound in database"""
        if not id:
            return

        self._bool = True
        self.id = id
        self.Config = config.getMainWindowConfig()
        componente = sql.getElement(id)
        self.formula = componente[1]
        self.nombre = componente[2]
        self.M = componente[3]
        self.SG = componente[124]
        self.Tc = unidades.Temperature(componente[4])
        self.Pc = unidades.Pressure(componente[5], "atm")
        self.Tb = unidades.Temperature(componente[131])
        self.Tf = unidades.Temperature(componente[132])
        if componente[125] != 0:
            self.f_acent = componente[125]
        elif self.Pc and self.Tc and self.Tb:
            self.f_acent = self._f_acent()
        else:
            self.f_acent = 0
        if componente[6] != 0:
            self.Vc = unidades.SpecificVolume(componente[6])
        elif self.f_acent != 0 and self.Tc != 0 and self.Pc != 0:
            self.Vc = Vc_Riedel(self.Tc, self.Pc, self.f_acent, self.M)
        else:
            self.Vc = 0
        if self.Tc:
            self.Zc = self.Pc*self.Vc/R/self.Tc
        else:
            self.Zc = 0
        if componente[7] != 0:
            self.API = componente[7]
        elif componente[124] != 0:
            self.API = 141.5/componente[124]-131.5
        else:
            self.API = 0
        self.cp = componente[8:14]

        # Parametric parameters
        self.antoine = componente[14:17]
        self.antoine += componente[152:156]
        self.wagner = componente[156:160]
        self._parametricMu = componente[21:23]
        self._parametricSigma = componente[23:25]
        self.henry = componente[17:21]

        # DIPPR parameters
        self._dipprRhoS = componente[25:33]
        self._dipprRhoL = componente[33:41]
        self._dipprPv = componente[41:49]
        self._dipprHv = componente[49:57]
        self._dipprCpS = componente[57:65]
        self._dipprCpL = componente[65:73]
        self._dipprCpG = componente[73:81]
        self._dipprMuL = componente[81:89]
        self._dipprMuG = componente[89:97]
        self._dipprKL = componente[97:105]
        self._dipprKG = componente[105:113]
        self._dipprSigma = componente[113:121]

        self.momento_dipolar = unidades.DipoleMoment(componente[121])
        if componente[123] != 0.0:
            self.rackett = componente[123]
        else:
            self.rackett = Rackett(self.f_acent)
        if componente[122] != 0.0:
            self.Vliq = componente[122]
#        elif self.Pc!=0 and self.Tc>298.15 :
#            self.Vliq=self.Volumen_Liquido_Constante()
        else:
            self.Vliq = 0

        if componente[126] != 0.0:
            self.parametro_solubilidad = unidades.SolubilityParameter(componente[126])
        else:
            self.parametro_solubilidad = self.Parametro_Solubilidad()
        self.Kw = componente[127]
        self.MSRK = componente[128:130]
        if componente[130] != 0.0:
            self.stiehl = componente[130]
        else:
            self.stiehl = 0
        # FIXME: No esta bien
#            self.stiehl=self.Stiehl_Polar_factor()
        self.CASNumber = componente[133]
        self.formula_alternativa = componente[134]
        self.UNIFAC = eval(componente[135])
        self.diametro_molecular = componente[136]
        self.ek = componente[137]
        self.UNIQUAC_area = componente[138]
        self.UNIQUAC_volumen = componente[139]
        if componente[140] == 0.0:
            self.f_acent_mod = componente[125]
        else:
            self.f_acent_mod = componente[140]
        self.calor_formacion = unidades.Enthalpy(componente[141]/self.M)
        self.energia_formacion = unidades.Enthalpy(componente[142]/self.M)
        self.wilson = componente[143]
        self.calor_combustion_neto = unidades.Enthalpy(componente[144]/self.M)
        self.calor_combustion_bruto = unidades.Enthalpy(componente[145]/self.M)
        self.nombre_alternativo = componente[146]
        self.V_char = componente[147]
        self.calor_formacion_solido = componente[148]
        self.energia_formacion_solido = componente[149]
        self.parametro_polar = componente[150]
        self.smile = componente[151]
        if self.smile != "" and os.environ["oasa"] == "True":
            import oasa
            mol = oasa.smiles.text_to_mol(self.smile, calc_coords=40)
            self.imageFile = tempfile.NamedTemporaryFile("w+r", suffix=".svg")
            oasa.svg_out.mol_to_svg(mol, self.imageFile.name)

        # TODO: Add branched property in database for each compound
        # For the moment define the branched property manually
        lineal = [1, 2, 3, 4, 6, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                  21, 90, 91, 92, 1763, 1765, 1767, 1738, 1769, 1770, 1844,
                  1741, 1842, 1771, 1742, 22, 23, 24, 25, 26, 28, 29, 30, 31,
                  35, 56, 57, 58]
        if self.id in lineal:
            self.branched = False
        else:
            self.branched = True

        # TODO: Añadir las contribuciones de grupos de parachor a la base de datos
        self.parachor = []
        # TODO: Añadir los tipos de cada elemento
        self.hidrocarburo = True
        self.Van_Veltzen = [] #Quizá se pueda derivar de otro grupos, UNIFAC o similar
        # TODO: Añadir caraceristicas químicas del componente
        self.isCiclico = False
        self.isHidrocarburo = True
        self.isLineal = True
        self.isAlcohol = False

        # TODO: Añadir parametros S1,S2 a la base de datos, API databook, pag 823
        self.SRKGraboski = [0, 0]

        # TODO: Añadir parámetros, archivo /media/datos/Biblioteca/archivos/Melhem, Almeida - A data Bank of Parameters for the Attractive-Aznar Telles.pdf
        self.Melhem = [0, 0]          #Alcoholes en archivo de abajo
        self.Almeida = [0, 0]

        # TODO: Añadir parámetros, archivo /media/datos/Biblioteca/archivos/alfas.pdf
        self.Mathias = 0
        self.MathiasCopeman = [0, 0, 0]
        self.Adachi = [0, 0]
        self.Andoulakis = [0, 0, 0]
        self.Yu_Lu = [0, 0, 0]

        # Desglosar formula en elementos y átomos de cada elemento
        decmp = atomic_decomposition(self.formula)
        self.C = decmp.get("C", 0)
        self.H = decmp.get("H", 0)
        self.O = decmp.get("O", 0)
        self.N = decmp.get("N", 0)
        self.S = decmp.get("S", 0)

        if self.C and self.H:
            self.HC = self.H/self.C

    def tr(self, T):
        return T/self.Tc

    def pr(self, P):
        return P/self.Pc

    def __bool__(self):
        return self._bool

    # Calculation of undefined properties of compound
    def _f_acent(self):
        """Acentric factor calculation in compounds with undefined property"""
        method = self.Config.getint("Transport", "f_acent")
        if method == 0:
            return facent_LeeKesler(self.Tb, self.Tc, self.Pc)
        elif method == 1:
            return prop_Edmister(Tb=self.Tb, Tc=self.Tc, Pc=self.Pc)
        elif method == 2:
            Pvr = self.Pv(0.7*self.Tc)/self.Pc
            return facent_AmbroseWalton(Pvr)

    def Volumen_Liquido_Constante(self):
        V=1/self.RhoL_Rackett(self.Tc)
        V=R_atml*1000*self.Tc/self.Pc.atm*self.rackett**(1+(1-self.tr(298.15))**(2.0/7)) #cm3/mol
        return V/(5.7+1611/self.Tc) #cm3/mol

    def Parametro_Solubilidad(self):
        """Método de cálculo del parametro de solubilidad, API databook pag 812"""
        if self._dipprHv==(0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0):
            Par=0
        else:
            DH=self.Hv_DIPPR(298.15).calg
            densidad = DIPPR("rhoL", T, self._dipprRhoL, M=self.M)
            Par=sqrt(DH-R_cal*298.15/self.M*densidad)
        return unidades.SolubilityParameter(Par, "calcc")

    def Stiehl_Polar_factor(self):
        f0=5.92714-6.09648/0.6-1.28862*log(0.6)+0.169347*0.6**6
        f1=15.2518-15.6875/0.6-13.4721*log(0.6)+0.43577*0.6**6
        pv=exp(f0*0.6+self.f_acent*f1*0.6)
        return log10(self.pr(P)/pv)

    # Ideal properties
    def _Cpo(self, T):
        """Ideal gas specific heat calculation procedure from polinomial
        coefficient in database in the form [A,B,C,D,E,F]
        Explained in procedure 7A1.1, pag 543

        .. math::
            Cp = A + BT + CT^2 + DT^3 + ET^4 + FT^5

        Parameters
        ----------
        T : float
            Temperature, [K]

        Notes
        -----
        The units in the calculate cp is in cal/mol·K

        References
        ----------
        .. [5] API. Technical Data book: Petroleum Refining 6th Edition
        """
        A, B, C, D, E, F = self.cp
        cp = A + B*T + C*T**2 + D*T**3 + E*T**4 + F*T**5
        return unidades.SpecificHeat(cp/self.M, "calgK")

    def _Ho(self, T):
        """Ideal gas enthalpy calculation from polinomial coefficient of
        specific heat saved in database
        Coefficient in database are in the form [A,B,C,D,E,F]
        Explained in procedure 7A1.1, pag 543

        .. math::
            Ho = BT + C/2T^2 + D/3T^3 + E/4T^4 + F/5T^5

        Parameters
        ----------
        T : float
            Temperature, [K]

        Notes
        -----
        The units in the calculate ideal enthalpy are in cal/mol·K, the
        reference state is set to T=298.15K

        References
        ----------
        .. [5] API. Technical Data book: Petroleum Refining 6th Edition
        """
        To = 298.15
        A, B, C, D, E, F = self.cp
        H = B*T + C/2*T**2 + D/3*T**3 + E/4*T**4 + F/5*T**5
        Ho = B*To + C/2*To**2 + D/3*To**3 + E/4*To**4 + F/5*To**5
        return unidades.Enthalpy((H-Ho)/self.M, "calg")

    def _so(self, T):
        """Ideal gas entropy calculation from polinomial coefficient of
        specific heat saved in database
        Coefficient in database are in the form [A,B,C,D,E,F]
        Explained in procedure 7A1.1, pag 543

        .. math::
            So = A\lnT + BT + C/2T^2 + D/3T^3 + E/4T^4 + F/5T^5

        Parameters
        ----------
        T : float
            Temperature, [K]

        Notes
        -----
        The units in the calculate ideal enthalpy are in cal/mol·K, the
        reference state is set to T=298.15K

        References
        ----------
        .. [5] API. Technical Data book: Petroleum Refining 6th Edition
        """
        A, B, C, D, E, F = self.cp
        so = A*log(T) + B*T + C/2*T**2 + D/3*T**3 + E/4*T**4 + F/5*T**5
        return unidades.SpecificHeat(so/self.M, "calgK")

    # Physical properties
    def RhoS(self, T):
        """Calculate the density of solid phase using the DIPPR equations"""
        return DIPPR("rhoS", T, self._dipprRhoS, M=self.M)

    def RhoL(self, T, P):
        """Calculate the density of liquid phase using any of available
        correlation"""
        rhoL = self.Config.getint("Transport", "RhoL")
        corr = self.Config.getint("Transport", "Corr_RhoL")
        if P < 1013250:
            if rhoL == 0 and self._dipprRhoL and \
                    self._dipprRhoL[6] <= T <= self._dipprRhoL[7]:
                return DIPPR("rhoL", T, self._dipprRhoL, M=self.M)
            elif rhoL == 1 and self.rackett != 0 and T < self.Tc:
                return RhoL_Rackett(T, self.Tc, self.Pc, self.rackett, self.M)
            elif rhoL == 2 and self.Vliq != 0:
                return self.RhoL_Cavett(T)
            elif rhoL == 3:
                if self.f_acent_mod != 0:
                    w = self.f_acent_mod
                else:
                    w = self.f_acent
                return RhoL_Costald(T, self.Tc, w, self.Vc)
            else:
                if self._dipprRhoL and \
                        self._dipprRhoL[6] <= T <= self._dipprRhoL[7]:
                    return DIPPR("rhoL", T, self._dipprRhoL, M=self.M)
                elif self.rackett != 0 and T < self.Tc:
                    return RhoL_Rackett(
                        T, self.Tc, self.Pc, self.rackett, self.M)
                elif self.Vliq != 0:
                    return self.RhoL_Cavett(T)
                else:
                    if self.f_acent_mod != 0:
                        w = self.f_acent_mod
                    else:
                        w = self.f_acent
                    return RhoL_Costald(T, self.Tc, w, self.Vc)
        else:
            if corr == 0:
                if self.f_acent_mod:
                    w = self.f_acent_mod
                else:
                    w = self.f_acent
                rhos = self.RhoL(T, 101325)
                Ps = self.Pv(T)
                return RhoL_ThomsonBrobstHankinson(
                    T, P, self.Tc, self.Pc, w, Ps, rhos)
            else:
                return self.RhoL_API(T, P)

    def RhoL_Cavett(self, T):
        """Método alternativo para calcular la densidad de liquidos haciendo uso
        de la ecuación de Cavett:       V = Vol_Con * (5.7 + 3Tr)
        donde   V: volumen de liquido en cc/mol
                Vol_Con es una constante para cada compuesto, situada en la base de datos en el puesto veinticinco
                Tr es la temperatura reducida
                Densidad obtenida en g/l"""
        return unidades.Density(1/(self.Vliq*(5.7+3*self.tr(T)))*1000*self.M)

    def RhoL_API(self, T, P):
        """Método alternativo para calcular la densidad de líquidos haciendo uso del método API
        Dens2=dens1*C2/C1
        donde:  C2,C1: valor empirico que se calcula a partir de una tabla de valores, C2 a T y P dadas, C1 a 60F y 1 atm
                dens1: densidad a 60ºF y 1 atm, situada en la base de datos en el puesto veintisiete
                densidad obtenida en mol/l"""
        pr=P/self.Pc
        A02=1.6368-0.04615*pr+2.1138e-3*pr**2-0.7845e-5*pr**3-0.6923e-6*pr**4
        A12=-1.9693-0.21874*pr-8.0028e-3*pr**2-8.2328e-5*pr**3+5.2604e-6*pr**4
        A22=2.4638-0.36461*pr-12.8763e-3*pr**2+14.8059e-5*pr**3-8.6895e-6*pr**4
        A32=-1.5841-0.25136*pr-11.3805e-3*pr**2+9.5672e-5*pr**3+2.1812e-6*pr**4
        C2=A02+A12*self.tr(T)+A22*self.tr(T)**2+A32*self.tr(T)**3
        A01=1.6368-0.04615*self.pr(1)+2.1138e-3*self.pr(1)**2-0.7845e-5*self.pr(1)**3-0.6923e-6*self.pr(1)**4
        A11=-1.9693-0.21874*self.pr(1)-8.0028e-3*self.pr(1)**2-8.2328e-5*self.pr(1)**3+5.2604e-6*self.pr(1)**4
        A21=2.4638-0.36461*self.pr(1)-12.8763e-3*self.pr(1)**2+14.8059e-5*self.pr(1)**3-8.6895e-6*self.pr(1)**4
        A31=-1.5841-0.25136*self.pr(1)-11.3805e-3*self.pr(1)**2+9.5672e-5*self.pr(1)**3+2.1812e-6*self.pr(1)**4
        t1=unidades.Temperature(60, "F")
        C1=A01+A11*self.tr(t1)+A21*self.tr(t1)**2+A31*self.tr(t1)**3
        d2=self.SG*1000*C2/C1
        #FIXME: C2 no sale bien, sale muchas veces negativo, repaso, repaso, pero no se porqué

        return unidades.Density(d2, "gl")

    def Pv(self, T):
        """Vapor pressure calculation procedure using the method defined in
        preferences"""
        method = self.Config.getint("Transport", "Pv")
        if method == 0 and self._dipprPv and \
                self._dipprPv[6] <= T <= self._dipprPv[7]:
            return DIPPR("Pv", T, self._dipprPv)
        elif method == 1 and self.wagner:
            return Pv_Wagner(T, self.Tc, self.Pc, self.wagner)
        elif method == 2 and self.antoine:
            return Pv_Antoine(T, self.antoine)
        elif method == 3 and self.Pc and self.Tc and self.f_acent:
            return Pv_AmbroseWalton(T, self.Tc, self.Pc, self.f_acent)
        elif method == 4 and self.Pc and self.Tc and self.f_acent:
            return Pv_Lee_Kesler(T, self.Tc, self.Pc, self.f_acent)
        elif method == 5 and self.Pc and self.Tc and self.Tb:
            return Pv_Riedel(T, self.Tc, self.Pc, self.Tb)
        elif method == 6 and self.Pc and self.Tc and self.Tb:
            return Pv_Sanjari(T, self.Tc, self.Pc, self.f_acent)
        elif method == 7 and self.Kw and self.Tb:
            return Pv_MaxwellBonnel(T, self.Tb, self.Kw)
        else:
            if self._dipprPv and self._dipprPv[6] <= T <= self._dipprPv[7]:
                return DIPPR("Pv", T, self._dipprPv)
            elif self.wagner:
                return Pv_Wagner(T, self.Tc, self.Pc, self.wagner)
            elif self.antoine:
                return Pv_Antoine(T, self.antoine)
            elif self.Pc and self.Tc and self.f_acent:
                return Pv_AmbroseWalton(T, self.Tc, self.Pc, self.f_acent)
            elif self.Pc and self.Tc and self.Tb:
                return Pv_Riedel(T, self.Tc, self.Pc, self.Tb)
            elif self.Kw and self.Tb:
                return Pv_MaxwellBonnel(T, self.Tb, self.Kw)

    def ThCond_Liquido(self, T, P):
        """Liquid thermal conductivity procedure using the method defined in
        preferences, use the decision diagram in [5]_ Figure 12-0.2 pag 1135"""
        ThCondL = self.Config.getint("Transport", "ThCondL")
        corr = self.Config.getint("Transport", "Corr_ThCondL")
        p = unidades.Pressure(P, "atm")
        if p.psi < 500:
            if ThCondL == 0 and self._dipprKL and \
                    self._dipprKL[6] <= T <= self._dipprKL[7]:
                return DIPPR("kL", T, self._dipprKL)
            elif ThCondL == 1 and T < self.Tc:
                return self.ThCond_Liquido_Pachaiyappan(T)
            else:
                # print "Warning: Thermal conductivity of %s out of range" % self.nombre
                return DIPPR("kL", self.i_dipprKL[7], self._dipprKL)
        else:
            if corr == 0:
                return self.ThCond_Liquido_Lenoir(T, P)
            else:
                return self.ThCond_Liquido_Kanitkar_Thodos(T, P)

    def ThCond_Liquido_Lenoir(self, T, P, ko=[]):
        """Método alternativo para el cálculo de la conductividad de líquidos a alta presión, API procedure 12A4.1, pag 1156
        Opcionalmente puede aceptar otro parametros ko que indica un valor experimental de la conductividad térmica (WmK, así como la temperatura y presión a la que se da, en un array [k,T,P]"""
        Tr=self.tr(T)
        if ko==[]:
            k1=self.ThCond_Liquido(T, 1)
            C1=17.77+0.065*self.pr(1)-7.764*Tr-2.065*Tr**2/exp(0.2*self.pr(1))
        else:
            k1=ko[0]
            C1=17.77+0.065*self.pr(ko[2])-7.764*self.tr(ko[1])-2.065*self.tr(ko[1])**2/exp(0.2*self.pr(ko[2]))
        C2=17.77+0.065*self.pr(P)-7.764*Tr-2.065*Tr**2/exp(0.2*self.pr(P))
        k=k1*C2/C1
        return unidades.ThermalConductivity(k)


    def ThCond_Gas(self, T, P):
        """Procedimiento que define el método más apropiado para el cálculo de la conductividad térmica del líquido, pag 1136"""
        ThCondG=self.Config.getint("Transport","ThCondG")
        p=unidades.Pressure(P)
        if p.psi<50:
            if ThCondG==0 and self._dipprKG and self._dipprKG[6]<=T<=self._dipprKG[7]:
                return self.DIPPR("kG", T, self._dipprKG)
            else:
                return self.ThCond_Gas_Misic_Thodos(T)
        elif self.id in [1, 46, 47, 48, 50, 51, 111]:
            return self.ThCond_Gas_Nonhidrocarbon(T, P)
        else:
            # TODO: fix crooks methods with lost lee_kesler library
            return self.ThCond_Gas(T, 101325.)
#            return self.ThCond_Gas_Crooks(T, P)

    def ThCond_Gas_Misic_Thodos(self, T):
        """Método alternativo para el cálculo de la conductividad térmica de gases a baja presión <4 atm, API Procedure 12B1.2 pag.1162"""
        l=(self.Tc.R)**(1./6)*self.M**0.5*(1/self.Pc.atm)**(2./3)
        cp=self.Cp_Gas_DIPPR(T)
        #TODO: Cuando se añada alguna propiedad en la base de datos que defina la naturaleza cíclica de los componentes se podrá generalizar este metodo a una temperatura reducida menor de 1 para compuestos cíclicos e hidrógeno, de momento todos se calculan por el metodo de tr mayor de 1.
#        if self.tr(T)<1:
#            k=1.188e-3*self.tr(T)*cp.BtulbF/l
#        else:
        k=2.67e-4*(14.52*self.tr(T)-5.14)**(2.0/3)*cp.BtulbF*self.M/l
        return unidades.ThermalConductivity(k, "BtuhftF")

    def ThCond_Gas_Crooks(self, T, P):
        """Método alternativo para el cálculo de la conductividad térmica de gases a alta presión >4 atm, API Procedure 12B4.1 pag.1170"""
        Tr=self.tr(T)
        Pr=self.pr(P)
        k=self.ThCond_Gas(T, 1)

        A=-0.0617*exp(1.91/Tr**9)
        B=2.29*exp(1.34/Tr**16)
        k_prima=1+(4.18/Tr**4+0.537*Pr/Tr**2)*(1-exp(A*Pr**B))+0.510*Pr/Tr**3*exp(A*Pr**B)
        # TODO: particularizar a compuestos cíclicos y no cíclicos
        k_prima2=1+1./Tr**5*Pr**4/(2.44*Tr**20+Pr**4)+0.012*Pr/Tr
        Cv=self.Cv_Lee_Kesler(T, P, 1)
        Cv0, Cvh, vr0, vrh=eos.Lee_Kesler_lib_Cp(Tr, Pr)
        Cv_prima=4.965-R_Btu*(1+Cv0)
        Cv_prima2=Cv.BtulbF-Cv_prima
        return unidades.ThermalConductivity(k*(Cv_prima/Cv.BtulbF*k_prima+Cv_prima2/Cv.BtulbF*k_prima2))

    def ThCond_Gas_Nonhidrocarbon(self, T, P):
        """Método de cálculo de la conductividad térmica de gases no hidrocarburos a alta presión, API procedure 12C1.1, pag 1174
        """
        t=unidades.Temperature(T)
        p=unidades.Pressure(P, "atm")
        if self.id==1:
            par=[4.681e-3, 2.e-4, -3.6e-8, 0.0, 0.0, 0.0, 1.7e-3]
        elif self.id==46:
            par=[4.561e-3, 1.61e-5, 0.0, 2.56e-9, 5.299e-3, 2.47e-3, 0.0]
        elif self.id==47:
            par=[5.95e-4, 1.71e-5, 0.0, -2.10e-8, 5.869e-3, 6.995e-3, 0.0]
        elif self.id==48:
            par=[1.757e-3, 1.55e-5, 0.0, 2.08e-8, 5.751e-3, 5.6e-3, 0.0]
        elif self.id==50:
            par=[-1.51e-3, 2.25e-5, 3.32e-10, 0.0, 0.0, 0.0, 0.0]
        elif self.id==51:
            par=[2.5826e-2, 1.35e-5, 0.0, -4.4e-7, 1.026e-2, -2.631e-2, 0.0]
        elif self.id==111:
            par=[-1.02e-3, 1.35e-5, 4.17e-9, 0.0, 0.0, 0.0, 0.0]

        k=par[0]+par[1]*t.R+par[2]*t.R**2+par[3]*p.psi+par[4]*p.psi/t.R**1.2+par[5]/(0.4*p.psi-0.001*t.R)**0.015+par[6]*log(p.psi)
        return unidades.ThermalConductivity(k, "BtuhftF")


    def Mu_Gas(self, T, P):
        """Procedimiento que define el método más apropiado para el cálculo de la viscosidad del líquido, pag 1026"""
        MuG = self.Config.getint("Transport", "MuG")
        if P/self.Pc < 0.6:
            if MuG == 0 and self._dipprMuG and \
                    self._dipprMuG[6] <= T <= self._dipprMuG[7]:
                return DIPPR("muG", T, self._dipprMuG)
            elif MuG == 1:
                return self.Mu_Gas_Chapman_Enskog(T)
            else:
                return MuG_StielThodos(T, self.Tc, self.Pc, self.M)
        else:
            if self.hidrocarburo:
                return self.Mu_Gas_Eakin_Ellingtong(T, P)
            else:
                return self.Mu_Gas_Carr(T, P)

    def Mu_Gas_Chapman_Enskog(self,T):
        """Método alternativo para calcular la viscosidad de gases (a baja presión):
        μg = 5/16(πMRT)^0.5/πO²Ω = 22.69*(MT)^0.5/O²Ω = Mp*P
        where:        M  = peso molecular
                      R  = constante del gas ideal
                      T  = temperatura
                      P  = presión
                      O  = diametro molécular
                      Ωv = colision integral
                      Mp = momento dipolar
        ref chemcad pag 81
        ref Properties gases and liquids pag. 470
        """
        #Metodo de Chung (pouling pag 473)
        if not self.diametro_molecular:
            diametro_molecular=0.809*self.Vc*self.M**(1./3)
        else: diametro_molecular=self.diametro_molecular
        if not self.ek:
            ek=self.Tc/1.2593
        else: ek=self.ek

        T_=T/ek
        if self.parametro_polar: #Polar, colisión integral de Brokaw (pouling pag 640)
            omega=1.03036/T_**0.15610+0.193/exp(0.47635*T_)+1.03587/exp(1.52996*T_)+1.76474/exp(3.89411*T_)+0.19*self.parametro_polar**2/T_
        else: #No polar, colisión integral de Neufeld
            omega=1.16145/T_**0.14874+0.52487/exp(0.7732*T_)+2.16178/exp(2.43787*T_)
        return unidades.Viscosity(26.69*(self.M*T)**0.5/diametro_molecular**2/omega, "microP")

    def Mu_Gas_Jossi(self, T, P, muo=0):
        """Método de cálculo de la viscosidad de hidrocarburos gaseosos pesados a alta presión,
        Jossi, J. A., Stiel, L. I., and Thodos, G., "The Viscosity of Pure Substances in the Dense Gaseous and Liquid Phases," American Institute of Chemical Engineers Journal, Vol. 8, 1962, pp. 59-63."""
        if muo==0:
            muo=self.Mu_Gas(T, 1)
        x=self.Tc**(1.0/6)/self.M**0.5/self.Pc.atm**(2.0/3)
        rhor=self.RhoG_Lee_Kesler(T, P)*self.Vc*self.M
        mu=((0.1023+0.023367*rhor+0.058533*rhor**2-0.040758*rhor**3+0.0093324*rho**4)**4-1e-4)*x+muo
        return unidades.Viscosity(mu, "cP")

    def Mu_Gas_Eakin_Ellingtong(self, T, P, muo=0):
        """Método de cálculo de la viscosidad de hidrocarburos gaseosos a alta presión, API procedure 11B4.1, pag 1107"""
        if muo==0:
            muo=self.Mu_Gas(T, 101325)
        else:
            muo=unidades.Viscosity(muo)
        x=self.Tc**(1.0/6)/self.M**0.5/self.Pc.atm**(2.0/3)
        rhor=self.RhoG_Lee_Kesler(T, P)*self.Vc*self.M
        mu=muo.cP+10.8e-5*(exp(1.439*rhor)-exp(-1.11*rhor**1.858))/x
        return unidades.Viscosity(mu, "cP")

    def Mu_Gas_Carr(self, T, P, muo=0):
        """Método de cálculo de la viscosidad de no hidrocarburos gaseosos a alta presión, API procedure 11C1.2, pag 1113"""
        Tr=self.tr(T)
        Pr=self.pr(P)
        if muo==0:
            muo=self.Mu_Gas(T, 101325)

        A1=83.8970*Tr**0.0105+0.6030*Tr**-0.0822+0.9017*Tr**-0.12-85.3080
        A2=1.514*Tr**-11.3036+0.3018*Tr**-0.6856+2.0636*Tr**-2.7611
        k=A1*1.5071*Pr**-0.4487+A2*(11.4789*Pr**0.2606-12.6843*Pr**0.1773+1.6953*Pr**-0.1052)
        return unidades.Viscosity(muo*k)

    def Mu_Gas_Stiel_Thodos(self, T, P, muo):
        """Método de cálculo de la viscosidad de gases polares a alta presión,
        Stiel, L. I. and Thodos, G., "The Viscosity of Polar Substances in the Dense Gaseous and Liquid Regions," American Institute of Chemical Engineers Journal, Vol. 10, No. 2, 1964, pp. 275-277."""
        if muo==0:
            muo=self.Mu_Gas(T, 1)
        rhor=self.RhoG_Lee_Kesler(T, P)*self.Vc*self.M
        x=self.Tc**(1.0/6)/self.M**0.5/self.Pc.atm**(2.0/3)

        if rhor<=0.1:
            mu=1.626e-4*rhor**1.111/x+muo
        elif 0.1<rhor<=0.9:
            mu=6.07e-6*(9.045*rhor+0.63)**1.739/x+muo
        else:
            mu=10**(4-10**(0.6239-0.1005*rhor))/x/1e4+muo
        return unidades.Viscosity(mu, "cP")


    def Mu_Liquido(self, T, P):
        """Procedimiento que define el método más apropiado para el cálculo de la viscosidad del líquido, pag 1026"""
        #Comparacion de métodos: pag 405 Vismanath
        MuL = self.Config.getint("Transport", "MuL")
        corr = self.Config.getint("Transport", "Corr_MuL")

        if P < 1013250:
            if MuL == 0 and self._dipprMuL and \
                    self._dipprMuL[6] <= T <= self._dipprMuL[7]:
                return DIPPR("muL", T, self._dipprMuL)
            elif MuL == 1 and self._parametricMu:
                return MuL_Parametric(T, self._parametricMu)
            elif MuL == 2 and self.Tc and self.Pc and self.f_acent:
                return MuL_LetsouStiel(
                        T, self.M, self.Tc, self.Pc, self.f_acent)
            elif MuL==3 and self.Van_Veltzen:
                return self.Mu_Liquido_Van_Veltzen(T)
            else:
                if self._dipprMuL and \
                        self._dipprMuL[6] <= T <= self._dipprMuL[7]:
                    return DIPPR("muL", T, self._dipprMuL)
                elif self._parametricMu:
                    return MuL_Parametric(T, self._parametricMu)
                elif self.Tc and self.Pc and self.f_acent and self.M:
                    return MuL_LetsouStiel(
                            T, self.M, self.Tc, self.Pc, self.f_acent)
                elif self.Van_Veltzen:
                    return self.Mu_Liquido_Van_Veltzen(T)
        else:
            if corr==0 and self.Tb<650:
            #En realidad el criterio de corte es los hidrocarburos de menos de 20 átomos de carbono (hidrocarburos de bajo peso molecular), pero aprovechando que la temperatura de ebullición es proporcional al peso molecular podemos usar esta
                return self.Mu_Liquido_Graboski_Braun(T, P)
            elif corr==1:
                return self.Mu_Liquido_Kouzel(T, P)
            elif corr==2 and self.Pc and self.f_acent:
                return self.Mu_Liquido_Lucas(T, P)
            else:
                if self.Tb<650:
                    return self.Mu_Liquido_Graboski_Braun(T, P)
                elif self.Pc and self.f_acent:
                    return self.Mu_Liquido_Lucas(T, P)
                else:
                    return self.Mu_Liquido_Kouzel(T, P)

    def Mu_Liquido_Van_Veltzen(self, T):
        """Método alternativo para calcular la viscosidad de líquidos haciendo uso de la contribución de los grupos moleculares, API procedure 11A2.3, pag 1048"""
        #TODO: Método dificil de implementar debido a que se tiene que calcular la contribución de grupos

    def Mu_Liquido_Lucas(self, T, P, muo=0):
        """Método de cálculo de la viscosidad de líquidos a alta presión
        ref Viswanath, Viscosidad de líquidos, pag 145
        ref Poling, The Properties of Gases and Liquids, pag 521"""
        P=unidades.Pressure(P, "atm")
        Tr=self.tr(T)
        if muo==0:
            muo=self.Mu_Liquido(T, 1)
        else:
            muo=unidades.Viscosity(muo)
        Pvp=self.Pv(T)
        DPr=(P.atm-Pvp.atm)/self.Pc.atm
        A=0.9991-(4.674e-4/(1.0523*Tr**-0.03877-1.0513))
        C=-0.07921+2.1616*Tr-13.4040*Tr**2+44.1706*Tr**3-84.8291*Tr**4+96.1209*Tr**5-59.8127*Tr**6+15.6719*Tr**7
        D=0.3257/(1.0039-Tr**2.573)**0.2906-0.2086
        mu=muo*(1+D*(DPr/2.118)**A)/(1+C*self.f_acent*DPr)
        return unidades.Viscosity(mu)

    def Mu_Liquido_Graboski_Braun(self, T, P):
        """Método alternativo para el cálculo de la viscosidad en líquidos de bajo peso molécular a altas presiones, API procedure 11A5.1 (pag. 1074)"""
        Pr=self.pr(P)
        Tr=self.tr(T)

        A1=3.0294*Tr**9.0740+0.0032*Tr**10.9399-0.3689
        A2=-0.038*Tr**-7.2309+0.0229*Tr**11.7631+0.5781
        A3=-0.1415*Tr**27.2842+0.0778*Tr**-4.3406+0.0014
        A4=0.0028*Tr**69.4404-0.0042*Tr**3.3586+0.0062
        A5=0.0107*Tr**-7.4626-85.8276*Tr**0.1392+87.3164
        mur0=A1*log10(Pr)+A2*log10(Pr)**2+A3*Pr+A4*Pr**2+A5

        if Pr <= 0.75:
            B1=-0.2462*Tr**0.0484-0.7275*log(Tr)-0.0588*Tr+0.0079
            B2=-0.3199*Tr**17.0626-0.0695*log(Tr)+0.1267*Tr-0.0101
            B3=4.7217*Tr**-1.9831+19.2008*Tr**-1.7595+65.5728*log(Tr)+0.6110*Tr-19.1590
        else:
            B1=-0.0214*Tr**0.0484-0.1827*log(Tr)-0.0183*Tr+0.0090
            B2=-0.3588*Tr**5.0537-0.1321*log(Tr)+0.0204*Tr-0.0075
            B3=3.7266*Tr**-2.5689+52.1358*Tr**0.3514-13.0750*log(Tr)+0.6358*Tr-56.6687
        mur1=B1*Pr+B2*log(Pr)+B3

        mur=mur0+self.f_acent*mur1
        #TODO: Para implementar correctamente este método hay que añadir a la base de datos los datos de la viscosidad en el punto crítico. De momento usaremos el procedimiento DIPPR como alternativa para obtener un valor de viscosidad en las condiciones estandart y una minitabla para los elementos que aparecen en la tabla 11A5.4 del API Databook
        if 1<self.id<=21 and self.id!=7:
            muc=[0, 0, 0.014, 0.02, 0.0237, 0.027, 0.0245, 0, 0.0255, 0.0350, 0.0264, 0.0273, 0.0282, 0.0291, 0.0305, 0.0309, 0.0315, 0.0328, 0.0337, 0.0348, 0.0355, 0.0362][self.id]
        elif self.id==90:
            muc=0.0370
        elif self.id==91:
            muc=0.0375
        elif self.id==92:
            muc=0.0388
        else:
#            muc=self.Mu_critica().cP
            To=298.15
            Po=101325
            muo=self.Mu_Liquido(To, 101325)

            Pr=self.pr(Po)
            Tr=self.tr(To)

            A1=3.0294*Tr**9.0740+0.0032*Tr**10.9399-0.3689
            A2=-0.038*Tr**-7.2309+0.0229*Tr**11.7631+0.5781
            A3=-0.1415*Tr**27.2842+0.0778*Tr**-4.3406+0.0014
            A4=0.0028*Tr**69.4404-0.0042*Tr**3.3586+0.0062
            A5=0.0107*Tr**-7.4626-85.8276*Tr**0.1392+87.3164
            muor0=A1*log10(Pr)+A2*log10(Pr)**2+A3*Pr+A4*Pr**2+A5

            if Pr <= 0.75:
                B1=-0.2462*Tr**0.0484-0.7275*log(Tr)-0.0588*Tr+0.0079
                B2=-0.3199*Tr**17.0626-0.0695*log(Tr)+0.1267*Tr-0.0101
                B3=4.7217*Tr**-1.9831+19.2008*Tr**-1.7595+65.5728*log(Tr)+0.6110*Tr-19.1590
            else:
                B1=-0.0214*Tr**0.0484-0.1827*log(Tr)-0.0183*Tr+0.0090
                B2=-0.3588*Tr**5.0537-0.1321*log(Tr)+0.0204*Tr-0.0075
                B3=3.7266*Tr**-2.5689+52.1358*Tr**0.3514-13.0750*log(Tr)+0.6358*Tr-56.6687
            muor1=B1*Pr+B2*log(Pr)+B3
            muor=muor0+self.f_acent*muor1
            muc=muo.cP/muor

        return unidades.Viscosity(mur*muc, "cP")

    def Mu_Liquido_Kouzel(self, T, P, mua=0):
        """Método alternativo para el cálculo de la viscosidad en líquidos de alto peso molécular a altas presiones, API procedure 11A5.5 (pag. 1081)
        como parámetro opcional se puede indicar la viscosidad a presión atmosferica a tempratura T"""
        p=unidades.Pressure(P, "atm")
        if mua==0:
            mua=self.Mu_Liquido(T, 1)
        mup=mua*10**(p.psig/1000*(-0.0102+0.04042*mua.cP**0.181))
        return unidades.Viscosity(mup)

    def Mu_critica(self):
        """Procedimiento que define la viscosidad crítica si no viene en la base de datos
        Eq 8.9 Riazi - Characterization and Properties of Petroleum fractions, pag 348"""
        x=self.Tc**(1.0/6)/self.M**0.5/self.Pc.atm**(2.0/3)
        return unidades.Viscosity(7.7e-4/x, "cP")

    def Tension(self, T):
        """Procedimiento que define el método más apropiado para el cálculo de la tensión superficial"""
        tension = self.Config.getint("Transport", "Tension")
        if tension == 0 and self._dipprSigma and \
                self._dipprSigma[6] <= T <= self._dipprSigma[7]:
            return DIPPR("sigma", T, self._dipprSigma)
        elif tension == 1 and self._parametricSigma:
            return Tension_Parametric(T, self._parametricSigma, self.Tc)
        elif tension==2 and self.parachor:
            return self.Tension_Parachor(T)
        elif tension==3:
            return self.Tension_MIller(T)
        elif tension==4 and self.stiehl:
            return self.Tension_Hakim(T)
        elif tension==5 and self.Kw:
            return self.Tension_Hydrocarbon(T)
        else:
            if self._dipprSigma and \
                    self._dipprSigma[6] <= T <= self._dipprSigma[7]:
                return DIPPR("sigma", T, self._dipprSigma)
            elif self._parametricSigma:
                return Tension_Parametric(T, self._parametricSigma, self.Tc)
            elif self.parachor:
                return self.Tension_Parachor(T)
            elif self.stiehl:
                return self.Tension_Hakim(T)
            elif self.Kw:
                return self.Tension_Hydrocarbon(T)
            else:
                return self.Tension_MIller(T)

    def Tension_Hakim(self, T):
        """Método alternativo para el cálculo de la tensión superficial de líquidos.
        ref Chemcad pag 85
        ref Properties of gases and liquids pag 693 y sig."""
        Qp=0.1574+0.385*self.f_acent-1.769*self.stiehl-13.69*self.stiehl**2-0.510*self.f_acent**2+1.298*self.f_acent*self.stiehl
        m=1.21+0.5385*self.f_acent-14.61*self.stiehl-32.07*self.stiehl**2-1.656*self.f_acent**2+22.03*self.f_acent*self.stiehl
        return unidades.Tension(self.Pc.atm**0.67*self.Tc**0.33*Qp*((1-self.tr(T))/0.4)**m, "dyncm")

    def Tension_Block_Bird(self, T):
        """Método alternativo para el cálculo de la tensión superficial de líquidos.
        ref Eq.8.88 Riazi-Characterization of petroleum fractions pag 373"""
        Tbr=self.Tb/self.Tc
        Q=0.1196*(1+Tbr*log(self.Pc.atm)/(1-Tbr))-0.279
        return unidades.Tension(self.Pc.atm**0.67*self.Tc**0.33*Q*(1-self.tr(T))**(11./9), "dyncm")

    def Tension_Miqueu(self, T):
        """Método alternativo para el cálculo de la tensión superficial de líquidos.
        Miqueu, C., Broseta, D., Satherley, J., Mendiboure, B., Lachiase, J., and Graciaa, A., "An Extended Scaled Equation for the Temperature Dependence of the Surface Tension of Pure Compounds Inferred From an Analysis of Experimental Data," Fluid Phase Equilibria, Vol. 172, 2000, pp. 169-182."""
        t=1-self.tr(T)
        return unidades.Tension(Bolzmann*1e7*self.Tc*(Avogadro/self.Vc.ccg)**(2./3)*(4.35+4.14*self.f_acent)*(1+0.19*t**0.5-0.25*t)*t**1.26, "dyncm")

    def Tension_Hydrocarbon(self, T):
        """Método alternativo para el cálculo de la tensión superficial de líquidos"""
        t=unidades.Temperature(T)
        return unidades.Tension(673.7/self.Kw*((self.Tc.R-t.R)/self.Tc.R)**1.232, "dyncm")

    def Tension_MIller(self, T):
        """Método alternativo para el cálculo de la tensión superficial de líquidos"""
        Q=0.1207*(1+self.tr(self.Tb)*log(self.Pc.atm)/(1-self.tr(self.Tb)))-0.281
        return unidades.Tension(self.Pc.atm**0.67*self.Tc**0.33*Q*(1-self.tr(T))**(11.0/9), "dyncm")

    def Tension_Parachor(self, T, parachor):
        """Método alternativo para el cálculo de la tensión superficial de líquidos haciendo uso de las contribuciones de grupos
        API procedure 10A1.4, pag 987"""
        #TODO: Mientras no sepa como automatizar el cálculo de las contribuciones de grupo, habra que indicarlo como parámetro
        rhoL = DIPPR("rhoL", T, self._dipprRhoL, M=self.M)/1000
        rhoG=self.RhoG_Lee_Kesler(T, 1)*self.M/1000
        sigma=(parachor/self.M*(rhoL-rhoG))**4
        return unidades.Tension(sigma, "dyncm")





    def Hv_DIPPR(self,T):
        """Cálculo del calor de vaporización usando las ecuaciones DIPPR
        Los parámetros se encuentran en la base de datos en decimosexta posición
        Calor de vaporización obtenido en (J/kmol)"""
        return DIPPR("Hv", T, self._dipprHv, self.M)

    def Hv_Lee_Kesler(self, T):
        """Método alternativo para el cálculo del calor de vaporización haciendo uso de las propiedades críticas
        Procedure API 7C1.16 Pag.680
        Valor en J/mol"""
        Pv=self.Pv_DIPPR(T)
        Tr=T/self.Tc
        Pr=Pv/self.Pc
        H_adimensional_vapor=eos.Lee_Kesler_Entalpia_lib(Tr, Pr, self.f_acent, 1)
        H_adimensional_liquido=eos.Lee_Kesler_Entalpia_lib(Tr, Pr, self.f_acent, 0)
        return unidades.Enthalpy(R*self.Tc/self.M*(H_adimensional_vapor-H_adimensional_liquido), "Jg")

    def Cp_Solido_DIPPR(self,T):
        """Cálculo de la capacidad calorifica del solido usando las ecuaciones DIPPR
        Los parámetros se encuentran en la base de datos en decimoseptima posición
        Capacidad calorifica obtenida en (J/kmol·K)"""
        return DIPPR("cpS", T, self._dipprCpS, self.M)

    def Cp_Liquido_DIPPR(self,T):
        """Cálculo de la capacidad calorifica del liquido usando las ecuaciones DIPPR
        Los parámetros se encuentran en la base de datos en decimoctava posición
        Capacidad calorifica obtenida en (J/kmol·K)"""
        return DIPPR("cpL", T, self._dipprCpL, self.M)

    def Cp_Hadden(self, T):
        """Método alternativo para el cálculo de la capacidad calorífica en líquidos por debajo de su punto de ebullición
        Procedure API 7D1.9 Pag.696"""
        #TODO: No fácil de implementar al depender del elemento

    def Cp_Gas_DIPPR(self,T):
        """Cálculo de la capacidad calorifica del gas ideal usando las ecuaciones DIPPR
        Los parámetros se encuentran en la base de datos en decimonovena posición
        Capacidad calorifica obtenida en (J/kmol·K)"""
        return DIPPR("cpG", T, self._dipprCpG, self.M)

    def Cp_Lee_Kesler(self, T, P, fase=None):
        """Método alternativo para el cálculo de la capacidad calorífica
        Procedure API 7D3.6 Pag.711"""
        Tr=self.tr(T)
        if fase==None:
            fase=self.Fase(T, P)
        Cv0, Cvh, vr0, vrh=eos.Lee_Kesler_lib_Cp(T, P, fase)

        B=0.1181193-0.265728/Tr-0.154790/Tr**2-0.030323/Tr**3
        C=0.0236744-0.0186984/Tr
        D=0.155488e-4+0.623689e-4/Tr
        dpdt_0=1/vr0*(1+(0.1181193+0.154790/Tr**2+2*0.030323/Tr**3)/vr0+0.0236744/vr0**2+0.155488e-4/vr0**5-2*0.042724/Tr**3/vr0**2*((0.65392+0.060167/vr0**2)*exp(-0.060167/vr0**2)))
        dpdv_0=-Tr/vr0**2*(1+2*B/vr0+3*C/vr0**2+6*D/vr0**5+0.042724/Tr**3/vr0**2*(3*0.65392+(5-2*(0.65392+0.060167/vr0**2))*0.060167/vr0**2)*exp(-0.060167/vr0**2))
        Cp0=1+Tr*dpdt_0**2/dpdv_0+Cv0

        B=0.2026579-0.331511/Tr-0.027655/Tr**2-0.203488/Tr**3
        C=0.0313385-0.0503618/Tr+0.016901/Tr**3
        D=0.48736e-4+0.0740336e-4/Tr
        dpdt_h=1/vrh*(1+(0.2026579+0.027655/Tr**2+2*0.203488/Tr**3)/vrh+(0.0313385-2*0.016901/Tr**3)/vrh**2+0.48736e-4/vrh**5-2*0.041577/Tr**3/vrh**2*((1.226+0.03754/vrh**2)*exp(-0.03754/vrh**2)))
        dpdv_h=-Tr/vrh**2*(1+2*B/vrh+3*C/vrh**2+6*D/vrh**5+0.041577/Tr**3/vrh**2*(3*1.226+(5-2*(1.226+0.03754/vrh**2))*0.03754/vrh**2)*exp(-0.03754/vrh**2))
        Cph=1+Tr*dpdt_h**2/dpdv_h+Cvh

        Cp_adimensional=Cp0+self.f_acent/factor_acentrico_octano*(Cph-Cp0)
        return unidades.SpecificHeat(self._Cpo(T).JgK-R/self.M*Cp_adimensional, "JgK")

    def Cv_Lee_Kesler(self, T, P, fase=None):
        """Método de cálculo de la capacidad calorífica a volumen constante
        Procedure API 7E1.6 Pag.726"""
        #FIXME: No sale, un factor de 100 tengo que añadir no sé de donde
        Pr=P/self.Pc
        Tr=T/self.Tc
        if fase==None:
            fase=self.Fase(T, P)
        Cpo=self._Cpo(T)
        Cv0, Cvh, vr0, vrh=eos.Lee_Kesler_lib_Cp(Tr, Pr, fase)
        Cv_adimensional=Cv0+self.f_acent/factor_acentrico_octano*(Cvh-Cv0)
        return unidades.SpecificHeat(100*(Cpo.JgK-R/self.M*(1+Cv_adimensional)), "JgK")


    def Cp_Cv_Lee_Kesler(self, T, P):
        """Método de cálculo de la capacidad calorífica a volumen constante
        Procedure API 7E1.6 Pag.726"""
        Cv=self.Cv_Lee_Kesler(T, P.atm)
        Cp=self.Cp_Lee_Kesler(T, P.atm)
#        print Cp.BtulbF, Cv
        return Cp/Cv


    def Fugacidad_Lee_Kesler(self, T, P):
        """Método de cálculo de la fugacidad
        Procedure API 7G1.8 Pag.752"""
        Tr=T/self.Tc
        Pr=P/self.Pc
        f=eos.Lee_Kesler_Fugacidad_lib(Tr, Pr, self.f_acent, self.Fase(T, P))
        return unidades.Pressure(P*exp(f), "atm")

    def Entropia_Lee_Kesler(self, T, P):
        """Método de cálculo de la entropia
        Procedure API 7F1.7 Pag.739"""
        Tr=T/self.Tc
        Pr=P/self.Pc
        S0=self._so(T)
        H_adimensional=eos.Lee_Kesler_Entalpia_lib(Tr, Pr, self.f_acent, self.Fase(T, P.atm))
        f=eos.Lee_Kesler_Fugacidad_lib(Tr, Pr, self.f_acent, self.Fase(T, P.atm))
        S=H_adimensional+f+log(P/101325)

        return unidades.SpecificHeat(S0.JgK-R*S/self.M, "JgK")

    def constante_Henry(self,T,parameters=None):
        """constante H obtenida en psia por unidad de fracción molar del gas
        lnH = A/T + B∗lnT + C∗T + D
        Solo disponible para algunos compuestos:
        Hydrogen, Helium, Argon, Neon, Krypton, Xenon, Oxygen, Nitrogen,
        Hydrogen sulfide, Carbon monoxide, Carbon dioxide, Sulfur dioxide,
        Nitrous oxide, Chlorine,Bromine, Iodine, Methane, Ethane, Propane,
        Ethylene, Ammonia.
        API procedure 9A7.1, pag 927
        Los parametros de la ecuación se encuentran en la base de datos
        en forma de lista en la posición décima
        """
        if parameters==None:
            parameters=self.henry
        t=unidades.Temperature(T)
        return exp(parameters[0]/t.R+parameters[1]*log(t.R)+parameters[2]*t.R+parameters[3])


    def Fase(self, T, P):
        """Método que calcula el estado en el que se encuentra la sustancia"""
        Pv=self.Pv(T).atm
        if Pv>P:
            return 1
        else:
            return 0







    def RhoG_Lee_Kesler(self, T, P):
        a, b=eos.SRK_lib(self, T)
        Z_srk=eos.Z_Cubic_EoS(T, P, b, a, b, 0, b)
        Vvo=Z_srk[0]*R_atml*T/P

        vr0v, vrhv, vr0l, vrhl=eos.Lee_Kesler_lib(T/self.Tc, P/self.Pc.atm, fase=1, Vvo=Vvo)
        z0v=P/self.Pc.atm*vr0v/T*self.Tc
        zhv=P/self.Pc.atm*vrhv/T*self.Tc
        z=z0v+self.f_acent/factor_acentrico_octano*(zhv-z0v)
        return P/z/R_atml/T


if __name__ == '__main__':
    import doctest
    doctest.testmod()

#    etilbenceno=Componente(45)
#    t=unidades.Temperature(180, "F")
#    print "DIPPR: ", etilbenceno.Tension_DIPPR(t).dyncm
#    print "Paramétrica: ", etilbenceno.Tension_Parametrica(t).dyncm
#    print "Hakim: ", etilbenceno.Tension_Hakim(t).dyncm
#    print "Miller: ", etilbenceno.Tension_MIller(t).dyncm
#    print "Hydrocarbon: ", etilbenceno.Tension_Hydrocarbon(t).dyncm
#    print "Parachor: ", etilbenceno.Tension_Parachor(t, 285.1).dyncm
#    print "Miqueu: ", etilbenceno.Tension_Miqueu(t).dyncm
#    print "Block Bird: ", etilbenceno.Tension_Block_Bird(t).dyncm

#    ipentano=Componente(7)
#    t=unidades.Temperature(212, "F")
#    print "DIPPR: ", ipentano.ThCond_Gas_DIPPR(t).BtuhftF
#    print "Misic-Thodos: ", ipentano.ThCond_Gas_Misic_Thodos(t).BtuhftF


#    heptano=Componente(11)
#    t=unidades.Temperature(572, "F")
#    p=unidades.Pressure(1450, "psi")
#    print "Crooks: ", heptano.ThCond_Gas_Crooks(t, p.atm).BtuhftF

#    oxigeno=Componente(47)
#    t=unidades.Temperature(984.6, "R")
#    p=unidades.Pressure(6075, "psi")
#    print "Nonhidrocarbon: ", oxigeno.ThCond_Gas_Nonhidrocarbon(t, p.atm).BtuhftF

#    butilbenceno=Componente(78)
#    t=unidades.Temperature(140, "F")
#    print "DIPPR: ", butilbenceno.ThCond_Liquido_DIPPR(t).BtuhftF
#    print "Pachaiyappan: ", butilbenceno.ThCond_Liquido_Pachaiyappan(t).BtuhftF

#    heptano=Componente(11)
#    t=unidades.Temperature(320, "F")
#    print "Kanitkar Thodos: ", heptano.ThCond_Liquido_Kanitkar_Thodos(t, 197.4).BtuhftF
#    print "Lenoir: ", heptano.ThCond_Liquido_Lenoir(t, 197.4).BtuhftF

#    decano=Componente(14)
#    t=unidades.Temperature(104, "F")
#    print "DIPPR: ", decano.Mu_Liquido_DIPPR(t).cP
#    print "Paramétrico: ", decano.Mu_Liquido_Parametrica(t).cP
#    print "Letsou Steil: ", decano.Mu_Liquido_Letsou_Steil(t).cP
#
#    pentano=Componente(50)
#    t=unidades.Temperature(200, "F")
#    p=unidades.Pressure(3000, "psi")
#    print "Graboski Broun: ", pentano.Mu_Liquido_Graboski_Braun(t, p.atm).cP
#    print "Lucas: ", pentano.Mu_Liquido_Lucas(t, p.atm).cP

#    tetralin=Componente(376)
#    t=unidades.Temperature(302, "F")
#    print "DIPPR:  %0.5f" % tetralin.Pv_DIPPR(t).psi
#    print "Antoine:   %0.5f" % tetralin.Pv_Antoine(t).psi
#    print "Lee-Kesler:  %0.5f" % tetralin.Pv_Lee_Kesler(t).psi
#    print "Maxwell-Bonnel: %0.5f" % tetralin.Pv_Maxwell_Bonnel(t).psi
#    print "Wagner: %0.5f" % tetralin.Pv_Wagner(t).psi

#    propano=Componente(4)
#    t=unidades.Temperature(30, "F")
#    print "DIPPR: ", propano.RhoL_DIPPR(t).gml
#    print "Rackett: ", propano.RhoL_Rackett(t).gml
#    print "Cavett: ", propano.RhoL_Cavett(t).gml
#    print "Costald: ", propano.RhoL_Costald(t).gml

#    octano=Componente(12)
#    t=unidades.Temperature(212, "F")
#    p=unidades.Pressure(4410, "psi")
#    print "Thomson Brobst Hankinson: ", octano.RhoL_Thomson_Brobst_Hankinson(t, p.atm).kgl
#    print "API: ", octano.RhoL_API(t, p.atm).kgl


#    ciclohexano=Componente(38)      #ej pag 637
#    t=unidades.Temperature(300, "F")
#    p=unidades.Pressure(1000, "psi")
#    print ciclohexano.Z_SRK(t, p.atm)
#    print ciclohexano.Lee_Kesler_Entalpia(t, p.atm).Btulb
#    print ciclohexano.Entropia(t, p.atm).BtulbF*1.8

#    print ciclohexano.Hv_Lee_Kesler(422.04), ciclohexano.Calor_vaporizacion(422.04)
#    print ciclohexano.Cp_Lee_Kesler(422.04, 68.046), ciclohexano.Cv_Lee_kesler(422.04, 68.046)


#    isobutano=Componente(5)
#    t=unidades.Temperature(370, "F")
#    p=unidades.Pressure(4000, "psi")
#    print isobutano.Lee_Kesler_Fugacidad(t, p.atm).psi #Ej pag 745
#    t=unidades.Temperature(475, "F")
#    print isobutano.Lee_Kesler_Entropia(t, p.atm).BtulbF #Ej pag 733

#    print "     SRK    Lee_Kesler    BWRS"
#    print "Z  %5.4f   %7.4f   %5.4f" % (isobutano.Z_SRK(t, p.atm), isobutano.Z_Lee_Kesler(t, p.atm), isobutano.Z_BWRS(t, p.atm))
#    print isobutano.RhoG_Lee_Kesler(t, p.atm)
#    print isobutano.RhoG_SRK(t, p.atm)
#    print isobutano.RhoG_BWRS(t, p.atm)
#    print isobutano.Entalpia_SRK(t, p.atm)

#    buteno=Componente(24)
#    print buteno.f_acent
#    print buteno.factor_acentrico()


#    butano=Componente(6)
#    T=unidades.Temperature(200, "F")
#    print unidades.Enthalpy(butano.Entalpia_formacion(T)).Btulb



#    decano=Componente(14)
#    t=unidades.Temperature(104, "F")
#    print "DIPPR: ", decano.Mu_Liquido_DIPPR(t).cP
#    print "Paramétrico: ", decano.Mu_Liquido_Parametrica(t).cP
#    print "Letsou Steil: ", decano.Mu_Liquido_Letsou_Steil(t).cP

#    pentano=Componente(8)
#    t=unidades.Temperature(200, "F")
#    p=unidades.Pressure(3000, "psi")
#    print "Graboski Broun: ", pentano.Mu_Liquido_Graboski_Braun(t, p.atm).cP
#    print "Lucas: ", pentano.Mu_Liquido_Lucas(t, p.atm).cP

#    metilciclohexano=Componente(39)
#    t=unidades.Temperature(300, "K")
#    p=unidades.Pressure(500, "bar")
#    muo=unidades.Viscosity(0.68, "cP")
#    print "Graboski Broun: ", metilciclohexano.Mu_Liquido_Graboski_Braun(t, p.atm).cP
#    print "Lucas: ", metilciclohexano.Mu_Liquido_Lucas(t, p.atm).cP


#    decano=Componente(14)
#    print decano.MuL_Kouzel(unidades.Temperature(120, "F"), unidades.Pressure(9940, "psi").atm, unidades.Viscosity(52.7, "cP")).cP

#    propano=Componente(4)
#    t=unidades.Temperature(176, "F")
#    print propano.MuG_Thodos(t).cP

#    metano=Componente(2)
#    t=unidades.Temperature(543, "F")
#    print metano.Mu_Gas(t, 1).cP
#    print metano.Mu_Gas_Thodos(t).cP
#    print metano.Mu_Gas_Eakin_Ellingtong(t, 50).cP

#    nitrogeno=Componente(46)
#    t=unidades.Temperature(-58, "F")
#    p=unidades.Pressure(1677, "psi")
#    print nitrogeno.Mu_Gas(t, p).cP
#    print nitrogeno.Mu_Gas_Carr(t, p.atm).cP


#    from pylab import arange, plot, show
#    nonano=Componente(13)
#    t=linspace(0.3,1, 10)
#    C1=[]
#    C2=[]
#    C3=[]
#    C5=[]
#    C10=[]
#    C30=[]
#    for i in t:
#        C1.append(nonano.C_API(i, 1))
#        C2.append(nonano.C_API(i, 2))
#        C3.append(nonano.C_API(i, 3))
#        C5.append(nonano.C_API(i, 5))
#        C10.append(nonano.C_API(i, 10))
#        C30.append(nonano.C_API(i, 30))
#        #        C3.append(nonano.RhoL_API(i, nonano.Pc*3))
##        C5.append(nonano.RhoL_API(i, nonano.Pc*5))
##        C10.append(nonano.RhoL_API(i, nonano.Pc*10))
##        C30.append(nonano.RhoL_API(i, nonano.Pc*30))
#
#    plot(t, C1, t, C2)
#    show()


#    agua=Componente(62)
#    print agua.composicion_molecular
#    oxigeno=Componente(47)
#    print oxigeno.composicion_molecular
#    benceno=Componente(40)
#    print benceno.composicion_molecular
#    cfc=Componente(241)
#    print cfc.composicion_molecular

#    i_pentano=Componente(7)
#    t=unidades.Temperature(68, "F")
#    print i_pentano.Solubilidad_agua(t)
#    benceno=Componente(40)
#    t=unidades.Temperature(104, "F")
#    print benceno.Solubilidad_agua(t)

#    hexano=Componente(10)
#    t=unidades.Temperature(212, "F")
#    print hexano.Solubilidad_en_agua(t)

#    fluorene=Componente(197)
#    t=unidades.Temperature(122, "F")
#    print fluorene.Solubilidad_en_agua(t)

#    sulfuro=Componente(50)
#    t=unidades.Temperature(77, "F")
#    print sulfuro.Solubilidad_Henry(t, 1)

#    agua=Componente(62)
#    T=unidades.Temperature(100, "C")
#
#    print "Liquid Thermal Conductivity: ", agua.ThCond_Liquido(T, 1), "W/mK"
#    print "Liquid viscosity: ", agua.Mu_Liquido(T, 1), "Pa·s"
#    print "Liquid surface tension: ", agua.Tension(T), "N/m"
#    print "Gas Thermal Conductivity: ", agua.ThCond_Gas(T, 1), "W/mK"
#    print "Gas viscosity: ", agua.Mu_Gas(T, 1), "Pa·s"
#
#    print "Vapor pressure: ", agua.Pv(T).atm, "atm"

#    propeno=Componente(23)
#    t=unidades.Temperature(302, "F")
#    p=unidades.Pressure(2290, "psi")
#    print propeno.Cp_Lee_Kesler(t, p.atm).BtulbF
#    print propeno.Cp_Cv_ratio(t, p.atm)

#    SO2=Componente(51)
#    t=unidades.Temperature(300, "C")
#    print SO2.Mu_Gas_Chapman_Enskog(t, 1).microP

#    from pylab import arange, plot, show
#    nonano=Componente(13)
#    p=arange(0.2*nonano.Pc,5*nonano.Pc,1)
#    C1=[]
#    C2=[]
#    C3=[]
#    C5=[]
#    C10=[]
#    C30=[]#    for i in p:
#        C1.append(nonano.pr(i)*nonano.Lee_Kesler(nonano.Tc*1, i)[0]/1)
#        C11.append(nonano.Lee_Kesler(nonano.Tc*1.1, i))
#        C12.append(nonano.Lee_Kesler(nonano.Tc*1.2, i))
#        C13.append(nonano.Lee_Kesler(nonano.Tc*1.3, i))
#        C15.append(nonano.Lee_Kesler(nonano.Tc*1.5, i))
#        C17.append(nonano.Lee_Kesler(nonano.Tc*1.7, i))
#        C2.append(nonano.pr(i)*nonano.Lee_Kesler(nonano.Tc*2, i)[0]/2)
#        C25.append(nonano.Lee_Kesler(nonano.Tc*2.5, i))
#        C3.append(nonano.Lee_Kesler(nonano.Tc*3, i))
#        C4.append(nonano.Lee_Kesler(nonano.Tc*4, i))
#    plot(p/nonano.Pc, C1)
#    show()

#    Hidrogeno=Componente(1)
#    print unidades.Temperature(Hidrogeno.Tc).R
#    print unidades.Pressure(Hidrogeno.Pc, "atm").psi
#    print Hidrogeno.f_acent

#    agua=Componente(62)
#    print agua.SRK_Z(298.15, 1)
#    print agua.SRK_RhoG(298.15, 1).kgm3
#    print agua.SRK_Entalpia(298.15, 1).MJkg

#    agua=Componente(62)
#    t=400
#    print agua.BWRS_Z(298.15, 1)
#    print agua.van_Waals_Z(t, 1), agua.PR_Z(t, 1), agua.RK_Z(t, 1), agua.HPW_Z(t, 1, -0.5)
#    print agua.RK_Z(t, 1), agua.Wilson_Z(t, 1), agua.SRK_Z(t, 1)
#    print agua.BWRS_RhoG(298.15, 1).kgm3
#    print agua.BWRS_Entalpia(298.15, 1).MJkg

#    print agua.PR_V(298.15, 1)
#    print agua.PR_RhoG(298.15, 1)
#    print agua.PR_Entalpia(t, 1).MJkg, agua.Lee_Kesler_Entalpia(t, 1).MJkg, agua.iapws_Entalpia(t, 1).MJkg
#    print agua.Lee_Kesler_Z(t, 1), agua.SRK_Z(t, 1)

#    print agua.Cp_Gas_DIPPR(400), iapws_Cp(400, 1), agua.Cp_ideal(400)
#    print agua.Hv_Lee_Kesler(t).MJkg, agua.Hv_DIPPR(t).MJkg
#    print agua.Lee_Kesler_Entalpia(t, 1).MJkg, agua.iapws_Entalpia(t, 1).MJkg, agua.Entalpia_ideal(t).MJkg
#    print agua.Cp_Lee_Kesler(t, 1).JkgK*agua.M, agua.iapws_Cp(t, 1).JkgK*agua.M
#    print agua.Lee_Kesler_Entropia(t, 1).JkgK, agua.iapws_Entropia(t, 1).JkgK

#    agua=Componente(62)
#    from scipy import arange
#    from pylab import plot, grid, show
#    d=arange(270, 500, 10.)
#    y=[]
#    y2=[]
#    y3=[]
#    delta=[]
#    for i in d:
#        y.append(agua.Lee_Kesler_Entalpia(i, 1))
#        y2.append(agua.TB_Entalpia(i, 1))
#        y3.append(agua.iapws_Entalpia(i, 1))
##        delta.append(y3[-1]-y2[-1])
#    plot(d, y, d, y2, d, y3)
#    grid(True)
#    show()
#


#    sulfuro=Componente(50)
#    t=300
#    p=1
#    print sulfuro.H2S_V(t, p).ccg*sulfuro.M
#    print sulfuro.H2S_RhoG(t, p).gcc
#    print sulfuro.H2S_Z(t, p), sulfuro.TB_Z(t, p)
#    print sulfuro.H2S_Fugacidad(t, p)
#    print sulfuro.H2S_Entalpia(t, p).Jg*sulfuro.M

#    agua=Componente(62)
#    t=273
#    p=1
#    print agua.TB_Fugacidad(t, p), agua.Lee_Kesler_Fugacidad(t, p)
#    print agua.TB_U_exceso(t, p), agua.TB_H_exceso(t, p), agua.TB_S_exceso(t, p), agua.TB_Cv_exceso(t, p)
#    print agua.TB_Entalpia(t, p).MJkg, agua.Lee_Kesler_Entalpia(t, p).MJkg, agua.iapws_Entalpia(t, p).MJkg
#    print agua.TB_Joule_Thomson(t, p)

#    solido=Componente(533)
#    print solido.PT_lib(300)
#
#    Hexano=Componente(10)
#    print Hexano.Mu_Liquido(340, 1)




#    agua=Componente(62)
#    print [agua.Tension_Parametrica(t) for t in range(300, 350, 10)]
#    print agua.RhoL_Tait_Costald(300, 1)
#    print agua.Tc, agua.Pc.bar, agua.f_acent


