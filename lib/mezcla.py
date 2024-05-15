#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>.


This module implement mixture component properties

:class:`Mezcla`: The main class with all integrated functionality. Use the
properties in database and calculate state properties with the methods chosen
in configuration

Liquid density calculation methods:
    * :func:`RhoL_RackettMix`
    * :func:`RhoL_CostaldMix`
    * :func:`RhoL_AaltoKeskinenMix`
    * :func:`RhoL_TaitCostaldMix`
    * :func:`RhoL_NasrifarMix`
    * :func:`RhoL_APIMix`

Liquid viscosity calculation methods:
    * :func:`MuL_KendallMonroe`
    * :func:`MuL_Chemcad`

Gas viscosity calculation methods:
    * :func:`MuG_Reichenberg`
    * :func:`MuG_Lucas`
    * :func:`MuG_Chung`
    * :func:`MuG_Wilke`
    * :func:`MuG_Herning`
    * :func:`MuG_P_Chung`
    * :func:`MuG_TRAPP`
    * :func:`MuG_DeanStielMix`
    * :func:`MuG_APIMix`

Liquid thermal conductivity calculation methods:
    * :func:`ThL_Li`
    * :func:`ThL_Power`

Gas thermal conductivity calculation methods:
    * :func:`ThG_MasonSaxena`
    * :func:`ThG_LindsayBromley`
    * :func:`ThG_Chung`
    * :func:`ThG_StielThodosYorizane`
    * :func:`ThG_TRAPP`
    * :func:`ThG_P_Chung`

Surface tension calculation methods:
    * :func:`Tension`

Mixture mass definition:
    * :func:`mix_unitmassflow`
    * :func:`mix_unitmolarflow`
    * :func:`mix_massflow_massfraction`
    * :func:`mix_massflow_molarfraction`
    * :func:`mix_molarflow_massfraction`
    * :func:`mix_molarflow_molarfraction`
"""


from math import exp, pi

from numpy.lib.scimath import log10, log
from numpy.linalg import solve
from scipy.constants import R

from lib.compuestos import (Componente, RhoL_Costald, RhoL_AaltoKeskinen,
                            RhoL_TaitCostald, RhoL_Nasrifar, MuG_DeanStiel,
                            MuG_API, ThG_StielThodos)
from lib.physics import R_atml, Collision_Neufeld
from lib import unidades, config
from lib.utilities import refDoc


__doi__ = {
    1:
        {"autor": "Wilke, C.R.",
         "title": "A Viscosity Equation for Gas Mixtures",
         "ref": "J. Chem. Phys. 18(4) (1950) 517-519",
         "doi": "10.1063/1.1747673"},
    2:
        {"autor": "API",
         "title": "Technical Data book: Petroleum Refining 6th Edition",
         "ref": "",
         "doi": ""},
    3:
        {"autor": "Poling, B.E, Prausnitz, J.M, O'Connell, J.P",
         "title": "The Properties of Gases and Liquids 5th Edition",
         "ref": "McGraw-Hill, New York, 2001",
         "doi": ""},
    4:
        {"autor": "Li, C.C.",
         "title": "Thermal Conductivity of Liquid Mixtures",
         "ref": "AIChE Journal 22(5) (1976) 927-930",
         "doi": "10.1002/aic.690220520"},
    5:
        {"autor": "Lindsay, A.L., Bromley, L.A.",
         "title": "Thermal Conductivity of Gas Mixtures",
         "ref": "Ind. & Eng. Chem. 42(8) (1950) 1508-1511",
         "doi": "10.1021/ie50488a017"},
    6:
        {"autor": "Mason, E.A., Saxena, S.C.",
         "title": "Approximate Formula for the Thermal Conductivity of Gas "
                  "Mixtures",
         "ref": "Fhys. Fluids 1(5) (1958) 361-369",
         "doi": "10.1063/1.1724352"},
    7:
        {"autor": "Yorizane, M., Yoshiumra, S., Masuoka, H., Yoshida, H.",
         "title": "Thermal Conductivities of Binary Gas Mixtures at High "
                  "Pressures: N2-O2, N2-Ar, CO2-Ar, CO2-CH4",
         "ref": "Ind. Eng. Chem. Fundam. 22(4) (1983) 458-462",
         "doi": "10.1021/i100012a018"},
    8:
        {"autor": "Livingston, J.k Morgan, R., Griggs, M.A.",
         "title": "The Properties of Mixed Liquids III. The Law of Mixtures I",
         "ref": "J. Am. Chem. Soc. 39 (1917) 2261-2275",
         "doi": "10.1021/ja02256a002"},
    9:
        {"autor": "Spencer, C.F., Danner, R.P.",
         "title": "Prediction of Bubble-Point Density of Mixtures",
         "ref": "J. Chem. Eng. Data 18(2) (1973) 230-234",
         "doi": "10.1021/je60057a007"},
    10:
        {"autor": "Hankinson, R.W., Thomson, G.H.",
         "title": "A New Correlation for Saturated Densities of Liquids and "
                  "Their Mixtures",
         "ref": "AIChE Journal 25(4) (1979) 653-663",
         "doi": "10.1002/aic.690250412"},
    11:
        {"autor": "Aalto, M., Keskinen, K.I., Aittamaa, J., Liukkonen, S.",
         "title": "An Improved Correlation for Compressed Liquid Densities of "
                  "Hydrocarbons. Part 2. Mixtures",
         "ref": "Fluid Phase Equilibria 114 (1996) 21-35",
         "doi": "10.1016/0378-3812(95)02824-2"},
    12:
        {"autor": "Thomson, G.H., Brobst, K.R., Hankinson, R.W.",
         "title": "An Improved Correlation for Densities of Compressed Liquids"
                  " and Liquid Mixtures",
         "ref": "AIChE Journal 28(4) (1982): 671-76",
         "doi": "10.1002/aic.690280420"},
    13:
        {"autor": "Nasrifar, K., Ayatollahi, S., Moshfeghian, M.",
         "title": "A Compressed Liquid Density Correlation",
         "ref": "Fluid Phase Equilibria 168 (2000) 149-163",
         "doi": "10.1016/s0378-3812(99)00336-2"},
    14:
        {"autor": "Rea, H.E., Spencer, C.F., Danner, R.P.",
         "title": "Effect of Pressure and Temperature on the Liquid Densities "
                  "of Pure Hydrocarbons",
         "ref": "J. Chem. Eng. Data 18(2) (1973) 227-230",
         "doi": "10.1021/je60057a003"},
    15:
        {"autor": "Chung, T.H., Ajlan, M., Lee, L.L., Starling, K.E.",
         "title": "Generalized Multiparameter Correlation for Nonpolar and "
                  "Polar Fluid Transport Properties",
         "ref": "Ind. Eng. Chem. Res. 27(4) (1988) 671-679",
         "doi": "10.1021/ie00076a024"},
    16:
        {"autor": "Younglove, B.A., Ely, J.F.",
         "title": "Thermophysical Properties of Fluids. II. Methane, Ethane, "
                  "Propane, Isobutane, and Normal Butane",
         "ref": "J. Phys. Chem. Ref. Data 16(4) (1987) 577-798",
         "doi": "10.1063/1.555785"},
    17:
        {"autor": "Ely, J.F.",
         "title": "An Enskog Correction for Size and Mass Difference Effects "
                  "in Mixture Viscosity Prediction",
         "ref": "J. Res. Natl. Bur. Stand. 86(6) (1981) 597-604",
         "doi": "10.6028/jres.086.028"},
    18:
        {"autor": "Chueh, P.L., Prausnitz, J.M.",
         "title": "Vapor-Liquid Equilibria at High Pressures: Calculation of "
                  "Critical Temperatures, Volumes and Pressures of Nonpolar "
                  "Mixtures",
         "ref": "AIChE Journal 13(6) (1967) 1107-1113",
         "doi": "10.1002/aic.690130613"},
    19:
        {"autor": "Dean, D.E., Stiel, L.I.",
         "title": "The Viscosity of Nonpolar Gas Mixtures at Moderate and High"
                  " Pressures",
         "ref": "AIChE Journal 11(3) (1965) 526-532 ",
         "doi": "10.1002/aic.690110330"},
    20:
        {"autor": "Kendall, J., Monroe, P.",
         "title": "The Viscosity of Liquids II. The Viscosity-Composition "
                  "Curve for Ideal Liquid Mixtures",
         "ref": "J. Am. Chem. Soc. 39(9) (1917) 1787-1802",
         "doi": "10.1021/ja02254a001"},
    21:
        {"autor": "Ely, J.F., Hanley, H.J.M.",
         "title": "A Computer Program for the Prediction of Viscosity and "
                  "Thermal Condcutivity in Hydrocarbon Mixtures",
         "ref": "NBS Technical Note 1039 (1981)",
         "doi": ""},
}


def mix_unitmassflow(unitMassFlow, cmps):
    """Calculate mixture composition properties with known unitMassFlow"""
    massFlow = sum(unitMassFlow)
    unitMolarFlow = [mass/cmp.M for mass, cmp in zip(unitMassFlow, cmps)]
    molarFlow = sum(unitMolarFlow)
    molarFraction = [mi/molarFlow for mi in unitMolarFlow]
    massFraction = [mi/massFlow for mi in unitMassFlow]

    kw = {}
    kw["unitMassFlow"] = unitMassFlow
    kw["unitMolarFlow"] = unitMolarFlow
    kw["molarFlow"] = molarFlow
    kw["massFlow"] = massFlow
    kw["molarFraction"] = molarFraction
    kw["massFraction"] = massFraction
    return kw


def mix_unitmolarflow(unitMolarFlow, cmps):
    """Calculate mixture composition properties with known unitMolarFlow"""
    molarFlow = sum(unitMolarFlow)
    unitMassFlow = [mol*cmp.M for mol, cmp in zip(unitMolarFlow, cmps)]
    massFlow = sum(unitMassFlow)
    molarFraction = [mi/molarFlow for mi in unitMolarFlow]
    massFraction = [mi/massFlow for mi in unitMassFlow]

    kw = {}
    kw["unitMassFlow"] = unitMassFlow
    kw["unitMolarFlow"] = unitMolarFlow
    kw["molarFlow"] = molarFlow
    kw["massFlow"] = massFlow
    kw["molarFraction"] = molarFraction
    kw["massFraction"] = massFraction
    return kw


def mix_massflow_molarfraction(massFlow, molarFraction, cmps):
    """Calculate mixture composition properties with known massFlow and
    molarFraction"""
    pesos = [x*cmp.M for x, cmp in zip(molarFraction, cmps)]
    M = sum(pesos)
    molarFlow = massFlow/M
    massFraction = [peso/M for peso in pesos]
    unitMassFlow = [x*massFlow for x in massFraction]
    unitMolarFlow = [x*molarFlow for x in molarFraction]

    kw = {}
    kw["unitMassFlow"] = unitMassFlow
    kw["unitMolarFlow"] = unitMolarFlow
    kw["molarFlow"] = molarFlow
    kw["massFlow"] = massFlow
    kw["molarFraction"] = molarFraction
    kw["massFraction"] = massFraction
    return kw


def mix_massflow_massfraction(massFlow, massFraction, cmps):
    """Calculate mixture composition properties with known massFlow and
    massFraction"""
    unitMassFlow = [x*massFlow for x in massFraction]
    unitMolarFlow = [mass/cmp.M for mass, cmp in zip(unitMassFlow, cmps)]
    molarFlow = sum(unitMolarFlow)
    molarFraction = [mi/molarFlow for mi in unitMolarFlow]

    kw = {}
    kw["unitMassFlow"] = unitMassFlow
    kw["unitMolarFlow"] = unitMolarFlow
    kw["molarFlow"] = molarFlow
    kw["massFlow"] = massFlow
    kw["molarFraction"] = molarFraction
    kw["massFraction"] = massFraction
    return kw


def mix_molarflow_molarfraction(molarFlow, molarFraction, cmps):
    """Calculate mixture composition properties with known molarFlow and
    molarFraction"""
    unitMolarFlow = [x*molarFlow for x in molarFraction]
    unitMassFlow = [mol*cmp.M for mol, cmp in zip(unitMolarFlow, cmps)]
    massFlow = sum(unitMassFlow)
    massFraction = [mi/massFlow for mi in unitMassFlow]

    kw = {}
    kw["unitMassFlow"] = unitMassFlow
    kw["unitMolarFlow"] = unitMolarFlow
    kw["molarFlow"] = molarFlow
    kw["massFlow"] = massFlow
    kw["molarFraction"] = molarFraction
    kw["massFraction"] = massFraction
    return kw


def mix_molarflow_massfraction(molarFlow, massFraction, cmps):
    """Calculate mixture composition properties with known molarFlow and
    massFraction"""
    moles = [x/cmp.M for x, cmp in zip(massFraction, cmps)]
    M = sum(moles)
    massFlow = molarFlow*M
    molarFraction = [mol/molarFlow for mol in moles]
    unitMassFlow = [x*massFlow for x in massFraction]
    unitMolarFlow = [x*molarFlow for x in molarFraction]

    kw = {}
    kw["unitMassFlow"] = unitMassFlow
    kw["unitMolarFlow"] = unitMolarFlow
    kw["molarFlow"] = molarFlow
    kw["massFlow"] = massFlow
    kw["molarFraction"] = molarFraction
    kw["massFraction"] = massFraction
    return kw


@refDoc(__doi__, [18, 2])
def Vc_ChuehPrausnitz(xi, Vci, Mi, hydrocarbon=None):
    r"""Calculates critic volume of a mixture using the Chueh-Prausnitz
    correlation, also referenced in API procedure 4B3.1 pag 314

    .. math::
        V_{cm} = \sum_i^n\phi_iV_{ci}+\sum_i^n\sum_j^n\phi_i\phi_j\upsilon_{ij}

    .. math::
        \phi_j = \frac{x_jV_{cj}^{2/3}}{\sum_{i=1}^nx_iV_{ci}^{2/3}}

    .. math::
        \upsilon_{ij} = \frac{V_{ij}\left(V_{ci}+V_{cj}\right)}{2}

    .. math::
        V_{ij} = -1.4684\eta_{ij}+C

    .. math::
        \eta_{ij} = \left|\frac{V_{ci}-V_{cj}}{V_{ci}+V_{cj}}\right|

    Parameters
    ----------
    xi : list
        Mole fractions of components, [-]
    Vci : list
        Critical volume of components, [m³/kg]
    Mi : list
        Molecular weight of components, [g/mol]
    hydrocarbon : list, optional
        Hydrocarbon flag of components, default True for all components

    Returns
    -------
    Vcm : float
        Critical volume of mixture, [m³/kg]

    Examples
    --------
    Example from [2]_; 63% nC4 37% nC7

    >>> Vc1 = unidades.SpecificVolume(0.0704, "ft3lb")
    >>> Vc2 = unidades.SpecificVolume(0.0691, "ft3lb")
    >>> Mi = [58.12, 100.2]
    >>> Mm = Mi[0]*0.63+Mi[1]*0.37
    >>> "%0.2f" % (Vc_ChuehPrausnitz([0.63, 0.37], [Vc1, Vc2], Mi).ft3lb*Mm)
    '4.35'
    """
    # Define default C parameters:
    if hydrocarbon is None:
        hydrocarbon = [True]*len(xi)

    # Convert critical volumes to molar base
    Vci = [Vc*M for Vc, M in zip(Vci, Mi)]

    Mm = sum([M*x for M, x in zip(Mi, xi)])

    Vij = []
    for Vc_i, C_i in zip(Vci, hydrocarbon):
        Viji = []
        for Vc_j, C_j in zip(Vci, hydrocarbon):
            if C_i and C_j:
                C = 0
            else:
                C = 0.1559
            mu = abs((Vc_i-Vc_j)/(Vc_i+Vc_j))
            Viji.append(-1.4684*mu + C)
        Vij.append(Viji)

    nuij = []
    for Viji, Vc_i in zip(Vij, Vci):
        nuiji = []
        for V, Vc_j in zip(Viji, Vci):
            nuiji.append(V*(Vc_i+Vc_j)/2)
        nuij.append(nuiji)

    # Eq 2
    phii = []
    for x_j, Vc_j in zip(xi, Vci):
        suma = sum([x_i*Vc_i**(2/3) for x_i, Vc_i in zip(xi, Vci)])
        phii.append(x_j*Vc_j**(2/3)/suma)

    # Eq 4 generalized
    sum1 = sum([phi*Vc for phi, Vc in zip(phii, Vci)])
    sum2 = 0
    for phi_i, nuiji in zip(phii, nuij):
        for phi_j, nu in zip(phii, nuiji):
            sum2 += phi_i*phi_j*nu
    Vcm = sum1 + sum2

    return unidades.SpecificVolume(Vcm/Mm)


# Liquid density correlations
@refDoc(__doi__, [9, 2, 3])
def RhoL_RackettMix(T, xi, Tci, Pci, Vci, Zrai, Mi):
    r"""Calculates saturated liquid densities of muxteres using the
    modified Rackett equation by Spencer-Danner, also referenced in API
    procedure 6A3.1 pag 479

    .. math::
        \frac{1}{\rho_{bp}} = R\left(\sum_i x_i \frac{T_{ci}}{P_{ci}}\right)
        Z_{RAm}^{1+\left(1+T_r\right)^{2/7}}

    .. math::
        Z_{RAm} = \sum_i x_iZ_{RAi}

    .. math::
        T_{cm} = \sum_i \phi_iT_{ci}

    .. math::
        \phi_i = \frac{x_iV_{ci}}{\sum_i x_iV_{ci}}

    .. math::
        1-k_{ij} = \frac{8\left(V_{ci}V{cj}\right)^{1/2}}
        {\left(V_{ci}^{1/3}+V_{cj}^{1/3}\right)^3}

    .. math::
        T_{cij} = \left(1+k{ij}\right)\left(T_{ci}T_{cj}\right)^{1/2}

    Parameters
    ----------
    T : float
        Temperature, [K]
    xi : list
        Mole fractions of components, [-]
    Tci : list
        Critical temperature of components, [K]
    Pci : list
        Critical pressure of components, [Pa]
    Vci : list
        Critical volume of components, [m³/kg]
    Zrai : list
        Racket constant of components, [-]
    Mi : list
        Molecular weight of components, [g/mol]

    Returns
    -------
    rho : float
        Bubble point liquid density at T, [kg/m³]

    Examples
    --------
    Example from [2]_; 58.71% ethane, 41.29% heptane at 91ºF

    >>> T = unidades.Temperature(91, "F")
    >>> xi = [0.5871, 0.4129]
    >>> Tc1 = unidades.Temperature(89.92, "F")
    >>> Tc2 = unidades.Temperature(512.7, "F")
    >>> Pc1 = unidades.Pressure(706.5, "psi")
    >>> Pc2 = unidades.Pressure(396.8, "psi")
    >>> Vc1 = unidades.SpecificVolume(0.0788, "ft3lb")
    >>> Vc2 = unidades.SpecificVolume(0.0691, "ft3lb")
    >>> Zrai = [0.2819, 0.261]
    >>> Mi = [30.07, 100.205]
    >>> args = (T, xi, [Tc1, Tc2], [Pc1, Pc2], [Vc1, Vc2], Zrai, Mi)
    >>> "%0.2f" % RhoL_RackettMix(*args).kgl
    '0.56'

    Example 5-3 from [3]_; 70% ethane, 30% nC10 at 344.26K

    >>> xi = [0.7, 0.3]
    >>> Tci = [305.32, 617.7]
    >>> Pci = [48.72e5, 21.1e5]
    >>> Vci = [145.5/1000, 624/1000]
    >>> Zrai = [0.282, 0.247]
    >>> Mi = [1, 1]
    >>> args = (344.26, [0.7, 0.3], Tci, Pci, Vci, Zrai, Mi)
    >>> "%0.1f" % (1/RhoL_RackettMix(*args).gcc)
    '120.0'
    """
    # Convert critical volumes to molar base
    Vci = [Vc*M for Vc, M in zip(Vci, Mi)]

    # Eq 11
    Zram = 0
    for x, Zra in zip(xi, Zrai):
        Zram += x*Zra

    # Eq 13
    suma = 0
    for x, Vc in zip(xi, Vci):
        suma += x*Vc
    phi = []
    for x, Vc in zip(xi, Vci):
        phi.append(x*Vc/suma)

    # Eq 16
    kij = []
    for Vc_i in Vci:
        kiji = []
        for Vc_j in Vci:
            kiji.append(1-8*((Vc_i**(1/3)*Vc_j**(1/3))**0.5/(
                Vc_i**(1/3)+Vc_j**(1/3)))**3)
        kij.append(kiji)

    # Eq 15
    Tcij = []
    for Tc_i, kiji in zip(Tci, kij):
        Tciji = []
        for Tc_j, k in zip(Tci, kiji):
            Tciji.append((Tc_i*Tc_j)**0.5*(1-k))
        Tcij.append(Tciji)

    # Eq 14
    Tcm = 0
    for phi_i, Tciji in zip(phi, Tcij):
        for phi_j, Tc in zip(phi, Tciji):
            Tcm += phi_i*phi_j*Tc

    # Eq 10
    suma = 0
    for x, Tc, Pc in zip(xi, Tci, Pci):
        suma += x*Tc/Pc*101325
    Tr = T/Tcm
    V = R_atml*suma*Zram**(1+(1-Tr)**(2/7))

    Mm = sum([M*x for M, x in zip(Mi, xi)])
    return unidades.Density(Mm/V)


@refDoc(__doi__, [10, 2, 3])
def RhoL_CostaldMix(T, xi, Tci, wi, Vci, Mi):
    r"""Calculates saturated liquid densities of pure components using the
    Corresponding STAtes Liquid Density (COSTALD) method, developed by
    Hankinson and Thomson, referenced too in API procedure 6A3.2 pag. 482

    .. math::
        T_{cm} = \frac{\sum_i\sum_jx_ix_j\left(V_i^oT_{ci}V_j^oT_{cj}\right)
        ^{1/2}}{V_m^o}

    .. math::
        V_m^o = \frac{\sum_i xiV_i^o + 3 \sum_i x_iV_i^{o^{2/3}}
        \sum_i x_iV_i^{o^{1/3}}}{4}

    .. math::
        \omega_m = \sum_i x_i\omega_{SRKi}

    Parameters
    ----------
    T : float
        Temperature [K]
    xi : list
        Mole fractions of components, [-]
    Tci : list
        Critical temperature of components, [K]
    wi : list
        Acentric factor optimized to SRK of components, [-]
    Vci : list
        Characteristic volume of components, [m³/kg]
    Mi : list
        Molecular weight of components, [g/mol]

    Returns
    -------
    rho : float
        Bubble point liquid density at T, [kg/m³]

    Examples
    --------
    Example from [2]_; 20% methane, 80% nC10 at 160ºF

    >>> T = unidades.Temperature(160, "F")
    >>> Tc1 = unidades.Temperature(-116.67, "F")
    >>> Tc2 = unidades.Temperature(652, "F")
    >>> Vc1 = unidades.SpecificVolume(1.592/16.04, "ft3lb")
    >>> Vc2 = unidades.SpecificVolume(9.919/142.28, "ft3lb")
    >>> Mi = [16.04, 142.28]
    >>> args = (T, [0.2, 0.8], [Tc1, Tc2], [0.0074, 0.4916], [Vc1, Vc2], Mi)
    >>> "%0.3f" % RhoL_CostaldMix(*args).kgl
    '0.667'

    Example 5-3 from [3]_; 70% ethane, 30% nC10 at 344.26K

    >>> xi = [0.7, 0.3]
    >>> Tci = [305.32, 617.7]
    >>> Vci = [145.5/1000, 624/1000]
    >>> wi = [0.099, 0.491]
    >>> Mi = [1, 1]
    >>> args = (344.26, [0.7, 0.3], Tci, wi, Vci, Mi)
    >>> "%0.1f" % (1/RhoL_CostaldMix(*args).gcc)
    '119.5'
    """
    # Convert critical volumes to molar base
    Vci = [Vc*M for Vc, M in zip(Vci, Mi)]

    # Apply mixing rules
    # Eq 24
    wm = 0
    for x, w in zip(xi, wi):
        wm += x*w

    # Eq 21
    sum1 = 0
    sum2 = 0
    sum3 = 0
    for x, Vc, Tc in zip(xi, Vci, Tci):
        sum1 += x*Vc
        sum2 += x*Vc**(2/3)
        sum3 += x*Vc**(1/3)
    Vcm = (sum1+3*sum2*sum3)/4

    # Eq 19 & 21
    Tcm = 0
    for x_i, Vc_i, Tc_i in zip(xi, Vci, Tci):
        for x_j, Vc_j, Tc_j in zip(xi, Vci, Tci):
            Tcm += x_i*x_j*(Vc_i*Tc_i*Vc_j*Tc_j)**0.5
    Tcm /= Vcm

    Mm = sum([M*x for M, x in zip(Mi, xi)])

    # Apply the pure component procedure with the mixing parameters
    rho = RhoL_Costald(T, Tcm, wm, Vcm)
    return unidades.Density(rho*Mm)


@refDoc(__doi__, [11, 3])
def RhoL_AaltoKeskinenMix(T, P, xi, Tci, Pci, Vci, wi, Mi, rhos):
    r"""Calculates compressed-liquid density of a mixture, using the
    Aalto-Keskinen modification of Chang-Zhao correlation

    .. math::
        T_{cm} = \frac{\sum_i\sum_jx_ix_j\left(V_i^oT_{ci}V_j^oT_{cj}\right)
        ^{1/2}}{V_m^o}

    .. math::
        V_m^o = \frac{\sum_i xiV_i^o + 3 \sum_i x_iV_i^{o^{2/3}}
        \sum_i x_iV_i^{o^{1/3}}}{4}

    .. math::
        P_{cm} = \frac{\left(0.291-0.08\omega_{SRKm}\right)RT_{cm}}{V_{cm}}

    .. math::
        \omega_{SRKm} = \left(\sum_i x_i\omega_{SRKi}\right)^2

    Parameters
    ----------
    T : float
        Temperature, [K]
    P : float
        Pressure, [Pa]
    xi : list
        Mole fractions of components, [-]
    Tci : list
        Critical temperature of components, [K]
    Pci : list
        Critical pressure of components, [Pa]
    Vci : list
        Critical volume of components, [m³/kg]
    wi : list
        Acentric factor (SRK optimized) of components, [-]
    Mi : list
        Molecular weight of components, [g/mol]
    rhos : float
        Boiling point liquid density, [kg/m³]

    Returns
    -------
    rho : float
        High-pressure liquid density, [kg/m³]

    Examples
    --------
    Example 5-4 from [3]_; 70% ethane, 30% nC10 at 344.26K and 10000psi

    >>> P = unidades.Pressure(10000, "psi")
    >>> xi = [0.7, 0.3]
    >>> Tci = [305.32, 617.7]
    >>> Pci = [48.72e5, 21.1e5]
    >>> Vci = [145.5/1000, 624/1000]
    >>> wi = [0.099, 0.491]
    >>> Mi = [1, 1]
    >>> args = (344.26, P, xi, Tci, Pci, Vci, wi, Mi, 1/116.43*1000)
    >>> "%0.2f" % (1/RhoL_AaltoKeskinenMix(*args).gcc)
    '99.05'
    """
    # Convert critical volumes to molar base
    Vci = [Vc*M for Vc, M in zip(Vci, Mi)]

    # Apply mixing rules
    # Eq A5
    wm = 0
    for x, w in zip(xi, wi):
        # Avoid problem with square root of negative number, hydrogen, methane
        if w < 0:
            w = 0
        wm += x*w**0.5
    wm *= wm

    # Eq 2
    sum1 = 0
    sum2 = 0
    sum3 = 0
    for x, Vc, Tc in zip(xi, Vci, Tci):
        sum1 += x*Vc
        sum2 += x*Vc**(2/3)
        sum3 += x*Vc**(1/3)
    Vcm = (sum1+3*sum2*sum3)/4

    # Eq 1
    Tcm = 0
    for x_i, Vc_i, Tc_i in zip(xi, Vci, Tci):
        for x_j, Vc_j, Tc_j in zip(xi, Vci, Tci):
            Tcm += x_i*x_j*(Vc_i*Tc_i*Vc_j*Tc_j)**0.5
    Tcm /= Vcm

    # Eq 4
    Pcm = (0.291-0.08*wm)*R*Tcm/Vcm*1000

    Ps = _Pv(T, Tcm, Pcm, wm)

    # Apply the pure component procedure with the mixing parameters
    rho = RhoL_AaltoKeskinen(T, P, Tcm, Pcm, wm, Ps, rhos)
    return unidades.Density(rho)


@refDoc(__doi__, [12, 2])
def RhoL_TaitCostaldMix(T, P, xi, Tci, Vci, wi, Mi, rhos):
    r"""Calculates compressed-liquid density of a mixture, using the
    Thomson-Brobst-Hankinson modification of Chang-Zhao correlation adapted to
    mixtures, also referenced in API procedure 6A3.4, pag 489

    .. math::
        T_{cm} = \frac{\sum_i\sum_jx_ix_j\left(V_i^oT_{ci}V_j^oT_{cj}\right)
        ^{1/2}}{V_m^o}

    .. math::
        V_m^o = \frac{\sum_i xiV_i^o + 3 \sum_i x_iV_i^{o^{2/3}}
        \sum_i x_iV_i^{o^{1/3}}}{4}

    .. math::
        P_{cm} = \frac{\left(0.291-0.08\omega_{SRKm}\right)RT_{cm}}{V_{cm}}

    .. math::
        \omega_{SRKm} = \left(\sum_i x_i\omega_{SRKi}\right)^2

    .. math::
        \frac{P_s}{P_{cm}} = P_m^{(0)}+\omega_{SRKm}P_m^{(1)}

    .. math::
        P_m^{(0)} = 5.8031817 \log T_{rm}+0.07608141\alpha

    .. math::
        P_m^{(1)} = 4.86601\beta

    .. math::
        \alpha = 35 - \frac{36}{T_{rm}}-96.736\log T_{rm}+T_{rm}^6

    .. math::
        \beta = \log T_{rm}+0.03721754\alpha

    Parameters
    ----------
    T : float
        Temperature, [K]
    P : float
        Pressure, [Pa]
    xi : list
        Mole fractions of components, [-]
    Tci : list
        Critical temperature of components, [K]
    Vci : list
        Critical volume of components, [m³/kg]
    wi : list
        Acentric factor (SRK optimized) of components, [-]
    Mi : list
        Molecular weight of components, [g/mol]
    rhos : float
        Boiling point liquid density, [kg/m³]

    Returns
    -------
    rho : float
        High-pressure liquid density, [kg/m³]

    Examples
    --------
    Example from [2]_; 20% ethane, 80% nC10 at 166F and 3000psi

    >>> T = unidades.Temperature(160, "F")
    >>> P = unidades.Pressure(3000, "psi")
    >>> xi = [0.2, 0.8]
    >>> Tc1 = unidades.Temperature(89.92, "F")
    >>> Tc2 = unidades.Temperature(652, "F")
    >>> Tci = [Tc1, Tc2]
    >>> Mi = [30.07, 142.286]
    >>> Vc1 = unidades.SpecificVolume(2.335/Mi[0], "ft3lb")
    >>> Vc2 = unidades.SpecificVolume(9.919/Mi[1], "ft3lb")
    >>> Vci = [Vc1, Vc2]
    >>> wi = [0.0983, 0.4916]
    >>>
    >>> Vs = unidades.SpecificVolume(2.8532/119.8428, "ft3lb")
    >>> args = (T, P, xi, Tci, Vci, wi, Mi, 1/Vs)
    >>> "%0.3f" % RhoL_TaitCostaldMix(*args).gcc
    '0.698'
    """
    # Convert critical volumes to molar base
    Vci = [Vc*M for Vc, M in zip(Vci, Mi)]

    # Apply mixing rules
    # Eq 16
    wm = 0
    for x, w in zip(xi, wi):
        wm += x*w

    # Eq 15
    sum1 = 0
    sum2 = 0
    sum3 = 0
    for x, Vc, Tc in zip(xi, Vci, Tci):
        sum1 += x*Vc
        sum2 += x*Vc**(2/3)
        sum3 += x*Vc**(1/3)
    Vcm = (sum1+3*sum2*sum3)/4

    # Eq 13
    Tcm = 0
    for x_i, Vc_i, Tc_i in zip(xi, Vci, Tci):
        for x_j, Vc_j, Tc_j in zip(xi, Vci, Tci):
            Tcm += x_i*x_j*(Vc_i*Tc_i*Vc_j*Tc_j)**0.5
    Tcm /= Vcm

    # Eq 17
    Pcm = (0.291-0.08*wm)*R*Tcm/Vcm*1000

    # Saturation Pressure
    Ps = _Pv(T, Tcm, Pcm, wm)

    # Apply the pure component procedure with the mixing parameters
    rho = RhoL_TaitCostald(T, P, Tcm, Pcm, wm, Ps, rhos)
    return unidades.Density(rho)


@refDoc(__doi__, [13])
def RhoL_NasrifarMix(T, P, xi, Tci, Vci, wi, Mi, rhos):
    r"""Calculates compressed-liquid density of a mixture, using the
    Nasrifar correlation

    .. math::
        T_{cm} = \frac{\sum_i\sum_jx_ix_j\left(V_i^oT_{ci}V_j^oT_{cj}\right)
        ^{1/2}}{V_m^o}

    .. math::
        V_m^o = \frac{\sum_i xiV_i^o + 3 \sum_i x_iV_i^{o^{2/3}}
        \sum_i x_iV_i^{o^{1/3}}}{4}

    .. math::
        P_{cm} = \frac{\left(0.291-0.08\omega_{SRKm}\right)RT_{cm}}{V_{cm}}

    .. math::
        \omega_{SRKm} = \sum_i x_i\omega_{SRKi}

    Parameters
    ----------
    T : float
        Temperature, [K]
    P : float
        Pressure, [Pa]
    xi : list
        Mole fractions of components, [-]
    Tci : list
        Critical temperature of components, [K]
    Vci : list
        Critical volume of components, [m³/kg]
    wi : list
        Acentric factor (SRK optimized) of components, [-]
    Mi : list
        Molecular weight of components, [g/mol]
    rhos : float
        Boiling point liquid density, [kg/m³]

    Returns
    -------
    rho : float
        High-pressure liquid density, [kg/m³]
    """
    # Convert critical volumes to molar base
    Vci = [Vc*M for Vc, M in zip(Vci, Mi)]

    # Apply mixing rules
    # Eq 25
    wm = 0
    for x, w in zip(xi, wi):
        wm += x*w

    # Eq 24
    sum1 = 0
    sum2 = 0
    sum3 = 0
    for x, Vc, Tc in zip(xi, Vci, Tci):
        sum1 += x*Vc
        sum2 += x*Vc**(2/3)
        sum3 += x*Vc**(1/3)
    Vcm = (sum1+3*sum2*sum3)/4

    # Eq 22
    Tcm = 0
    for x_i, Vc_i, Tc_i in zip(xi, Vci, Tci):
        for x_j, Vc_j, Tc_j in zip(xi, Vci, Tci):
            Tcm += x_i*x_j*(Vc_i*Tc_i*Vc_j*Tc_j)**0.5
    Tcm /= Vcm

    Mm = 0
    for x, M in zip(xi, Mi):
        Mm += x*M

    # Eq 26
    Pcm = (0.291-0.08*wm)*R*Tcm/Vcm*1000

    # Saturation Pressure
    Ps = _Pv(T, Tcm, Pcm, wm)

    # Apply the pure component procedure with the mixing parameters
    rho = RhoL_Nasrifar(T, P, Tcm, Pcm, wm, Mm, Ps, rhos)
    return unidades.Density(rho)


@refDoc(__doi__, [14, 2])
def RhoL_APIMix(T, P, xi, Tci, Pci, rhos, To=None, Po=None):
    r"""Calculates compressed-liquid density, using the analytical expression
    of Lu Chart referenced in API procedure 6A2.22

    .. math::
        \rho_2 = \rho_1\frac{C_2}{C_1}

    Parameters
    ----------
    T : float
        Temperature, [K]
    P : float
        Pressure, [Pa]
    xi : list
        Mole fractions of components, [-]
    Tci : list
        Critical temperature of components, [K]
    Pci : list
        Critical pressure of components, [Pa]
    rhos : float
        Liquid density at 60ºF, [kg/m^3]
    To : float, optional
        Reference temperature with known density, [K]
    Po : float, optional
        Reference pressure with known density, [Pa]

    Returns
    -------
    rho : float
        High-pressure liquid density, [kg/m^3]

    Examples
    --------
    Example A from [2]; 60.52% Ethylene, 39.48% n-heptane at 162.7ºF and 900psi

    >>> T = unidades.Temperature(162.7, "F")
    >>> P = unidades.Pressure(900, "psi")
    >>> x = [0.6052, 0.3948]
    >>> Tc = unidades.Temperature(48.58, "F")
    >>> Tc2 = unidades.Temperature(512.7, "F")
    >>> Pc = unidades.Pressure(729.8, "psi")
    >>> Pc2 = unidades.Pressure(396.8, "psi")
    >>> rs = unidades.Density(37.55, "lbft3")
    >>> To = unidades.Temperature(49, "F")
    >>> Po = unidades.Pressure(400, "psi")
    >>> "%0.1f" % RhoL_APIMix(T, P, x, [Tc, Tc2], [Pc, Pc2], rs, To, Po).lbft3
    '31.5'

    Example C from [2]; 20% methane, 80% n-C10 at 160ºF and 3000psi

    >>> T = unidades.Temperature(162.7, "F")
    >>> P = unidades.Pressure(900, "psi")
    >>> x = [0.2, 0.8]
    >>> Tc1 = unidades.Temperature(-116.63, "F")
    >>> Tc2 = unidades.Temperature(652.1, "F")
    >>> Pc1 = unidades.Pressure(667.8, "psi")
    >>> Pc2 = unidades.Pressure(304, "psi")
    >>> rs = unidades.Density(41.65, "lbft3")
    >>> "%0.1f" % RhoL_APIMix(T, P, x, [Tc1, Tc2], [Pc1, Pc2], rs).lbft3
    '42.6'
    """
    # Define reference state
    if To is None:
        To = T
    if Po is None:
        Po = 101325

    # Calculate pseudocritical properties
    Tc = 0
    Pc = 0
    for x, Tc_i, Pc_i in zip(xi, Tci, Pci):
        Tc += x*Tc_i
        Pc += x*Pc_i

    def C(Tr, Pr):
        """Polinomial ajust of Lu chart, referenced in [14]"""
        A0 = 1.6368-0.04615*Pr+2.1138e-3*Pr**2-0.7845e-5*Pr**3-0.6923e-6*Pr**4
        A1 = -1.9693+0.21874*Pr-8.0028e-3*Pr**2-8.2328e-5*Pr**3+5.2604e-6*Pr**4
        A2 = 2.4638-.36461*Pr+12.8763e-3*Pr**2+14.8059e-5*Pr**3-8.6895e-6*Pr**4
        A3 = -1.5841+.25136*Pr-11.3805e-3*Pr**2+9.5672e-5*Pr**3+2.1812e-6*Pr**4
        C = A0 + A1*Tr + A2*Tr**2 + A3*Tr**3
        return C

    C1 = C(To/Tc, Po/Pc)
    C2 = C(T/Tc, P/Pc)

    # Eq 1
    d2 = rhos*C2/C1
    return unidades.Density(d2)


@refDoc(__doi__, [11, 12])
def _Pv(T, Tcm, Pcm, wm):
    """Pseudo vapor presure from pseudocritical properties to use in
    compressed liquid density correlations

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tcm : float
        Pseudo critical temperature of mixture, [K]
    Pcm : float
        Pseudo critical pressure of mixture, [Pa]
    wm : float
        Acentric factor of mixture, [-]

    Returns
    -------
    Pbp : float
        Boiling point pressure, [Pa]
    """
    Trm = T/Tcm
    alpha = 35 - 36/Trm - 96.736*log10(Trm) + Trm**6                   # Eq 22
    beta = log10(Trm) + 0.03721754*alpha                               # Eq 23
    Pro = 5.8031817*log10(Trm)+0.07608141*alpha                        # Eq 20

    # Using Aalto modificated equation Eq 8
    Pr1 = 4.86601*beta+0.03721754*alpha                                # Eq 21

    Pbp = Pcm*10**(Pro+wm*Pr1)
    return unidades.Pressure(Pbp)


# Liquid viscosity correlations
@refDoc(__doi__, [20, 2, 3])
def MuL_KendallMonroe(xi, mui):
    r"""Calculate viscosity of liquid mixtures using the Kendall-Monroe method,
    also referenced in API procedure 11A3.1, pag 1051

    .. math::
        \mu_m = \left(\sum_{i=1}^nx_i\mu_i^{1/3}\right)^3

    Parameters
    ----------
    xi : list
        Mole fractions of components, [-]
    mui : list
        Viscosities of components, [Pa·s]

    Returns
    -------
    mu : float
        Viscosity of mixture, [Pa·s]

    Examples
    --------
    Example A from [2]_; 29.57% nC16, 35.86% benzene, 34.57% nC6 at 77ºF

    >>> x = [0.2957, 0.3586, 0.3457]
    >>> mu = [3.03e-3, 0.6e-3, 0.3e-3]
    >>> "%0.2f" % MuL_KendallMonroe(x, mu).cP
    '0.89'

    Example B from [2]_; 25% nC3, 50% nC5, 25% cycloC6 at 160ºF
    >>> x = [0.25, 0.5, 0.25]
    >>> mu = [0.109e-3, 0.218e-3, 0.63e-3]
    >>> "%0.3f" % MuL_KendallMonroe(x, mu).cP
    '0.256'
    """
    mum = 0
    for x, mu in zip(xi, mui):
        mum += x*mu**(1/3)
    return unidades.Viscosity(mum**3)


def MuL_Chemcad(xi, Mi, mui):
    r"""Calculate viscosity of liquid mixtures using the mixing rules of
    Arrhenius and with modification used in CHEMCAD®

    .. math::
        \mu_m = e^A

    .. math::
        A = \sum_i \frac{x_iM_i}{M_m}\log\mu_i if M > 2000

    .. math::
        A = \sum_i x_i\log\mu_i if M < 2000

    Parameters
    ----------
    xi : list
        Mole fractions of components, [-]
    Mi : list
        Molecular weight of components, [g/mol]
    mui : list
        Viscosities of components, [Pa·s]

    Returns
    -------
    mu : float
        Viscosity of mixture, [Pa·s]

    Notes
    -----
    I can't find any reference for this mixing rule. In some place is refered
    as the Arrhenius mixing rules, but only for the los molecular weight form.
    """
    Mm = sum([M*x for M, x in zip(Mi, xi)])

    if max(Mi) > 2000:
        A = 0
        for x, M, mu in zip(xi, Mi, mui):
            A += x*M/Mm*log(mu)
    else:
        A = 0
        for x, mu in zip(xi, mui):
            A += x*log(mu)

    return unidades.Viscosity(exp(A))


# Gas viscosity correlations
@refDoc(__doi__, [3])
def MuG_Reichenberg(T, xi, Tci, Pci, Mi, mui, Di):
    r"""Calculate viscosity of gas mixtures using the Reichenberg method as
    explain in [3]_

    .. math::
        \eta_m = \sum _{i=1}^n K_i\left(1+2\sum_{j=1}^{i-1}H_{ij}K_j +
        \sum_{j=1≠i}^n \sum_{k=1≠i}^nH_{ij}H_{ik}K_jK_k\right)

    .. math::
        K_i = \frac{x_i\eta_i}{x_i+\eta_i \sum_{k=1≠i}^n x_kH_{ik}
        \left[3+\left(2M_k/M_i\right)\right]}

    .. math::
        H_{ij} = \left[\frac{M_iM_j}{32\left(M_i+M_j\right)^3}\right]^{1/2}
        \left(C_i+C_j\right)^2 \frac{\left[1+0.36T_{rij}\left(T_{rij}-1
        \right)\right]^{1/6}F_{Rij}}{T_{rij}^{1/2}}

    .. math::
        F_{Ri} = \frac{T_{ri}^{3.5}+\left(10\mu_{ri}\right)^7}
        {T_{ri}^{3.5}\left[1+\left(10\mu_{ri}\right)^7\right]}

    .. math::
        C_i = \frac{M_i^{1/4}}{\left(\eta_iU_i\right)^{1/2}}

    .. math::
        U_i = \frac{\left[1+0.36T_{ri}\left(T_{ri}-1\right)\right]^{1/6}F_{Ri}}
        {T_{ri}^{1/2}}

    .. math::
        T_{rij} = \frac{T}{\left(T_{ci}T_{cj}\right)^{1/2}}

    .. math::
        \mu_{rij} = \left(\mu_{ri}\mu_{rj}\right)^{1/2}

    .. math::
        \mu_r = 52.46\frac{\mu^2P_c}{T_c^2}


    Parameters
    ----------
    T : float
        Temperature, [K]
    xi : list
        Mole fractions of components, [-]
    Tci : list
        Critical temperature of components, [K]
    Pci : list
        Critical pressure of components, [Pa]
    Mi : list
        Molecular weights of components, [g/mol]
    mui : list
        Viscosities of components, [Pa·s]
    Di : list
        Dipole moment of components, [Debye]

    Returns
    -------
    mu : float
        Viscosity of mixture, [Pa·s]

    Examples
    --------
    Example 9-4 from [3]_; 28.6% N2, 71.4% R22 at 50ºC

    >>> T = unidades.Temperature(50, "C")
    >>> x = [0.286, 0.714]
    >>> Tc = [126.2, 369.28]
    >>> Pc = [33.98e5, 49.86e5]
    >>> M = [28.014, 86.468]
    >>> mu = [188e-7, 134e-7]
    >>> D = [0, 1.4]
    >>> "%0.1f" % MuG_Reichenberg(T, x, Tc, Pc, M, mu, D).microP
    '146.2'
    """
    # Calculate reduced temperatures
    Tri = [T/Tc for Tc in Tci]
    Trij = []
    for Tc_i in Tci:
        Triji = []
        for Tc_j in Tci:
            Triji.append(T/(Tc_i*Tc_j)**0.5)
        Trij.append(Triji)

    # Calculate reduced viscosity
    muri = [52.46*D**2*Pc*1e-5/Tc**2 for D, Pc, Tc in zip(Di, Pci, Tci)]
    murij = []
    for mur_i in muri:
        muriji = []
        for mur_j in muri:
            muriji.append((mur_i*mur_j)**0.5)
        murij.append(muriji)

    # Polar correction, Eq 9-5.5
    Fri = []
    for Tr, mur in zip(Tri, muri):
        Fri.append((Tr**3.5+(10*mur)**7)/Tr**3.5/(1+(10*mur)**7))
    Frij = []
    for Triji, muriji in zip(Trij, murij):
        Friji = []
        for Tr, mur in zip(Triji, muriji):
            Friji.append((Tr**3.5+(10*mur)**7)/Tr**3.5/(1+(10*mur)**7))
        Frij.append(Friji)

    # Eq 9-5.3
    Ui = []
    for Tr, Fr in zip(Tri, Fri):
        Ui.append((1+0.36*Tr*(Tr-1))**(1/6)*Fr/Tr**0.5)

    # Eq 9-5.4
    Ci = [M**0.25/(mu*U)**0.5 for M, mu, U in zip(Mi, mui, Ui)]

    # Eq 9-5.6
    Hij = []
    for M_i, C_i, Triji, Friji in zip(Mi, Ci, Trij, Frij):
        Hiji = []
        for M_j, C_j, Tr, Fr in zip(Mi, Ci, Triji, Friji):
            Hiji.append((M_i*M_j/32/(M_i+M_j)**3)**0.5*(C_i+C_j)**2*(
                        1+0.36*Tr*(Tr-1))**(1/6)*Fr/Tr**0.5)
        Hij.append(Hiji)

    # Eq 9-5.2
    sumai = []
    for i, M_i in enumerate(Mi):
        sumaij = 0
        for j, M_j in enumerate(Mi):
            if i != j:
                sumaij += xi[j]*Hij[i][j]*(3+2*M_j/M_i)
        sumai.append(sumaij)
    Ki = [x*mu/(x+mu*suma) for x, mu, suma in zip(xi, mui, sumai)]

    # Eq 9-5.1
    mu = 0
    for i, K_i in enumerate(Ki):
        sum1 = 0
        for H, K_j in zip(Hij[i], Ki[:i]):
            sum1 += H*K_j

        sum2 = 0
        for j, K_j in enumerate(Ki):
            for k, K_k in enumerate(Ki):
                if j != i and k != i:
                    sum2 += Hij[i][j]*Hij[i][k]*K_j*K_k

        mu += K_i*(1+2*sum1+sum2)

    return unidades.Viscosity(mu)


@refDoc(__doi__, [1, 2, 3])
def MuG_Wilke(xi, Mi, mui):
    r"""Calculate viscosity of gas mixtures using the Wilke mixing rules, also
    referenced in API procedure 11B2.1, pag 1102

    .. math::
        \mu_{m} = \sum_{i=1}^n \frac{\mu_i}{1+\frac{1}{x_i}\sum_{j=1}
        ^n x_j \phi_{ij}}

    .. math::
        \phi_{ij} = \frac{\left[1+\left(\mu_i/\mu_j\right)^{0.5}(M_j/M_i)
        ^{0.25}\right]^2}{4/\sqrt{2}\left[1+\left(M_i/M_j\right)\right]
        ^{0.5}}

    Parameters
    ----------
    xi : list
        Mole fractions of components, [-]
    Mi : list
        Molecular weights of components, [g/mol]
    mui : list
        Viscosities of components, [Pa·s]

    Returns
    -------
    mu : float
        Viscosity of mixture, [Pa·s]

    Examples
    --------
    Selected point in Table II from [1]_
    CO2 3.5%, O2 0.3%, CO 27.3%, H2 14.4%, CH4 3.7%, N2 50%, C2 0.8%

    >>> from numpy import r_
    >>> from lib.mEoS import CO2, O2, CO, H2, CH4, N2, C2
    >>> Mi = [CO2.M, O2.M, CO.M, H2.M, CH4.M, N2.M, C2.M]
    >>> xi = [0.035, 0.003, 0.273, 0.144, 0.037, 0.5, 0.008]
    >>> mui = r_[147.2, 201.9, 174.9, 87.55, 108.7, 178.1, 90.9]
    >>> "%0.7f" % MuG_Wilke(xi, Mi, mui*1e-9).dynscm2
    '0.0001711'

    Example A from [2]_; 58.18% H2, 41.82% propane at 77ºF and 14.7 psi

    >>> mu = MuG_Wilke([0.5818, 0.4182], [2.02, 44.1], [8.91e-6, 8.22e-6])
    >>> "%0.4f" % mu.cP
    '0.0092'

    Example B from [2]_ 95.6% CH4, 3.6% C2, 0.5% C3, 0.3% N2

    >>> x = [0.956, 0.036, 0.005, 0.003]
    >>> Mi = [16.04, 30.07, 44.1, 28.01]
    >>> mui = [1.125e-5, 9.5e-6, 8.4e-6, 1.79e-5]
    >>> "%0.5f" % MuG_Wilke(x, Mi, mui).cP
    '0.01117'

    Example 9-5 from [3]_, methane 69.7%, n-butane 30.3%

    >>> mui = [109.4e-7, 72.74e-7]
    >>> "%0.2f" % MuG_Wilke([0.697, 0.303], [16.043, 58.123], mui).microP
    '92.25'
    """
    # Eq 4
    kij = []
    for i, mi in enumerate(Mi):
        kiji = []
        for j, mj in enumerate(Mi):
            if i == j:
                kiji.append(0)
            else:
                kiji.append((1+(mui[i]/mui[j])**0.5*(mj/mi)**0.25)**2 /
                            8**0.5/(1+mi/mj)**0.5)
        kij.append(kiji)

    # Eq 13
    # Precalculate of internal sum
    suma = []
    for i, Xi in enumerate(xi):
        sumai = 0
        if Xi != 0:
            for j, Xj in enumerate(xi):
                sumai += kij[i][j]*Xj/Xi
        suma.append(sumai)

    mu = 0
    for mu_i, sumai in zip(mui, suma):
        mu += mu_i/(1+sumai)
    return unidades.Viscosity(mu)


@refDoc(__doi__, [3])
def MuG_Herning(xi, Mi, mui):
    r"""Calculate viscosity of gas mixtures using the Herning-Zipperer mixing
    rules, as explain in [3]_

    .. math::
        \mu_{m} = \sum_{i=1}^n \frac{\mu_i}{1+\frac{1}{x_i}\sum_{j=1}
        ^n x_j \phi_{ij}}

    .. math::
        \phi_{ij} = \left(\frac{M_j}{M_i}\right)^{1/2}

    Parameters
    ----------
    xi : list
        Mole fractions of components, [-]
    Mi : list
        Molecular weights of components, [g/mol]
    mui : list
        Viscosities of components, [Pa·s]

    Returns
    -------
    mu : float
        Viscosity of mixture, [Pa·s]

    Examples
    --------
    Example 9-6 from [3]_, methane 69.7%, n-butane 30.3%

    >>> mui = [109.4e-7, 72.74e-7]
    >>> "%0.1f" % MuG_Herning([0.697, 0.303], [16.043, 58.123], mui).microP
    '92.8'
    """
    kij = []
    for i, mi in enumerate(Mi):
        kiji = []
        for j, mj in enumerate(Mi):
            if i == j:
                kiji.append(0)
            else:
                kiji.append((mj/mi)**0.5)
        kij.append(kiji)

    suma = []
    for i, Xi in enumerate(xi):
        sumai = 0
        if Xi != 0:
            for j, Xj in enumerate(xi):
                sumai += kij[i][j]*Xj/Xi
        suma.append(sumai)

    mu = 0
    for mu_i, sumai in zip(mui, suma):
        mu += mu_i/(1+sumai)
    return unidades.Viscosity(mu)


@refDoc(__doi__, [3])
def MuG_Lucas(T, P, xi, Tci, Pci, Vci, Zci, Mi, Di):
    r"""Calculate the viscosity of a gas mixture using the Lucas mixing rules
    as explain in [3]_.

    .. math::
        T_{cm} = \sum_i x_iT_{ci}

    .. math::
        P_{cm} = RT_{cm}\frac{\sum_i x_iZ_{ci}}{\sum_i x_iV_{ci}}

    .. math::
        M_m = \sum_i x_iM_i

    .. math::
        F_{Pm}^o = \sum_i x_iF_{Pi}^o

    .. math::
        F_{Qm}^o = \left(\sum_i x_iF_{Qi}^o\right)A

    .. math::
        A = 1-0.01\left(\frac{M_H}{M_L}\right)^{0.87}

    Parameters
    ----------
    T : float
        Temperature, [K]
    P : float
        Pressure, [Pa]
    xi : list
        Mole fractions of components, [-]
    Tci : list
        Critical temperature of components, [K]
    Pci : list
        Critical pressure of components, [Pa]
    Vci : list
        Critical volume of components, [m³/kg]
    Zci : list
        Critical compressibility factor of components, [-]
    Mi : list
        Molecular weights of components, [g/mol]
    Di : list
        Dipole moment of components, [Debye]

    Returns
    -------
    mu : float
        Viscosity of gas, [Pa·s]

    Examples
    --------
    Example 9-7 in [3]_, 67.7% NH3 33.3% H2 at 33ºC

    >>> Tc = [405.5, 33.2]
    >>> Pc = [113.5e5, 13e5]
    >>> Vc =[72.5/17.031/1000, 64.3/2.016/1000]
    >>> Zc = [0.244, 0.306]
    >>> M = [17.031, 2.0158]
    >>> D = [1.47, 0]
    >>> mu = MuG_Lucas(33+273.15, 101325, [0.677, 0.323], Tc, Pc, Vc, Zc, M, D)
    >>> "%0.1f" % mu.microP
    '116.3'
    """
    # Use critical volume in molar base
    Vci = [Vc*M*1000 for Vc, M in zip(Vci, Mi)]

    # Calculate critical properties of mixture
    Tcm = sum([x*Tc for x, Tc in zip(xi, Tci)])
    Zcm = sum([x*Zc for x, Zc in zip(xi, Zci)])
    Vcm = sum([x*Vc for x, Vc in zip(xi, Vci)])
    Mm = sum([x*M for x, M in zip(xi, Mi)])

    Pcm = R*Tcm*Zcm/Vcm*1e6

    Trm = T/Tcm
    Prm = P/Pcm

    Pcm_bar = Pcm*1e-5
    X = 0.176*Tcm**(1/6)/Mm**0.5/Pcm_bar**(2/3)

    # Polarity and quantum effects correction factors
    Fpoi = []
    Fqoi = []
    for Tc, Pc, Zc, M, D in zip(Tci, Pci, Zci, Mi, Di):
        Tr = T/Tc
        Pc_bar = Pc*1e-5
        mur = 52.46*D**2*Pc_bar/Tc**2
        if mur < 0.022:
            Fpoi.append(1)
        elif mur < 0.075:
            Fpoi.append(1 + 30.55*(0.292-Zc)**1.72)
        else:
            Fpoi.append(1 + 30.55*(0.292-Zc)**1.72*abs(0.96+0.1*(Tr-0.7)))

        if Tr < 12:
            sign = -1
        else:
            sign = 1
        if M == 2.0158:
            Q = 0.76  # Hydrogen
            Fqoi.append(1.22*Q**0.15*(1+0.00385*((Tr-12)**2)**(1/M)*sign))
        elif M == 4.0026:
            Q = 1.38  # Helium
            Fqoi.append(1.22*Q**0.15*(1+0.00385*((Tr-12)**2)**(1/M)*sign))
        else:
            Fqoi.append(1)

    Fpom = sum([x*Fpo for x, Fpo in zip(xi, Fpoi)])
    Fqom = sum([x*Fqo for x, Fqo in zip(xi, Fqoi)])

    # Calculate A factor
    Mh = max(Mi)
    Ml = min(Mi)
    xh = xi[Mi.index(Mh)]
    if Mh/Ml > 9 and 0.05 < xh < 0.7:
        A = 1 - 0.01*(Mh/Ml)**0.87
    else:
        A = 1
    Fqom *= A

    Z1 = Fpom*Fqom*(0.807*Trm**0.618 - 0.357*exp(-0.449*Trm) +
                    0.34*exp(-4.058*Trm) + 0.018)

    if Prm < 0.6:
        # Low pressure correlation
        mu = Z1/X

    else:
        # High pressure correlation
        if Trm <= 1:
            alfa = 3.262 + 14.98*Prm**5.508
            beta = 1.39 + 14.98*Prm
            Z2 = 0.6 + 0.76*Prm**alfa + (6.99*Prm**beta-0.6)*(1-Trm)
        else:
            a = 1.245e-3/Trm*exp(5.1726*Trm**-0.3286)
            b = a*(1.6553*Trm-1.2723)
            c = 0.4489/Trm*exp(3.0578*Trm**-37.7332)
            d = 1.7368/Trm*exp(2.231*Trm**-7.6351)
            e = 1.3088
            f = 0.9425*exp(-0.1853*Trm**0.4489)
            Z2 = Z1*(1+a*Prm**e/(b*Prm**f+1/(1+c*Prm**d)))

        Y = Z2/Z1
        Fp = (1+(Fpom-1)/Y**3)/Fpom
        Fq = (1+(Fqom-1)*(1/Y-0.007*log(Y)**4))/Fqom
        mu = Z2*Fp*Fq/X

    return unidades.Viscosity(mu, "microP")


@refDoc(__doi__, [15, 3])
def MuG_Chung(T, xi, Tci, Vci, Mi, wi, Di, ki):
    r"""Calculate the viscosity of a gas mixture using the Chung correlation

    .. math::
        \mu=40.785\frac{F_{cm}\left(M_mT\right)^{1/2}}{V_{cm}^{2/3}\Omega_v}

    .. math::
        \sigma_i = 0.809V_{ci}^{1/3}

    .. math::
        \sigma_{ij} = \xi\left(\sigma_i\sigma_j\right)^{1/2}

    .. math::
        \sigma_m^3 = \sum_i \sum_j x_ix_j\sigma_{ij}^3

    .. math::
        \frac{\epsilon_i}{k} = \frac{T_{ci}}{1.2593}

    .. math::
        \frac{\epsilon_{ij}}{k} = \zeta\left(\frac{\epsilon_i}{k}
        \frac{\epsilon_j}{k}\right)^{1/2}

    .. math::
        \left(\frac{\epsilon}{k}\right)_m = \frac{\sum_i \sum_j x_ix_j
        \left(\epsilon_{ij}/k\right)\sigma_{ij}^3}{\sigma_m^3}

    .. math::
        \omega_{ij} = \frac{\omega_i+\omega_j}{2}

    .. math::
        \kappa_{ij} = \left(\kappa_i+\kappa_j\right)^{1/2}

    .. math::
        M_{ij} = \frac{2M_iM_j}{M_i+M_j}

    .. math::
        M_m = \left[\frac{\sum_i\sum_jx_ix_j\left(\epsilon_{ij}/k\right)\sigma_
        {ij}^3M_{ij}^{1/2}}{\left(\epsilon/k\right)_m\sigma_m^3}\right]^2

    .. math::
        \omega_m =\frac{\sum_i\sum_jx_ix_j\omega_{ij}\sigma_{ij}^3}{\sigma_m^3}

    .. math::
        \mu_m^4 = \sigma_m^3\sum_i\sum_j\frac{x_ix_j\mu_i^2\mu_j^2}
        {\sigma_{ij}^3}

    .. math::
        \kappa_m = \sum_i \sum_j x_ix_j\kappa_{ij}

    .. math::
        F_{cm} = 1-0.2756\omega_m+0.059035\mu_{rm}^4+\kappa_m

    .. math::
        mu_{rm} = \frac{131.3\mu}{\left(V_{cm}T_{cm}\right)^{1/2}}

    .. math::
        T_{cm} = 1.2593\left(\frac{\epsilon}{k}\right)_m

    .. math::
        V_{cm} = \left(\frac{\sigma_m}{0.809}\right)^3

    Parameters
    ----------
    T : float
        Temperature, [K]
    xi : list
        Mole fractions of components, [-]
    Tci : list
        Critical temperature of components, [K]
    Vci : list
        Critical volume of components, [m³/kg]
    Mi : list
        Molecular weights of components, [g/mol]
    wi : list
        Acentric factor of components, [-]
    Di : list
        Dipole moment of components, [Debye]
    ki : list
        Correction factor for polar substances, [-]

    Returns
    -------
    mu : float
        Viscosity of gas, [Pa·s]

    Notes
    -----
    The method use binary interaction parameters for disimilar molecules, polar
    of hidrogen-bonding substances. That parameters must be calculated from
    experimental data. In this implementation this paramter set equal to unity

    Examples
    --------
    Example 9-8 in [3]_, 20.4% H2S 79.6% ethyl ether at 331K

    >>> x = [0.204, 0.796]
    >>> Tc = [373.4, 466.7]
    >>> Vc = [98/34.082/1000, 280/74.123/1000]
    >>> M = [34.082, 74.123]
    >>> w = [0.09, 0.281]
    >>> mu = [0.9, 1.3]
    >>> k = [0, 0]
    >>> "%0.1f" % MuG_Chung(331, x, Tc, Vc, M, w, mu, k).microP
    '87.6'
    """
    # Use critical volume in molar base
    Vci = [Vc*M*1000 for Vc, M in zip(Vci, Mi)]

    sigmai = [0.809*Vc**(1/3) for Vc in Vci]                            # Eq 4
    eki = [Tc/1.2593 for Tc in Tci]                                     # Eq 5

    # Eq 23
    sigmaij = []
    for s_i in sigmai:
        sigmaiji = []
        for s_j in sigmai:
            sigmaiji.append((s_i*s_j)**0.5)
        sigmaij.append(sigmaiji)

    # Eq 24
    ekij = []
    for ek_i in eki:
        ekiji = []
        for ek_j in eki:
            ekiji.append((ek_i*ek_j)**0.5)
        ekij.append(ekiji)

    # Eq 25
    wij = []
    for w_i in wi:
        wiji = []
        for w_j in wi:
            wiji.append((w_i+w_j)/2)
        wij.append(wiji)

    # Eq 26
    Mij = []
    for M_i in Mi:
        Miji = []
        for M_j in Mi:
            Miji.append(2*M_i*M_j/(M_i+M_j))
        Mij.append(Miji)

    # Eq 27
    kij = []
    for k_i in ki:
        kiji = []
        for k_j in ki:
            kiji.append((k_i*k_j)**0.5)
        kij.append(kiji)

    # Eq 14
    sm = 0
    for x_i, sigmaiji in zip(xi, sigmaij):
        for x_j, sij in zip(xi, sigmaiji):
            sm += x_i*x_j*sij**3
    sm = sm**(1/3)

    # Eq 15
    ekm = 0
    for x_i, sigmaiji, ekiji in zip(xi, sigmaij, ekij):
        for x_j, sij, ek in zip(xi, sigmaiji, ekiji):
            ekm += x_i*x_j*ek*sij**3
    ekm /= sm**3

    # Eq 18
    wm = 0
    for x_i, sigmaiji, wiji in zip(xi, sigmaij, wij):
        for x_j, sij, w in zip(xi, sigmaiji, wiji):
            wm += x_i*x_j*w*sij**3
    wm /= sm**3

    # Eq 19
    Mm = 0
    for x_i, sigmaiji, ekiji, Miji in zip(xi, sigmaij, ekij, Mij):
        for x_j, s, ek, M in zip(xi, sigmaiji, ekiji, Miji):
            Mm += x_i*x_j*ek*s**2*M**0.5
    Mm = (Mm/(ekm*sm**2))**2

    # Eq 20
    Dm = 0
    for x_i, sigmaiji, ekiji, D_i in zip(xi, sigmaij, ekij, Di):
        for x_j, s, ek, D_j in zip(xi, sigmaiji, ekiji, Di):
            Dm += x_i*x_j*(D_i*D_j)**2/ek/s**3
    Dm = (Dm*ekm*sm**3)**0.25

    # Eq 21
    km = 0
    for x_i, kiji in zip(xi, kij):
        for x_j, k in zip(xi, kiji):
            km += x_i*x_j*k

    Vcm = (sm/0.809)**3                                                # Eq 16
    Tcm = 1.2593*ekm                                                   # Eq 17
    murm = 131.3*Dm/(Vcm*Tcm)**0.5                                     # Eq 22

    T_ = T/ekm
    omega = Collision_Neufeld(T_)

    Fcm = 1 - 0.2756*wm + 0.059035*murm**4 + km                         # Eq 7
    mu = 40.785*Fcm*Mm**0.5*T**0.5/Vcm**(2/3)/omega                     # Eq 6
    return unidades.Viscosity(mu, "microP")


@refDoc(__doi__, [15, 1])
def MuG_P_Chung(T, xi, Tci, Vci, Mi, wi, Di, ki, rho, muo):
    r"""Calculate the viscosity of a compressed gas using the Chung correlation

    .. math::
        \mu=40.785\frac{F_c\left(MT\right)^{1/2}}{V_c^{2/3}\Omega_v}

    Parameters
    ----------
    T : float
        Temperature, [K]
    xi : list
        Mole fractions of components, [-]
    Tci : list
        Critical temperature of components, [K]
    Vci : list
        Critical volume of components, [m³/kg]
    Mi : list
        Molecular weights of components, [g/mol]
    wi : list
        Acentric factor of components, [-]
    Di : list
        Dipole moment of components, [Debye]
    ki : list
        Correction factor for polar substances, [-]
    rho : float
        Density, [kg/m³]
    muo : float
        Viscosity of low-pressure gas, [Pa·s]

    Returns
    -------
    mu : float
        Viscosity of gas mixture, [Pa·s]
    """
    # Use critical volume in molar base
    Vci = [Vc*M*1000 for Vc, M in zip(Vci, Mi)]

    sigmai = [0.809*Vc**(1/3) for Vc in Vci]                            # Eq 4
    eki = [Tc/1.2593 for Tc in Tci]                                     # Eq 5

    # Eq 23
    sigmaij = []
    for s_i in sigmai:
        sigmaiji = []
        for s_j in sigmai:
            sigmaiji.append((s_i*s_j)**0.5)
        sigmaij.append(sigmaiji)

    # Eq 24
    ekij = []
    for ek_i in eki:
        ekiji = []
        for ek_j in eki:
            ekiji.append((ek_i*ek_j)**0.5)
        ekij.append(ekiji)

    # Eq 25
    wij = []
    for w_i in wi:
        wiji = []
        for w_j in wi:
            wiji.append((w_i+w_j)/2)
        wij.append(wiji)

    # Eq 26
    Mij = []
    for M_i in Mi:
        Miji = []
        for M_j in Mi:
            Miji.append(2*M_i*M_j/(M_i+M_j))
        Mij.append(Miji)

    # Eq 27
    kij = []
    for k_i in ki:
        kiji = []
        for k_j in ki:
            kiji.append((k_i*k_j)**0.5)
        kij.append(kiji)

    # Eq 14
    sm = 0
    for x_i, sigmaiji in zip(xi, sigmaij):
        for x_j, sij in zip(xi, sigmaiji):
            sm += x_i*x_j*sij**3
    sm = sm**(1/3)

    # Eq 15
    ekm = 0
    for x_i, sigmaiji, ekiji in zip(xi, sigmaij, ekij):
        for x_j, sij, ek in zip(xi, sigmaiji, ekiji):
            ekm += x_i*x_j*ek*sij**3
    ekm /= sm**3

    # Eq 18
    wm = 0
    for x_i, sigmaiji, wiji in zip(xi, sigmaij, wij):
        for x_j, sij, w in zip(xi, sigmaiji, wiji):
            wm += x_i*x_j*w*sij**3
    wm /= sm**3

    # Eq 19
    Mm = 0
    for x_i, sigmaiji, ekiji, Miji in zip(xi, sigmaij, ekij, Mij):
        for x_j, s, ek, M in zip(xi, sigmaiji, ekiji, Miji):
            Mm += x_i*x_j*ek*s**2*M**0.5
    Mm = (Mm/(ekm*sm**2))**2

    # Eq 20
    Dm = 0
    for x_i, sigmaiji, ekiji, D_i in zip(xi, sigmaij, ekij, Di):
        for x_j, s, ek, D_j in zip(xi, sigmaiji, ekiji, Di):
            Dm += x_i*x_j*(D_i*D_j)**2/ek/s**3
    Dm = (Dm*ekm*sm**3)**0.25

    # Eq 21
    km = 0
    for x_i, kiji in zip(xi, kij):
        for x_j, k in zip(xi, kiji):
            km += x_i*x_j*k

    rho = rho/Mm/1000
    Vcm = (sm/0.809)**3                                                # Eq 16
    Tcm = 1.2593*ekm                                                   # Eq 17
    murm = 131.3*Dm/(Vcm*Tcm)**0.5                                     # Eq 22

    T_ = T/ekm

    # Table II
    dat = [
        (6.32402, 50.4119, -51.6801, 1189.02),
        (0.12102e-2, -0.11536e-2, -0.62571e-2, 0.37283e-1),
        (5.28346, 254.209, -168.481, 3898.27),
        (6.62263, 38.09570, -8.46414, 31.4178),
        (19.74540, 7.63034, -14.35440, 31.5267),
        (-1.89992, -12.53670, 4.98529, -18.1507),
        (24.27450, 3.44945, -11.29130, 69.3466),
        (0.79716, 1.11764, 0.12348e-1, -4.11661),
        (-0.23816, 0.67695e-1, -0.81630, 4.02528),
        (0.68629e-1, 0.34793, 0.59256, -0.72663)]

    # Eq 11
    A = []
    for ao, a1, a2, a3 in dat:
        A.append(ao + a1*wm + a2*murm**4 + a3*km)
    A1, A2, A3, A4, A5, A6, A7, A8, A9, A10 = A

    Y = rho*Vcm/6
    G1 = (1-0.5*Y)/(1-Y)**3
    G2 = (A1*((1-exp(-A4*Y))/Y)+A2*G1*exp(A5*Y)+A3*G1)/(A1*A4+A2+A3)

    muk = muo*(1/G2 + A6*Y)
    mup = (36.344e-6*(Mm*Tcm)**0.5/Vcm**(2/3))*A7*Y**2*G2*exp(
        A8+A9/T_+A10/T_**2)

    return unidades.Viscosity(muk+mup, "P")


@refDoc(__doi__, [1, 21, 16, 17])
def MuG_TRAPP(T, P, xi, Tci, Vci, Zci, Mi, wi, rho, muo):
    r"""Calculate the viscosity of a compressed gas using the TRAPP (TRAnsport
    Property Prediction) method.

    .. math::
        \eta_m - \eta_m^o = F_{\eta m}\left(\eta^R-\eta^{Ro}\right) +
        \Delta\eta^{ENSKOG}

    .. math::
        h_m = \sum_i\sum_jx_ix_jh_{ij}

    .. math::
        f_mh_m = \sum_i\sum_jx_ix_jf_{ij}h_{ij}

    .. math::
        h_{ij} = \frac{\left(h_i^{1/3}+h_j^{1/3}\right)^3}{8}

    .. math::
        f_{ij} = \left(f_if_j\right)^{1/2}

    .. math::
        T_o = T/f_m

    .. math::
        \rho_o = \rho h_m

    .. math::
        F_{\eta m} = \frac{1}{44.094^{1/2}h_m^2} \sum_i\sum_j x_ix_j
        \left(f_{ij}M_{ij}\right)^{1/2}h_{ij}^{4/3}

    .. math::
        M_{ij} = \frac{2M_iM_j}{M_i+M_j}

    .. math::
        \Delta\eta^{ENSKOG} = \eta_m^{ENSKOG} - \eta_x^{ENSKOG}

    .. math::
        \eta_m^{ENSKOG} = \sum_i \beta_iY_i + \sum_i\sum_jx_ix_j\sigma_{ij}^6
        \eta_{ij}^og_{ij}

    .. math::
        \sigma_i = 4.771h_i^{1/3}

    .. math::
        \sigma_{ij} = \frac{\sigma_i+\sigma_j}{2}

    .. math::
        g_{ij} = \frac{1}{1-\xi}+\frac{3\xi}{\left(1-\xi\right)^2}\Theta_{ij}+
        \frac{2\xi^2}{\left(1-\xi\right)^3}\Theta_{ij}^2

    .. math::
        \Theta_{ij} = \frac{\sigma_i\sigma_j}{2\sigma_{ij}}
        \frac{\sum_kx_k\sigma_k^2}{\sum_k x_k\sigma_k^3}

    .. math::
        \xi = 6.023e-4\frac{\pi}{6}\rho\sum_i x_i\sigma_i^3

    .. math::
        Y_i = x_i\left[1+\frac{8\pi}{15}6.023e-4\rho\sum_jx_j\left(\frac{M_j}
        {M_i+M_j}\right)\sigma_{ij}^3g_{ij}\right]

    .. math::
        \sum_i B_{ij}\beta_j = Y_i

    .. math::
        B_{ij} = 2\sum_kx_ix_k\frac{g_{ik}}{\eta_{ij}^o}\left(\frac{M_k}
        {M_i+M_k}\right)^2\left[\left(1+\frac{5}{3}\frac{M_i}{M_k}\right)
        \delta_{ij}-\frac{2}{3}\frac{M_i}{M_k}\delta_{jk}\right]

    .. math::
        \sigma_x = \left(\sum_i\sum_jx_ix_j\sigma_{ij}^3\right)^{1/3}

    .. math::
        M_x = \frac{\left(\sum_i\sum_jx_ix_jM_{ij}^{1/2}\sigma_{ij}^4\right)^2}
        {\sigma_x^8}

    Parameters
    ----------
    T : float
        Temperature, [K]
    xi : list
        Mole fractions of components, [-]
    Tci : list
        Critical temperature of components, [K]
    Vci : list
        Critical volume of components, [m³/kg]
    Zci : list
        Critical compressibility factor of components, [K]
    Mi : list
        Molecular weights of components, [g/mol]
    wi : list
        Acentric factor of components, [-]
    rho : float
        Density, [kg/m³]
    muo : float
        Viscosity of low-pressure gas, [Pa·s]

    Returns
    -------
    mu : float
        Viscosity of gas, [Pa·s]

    Examples
    --------
    Example 9-14 in [3]_, 80% methane, 20% nC10 at 377.6K and 413.7bar

    >>> x = [0.8, 0.2]
    >>> M = [16.043, 142.285]
    >>> Tc = [190.56, 617.7]
    >>> Vc = [98.6/16.043/1000, 624/142.285/1000]
    >>> Zc = [0.286, 0.256]
    >>> w = [0.011, 0.49]
    >>> rho = 1/243.8*58.123*1000
    >>> mu = MuG_TRAPP(377.6, 413.7e5, x, Tc, Vc, Zc, M, w, 448.4, 10.2e-6)
    >>> "%0.1f" % mu.muPas
    '82.6'
    """
    # Reference fluid properties, propane
    TcR = 369.83
    rhocR = 1/200  # mol/cm³
    ZcR = 0.276
    wR = 0.152

    # Convert volume to molar base
    Vci = [Vc*M*1000 for Vc, M in zip(Vci, Mi)]
    Mm = sum([x*M for x, M in zip(xi, Mi)])
    rho = rho/Mm

    # Calculate shape factor for mixture
    fi = []
    hi = []
    for Tc, Vc, Zc, w in zip(Tci, Vci, Zci, wi):
        fi.append(Tc/TcR*(1+(w-wR)*(0.05203-0.7498*log(T/Tc))))
        hi.append(rhocR*Vc*ZcR/Zc*(1-(w-wR)*(0.1436-0.2822*log(T/Tc))))

    # Eq 9-7.5
    fij = []
    for f_i in fi:
        fiji = []
        for f_j in fi:
            fiji.append((f_i*f_j)**0.5)
        fij.append(fiji)

    # Eq 9-7.4
    hij = []
    for h_i in hi:
        hiji = []
        for h_j in hi:
            hiji.append((h_i**(1/3)+h_j**(1/3))**3/8)
        hij.append(hiji)

    # Eq 9-7.2
    hm = 0
    for x_i, hiji in zip(xi, hij):
        for x_j, h in zip(xi, hiji):
            hm += x_i*x_j*h

    # Eq 9-7.3
    fm = 0
    for x_i, hiji, fiji in zip(xi, hij, fij):
        for x_j, h, f in zip(xi, hiji, fiji):
            fm += x_i*x_j*f*h/hm

    # Eq 9-7.9
    Mij = []
    for M_i in Mi:
        Miji = []
        for M_j in Mi:
            Miji.append(2*M_i*M_j/(M_i+M_j))
        Mij.append(Miji)

    To = T/fm                                                       # Eq 9-7.6
    rho0 = rho/1000*hm                                              # Eq 9-7.7

    # Eq 9-7.8
    suma = 0
    for x_i, fiji, Miji, hiji in zip(xi, fij, Mij, hij):
        for x_j, f, M, h in zip(xi, fiji, Miji, hiji):
            suma += x_i*x_j*(f*M)**0.5*h**(4/3)
    Fnm = 44.094**-0.5/hm**2*suma

    # Calculation of reference residual viscosity
    # Coefficients in [16]_, pag 796
    # Density are in mol/dm³
    rho0 *= 1000
    rhocR *= 1000
    G = -14.113294896 + 968.22940153/To
    H = rho0**0.5*(rho0-rhocR)/rhocR
    G2 = 13.686545032 - 12511.628378/To**1.5
    G3 = 0.0168910864 + 43.527109444/To + 7659.4543472/To**2
    F = G + G2*rho0**0.1+G3*H
    muR = exp(F)-exp(G)

    # Calculate of Δη
    sigmai = [4.771*h**(1/3) for h in hi]
    sigmaij = []
    for s_i in sigmai:
        sigmaiji = []
        for s_j in sigmai:
            sigmaiji.append((s_i+s_j)/2)
        sigmaij.append(sigmaiji)

    # Eq 9-7.16
    X = 6.023e-4*pi/6*rho*sum([x*s**3 for x, s in zip(xi, sigmai)])

    # Eq 9-7.15
    sum2 = sum([x*s**2 for x, s in zip(xi, sigmai)])
    sum3 = sum([x*s**3 for x, s in zip(xi, sigmai)])
    titaij = []
    for s_i, siji in zip(sigmai, sigmaij):
        titaiji = []
        for s_j, sij in zip(sigmai, siji):
            titaiji.append(s_i*s_j/2/sij*sum2/sum3)
        titaij.append(titaiji)

    # Eq 9-7.14
    gij = []
    for titaiji in titaij:
        giji = []
        for tij in titaiji:
            giji.append(1/(1-X)+3*X/(1-X)**2*tij+2*X**2/(1-X)**3*tij**2)
        gij.append(giji)

    # Eq 9-3.8
    muij = []
    for miji, siji in zip(Mij, sigmaij):
        muiji = []
        for m, s in zip(miji, siji):
            muiji.append(2.669*(m*T)**0.5/s**2)
        muij.append(muiji)

    # Eq 9-7.19
    Bij = []
    for i, M_i in enumerate(Mi):
        Biji = []
        for j, M_j in enumerate(Mi):
            B = 0
            for k, M_k in enumerate(Mi):
                if i == j:
                    dij = 1
                else:
                    dij = 0
                if j == k:
                    djk = 1
                else:
                    djk = 0
                B += xi[i]*xi[k]*gij[i][k]/muij[i][k]*(M_k/(M_i+M_k))**2*(
                    (1+5/3*M_i/M_k)*dij - 2/3*M_i/M_k*djk)
            B *= 2e-1
            Biji.append(B)
        Bij.append(Biji)

    # Eq 9-7.17
    Yi = []
    for x_i, M_i, siji, giji in zip(xi, Mi, sigmaij, gij):
        suma = 0
        for x_j, M_j, s, g in zip(xi, Mi, siji, giji):
            suma += x_j*M_j/(M_i+M_j)*s**3*g
        Yi.append(x_i*(1+8*pi/15*6.023e-4*rho*suma))

    # Eq 9-7.18
    betai = solve(Bij, Yi)

    # Eq 9-7.11
    sum1 = sum([b*Y for b, Y in zip(betai, Yi)])
    sum2 = 0
    for x_i, siji, muiji, giji in zip(xi, sigmaij, muij, gij):
        for x_j, s, mu, g in zip(xi, siji, muiji, giji):
            sum2 += x_i*x_j*s**6*mu*g
    alfa = 48/25/pi*(2*pi/3*6.023e-4)**2
    eta_m = sum1 + alfa*10*rho**2*sum2

    # ηx for pure hypothethical fluid
    # Eq 9-7.20
    sigmax = 0
    for x_i, siji in zip(xi, sigmaij):
        for x_j, s in zip(xi, siji):
            sigmax += x_i*x_j*s**3
    sigmax = sigmax**(1/3)

    # Eq 9-7.21
    Mx = 0
    for x_i, siji, Miji in zip(xi, sigmaij, Mij):
        for x_j, s, M in zip(xi, siji, Miji):
            Mx += x_i*x_j*M**0.5*s**4
    Mx = Mx**2/sigmax**8

    X = 6.023e-4*pi/6*rho*sigmax**3
    gxx = 1/(1-X)+3*X/(1-X)**2*0.5+2*X**2/(1-X)**3*0.25
    Yx = 1+8*pi/15*6.023e-4*rho*sigmax**3*gxx/2
    mux = 26.69*(Mx*T)**0.5/sigmax**2
    Bxx = gxx/mux
    betax = Yx/Bxx
    eta_x = betax*Yx+alfa*rho**2*sigmax**6*mux*gxx

    # Eq 9-7.10
    # The 0.1 factor because this values are im μP, to convert to μPa·s
    Dmu = (eta_m-eta_x)*0.1
    return unidades.Viscosity(Fnm*muR + Dmu + muo*1e6, "muPas")


@refDoc(__doi__, [19, 2])
def MuG_DeanStielMix(xi, Tci, Pci, Mi, rhoc, rho, muo):
    r"""Calculate the viscosity of a compressed gas using the Dean-Stiel
    correlation, also referenced in API databook Procedure 11B4.1, pag 1107

    .. math::
        \left(\mu-\mu_o\right)\xi=10.8x10^{-5}\left[exp\left(1.439\rho_r\right)
        -exp\left(-1.11\rho_r^{1.858}\right)\right]

    .. math::
        \xi = \frac{T_{pc}^{1/6}}{M_m^{1/2}P_{pc}^{2/3}

    .. math::
        T_{pc} = \sum_i x_iT_{ci}

    .. math::
        P_{pc} = \sum_i x_iP_{ci}

    .. math::
        M_m = \sum_i x_iM_i

    Parameters
    ----------
    xi : list
        Mole fractions of components, [-]
    Tci : float
        Critical temperature of components, [K]
    Pci : float
        Critical pressure of components, [Pa]
    Mi : list
        Molecular weight of components, [g/mol]
    rhoc : float
        Critical Density of mixture, [kg/m³]
    rho : float
        Density, [kg/m³]
    muo : float
        Viscosity of low-pressure gas, [Pa·s]

    Returns
    -------
    mu : float
        Viscosity of gas, [Pa·s]

    Examples
    --------
    Example in [2]_, 60% Methane 40% propane at 1500psi and 257ºF

    >>> Tc1 = unidades.Temperature(-116.66, "F")
    >>> Tc2 = unidades.Temperature(206.02, "F")
    >>> Pc1 = unidades.Pressure(667.04, "psi")
    >>> Pc2 = unidades.Pressure(616.13, "psi")
    >>> x = [0.6, 0.4]
    >>> args = (x, [Tc1, Tc2], [Pc1, Pc2], [16.04, 44.1], 1, 0.5283, 123e-7)
    >>> "%0.4f" % MuG_DeanStielMix(*args).cP
    '0.0163'
    """
    # Calculate pseudo critical properties
    Tpc = sum([x*Tc for x, Tc in zip(xi, Tci)])                         # Eq 5
    Ppc = sum([x*Pc for x, Pc in zip(xi, Pci)])                         # Eq 6
    M = sum([x*M_i for x, M_i in zip(xi, Mi)])

    return MuG_DeanStiel(Tpc, Ppc, rhoc, M, rho, muo)


@refDoc(__doi__, [2])
def MuG_APIMix(T, P, xi, Tci, Pci, muo):
    r"""Calculate the viscosity of nonhydrocarbon gases at high pressure using
    the linearization of Carr figure as give in API Databook procedure 11C1.2,
    pag 1113

    .. math::
        \frac{\mu}{\mu_o}=A_1hP_r^f + A_2\left(kP_r^l+mP_r^n+pP_r^q\right)

    Parameters
    ----------
    T : float
        Temperature, [K]
    P : float
        Pressure, [Pa]
    xi : list
        Mole fractions of components, [-]
    Tci : float
        Critical temperature of components, [K]
    Pci : float
        Critical pressure of components, [Pa]
    muo : float
        Viscosity of low-pressure gas, [Pa·s]

    Returns
    -------
    mu : float
        Viscosity of gas, [Pa·s]

    Notes
    -----
    This method is recomended for gaseous nonhydrocarbons at high pressure,
    although this method is also applicable for hydrocarbons.

    Examples
    --------
    Example B in [2]_, 60% Methane 40% propane at 1500psi and 257ºF

    >>> T = unidades.Temperature(257, "F")
    >>> P = unidades.Pressure(1500, "psi")
    >>> Tc1 = unidades.Temperature(-116.66, "F")
    >>> Tc2 = unidades.Temperature(206.02, "F")
    >>> Pc1 = unidades.Pressure(667.04, "psi")
    >>> Pc2 = unidades.Pressure(616.13, "psi")
    >>> x = [0.6, 0.4]
    >>> "%0.4f" % MuG_APIMix(T, P, x, [Tc1, Tc2], [Pc1, Pc2], 123e-7).cP
    '0.0173'
    """
    # Calculate pseudo critical properties
    Tpc = sum([x*Tc for x, Tc in zip(xi, Tci)])                         # Eq 5
    Ppc = sum([x*Pc for x, Pc in zip(xi, Pci)])                         # Eq 6

    return MuG_API(T, P, Tpc, Ppc, muo)


# Liquid thermal conductivity correlations
@refDoc(__doi__, [4, 2])
def ThL_Li(xi, Vi, Mi, ki):
    r"""Calculate thermal conductiviy of liquid nmixtures using the Li method,
    also referenced in API procedure 12A2.1, pag 1145

    .. math::
        \lambda_{m} = \sum_{i} \sum_{j} \phi_i \phi_j k_{ij}

    .. math::
        k_{ij} = 2\left(\lambda_i^{-1}+\lambda_j^{-1}\right)^{-1}

    .. math::
        \phi_{i} = \frac{x_iV_i}{\sum_j x_jV_j}

    Parameters
    ----------
    xi : list
        Mole fractions of components, [-]
    Vi : list
        Specific volume of components, [m³/kg]
    Mi : list
        Molecular weight of components, [g/mol]
    ki : list
        Thermal conductivities of components, [W/m·K]

    Returns
    -------
    k : float
        Thermal conductivities of mixture, [W/m·K]

    Examples
    --------
    Example from [2]_ 68% nC7, 32% CycloC5 at 32F and 1atm

    >>> V1 = unidades.MolarVolume(2.285, "ft3lbmol")
    >>> V2 = unidades.MolarVolume(1.473, "ft3lbmol")
    >>> k1 = unidades.ThermalConductivity(0.07639, "BtuhftF")
    >>> k2 = unidades.ThermalConductivity(0.08130, "BtuhftF")
    >>> "%0.5f" % ThL_Li([0.68, 0.32], [V1, V2], [1, 1], [k1, k2]).BtuhftF
    '0.07751'
    """
    # Use critical volume in molar base
    Vi = [V*M*1000 for V, M in zip(Vi, Mi)]

    # Calculation of binary thermal conductivity pair, Eq 2
    kij = []
    for k_i in ki:
        kiji = []
        for k_j in ki:
            kiji.append(2/(1/k_i+1/k_j))
        kij.append(kiji)

    # Calculation of volume fraction, Eq 3
    suma = 0
    for x, V in zip(xi, Vi):
        suma += x*V
    phi = []
    for x, V in zip(xi, Vi):
        phi.append(x*V/suma)

    # Calculation of misture thermal conductivity, Eq 1
    k = 0
    for phi_i, ki in zip(phi, kij):
        for phi_j, k_ij in zip(phi, ki):
            k += phi_i*phi_j*k_ij

    return unidades.ThermalConductivity(k)


@refDoc(__doi__, [3])
def ThL_Power(wi, ki):
    r"""Calculate thermal conductiviy of liquid nmixtures using the power law
    method, as referenced in [3]_

    .. math::
        \lambda_{m} = \frac{1}{\sqrt{\sum_{i} w_i/\lambda_i^2}}

    Parameters
    ----------
    wi : list
        Weight fractions of components, [-]
    ki : list
        Thermal conductivities of components, [W/m·K]

    Returns
    -------
    k : float
        Thermal conductivities of mixture, [W/m·K]

    Notes
    -----
    This method shold not be used if water is in the mixture or if pure
    component thermal conductivities are very different, ki/kj < 2
    """
    km = 0
    for w, k in zip(wi, ki):
        km += w/k**2

    km = km**-0.5
    return unidades.ThermalConductivity(km)


# Gas thermal conductivity correlations
@refDoc(__doi__, [5, 2])
def ThG_LindsayBromley(T, xi, Mi, Tbi, mui, ki):
    r"""Calculate thermal conductiviy of gas mixtures using the Lindsay-Bromley
    method, also referenced in API procedure 12B2.1, pag 1164

    .. math::
        k = \sum_{i} \frac{k_i}{\frac{1}{x_i}\sum x_i A_{ij}}

    .. math::
        A_{ij} = \frac{1}{4} \left\{ 1 + \left[\frac{\mu_i}{\mu_j}
        \left(\frac{M_j}{M_i}\right)^{0.75} \left(\frac{1+S_i/T}{1+S_j/T}
        \right)\right]^{0.5}\right\}^2\left(\frac{1+S_{ij}/T}{1+S_i/T}\right)

    .. math::
        S_{ij} = (S_i S_j)^{0.5}

    Parameters
    ----------
    T : float
        Temperature, [K]
    xi : list
        Mole fractions of components, [-]
    Mi : float
        Molecular weights of components, [g/mol]
    Tbi : float
        Boiling points of components, [K]
    mui : float
        Gas viscosities of components, [Pa·s]
    ki : list
        Thermal conductivities of components, [W/m·K]

    Returns
    -------
    k : float
        Thermal conductivities of mixture, [W/m·K]

    Examples
    --------
    Example from [2]_ 29.96% nC5, 70.04% nC6 at 212F and 1atm

    >>> T = unidades.Temperature(212, "F")
    >>> x = [0.2996, 0.7004]
    >>> M = [72.15, 86.18]
    >>> Tb1 = unidades.Temperature(96.93, "F")
    >>> Tb2 = unidades.Temperature(155.71, "F")
    >>> mu1 = unidades.Viscosity(0.008631, "cP")
    >>> mu2 = unidades.Viscosity(0.008129, "cP")
    >>> k1 = unidades.ThermalConductivity(0.01280, "BtuhftF")
    >>> k2 = unidades.ThermalConductivity(0.01165, "BtuhftF")
    >>> k = ThG_LindsayBromley(T, x, M, [Tb1, Tb2], [mu1, mu2], [k1, k2])
    >>> "%0.5f" % k.BtuhftF
    '0.01197'
    """
    # Calculation of Sutherland constants, Eq 14
    S = []
    for Tb, M in zip(Tbi, Mi):
        if M == 2.0158 or M == 4.0026:
            # Hydrogen or helium case
            S.append(79)
        else:
            S.append(1.5*Tb)

    # Geometric mean of collision Sutherland constants, Eq 15
    Sij = []
    for Si in S:
        Siji = []
        for Sj in S:
            Siji.append((Si*Sj)**0.5)
        Sij.append(Siji)

    # Eq 12
    Aij = []
    for mu_i, M_i, Si, Siji in zip(mui, Mi, S, Sij):
        Aiji = []
        for mu_j, M_j, Sj, S_ij in zip(mui, Mi, S, Siji):
            Aiji.append(0.25*(1+(mu_i/mu_j*(M_j/M_i)**0.75*(1+Si/T)/(
                1+Sj/T))**0.5)**2 * (1+S_ij/T)/(1+Si/T))
        Aij.append(Aiji)

    # Calculate thermal conductivity, Eq 11
    sumaj = []
    for Aiji in Aij:
        suma = 0
        for xj, A_ij in zip(xi, Aiji):
            suma += A_ij*xj
        sumaj.append(suma)
    k = 0
    for x_i, k_i, si in zip(xi, ki, sumaj):
        k += k_i*x_i/si
    return unidades.ThermalConductivity(k)


@refDoc(__doi__, [6, 3])
def ThG_MasonSaxena(xi, Mi, mui, ki):
    r"""Calculate thermal conductiviy of gas mixtures using the Mason-Saxena
    method

    .. math::
        k = \sum_{i} \frac{k_i}{\frac{1}{x_i}\sum x_i A_{ij}}

    .. math::
        A_{ij} = \frac{\epsilon \left[1+\left(\lambda_{tri}/\lambda_{trj}
        \right)^{1/2} \left(M_i/M_j\right)^{1/4}\right]^2}
        {\left[8\left(1+M_i/M_j\right)\right]^{1/2}}

    .. math::
        \frac{\lambda_{tri}}{\lambda_{trj}}=\frac{\mu_i}{\mu_j}\frac{M_i}{M_j}

    Parameters
    ----------
    xi : list
        Mole fractions of components, [-]
    Mi : float
        Molecular weights of components, [g/mol]
    mui : float
        Gas viscosities of components, [Pa·s]
    ki : list
        Thermal conductivities of components, [W/m·K]

    Returns
    -------
    k : float
        Thermal conductivities of mixture, [W/m·K]

    Examples
    --------
    Example 10-5 from [3]_; 25% benzene, 75% Argon at 100.6ºC and 1bar

    >>> xi = [0.25, 0.75]
    >>> Mi = [78.114, 39.948]
    >>> mui = [92.5e-7, 271e-7]
    >>> ki = [0.0166, 0.0214]
    >>> "%0.4f" % ThG_MasonSaxena(xi, Mi, mui, ki)
    '0.0184'
    """
    # Aij coefficient with ε=1 as explain in [3]_, Eq 21
    Aij = []
    for mu_i, M_i in zip(mui, Mi):
        Aiji = []
        for mu_j, M_j in zip(mui, Mi):
            # Monatomic value of thermal conductivity ratio, Eq 22
            lt_ij = mu_i*M_j/mu_j/M_i
            Aiji.append((1+lt_ij**0.5*(M_i/M_j)**0.25)**2/(8*(1+M_i/M_j))**0.5)
        Aij.append(Aiji)

    # Calculate thermal conductivity, Eq 20
    sumaj = []
    for Aiji in Aij:
        suma = 0
        for xj, A_ij in zip(xi, Aiji):
            suma += A_ij*xj
        sumaj.append(suma)
    k = 0
    for x_i, k_i, si in zip(xi, ki, sumaj):
        k += k_i*x_i/si
    return unidades.ThermalConductivity(k)


@refDoc(__doi__, [15, 1])
def ThG_Chung(T, xi, Tci, Vci, Mi, wi, Cvi, mu):
    r"""Calculate thermal conductivity of gas mixtures at low pressure using
    the Chung correlation

    .. math::
        \lambda_o = \frac{7.452\mu_o\Psi}{M}

    .. math::
        \Psi = 1 + \alpha \left\{[0.215+0.28288\alpha-1.061\beta+0.26665Z]/
        [0.6366+\beta Z + 1.061 \alpha \beta]\right\}

    .. math::
        \alpha = \frac{C_v}{R}-1.5

    .. math::
        \beta = 0.7862-0.7109\omega + 1.3168\omega^2

    .. math::
        Z=2+10.5T_r^2

    Parameters
    ----------
    T : float
        Temperature, [K]
    xi : list
        Mole fractions of components, [-]
    Tci : list
        Critical temperature of compounds, [K]
    Vci : list
        Critical volume of compounds, [m³/kg]
    Mi : list
        Molecular weight of components, [g/mol]
    wi :  list
        Acentric factor of compounds, [-]
    Cvi : list
        Ideal gas heat capacity at constant volume of components, [J/kg·K]
    mu : float
        Gas viscosity [Pa·s]

    Returns
    -------
    k : float
        Thermal conductivity, [W/m/k]

    Examples
    --------
    Example 10-5 from [3]_; 25% benzene, 75% Argon at 100.6ºC and 1bar

    >>> T = unidades.Temperature(100.6, "C")
    >>> xi = [0.25, 0.75]
    >>> Tci = [562.05, 150.86]
    >>> Vci = [256/78.114/1000, 74.57/39.948/1000]
    >>> Mi = [78.114, 39.948]
    >>> wi = [0.21, -0.002]
    >>> Cvi = [96.2/78.114*1000, 12.5/39.948*1000]
    >>> mui = [92.5e-7, 271e-7]
    >>> ki = [0, 0]
    >>> mu = MuG_Chung(T, xi, Tci, Vci, Mi, wi, mui, ki)
    >>> "%0.1f" % mu.microP
    '183.3'
    >>> "%0.4f" % ThG_Chung(T, xi, Tci, Vci, Mi, wi, Cvi, mu)
    '0.0222'
    """
    # Molar values
    Vci = [Vc*M*1000 for Vc, M in zip(Vci, Mi)]
    Cvi = [Cv*M/1000 for Cv, M in zip(Cvi, Mi)]
    Cvm = sum([x*Cv for x, Cv in zip(xi, Cvi)])

    sigmai = [0.809*Vc**(1/3) for Vc in Vci]                            # Eq 4
    eki = [Tc/1.2593 for Tc in Tci]                                     # Eq 5

    # Eq 23
    sigmaij = []
    for s_i in sigmai:
        sigmaiji = []
        for s_j in sigmai:
            sigmaiji.append((s_i*s_j)**0.5)
        sigmaij.append(sigmaiji)

    # Eq 24
    ekij = []
    for ek_i in eki:
        ekiji = []
        for ek_j in eki:
            ekiji.append((ek_i*ek_j)**0.5)
        ekij.append(ekiji)

    # Eq 25
    wij = []
    for w_i in wi:
        wiji = []
        for w_j in wi:
            wiji.append((w_i+w_j)/2)
        wij.append(wiji)

    # Eq 26
    Mij = []
    for M_i in Mi:
        Miji = []
        for M_j in Mi:
            Miji.append(2*M_i*M_j/(M_i+M_j))
        Mij.append(Miji)

    # Eq 14
    sm = 0
    for x_i, sigmaiji in zip(xi, sigmaij):
        for x_j, sij in zip(xi, sigmaiji):
            sm += x_i*x_j*sij**3
    sm = sm**(1/3)

    # Eq 15
    ekm = 0
    for x_i, sigmaiji, ekiji in zip(xi, sigmaij, ekij):
        for x_j, sij, ek in zip(xi, sigmaiji, ekiji):
            ekm += x_i*x_j*ek*sij**3
    ekm /= sm**3

    # Eq 18
    wm = 0
    for x_i, sigmaiji, wiji in zip(xi, sigmaij, wij):
        for x_j, sij, w in zip(xi, sigmaiji, wiji):
            wm += x_i*x_j*w*sij**3
    wm /= sm**3

    # Eq 19
    Mm = 0
    for x_i, sigmaiji, ekiji, Miji in zip(xi, sigmaij, ekij, Mij):
        for x_j, s, ek, M in zip(xi, sigmaiji, ekiji, Miji):
            Mm += x_i*x_j*ek*s**2*M**0.5
    Mm = (Mm/(ekm*sm**2))**2

    Tcm = 1.2593*ekm
    Trm = T/Tcm

    alpha = Cvm/R - 1.5
    beta = 0.7862 - 0.7109*wm + 1.3168*wm**2
    Z = 2 + 10.5*Trm**2
    phi = 1 + alpha*((0.215 + 0.28288*alpha - 1.061*beta + 0.26665*Z)/(
        0.6366 + beta*Z + 1.061*alpha*beta))

    # Eq 9
    k = 7.452*mu*10/Mm*phi   # Viscosity in P
    return unidades.ThermalConductivity(k, "calscmK")


@refDoc(__doi__, [15, 1])
def ThG_P_Chung(T, xi, Tci, Vci, Mi, wi, Di, ki, rho, ko):
    r"""Calculate the thermal conductivity of a compressed gas mixture using
    the Chung correlation

    .. math::
        \lambda = \frac{31.2 \eta^o\Psi}{M}(1/G_2+B_6y)+qB_7y^2T_r^{1/2}G_2

    Parameters
    ----------
    T : float
        Temperature, [K]
    xi : list
        Mole fractions of components, [-]
    Tci : list
        Critical temperature of compounds, [K]
    Vci : list
        Critical volume of compounds, [m³/kg]
    Mi : list
        Molecular weight of components, [g/mol]
    wi :  list
        Acentric factor of compounds, [-]
    Di : list
        Dipole moment of compounds, [Debye]
    ki : list
        Corection factor for polar substances, [-]
    rho : float
        Density, [kg/m³]
    ko : float
        Low-pressure gas thermal conductivity[Pa*S]

    Returns
    -------
    k : float
        High-pressure gas thermal conductivity [W/m/k]

    Examples
    --------
    Example 10-7 from [3]_; 75.5% methane, 24.5% CO2 at 370.8K and 174.8bar

    >>> Tc = [190.56, 304.12]
    >>> Vc = [98.6/16.043/1000, 94.07/44.01/1000]
    >>> M = [16.043, 44.01]
    >>> w = [0.011, 0.225]
    >>> x = [0.755, 0.245]
    >>> D = [0, 0]
    >>> k = [0, 0]
    >>> args = (370.8, x, Tc, Vc, M, w, D, k, 1/159*22.9*1000, 0.0377)
    >>> "%0.3f" % ThG_P_Chung(*args)
    '0.058'
    """
    # Thermal conductivity in procedure in cal/s·cm·K
    ko = unidades.ThermalConductivity(ko).calscmK

    # Use critical volume in molar base
    Vci = [Vc*M*1000 for Vc, M in zip(Vci, Mi)]

    sigmai = [0.809*Vc**(1/3) for Vc in Vci]                            # Eq 4
    eki = [Tc/1.2593 for Tc in Tci]                                     # Eq 5

    # Eq 23
    sigmaij = []
    for s_i in sigmai:
        sigmaiji = []
        for s_j in sigmai:
            sigmaiji.append((s_i*s_j)**0.5)
        sigmaij.append(sigmaiji)

    # Eq 24
    ekij = []
    for ek_i in eki:
        ekiji = []
        for ek_j in eki:
            ekiji.append((ek_i*ek_j)**0.5)
        ekij.append(ekiji)

    # Eq 25
    wij = []
    for w_i in wi:
        wiji = []
        for w_j in wi:
            wiji.append((w_i+w_j)/2)
        wij.append(wiji)

    # Eq 26
    Mij = []
    for M_i in Mi:
        Miji = []
        for M_j in Mi:
            Miji.append(2*M_i*M_j/(M_i+M_j))
        Mij.append(Miji)

    # Eq 27
    kij = []
    for k_i in ki:
        kiji = []
        for k_j in ki:
            kiji.append((k_i*k_j)**0.5)
        kij.append(kiji)

    # Eq 14
    sm = 0
    for x_i, sigmaiji in zip(xi, sigmaij):
        for x_j, sij in zip(xi, sigmaiji):
            sm += x_i*x_j*sij**3
    sm = sm**(1/3)

    # Eq 15
    ekm = 0
    for x_i, sigmaiji, ekiji in zip(xi, sigmaij, ekij):
        for x_j, sij, ek in zip(xi, sigmaiji, ekiji):
            ekm += x_i*x_j*ek*sij**3
    ekm /= sm**3

    # Eq 18
    wm = 0
    for x_i, sigmaiji, wiji in zip(xi, sigmaij, wij):
        for x_j, sij, w in zip(xi, sigmaiji, wiji):
            wm += x_i*x_j*w*sij**3
    wm /= sm**3

    # Eq 19
    Mm = 0
    for x_i, sigmaiji, ekiji, Miji in zip(xi, sigmaij, ekij, Mij):
        for x_j, s, ek, M in zip(xi, sigmaiji, ekiji, Miji):
            Mm += x_i*x_j*ek*s**2*M**0.5
    Mm = (Mm/(ekm*sm**2))**2

    # Eq 20
    Dm = 0
    for x_i, sigmaiji, ekiji, D_i in zip(xi, sigmaij, ekij, Di):
        for x_j, s, ek, D_j in zip(xi, sigmaiji, ekiji, Di):
            Dm += x_i*x_j*(D_i*D_j)**2/ek/s**3
    Dm = (Dm*ekm*sm**3)**0.25

    # Eq 21
    km = 0
    for x_i, kiji in zip(xi, kij):
        for x_j, k in zip(xi, kiji):
            km += x_i*x_j*k

    rho = rho/Mm/1000
    Vcm = (sm/0.809)**3                                                # Eq 16
    Tcm = 1.2593*ekm                                                   # Eq 17
    murm = 131.3*Dm/(Vcm*Tcm)**0.5                                     # Eq 22

    T_ = T/ekm

    # Table IV
    dat = [
        (2.41657, 0.74824, -0.91858, 121.72100),
        (-0.50924, -1.50936, -49.99120, 69.98340),
        (6.61069, 5.62073, 64.75990, 27.03890),
        (14.54250, -8.91387, -5.63794, 74.34350),
        (0.79274, 0.82019, -0.69369, 6.31734),
        (-5.86340, 12.80050, 9.58926, -65.52920),
        (81.17100, 114.15800, -60.84100, 466.77500)]

    # Eq 13
    B = []
    for bo, b1, b2, b3 in dat:
        B.append(bo + b1*wm + b2*murm**4 + b3*km)
    B1, B2, B3, B4, B5, B6, B7 = B

    Y = rho*Vcm/6
    G1 = (1-0.5*Y)/(1-Y)**3
    H2 = (B1*((1-exp(-B4*Y))/Y)+B2*G1*exp(B5*Y)+B3*G1)/(B1*B4+B2+B3)

    # Eq 12
    kk = ko*(1/H2 + B6*Y)
    kp = (3.039e-4*(Tcm/Mm)**0.5/Vcm**(2/3))*B7*Y**2*H2*T_**0.5

    return unidades.ThermalConductivity(kk+kp, "calscmK")


@refDoc(__doi__, [7, 3])
def ThG_StielThodosYorizane(T, xi, Tci, Pci, Vci, wi, Mi, V, ko):
    r"""Calculate thermal conductiviy of gas mixtures at high pressure using
    the Stiel-Thodos correlation for pure compound using the mixing rules
    defined by Yorizane et al.

    .. math::
        T_{cm} = \frac{\sum_i \sum_j x_ix_jV_{cij}T_{cij}}{V_{cm}}

    .. math::
        V_{cm} = \sum_i \sum_j x_ix_jV_{cij}

    .. math::
        \omega_m = \sum_i x_i\omega_i

    .. math::
        Z_{cm} = 0.291-0.08\omega_m

    .. math::
        P_{cm} = \frac{Z_{cm}RT_{cm}}{V_{cm}}

    .. math::
        M_m = \sum_i x_iM_i

    .. math::
        T_{cij} = \left(T_{ci}T_{cj}\right)^{1/2}

    .. math::
        V_{cij} = \frac{\left(V_{ci}^{1/3}+V_{cj}^{1/3}\right)^3}{8}

    Parameters
    ----------
    T : float
        Temperature, [K]
    xi : list
        Mole fractions of components, [-]
    Tci : list
        Critical temperature of compounds, [K]
    Vci : list
        Critical volume of compounds, [m³/kg]
    Zci : list
        Critical pressure of compounds, [-]
    wi :  list
        Acentric factor of compounds, [-]
    Mi : list
        Molecular weight of compounds, [g/mol]
    V : float
        Specific volume, [m³/kg]
    ko : list
        Thermal conductivities of mixture at low pressure, [W/m·K]

    Returns
    -------
    k : float
        Thermal conductivities of mixture, [W/m·K]

    Examples
    --------
    Example 10-6 from [3]_; 75.5% methane, 24.5% CO2 at 370.8K and 174.8bar

    >>> Tc = [190.56, 304.12]
    >>> Pc = [45.99e5, 73.74e5]
    >>> Vc = [98.6/16.043/1000, 94.07/44.01/1000]
    >>> M = [16.043, 44.01]
    >>> w = [0.011, 0.225]
    >>> x = [0.755, 0.245]
    >>> args = (370.8, x, Tc, Pc, Vc, w, M, 159/22.9/1000, 0.0377)
    >>> "%0.4f" % ThG_StielThodosYorizane(*args)
    '0.0527'
    """
    # Use critical volume in molar base
    Vci = [Vc*M*1000 for Vc, M in zip(Vci, Mi)]

    # Eq 8; missing rules for critical properties
    wm = sum([x*w for x, w in zip(xi, wi)])
    Mm = sum([x*M for x, M in zip(xi, Mi)])
    Zcm = 0.291-0.08*wm

    Vcij = []
    for Vc_i in Vci:
        Vciji = []
        for Vc_j in Vci:
            Vciji.append((Vc_i**(1/3)+Vc_j**(1/3))**3/8)
        Vcij.append(Vciji)

    Vcm = 0
    for x_i, Vciji in zip(xi, Vcij):
        for x_j, Vc in zip(xi, Vciji):
            Vcm += x_i*x_j*Vc

    Tcij = []
    for Tc_i in Tci:
        Tciji = []
        for Tc_j in Tci:
            Tciji.append((Tc_i*Tc_j)**0.5)
        Tcij.append(Tciji)

    Tcm = 0
    for x_i, Vciji, Tciji in zip(xi, Vcij, Tcij):
        for x_j, Vc, Tc in zip(xi, Vciji, Tciji):
            Tcm += x_i*x_j*Vc*Tc/Vcm

    Pcm = Zcm*R*Tcm/Vcm*1e6
    Vcm = Vcm/Mm/1000

    km = ThG_StielThodos(T, Tcm, Pcm, Vcm, Mm, V, ko)
    return unidades.ThermalConductivity(km)


@refDoc(__doi__, [1, 21])
def ThG_TRAPP(T, xi, Tci, Vci, Zci, wi, Mi, rho, ko):
    r"""Calculate the thermal conductivity of gas mixtures at high pressure
    using the TRAPP (TRAnsport Property Prediction) method.

    .. math::
        \lambda_m = \lambda_m^o+F_{\lambda m}X_{\lambda m}
        \left(\lambda^R-\lambda^{Ro}\right)

    .. math::
        h_m = \sum_i \sum_j x_ix_jh_{ij}

    .. math::
        f_m = \frac{\sum_i \sum_j x_ix_jf_{ij}h_{ij}}{h_m}

    .. math::
        h_{ij}=\frac{\left(h_i^{1/3}+h_j^{1/3}\right)^3}{8}

    .. math::
        f_{ij} = \left(f_if_j\right)^{1/2}

    .. math::
        T_o = T/f_m

    .. math::
        \rho_o = \rho h_m

    .. math::
        F_{\lambda m} = \frac{44.094^{1/2}}{h_m^2} \sum_i \sum_j x_ix_j
        \left(\frac{f_{ij}}{M_{ij}}\right)^{1/2}h_{ij}^{4/3}

    .. math::
        M_{ij} = \left(\frac{1}{2M_i}+\frac{1}{2M_j}\right)^{-1}

    .. math::
        X_{\lambda m} = \left[1+\frac{2.1866\left(\omega_m-\omega^R\right)}
        {1-0.505\left(\omega_m-\omega^R\right)}\right]^{1/2}

    .. math::
        \omega_m = \sum_i x_i\omega_i

    Parameters
    ----------
    T : float
        Temperature, [K]
    xi : list
        Mole fractions of components, [-]
    Tci : list
        Critical temperature of compounds, [K]
    Vci : list
        Critical volume of compounds, [m³/kg]
    Zci : list
        Critical pressure of compounds, [Pa]
    wi :  list
        Acentric factor of compounds, [-]
    Mi : list
        Molecular weight of compounds, [g/mol]
    rho : float
        Density, [kg/m3]
    ko : float
        Low-pressure gas thermal conductivity, [Pa*S]

    Returns
    -------
    k : float
        High-pressure gas thermal conductivity [W/m/k]

    Examples
    --------
    Example 10-8 from [3]_; 75.5% methane, 24.5% CO2 at 370.8K and 174.8bar

    >>> Tc = [190.56, 304.12]
    >>> Vc = [98.6/16.043/1000, 94.07/44.01/1000]
    >>> Zc = [0.286, 0.274]
    >>> M = [16.043, 44.01]
    >>> w = [0.011, 0.225]
    >>> x = [0.755, 0.245]
    >>> args = (370.8, x, Tc, Vc, Zc, w, M, 1/159*22.9*1000, 0.0377)
    >>> "%0.4f" % ThG_TRAPP(*args)
    '0.0549'
    """
    # Reference fluid properties, propane
    TcR = 369.83
    rhocR = 1/200  # mol/cm³
    ZcR = 0.276
    wR = 0.152

    # Convert volume to molar base
    Vci = [Vc*M*1000 for Vc, M in zip(Vci, Mi)]
    Mm = sum([x*M for x, M in zip(xi, Mi)])
    rho = rho/Mm/1000

    # Calculate shape factor for mixture
    fi = []
    hi = []
    for Tc, Vc, Zc, w in zip(Tci, Vci, Zci, wi):
        fi.append(Tc/TcR*(1+(w-wR)*(0.05203-0.7498*log(T/Tc))))
        hi.append(rhocR*Vc*ZcR/Zc*(1-(w-wR)*(0.1436-0.2822*log(T/Tc))))

    fij = []
    for f_i in fi:
        fiji = []
        for f_j in fi:
            fiji.append((f_i*f_j)**0.5)
        fij.append(fiji)

    hij = []
    for h_i in hi:
        hiji = []
        for h_j in hi:
            hiji.append((h_i**(1/3)+h_j**(1/3))**3/8)
        hij.append(hiji)

    hm = 0
    for x_i, hiji in zip(xi, hij):
        for x_j, h in zip(xi, hiji):
            hm += x_i*x_j*h

    fm = 0
    for x_i, hiji, fiji in zip(xi, hij, fij):
        for x_j, h, f in zip(xi, hiji, fiji):
            fm += x_i*x_j*f*h/hm

    To = T/fm
    rho0 = rho*hm

    Mij = []
    for M_i in Mi:
        Miji = []
        for M_j in Mi:
            Miji.append(1/(1/2/M_i+1/2/M_j))
        Mij.append(Miji)

    suma = 0
    for x_i, fiji, Miji, hiji in zip(xi, fij, Mij, hij):
        for x_j, f, M, h in zip(xi, fiji, Miji, hiji):
            suma += x_i*x_j*(f/M)**0.5*h**(4/3)
    Flm = 44.094**0.5/hm**2*suma

    wm = sum([x*w for x, w in zip(xi, wi)])
    Xlm = (1+2.1866*(wm-wR)/(1-0.505*(wm-wR)))**0.5

    # Coefficients in [53]_, pag 796
    # Density are in mol/dm³
    rho0 *= 1000
    rhocR *= 1000

    rhorR = rho0/rhocR
    TrR = To/TcR
    lR = 15.2583985944*rhorR + 5.29917319127*rhorR**3 + \
        (-3.05330414748+0.450477583739/TrR)*rhorR**4 + \
        (1.03144050679-0.185480417707/TrR)*rhorR**5

    return unidades.ThermalConductivity(Flm*Xlm*lR*1e-3 + ko)


@refDoc(__doi__, [8, 2])
def Tension(xi, sigmai):
    r"""Calculate surface tension of liquid nmixtures using the Morgan-Griggs
    law of mixtures method, also referenced in API procedure 10A2.1, pag 991

    .. math::
        \sigma_{m} = \sum_{i} x_i \sigma_i

    Parameters
    ----------
    xi : list
        Mole fractions of components, [-]
    sigmai : list
        Surface tension of components, [N/m]

    Returns
    -------
    sigma : float
        Surface tension of mixture, [N/m]

    Examples
    --------
    Example from [2]_ 37.9% benzene, 62.1% CycloC6 at 77F and 1atm

    >>> s1 = unidades.Tension(28.2, "dyncm")
    >>> s2 = unidades.Tension(24.3, "dyncm")
    >>> "%0.1f" % Tension([0.379, 0.621], [s1, s2]).dyncm
    '25.8'
    """
    sigma = 0
    for x, s in zip(xi, sigmai):
        sigma += x*s
    return unidades.Tension(sigma)


class Mezcla(config.Entity):
    """
    Class to model mixture calculation, components, physics properties, mixing
    rules, entity save/load layer.

    Parameters
    ----------
    tipo : int
        kind of mix definition:

            * 0 : Undefined
            * 1 : Unitary Massflow
            * 2 : Unitary Molarflow
            * 3 : Mass flow and molar fractions
            * 4 : Mass flow and mass fractions
            * 5 : Molar flow and molar fractions
            * 6 : Molar flow and mass fractions

    ids : list
        Index of component in database, [-]
    customCmp : list
        List with additional component defined out of main database
    fraccionMolar: list
        Molar fraccion list of compounds, [-]
    fraccionMasica: list
        Mass fraccion list of compounds, [-]
    caudalMasico: float
        Total mass flow, [kg/h)
    caudalMolar: float
        Total molar flow, [kmol/h)
    caudalUnitarioMasico: list
        Mass flow for each compund, [kg/h]
    caudalUnitarioMolar: list
        Molar flow for each compund, [kmol/h]

    Notes
    -----
    Additionally can define custom calculation method with the parameters:

        * *ids*: List with component index of mixture
        * *rhoLMix*: Liquid density correlation index
        * *Corr_RhoLMix*: Compressed liquid density correlation index
        * *MuLMix*: Liquid viscosity correlation index
        * *MuGMix*: Gas viscosity correlation index
        * *corr_MuGMix*: Compressed gas viscosity correlation index
        * *ThCondLMix*: Liquid thermal conductivity correlation index
        * *ThCondGMix*: Gas thermal conductivity correlation index
        * *corr_ThCondGMix*: Compressed gas thermal conductivity correlation

    These options overwrite the project configuration and the user
    configuration, for now only in API usage. Not custom stream property
    definition in main program

    Examples
    --------
    This are several ejemples of usage of this class with several configuration
    definition, obviously not all correlation return valid values.

    Liquid density methods: Example from [2]_;
    58.71% ethane, 41.29% heptane at 91ºF

    >>> xi = [0.5871, 0.4129]
    >>> c0 = Mezcla(2, ids=[3, 11], caudalUnitarioMolar=xi, RhoLMix=0)
    >>> c1 = Mezcla(2, ids=[3, 11], caudalUnitarioMolar=xi, RhoLMix=1)
    >>> args = (unidades.Temperature(91, "F"), 101325)
    >>> "%0.2f %0.2f" % (c0.RhoL(*args).kgl, c1.RhoL(*args).kgl)
    '0.57 0.57'

    Example from [3]_; 20% ethane, 80% nC10 at 166F and 3000psi

    >>> c0 = Mezcla(2, ids=[3, 14], caudalUnitarioMolar=[2, 8], RhoLPMix=0)
    >>> c1 = Mezcla(2, ids=[3, 14], caudalUnitarioMolar=[2, 8], RhoLPMix=1)
    >>> c2 = Mezcla(2, ids=[3, 14], caudalUnitarioMolar=[2, 8], RhoLPMix=2)
    >>> args = (unidades.Temperature(160, "F"), unidades.Pressure(3000, "psi"))
    >>> "%0.3f %0.3f" % (c0.RhoL(*args).kgl, c1.RhoL(*args).kgl)
    '0.678 0.681'
    >>> "%0.3f" % c2.RhoL(*args).kgl
    '0.687'

    Gas viscosity methods: Example 9-4 from [3]_; 28.6% N2, 71.4% R22 at 50ºC

    >>> c0 = Mezcla(2, ids=[46, 220], caudalUnitarioMolar=[286, 714], MuGMix=0)
    >>> c1 = Mezcla(2, ids=[46, 220], caudalUnitarioMolar=[286, 714], MuGMix=1)
    >>> c2 = Mezcla(2, ids=[46, 220], caudalUnitarioMolar=[286, 714], MuGMix=2)
    >>> c3 = Mezcla(2, ids=[46, 220], caudalUnitarioMolar=[286, 714], MuGMix=3)
    >>> c4 = Mezcla(2, ids=[46, 220], caudalUnitarioMolar=[286, 714], MuGMix=4)
    >>> args = (unidades.Temperature(50, "C"), 101325, 0)
    >>> "%0.2f %0.2f" % (c0.Mu_Gas(*args).microP, c1.Mu_Gas(*args).microP)
    '151.22 160.94'
    >>> "%0.2f %0.2f" % (c2.Mu_Gas(*args).microP, c3.Mu_Gas(*args).microP)
    '156.18 148.89'
    >>> "%0.2f" % c4.Mu_Gas(*args).microP
    '148.66'

    Example 9-14 in [3]_, 80% methane, 20% nC10 at 377.6K and 413.7bar

    >>> c0 = Mezcla(2, ids=[2, 14], caudalUnitarioMolar=[8, 2], MuGPMix=0)
    >>> c1 = Mezcla(2, ids=[2, 14], caudalUnitarioMolar=[8, 2], MuGPMix=1)
    >>> c2 = Mezcla(2, ids=[2, 14], caudalUnitarioMolar=[8, 2], MuGPMix=2)
    >>> c3 = Mezcla(2, ids=[2, 14], caudalUnitarioMolar=[8, 2], MuGPMix=3)
    >>> c4 = Mezcla(2, ids=[2, 14], caudalUnitarioMolar=[8, 2], MuGPMix=4)
    >>> args = (377.6, 413.7e5, 448.4)

    >>> "%0.2f %0.2f" % (c0.Mu_Gas(*args).muPas, c1.Mu_Gas(*args).muPas)
    '55.66 112.22'
    >>> "%0.2f %0.2f" % (c2.Mu_Gas(*args).muPas, c3.Mu_Gas(*args).muPas)
    '83.37 40.78'
    >>> "%0.2f" % c4.Mu_Gas(*args).muPas
    '39.91'

    Liquid viscosity methods: Example A from [2]_;
    29.57% nC16, 35.86% benzene, 34.57% nC6 at 77ºF

    >>> x = [0.2957, 0.3586, 0.3457]
    >>> c0 = Mezcla(2, ids=[20, 40, 10], caudalUnitarioMolar=x, MuLMix=0)
    >>> c1 = Mezcla(2, ids=[20, 40, 10], caudalUnitarioMolar=x, MuLMix=1)
    >>> args = (unidades.Temperature(77, "F"), 101325)
    >>> "%0.2f %0.2f" % (c0.Mu_Liquido(*args).cP, c1.Mu_Liquido(*args).cP)
    '0.88 0.76'

    Liquid thermal conductivity methods:
    Example from [2]_ 68% nC7, 32% CycloC5 at 32F and 1atm
    >>> x = [0.68, 0.32]
    >>> c0 = Mezcla(2, ids=[11, 36], caudalUnitarioMolar=x, ThCondLMix=0)
    >>> c1 = Mezcla(2, ids=[11, 36], caudalUnitarioMolar=x, ThCondLMix=1)
    >>> args = (unidades.Temperature(32, "F"), 101325, 0)
    >>> "%0.3f %0.3f" % (c0.ThCond_Liquido(*args), c1.ThCond_Liquido(*args))
    '0.132 0.132'

    Gas thermal conductivity methods:
    Example 10-5 from [3]_; 25% benzene, 75% Argon at 100.6ºC and 1bar

    >>> x = [0.25, 0.75]
    >>> c0 = Mezcla(2, ids=[40, 98], caudalUnitarioMolar=x, ThCondGMix=0)
    >>> c1 = Mezcla(2, ids=[40, 98], caudalUnitarioMolar=x, ThCondGMix=1)
    >>> c2 = Mezcla(2, ids=[40, 98], caudalUnitarioMolar=x, ThCondGMix=2)
    >>> args = (unidades.Temperature(100.6, "C"), 1e5, 0)
    >>> "%0.4f %0.4f" % (c0.ThCond_Gas(*args), c1.ThCond_Gas(*args))
    '0.0187 0.0197'
    >>> "%0.4f" % c2.ThCond_Gas(*args)
    '0.0224'

    Example 10-6 from [3]_; 75.5% methane, 24.5% CO2 at 370.8K and 174.8bar
    Here the Chung method is the best option to get the low pressure thermal
    conductivity of mixture

    >>> x = [0.755, 0.245]
    >>> kw = {"caudalUnitarioMolar": x, "ThCondGMix": 2}
    >>> c0 = Mezcla(2, ids=[2, 49], ThCondGPMix=0, **kw)
    >>> c1 = Mezcla(2, ids=[2, 49], ThCondGPMix=1, **kw)
    >>> c2 = Mezcla(2, ids=[2, 49], ThCondGPMix=2, **kw)
    >>> args = (370.8, 174.8e5, 1/159*22.9*1000)
    >>> "%0.4f %0.4f" % (c0.ThCond_Gas(*args), c1.ThCond_Gas(*args))
    '0.0526 0.0548'
    >>> "%0.4f" % c2.ThCond_Gas(*args)
    '0.0580'
    """
    kwargs = {"ids": [],
              "customCmp": [],
              "caudalMasico": 0.0,
              "caudalMolar": 0.0,
              "caudalUnitarioMolar": [],
              "caudalUnitarioMasico": [],
              "fraccionMolar": [],
              "fraccionMasica": [],

              "RhoLMix": None,
              "RhoLPMix": None,
              "MuLMix": None,
              "MuGMix": None,
              "MuGPMix": None,
              "ThCondLMix": None,
              "ThCondGMix": None,
              "ThCondGPMix": None
              }

    METHODS_RhoL = ["Rackett", "COSTALD"]
    METHODS_RhoLP = ["Aalto-Keskinen (1996)", "Tait-COSTALD (1982",
                     "Nasrifar (2000)", "API"]

    METHODS_MuG = ["Reichenberg (1975)", "Lucas (1984)", "Chung (1988)",
                   "Wilke (1950)", "Herning-Zipperer (1936)"]
    METHODS_MuGP = ["Lucas (1984)", "Chung (1988)", "TRAPP (1996)",
                    "Dean-Stiel (1965)", "API"]
    METHODS_MuL = ["Kendall-Monroe", "Arrhenius"]

    METHODS_ThG = ["Mason-Saxena (1958)", "Lindsay-Bromley (1950)",
                   "Chung (1988)"]
    METHODS_ThGP = ["Stiel-Thodos-Yorizane (1983)", "TRAPP", "Chung (1988)"]
    METHODS_ThL = ["Li (1976)", "Power Law"]

    def __init__(self, tipo=0, **kwargs):
        if tipo == 0:
            self._bool = False
            return

        self._bool = True

        self.kwargs = Mezcla.kwargs.copy()
        self.kwargs.update(kwargs)
        self.Config = config.getMainWindowConfig()
        if self.kwargs["ids"]:
            self.ids = self.kwargs.get("ids")
        elif self.kwargs["customCmp"]:
            # Initialice ids variable to avoid acumulate in several instances
            self.ids = []
        else:
            txt = self.Config.get("Components", "Components")
            if isinstance(txt, str):
                self.ids = eval(txt)
            else:
                self.ids = txt
        self.componente = [Componente(int(i), **kwargs) for i in self.ids]
        for cmp in self.kwargs["customCmp"]:
            self.componente.append(cmp)
            self.ids.append(0)

        fraccionMolar = self.kwargs.get("fraccionMolar", None)
        fraccionMasica = self.kwargs.get("fraccionMasica", None)
        caudalMasico = self.kwargs.get("caudalMasico", None)
        caudalMolar = self.kwargs.get("caudalMolar", None)
        caudalUnitarioMasico = self.kwargs.get("caudalUnitarioMasico", None)
        caudalUnitarioMolar = self.kwargs.get("caudalUnitarioMolar", None)

        # normalizce fractions to unit
        if fraccionMolar:
            suma = float(sum(fraccionMolar))
            fraccionMolar = [x/suma for x in fraccionMolar]
        if fraccionMasica:
            suma = float(sum(fraccionMasica))
            fraccionMasica = [x/suma for x in fraccionMasica]

        # calculate all concentration units
        if tipo == 1:
            kw = mix_unitmassflow(caudalUnitarioMasico, self.componente)
        elif tipo == 2:
            kw = mix_unitmolarflow(caudalUnitarioMolar, self.componente)
        elif tipo == 3:
            kw = mix_massflow_molarfraction(
                caudalMasico, fraccionMolar, self.componente)
        elif tipo == 4:
            kw = mix_massflow_massfraction(
                caudalMasico, fraccionMasica, self.componente)
        elif tipo == 5:
            kw = mix_molarflow_molarfraction(
                caudalMolar, fraccionMolar, self.componente)
        elif tipo == 6:
            kw = mix_molarflow_massfraction(
                caudalMolar, fraccionMasica, self.componente)

        caudalMolar = kw["molarFlow"]
        caudalMasico = kw["massFlow"]
        caudalUnitarioMasico = kw["unitMassFlow"]
        caudalUnitarioMolar = kw["unitMolarFlow"]
        fraccionMolar = kw["molarFraction"]
        fraccionMasica = kw["massFraction"]

#        # Clean component with null composition
#        self.zeros = []
#        for i, x in enumerate(fraccionMolar):
#            if not x:
#                self.zeros.append(i)
#
#        for i in self.zeros[::-1]:
#            del self.ids[i]
#            del self.componente[i]
#            del fraccionMolar[i]
#            del fraccionMasica[i]
#            del caudalUnitarioMasico[i]
#            del caudalUnitarioMolar[i]

        self.fraccion = [unidades.Dimensionless(f) for f in fraccionMolar]
        self.caudalmasico = unidades.MassFlow(caudalMasico)
        self.caudalmolar = unidades.MolarFlow(caudalMolar)
        self.fraccion_masica = [unidades.Dimensionless(f)
                                for f in fraccionMasica]
        self.caudalunitariomasico = [unidades.MassFlow(q)
                                     for q in caudalUnitarioMasico]
        self.caudalunitariomolar = [unidades.MolarFlow(q)
                                    for q in caudalUnitarioMolar]
        self.M = unidades.Dimensionless(caudalMasico/caudalMolar)

        if tipo == 0:
            self._bool = False
            self.status = 0
            return
        else:
            self._bool = True
            self.status = 1

        # Calculate critic temperature, API procedure 4B1.1 pag 304
        V = sum([xi*cmp.Vc for xi, cmp in zip(self.fraccion, self.componente)])
        k = [xi*cmp.Vc/V for xi, cmp in zip(self.fraccion, self.componente)]
        Tcm = sum([ki*cmp.Tc for ki, cmp in zip(k, self.componente)])
        self.Tc = unidades.Temperature(Tcm)

        # Calculate pseudocritic temperature
        t = sum([xi*cmp.Tc for xi, cmp in zip(self.fraccion, self.componente)])
        self.tpc = unidades.Temperature(t)

        # Calculate pseudocritic pressure
        p = sum([xi*cmp.Pc for xi, cmp in zip(self.fraccion, self.componente)])
        self.ppc = unidades.Pressure(p)

        # Calculate critic pressure, API procedure 4B2.1 pag 307
        sumaw = 0
        for xi, cmp in zip(self.fraccion, self.componente):
            sumaw += xi*cmp.f_acent
        pc = self.ppc+self.ppc*(5.808+4.93*sumaw)*(self.Tc-self.tpc)/self.tpc
        self.Pc = unidades.Pressure(pc)

        # Calculate acentric factor, API procedure 6B2.2-6 pag 523
        self.f_acent = sum([xi*cmp.f_acent for xi, cmp in
                            zip(self.fraccion, self.componente)])
        self.f_acent_mod = sum([xi*cmp.f_acent_mod for xi, cmp in
                                zip(self.fraccion, self.componente)])

        # Calculate critic volume
        Vci = self._arraylize("Vc")
        Mi = self._arraylize("M")
        hc = self._arraylize("isHydrocarbon")
        self.Vc = Vc_ChuehPrausnitz(self.fraccion, Vci, Mi, hydrocarbon=hc)

        tb = [xi*cmp.Tb for xi, cmp in zip(self.fraccion, self.componente)]
        self.Tb = unidades.Temperature(sum(tb))
        self.SG = sum([xi*cmp.SG for xi, cmp in
                       zip(self.fraccion, self.componente)])

    def __call__(self):
        pass

    def _arraylize(self, prop, unit=None):
        """Get the compounds property prop as list
        prop: a string code with the property to return
            f_acent, M, Vc, Tc,...
        """
        array = []
        for cmp in self.componente:
            value = cmp.__getattribute__(prop)
            if unit:
                value = value.__getattribute__(unit)
            array.append(value)
        return array

    @refDoc(__doi__, [2], tab=8)
    def _Ho(self, T):
        r"""Ideal gas enthalpy, referenced in API procedure 7B4.1, pag 645

        .. math::
            H_m^o = \sum_i x_wH_i^o

        Parameters
        ----------
        T : float
            Temperature, [K]
        """
        h = 0
        for xw, cmp in zip(self.fraccion_masica, self.componente):
            h += xw*cmp._Ho(T)
        return unidades.Enthalpy(h)

    @refDoc(__doi__, [2], tab=8)
    def _so(self, T):
        r"""
        Ideal gas entropy, referenced in API procedure 7F2.1, pag 741

        .. math::
            S_m^o = \sum_i x_wS_i^o - \frac{R}{M} x_i\lnx_i

        Parameters
        ----------
        T : float
            Temperature, [K]
        """
        s = 0
        for x, xw, cmp in zip(
                self.fraccion, self.fraccion_masica, self.componente):
            s += xw*cmp._So(T) + R/cmp.M*x*log(x)
        return unidades.SpecificHeat(s)

    def Cp_Gas(self, T, P):
        """Calculate specific heat from gas, API procedure 7D4.1, pag 714"""
        Cp = 0
        for xi, cmp in zip(self.fraccion_masica, self.componente):
            Cp += xi*cmp.Cp_Gas_DIPPR(T)
        return unidades.SpecificHeat(Cp)

    def Cp_Liquido(self, T):
        """Calculate specific heat from liquid, API procedure 7D1.9, pag 714"""
        Cp = 0
        for xi, cmp in zip(self.fraccion_masica, self.componente):
            Cp += xi*cmp.Cp_Liquido(T)
        return unidades.SpecificHeat(Cp)

    def RhoL(self, T, P):
        """Calculate the density of liquid phase using any of available
        correlation"""
        method = self.kwargs["RhoLMix"]
        if method is None or method >= len(Mezcla.METHODS_RhoL):
            method = self.Config.getint("Transport", "RhoLMix")
        Pcorr = self.kwargs["RhoLPMix"]
        if Pcorr is None or method >= len(Mezcla.METHODS_RhoLP):
            Pcorr = self.Config.getint("Transport", "Corr_RhoLMix")

        # Calculate of low pressure viscosity
        if method == 0:
            Tci = self._arraylize("Tc")
            Pci = self._arraylize("Pc")
            Vci = self._arraylize("Vc")
            Zrai = self._arraylize("rackett")
            Mi = self._arraylize("M")
            rhos = RhoL_RackettMix(T, self.fraccion, Tci, Pci, Vci, Zrai, Mi)
        elif method == 1:
            Tci = self._arraylize("Tc")
            Vci = self._arraylize("Vc")
            wi = self._arraylize("f_acent")
            Mi = self._arraylize("M")
            rhos = RhoL_CostaldMix(T, self.fraccion, Tci, wi, Vci, Mi)

        # Add correction factor for high pressure
        if P < 1e6:
            rho = rhos
        elif Pcorr == 0:
            Tci = self._arraylize("Tc")
            Pci = self._arraylize("Pc")
            Vci = self._arraylize("Vc")
            Mi = self._arraylize("M")
            wi = self._arraylize("f_acent")
            Mi = self._arraylize("M")
            rho = RhoL_AaltoKeskinenMix(
                T, P, self.fraccion, Tci, Pci, Vci, wi, Mi, rhos)
        elif Pcorr == 1:
            Tci = self._arraylize("Tc")
            Vci = self._arraylize("Vc")
            Mi = self._arraylize("M")
            wi = self._arraylize("f_acent")
            Mi = self._arraylize("M")
            rho = RhoL_TaitCostaldMix(
                T, P, self.fraccion, Tci, Vci, wi, Mi, rhos)
        elif Pcorr == 2:
            Tci = self._arraylize("Tc")
            Vci = self._arraylize("Vc")
            Mi = self._arraylize("M")
            wi = self._arraylize("f_acent")
            Mi = self._arraylize("M")
            rho = RhoL_NasrifarMix(T, P, self.fraccion, Tci, Vci, wi, Mi, rhos)
        elif Pcorr == 3:
            Tci = self._arraylize("Tc")
            Mi = self._arraylize("M")
            wi = self._arraylize("f_acent")
            Mi = self._arraylize("M")
            rho = RhoL_APIMix(T, P, self.fraccion, Tci, Pci, rhos)

        return rho

    def Mu_Gas(self, T, P, rho):
        """General method for calculate viscosity of gas"""
        method = self.kwargs["MuGMix"]
        if method is None or method >= len(Mezcla.METHODS_MuG):
            method = self.Config.getint("Transport", "MuGMix")
        Pcorr = self.kwargs["MuGPMix"]
        if Pcorr is None or method >= len(Mezcla.METHODS_MuGP):
            Pcorr = self.Config.getint("Transport", "Corr_MuGMix")

        # Calculate of low pressure viscosity
        if method == 0:
            Tci = self._arraylize("Tc")
            Pci = self._arraylize("Pc")
            Mi = self._arraylize("M")
            Mi = self._arraylize("M")
            Di = self._arraylize("dipole", "Debye")
            mui = [cmp.Mu_Gas(T, 101325, rho) for cmp in self.componente]
            muo = MuG_Reichenberg(T, self.fraccion, Tci, Pci, Mi, mui, Di)
        elif method == 1:
            Tci = self._arraylize("Tc")
            Pci = self._arraylize("Pc")
            Vci = self._arraylize("Vc")
            Zci = self._arraylize("Zc")
            Mi = self._arraylize("M")
            Di = self._arraylize("dipole", "Debye")
            muo = MuG_Lucas(
                T, 101325, self.fraccion, Tci, Pci, Vci, Zci, Mi, Di)
        elif method == 2:
            Tci = self._arraylize("Tc")
            Vci = self._arraylize("Vc")
            Mi = self._arraylize("M")
            wi = self._arraylize("f_acent")
            Di = self._arraylize("dipole", "Debye")
            ki = [cmp._K_Chung() for cmp in self.componente]
            muo = MuG_Chung(T, self.fraccion, Tci, Vci, Mi, wi, Di, ki)
        elif method == 3:
            Mi = self._arraylize("M")
            mui = [cmp.Mu_Gas(T, 101325, None) for cmp in self.componente]
            muo = MuG_Wilke(self.fraccion, Mi, mui)
        elif method == 4:
            Mi = self._arraylize("M")
            mui = [cmp.Mu_Gas(T, 101325, None) for cmp in self.componente]
            muo = MuG_Herning(self.fraccion, Mi, mui)

        # Add correction factor for high pressure
        if P < 1e6:
            mu = muo
        elif Pcorr == 0:
            Tci = self._arraylize("Tc")
            Pci = self._arraylize("Pc")
            Vci = self._arraylize("Vc")
            Zci = self._arraylize("Zc")
            Mi = self._arraylize("M")
            Di = self._arraylize("dipole", "Debye")
            mu = MuG_Lucas(T, P, self.fraccion, Tci, Pci, Vci, Zci, Mi, Di)
        elif Pcorr == 1:
            Tci = self._arraylize("Tc")
            Vci = self._arraylize("Vc")
            Mi = self._arraylize("M")
            wi = self._arraylize("f_acent")
            Di = self._arraylize("dipole", "Debye")
            ki = [cmp._K_Chung() for cmp in self.componente]
            mu = MuG_P_Chung(
                T, self.fraccion, Tci, Vci, Mi, wi, Di, ki, rho, muo)
        elif Pcorr == 2:
            Tci = self._arraylize("Tc")
            Vci = self._arraylize("Vc")
            Zci = self._arraylize("Zc")
            Mi = self._arraylize("M")
            wi = self._arraylize("f_acent")
            mu = MuG_TRAPP(
                T, P, self.fraccion, Tci, Vci, Zci, Mi, wi, rho, muo)
        elif Pcorr == 3:
            Tci = self._arraylize("Tc")
            Pci = self._arraylize("Pc")
            Vci = self._arraylize("Vc")
            Mi = self._arraylize("M")
            rhoc = 1/Vc_ChuehPrausnitz(self.fraccion, Vci, Mi)
            mu = MuG_DeanStielMix(self.fraccion, Tci, Pci, Mi, rhoc, rho, muo)
        elif Pcorr == 4:
            Tci = self._arraylize("Tc")
            Pci = self._arraylize("Pc")
            mu = MuG_APIMix(T, P, self.fraccion, Tci, Pci, muo)

        return mu

    def Mu_Liquido(self, T, P):
        """General method for calculate viscosity of Liquid"""
        method = self.kwargs["MuLMix"]
        if method is None or method >= len(Mezcla.METHODS_MuL):
            method = self.Config.getint("Transport", "MuLMix")

        if method == 0:
            mui = [cmp.Mu_Liquido(T, P) for cmp in self.componente]
            mu = MuL_KendallMonroe(self.fraccion, mui)
        elif method == 1:
            Mi = self._arraylize("M")
            mui = [cmp.Mu_Liquido(T, P) for cmp in self.componente]
            mu = MuL_Chemcad(self.fraccion, Mi, mui)

        return mu

    def Tension(self, T):
        """General method for calculate surface tension"""
        sigmai = [cmp.Tension(T) for cmp in self.componente]
        tension = Tension(self.fraccion, sigmai)
        return unidades.Tension(tension)

    def ThCond_Liquido(self, T, P, rho):
        """General method for calculate thermal conductivity of liquid"""
        method = self.kwargs["ThCondLMix"]
        if method is None or method >= len(Mezcla.METHODS_ThL):
            method = self.Config.getint("Transport", "ThCondLMix")

        if method == 0:
            Vi = [1/cmp.RhoL(T, P) for cmp in self.componente]
            Mi = self._arraylize("M")
            ki = [cmp.ThCond_Liquido(T, P, rho) for cmp in self.componente]
            k = ThL_Li(self.fraccion, Vi, Mi, ki)
        elif method == 1:
            ki = [cmp.ThCond_Liquido(T, P, rho) for cmp in self.componente]
            k = ThL_Power(self.fraccion_masica, ki)

        return k

    def ThCond_Gas(self, T, P, rho):
        """General method for calculate thermal conductivity of gas"""
        method = self.kwargs["ThCondGMix"]
        if method is None or method >= len(Mezcla.METHODS_ThG):
            method = self.Config.getint("Transport", "ThCondGMix")
        Pcorr = self.kwargs["ThCondGPMix"]
        if Pcorr is None or method >= len(Mezcla.METHODS_ThGP):
            Pcorr = self.Config.getint("Transport", "Corr_ThCondGMix")

        # Calculate of low pressure viscosity
        if method == 0:
            Mi = self._arraylize("M")
            mui = [cmp.Mu_Gas(T, 101325, rho) for cmp in self.componente]
            ki = [cmp.ThCond_Gas(T, P, rho) for cmp in self.componente]
            ko = ThG_MasonSaxena(self.fraccion, Mi, mui, ki)
        elif method == 1:
            Mi = self._arraylize("M")
            Tbi = self._arraylize("Tb")
            mui = [cmp.Mu_Gas(T, 101325, rho) for cmp in self.componente]
            ki = [cmp.ThCond_Gas(T, P, rho) for cmp in self.componente]
            ko = ThG_LindsayBromley(T, self.fraccion, Mi, Tbi, mui, ki)
        elif method == 2:
            Tci = self._arraylize("Tc")
            Vci = self._arraylize("Vc")
            Mi = self._arraylize("M")
            wi = self._arraylize("f_acent")
            Cvi = [cmp.Cv(T) for cmp in self.componente]
            Di = self._arraylize("dipole", "Debye")
            ki = [cmp._K_Chung() for cmp in self.componente]
            mu = MuG_Chung(T, self.fraccion, Tci, Vci, Mi, wi, Di, ki)
            ko = ThG_Chung(T, self.fraccion, Tci, Vci, Mi, wi, Cvi, mu)

        # Add correction factor for high pressure
        if P < 1e6:
            k = ko
        elif Pcorr == 0:
            Tci = self._arraylize("Tc")
            Pci = self._arraylize("Pc")
            Vci = self._arraylize("Vc")
            wi = self._arraylize("f_acent")
            Mi = self._arraylize("M")
            k = ThG_StielThodosYorizane(
                T, self.fraccion, Tci, Pci, Vci, wi, Mi, 1/rho, ko)
        elif Pcorr == 1:
            Tci = self._arraylize("Tc")
            Vci = self._arraylize("Vc")
            Zci = self._arraylize("Zc")
            wi = self._arraylize("f_acent")
            Mi = self._arraylize("M")
            k = ThG_TRAPP(T, self.fraccion, Tci, Vci, Zci, wi, Mi, rho, ko)
        elif Pcorr == 2:
            Tci = self._arraylize("Tc")
            Vci = self._arraylize("Vc")
            Mi = self._arraylize("M")
            wi = self._arraylize("f_acent")
            Di = self._arraylize("dipole", "Debye")
            ki = [cmp._K_Chung() for cmp in self.componente]
            k = ThG_P_Chung(
                T, self.fraccion, Tci, Vci, Mi, wi, Di, ki, rho, ko)

        return k

    # Entity related functionality
    def writeStatetoJSON(self, state):
        mezcla = {}
        if self._bool:
            mezcla["ids"] = self.ids
            mezcla["fraction"] = self.fraccion
            mezcla["massFraction"] = self.fraccion_masica
            mezcla["massUnitFlow"] = self.caudalunitariomasico
            mezcla["molarUnitFlow"] = self.caudalunitariomolar
            mezcla["massFlow"] = self.caudalmasico
            mezcla["molarFlow"] = self.caudalmolar

            mezcla["M"] = self.M
            mezcla["Tc"] = self.Tc
            mezcla["tpc"] = self.tpc
            mezcla["ppc"] = self.ppc
            mezcla["Pc"] = self.Pc
            mezcla["w"] = self.f_acent
            mezcla["wm"] = self.f_acent_mod
            mezcla["Vc"] = self.Vc
            mezcla["Tb"] = self.Tb
            mezcla["SG"] = self.SG
        state["mezcla"] = mezcla

    def readStatefromJSON(self, mezcla):
        if mezcla:
            self._bool = True
            self.ids = mezcla["ids"]
            self.componente = [Componente(int(i)) for i in self.ids]
            self.fraccion = [
                unidades.Dimensionless(x) for x in mezcla["fraction"]]
            self.fraccion_masica = [
                unidades.Dimensionless(x) for x in mezcla["massFraction"]]
            self.caudalunitariomasico = [
                unidades.MassFlow(x) for x in mezcla["massUnitFlow"]]
            self.caudalunitariomolar = [
                unidades.MolarFlow(x) for x in mezcla["molarUnitFlow"]]
            self.caudalmasico = unidades.MassFlow(mezcla["massFlow"])
            self.caudalmolar = unidades.MolarFlow(mezcla["molarFlow"])

            self.M = unidades.Dimensionless(mezcla["M"])
            self.Tc = unidades.Temperature(mezcla["Tc"])
            self.tpc = unidades.Temperature(mezcla["tpc"])
            self.ppc = unidades.Pressure(mezcla["ppc"])
            self.Pc = unidades.Pressure(mezcla["Pc"])
            self.f_acent = unidades.Dimensionless(mezcla["w"])
            self.f_acent_mod = unidades.Dimensionless(mezcla["wm"])
            self.Vc = unidades.SpecificVolume(mezcla["Vc"])
            self.Tb = unidades.Temperature(mezcla["Tb"])
            self.SG = unidades.Dimensionless(mezcla["SG"])

    def recallZeros(self, lista, val=0):
        """Method to return any list with null component added"""
        l = lista[:]
        for i in self.zeros:
            l.insert(i, val)
        return l

# TODO:
# ThL_Rowley
# Parachor method for tension


if __name__ == '__main__':
    # T = unidades.Temperature(400, "F")
    # mezcla = Mezcla(1, ids=[4, 40], caudalUnitarioMasico=[26.92, 73.08])
    # print(mezcla.Tension(T))

    mezcla = Mezcla(2, ids=[1, 2, 40, 41], caudalUnitarioMolar=[
        0.004397674808848511, 0.057137022156057246, 0.7079892468704139, 0.23047605616468061])
    P = unidades.Pressure(485, "psi")
    T = unidades.Temperature(100, "F")
    print(mezcla.RhoL(T, P))
