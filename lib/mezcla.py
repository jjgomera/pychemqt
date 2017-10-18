#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>."""


###############################################################################
# Module with mixture definition
#   -Mixture definition procedures:
#      _mix_from_unitmassflow
#      _mix_from_unitmolarflow
#      _mix_from_massflow_and_molarfraction
#      _mix_from_massflow_and_massfraction
#      _mix_from_molarflow_and_molarfraction
#      _mix_from_molarflow_and_massfraction

#   -Mezcla: Mixture related calculation
###############################################################################


# TODO:
# ThG_Chung
# ThL_Rowley


from scipy import roots, log, sqrt, log10, exp, sin, zeros

from lib.compuestos import (Componente, RhoL_Costald, RhoL_AaltoKeskinen,
                            ThG_StielThodos)
from lib.physics import R_atml, R
from lib import unidades, config
from lib.elemental import Elemental

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
         "doi": "10.1002/aic.690220520 "},
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
        {"autor": "",
         "title": "",
         "ref": "",
         "doi": ""},
}


def _mix_from_unitmassflow(unitMassFlow, cmps):
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


def _mix_from_unitmolarflow(unitMolarFlow, cmps):
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


def _mix_from_massflow_and_molarfraction(massFlow, molarFraction, cmps):
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


def _mix_from_massflow_and_massfraction(massFlow, massFraction, cmps):
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


def _mix_from_molarflow_and_molarfraction(molarFlow, molarFraction, cmps):
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


def _mix_from_molarflow_and_massfraction(molarFlow, massFraction, cmps):
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


# Liquid density correlations
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

    References
    ----------
    .. [9] Spencer, C.F., Danner, R.P. Prediction of Bubble-Point Density of
        Mixtures. J. Chem. Eng. Data 18(2) (1973) 230-234
    .. [2] API. Technical Data book: Petroleum Refining 6th Edition
    .. [3] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
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

    Mm = 0
    for x, M in zip(xi, Mi):
        Mm += x*M
    return unidades.Density(Mm/V)


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


    References
    ----------
    .. [10] Hankinson, R.W., Thomson, G.H. A New Correlation for Saturated
        Densities of Liquids and Their Mixtures. AIChE Journal 25(4) (1979)
        653-663
    .. [2] API. Technical Data book: Petroleum Refining 6th Edition
    .. [3] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
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

    Mm = 0
    for x, M in zip(xi, Mi):
        Mm += x*M

    # Apply the pure component procedure with the mixing parameters
    rho = RhoL_Costald(T, Tcm, wm, Vcm)
    return unidades.Density(rho*Mm)


def RhoL_AaltoKeskinenMix(T, P, xi, Tci, Pci, Vci, wi, Mi, Ps, rhos):
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
    Ps : float
        Saturation pressure, [Pa]
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
    >>> args = (344.26, P, xi, Tci, Pci, Vci, wi, Mi, 9.2e5, 1/116.43*1000)
    >>> "%0.2f" % (1/RhoL_AaltoKeskinenMix(*args).gcc)
    '99.05'

    References
    ----------
    .. [11] Aalto, M., Keskinen, K.I., Aittamaa, J., Liukkonen, S. An Improved
        Correlation for Compressed Liquid Densities of Hydrocarbons. Part 2.
        Mixtures. Fluid Phase Equilibria 114 (1996) 21-35
    .. [3] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    """
    # Convert critical volumes to molar base
    Vci = [Vc*M for Vc, M in zip(Vci, Mi)]

    # Apply mixing rules
    # Eq A5
    wm = 0
    for x, w in zip(xi, wi):
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

    Mm = 0
    for x, M in zip(xi, Mi):
        Mm += x*M

    # Eq 4
    Pcm = (0.291-0.08*wm)*R*Tcm/Vcm*1000

    # Apply the pure component procedure with the mixing parameters
    rho = RhoL_AaltoKeskinen(T, P, Tcm, Pcm, wm, Ps, rhos)
    return unidades.Density(rho/Mm)


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

    >>> from scipy import r_
    >>> from lib.mEoS import CO2, O2, CO, H2, CH4, N2, C2
    >>> Mi = [CO2.M, O2.M, CO.M, H2.M, CH4.M, N2.M, C2.M]
    >>> xi = [0.035, 0.003, 0.273, 0.144, 0.037, 0.5, 0.008]
    >>> mui = r_[147.2, 201.9, 174.9, 87.55, 108.7, 178.1, 90.9]
    >>> "%0.7f" % MuG_Wilke(xi, Mi, mui*1e-9).dynscm2
    '0.0001711'

    Example A from [5]_; 58.18% H2, 41.82% propane at 77ºF and 14.7 psi

    >>> mu = MuG_Wilke([0.5818, 0.4182], [2.02, 44.1], [8.91e-6, 8.22e-6])
    >>> "%0.4f" % mu.cP
    '0.0092'

    Example B from [5]_ 95.6% CH4, 3.6% C2, 0.5% C3, 0.3% N2

    >>> x = [0.956, 0.036, 0.005, 0.003]
    >>> Mi = [16.04, 30.07, 44.1, 28.01]
    >>> mui = [1.125e-5, 9.5e-6, 8.4e-6, 1.79e-5]
    >>> "%0.5f" % MuG_Wilke(x, Mi, mui).cP
    '0.01117'

    Example 9-5 from [3]_, methane 69.7%, n-butane 30.3%

    >>> mui = [109.4e-7, 72.74e-7]
    >>> "%0.2f" % MuG_Wilke([0.697, 0.303], [16.043, 58.123], mui).microP
    '92.25'

    References
    ----------
    .. [1] Wilke, C.R. A Viscosity Equation for Gas Mixtures. J. Chem. Phys.
        18(4) (1950) 517-519
    .. [2] API. Technical Data book: Petroleum Refining 6th Edition
    .. [3] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
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


# Liquid thermal conductivity correlations
def ThL_Li(xi, Vi, ki):
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
        Molar volumes of components, [m³/kmol]
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
    >>> "%0.5f" % ThL_Li([0.68, 0.32], [V1, V2], [k1, k2]).BtuhftF
    '0.07751'

    References
    ----------
    .. [4] Li, C.C. Thermal Conductivity of Liquid Mixtures. AIChE Journal
        22(5) (1976) 927-930
    .. [2] API. Technical Data book: Petroleum Refining 6th Edition
    """
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

    References
    ----------
    .. [3] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    """
    km = 0
    for w, k in zip(wi, ki):
        km += w/k**2

    km = km**-0.5
    return unidades.ThermalConductivity(km)




# Gas thermal conductivity correlations
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

    References
    ----------
    .. [5] Lindsay, A.L., Bromley, L.A. Thermal Conductivity of Gas Mixtures.
        Ind. & Eng. Chem. 42(8) (1950) 1508-1511
    .. [2] API. Technical Data book: Petroleum Refining 6th Edition
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

    References
    ----------
    .. [6] Mason, E.A., Saxena, S.C. Approximate Formula for the Thermal
        Conductivity of Gas Mixtures. Fhys. Fluids 1(5) (1958) 361-369
    .. [3] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
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

    References
    ----------
    .. [7] Yorizane, M., Yoshiumra, S., Masuoka, H., Yoshida, H. Thermal
        Conductivities of Binary Gas Mixtures at High Pressures: N2-O2, N2-Ar,
        CO2-Ar, CO2-CH4. Ind. Eng. Chem. Fundam. 22(4) (1983) 458-462
    .. [3] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
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
    '0.0550'

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
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
        hi.append(rhocR*Vc*ZcR/Zc*(1+(w-wR)*(0.1436-0.2882*log(T/Tc))))

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

    References
    ----------
    .. [8] Livingston, J.k Morgan, R., Griggs, M.A. The Properties of Mixed
        Liquids III. The Law of Mixtures I. J. Am. Chem. Soc. 39 (1917)
        2261-2275
    .. [2] API. Technical Data book: Petroleum Refining 6th Edition
    """
    sigma = 0
    for x, s in zip(xi, sigmai):
        sigma += x*s
    return unidades.Tension(sigma)




class Mezcla(config.Entity):
    """
    Class to model mixure calculation, component, physics properties, mix rules
    Parameter:
        tipo: kind of mix definition
            0 - Undefined
            1 - Unitary Massflow
            2 - Unitary Molarflow
            3 - Mass flow and molar fractions
            4 - Mass flow and mass fractions
            5 - Molar flow and molar fractions
            6 - Molar flow and mass fractions
        kwargs: any of this variable for define mixture
            fraccionMolar
            fraccionMasica
            caudalMasico
            caudalMolar
            caudalUnitarioMasico
            caudalUnitarioMolar
    """
    kwargs = {"caudalMasico": 0.0,
              "caudalMolar": 0.0,
              "caudalUnitarioMolar": [],
              "caudalUnitarioMasico": [],
              "fraccionMolar": [],
              "fraccionMasica": []}

    def __init__(self, tipo=0, **kwargs):
        if tipo == 0:
            self._bool = False
            return

        self._bool = True

        self.kwargs = Mezcla.kwargs.copy()
        self.kwargs.update(kwargs)
        if self.kwargs["ids"]:
            self.ids = self.kwargs.get("ids")
        else:
            Config = config.getMainWindowConfig()
            txt = Config.get("Components", "Components")
            if isinstance(txt, str):
                self.ids = eval(txt)
            else:
                self.ids = txt
        self.componente = [Componente(int(i)) for i in self.ids]
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
            kw = _mix_from_unitmassflow(caudalUnitarioMasico, self.componente)
        elif tipo == 2:
            kw = _mix_from_unitmolarflow(caudalUnitarioMolar, self.componente)
        elif tipo == 3:
            kw = _mix_from_massflow_and_molarfraction(
                caudalMasico, fraccionMolar, self.componente)
        elif tipo == 4:
            kw = _mix_from_massflow_and_massfraction(
                caudalMasico, fraccionMasica, self.componente)
        elif tipo == 5:
            kw = _mix_from_molarflow_and_molarfraction(
                caudalMolar, fraccionMolar, self.componente)
        elif tipo == 6:
            kw = _mix_from_molarflow_and_massfraction(
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

        mixing = [self.Mix_van_der_Waals, self.Mix_Stryjek_Vera,
                  self.Mix_Panagiotopoulos, self.Mix_Melhem]
        conf = config.getMainWindowConfig().getint("Thermo", "Mixing")
        self.Mixing_Rule = mixing[conf]

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

        # Calculate critic volume, API procedure 4B3.1 pag 314
        sumaxvc23 = sum([xi*cmp.Vc**(2./3) for xi, cmp in
                         zip(self.fraccion, self.componente)])
        k = [xi*cmp.Vc**(2./3)/sumaxvc23 for xi, cmp in
             zip(self.fraccion, self.componente)]

        # TODO: Calculate C value from component type.
        # For now it suppose all are hidrycarbon (C=0)
        C = 0

        V = [[-1.4684*abs((cmpi.Vc-cmpj.Vc)/(cmpi.Vc+cmpj.Vc))+C
              for cmpj in self.componente] for cmpi in self.componente]
        v = [[V[i][j]*(cmpi.Vc+cmpj.Vc)/2. for j, cmpj in enumerate(
            self.componente)] for i, cmpi in enumerate(self.componente)]
        suma1 = sum([ki*cmp.Vc for ki, cmp in zip(k, self.componente)])
        suma2 = sum([ki*kj*v[i][j] for j, kj in enumerate(k)
                     for i, ki in enumerate(k)])
        self.Vc = unidades.SpecificVolume((suma1+suma2)*self.M)

        tb = [xi*cmp.Tb for xi, cmp in zip(self.fraccion, self.componente)]
        self.Tb = unidades.Temperature(sum(tb))
        self.SG = sum([xi*cmp.SG for xi, cmp in
                       zip(self.fraccion, self.componente)])

    def __call__(self):
        pass

    def recallZeros(self, lista, val=0):
        """Method to return any list with null component added"""
        l = lista[:]
        for i in self.zeros:
            l.insert(i, val)
#        return l

    def tr(self, T):
        """reduced temperature"""
        return T/self.tpc

    def pr(self, P):
        """reduced pressure"""
        return P/self.ppc

    def Kij(self, T=0, EOS=None):
        """Calculate binary interaction matrix for component of mixture,
        use bip data from dat/bip directory
        Parameter:
            T: opcional temperatura for generalized method
            EOS: name of equation of state, bwrs, nrtl, pr, srk, uniq, wils
        API procedure 8D1.1 pag 819, equations pag 827"""
        # TODO: import data here from file and remove lib/bip
        if EOS:
            kij = []
            for i in self.ids:
                kijk = []
                for j in self.ids:
                    if i == j:
                        kijk.append(0)
                    else:
                        for indice in EOS:
                            if i in indice[0:2] and j in indice[0:2]:
                                kijk.append(indice[2])
                                break
                        else:
                            if i == 1 or j == 1:
                                kijk.append(1/(344.23*exp(-0.48586*T/databank.base_datos[1][3])+1))
                            elif i == 2 or j == 2:
                                kijk.append(0.014*(abs(self.componente[self.ids.index(i)].SolubilityParameter-self.componente[self.ids.index(j)].SolubilityParameter)))
                            elif i == 46 or j == 46:
                                kijk.append(0.0403*(abs(self.componente[self.ids.index(i)].SolubilityParameter-self.componente[self.ids.index(j)].SolubilityParameter)))
                            elif i == 48 or j == 48:
                                kijk.append(0)
                            elif i == 49 or j == 49:
                                kijk.append(0.1)
                            elif i == 50 or j == 50:
                                kijk.append(0.0316*(abs(self.componente[self.ids.index(i)].SolubilityParameter-self.componente[self.ids.index(j)].SolubilityParameter)))
                            else:
                                kijk.append(0)
                kij.append(kijk)
        else:
            kij = zeros((len(self.ids), len(self.ids)))
        return kij

    def _Critical_API(self):
        """Método de cálculo de las propiedades críticas, haciendo uso de la ecuación de estado de Soave-Redlich-Kwong Modificada Thorwart-Daubert, API procedure 4B4.1 pag 317"""

        def parameters(T):
            ai=[]
            bi=[]
            for componente in self.componente:
                a, b=eos.SRK_Thorwart_lib(componente, T)
                ai.append(a)
                bi.append(b)
            b=sum([fraccion*b for fraccion, b in zip(self.fraccion, bi)])

            k=self.Kij(srk)

            aij=[[(ai[i]*ai[j])**0.5*(1-k[i][j]) for j in range(len(self.componente))] for i in range(len(self.componente))]
            a=sum([fraccioni*fraccionj*aij[i][j] for j, fraccionj in enumerate(self.fraccion) for i, fraccioni in enumerate(self.fraccion)])
            a_=[2*sum([fraccion*a for fraccion, a in zip(self.fraccion, aij[i])]) for i in range(len(self.fraccion))]

            return ai, bi, b, k, aij, a, a_


        def q(T, V):
            """Subrutina de cálculo del determinante de Q, eq 4B4.1-5"""
            ai, bi, b, k, aij, a, a_=parameters(T)
            B1=[[2*a*bi[i]*bi[j]-b*(a_[i]*bi[j]+a_[j]*bi[i]) for j in range(len(self.fraccion))] for i in range(len(self.fraccion))]
            B2=[[-B1[i][j]-2*aij[i][j]*b**2 for j in range(len(self.fraccion))] for i in range(len(self.fraccion))]
            d=[]
            for i in range(len(self.componente)):
                d_=[]
                for j in range(len(self.componente)):
                    if i==j:
                        d_.append(1)
                    else: d_.append(0)
                d.append(d_)

            Q=[[R_atml*T*(d[i][j]/self.fraccion[i]+(bi[i]+bi[j])/(V-b)+bi[i]*bi[j]/(V-b)**2)+a*bi[i]*bi[j]/b/(V+b)**2+B1[i][j]/b**2/(V+b)+B2[i][j]/b**3*log((V+b)/V)  for j in range(len(self.componente))] for i in range(len(self.componente))]
            q=triu(Q)
            return q

        def C(T, V):
            """Subrutina de cálculo del parámetro C, eq 4B4.1-19"""
            Q=q(T, V)
            ai, bi, b, k, aij, a, a_=parameters(T)

            deltaN=[1]
            for i in range(len(self.componente)-1, 0, -1):
                deltaN.insert(0, -deltaN[0]*Q[i-1][i]/Q[i-1][i-1])
            suma=sum([n**2 for n in deltaN])**0.5
            deltaN=[n/suma for n in deltaN]

            h=[]
            for i in range(len(self.componente)):
                hi=[]
                for j in range(len(self.componente)):
                    hij=[]
                    for k in range(len(self.componente)):
                        if i==j and j==k: hij.append(1)
                        elif i!=j and i!=k and j!=k: hij.append(6)
                        else: hij.append(3)
                    hi.append(hij)
                h.append(hi)

            F=[[[b_i*b_j*b_k for b_k in bi] for b_j in bi] for b_i in bi]
            D=[[[b*a_[i]*bi[j]*bi[k]+a_[j]*bi[i]*bi[k]+a_[k]*bi[i]*bi[j]-3*F[i][j][k]*a for k in range(len(self.componente))] for j in range(len(self.componente))] for i in range(len(self.componente))]
            E=[[[2*D[i][j][k]-2*b**2*(aij[i][j]*bi[k]+aij[i][k]*bi[j]+aij[j][k]*bi[i]) for k in range(len(self.componente))] for j in range(len(self.componente))] for i in range(len(self.componente))]

            d=[]
            for i in range(len(self.componente)):
                d_=[]
                for j in range(len(self.componente)):
                    d__=[]
                    for k in range(len(self.componente)):
                        if i==j and i==k: d__.append(1)
                        else: d__.append(0)
                    d_.append(d__)
                d.append(d_)

            delta=[[[R_atml*T*(d[i][j][k]/self.fraccion[i]**2+(bi[i]*bi[j]+bi[j]*bi[k]+bi[k]*bi[i])/(V-b)**2+2*F[i][j][k]/(V-b)**3)-2*a*F[i][j][k]/b/(V+b)**3+D[i][j][k]/b**2/(V+b)**2+E[i][j][k]/b**3/(V+b)-E[i][j][k]/b**4*log((V+b)/V)  for k in range(len(self.componente))] for j in range(len(self.componente))] for i in range(len(self.componente))]
            sumaijk=sum([sum([sum([h[i][j][k]*delta[i][j][k]*deltaN[i]*deltaN[j]*deltaN[k] for k in range(len(self.componente))]) for j in range(len(self.componente))]) for i in range(len(self.componente))])

            return ((V-b)/2/b)**2*sumaijk

        To=1.5*self.tpc
        Vo=sum([fraccion*R_atml*componente.Tc/3/componente.Pc.atm for fraccion, componente in zip(self.fraccion, self.componente)])
        funcion = lambda par: (det(q(par[0], par[1])), C(par[0], par[1]))
        Tc, Vc=fsolve(funcion, [To, Vo])

        ai, bi, b, k, aij, a, a_=parameters(Tc)
        Pc=R_atml*Tc/(Vc-b)-a/(Vc*(Vc+b))
        VcCorr=sum([fraccion/(R_atml*componente.Tc/3/componente.Pc.atm-componente.Vc*componente.M) for fraccion, componente in zip(self.fraccion, self.componente)])

        return unidades.Temperature(Tc), unidades.Pressure(Pc, "atm"), unidades.SpecificVolume(VcCorr/self.M, "lg")

    # Mixing Rules
    def Mix_van_der_Waals(self, parameters, kij):
        """Miwing rules of van der Waals"""
        ai = parameters[0]
        bi = parameters[1:]
        b = [0]*len(bi)
        a = 0
        for i in range(len(self.componente)):
            for j in range(len(bi)):
                b[j] += self.fraccion[i]*bi[j][i]
            for j in range(len(self.componente)):
                a += self.fraccion[i]*self.fraccion[j]*(ai[i]*ai[j])**0.5*(1-kij[i][j])
        return tuple([a]+b)

    def Mix_Stryjek_Vera(self, parameters, kij):
        """Mixing rules of Stryjek and Vera (1986)"""
        ai = parameters[0]
        bi = parameters[1:]
        b = [0]*len(bi)
        a = 0
        for i in range(len(self.componente)):
            for j in range(len(bi)):
                b[j] += self.fraccion[i]*bi[j][i]
            for j in range(len(self.componente)):
                if kij[i][j] == 0 and kij[j][i] == 0:
                    k = 0.
                else:
                    k = kij[i][j]*kij[j][i]/(self.fraccion[i]*kij[i][j]+self.fraccion[j]*kij[j][i])
                a += self.fraccion[i]*self.fraccion[j]*(ai[i]*ai[j])**0.5*(1-k)
        return tuple([a]+b)

    def Mix_Panagiotopoulos(self, parameters, kij):
        """Mixing Rules of Panagiotopoulos (1985)"""
        ai = parameters[0]
        bi = parameters[1:]
        b = [0]*len(bi)
        a = 0
        for i in range(len(self.componente)):
            for j in range(len(bi)):
                b[j] += self.fraccion[i]*bi[j][i]
            for j in range(len(self.componente)):
                a += self.fraccion[i]*self.fraccion[j]*(ai[i]*ai[j])**0.5*(1-kij[i][j]+(kij[i][j]-kij[j][i])*self.fraccion[i])
        return tuple([a]+b)

    def Mix_Melhem(self, parameters, kij):
        """Mixing Rules of Melhem (1991)"""
        ai = parameters[0]
        bi = parameters[1:]
        b = [0]*len(bi)
        a = 0
        for i in range(len(self.componente)):
            for j in range(len(bi)):
                b[j] += self.fraccion[i]*bi[j][i]
            for j in range(len(self.componente)):
                a += self.fraccion[i]*self.fraccion[j]*(ai[i]*ai[j])**0.5*(1-kij[i][j]+(kij[i][j]-kij[j][i])*self.fraccion[i]/(self.fraccion[i]+self.fraccion[j]))
        return tuple([a]+b)

    def Lee_Kesler_Entalpia(self, T, P):
        """Método de cálculo de la entalpía haciendo uso de las propiedades críticas, método de Lee-Kesler
        Procedure API 7B3.7 Pag.643"""
        P = unidades.Pressure(P, "atm")
        H_adimensional = self.Lee_Kesler_Entalpia_lib(T, P.atm, self.Fase(T, P.atm))
        H_ideal = self._Ho(T)
        return unidades.Enthalpy(H_ideal.Jg-H_adimensional*R*self.Tc/self.M, "Jg")

    def _Ho(self, T):
        """Ideal enthalpy"""
        entalpia = 0
        for xi, cmp in zip(self.fraccion_masica, self.componente):
            entalpia += xi*cmp._Ho(T)
        return unidades.Enthalpy(entalpia)

    def Hv_DIPPR(self, T):
        entalpia = 0
        for xi, cmp in zip(self.fraccion, self.componente):
            if cmp.calor_vaporizacion[-1] < T:
                Hi = cmp.Hv_DIPPR(cmp.calor_vaporizacion[-1])
            elif cmp.calor_vaporizacion[-2] > T:
                Hi = cmp.Hv_DIPPR(cmp.calor_vaporizacion[-2])
            else:
                Hi = cmp.Hv_DIPPR(T)
            entalpia += xi*Hi
        return unidades.Enthalpy(entalpia/self.M, "Jkg")

    def Entalpia(self, T, P):
        """Método de cálculo de la entalpía, API procedure 7B4.1, pag 645"""
        Ho = self._Ho(T)
        # TODO:

    def Calor_vaporizacion(self):
        """Método de cálculo del calor de vaporización, API procedure 7C2.1, pag 682"""
        pass

    def Cp_Liquido(self, T):
        """Calculate specific heat from liquid, API procedure 7D2.1, pag 700"""
        Cp = 0
        for i, cmp in enumerate(self.componente):
            Cp += self.fraccion_masica[i]*cmp.Cp_Liquido_DIPPR(T)
        return unidades.SpecificHeat(Cp)

    def Cp_Gas(self, T, P):
        """Calculate specific heat from gas, API procedure 7D4.1, pag 714"""
        Cp = 0
        for i, cmp in enumerate(self.componente):
            Cp += self.fraccion_masica[i]*cmp.Cp_Gas_DIPPR(T)
        return unidades.SpecificHeat(Cp)
        # TODO: Añadir correción a alta presión

    def Cp_v_ratio(self):
        """Calculo de la relación de capacidades caloríficas
        Procedure API 7E2.1 Pag.728"""
        pass

    def Entropia(self, T):
        """Método de cálculo de la entropia
        Procedure API 7F2.1, Pag.741"""
        pass

    def flash_water(self, T, P):
        """Método de cálculo del equilibrio líquido-vapor en sistemas hidrocarburo-agua, API procedure 9A6.1, pag 922"""
        # TODO:
        #Volumen molar
        a = 1
        b = 1
        V = roots([1, -R_atml*T/P, -(b**2+R_atml*T/P*b-a/P), a*b/P])
        return V

    def SOUR_water(self, T):
        """Método de cálculo del equilibrio líquido-vapor en sistemas H2O-NH3-H2S, API procedure 9A7.3 pag 930, en chemcad modelo k: sour water"""

        t = Temperature(T)

        if 63 in self.ids:
            amoniaco = True
            indice_amoniaco = self.ids.index(63)
        else:
            amoniaco = False
        if 50 in self.ids:
            sulfuro = True
            indice_sulfuro = self.ids.index(50)
        else:
            sulfuro = False
        if sulfuro and amoniaco:
            ternario = True
        else:
            ternario = False

        if ternario:
            ppmnh3 = log10(self.fraccion_masica[indice_amoniaco]*1e6)
            ppmh2s = log10(self.fraccion_masica[indice_sulfuro]*1e6)
            ratio = self.fraccion[indice_amoniaco]/self.fraccion[indice_sulfuro]
            print((10**ppmh2s, ratio))

            def Z(ind, M=ratio):
                """Subrutina que nos devuelve el valor de Zm"""
                if ind == 1:
                    return 1/ratio**2
                elif ind == 2:
                    return log10(ratio)
                else:
                    return sin(ratio)/ratio

            # Table 9A7.4
            if t.R < 537.7:
                if ratio > 0.7:
                    parnh3 = [Z(1), 0.02998, 0.838848, -1.34384, -2.8802]
                else:
                    parnh3 = [Z(2), 0.13652, 0.074081, 1.12893, -3.8683]
                if ratio < 1.:
                    parh2s = [Z(3), -0.00562, 1.028528, 8.33185, -8.9175]
                elif ratio < 4.:
                    parh2s = [Z(2), 0.00227, 0.980646, -2.82752, -2.3105]
                else:
                    parh2s = [Z(2), 0.00138, 0.881842, -1.08346, -2.7605]
            elif t.R < 559.7:
                if ratio > 0.7:
                    parnh3 = [Z(1), 0.03894, 0.775221, -1.17860, -2.5972]
                else:
                    parnh3 = [Z(2), 0.06415, 0.659240, 1.24672, -4.5124]
                if ratio < 1.:
                    parh2s = [Z(3), -0.00487, 1.025371, 7.46512, -7.9480]
                elif ratio < 4.:
                    parh2s = [Z(2), 0.01018, 0.915687, -2.65953, -1.9435]
                else:
                    parh2s = [Z(2), 0.00668, 0.849756, -1.10875, -2.4133]
            elif t.R < 579.7:
                if ratio > 0.7:
                    parnh3 = [Z(1), 0.02393, 0.881656, -1.08692, -2.5668]
                else:
                    parnh3 = [Z(2), 0.09581, 0.420378, 1.30988, -3.7033]
                if ratio < 1.:
                    parh2s = [Z(3), 0.00308, 0.976347, 6.81079, -7.1893]
                elif ratio < 4.:
                    parh2s = [Z(2), 0.00562, 0.932481, -2.48535, -1.6691]
                else:
                    parh2s = [Z(2), 0.01577, 0.778382, -1.07217, -1.9946]
            elif t.R < 609.7:
                if ratio > 0.7:
                    parnh3 = [Z(1), 0.02580, 0.856844, -0.90996, -2.2738]
                else:
                    parnh3 = [Z(2), 0.07081, 0.625739, 1.23245, -3.5797]
                if ratio < 0.9:
                    parh2s = [Z(2), -0.02064, 1.090006, -0.66045, -1.0928]
                elif ratio < 2.9:
                    parh2s = [Z(2), 0.00144, 0.961073, -3.00981, -1.2917]
                else:
                    parh2s = [Z(2), 0.01306, 0.813499, -1.15162, -1.5845]
            elif t.R < 659.7:
                if ratio > 1.6:
                    parnh3 = [Z(2), 0.06548, 0.621031, 0.24106, -1.9006]
                elif ratio > 1.1:
                    parnh3 = [Z(2), 0.05838, 0.683411, 1.95035, -2.4470]
                else:
                    parnh3 = [Z(2), -0.07626, 1.462759, 2.57458, -3.5185]
                if ratio < 1.5:
                    parh2s = [Z(2), 0.00398, 0.997000, -1.95892, -1.1374]
                else:
                    parh2s = [Z(2), -0.00651, 1.037056, -1.24395, -1.4799]
            elif t.R < 679.7:
                if ratio > 1.5:
                    parnh3 = [Z(2), 0.05840, 0.647857, 0.28668, -1.8320]
                elif ratio > 1.:
                    parnh3 = [Z(2), 0.02447, 0.864904, 2.60250, -2.5186]
                else:
                    parnh3 = [Z(2), -0.08280, 1.479719, 2.26608, -3.3007]
                if ratio < 1.5:
                    parh2s = [Z(2), 0.00892, 0.957002, -1.84465, -0.9934]
                else:
                    parh2s = [Z(2), 0.00674, 0.910272, -1.16113, -1.0837]
            elif t.R < 699.7:
                if ratio > 1.6:
                    parnh3 = [Z(2), 0.04542, 0.756870, 0.18305, -1.8041]
                elif ratio > 1.:
                    parnh3 = [Z(2), 0.03659, 0.790904, 2.00739, -2.2410]
                else:
                    parnh3 = [Z(2), -0.07386, 1.445958, 1.79890, -3.1476]
                parh2s = [Z(2), 0.0, 1.002832, -1.40994, -0.9558]
            elif t.R < 719.7:
                if ratio > 1.6:
                    parnh3 = [Z(2), 0.04603, 0.728930, 0.17864, -1.5652]
                elif ratio > 1.:
                    parnh3 = [Z(2), 0.03493, 0.796292, 1.69221, -2.0369]
                else:
                    parnh3 = [Z(2), -0.07354, 1.437440, 1.73508, -2.8916]
                if ratio < 1.5:
                    parh2s = [Z(2), 0.00871, 0.964568, -1.44142, -0.7861]
                else:
                    parh2s = [Z(2), -0.01035, 1.054256, -1.22074, -0.9966]
            else:
                if ratio > 1.5:
                    parnh3 = [Z(2), 0.04484, 0.728061, 0.25132, -1.5543]
                elif ratio > 0.9:
                    parnh3 = [Z(2), 0.00193, 0.998553, 1.87924, -2.1925]
                else:
                    parnh3 = [Z(2), -0.06773, 1.391080, 1.68263, -2.6696]
                parh2s = [Z(2), 0.0, 1.000692, -1.23871, -0.7932]

            ppnh3 = 10**(parnh3[1]*ppmnh3**2+parnh3[2]*ppmnh3+parnh3[3]*parnh3[0]+parnh3[4])
            pph2s = 10**(parh2s[1]*ppmh2s**2+parh2s[2]*ppmh2s+parh2s[3]*parh2s[0]+parh2s[4])
            return Pressure(ppnh3, "mmHg"), Pressure(pph2s, "mmHg")

        else:
            if amoniaco:
                ppm = log10(self.fraccion_masica[indice_amoniaco]*1e6)
                pp = 10**(0.01129*ppm**2+0.9568*ppm+0.00719*t.R-6.83498)
            else:
                ppm = log10(self.fraccion_masica[indice_sulfuro]*1e6)
                pp = 10**(-0.00322*ppm**2+1.0318*ppm+0.00248*t.R-1.95729)

            return Pressure(pp, "mmHg")

    def SOUR_water_ph(self, T):
        """Método de cálculo del ph de disoluciones de NH3-H2S en agua, API procedure 9A8.1 pag 934"""
        t = Temperature(T)
        if 63 in self.ids:
            amoniaco = True
            indice_amoniaco = self.ids.index(63)
        else:
            amoniaco = False
        if 50 in self.ids:
            sulfuro = True
            indice_sulfuro = self.ids.index(50)
        else:
            sulfuro = False
        if sulfuro and amoniaco:
            ternario = True
        else:
            ternario = False

        if ternario:
            M = self.fraccion[indice_amoniaco]/self.fraccion[indice_sulfuro]
            pH = 10**(-0.0084*log10(log10(M))+0.1129*log10(M)-0.00032*T.R+1.0689)
        else:
            if amoniaco:
                ppm = log10(self.fraccion_masica[indice_amoniaco]*1e6)
                pH = 10**(0.480*ppm-0.0106*t.R-15.236)
            else:
                ppm = log10(self.fraccion_masica[indice_sulfuro]*1e6)
                pH = 10**(-0.505*ppm-0.0019*t.R+6.807)
        return pH

    def RhoL_Rackett(self, T):
        """Calculate saturated liquid density by Rackett method,
        procedure API 6A3.1 pag.479
        Value in mol/l"""

        # eq 6A3.1-2
        Zram = 0
        for i in range(len(self.componente)):
            Zram += self.fraccion[i]*self.componente[i].rackett

        # eq 6A3.1-5
        suma = 0
        for i in range(len(self.componente)):
            suma += self.fraccion[i]*self.componente[i].Vc
        fi = []
        for i in range(len(self.componente)):
            fi.append(self.fraccion[i]*self.componente[i].Vc/suma)

        # eq 6A3.1-7
        k = []
        for i in self.componente:
            ki = []
            for j in self.componente:
                ki.append(1-(sqrt(i.Vc**(1./3)*j.Vc**(1./3))*2 /
                             (i.Vc**(1./3)+j.Vc**(1./3)))**3)
            k.append(ki)

        # eq 6A3.1-6
        Tc = []
        for i in range(len(self.componente)):
            Tci = []
            for j in range(len(self.componente)):
                Tci.append(sqrt(self.componente[i].Tc*self.componente[j].Tc) *
                           (1-k[i][j]))
            Tc.append(Tci)

        # eq 6A3.1-4
        Tmc = 0
        for i in range(len(self.componente)):
            for j in range(len(self.componente)):
                Tmc += fi[i]*fi[j]*Tc[i][j]

        # eq 6A3.1-3
        Tr = T/Tmc

        # eq 6A3.1-1
        suma = 0
        for i in range(len(self.componente)):
            suma += self.fraccion[i]*self.componente[i].Tc/self.componente[i].Pc
        inv = R_atml*suma*Zram**(1+(1-Tr)**(2./7))
        return unidades.Density(1/inv*self.M, "gl")

    def _lib_Costald(self):
        """Library for saturated liquid density by Costald method,
        Value in mol/l"""
        # eq 6A3.1-2
        suma1 = 0
        suma2 = 0
        suma3 = 0
        for i in range(len(self.componente)):
            suma1 += self.fraccion[i]*self.componente[i].Vc
            suma2 += self.fraccion[i]*self.componente[i].Vc**(2./3)
            suma3 += self.fraccion[i]*self.componente[i].Vc**(1./3)
        Vm = (suma1+3*suma2*suma3)/4

        # eq 6A3.1-5
        suma = 0
        for i in range(len(self.componente)):
            for j in range(len(self.componente)):
                suma += self.fraccion[i]*self.fraccion[j]*sqrt(self.componente[i].Vc*self.componente[j].Vc*self.componente[i].Tc*self.componente[j].Tc)
        Tmc = suma/Vm
        return Vm, Tmc

    def RhoL_Costald(self, T):
        """Calculate saturated liquid density by Costald method,
        procedure API 6A3.2 pag.482
        Valor obtenido en mol/l"""

        Vm, Tmc = self._lib_Costald()
        Tr = T/Tmc    #eq 6A3.2-4
        Vr0 = 1-1.52816*(1-Tr)**(1./3)+1.43907*(1-Tr)**(2./3)-0.81446*(1-Tr)+0.190454*(1-Tr)**(4./3)        #eq 6A3.2-7
        Vrd = (-0.296123+0.386914*Tr-0.0427258*Tr**2-0.0480645*Tr**3)/(Tr-1.00001)    #eq 6A3.2-8

        # eq 6A3.2-1
        return unidades.Density(1/(Vm*Vr0*(1-self.f_acent_mod*Vrd))*self.M, "gl")

    def RhoL_Tait_Costald(self, T, P):
        """Calculate saturated liquid density by Tait-Costald method,
        procedure API 6A3.4 pag.489
        Presión dada en pascales
        Densidad obtenida en mol/l"""
        # FIXME: dont work
        densidad_bp = self.RhoL_Costald(T)
        Vm, Tmc = self._lib_Costald()
        Tr = T/Tmc
        zmc = 0.291-0.08*self.f_acent_mod
        pmc = zmc*R_atml*Tmc/Vm

        # eq 6A3.4-12
        alfa = 35.-36./Tr-96.736*log10(Tr)+Tr**6
        beta = log10(Tr)+0.03721754*alfa
        prm0 = 5.8031817*log10(Tr)+0.07608141*alfa
        prm1 = 4.86601*beta
        presion_bp = 10**(prm0+self.f_acent_mod*prm1)*pmc

        # eq 6A3.4-2
        C = 0.0861488+0.0344483*self.f_acent_mod
        e = exp(4.79594+0.250047*self.f_acent_mod+1.14188*self.f_acent_mod**2)
        B = pmc*(-1-9.070217*(1-Tr)**(1./3)+62.45326*(1-Tr)**(2./3)-135.1102*(1-Tr)+e*(1-Tr)**(4./3))
        return unidades.Density(densidad_bp/(1-C*log((B+P)/(B+presion_bp))))

    def RhoL_API(self, T, P):
        """Calculate liquid density, API procedure 6A3.3, pag 485"""
        Tcm = Pcm = 0
        for i in range(len(self.componente)):
            Tcm += self.fraccion[i]*self.componente[i].Tc
            Pcm += self.fraccion[i]*self.componente[i].Pc

        Pr = self.pr(P, Pcm)
        Pr0 = self.pr(1, Pcm)
        Tr = self.tr(T, Tcm)
        Tr0 = self.tr(288.71, Tcm)

        # FIXME: Add correction in gass at 1 atm and 60ºF
        suma1 = suma2 = 0
        for i in range(len(self.componente)):
            suma1 += self.fraccion[i]*self.componente[i].M
            suma2 += self.fraccion[i]*self.componente[i].M/self.componente[i].SG/1000
        rho1 = suma1/suma2

        A02 = 1.6368-0.04615*Pr+2.1138e-3*Pr**2-0.7845e-5*Pr**3-0.6923e-6*Pr**4
        A12 = -1.9693-0.21874*Pr-8.0028e-3*Pr**2-8.2328e-5*Pr**3+5.2604e-6*Pr**4
        A22 = 2.4638-0.36461*Pr-12.8763e-3*Pr**2+14.8059e-5*Pr**3-8.6895e-6*Pr**4
        A32 = -1.5841-0.25136*Pr-11.3805e-3*Pr**2+9.5672e-5*Pr**3+2.1812e-6*Pr**4
        C2 = A02+A12*Tr+A22*Tr**2+A32*Tr**3
        A01 = 1.6368-0.04615*Pr0+2.1138e-3*Pr0**2-0.7845e-5*Pr0**3-0.6923e-6*Pr0**4
        A11 = -1.9693-0.21874*Pr0-8.0028e-3*Pr0**2-8.2328e-5*Pr0**3+5.2604e-6*Pr0**4
        A21 = 2.4638-0.36461*Pr0-12.8763e-3*Pr0**2+14.8059e-5*Pr0**3-8.6895e-6*Pr0**4
        A31 = -1.5841-0.25136*Pr0-11.3805e-3*Pr0**2+9.5672e-5*Pr0**3+2.1812e-6*Pr0**4
        C1 = A01+A11*Tr0+A21*Tr0**2+A31*Tr0**3
        return unidades.Density(rho1*C2/C1)

    def Tension(self,T):
        """Calculate surface tension at low pressure,
        API procedure 10A2.1, pag 991"""
        tension = sum([xi*cmp.Tension(T) for xi, cmp in
                       zip(self.fraccion, self.componente)])
        return unidades.Tension(tension)

    def Tension_superficial_presion(self,T, parachor, fraccion_liquido, fraccion_vapor):
        """Calculate surface tension at high pressure,
        API procedure 10A2.2, pag 993
        Este procedimiento tiene un uso interno, ya que necesita como parámetros además las fracciones molares de los componentes en ambas fases, datos provenientes de un cálculo flash previo"""
        # TODO: parachor tiene que ser indicado como parámetro mientras no sepa como calcularlo para cada elemento
        mv = ml = suma = 0
        for i in range(len(self.componente)):
            mv += self.componente[i].M*fraccion_vapor[i]
            ml += self.componente[i].M*fraccion_liquido[i]
        rhoV = RhoG_Lee_Kesler(T, 1)
        suma = 0
        for i in range(len(self.componente)):
            suma += parachor[i]*(self.RhoL_Rackett(T)/ml*fraccion_liquido[i]-rhoV/mv*fraccion_vapor[i])

        return unidades.Tension(suma**4, "dyncm")


    def Tension_inferfacial_water(self, T):
        """Método de cálculo de la tensión interfacial entre agua e hidrocarburos, API procedure 10B1.3, pag 1007"""
        agua = Componente(62)
        sigma_w=agua.Tension_parametrica(T).dyncm
        sigma_h=self.Tension_superficial(T).dyncm
        return unidades.Tension(sigma_h+sigma_w-1.1*sqrt(sigma_h*sigma_w), "dyncm")


    def Mu_Liquido(self, T, P):
        """Calculate liquid viscosity, API procedure 11A3.1, pag 1051"""
        suma = 0
        for xi, cmp in zip(self.fraccion, self.componente):
            suma += xi*cmp.Mu_Liquido(T, P)**(1./3)
        return unidades.Viscosity(suma**3)

    def Mu_Gas_Stiel(self, T, P, rhoG=0, muo=0):
        """Calculate gas viscosity at high pressure, API procedure 11B4.1, pag 1107"""
        if muo == 0:
            muo = self.Mu_Gas_Wilke(T)
        if rhoG == 0:
            rhoG = P/self.Z/R_atml/T
        x = self.tpc**(1.0/6)/self.M**0.5/self.ppc**(2.0/3)
        rhor = rhoG*self.Vc/self.M
        mu = muo.cP+10.8e-5*(exp(1.439*rhor)-exp(-1.11*rhor**1.858))/x
        return unidades.Viscosity(mu, "cP")

    def Mu_Gas_Carr(self, T, P, muo=0):
        """Calculate gas viscosity ad high pressure, specify for hydrocarbon
        API procedure 11C1.2, pag 1113"""

        Tr = T/self.tpc
        Pr = P/self.ppc
        muo = self.Mu_Gas_Wilke(T)

        A1 = 83.8970*Tr**0.0105+0.6030*Tr**-0.0822+0.9017*Tr**-0.12-85.3080
        A2 = 1.514*Tr**-11.3036+0.3018*Tr**-0.6856+2.0636*Tr**-2.7611
        k = A1*1.5071*Pr**-0.4487+A2 * \
            (11.4789*Pr**0.2606-12.6843*Pr**0.1773+1.6953*Pr**-0.1052)
        return unidades.Viscosity(muo*k)

    def Mu_Gas(self, T, P, rho):
        """General method for calculate viscosity of gas"""
        if P < 2:
            Mi = [cmp.M for cmp in self.componente]
            mui = [cmp.Mu_Gas(T, 1, rho) for cmp in self.componente]
            return MuG_Wilke(self.fraccion, Mi, mui)
        else:
            return self.Mu_Gas_Stiel(T, P)
        # TODO: Add Carr method when it's availabe in database component type

    def ThCond_Liquido(self, T, P, rho):
        """Calculate liquid thermal conductivity, API procedure 12A2.1, pag 1145"""
        ki = []
        for cmp in self.componente:
            ki.append(cmp.ThCond_Liquido(T, P, rho))
        kij = []
        for i in range(len(self.componente)):
            kiji = []
            for j in range(len(self.componente)):
                kiji.append(2/(1/ki[i]+1/ki[j]))
            kij.append(kiji)
        xV = 0
        Vi = []
        for i in range(len(self.componente)):
            Vi.append(self.componente[i].M/self.componente[i].RhoL(T, P))
            xV += self.fraccion[i]*Vi[i]
        fi_bruto = []
        suma = 0
        for i in range(len(self.componente)):
            fi_bruto.append(self.fraccion[i]*Vi[i]/xV)
            suma += fi_bruto[i]
        fi = []
        for i in range(len(self.componente)):
            fi.append(fi_bruto[i]/suma)
        k = 0
        for i in range(len(self.componente)):
            for j in range(len(self.componente)):
                k += fi[i]*fi[j]*kij[i][j]
        return unidades.ThermalConductivity(k)

    def ThCond_Gas(self, T, P, rho):
        """Calculate gas thermal conductivity, API procedure 12A2.1, pag 1145"""
        ki = []
        S = []
        mu = []
        for cmp in self.componente:
            ki.append(cmp.ThCond_Gas(T, P, rho))
            if cmp.indice == 1:
                S.append(78.8888889)
            else:
                S.append(1.5*cmp.Tb)
            mu.append(cmp.Mu_Gas(T, P, rho))

        Sij = []
        for i in range(len(self.componente)):
            Siji = []
            for j in range(len(self.componente)):
                Siji.append(sqrt(S[i]*S[j]))
            Sij.append(Siji)

        Aij = []
        for i in range(len(self.componente)):
            Aiji = []
            for j in range(len(self.componente)):
                Aiji.append(0.25*(1+sqrt(mu[i]/mu[j]*(self.componente[j].M/self.componente[i].M)**0.75*(1+S[i]/T)/(1+S[j]/T)))**2*(1+Sij[i][j]/T)/(1+S[i]/T))
            Aij.append(Aiji)

        sumaj = []
        for i in range(len(self.componente)):
            suma = 0
            for j in range(len(self.componente)):
                suma += Aij[i][j]*self.fraccion[j]
            sumaj.append(suma)

        k = 0
        for i in range(len(self.componente)):
            k += ki[i]*self.fraccion[i]/sumaj[i]

        return unidades.ThermalConductivity(k)

    def Solubilidad_agua(self, T):
        """Método de cálculo de la solubilidad de agua en el componente, API procedure 9A1.1, 9A1.5, pag 897
        Solubilidad obtenida expresada en fracción molar"""
        # TODO: Rehacer para la última versión de API
        t = unidades.Temperature(T)
        modo1 = [4, 6, 5, 8, 7, 10, 55, 11, 541, 12, 82, 14, 38, 60, 23, 27,
                 35, 28, 40, 41, 45, 42, 43, 71, 861, 885, 178]
        if self.indice in modo1:
            a1 = [-6.158, -0.792, 10.289, -6.46, -6.361, -3.441, 11.511,
                  -0.082, 12.831, 0.386, 1.934, 0.998, -3.615, 0.579, 16.123,
                  3.489, 1.202, 2.041, -0.595, -1.271, -0.587, 0.382, -1.024,
                  20.439, 0.879, 0.905, -21.475]
            a2 = [-964, -2540, -6054, -572, -1007, -1306, -5918, -2428, -6281,
                  -2536, -2939, -2689, -1263, -2569, -7437, -3270, -2661,
                  -2489, -1591, -1386, -1631, -1896, -1705, -7147, -2192,
                  -2209, 3902]
            a3 = [0.8338e-2, 0.3597e-2, -0.464e-2, 0.7548e-2, 0.9174e-2,
                  0.4815e-2, -0.679e-2, 0.2379e-2, -0.789e-2, 0.1958e-2,
                  0.076e-2, 0.1416e-2, 0.4778e-2, 0.1654e-2, -1.052e-2,
                  -0.073e-2, 0.1262e-2, -0.017e-2, 0.198e-2, 0.2451e-2,
                  0.1955e-2, 0.1123e-2, 0.2584e-2, -1.8040e-2, 0.1007e-2,
                  0.1010e-2, 2.0184e-2]
            indice = modo1.index(self.indice)
            return 10**(a1[indice]+a2[indice]/t.R+a3[indice]*t.R)
        else:
            H = self.composicion_molecular[1][self.composicion_molecular[0].index("H")]*Elemental(1).atomic_mass
            C = self.composicion_molecular[1][self.composicion_molecular[0].index("C")]*Elemental(6).atomic_mass
            return 10**(-(4200*H/C+1050)*(1.8/t.R-0.0016))

    def Solubilidad_en_agua(self, T):
        """Método de cálculo de la solubilidad en agua del componente en condiciones de equilibrio liquido-vapor, API procedure 9A2.1, 9A2.3, pag 899, método renovado en API v7
        Solubilidad obtenida expresada en fracción molar"""
        t = unidades.Temperature(T)
        modo1 = [4, 6, 8, 10, 11, 12, 13, 36, 38, 39, 60, 389, 24, 35, 83, 40,
                 41, 45, 70, 377, 1486, 886, 71, 77, 862, 885, 178, 191, 1497]
        modo2 = [185, 193, 194, 197, 409, 200, 203, 204, 206, ]
        if self.indice in modo1:
            a1 = [-314.3825, -292.4658, -361.7959, -406.5873, -430.4197,
                  -450.7616, -469.842, -235.6808, -326.816, -356.3782,
                  -385.8427, -444.1288, -228.836, -291.413, -366.146,
                  -195.5673, -218.6392, -241.8211, -330.4825, -375.6332,
                  -420.7929, -465.8706, -262.5274, -248.0811, -311.3553,
                  -308.4157, -288.0113, -220.5168, -247.3787]
            a2 = [24225.156, 21452.23, 26167.45, 29388.83, 31018.136,
                  32355.695, 33782.08, 16843.48, 23264.01, 25331.92, 27399.83,
                  31535.66, 16995.95, 20436.66, 25673.53, 13544.694, 15119.661,
                  16694.626, 22994.478, 26144.406, 29294.352, 32444.28,
                  17838.556, 16349.596, 21748.495, 20014.725, 19278.918,
                  12603.309, 15072.32]
            a3 = [41.50287, 38.58148, 47.97436, 53.89582, 56.95927, 59.55451,
                  61.94, 30.89572, 43.298, 47.1467, 50.9954, 58.6928, 30.0380,
                  38.4871, 48.3494, 25.8585, 28.8653, 31.8721, 43.8994, 49.913,
                  55.9266, 61.9402, 34.67748, 32.77122, 41.11214, 40.74593,
                  38.54914, 29.35657, 32.71111]
            indice = modo1.index(self.indice)
            return exp(a1[indice]+a2[indice]/t.R+a3[indice]*log(t.R))
        elif self.indice in modo2:
            a1 = [-373.24613, -334.02054, -332.69462, -380.76757, -420.02594,
                  -350.23574, -0.84491, -380.93125, -7.53314]
            a2 = [21618.93, 17740.54, 17633.99, 20361.3, 21096.86, 17490.54,
                  -9069.58, 19440.39, -8120.29]
            a3 = [51.01156, 45.56672, 45.48078, 52.08392, 5753568, 47.9872, 0,
                  51.93974, 0]
            indice = modo2.index(self.indice)
            return exp(a1[indice]+a2[indice]/t.R+a3[indice]*log(t.R))

    def Solubilidad_en_agua_25(self):
        """Método de cálculo de la solubilidad en agua saturada a 25ºC, API procedure 9A2.6, pag 909"""
        #TODO: Díficil de implementar mientras no se añada a la base de datos alguna propiedad que indique la naturaleza química del compuesto

    def Solubilidad_Henry(self,T, P):
        """Solubilidad de gases en líquidos usando la ley de Henry
        Temperatura dada en ºR
        constante H obtenida en psia por unidad de fracción molar del gas
        lnH = A/T + B∗lnT + C∗T + D
        Solo disponible para algunos compuestos:
        Hydrogen, Helium, Argon, Neon, Krypton, Xenon, Oxygen, Nitrogen,
        Hydrogen sulfide, Carbon monoxide, Carbon dioxide, Sulfur dioxide,
        Nitrous oxide, Chlorine,Bromine, Iodine, Methane, Ethane, Propane,
        Ethylene, Ammonia.
        API procedure 9A7.1, pag 927
        Los parametros de la ecuación se encuentran en la base de datos
        en forma de lista en la posición décima
        Solubilidad obtenida en fracción molar"""
        t = unidades.Temperature(T)
        p = unidades.Pressure(P, "atm")
        H = exp(self.henry[0]/T.R+self.henry[1]*log(T.R)+self.henry[2]*T.R+self.henry[3])
        return p.psi/H

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
            self.fraccion = [unidades.Dimensionless(x) for x in mezcla["fraction"]]
            self.fraccion_masica = [unidades.Dimensionless(x) for x in mezcla["massFraction"]]
            self.caudalunitariomasico = [unidades.MassFlow(x) for x in mezcla["massUnitFlow"]]
            self.caudalunitariomolar = [unidades.MolarFlow(x) for x in mezcla["molarUnitFlow"]]
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


if __name__ == '__main__':

#    aire=Mezcla([46, 47, 49], [0.781, 0.209, 0.01])
#    print aire.RhoL_Rackett(20)
#    print aire.RhoL_Costald(20)[0]
#    print 1/aire.RhoG_Lee_Kesler(300, 1),  "l/mol"
#    print aire.RhoL_API(20, 1)
#    print aire.Viscosidad_liquido(20)
#    suma=0
#    for i in aire.fraccion_masica:
#        print i*0.9
#        suma+=i
#    print suma


#    ejemplo1=Mezcla([3, 11], [0.5871, 0.4129])
#    print Densidad(ejemplo1.densidad_Rackett(Temperature(91, "F"))).lbft3
#
#    ejemplo2=Mezcla([2, 14], [0.2, 0.8])
#    print Densidad(ejemplo2.densidad_Costald(Temperature(160, "F"))[0]).lbft3
#
#    ejemplo3=Mezcla([3, 14], [0.2, 0.8])
#    print Densidad(ejemplo3.densidad_Tait_Costald(Temperature(160, "F"), Pressure(3000, "psi").atm), "gl").lbft3

#    ejemplo5=Mezcla([1, 25, 49], [0.3, 0.415, 0.285], 300)
#    print Temperature(ejemplo5.Tc()).R
#    print ejemplo5.propiedades_criticas()
#    print ejemplo5.kij()


#    ejemplo6=Mezcla([45, 12], [0.3, 0.7])
#    print Pressure(ejemplo6.Pc(), "atm").psi

#    ejemplo7=Mezcla([6, 11], [0.63, 0.37])
#    print Volumen(ejemplo7.Vc()).ft3*lb

#    ejemplo8=Mezcla([2, 3], [0.1, 0.9])
#    print Temperature(ejemplo8.Tc()).R
#    print Pressure(ejemplo8.Pc(), "atm").psi
#    print Volumen(ejemplo8.Vc()).ft3*lb
#    print ejemplo8.kij()

#    ejemplo9=Mezcla([2, 3, 4, 6, 46], [0.868, 0.065, 0.025, 0.007, 0.034])
#    print ejemplo9.Tb().F
#    print ejemplo9.tc_gas().F

#    ejemplo10=Mezcla([133, 7], [7.56, 92.44])
#    print ejemplo10.Reid()

#    ejemplo11=Mezcla([3, 14], [0.2, 0.8])
#    t=Temperature(160, "F")
#    p=Pressure(3000, "psi")
#    print ejemplo11.RhoL_API(t, p.atm)

#    ejemplo12=Mezcla([3, 14], [0.2, 0.8])
#    print ejemplo12.flash_water(300, 1)

#    amoniaco=Mezcla([62, 63], [0.9894, 0.0106])
#    t=unidades.Temperature(248, "F")
#    print amoniaco.fraccion_masica[1]*1e6
#    print amoniaco.SOUR_water(t).psi

#    sour=Mezcla([62, 63, 50], [0.9796, 0.0173, 0.0031])
#    print sour.SOUR_water(t)[1].mmHg

#    sour=Mezcla([62, 63, 50], [0.97, 0.02, 0.01])
#    t=unidades.Temperature(80, "F")
#    print sour.SOUR_water_ph(t)


#    ejemplo13=Mezcla([40, 38], [0.379, 0.621])
#    print ejemplo13.Tension_superficial(unidades.Temperature(77, "F")).dyncm

#    ejemplo14=Mezcla([2, 4], [0, 1])
#    print ejemplo14.Tension_superficial_presion([72.6, 150.8], [0.418, 0.582], [0.788, 0.212])

#    ejemplo15=Mezcla([14], [1])
#    print ejemplo15.Tension_inferfacial_water(100+273.15).dyncm

#    ejemplo16=Mezcla([20, 40, 10], [0.2957, 0.3586, 0.3457])
#    print ejemplo16.Viscosidad_liquido(unidades.Temperature(77, "F")).cP

#    ejemplo17=Mezcla([4, 8, 38], [0.25, 0.5, 0.25])
#    print ejemplo17.Viscosidad_liquido(unidades.Temperature(160, "F")).cP
#    print ejemplo17.Viscosidad_gas(unidades.Temperature(160, "F"))

#    ejemplo18=Mezcla([1, 4], [0.5818, 0.4182])
#    t=unidades.Temperature(77, "F")
#    print ejemplo18.Mu_Gas_Wilke(t).cP

#    ejemplo=Mezcla([2, 3, 4, 46], [0.956, 0.036, 0.005, 0.003])
#    t=unidades.Temperature(85, "F")
#    print ejemplo.Mu_Gas_Wilke(t).cP

#    ejemplo19=Mezcla([2, 3, 4, 46], [95.6, 3.6, 0.5, 0.3])
#    print ejemplo19.Viscosidad_gas(unidades.Temperature(85, "F")).cP
#    print ejemplo19.Viscosidad_gas_presion(unidades.Temperature(85, "F"), 5).cP

#    ejemplo20=Mezcla([2, 4], [0.6, 0.4])
#    print ejemplo20.Viscosidad_gas(398.15).cP
#    print ejemplo20.Viscosidad_gas_presion(398.15, 102.07).cP
#    print ejemplo20.MuG_Carr(398.15, 102.07).cP

#    ejemplo21=Mezcla([11, 36], [0.68, 0.32])
#    t=unidades.Temperature(32, "F")
#    print ejemplo21.Conductividad_Termica_liquido(t, 1).BtuhftF

#    ejemplo22=Mezcla([8, 10], [0.2996, 0.7004])
#    t=unidades.Temperature(212, "F")
#    print ejemplo22.Conductividad_Termica_gas(t, 1).BtuhftF

#    ejemplo23=Mezcla([5, 10], [0.45, 0.55])
#    t=unidades.Temperature(200, "C")
#    print ejemplo23.RhoG_SRK(t, 1)
#    print ejemplo23.Z_SRK(t, 1)
#    print ejemplo23.Entalpia_SRK(t, 1)
#    print ejemplo23.RhoG_MSRK(t, 1)
#    print ejemplo23.Z_MSRK(t, 1)
#    print ejemplo23.Z_BWRS(t, 1)

#    ejemplo24=Mezcla([49, 4], [0.5, 0.5])
#    p=unidades.Pressure(140, "bar")
#    print ejemplo24.Z_RK(450, p.atm)
#    print ejemplo24.Z_SRK(450, p.atm)
#    print ejemplo24.Z_PR(450, p.atm)


#    ejemplo25=Mezcla([2, 3, 4, 50, 49, 46], [0.1, 0.3, 0.25, 0.05, 0.23, 0.07])
#    print ejemplo25.Z_PR(298.15, 1)
#    print ejemplo25.RhoG_PR(298.15, 1)
#    print ejemplo25.Fugacidad_PR(298.15, 1)
#    print ejemplo25.Mu_Gas_Wilke(298.15)
#    print ejemplo25.ThCond_Gas(298.15, 1)

#    ejemplo26=Mezcla([10, 38, 22, 61], [.3, 0.5, 0.05, 0.15])
#    t=unidades.Temperature(76, "C")
#    ejemplo26.Flash_SRK(320, 1)

#    ejemplo27=Mezcla( [46, 47], [0.781, 0.209])
#    t=unidades.Temperature(76, "C")
#    print ejemplo27.Flash_SRK(t, 1)

#    ejemplo28=Mezcla([8, 12], [39.2, 60.8]) #Ej pag 647
#    t=unidades.Temperature(590, "F")
#    p=unidades.Pressure(1400, "psi")
#    print ejemplo28.Lee_Kesler_Entalpia(t, p.atm)

#    t=unidades.Temperature(160, "F")
#    p=unidades.Pressure(3000, "psi")
#    mezcla=Mezcla([2, 14], [0.2, 0.8]) #Ej pag 482
#    print mezcla.RhoL_Costald(t).gml
#
#    mezcla=Mezcla([3, 14], [0.2, 0.8]) #Ej pag 490
#    print mezcla.RhoL_Tait_Costald(t, p.atm).gml

#    agua=Mezcla( [62], [1])
#    t=unidades.Temperature(76, "C")
#    print agua.Flash_SRK(t, 1)

#    ejemplo=Corriente(340, 1, 1, Mezcla([10, 38, 22, 61], [.3, 0.5, 0.05, 0.15]))
#    print  ejemplo.Entalpia.MJh

#    agua1=Corriente(298.15, 1, 1, Mezcla([62], [1.0]))
#    agua2=Corriente(350, 1, 1, Mezcla([62], [1.0]))
#    print agua2.Entalpia.MJkg-agua1.Entalpia.MJkg

#    ejemplo2=Mezcla([10, 38, 22, 61], [0.3045265648865772, 0.61740493041538613, 0.0010833278211912348, 0.076985176876845404])
#    print ejemplo2.RhoL_Tait_Costald(340, 1)

#    ejemplo=Mezcla([10, 38, 22, 61], [.3, 0.5, 0.05, 0.15])
#    print ejemplo.SRK_lib(300)

#    ejemplo=Mezcla([2, 4], [0.6, 0.4])
#    t=unidades.Temperature(257, "F")
#    p=unidades.Pressure(1500, "psi")
#    print ejemplo.Mu_Gas_Stiel(t, p.atm)

#    ejemplo=Mezcla([11, 36], [0.68, 0.32])
#    t=unidades.Temperature(32, "F")
#    print ejemplo.ThCond_Liquido(t, 1).BtuhftF

#    ejemplo=Mezcla([10, 38, 22, 61], [.3, 0.5, 0.05, 0.15])
#    t=unidades.Temperature(32, "F")
#    print ejemplo.ThCond_Liquido(t, 1).BtuhftF

#    ejemplo=Mezcla([40, 38], [0.552, 0.448]) #Ej. pag 701
#    t=unidades.Temperature(72.1, "F")
#    print ejemplo.Cp_liquido(t).BtulbF
#
#
#    ejemplo=Mezcla([6, 11], [0.63, 0.37])
#    print ejemplo.Vc.ft3lb/ejemplo.M

#    ejemplo=Mezcla([2, 3], [0.1, 0.9])
#    print ejemplo.Tc.R, ejemplo.Pc.psi, ejemplo.Vc.ft3lb/ejemplo.M


#    z=0.965
#    mez=Mezcla(tipo=3, fraccionMolar=[z, 1-z], caudalMasico=1.)
#    tb=mez.componente[0].Tb
#    print tb
#    corr=Corriente(T=tb, P=101325., mezcla=mez)
#    print corr.eos._Dew_T()

    mezcla=Mezcla(caudalMasico=0.01, ids=[10, 38, 22, 61], fraccionMolar=[.3, 0.5, 0.05, 0.15], tipo=3)
    corr=Corriente(T=300, P=101325., mezcla=mezcla)
    print(corr.x, corr.Q, corr.P, corr.caudalmasico)
