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


This module implement pure component properties

:class:`Componente`: The main class with all integrated functionality. Use the
properties in database and calculate state properties with the methods chossen
in configuration

Liquid density calculation methods:
    * :func:`DIPPR`
    * :func:`RhoL_Rackett`
    * :func:`RhoL_Cavett`
    * :func:`RhoL_Costald`
    * :func:`RhoL_YenWoods`
    * :func:`RhoL_YamadaGunn`
    * :func:`RhoL_Bhirud`
    * :func:`RhoL_Mchaweh`
    * :func:`RhoL_Riedel`
    * :func:`RhoL_ChuehPrausnitz`
    * :func:`RhoL_TaitCostald`
    * :func:`RhoL_ChangZhao`
    * :func:`RhoL_AaltoKeskinen`
    * :func:`RhoL_AaltoKeskinen2`
    * :func:`RhoL_Nasrifar`
    * :func:`RhoL_API`

Liquid viscosity calculation methods:
    * :func:`DIPPR`
    * :func:`MuL_Parametric`
    * :func:`MuL_LetsouStiel`
    * :func:`MuL_PrzedzieckiSridhar`
    * :func:`MuL_Lucas`
    * :func:`MuL_API`
    * :func:`MuL_Kouzel`

Gas viscosity calculation methods:
    * :func:`DIPPR`
    * :func:`MuG_ChapmanEnskog`
    * :func:`MuG_Chung`
    * :func:`MuG_Lucas`
    * :func:`MuG_StielThodos`
    * :func:`MuG_Gharagheizi`
    * :func:`MuG_YoonThodos`
    * :func:`MuG_P_Chung`
    * :func:`MuG_Brule`
    * :func:`MuG_Jossi`
    * :func:`MuG_TRAPP`
    * :func:`MuG_P_StielThodos`
    * :func:`MuG_Reichenberg`
    * :func:`MuG_DeanStiel`
    * :func:`MuG_API`

Gas thermal conductivity calculation methods:
    * :func:`DIPPR`
    * :func:`ThG_MisicThodos`
    * :func:`ThG_Chung`
    * :func:`ThG_Eucken`
    * :func:`ThG_EuckenMod`
    * :func:`ThG_RiaziFaghri`
    * :func:`ThG_StielThodos`
    * :func:`ThG_P_Chung`
    * :func:`ThG_TRAPP`

Liquid thermal conductivity calculation methods:
    * :func:`DIPPR`
    * :func:`ThL_Pachaiyappan`
    * :func:`ThL_SatoRiedel`
    * :func:`ThL_KanitkarThodos`
    * :func:`ThL_RiaziFaghri`
    * :func:`ThL_Gharagheizi`
    * :func:`ThL_LakshmiPrasad`
    * :func:`ThL_Nicola`
    * :func:`ThL_Lenoir`
    * :func:`ThL_Missenard`

Vapor pressure calculation methods:
    * :func:`DIPPR`
    * :func:`Pv_Wagner`
    * :func:`Pv_Antoine`
    * :func:`Pv_AmbroseWalton`
    * :func:`Pv_Lee_Kesler`
    * :func:`Pv_Riedel`
    * :func:`Pv_Sanjari`
    * :func:`Pv_MaxwellBonnel`

Surface tension calculation methods:
    * :func:`DIPPR`
    * :func:`Tension_Parametric`
    * :func:`Tension_BlockBird`
    * :func:`Tension_Pitzer`
    * :func:`Tension_ZuoStenby`
    * :func:`Tension_SastriRao`
    * :func:`Tension_Hakim`
    * :func:`Tension_Miqueu`

Acentric factor calculation methods:
    * :func:`facent_LeeKesler`
    * :func:`prop_Edmister`
    * :func:`facent_AmbroseWalton`

Others method:
    * :func:`Vc_Riedel`
    * :func:`Henry`
    * :func:`atomic_decomposition`
    * :func:`refrigerantCode`
'''


import math
import os
import re

from numpy.lib.scimath import log, log10
from numpy import exp, cosh, sinh, tanh, roots, absolute, array
from scipy.optimize import fsolve
from scipy.constants import R, Avogadro, Boltzmann
from scipy.interpolate import interp1d, interp2d

from lib.physics import R_atml, Collision_Neufeld
from lib import unidades, config, sql
from lib.utilities import refDoc


__doi__ = {
    1:
        {"autor": "Poling, B.E, Prausnitz, J.M, O'Connell, J.P",
         "title": "The Properties of Gases and Liquids 5th Edition",
         "ref": "McGraw-Hill, New York, 2001",
         "doi": ""},
    2:
        {"autor": "Ahmed, T.",
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
        {"autor": "Lee, B.I., Kesler, M.G.",
         "title": "A Generalized Thermodynamic Correlation Based on "
                  "Three-Parameter Corresponding States",
         "ref": "AIChE Journal 21(3) (1975) 510-527",
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
         "doi": "10.1016/0300-9467(92)80037-b"},
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
        {"autor": "Riedel, L.",
         "title": "Kritischer Koeffizient, Dichte des gesättigten Dampfes und "
                  "Verdampfungswärme: Untersuchungen über eine Erweiterung des"
                  " Theorems der übereinstimmenden Zustände. Teil III",
         "ref": "Chem. Ingr. Tech., 26(12) (1954) 679-683",
         "doi": "10.1002/cite.330261208"},
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
        {"autor": "Stiel, L.I., Thodos, G.",
         "title": "The Viscosity of Nonpolar Gases at Normal Pressures",
         "ref": "AIChE J. 7(4) (1961) 611-615",
         "doi": "10.1002/aic.690070416"},
    25:
        {"autor": "Misic, D., Thodos, G.",
         "title": "The Thermal Conductivity of Hydrocarbon Gases at Normal "
                  "Pressures",
         "ref": "AIChE Journal 7(2) (1961) 264-267",
         "doi": "10.1002/aic.690070219"},
    26:
        {"autor": "Yen, L.C., Woods, S.S.",
         "title": "A Generalized Equation for Computer Calculation of Liquid "
                  "Densities",
         "ref": "AIChE Journal 12(1) (1966) 95-99",
         "doi": "10.1002/aic.690120119"},
    27:
        {"autor": "Gunn, R.D., Yamada, T.",
         "title": "A Corresponding States Correlation of Saturated Liquid "
                  "Volumes",
         "ref": "AIChE Journal 17(6) (1971) 1341-1345",
         "doi": "10.1002/aic.690170613"},
    28:
        {"autor": "Bhirud, V.L.",
         "title": "Saturated Liquid Densities of Normal Fluids",
         "ref": "AIChE Journal 24(6) (1978) 1127-1131",
         "doi": "10.1002/aic.690240630"},
    29:
        {"autor": "Mchaweh, A., Alsaygh, A., Nasrifar, Kh., Moshfeghian, M.",
         "title": "A Simplified Method for Calculating Saturated Liquid "
                  "Densities",
         "ref": "Fluid Phase Equilibria 224 (2004) 157-167",
         "doi": "10.1016/j.fluid.2004.06.054"},
    30:
        {"autor": "Chang, C.H., Zhao, X.M.",
         "title": "A New Generalized Equation for Predicting Volume of "
                  "Compressed Liquids",
         "ref": "Fluid Phase Equilibria, 58 (1990) 231-238",
         "doi": "10.1016/0378-3812(90)85134-v"},
    31:
        {"autor": "Aalto, M., Keskinen, K.I., Aittamaa, J., Liukkonen, S.",
         "title": "An Improved Correlation for Compressed Liquid Densities of "
                  "Hydrocarbons. Part 1. Pure Compounds",
         "ref": "Fluid Phase Equilibria 114 (1996) 1-19",
         "doi": "10.1016/0378-3812(95)02822-6"},
    32:
        {"autor": "Riedel, L.",
         "title": "Die Flüssigkeitsdichte im Sättigungszustand. Untersuchungen"
                  " über eine Erweiterung des Theorems der übereinstimmenden "
                  "Zustände. Teil II.",
         "ref": "Chem. Eng. Tech. 26(5) (1954) 259-264",
         "doi": "10.1002/cite.330260504"},
    33:
        {"autor": "Rea, H.E., Spencer, C.F., Danner, R.P.",
         "title": "Effect of Pressure and Temperature on the Liquid Densities "
                  "of Pure Hydrocarbons",
         "ref": "J. Chem. Eng. Data 18(2) (1973) 227-230",
         "doi": "10.1021/je60057a003"},
    34:
        {"autor": "Chueh, P.L., Prausnitz, J.M.",
         "title": "Vapor-Liquid Equilibria at High Pressures: Calculation of "
                  "Partial Molar Volumes in Nonpolar Liquid Mixtures",
         "ref": "AIChE Journal 13(6) (1967) 1099-1107",
         "doi": "10.1002/aic.690130612"},
    35:
        {"autor": "Spencer, C.F., Danner, R.P.",
         "title": "Improved Equation for Prediction of Saturated Liquid "
                  "Density",
         "ref": "J. Chem. Eng. Data 17(2) (1972) 236-241",
         "doi": "10.1021/je60053a012"},
    36:
        {"autor": "Nasrifar, K., Ayatollahi, S., Moshfeghian, M.",
         "title": "A Compressed Liquid Density Correlation",
         "ref": "Fluid Phase Equilibria 168 (2000) 149-163",
         "doi": "10.1016/s0378-3812(99)00336-2"},
    37:
        {"autor": "Aalto, M., Keskinen, K.I.",
         "title": "Liquid Densities at High Pressures",
         "ref": "Fluid Phase Equilibria 166 (1999) 183-205",
         "doi": "10.1016/s0378-3812(99)00300-3"},
    38:
        {"autor": "Pal, A., Pope, G., Arai, Y., Carnahan, N., Kobayashi, R.",
         "title": "Experimental Pressure-Volume-Temperature Relations for "
                  "Saturated and Compressed Fluid Ethane",
         "ref": "J. Chem. Eng. Data 21(4) (1976) 394-397",
         "doi": "10.1021@je60071a008"},
    39:
        {"autor": "Brock, J.R., Bird, R.B.",
         "title": "Surface Tension and the Principle of Corresponding States",
         "ref": "AIChE Journal 1(2) (1955) 174-177",
         "doi": "10.1002/aic.690010208"},
    40:
        {"autor": "Miller, D.G., Thodos, G.",
         "title": "Reduced Frost-Kalkwarf Vapor Pressure Equation",
         "ref": "I&EC Fundamentals 2(1) (1963) 78-80",
         "doi": "10.1021/i160005a015"},
    41:
        {"autor": "Zuo, Y., Stenby, E.H.",
         "title": "Corresponding-States and Parachor Models for the "
                  "Calculation of Interfacial Tensions",
         "ref": "Can. J. Chem. Eng. 75(6) (1997) 1130-1137",
         "doi": "10.1002/cjce.5450750617"},
    42:
        {"autor": "Sastri, S.R.S., Rao, K.K.",
         "title": "A Simple Method to Predict Surface Tension of Organic "
                  "Liquids",
         "ref": "Chem. Eng. Journal 59(2) (1995) 181-186",
         "doi": "10.1016/0923-0467(94)02946-6"},
    43:
        {"autor": "Hakim, D.I., Steinberg, D., Stiel, L.I.",
         "title": "Generalized Relationship for the Surface Tension of Polar "
                  "Fluids",
         "ref": "I&EC Fundamentals 10(1) (1971) 174-75.",
         "doi": "10.1021/i160037a032"},
    44:
        {"autor": "Miqueu, C., Broseta, D., Satherley, J., Mendiboure, B., "
                  "Lachaise, J., Graciaa, A.",
         "title": "An Extended Scaled Equation for the Temperature Dependence "
                  "of the Surface Tension of Pure Compounds Inferred from an "
                  "Analysis of Experimental Data",
         "ref": "Fluid Phase Equilibria 172(2) (2000) 169-182",
         "doi": "10.1016/s0378-3812(00)00384-8"},
    45:
        {"autor": "Przedziecki, J.W., Sridhar, T.",
         "title": "Prediction of Liquid Viscosities",
         "ref": "AIChE Journal 31(2) (1985) 333-335",
         "doi": "10.1002/aic.690310225"},
    46:
        {"autor": "Lucas, K.",
         "title": "Die Druckabhängigheit der Viskosität von Flüssigkeiten, "
                  "eine Einfache Abschätzung",
         "ref": "Chem. Ing. Tech. 46(4) (1981) 959-960",
         "doi": "10.1002/cite.330531209"},
    47:
        {"autor": "Riazi, M. R.",
         "title": "Characterization and Properties of Petroleum Fractions.",
         "ref": "ASTM manual series MNL50, 2005",
         "doi": ""},
    48:
        {"autor": "Brokaw, R.S.",
         "title": "Predicting Transport Properties of Dilute Gases",
         "ref": "I&EC Process Design and Development 8(22) (1969) 240-253",
         "doi": "10.1021/i260030a015"},
    49:
        {"autor": "Chung, T.H., Ajlan, M., Lee, L.L., Starling, K.E.",
         "title": "Generalized Multiparameter Correlation for Nonpolar and "
                  "Polar Fluid Transport Properties",
         "ref": "Ind. Eng. Chem. Res. 27(4) (1988) 671-679",
         "doi": "10.1021/ie00076a024"},
    50:
        {"autor": "Chung, T.H., Lee, L.L., Starling, K.E.",
         "title": "Applications of Kinetic Gas Theories and Multiparameter "
                  "Correlation for Prediction of Dilute Gas Viscosity and "
                  "Thermal Conductivity",
         "ref": "Ind. Eng. Chem. Fundam. 23(1) (1984) 8-13",
         "doi": "10.1021/i100013a002"},
    51:
        {"autor": "Jossi, J.A., Stiel, L.I., Thodos, G.",
         "title": "The Viscosity of Pure Substances in the Dense Gaseous and "
                  "Liquid Phases",
         "ref": "AIChE Journal 8(1) (1962) 59-63",
         "doi": "10.1002/aic.690080116"},
    52:
        {"autor": "Stiel, L.I., Thodos, G.",
         "title": "The Viscosity of Polar Substances in the Dense Gaseous and "
                  "Liquid Regions",
         "ref": "AIChE Journal 10(2) (1964) 275-277",
         "doi": "10.1002/aic.690100229"},
    53:
        {"autor": "Younglove, B.A., Ely, J.F.",
         "title": "Thermophysical Properties of Fluids. II. Methane, Ethane, "
                  "Propane, Isobutane, and Normal Butane",
         "ref": "J. Phys. Chem. Ref. Data 16(4) (1987) 577-798",
         "doi": "10.1063/1.555785"},
    54:
        {"autor": "Brulé, M.R., Starling, K.E.",
         "title": "Thermophysical Properties of Complex Systems: Applications "
                  "of Multiproperty Analysis",
         "ref": "Ind. Eng. Chem. Process Dev. 23 (1984) 833-845",
         "doi": "10.1021/i200027a035"},
    55:
        {"autor": "Dean, D.E., Stiel, L.I.",
         "title": "The Viscosity of Nonpolar Gas Mixtures at Moderate and High"
                  " Pressures",
         "ref": "AIChE Journal 11(3) (1965) 526-532 ",
         "doi": "10.1002/aic.690110330"},
    56:
        {"autor": "Yoon, P., Thodos, G.",
         "title": "Viscosity of Nonpolar Gaseous Mixtures at Normal Pressures",
         "ref": "AIChE Journal 16(2) (1970) 300-304",
         "doi": "10.1002/aic.690160225"},
    57:
        {"autor": "Gharagheizi, F., Eslamimanesh, A., Sattari, M., Mohammadi, "
                  "A.H., Richon, D.",
         "title": "Corresponding States Method for Determination of the "
                  "Viscosity of Gases at Atmospheric Pressure",
         "ref": "I&EC Research 51(7) (2012) 3179-3185",
         "doi": "10.1021/ie202591f"},
    58:
        {"autor": "Stiel, L.I., Thodos, G.",
         "title": "The Thermal Conductivity of Nonpolar Substances in the "
                  "Dense Gaseous and Liquid Regions",
         "ref": "AIChE Journal 10(1) (1964) 26-30",
         "doi": "10.1002/aic.690100114"},
    59:
        {"autor": "Lenoir, J.M.",
         "title": "Effect of Pressure on Thermal Conductivity of Liquids",
         "ref": "Petroelum Refiner 36(8) 1508",
         "doi": ""},
    60:
        {"autor": "Edwards, T.J., Newman, J., Prausnitz, J.M.",
         "title": "Thermodynamics of Aqueous Solutions Containing Volatile "
                  "Weak Electrolytes",
         "ref": "AIChE Journal 21(2) (1975) 248-259",
         "doi": "10.1002/aic.690210205"},
    61:
        {"autor": "Ely, J.F., Hanley, H.J.M.",
         "title": "A Computer Program for the Prediction of Viscosity and "
                  "Thermal Condcutivity in Hydrocarbon Mixtures",
         "ref": "NBS Technical Note 1039 (1981)",
         "doi": ""},
    62:
        {"autor": "ASHRAE",
         "title": "Designation and Safety Classification of Refrigerants",
         "ref": "Standard 34-2010",
         "doi": ""}
}


def atomic_decomposition(cmp):
    """Procedure to decompose a molecular string representation in its atomic
    composition. Support both expanded and short formula

    Parameters
    ----------
    cmp : string
        Compound formula

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


@refDoc(__doi__, [62])
def refrigerantCode(cmp):
    """ASHRAE refrigerant code, calculate only the numbers, the letter with
    isomers definitions and the inorganic R6x definition are unsupported
    because there aré very compound specific

    Parameters
    ----------
    cmp : string
        Compound expanded formula

    Returns
    -------
    code : string
        ASHRAE refrigerant code

    Examples
    --------
    >>> refrigerantCode("CF3CF=CH2")
    'R1234'

    """
    dcomp = atomic_decomposition(cmp)
    double = cmp.count("=")
    cyclo = cmp.count("cyclo")

    code = "R"
    if cyclo > 0:
        code += "C"
    if double > 0:
        code += str(double)

    code += str(dcomp["C"]-1)
    code += str(dcomp.get("H", 0)+1)
    code += str(dcomp["F"])

    return code


def DIPPR(prop, T, args, Tc=None, M=None):
    r"""Procedure to implement the DIPPR equations valid to calculate several
    physical properties of compounds.

    Parameters
    ----------
    prop : string
        Property to calculate, any of:
        rhoS, rhoL, Hv, Pv, cpS, cpL, cpG, muL, muG, kL, kG, sigma
    T : float
        Temperature, [K]
    args : list
        Coefficients for DIPPR equation, [eq, A, B, C, D, E]
    Tc : float, optional
        Critical temperature, [K]
    M : float, optional
        Molecular weight, [g/mol]

    Notes
    -----
    The properties this method can calculate, and the units for the calculated
    properties are:
        * rhoS: Solid density, [kmol/m³]
        * rhoL: Liquid density, [kmol/m³]
        * Pv: Vapor pressure, [Pa]
        * Hv: Heat of vaporization, [J/kmol]
        * cpS: Solid heat capacity, [J/kmol·K]
        * cpL: Liquid heat capacity, [J/kmol·K]
        * cpG: Ideal gas heat capacity, [J/kmol·K]
        * muL: Liquid viscosity, [Pa·s]
        * muG: Vapor viscosity, [Pa·s]
        * kL: Liquid thermal conductivity, [W/m·K]
        * kG: Vapor thermal conductivity, [W/m·K]
        * sigma: Surface Tension, [N/m]

    The first element in args define the equation to use:
        * Eq 1: :math:`Y = A+BT+CT^2+DT^3+ET^4`
        * Eq 2: :math:`Y = exp(A+BT+Cln(T)+DT^E)`
        * Eq 3: :math:`Y = A*T^B/(1+CT+DT^2)`
        * Eq 4: :math:`Y = A+Bexp(-C/T^D)`
        * Eq 5: :math:`Y = A + BT + CT^3 + DT^8 + ET^9`
        * Eq 6: :math:`Y = A/(B^(1+(1-T/C)^D)`
        * Eq 7: :math:`Y = A*(1-Tr)^(B+CTr+DTr^2+ETr^3)`
        * Eq 8: :math:`Y = A+ B*((C/T)/sinh(C/T))^2 + D*((E/T)/cosh(E/T))^2`
        * Eq 9: :math:`Y = A^2/Tr+B-2ACTr-ADTr^2-C^2Tr^3/3-CDTr^4/2-D^2Tr^5/5`

    where:
        * T: Temperature, [K]
        * Tr: Reduced temperature T/Tc
        * A,B,C,D,E: Parameters of equation

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
        mul = 1/M
    elif "cp" in prop:
        unit = unidades.SpecificHeat
        mul = 1/M
    elif "mu" in prop:
        unit = unidades.Viscosity
    elif "k" in prop:
        unit = unidades.ThermalConductivity
    elif prop == "sigma":
        unit = unidades.Tension

    eq, A, B, C, D, E = args
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
@refDoc(__doi__, [21, 35, 5])
def RhoL_Rackett(T, Tc, Pc, Zra, M):
    r"""Calculates saturated liquid densities of pure components using the
    modified Rackett equation by Spencer-Danner

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
        Molecular weight, [g/mol]

    Returns
    -------
    rho : float
        Saturated liquid density at T, [kg/m³]

    Examples
    --------
    Example from 5_; propane at 30ºF, API procedure 6A2.13 pag.454

    >>> T = unidades.Temperature(30, "F")
    >>> Tc = unidades.Temperature(206.06, "F")
    >>> Pc = unidades.Pressure(616, "psi")
    >>> "%0.3f" % RhoL_Rackett(T, Tc, Pc, 0.2763, 44.1).kgl
    '0.531'
    """
    Pc_atm = Pc/101325
    Tr = T/Tc
    V = R_atml*Tc/Pc_atm*Zra**(1+(1-Tr)**(2/7))
    return unidades.Density(M/V)


@refDoc(__doi__, [20, 5])
def RhoL_Costald(T, Tc, w, Vc):
    r"""Calculates saturated liquid densities of pure components using the
    Corresponding STAtes Liquid Density (COSTALD) method, developed by
    Hankinson and Thomson, referenced too in API procedure 6A2.15 pag. 462

    .. math::
        \frac{V}{V^{o}}=V_{R}^{(0)}\left(1-\omega_{SRK}V_{R}^{(1)}\right)

    .. math::
        V_{R}^{(0)}=1+a\left(1-T_{R}\right)^{1/3}+b\left(1-T_{R}\right)^{2/3}
        +c\left(1-T_{R}\right)+d\left(1-T_{R}\right)^{4/3}

    .. math::
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
    Example 1 from 5_; propane at 30ºF

    >>> T = unidades.Temperature(30, "F")
    >>> Tc = unidades.Temperature(206.01, "F")
    >>> Vc = unidades.SpecificVolume(3.205/44.097, "ft3lb")
    >>> "%0.3f" % RhoL_Costald(T, Tc, 0.1532, Vc).kgl
    '0.530'
    """
    a = -1.52816
    b = 1.43907
    c = -0.81446
    d = 0.190454
    e = -0.296123
    f = 0.386914
    g = -0.0427258
    h = -0.0480645

    if T > Tc:
        Tr = 1
    else:
        Tr = T/Tc

    # Eq 17
    Vr0 = 1 + a*(1-Tr)**(1/3) + b*(1-Tr)**(2/3) + c*(1-Tr) + d*(1-Tr)**(4/3)
    Vr1 = (e + f*Tr + g*Tr**2 + h*Tr**3)/(Tr-1.00001)                  # Eq 18

    # TODO: Add V* to the database
    V = Vc*Vr0*(1-w*Vr1)                                               # Eq 16
    return unidades.Density(1/V)


def RhoL_Cavett(T, Tc, M, Vliq):
    r"""Calculates saturated liquid densities of pure components using the
    Cavett equation. Referenced in Chemcad Physical properties user guide

    .. math::
        \frac{1}{\rho} = V_{liq}\left(5.7+3*T_r\right)

    Vliq is the liquid volume constant saved in database for many compounds

    Parameters
    ----------
    T : float
        Temperature [K]
    Tc : float
        Critical temperature [K]
    M : float
        Molecular weight, [g/mol]
    Vliq : float
        Liquid mole volume constant, [cm³/g]

    Returns
    -------
    rho : float
        Saturated liquid density at T, [kg/m³]
    """
    Tr = T/Tc
    V = Vliq*(5.7+3*Tr)/M
    return unidades.Density(1/V, "gcc")


@refDoc(__doi__, [26])
def RhoL_YenWoods(T, Tc, Vc, Zc):
    r"""Calculates saturation liquid density using the Yen-Woods correlation

    .. math::
        \rho_s/\rho_c = 1 + A(1-T_r)^{1/3} + B(1-T_r)^{2/3} + D(1-T_r)^{4/3}

    .. math::
        A = 17.4425 - 214.578Z_c + 989.625Z_c^2 - 1522.06Z_c^3

    .. math::
        B = -3.28257 + 13.6377Z_c + 107.4844Z_c^2-384.211Z_c^3
        \text{ if } Zc \le 0.26

    .. math::
        B = 60.2091 - 402.063Z_c + 501Z_c^2 + 641Z_c^3
        \text{ if } Zc \ge 0.26

    .. math::
        D = 0.93-B

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tc : float
        Critical temperature, [K]
    Vc : float
        Critical volume, [m^3/mol]
    Zc : float
        Critical compressibility factor, [-]

    Returns
    -------
    rhos : float
        Liquid density, [kg/m³]
    """
    Tr = T/Tc
    A = 17.4425 - 214.578*Zc + 989.625*Zc**2 - 1522.06*Zc**3           # Eq 4
    if Zc <= 0.26:
        B = -3.28257 + 13.6377*Zc + 107.4844*Zc**2 - 384.211*Zc**3     # Eq 5A
    else:
        B = 60.2091 - 402.063*Zc + 501.0*Zc**2 + 641.0*Zc**3           # Eq 5B
    D = 0.93 - B                                                       # Eq 6

    # Eq 2
    rhos = (1 + A*(1-Tr)**(1/3) + B*(1-Tr)**(2/3) + D*(1-Tr)**(4/3))/Vc
    return unidades.Density(rhos)


@refDoc(__doi__, [23, 27])
def RhoL_YamadaGunn(T, Tc, Pc, w, M):
    r"""Calculates saturation liquid volume, using Gunn-Yamada correlation

    .. math::
        V/V_sc = V_R^{(0)}\left(1-\omega\delta\right)

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
    M : float
        Molecular weight, [g/mol]

    Returns
    -------
    rhos : float
        Liquid density, [kg/m³]

    Notes
    -----
    The equation is defined in [27]_ in volumen terms.
    """
    Tr = T/Tc

    if Tr < 0.8:
        # Eq 3
        Vr = 0.33593-0.33953*Tr+1.51941*Tr**2-2.02512*Tr**3+1.11422*Tr**4
    elif Tr < 1:
        # Eq 4
        Vr = 1+1.3*(1-Tr)**0.5*log10(1-Tr)-0.50879*(1-Tr)-0.91534*(1-Tr)**2
    elif Tr == 1:
        Vr = 1

    # Eq 5
    d = 0.29607-0.09045*Tr-0.04842*Tr**2

    Zsc = 0.292-0.0967*w                                                # Eq 7
    Vsc = Zsc*R*Tc/Pc

    V = Vsc*Vr*(1-w*d)                                                  # Eq 1
    return unidades.Density(M/V, "gm3")


@refDoc(__doi__, [28])
def RhoL_Bhirud(T, Tc, Pc, w, M):
    r"""Calculates saturation liquid density using the Bhirud correlation

    .. math::
        \ln \frac{P_c V_s}{RT} = \ln U^{(0)} + \omega\ln U^{(1)}

    .. math::
        \ln U^{(0)} = 1.39644 - 24.076T_r + 102.615T_r^2 - 255.719T_r^3
        + 355.805T_r^4 - 256.671T_r^5 + 75.1088T_r^6

    .. math::
        \ln U^{(1)} = 13.4412 - 135.7437T_r + 533.380T_r^2 - 1091.453T_r^3
        + 1231.43T_r^4 - 728.227T_r^5 + 176.737T_r^6

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
    M : float
        Molecular weight, [g/mol]

    Returns
    -------
    rhos : float
        Liquid density, [kg/m³]

    Notes
    -----
    Raise :class:`NotImplementedError` if Tr is > 1
    """
    Tr = T/Tc
    if Tr <= 0.98:
        lnU0 = 1.39644 - 24.076*Tr + 102.615*Tr**2 - 255.719*Tr**3 \
            + 355.805*Tr**4 - 256.671*Tr**5 + 75.1088*Tr**6            # Eq 9
        lnU1 = 13.4412 - 135.7437*Tr + 533.380*Tr**2 - 1091.453*Tr**3 \
            + 1231.43*Tr**4 - 728.227*Tr**5 + 176.737*Tr**6            # Eq 10
    elif Tr < 1:
        # Interpolation data from Table 1
        Trs_ = [0.98, 0.982, 0.984, 0.986, 0.988, 0.99, 0.992, 0.994, 0.996,
                0.998, 0.999, 1]
        lnU0_ = [-1.6198, -1.604, -1.59, -1.578, -1.564, -1.548, -1.533,
                 -1.515, -1.489, -1.454, -1.425, -1.243]
        lnU1_ = [-0.4626, -0.459, -0.451, -0.441, -0.428, -0.412, -0.392,
                 -0.367, -0.337, -0.302, -0.283, -0.2629]

        lnU0 = interp1d(Trs_, lnU0_, kind='cubic')(Tr)
        lnU1 = interp1d(Trs_, lnU1_, kind='cubic')(Tr)
    else:
        raise NotImplementedError("Input out of bound")

    # Eq 8
    U = exp(lnU0 + w*lnU1)
    Vs = U*R*T/Pc
    return unidades.Density(M/Vs, "gm3")


@refDoc(__doi__, [29])
def RhoL_Mchaweh(T, Tc, Vc, w, delta):
    r"""Calculates saturated liquid density using the Mchaweh correlation

    .. math::
        \rho_s = \rho_c\rho_o\left[1+\delta_{SRK}\left(\alpha_{SRK}-1
        \right)^{1/3}\right]

    .. math::
        \rho_o = 1+1.169\tau^{1/3}+1.818\tau^{2/3}-2.658\tau+2.161\tau^{4/3}

    .. math::
        \tau = 1-\frac{(T_r)}{\alpha_{SRK}}

    .. math::
        \alpha_{SRK} = \left[1 + m\left(1-\sqrt{T_r}\right)\right]^2

    .. math::
        m = 0.480 + 1.574\omega - 0.176\omega^2

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tc : float
        Critical temperature, [K]
    Vc : float
        Critical volume, [m^3/kg]
    w : float
        Acentric factor, [-]
    delta : float
        Correlation parameter, [-]

    Returns
    -------
    rhos : float
        Liquid density, [kg/m³]
    """
    Tr = T/Tc

    m = 0.480 + 1.574*w - 0.176*w**2                                   # Eq 14
    alpha = (1+m*(1-Tr**0.5))**2                                       # Eq 15
    tau = 1 - Tr/alpha                                                 # Eq 17

    rho0 = 1 + 1.169*tau**(1/3) + 1.818*tau**(2/3) - 2.658*tau + \
        2.161*tau**(4/3)                                               # Eq 16

    rhos = rho0/Vc*(1+delta*(alpha-1)**(1/3))           # Eq 15
    return unidades.Density(rhos)


Mchaweh_d = {28: 0.57510,
             24: 2.09626,
             35: 2.42033,
             83: 1.37209,
             146: -0.13014,
             130: -1.03473,
             140: 1.75281,
             65: 0.07747,
             63: 6.79976,
             98: -3.25962,
             40: 1.44502,
             343: 2.68664,
             # biphenyl : 4.33458,
             6: 0.13196,
             49: 1.54186,
             48: 0.80200,
             105: 1.30100,
             172: 4.46227,
             # chloroethane, R-160 : 6.09099,
             25: 1.12156,
             30: 5.68825,
             # cumene : 3.93381,
             38: -0.45300,
             325: 2.70939,
             36: -0.28797,
             69: 1.79061,
             14: 5.12504,
             3: -1.46429,
             22: -0.46558,
             45: 2.94376,
             59: 0.94253,
             162: 0.10621,
             208: -3.03979,
             50: 0.59860,
             1: -19.75170,
             5: 3.10327,
             27: 0.98728,
             7: 0.09888,
             61: 4.21712,
             145: -1.10422,
             # Krypton : 4.50251,
             2: -3.20525,
             117: 0.71947,
             39: -0.57591,
             37: 0.59017,
             160: 0.23343,
             107: -5.63803,
             11: 1.87922,
             10: 1.19300,
             46: -0.79463,
             91: -3.71806,
             13: 3.60629,
             # Deuterium : -5.33454,
             8: 0.62738,
             90: 9.08709,
             12: 3.64147,
             47: -2.70491,
             19: 6.93421,
             229: -2.62285,
             400: -1.76731,
             57: 10.90360,
             4: -0.25595,
             23: 1.47939,
             # R-134 : -8.66475,
             # R-245 : 1.16965,
             # R-32 : 6.66817,
             # R-40 : 7.86026,
             # R-500 : -2.47818,
             # R-502 : 0.54781,
             # R-C318 : 1.04585,
             51: 3.29577,
             100: -1.31517,
             41: 1.39931,
             26: 0.19733,
             # xenon : -1.27292,
             }


@refDoc(__doi__, [32])
def RhoL_Riedel(T, Tc, Vc, w):
    r"""Calculates saturation liquid density using the Riedel correlation

    .. math::
        \rho_s/\rho_c = 1 + 0.85\left(1-T_r\right) +
        \left(1.6916+0.984\omega\right)\left(1-T_r\right)^{1/3}

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tc : float
        Critical temperature, [K]
    Vc : float
        Critical volume, [m^3/mol]
    w : float
        Acentric factor, [-]

    Returns
    -------
    rhos : float
        Liquid density, [kg/m³]
    """
    Tr = T/Tc

    rhos = (1+0.85*(1-Tr)+(1.6916+0.984*w)*(1-Tr)**(1/3))/Vc
    return unidades.Density(rhos)


@refDoc(__doi__, [33])
def RhoL_ChuehPrausnitz(T, Tc, Vc, w):
    r"""Calculates saturation liquid density using the Chueh-Prausnitz
    correlation

    .. math::
        V_s/V_c = V_R^{(0)} + \omega V_R^{(1)} + \omega^2V_R^{(2)}

    .. math::
        V_R^{(i)} = a^{(i)} + b^{(i)}T_R + c^{(i)}T_R^2 + d^{(i)}T_R^3 +
        e^{(i)}/T_R + f^{(i)}\ln{1-T_R}

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tc : float
        Critical temperature, [K]
    Vc : float
        Critical volume, [m^3/mol]
    w : float
        Acentric factor, [-]

    Returns
    -------
    rhos : float
        Liquid density, [kg/m³]
    """
    # Table 1
    a = (0.11917, 0.98465, -0.55314)
    b = (0.009513, -1.60378, -0.15793)
    c = (0.21091, 1.82484, -1.01601)
    d = (-0.06922, -0.61432, 0.34095)
    e = (0.07480, -0.34546, 0.46795)
    f = (-0.084476, 0.087037, -0.239938)

    Tr = T/Tc

    # Eq 6
    v0 = a[0] + b[0]*Tr + c[0]*Tr**2 + d[0]*Tr**3 + e[0]/Tr + f[0]*log(1-Tr)
    v1 = a[1] + b[1]*Tr + c[1]*Tr**2 + d[1]*Tr**3 + e[1]/Tr + f[1]*log(1-Tr)
    v2 = a[2] + b[2]*Tr + c[2]*Tr**2 + d[2]*Tr**3 + e[2]/Tr + f[2]*log(1-Tr)

    # Eq 5
    Vr = v0 + w*v1 + w**2*v2
    return unidades.Density(1/Vr/Vc)


@refDoc(__doi__, [22, 5])
def RhoL_TaitCostald(T, P, Tc, Pc, w, Ps, rhos):
    r"""Calculates compressed-liquid density, using the Thomson-Brobst-
    Hankinson generalization of Tait equation, also referenced in API procedure
    6A2.23 pag. 477

    .. math::
        V = V_s\left(1-C\ln\frac{B + P}{B + P_s}\right)

    .. math::
        \frac{B}{P_c} = -1 + a\tau^{1/3} + b\tau^{2/3} + d\tau + e\tau^{4/3}

    .. math::
        e = \exp(f + g\omega_{SRK} + h \omega_{SRK}^2)

    .. math::
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
    w : float
        Acentric factor (SRK optimized), [-]
    rhos : float
        Saturation liquid volume, [kg/m^3]

    Returns
    -------
    rho : float
        High-pressure liquid density, [kg/m^3]

    Examples
    --------
    Example from 5_; n-octane at 212ºF and 4410 psi

    >>> T = unidades.Temperature(212, "F")
    >>> P = unidades.Pressure(4410, "psi")
    >>> Tc = unidades.Temperature(564.22, "F")
    >>> Pc = unidades.Pressure(360.6, "psi")
    >>> Ps = unidades.Pressure(6.74, "psi")
    >>> rs = RhoL_Rackett(T, Tc, Pc, 0.2569, 114.232)
    >>> "%0.3f" % (1/rs.lbft3*114.232)
    '2.874'
    >>> "%0.3f" % RhoL_TaitCostald(T, P, Tc, Pc, 0.3962, Ps, rs).kgl
    '0.676'
    """
    Tr = T/Tc
    C = 0.0861488 + 0.0344483*w                                         # Eq 8
    e = exp(4.79594 + 0.250047*w + 1.14188*w**2)                        # Eq 7
    B = Pc*(-1 - 9.070217*(1-Tr)**(1/3) + 62.45326*(1-Tr)**(2/3) -
            135.1102*(1-Tr) + e*(1-Tr)**(4/3))

    # Eq 5
    rho = rhos/(1-C*log((B+P)/(B+Ps)))
    return unidades.Density(rho, "gl")


@refDoc(__doi__, [30])
def RhoL_ChangZhao(T, P, Tc, Pc, w, Ps, rhos):
    r"""Calculates compressed-liquid density, using the Chang-Zhao correlation

    .. math::
        V = V_s\frac{AP_c + C^{\left(D-T_r\right)^B}\left(P-P_{vp}\right)}
        {AP_c + C\left(P-P_{vp}\right)}

    .. math::
        A=\sum_{i=0}^{5}a_{i}T_{r}^{i}

    .. math::
        B=\sum_{j=0}^{2}b_{j}\omega^{j}

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
    w : float
        Acentric factor (SRK optimized), [-]
    rhos : float
        Saturation liquid volume, [kg/m^3]

    Returns
    -------
    rho : float
        High-pressure liquid density, [kg/m^3]
    """
    Tr = T/Tc
    Pr = P/Pc
    Psr = Ps/Pc

    # Constants, Table 2
    a = (99.42, -78.68, -75.18, 41.49, 7.257)
    b = (0.38144, -0.30144)

    A = sum(ai*Tr**i for i, ai in enumerate(a))                         # Eq 8
    B = sum(bi*w**i for i, bi in enumerate(b))                          # Eq 7

    # Eq 5
    rho = rhos*(A+2.81*(Pr-Psr))/(A+2.81**(1.1-Tr)**B*(Pr-Psr))         # Eq 9
    return unidades.Density(rho)


@refDoc(__doi__, [31, 1])
def RhoL_AaltoKeskinen(T, P, Tc, Pc, w, Ps, rhos):
    r"""Calculates compressed-liquid density, using the Aalto-Keskinen
    modification of Chang-Zhao correlation

    .. math::
        V = V_s\frac{AP_c + C^{\left(D-T_r\right)^B}\left(P-P_{vp}\right)}
        {AP_c + C\left(P-P_{vp}\right)}

    .. math::
        A = a_0 + a_1T_r + a_2T_r^3 + a_3T_r^6 + a_4/T_r

    .. math::
        B = b_0 + \omega_SRKb_1

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
    w : float
        Acentric factor (SRK optimized), [-]
    Ps : float
        Saturation pressure, [Pa]
    rhos : float
        Saturation liquid density, [kg/m³]

    Returns
    -------
    rho : float
        High-pressure liquid density, [kg/m³]

    Notes
    -----
    This correlation improve the Hankinson-Boost-Thomson and Chung-Huang method
    in the region near to Tc.

    Examples
    --------
    Example 4-8 from 1_; ammonia at 400bar at 300K and 400K

    >>> P = unidades.Pressure(400, "bar")
    >>> Pc = unidades.Pressure(113.53, "bar")
    >>> Ps1 = unidades.Pressure(10.61, "bar")
    >>> rs1 = unidades.Density(1/28.38*17.031, "gcc")
    >>> r1 = RhoL_AaltoKeskinen(300, P, 405.4, Pc, 0.256, Ps1, rs1)
    >>> Ps2 = unidades.Pressure(102.97, "bar")
    >>> rs2 = unidades.Density(1/49.15*17.031, "gcc")
    >>> r2 = RhoL_AaltoKeskinen(400, P, 405.4, Pc, 0.256, Ps2, rs2)
    >>> "%0.2f %0.2f" % (1/r1.gcc*17.031, 1/r2.gcc*17.031)
    '27.19 35.60'

    Example 4-9 from 1_; m-cresol at 3000bar and 503.15K

    >>> P = unidades.Pressure(3000, "bar")
    >>> Pc = unidades.Pressure(45.6, "bar")
    >>> Ps = unidades.Pressure(1, "bar")
    >>> rs = unidades.Density(1/127.31*108.14, "gcc")
    >>> r = RhoL_AaltoKeskinen(503.15, P, 705.7, Pc, 0.452, Ps, rs)
    >>> "%0.2f" % (1/r.gcc*108.14)
    '112.97'
    """
    Tr = T/Tc
    Pr = P/Pc
    Psr = Ps/Pc

    # Constants, Table 4
    a = (-170.335, -28.5784, 124.809, -55.5393, 130.010)
    b = (0.164813, -0.0914427)
    C = math.e

    A = a[0] + a[1]*Tr + a[2]*Tr**3 + a[3]*Tr**6 + a[4]/Tr             # Eq 17
    B = b[0] + b[1]*w                                                  # Eq 18

    # Eq 5
    rho = rhos*(A+C*(Pr-Psr))/(A+C**(1.00588-Tr)**B*(Pr-Psr))          # Eq 14
    return unidades.Density(rho)


@refDoc(__doi__, [37, 38])
def RhoL_AaltoKeskinen2(T, P, Tc, Pc, w, Ps, rhos):
    r"""Calculates compressed-liquid density, using the Aalto-Keskinen
    modification of Chang-Zhao correlation extended to a more high pressure
    range

    .. math::
        V = V_s\frac{AP_c + C^{\left(D-T_r\right)^B}\left(P-P_{vp}\right)^E}
        {AP_c + C\left(P-P_{vp}\right)^E}

    .. math::
        A = a_0 + a_1T_r + a_2T_r^3 + a_3T_r^6 + a_4/T_r

    .. math::
        B = b_0 + \frac{b_1}{b_2+\omega_SRK}

    .. math::
        C = c_1\left(1-T_r\right)^{c_2}+\left(1-\left(1-T_r\right)^{c_2}\right)
        \exp\left(c_3+c_4\left(P-P_s\right)\right)

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
    w : float
        Acentric factor (SRK optimized), [-]
    Ps : float
        Saturation pressure, [Pa]
    rhos : float
        Saturation liquid density, [kg/m³]

    Returns
    -------
    rho : float
        High-pressure liquid density, [kg/m³]

    Notes
    -----
    This correlation increase the high pressure range of previous Aalto
    Correlation

    Examples
    --------
    Selected values from experimental data for ethane used in correlation
    development, [38]_

    >>> from lib.mEoS import C2
    >>> et = C2()
    >>> Ps = et._Vapor_Pressure(293.608)
    >>> rhos = et._Liquid_Density(293.608)
    >>> args = (293.608, 71.4671e6, C2.Tc, C2.Pc, C2.f_acent, Ps, rhos)
    >>> "%0.2f" % RhoL_AaltoKeskinen2(*args).gcc
    '0.50'

    >>> Ps = et._Vapor_Pressure(281.789)
    >>> rhos = et._Liquid_Density(281.789)
    >>> args = (281.789, 8.4e6, C2.Tc, C2.Pc, C2.f_acent, Ps, rhos)
    >>> "%0.2f" % RhoL_AaltoKeskinen2(*args).gcc
    '0.41'
    """
    Tr = T/Tc
    Pr = P/Pc
    Psr = Ps/Pc

    # Constants, Table 2
    a = (482.85416, -1154.2977, 790.09727, -212.14413, 93.4904)
    b = (0.0264002, 0.42711522, 0.5)
    c = (9.2892236, 2.5103968, 0.59397220, 0.0010895002)
    E = 0.80329503

    A = a[0] + a[1]*Tr + a[2]*Tr**3 + a[3]*Tr**6 + a[4]/Tr             # Eq 5
    B = b[0] + b[1]/(b[2]+w)                                           # Eq 6
    C = c[0]*(1-Tr)**c[1] + (1-(1-Tr)**c[1])*exp(c[2]+c[3]*(Pr-Psr))   # Eq 7

    # Eq 4
    rho = rhos*(A+C*(Pr-Psr)**E)/(A+C**(1.00001-Tr)**B*(Pr-Psr)**E)    # Eq 14
    return unidades.Density(rho)


@refDoc(__doi__, [36])
def RhoL_Nasrifar(T, P, Tc, Pc, w, M, Ps, rhos):
    r"""Calculates compressed liquid density using the Nasrifar correlation

    .. math::
        \frac{v-v_{s}}{v_{\infty}-v_{s}}=C\Psi

    .. math::
        \Psi=\frac{J+L\left(P_{r}-P_{rs}\right)+M\left(P_{r}-P_{rs}\right)^{3}}
        {F+G\left(P_{r}-P_{rs}\right)+I\left(P_{r}-P_{rs}\right)^{3}}

    .. math::
        J=j_{0}+j_{1}\left(1-T_{r}\right)^{1/3}+j_{2}\left(1-T_{r}\right)^{2/3}

    .. math::
        F=f_{0}\left(1-T_{r}\right)

    .. math::
        C=c_{0}+c_{1}\omega_{SRK}

    .. math::
        v_{\infty}=\varOmega\frac{RT_{c}}{P_{c}}

    .. math::
        \varOmega=\varOmega_{0}+\varOmega_{1}\omega_{SRK}

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
    M : float
        Molecular weight, [g/mol]
    Ps : float
        Saturation pressure, [Pa]
    rhos : float
        Saturation liquid density, [kg/m^3]

    Returns
    -------
    rho : float
        Liquid density, [kg/m³]
    """
    # Table 1
    j = (1.3168e-3, 3.4448e-2, 5.4131e-2)
    L = 9.6840e-2
    m = 8.6761e-6
    f = 48.8756
    G = 0.7185
    i = 3.4031e-5
    c = (5.5526, -2.7659)
    om = (7.9019e-2, -2.8431e-2)

    Tr = T/Tc
    Pr = P/Pc
    Prs = Ps/Pc
    vs = 1/rhos

    J = j[0] + j[1]*(1-Tr)**(1/3) + j[2]*(1-Tr)**(2/3)                 # Eq 17
    F = f*(1-Tr)                                                       # Eq 18
    C = c[0] + c[1]*w                                                  # Eq 19
    OM = om[0] + om[1]*w                                               # Eq 21

    # R value change to mass base in unit m³Pa/kgK
    v_inf = OM*R*1000/M*Tc/Pc                                          # Eq 18

    phi = (J+L*(Pr-Prs)+m*(Pr-Prs)**3)/(F+G*(Pr-Prs)+i*(Pr-Prs)**3)    # Eq 16

    # Eq 9
    v = C*tanh(phi)*(v_inf-vs)+vs
    return unidades.Density(1/v)


@refDoc(__doi__, [33, 5])
def RhoL_API(T, P, Tc, Pc, SG, rhos):
    r"""Calculates compressed-liquid density, using the analytical expression
    of Lu Chart referenced in API procedure 6A2.22

    .. math::
        \rho_2 = \rho_1\frac{C_2}{SG}

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
    SG : float
        Specific gravity at 60ºF, [-]
    rhos : float
        Liquid density at 60ºF, [kg/m^3]

    Returns
    -------
    rho : float
        High-pressure liquid density, [kg/m^3]

    Examples
    --------
    Example from 5_; n-nonane at 220ºF and 1000psi, the API databook use the
    original Lu Chart so the result don't have to be exact

    >>> T = unidades.Temperature(220, "F")
    >>> P = unidades.Pressure(1000, "psi")
    >>> Tc = unidades.Temperature(610.7, "F")
    >>> Pc = unidades.Pressure(332, "psi")
    >>> rs = unidades.Density(44.94, "lbft3")
    >>> "%0.1f" % RhoL_API(T, P, Tc, Pc, 1.077, rs).lbft3
    '41.7'
    """
    Pr = P/Pc
    Tr = T/Tc
    A0 = 1.6368-0.04615*Pr+2.1138e-3*Pr**2-0.7845e-5*Pr**3-0.6923e-6*Pr**4
    A1 = -1.9693+0.21874*Pr-8.0028e-3*Pr**2-8.2328e-5*Pr**3+5.2604e-6*Pr**4
    A2 = 2.4638-0.36461*Pr+12.8763e-3*Pr**2+14.8059e-5*Pr**3-8.6895e-6*Pr**4
    A3 = -1.5841+0.25136*Pr-11.3805e-3*Pr**2+9.5672e-5*Pr**3+2.1812e-6*Pr**4
    C2 = A0 + A1*Tr + A2*Tr**2 + A3*Tr**3

    # Eq 1
    d2 = rhos*C2/SG
    return unidades.Density(d2)


# Vapor pressure correlations
@refDoc(__doi__, [1, 3])
def Pv_Antoine(T, args, Tc=None, base=math.e, Punit="mmHg"):
    r"""Vapor Pressure calculation procedure using the Antoine equation

    .. math::
        \log_{\text{base}} P^{\text{sat}} = A - \frac{B}{T+C}

    The method implement too the extended Antoine Equation

    .. math::
        \log_{10} P^{sat} = A - \frac{B}{T + C} + 0.43429x^n + Ex^8 + Fx^{12}

    .. math::
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
    Example 7-1 in 1_, furan at 309.429 K

    >>> P = Pv_Antoine(309.429, (4.1199, 1070.2, -44.32), base=10, Punit="bar")
    >>> "%0.4f" % P.bar
    '1.2108'
    """
    # Clean args input of null values
    if len(args) > 3 and args[3] is None:
        args = args[:3]

    if len(args) == 3:
        A, B, C = args
        Pv = base**(A-B/(T+C))
    elif len(args) == 7 and args[3] is not None and Tc is not None:
        A, B, C, n, E, F, to = args
        x = (T-to)/Tc
        if x <= 0:
            Pv = base**(A-B/(T+C))
        else:
            Pv = base**(A - B/(T+C) + 0.43429*x**n + E*x**8 + F*x**12)
    return unidades.Pressure(Pv, Punit)


@refDoc(__doi__, [4, 2])
def Pv_Lee_Kesler(T, Tc, Pc, w):
    r"""Calculates vapor pressure of a fluid using the Lee-Kesler correlation

    The vapor pressure is given by:

    .. math::
        \ln P_r = f^{(0)} + \omega f^{(1)}

    .. math::
        f^{(0)} = 5.92714-\frac{6.09648}{T_r}-1.28862\ln T_r + 0.169347T_r^6

    .. math::
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
    Example 1.2 from 4_; propane at 80ºF

    >>> T = unidades.Temperature(80, "F")
    >>> Tc = unidades.Temperature(666.01, "R")
    >>> Pc = unidades.Pressure(616.3, "psi")
    >>> "%0.0f" % Pv_Lee_Kesler(T, Tc, Pc, 0.1522).psi
    '144'
    """
    # Eq 17, pag 525
    Tr = T/Tc
    f0 = 5.92714 - 6.09648/Tr - 1.28862*log(Tr) + 0.169347*Tr**6
    f1 = 15.2518 - 15.6875/Tr - 13.4721*log(Tr) + 0.43577*Tr**6
    return unidades.Pressure(exp(f0 + w*f1)*Pc)


@refDoc(__doi__, [6, 1, 5, 7])
def Pv_Wagner(T, Tc, Pc, args):
    r"""Calculates vapor pressure of a fluid using the Wagner correlation

    .. math::
        \ln P^{v}= \ln P_c + \frac{a\tau + b \tau^{1.5} + c\tau^{3}
        + d\tau^6} {T_r}

    .. math::
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
    method implement the origintal form of Wagner as in 6_, with the
    parameters from McGarry. API use other same different form.
    """
    a, b, c, d = args
    Tr = T/Tc
    tau = 1-Tr
    Pv = Pc/Tr*exp(a*tau + b*tau**1.5 + c*tau**3 + d*tau**6)            # Eq 14
    return unidades.Pressure(Pv)


@refDoc(__doi__, [8, 1])
def Pv_AmbroseWalton(T, Tc, Pc, w):
    r"""Calculates vapor pressure of a fluid using the Ambrose-Walton
    corresponding-states correlation

    .. math::
        \ln P_r=f^{(0)}+\omega f^{(1)}+\omega^2f^{(2)}

    .. math::
        f^{(0)}=\frac{-5.97616\tau + 1.29874\tau^{1.5}- 0.60394\tau^{2.5}
        -1.06841\tau^5}{T_r}

    .. math::
        f^{(1)}=\frac{-5.03365\tau + 1.11505\tau^{1.5}- 5.41217\tau^{2.5}
        -7.46628\tau^5}{T_r}

    .. math::
        f^{(2)}=\frac{-0.64771\tau + 2.41539\tau^{1.5}- 4.26979\tau^{2.5}
        +3.25259\tau^5}{T_r}

    .. math::
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
    Example 7-3 from 2_; ethylbenzene at 347.25 K.

    >>> "%0.4f" % Pv_AmbroseWalton(347.25, 617.15, 36.09E5, 0.304).bar
    '0.1328'
    >>> "%0.3f" % Pv_AmbroseWalton(460, 617.15, 36.09E5, 0.304).bar
    '3.325'
    """
    Tr = T/Tc
    t = 1 - T/Tc
    f0 = (-5.97616*t + 1.29874*t**1.5 - 0.60394*t**2.5 - 1.06841*t**5)/Tr
    f1 = (-5.03365*t + 1.11505*t**1.5 - 5.41217*t**2.5 - 7.46628*t**5)/Tr
    f2 = (-0.64771*t + 2.41539*t**1.5 - 4.26979*t**2.5 + 3.25259*t**5)/Tr

    # Eq 8
    Pv = Pc*exp(f0 + w*f1 + w**2*f2)
    return unidades.Pressure(Pv)


@refDoc(__doi__, [1])
def Pv_Riedel(T, Tc, Pc, Tb):
    r"""Calculate vapor pressure of a fluid using the Rieel
    corresponding-states correlation

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
    Example 7-4 from 2_; ethylbenzene

    >>> "%0.3f" % Pv_Riedel(347.25, 617.15, 36.09E5, 409.36).bar
    '0.131'
    >>> "%0.2f" % Pv_Riedel(460, 617.15, 36.09E5, 409.36).bar
    '3.35'
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


@refDoc(__doi__, [5])
def Pv_MaxwellBonnel(T, Tb, Kw):
    r"""Calculates vapor pressure of a fluid using the Maxell-Bonnel
    correlation as explain in 5_, procedure 5A1.18, Pag. 394

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
    Example in 5_, tetralin at 302ºF

    >>> T = unidades.Temperature(302, "F")
    >>> Tb = unidades.Temperature(405.7, "F")
    >>> Pv = Pv_MaxwellBonnel(T, Tb, 9.78)
    >>> "%0.1f" % Pv.psi
    '3.1'
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


@refDoc(__doi__, [9])
def Pv_Sanjari(T, Tc, Pc, w):
    r"""Calculates vapor pressure of a fluid using the Sanjari correlation
    pressure, and acentric factor.

    The vapor pressure of a chemical at `T` is given by:

    .. math::
        P_{v} = P_c\exp(f^{(0)} + \omega f^{(1)} + \omega^2f^{(2)})

    .. math::
        f^{(0)} = a_1 + \frac{a_2}{T_r} + a_3\ln T_r + a_4 T_r^{1.9}

    .. math::
        f^{(1)} = a_5 + \frac{a_6}{T_r} + a_7\ln T_r + a_8 T_r^{1.9}

    .. math::
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


# Liquid viscosity correlations
def MuL_Parametric(T, args):
    r"""Calculates liquid viscosity using a paremtric equation

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


@refDoc(__doi__, [10])
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

    Examples
    --------
    Example 9.19 from 1_ 4Ed; propanol at 433.2K

    >>> Vc = 316/92.14/1000
    >>> "%0.3f" % MuL_LetsouStiel(433.2, 60.10, 536.8, 51.7e5, 0.623).cP
    '0.171'
    """
    Pc_atm = unidades.Pressure(Pc).atm
    Tr = T/Tc

    x = Tc**(1/6)/M**0.5/Pc_atm**(2/3)
    x0 = 0.015178 - 0.021351*Tr + 0.007503*Tr**2
    x1 = 0.042559 - 0.07675*Tr + 0.034007*Tr**2
    mu = (x0+w*x1)/x
    return unidades.Viscosity(mu, "cP")


@refDoc(__doi__, [45, 1])
def MuL_PrzedzieckiSridhar(T, Tc, Pc, Vc, w, M, Tf, Vr=None, Tv=None):
    r"""Calculates the viscosity of a liquid using the Przezdziecki-Sridhar
    correlation

    .. math::
        \frac{1}{\mu} = B \left(\frac{V-V_o}{V_o}\right)

    .. math::
        B = \frac{0.33V_c}{f_1}-1.12

    .. math::
        f_1 = 4.27+0.032M_w-0.077P_c+0.014T_f-3.82\frac{T_f}{T_c}

    .. math::
        V_o = 0.0085T_c\omega-2.02+\frac{V_{m}}{0.342(T_f/T_c)+0.894}

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tc : float
        Critical temperature, [K]
    Pc : float
        Critical pressure, [Pa]
    Vc : float
        Critical volume, [m³/kg]
    w : float
        Acentric factor, [-]
    M : float
        Molecular weight, [g/mol]
    Tf : float
        Melting point, [K]
    Vr : float, optional
        Liquid volume at Tr, [m³/kg]
    Tv : float, optional
        Temperature of known volume, [K]

    Returns
    -------
    mu : float
        Viscosity, [Pa·s]

    Notes
    -----
    The refernce volume is a optional volume, it use the critical point as
    default values

    Examples
    --------
    Example 9.196 from 1_; toluene at 383K

    >>> Vc = 316/92.14/1000
    >>> V = 106.87/92.14/1000
    >>> args = (383, 591.75, 41.08e5, Vc, 0.264, 92.14, 178, V, 298.15)
    >>> "%0.3f" % MuL_PrzedzieckiSridhar(*args).cP
    '0.223'
    """
    # Define default value por reference state:
    if Vr is None:
        Vr = Vc
        Tv = Tc

    # Pa in atm
    Pc = Pc/101325
    # Volume to mol base
    Vc = Vc*M*1000
    Vr = Vr*M*1000

    def f(Tr):
        G = 0.29607 - 0.09045*Tr - 0.04842*Tr**2                       # Eq 11

        # Eq 12
        Vr = .33593 - .33953*Tr + 1.51941*Tr**2 - 2.02512*Tr**3 + 1.11422*Tr**4
        f2 = Vr*(1-w*G)                                                # Eq 9
        return f2

    f2 = f(T/Tc)
    f2R = f(Tv/Tc)
    f2m = f(Tf/Tc)
    V = f2/f2R*Vr
    Vm = f2m/f2R*Vr

    # Eq 12
    Vo = 0.0085*Tc*w - 2.02 + Vm/(0.342*Tf/Tc + 0.894)

    f1 = 4.27 + 0.032*M - 0.077*Pc + 0.014*Tf - 3.82*Tf/Tc             # Eq 6
    B = 0.33*Vc/f1-1.12                                                # Eq 5

    # Eq 4
    mu = Vo/B/(V-Vo)
    return unidades.Viscosity(mu, "cP")


@refDoc(__doi__, [46, 1])
def MuL_Lucas(T, P, Tc, Pc, w, Ps, mus):
    r"""Calculate the viscosity of liquid at high pressure using the Lucas
    correlation

    .. math::
        \eta\left(T,P\right)=\eta_{S}\left(T\right)F_{p}\left(T_{r},P_{r},
                \omega\right)

    .. math::
        F_{p}\left(T_{r},P_{r},\omega\right)=\frac{F_{p}^{ref}\left(T_{r},P_{r}
          \right)}{1+F_{s}\left(T_{r},\omega\right)\left(P_{r}-P_{sr}\right)}

    .. math::
        F_{p}^{ref}\left(T_{r},P_{r}\right)=1+f_{2}\left(T_{r}\right)\left[
          \left(P_{r}-P_{sr}\right)/2.11824066\right]^{f_{1}\left(T_{r}\right)}

    .. math::
        f_{1}\left(T_{r}\right)=0.9990614-\frac{0.00046739}{1.052278T_{r}^
          {-0.03876963}-1.05134195}

    .. math::
        f_{2}\left(T_{r}\right)=-0.20863153+\frac{0.32569953}{\left(
          1.00383978-T_{r}^{2.57327058}\right)^{0.29063299}}

    .. math::
        f_{s}\left(T_{r},\omega\right)=\omega\left(-0.079206+2.161577T_{r}-
          13.403985T_{r}^{2}+44.170595T_{r}^{3}-84.829114T_{r}^{4}+
          96.120856T_{r}^{5}-59.812675T_{r}^{6}+15.671878T_{r}^{7}\right)

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
    w : float
        Acentric factor, [-]
    Ps : float
        Saturation pressure, [Pa]
    mus: float
        Viscosity of saturated liquid, [Pa·s]

    Returns
    -------
    mu : float
        Viscosity at high pressure, [Pa·s]

    Examples
    --------
    Example 9.15 from 1_, methylcyclohexane at 300K 500 bar

    >>> "%0.2f" % MuL_Lucas(300, 500e5, 572.19, 34.7e5, 0.236, 0, 0.00068).cP
    '1.07'

    Selected value from Table 1 in [46]_, hydrogen

    >>> from lib.mEoS import pH2
    >>> y = pH2.younglove
    >>> T = 0.904*y["Tc"]
    >>> P = 7.71*y["Pc"]*1e3
    >>> Ps = pH2()._Vapor_Pressure(T)
    >>> "%0.2f" % MuL_Lucas(T, P, y["Tc"], y["Pc"]*1e3, pH2.f_acent, Ps, 1)
    '1.00'
    """
    Tr = T/Tc
    # Discard point in gas phase
    if P < Ps:
        dPr = 0
    else:
        dPr = P/Pc-Ps/Pc

    f1 = 0.9990614 - 4.6739e-4/(1.052278*Tr**-0.03876963 - 1.05134195)  # Eq 4
    # Eq 5
    f2 = -0.20863153 + 0.32569953/(1.00383978 - Tr**2.57327058)**0.29063299
    Fs = w * (-0.079206 + 2.161577*Tr - 13.403985*Tr**2 + 44.170595*Tr**3
              - 84.829114*Tr**4 + 96.120856*Tr**5 - 59.812675*Tr**6
              + 15.671878*Tr**7)                                        # Eq 6
    Fpr = 1 + f2*(dPr/2.11824066)**f1                                   # Eq 3
    Fp = Fpr/(1+Fs*dPr)                                                 # Eq 2

    mu = mus*Fp                                                         # Eq 1
    return unidades.Viscosity(mu)


@refDoc(__doi__, [5])
def MuL_API(T, P, Tc, Pc, w, muc):
    """Calculate the viscosity of liquid at high pressure using the API
    correlation, API procedure 11A5.1, pag 1074.

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
    w : float
        Acentric factor, [-]
    muc: float
        Viscosity of critical point, [Pa·s]

    Returns
    -------
    mu : float
        Viscosity at high pressure, [Pa·s]

    Notes
    -----
    Procedure valid for low-molecular weight hydrocarbons at high pressure. A
    compound with lower than 20 carbon atoms are valid for this method.

    Examples
    --------
    Example from 5_; pentane at 200ºF and 3000 psi
    The reference has a typo in mur0 calculation, the correct value is 6.20
    (not 5.20) and so the calculated viscosity is 0.171 cP, near to the
    experimental value of 0.166 cP.

    >>> T = unidades.Temperature(200, "F")
    >>> Tc = unidades.Temperature(385.7, "F")
    >>> P = unidades.Pressure(3000, "psi")
    >>> Pc = unidades.Pressure(488.8, "psi")
    >>> "%0.3f" % MuL_API(T, P, Tc, Pc, 0.2515, 2.55e-5).cP
    '0.171'
    """
    Tr = T/Tc
    Pr = P/Pc

    # Eqs 11A5.1-3
    A1 = 3.0294*Tr**9.0740 + 0.0032*Tr**10.9399 - 0.3689
    A2 = -0.038*Tr**-7.2309 + 0.0229*Tr**11.7631 + 0.5781
    A3 = -0.1415*Tr**27.2842 + 0.0778*Tr**-4.3406 + 0.0014
    A4 = 0.0028*Tr**69.4404 - 0.0042*Tr**3.3586 + 0.0062
    A5 = 0.0107*Tr**-7.4626 - 85.8276*Tr**0.1392 + 87.3164

    # Eq 11A5.1-2
    mur0 = A1*log10(Pr) + A2*log10(Pr)**2 + A3*Pr + A4*Pr**2 + A5

    # Eqs 11A5.1-5
    if Pr <= 0.75:
        B1 = -0.2462*Tr**0.0484 - 0.7275*log(Tr) - 0.0588*Tr + 0.0079
        B2 = -0.3199*Tr**17.0626 - 0.0695*log(Tr) + 0.1267*Tr - 0.0101
        B3 = 4.7217*Tr**-1.9831 + 19.2008*Tr**-1.7595 + 65.5728*log(Tr) + \
            0.6110*Tr-19.1590
    else:
        B1 = -0.0214*Tr**0.0484 - 0.1827*log(Tr)-0.0183*Tr + 0.0090
        B2 = -0.3588*Tr**5.0537 - 0.1321*log(Tr)+0.0204*Tr - 0.0075
        B3 = 3.7266*Tr**-2.5689 + 52.1358*Tr**0.3514 - 13.0750*log(Tr) + \
            0.6358*Tr-56.6687

    # Eqs 11A5.1-4
    mur1 = B1*Pr + B2*log(Pr) + B3

    # Eqs 11A5.1-1
    mur = mur0 + w*mur1
    return unidades.Viscosity(mur*muc)


@refDoc(__doi__, [5])
def MuL_Kouzel(T, P, muo):
    """Calculate the viscosity of liquid at high pressure using the API
    correlation, API procedure 11A5.5, pag 1081.

    Parameters
    ----------
    T : float
        Temperature, [K]
    P : float
        Pressure, [Pa]
    muo: float
        Viscosity of atmospheric pressure, [Pa·s]

    Returns
    -------
    mu : float
        Viscosity at high pressure, [Pa·s]

    Notes
    -----
    Procedure valid for high-molecular weight hydrocarbons at high pressure. A
    compound with more than 20 carbon atoms are valid for this method.

    Examples
    --------
    Example from 5_; lubricating oil at 120.2ºF and 9940 psi

    >>> T = unidades.Temperature(120.2, "F")
    >>> P = unidades.Pressure(9940, "psi")
    >>> "%0.1f" % MuL_Kouzel(T, P, 0.0527).cP
    '277.2'
    """
    # Unit conversion
    psig = unidades.Pressure(P).psig
    muocp = unidades.Viscosity(muo).cP

    # Eq 11A5.5-1
    mu = muocp*10**(psig/1000*(-0.0102+0.04042*muocp**0.181))

    return unidades.Viscosity(mu, "cP")


# Vas viscosity correlations
@refDoc(__doi__, [1])
def MuG_ChapmanEnskog(T, M, sigma, omega):
    r"""Calculate the viscosity of a gas using the Chapman-Enskog correlation

    .. math::
        \mu=\frac{26.69\left(MT\right)^{1/2}}{\sigma^2\Omega_v}

    Parameters
    ----------
    T : float
        Temperature, [K]
    M : float
        Molecular weight, [g/mol]
    sigma : float
        hard sphere diameter, [Å]
    omega : float
        Collision integral, [-]

    Returns
    -------
    mu : float
        Viscosity of gas, [Pa·s]
    """
    mu = 26.69*M**0.5*T**0.5/sigma**2/omega
    return unidades.Viscosity(mu, "microP")


@refDoc(__doi__, [24, 5])
def MuG_StielThodos(T, Tc, Pc, M):
    r"""Calculate the viscosity of a gas using the Stiel-Thodos correlation,
    also referenced in API procedure 11B1.3, pag 1099

    .. math::
        \mu=N/\xi

    .. math::
        \xi=\frac{T_{c}^{1/6}}{M^{1/2}P_{c}^{2/3}}

    .. math::
        N=3.4e^{-4}T_{r}^{0.94}   for Tr ≤ 1.5

    .. math::
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
    Example A in 5_, Propane at 176ºF

    >>> T = unidades.Temperature(176, "F")
    >>> Tc = unidades.Temperature(206, "F")
    >>> Pc = unidades.Pressure(616, "psi")
    >>> "%0.4f" % MuG_StielThodos(T, Tc, Pc, 44.1).cP
    '0.0100'

    Example B in 5_, Methane at 543ºF

    >>> T = unidades.Temperature(543, "F")
    >>> Tc = unidades.Temperature(-116.67, "F")
    >>> Pc = unidades.Pressure(667, "psi")
    >>> "%0.4f" % MuG_StielThodos(T, Tc, Pc, 16.04).cP
    '0.0176'
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


@refDoc(__doi__, [57])
def MuG_Gharagheizi(T, Tc, Pc, M):
    r"""Calculates the viscosity of a gas using the Gharagheizi et al.
    correlation

    .. math::
        \mu = 10^{-5} P_cT_r + \left(0.091-\frac{0.477}{M}\right)T +
        M \left(10^{-5}P_c-\frac{8M^2}{T^2}\right)
        \left(\frac{10.7639}{T_c}-\frac{4.1929}{T}\right)

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
    Methane at 120K

    >>> "%0.6e" % MuG_Gharagheizi(120, 190.564, 45.99e5, 16.04246)
    '5.215762e-06'

    1-Octanol at 120K

    >>> "%0.6e" % MuG_Gharagheizi(468.35, 652.5, 27.77e5, 130.22792)
    '8.751141e-06'
    """
    Tr = T/Tc

    # Eq 4
    mu = 1e-5*Pc*Tr + (0.091-0.477/M)*T + \
        M*(1e-5*Pc-8*M**2/T**2)*(10.7639/Tc - 4.1929/T)
    return unidades.Viscosity(mu*1e-7)


@refDoc(__doi__, [56])
def MuG_YoonThodos(T, Tc, Pc, M):
    r"""Calculates the viscosity of a gas using an Yoon-Thodos correlation

    .. math::
        \eta^o\xi = 46.1T_r^{0.618}-20.4\exp(-0.449T_r)+19.4\exp(-4.058T_r)+1

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

    Notes
    -----
    This method is valid only for nonpolar gases
    """
    Pc_atm = Pc/101325

    Tr = T/Tc
    x = Tc**(1/6)/M**0.5/Pc_atm**(2/3)

    if M == 2.0158:
        # Eq 3, Hydrogen case
        mur = 47.65*Tr**0.657 - 20*exp(-0.858*Tr) + 19*exp(-3.995*Tr) + 1
    elif M == 4.0026:
        # Eq 4, Helium case
        mur = 52.57*Tr**0.656 - 18.9*exp(-1.144*Tr) + 17.9*exp(-5.182*Tr) + 1
    else:
        # Eq 2, General case for nonpolar gases
        mur = 46.1*Tr**0.618 - 20.4*exp(-0.449*Tr) + 19.4*exp(-4.058*Tr) + 1

    return unidades.Viscosity(mur*1e-5/x, "cP")


@refDoc(__doi__, [49, 50, 1])
def MuG_Chung(T, Tc, Vc, M, w, D, k=0):
    r"""Calculate the viscosity of a gas using the Chung et al. correlation

    .. math::
        \mu=40.785\frac{F_c\left(MT\right)^{1/2}}{V_c^{2/3}\Omega_v}

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
    w : float
        Acentric factor, [-]
    D : float
        Dipole moment, [Debye]
    k : float, optional
        Corection factor for polar substances, [-]

    Returns
    -------
    mu : float
        Viscosity of gas, [Pa·s]

    Examples
    --------
    Example 9-1 in 1_, SO2 at 300ºC

    >>> T = unidades.Temperature(300, "C")
    >>> "%0.1f" % MuG_Chung(T, 430.8, 122e-3/64.065, 64.065, 0.257, 1.6).microP
    '245.5'
    """
    # Vc in molar base
    Vc = Vc*M*1000

    T_ = 1.2593*T/Tc
    omega = Collision_Neufeld(T_)
    mur = 131.3*D/(Vc*Tc)**0.5                                          # Eq 8
    Fc = 1 - 0.2756*w + 0.059035*mur**4 + k                             # Eq 7
    mu = 40.785*Fc*M**0.5*T**0.5/Vc**(2/3)/omega                        # Eq 6
    return unidades.Viscosity(mu, "microP")


@refDoc(__doi__, [49, 50, 1])
def MuG_P_Chung(T, Tc, Vc, M, w, D, k, rho, muo):
    r"""Calculate the viscosity of a compressed gas using the Chung correlation

    .. math::
        \mu=40.785\frac{F_c\left(MT\right)^{1/2}}{V_c^{2/3}\Omega_v}

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tc : float
        Critical temperature, [K]
    Vc : float
        Critical volume, [m³/kg]
    M : float
        Molecular weight, [g/mol]
    w : float
        Acentric factor, [-]
    D : float
        Dipole moment, [Debye]
    k : float, optional
        Corection factor for polar substances, [-]
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
    Example 9-12 in 1_, ammonia at 520K and 600bar

    >>> Vc = 72.4/17.031/1e3
    >>> rho = 1/48.2*17.031*1e3
    >>> mu = MuG_P_Chung(520, 405.5, Vc, 17.031, 0.256, 1.47, 0, rho, 182e-7)
    >>> "%0.0f" % mu.microP
    '455'
    """
    # units in molar base
    Vc = Vc*M*1000
    rho = rho/M/1000

    T_ = 1.2593*T/Tc
    mur = 131.3*D/(Vc*Tc)**0.5                                          # Eq 8

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
        A.append(ao + a1*w + a2*mur**4 + a3*k)
    A1, A2, A3, A4, A5, A6, A7, A8, A9, A10 = A

    Y = rho*Vc/6
    G1 = (1-0.5*Y)/(1-Y)**3
    G2 = (A1*((1-exp(-A4*Y))/Y)+A2*G1*exp(A5*Y)+A3*G1)/(A1*A4+A2+A3)

    muk = 10*muo*(1/G2 + A6*Y)
    mup = (36.344e-6*(M*Tc)**0.5/Vc**(2/3))*A7*Y**2*G2*exp(A8+A9/T_+A10/T_**2)

    return unidades.Viscosity(muk+mup, "P")


@refDoc(__doi__, [1])
def MuG_Reichenberg(T, P, Tc, Pc, Vc, M, D, muo):
    r"""Calculate the viscosity of a compressed gas using the Reichenberg
    correlation as explain in 1_

    .. math::
        \frac{\mu}{\mu^o}=1+Q\frac{AP_r^{3/2}}{BP_r+\left(1+CP_r^D\right)^{-1}}

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tc : float
        Critical temperature, [K]
    Pc : float
        Critical pressure, [Pa]
    Vc : float
        Critical volume, [m³/kg]
    M : float
        Molecular weight, [g/mol]
    D : float
        Dipole moment, [Debye]
    muo : float
        Viscosity of low-pressure gas, [Pa·s]

    Returns
    -------
    mu : float
        Viscosity of gas, [Pa·s]

    Examples
    --------
    Example 9-9 in 1_, n-pentane at 500K and 101bar

    >>> mu = MuG_Reichenberg(500, 101e5, 469.7, 33.7e5, 0, 0, 0, 114e-7)
    >>> "%0.0f" % mu.microP
    '520'
    """
    # units in molar base
    if D:
        Vc = Vc*M/1000
        mu_r = 131.3*D/(Vc*Tc)**0.5
        Q = 1-5.655*mu_r
    else:
        Q = 1

    Tr = T/Tc
    Pr = P/Pc

    A = 1.9824e-3/Tr*exp(5.2683*Tr**-0.5767)
    B = A*(1.6552*Tr-1.276)
    C = 1.319/Tr*exp(3.7035*Tr**-79.8678)
    D = 2.9496/Tr*exp(2.9190*Tr**-16.6169)

    mur = 1+Q*A*Pr**1.5/(B*Pr+1/(1+C*Pr**D))
    return unidades.Viscosity(mur*muo)


@refDoc(__doi__, [1])
def MuG_Lucas(T, P, Tc, Pc, Zc, M, D):
    """Calculate the viscosity of a gas using the Lucas correlation
    as explain in 1_. This method can calculate the viscosity at any pressure

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
    Zc : float
        Critical compressibility factor, [-]
    M : float
        Molecular weight, [g/mol]
    D : float
        Dipole moment, [Debye]

    Returns
    -------
    mu : float
        Viscosity of gas, [Pa·s]

    Examples
    --------
    Example 9-2 in 1_, methanol at 550K and 1bar

    >>> mu = MuG_Lucas(550, 1e5, 512.64, 80.97e5, 0.224, 32.042, 1.7)
    >>> "%0.0f" % mu.microP
    '178'

    Example 9-10 in 1_, ammonia at 420K and 300bar

    >>> muo = MuG_Lucas(420, 1e5, 405.5, 113.53e5, 0.244, 17.031, 1.47)
    >>> "%0.0f" % muo.microP
    '147'
    >>> mu = MuG_Lucas(420, 3e7, 405.5, 113.53e5, 0.244, 17.031, 1.47)
    >>> "%0.0f" % mu.microP
    '603'
    """
    Tr = T/Tc
    Pr = P/Pc

    Pc_bar = Pc*1e-5
    xi = 0.176*Tc**(1/6)/M**0.5/Pc_bar**(2/3)

    # Polarity and quantum effects correction factors
    mur = 52.46*D**2*Pc_bar/Tc**2
    if mur < 0.022:
        Fpo = 1
    elif mur < 0.075:
        Fpo = 1 + 30.55*(0.292-Zc)**1.72
    else:
        Fpo = 1 + 30.55*(0.292-Zc)**1.72*abs(0.96+0.1*(Tr-0.7))

    if Tr < 12:
        sign = -1
    else:
        sign = 1

    if M == 2.0158:
        Q = 0.76  # Hydrogen
        Fqo = 1.22*Q**0.15*(1+0.00385*((Tr-12)**2)**(1/M)*sign)
    elif M == 4.0026:
        Q = 1.38  # Helium
        Fqo = 1.22*Q**0.15*(1+0.00385*((Tr-12)**2)**(1/M)*sign)
    else:
        Fqo = 1

    Z1 = Fpo*Fqo*(0.807*Tr**0.618 - 0.357*exp(-0.449*Tr)
                  + 0.34*exp(-4.058*Tr) + 0.018)

    if Pr < 0.6:
        # Low pressure correlation
        mu = Z1/xi

    else:
        # High pressure correlation

        if Tr <= 1:
            alfa = 3.262 + 14.98*Pr**5.508
            beta = 1.39 + 14.98*Pr
            Z2 = 0.6 + 0.76*Pr**alfa + (6.99*Pr**beta-0.6)*(1-Tr)
        else:
            a = 1.245e-3/Tr*exp(5.1726*Tr**-0.3286)
            b = a*(1.6553*Tr-1.2723)
            c = 0.4489/Tr*exp(3.0578*Tr**-37.7332)
            d = 1.7368/Tr*exp(2.231*Tr**-7.6351)
            e = 1.3088
            f = 0.9425*exp(-0.1853*Tr**0.4489)
            Z2 = Z1*(1+a*Pr**e/(b*Pr**f+1/(1+c*Pr**d)))

        Y = Z2/Z1
        Fp = (1+(Fpo-1)/Y**3)/Fpo
        Fq = (1+(Fqo-1)*(1/Y-0.007*log(Y)**4))/Fqo
        mu = Z2*Fp*Fq/xi

    return unidades.Viscosity(mu, "microP")


@refDoc(__doi__, [51])
def MuG_Jossi(Tc, Pc, rhoc, M, rho, muo):
    r"""Calculate the viscosity of a compressed gas using the Jossi correlation

    .. math::
        \left[\left(\mu-\mu^o\right)\xi_T+1\right]^{1/4}=1.023+0.23364\rho_r+
        0.58533\rho_r^2-0.40758\rho_r^3+0.093324\rho_r^4

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tc : float
        Critical temperature, [K]
    Pc : float
        Critical pressure, [Pa]
    rhoc : float
        Critical density, [kg/m3]
    M : float
        Molecular weight, [g/mol]
    rho : float
        Density, [kg/m3]
    muo : float
        Viscosity of low-pressure gas, [Pa·s]

    Returns
    -------
    mu : float
        Viscosity of gas, [Pa·s]

    Notes
    -----
    This method is valid only for non polar substances, the paper give
    alternate equations for hydrogen, water and ammonia but there isn't a
    general correlation for polar compounds.

    Examples
    --------
    Example 9-11 in 1_, isobutane at 500K and 100bar

    >>> rhoc = 1/262.7*58.123*1000
    >>> rho = 1/243.8*58.123*1000
    >>> "%0.0f" % MuG_Jossi(407.85, 36.4e5, rhoc, 58.123, rho, 120e-7).microP
    '275'
    """
    Pc_atm = Pc/101325

    x = Tc**(1/6)/M**0.5/Pc_atm**(2/3)
    rhor = rho/rhoc

    # Eq 8
    mur = 0.1023 + 0.023364*rhor + 0.058533*rhor**2 - 0.040758*rhor**3 + \
        0.0093324*rhor**4
    mu = (mur**4-1e-4)/x+muo*1e3
    return unidades.Viscosity(mu, "cP")


@refDoc(__doi__, [52])
def MuG_P_StielThodos(Tc, Pc, rhoc, M, rho, muo):
    r"""Calculate the viscosity of a compressed gas using the Stiel-Thodos
    correlation. This method is valid for polar substances.

    .. math::
        \left(\mu-\mu^o\right)\xi=1.656\rho_r^{1.111}, \rho_r ≤ 0.1

    .. math::
        \left(\mu-\mu^o\right)\xi=0.0607\left(9.045\rho_r+0.63\right)^{1.739},
        0.1 ≤ \rho_r ≤ 0.9

    .. math::
        log\left[4-log\left(\left(\mu-\mu^o\right)\xi\right)\right]=
        0.6439-0.1005\rho_r-\Delta, 0.9 ≤ \rho_r ≤2.6

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tc : float
        Critical temperature, [K]
    Pc : float
        Critical pressure, [Pa]
    rhoc : float
        Critical density, [kg/m3]
    M : float
        Molecular weight, [g/mol]
    rho : float
        Density, [kg/m3]
    muo : float
        Viscosity of low-pressure gas, [Pa·s]

    Returns
    -------
    mu : float
        Viscosity of gas, [Pa·s]
    """
    Pc_atm = Pc/101325
    x = Tc**(1/6)/M**0.5/Pc_atm**(2/3)
    rhor = rho/rhoc

    if rhor <= 0.1:
        mur = 1.656e-5*rhor**1.111                                       # Eq 4
    elif 0.1 < rhor <= 0.9:
        mur = 0.607e-5*(9.045*rhor+0.63)**1.739                          # Eq 5
    else:
        if rhor < 2.2:
            D = 0
        else:
            D = 4.75e-4*(rhor**3-10.65)**2
        mur = 10**(4-10**(0.6439-0.1005*rhor-D))                         # Eq 6
    return unidades.Viscosity(mur/x+muo*1e7, "microP")


@refDoc(__doi__, [1, 61, 53])
def MuG_TRAPP(T, Tc, Vc, Zc, M, w, rho, muo):
    """Calculate the viscosity of a compressed gas using the TRAPP (TRAnsport
    Property Prediction) method.

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tc : float
        Critical temperature, [K]
    Zc : float
        Critical pressure, [Pa]
    rhoc : float
        Critical density, [kg/m3]
    M : float
        Molecular weight, [g/mol]
    rho : float
        Density, [kg/m3]
    muo : float
        Viscosity of low-pressure gas, [Pa·s]

    Returns
    -------
    mu : float
        Viscosity of gas, [Pa·s]

    Examples
    --------
    Example 9-13 in 1_, isobutane at 500K and 100bar

    >>> Vc = 259/58.124/1000
    >>> rho = 1/243.8*58.123*1000
    >>> mu = MuG_TRAPP(500, 407.85, Vc, 0.278, 58.124, 0.186, rho, 120e-7)
    >>> "%0.0f" % mu.microP
    '267'
    """
    # Reference fluid properties, propane
    TcR = 369.83
    rhocR = 1/200  # mol/cm³
    ZcR = 0.276
    wR = 0.152

    # Convert volume to molar base
    Vc = Vc*M*1000
    rho = rho/M/1000
    rhoc = 1/Vc

    Tr = T/Tc

    f = Tc/TcR*(1+(w-wR)*(0.05203-0.7498*log(Tr)))
    h = rhocR/rhoc*ZcR/Zc*(1-(w-wR)*(0.1436-0.2822*log(Tr)))

    To = T/f
    rho0 = rho*h
    Fn = (M/44.094*f)**0.5/h**(2/3)

    # Coefficients in [53]_, pag 796
    # Density are in mol/dm³
    rho0 *= 1000
    rhocR *= 1000

    G = -14.113294896 + 968.22940153/To
    H = rho0**0.5*(rho0-rhocR)/rhocR
    G2 = 13.686545032 - 12511.628378/To**1.5
    G3 = 0.0168910864 + 43.527109444/To + 7659.4543472/To**2
    F = G + G2*rho0**0.1+G3*H
    muR = exp(F)-exp(G)
    return unidades.Viscosity(Fn*muR*1e-6 + muo)


@refDoc(__doi__, [54])
def MuG_Brule(T, Tc, Vc, M, w, rho, muo):
    r"""Calculate the viscosity of a compressed gas using the Chung correlation

    .. math::
        \mu=40.785\frac{F_c\left(MT\right)^{1/2}}{V_c^{2/3}\Omega_v}

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tc : float
        Critical temperature, [K]
    Vc : float
        Critical volume, [m³/kg]
    M : float
        Molecular weight, [g/mol]
    w : float
        Acentric factor, [-]
    rho : float
        Density, [kg/m³]
    muo : float
        Viscosity of low-pressure gas, [Pa·s]

    Returns
    -------
    mu : float
        Viscosity of gas, [Pa·s]
    """
    # units in molar base
    Vc = Vc*M*1000
    rho = rho/M/1000

    T_ = 1.2593*T/Tc

    # Table I
    dat = {
        0: (1., -0.2756),
        1: (17.4499, 34.0631),
        2: (-0.961125e-3, 0.723459e-2),
        3: (51.0431, 169.460),
        4: (-0.605917, 71.1743),
        5: (21.3818, -2.11014),
        6: (4.66798, -39.9408),
        7: (3.76241, 56.6234),
        8: (1.00377, 3.13962),
        9: (-0.7774233e-1, -3.58446),
        10: (0.317523, 1.15995)}

    E = []
    for i in range(1, 11):
        ai, bi = dat[i]
        E.append(ai + bi*w)
    E1, E2, E3, E4, E5, E6, E7, E8, E9, E10 = E

    Y = rho*Vc/6
    G1 = (1-0.5*Y)/(1-Y)**3
    G2 = (E1*((1-exp(-E4*Y))/Y)+E2*G1*exp(E5*Y)+E3*G1)/(E1*E4+E2+E3)

    muk = muo*(1/G2 + E6*Y)
    mup = (36.344e-6*(M*Tc)**0.5/Vc**(2/3))*E7*Y**2*G2*exp(E8+E9/T_+E10/T_**2)

    return unidades.Viscosity(muk+mup, "P")


@refDoc(__doi__, [55, 5])
def MuG_DeanStiel(Tc, Pc, rhoc, M, rho, muo):
    r"""Calculate the viscosity of a compressed gas using the Dean-Stiel
    correlation, also referenced in API databook Procedure 11B4.1, pag 1107

    .. math::
        \left(\mu-\mu_o\right)\xi=10.8x10^{-5}\left[exp\left(1.439\rho_r\right)
        -exp\left(-1.11\rho_r^{1.858}\right)\right]

    Parameters
    ----------
    Tc : float
        Critical temperature, [K]
    Pc : float
        Critical pressure, [Pa]
    rhoc : float
        Critical density, [kg/m³]
    M : float
        Molecular weight, [g/mol]
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
    Example in 5_, mixture at 1500psi and 257ºF

    >>> Tc = unidades.Temperature(472.09, "R")
    >>> Pc = unidades.Pressure(646.68, "psi")
    >>> "%0.4f" % MuG_DeanStiel(Tc, Pc, 1, 27.264, 0.5283, 123e-7).cP
    '0.0163'
    """
    Pc_atm = Pc/101325

    x = Tc**(1/6)/M**0.5/Pc_atm**(2.0/3)
    rhor = rho/rhoc

    # Eq 13
    mur = 10.8e-5*(exp(1.439*rhor)-exp(-1.11*rhor**1.858))
    return unidades.Viscosity(muo*1e3 + mur/x, "cP")


@refDoc(__doi__, [5])
def MuG_API(T, P, Tc, Pc, muo):
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
    Tc : float
        Critical temperature, [K]
    Pc : float
        Critical pressure, [Pa]
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
    Example in 5_, nitrogen at -58ºF and 1677psi

    >>> T = unidades.Temperature(-58, "F")
    >>> Tc = unidades.Temperature(-232.5, "F")
    >>> P = unidades.Pressure(1677, "psi")
    >>> Pc = unidades.Pressure(493.1, "psi")
    >>> "%0.4f" % MuG_API(T, P, Tc, Pc, 1.44e-5).cP
    '0.0203'
    """
    Tr = T/Tc
    Pr = P/Pc

    A1 = 83.8970*Tr**0.0105 + 0.603*Tr**-0.0822 + 0.9017*Tr**-0.12 - 85.308
    A2 = 1.514*Tr**-11.3036 + 0.3018*Tr**-0.6856 + 2.0636*Tr**-2.7611
    mur = A1*1.5071*Pr**-0.4487 + A2*(
        11.4789*Pr**0.2606 - 12.6843*Pr**0.1773 + 1.6953*Pr**-0.1052)
    return unidades.Viscosity(muo*mur)


# Liquid thermal conductivity correlations
@refDoc(__doi__, [11])
def ThL_RiaziFaghri(T, Tb, SG):
    r"""Calculates thermal conductivity of liquid hydrocarbon at low pressure
    using the Riazi-Faghri correlation.

    .. math::
        \kappa = aT_{b}^{b}SG^{c}

    .. math::
        a = \exp\left(-4.5093-0.6844t-0.1305t^{2}\right)

    .. math::
        b = 0.3003+0.0918t+0.0195t^{2}

    .. math::
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
    """
    # Convert input Tb in Kelvin to Fahrenheit to use in the correlation
    Tb_R = unidades.K2R(Tb)
    t = unidades.K2F(T)/100
    A = exp(-4.5093-0.6844*t-0.1305*t**2)                              # Eq 9
    B = 0.3003+0.0918*t+0.0195*t**2                                    # Eq 10
    C = 0.1029+0.0894*t+0.0292*t**2                                    # Eq 11
    k = A*Tb_R**B*SG**C                                                # Eq 7
    return unidades.ThermalConductivity(k, "BtuhftF")


@refDoc(__doi__, [12])
def ThL_Gharagheizi(T, Pc, Tb, M, w):
    r"""Calculates the thermal conductivity of liquid using the Gharagheizi
    correlation.

    .. math::
        \kappa = 10^{-4}\left(10\omega+2P_c-2T+4+1.908\left(T_b+\frac{1.009B^2}
        {M^2}\right)+\frac{3.9287M^4}{B^4}+\frac{A}{B^8}\right)

    .. math::
        A = 3.8588M^8\left(1.0045B+6.5152M-8.9756\right)

    .. math::
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
    """

    Pc_bar = Pc*1e-5
    B = 16.0407*M+2*Tb-27.9074                                          # Eq 6
    A = 3.8588*M**8*(1.0045*B + 6.5152*M - 8.9756)                      # Eq 5
    k = 1e-4*(10*w + 2*Pc_bar - 2*T + 4 + 1.908*(Tb+1.009*B**2/M**2)
              + 3.9287*M**4/B**4 + A/B**8)                              # Eq 4
    return unidades.ThermalConductivity(k)


@refDoc(__doi__, [13])
def ThL_LakshmiPrasad(T, M):
    r"""Calculates the thermal conductivity of liquid using the Lakshmi-Prasad
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
        Thermal conductivity, [W/m·K]
    """
    # Eq 5
    k = 0.0655 + (1.3855 - 0.00197*T)/M**0.5 - 0.00005*T
    return unidades.ThermalConductivity(k)


@refDoc(__doi__, [14])
def ThL_Nicola(T, M, Tc, Pc, w, mu=None):
    r"""Calculates the thermal conductivity of liquid using the Nicola
    correlation.

    .. math::
        \frac{\lambda}{\lambda_o} = aT_r + bPc + c\omega +
        \left(\frac{e}{M}\right)^{d}

    .. math::
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
    """
    Pc_bar = unidades.Pressure(Pc).bar
    if mu:
        # Eq 4 using dipole moment of compound
        k = 0.6542*(-0.2034*T/Tc+0.0013*Pc_bar+0.1714*w+(1/M)**0.3539-0.007*mu)
    else:
        # Eq 3
        k = 0.5147*(-0.2537*T/Tc+0.0017*Pc_bar+0.1501*w+(1/M)**-0.2999)
    return unidades.ThermalConductivity(k)


@refDoc(__doi__, [1])
def ThL_SatoRiedel(T, Tc, M, Tb):
    r"""Calculate the thermal conductivity of a liquid using the Sato-Riedel
    correlation, as explain in [1]_.

    .. math::
        k = \frac{1.1053152}{\sqrt{MW}}\frac{3+20(1-T_r)^{2/3}}
        {3+20(1-T_{br})^{2/3}}

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tc : float
        Critical temperature, [K]
    M : float
        Molecular weight, [g/mol]
    Tb : float
        Boiling temperature, [K]

    Returns
    -------
    k : float
        Liquid thermal conductivity, [W/m·k]

    Examples
    --------
    Example 10-10 from 1_ (4th Edition); CCl4 at 293K

    >>> "%0.3f" % ThL_SatoRiedel(293, 556.4, 153.823, 349.9)
    '0.101'
    """
    Tr = T/Tc
    Tbr = Tb/Tc
    k = 1.1053152/M**0.5*(3+20*(1-Tr)**(2/3))/(3+20*(1-Tbr)**(2/3))
    return unidades.ThermalConductivity(k)


@refDoc(__doi__, [15, 5])
def ThL_Pachaiyappan(T, Tc, M, rho, branched=True):
    r"""Calculates the thermal conductivity of liquid using the Pachaiyappan
    correlation as explain in 5_, procedure 12A1.2, pag 1141

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
    Example in 5_, n-butylbenzene at 140ºF

    >>> T = unidades.Temperature(140, "F")
    >>> Tc = unidades.Temperature(729.32, "F")
    >>> rho = unidades.Density(53.76, "lbft3")
    >>> k = ThL_Pachaiyappan(T, Tc, 134.22, rho)
    >>> "%0.4f" % k.BtuhftF
    '0.0673'
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


@refDoc(__doi__, [16, 5])
def ThL_KanitkarThodos(T, P, Tc, Pc, Vc, M, rho):
    r"""Calculates the thermal conductivity of liquid using the Kanitkar-Thodos
    correlation as explain in 5_, procedure 12A1.3, pag 1143

    .. math::
        \kappa\lambda = -1.884e-6P_r^2 + 1.442e-3P_r +
        \alpha\exp\left(\beta\rho_r\right)

    .. math::
        \alpha = \frac{7.137e-3}{\beta^{3.322}}

    .. math::
        \beta = 0.4 + \frac{0.986}{\exp{0.58\lambda}}

    .. math::
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
        Critical specific volume, [m³/kg]
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
    Example in 5_, n-heptane at 320ºF and 197.4 atm

    >>> T = unidades.Temperature(320, "F")
    >>> P = unidades.Pressure(197.4, "atm")
    >>> Tc = unidades.Temperature(512.69, "F")
    >>> Pc = unidades.Pressure(397.41, "psi")
    >>> Vc = unidades.SpecificVolume(0.0684, "ft3lb")
    >>> rho = unidades.Density(37.93, "lbft3")
    >>> k = ThL_KanitkarThodos(T, P, Tc, Pc, Vc, 100.2, rho)
    >>> "%0.5f" % k.BtuhftF
    '0.06957'
    """

    Pc_atm = Pc/101325
    Tc_R = unidades.K2R(Tc)
    Pr = P/Pc
    rhor = rho*Vc

    li = Tc_R**(1/6)*M**0.5/Pc_atm**(2/3)
    b = 0.4 + 0.986/exp(0.58*li)                                        # Eq 8
    alfa = 7.137e-3/b**3.322                                            # Eq 7
    k = (-1.884e-6*Pr**2+1.442e-3*Pr+alfa*exp(b*rhor))/li               # Eq 6
    return unidades.ThermalConductivity(k, "BtuhftF")


@refDoc(__doi__, [59, 1, 5])
def ThL_Lenoir(T, P, Tc, Pc, ko, To=None, Po=None):
    r"""Calculates the thermal conductivity of liquid using the Lenoir
    correlation as explain in 5_, procedure 12A4.1, pag 1156

    .. math::
        k_2 = k_1\frac{C_2}{C_1}

    .. math::
        C = 17.77+0.065*P_r-7.764*T_r-\frac{2.065T_r^2}{exp(0.2P_r}

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
    ko : float
        Reference thermal conductivity [W/m·k]
    To : float, optional
        Temperature with known thermal conductivity, default T, [K]
    Po : float, optional
        Pressure with known thermal conductivity, default 101325, [Pa]

    Returns
    -------
    k : float
        Thermal conductivity [W/m·k]

    Notes
    -----
    Raise :class:`NotImplementedError` if input pair isn't in limit:

        * 0.4 ≤ Tr ≤ 0.8
        * P ≥ 500 psi

    Examples
    --------
    Example in 5_, toluene at 87.5ºF and 22.044psi

    >>> T = unidades.Temperature(87.5, "F")
    >>> P = unidades.Pressure(22044, "psi")
    >>> Tc = unidades.Temperature(1065.22, "R")
    >>> Pc = unidades.Pressure(595.9, "psi")
    >>> Pref = unidades.Pressure(14.696, "psi")
    >>> ko = unidades.ThermalConductivity(0.07425, "BtuhftF")
    >>> k = ThL_Lenoir(T, P, Tc, Pc, ko, T, Pref)
    >>> "%0.5f" % k.BtuhftF
    '0.09074'

    Example 10-12 in 1_, NO2 at 311K and 276bar

    >>> "%0.2f" % ThL_Lenoir(311, 276e5, 431.35, 101.33e5, 0.124, 311, 2.1e5)
    '0.13'
    """
    # Definiton of optional values
    if To is None:
        To = T
    if Po is None:
        Po = 101325

    Tr2 = T/Tc
    Pr2 = P/Pc
    Tr1 = To/Tc
    Pr1 = Po/Pc

    Pmin = unidades.Pressure(500, "psi")
    if Tr2 < 0.4 or Tr2 > 0.8 or P < Pmin:
        raise NotImplementedError("Input out of bound")

    C1 = 17.77 + 0.065*Pr1 - 7.764*Tr1 - 2.054*Tr1**2/exp(0.2*Pr1)
    C2 = 17.77 + 0.065*Pr2 - 7.764*Tr2 - 2.054*Tr2**2/exp(0.2*Pr2)
    return unidades.ThermalConductivity(ko*C2/C1)


@refDoc(__doi__, [1])
def ThL_Missenard(T, P, Tc, Pc, ko):
    r"""Calculates the thermal conductivity of liquid using the Missenard
    correlation, as explain in 1_

    .. math::
        \frac{\lambda}{\lambda_o} = 1+QP_r^{0.7}

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
    ko : float
        Thermal conductivity at low pressure, [W/m·k]

    Returns
    -------
    k : float
        Thermal conductivity, [W/m/K]

    Notes
    -----
    Raise :class:`NotImplementedError` if input pair isn't in limit:

        * 0.5 ≤ Tr ≤ 0.8
        * 1 ≤ Pr ≤ 200

    Examples
    --------
    Example 10-13 from 1_; toluene at 6330bar and 304K

    >>> "%0.3f" % ThL_Missenard(304, 6330e5, 591.75, 41.08e5, 0.129)
    '0.220'
    """
    Tr = T/Tc
    Pr = P/Pc

    # Check range of validity
    if Tr < 0.5 or Tr > 0.8 or Pr < 1 or Pr > 200:
        raise NotImplementedError("Input out of bound")

    # Interpolate over table to get the Q parameter
    Tri = [0.8, 0.7, 0.6, 0.5]
    Pri = [1, 5, 10, 50, 100, 200]
    Qi = array([[0.036, 0.038, 0.038, 0.038, 0.038, 0.038],
                [0.018, 0.025, 0.027, 0.031, 0.032, 0.032],
                [0.015, 0.020, 0.022, 0.024, 0.025, 0.025],
                [0.012, 0.0165, 0.017, 0.019, 0.020, 0.020]])
    f_Q = interp2d(Pri, Tri, Qi)
    Q = f_Q(Pr, Tr)

    k = ko*(1 + Q*Pr**0.7)
    return unidades.ThermalConductivity(k)


# Gas Thermal conductivity
@refDoc(__doi__, [25, 5])
def ThG_MisicThodos(T, Tc, Pc, M, Cp):
    r"""Calculates thermal conductivity of gas hydrocarbon at low pressure
    using the Misic-Thodos correlation, also referenced in API Procedure
    12B1.2 pag.1162

    .. math::
        \kappa = 1.188e^{-3}\frac{T_rC_p}{\lambda}   for Tr ≤ 1

    .. math::
        \kappa = 2.67e^{-4}\left(14.52T_r-5.14\right)^{2/3}
        \frac{C_p}{\lambda}   for Tr ≤ 1

    .. math::
        \lambda=\frac{T_{c}^{1/6}}{M^{1/2}P_{c}^{2/3}}

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
    Cp : float
        Isobaric heat capacity, [cal/gK]

    Returns
    -------
    k : float
        Thermal conductivity, [cal/scmK]

    Notes
    -----
    Range of validity:
        P ≤ 50 psi

    Examples
    --------
    Example in 5_, 2-methylbutane at 212ºF and 1 atm

    >>> T = unidades.Temperature(212, "F")
    >>> Tc = unidades.Temperature(369.1, "F")
    >>> Pc = unidades.Pressure(498.38, "psi")
    >>> cp = unidades.SpecificHeat(34.49/72.15, "BtulbF")
    >>> "%0.3f" % ThG_MisicThodos(T, Tc, Pc, 72.15, cp).BtuhftF
    '0.013'
    """
    Pc_atm = Pc/101325
    Tr = T/Tc
    cp = unidades.MolarSpecificHeat(Cp*M).calmolK

    li = Tc**(1/6)*M**0.5/Pc_atm**(2/3)
    if Tr < 1:
        k = 0.445e-5*Tr*cp/li                                           # Eq 6
    else:
        k = 1e-6*(14.52*Tr-5.14)**(2/3)*cp/li                           # Eq 7
    return unidades.ThermalConductivity(k, "calscmK")


@refDoc(__doi__, [11])
def ThG_RiaziFaghri(T, Tb, SG):
    r"""Calculates thermal conductivity of gas hydrocarbon at low pressure
    using the Riazi-Faghri correlation.

    .. math::
        \kappa = aT_{b}^{b}SG^{c}

    .. math::
        a = \exp\left(-4.5093-0.6844t-0.1305t^{2}\right)

    .. math::
        b = 0.3003+0.0918t+0.0195t^{2}

    .. math::
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

        * 150ºF ≤ T ≤ 550ºF
        * 0.65 ≤ SG ≤ 0.9
    """
    # Convert input Tb in Kelvin to Fahrenheit to use in the correlation
    Tb_R = unidades.K2R(Tb)
    t = unidades.K2F(T)/100
    A = exp(21.78-8.07986*t+1.12981*t**2-0.05309*t**3)                 # Eq 16
    B = -4.13948+1.29924*t-0.17813*t**2+0.00833*t**3                   # Eq 17
    C = 0.19876-0.0313*t-0.00567*t**2                                  # Eq 18
    k = A*Tb_R**B*SG**C                                                # Eq 7
    return unidades.ThermalConductivity(k, "BtuhftF")


@refDoc(__doi__, [1])
def ThG_Eucken(M, Cv, mu):
    r"""Calculates thermal conductivity of gas al low pressure using the Eucken
    correlation as explain in 1_

    .. math::
        \frac{\lambda M}{\eta C_v} = 1 + \frac{9/4}{C_v/R}

    Parameters
    ----------
    M : float
        Molecular weight, [g/mol]
    Cv : float
        Ideal gas heat capacity at constant volume, [J/kg·K]
    mu : float
        Gas viscosity [Pa·s]

    Returns
    -------
    k : float
        Thermal conductivity, [W/m·K]

    Examples
    --------
    Example 10-1 from 1_; 2-methylbutane at 100ºC and 1bar

    >>> cv_mass = 135.8/72.151*1000
    >>> "%0.4f" % ThG_Eucken(72.151, cv_mass, 8.72e-6)
    '0.0187'
    """
    Cvm = Cv*M/1000
    M = M/1000.
    k = (1 + 9/4/(Cvm/R))*mu*Cvm/M
    return unidades.ThermalConductivity(k)


@refDoc(__doi__, [1])
def ThG_EuckenMod(M, Cv, mu):
    r"""Calculates thermal conductivity of gas al low pressure using the
    modified Eucken correlation as explain in 1_

    .. math::
        \frac{\lambda M}{\eta C_v} = 1.32 + \frac{1.77}{C_v/R}

    Parameters
    ----------
    M : float
        Molecular weight, [g/mol]
    Cv : float
        Ideal gas heat capacity at constant volume, [J/kg·K]
    mu : float
        Gas viscosity [Pa·s]

    Returns
    -------
    k : float
        Thermal conductivity, [W/m/K]

    Examples
    --------
    Example 10-1 from 1_; 2-methylbutane at 100ºC and 1bar

    >>> "%0.4f" % ThG_EuckenMod(72.151, 135.8/72.151*1000, 8.72e-6)
    '0.0234'
    """
    Cvm = Cv*M/1000
    M = M/1000.
    k = (1.32 + 1.77/(Cvm/R))*mu*Cvm/M
    return unidades.ThermalConductivity(k)


@refDoc(__doi__, [49, 1])
def ThG_Chung(T, Tc, M, w, Cv, mu):
    r"""Calculate thermal conductivity of gas at low pressure using the Chung
    correlation

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
    M : float
        Molecular weight, [g/mol]
    Tc : float
        Critical temperature, [K]
    w : float
        Acentric factor, [-]
    Cv : float
        Ideal gas heat capacity at constant volume, [J/kg·K]
    mu : float
        Gas viscosity [Pa·s]

    Returns
    -------
    k : float
        Thermal conductivity, [W/m·K]

    Examples
    --------
    Example 10-1 from 1_; 2-methylbutane at 100ºC and 1bar

    >>> cv_mass = 135.8/72.151*1000
    >>> "%0.4f" % ThG_Chung(373.15, 460.39, 72.151, 0.272, cv_mass, 8.72e-6)
    '0.0229'
    """
    Tr = T/Tc
    Cvm = Cv*M/1000

    alpha = Cvm/R - 1.5
    beta = 0.7862 - 0.7109*w + 1.3168*w**2
    Z = 2 + 10.5*Tr**2
    phi = 1 + alpha*((0.215 + 0.28288*alpha - 1.061*beta + 0.26665*Z)/(
        0.6366 + beta*Z + 1.061*alpha*beta))

    # Eq 9
    k = 7.452*mu*10/M*phi   # Viscosity in P
    return unidades.ThermalConductivity(k, "calscmK")


@refDoc(__doi__, [5])
def ThG_NonHydrocarbon(T, P, id):
    r"""Calculates thermal conductivity of selected nonhydrocarbon, referenced
    in API procedure 12C1.1, pag 1174

    .. math::
        \kappa = A + BT + CT^2 + DP + E\frac{P}{T^{1.2}} + \frac{F}{
        \left(0.4P-0.001T\right)^{0.015}} + G\ln{P}

    Parameters
    ----------
    T : float
        Temperature, [K]
    P : float
        Pressure, [Pa]
    id : int
        Index of compound in database

    Returns
    -------
    k : float
        Thermal conductivity, [Btu/hftºF]

    Notes
    -----
    This method calculate the thermal conductivity of selected nonhydrocarbon
    gases, the available compounds are:
        * 1   -   Hydrogen
        * 46  -   Nitrogen
        * 47  -   Oxygen
        * 48  -   Carbon Monoxide
        * 50  -   Hydrogen Sulfide
        * 51  -   Sulfur dioxide
        * 111 -   Sulfur trioxide

    The range of validity of relation depends of compounds, it's checked in
    procedure and raise a :class:`NotImplementedError` when inputs are out of
    bound or the id of compound isn't supported:

        * N2, CO     - 150ºR ≤ T ≤ 2460ºR, 15psi ≤ P ≤ 10000psi
        * O2         - 150ºR ≤ T ≤ 2460ºR, 15psi ≤ P ≤ 15000psi
        * H2         - 260ºR ≤ T ≤ 2260ºR, 15psi ≤ P ≤ 10000psi
        * SO2        - 960ºR ≤ T ≤ 2460ºR, 15psi ≤ P ≤ 10000psi
        * H2S, SO3   - 460ºR ≤ T ≤ 2460ºR, P atmospheric

    Examples
    --------
    Example from 5_; oxygen at 984.67ºR and 6075psi

    >>> T = unidades.Temperature(984.67, "R")
    >>> P = unidades.Pressure(6075, "psi")
    >>> "%0.5f" % ThG_NonHydrocarbon(T, P, 47).BtuhftF
    '0.03265'
    """
    # Table 12C1.2
    dat = {
        1: (4.681e-3, 2.e-4, -3.6e-8, 0.0, 0.0, 0.0, 1.7e-3),
        46: (4.561e-3, 1.61e-5, 0.0, 2.56e-9, 5.299e-3, 2.47e-3, 0.0),
        47: (5.95e-4, 1.71e-5, 0.0, -2.10e-8, 5.869e-3, 6.995e-3, 0.0),
        48: (1.757e-3, 1.55e-5, 0.0, 2.08e-8, 5.751e-3, 5.6e-3, 0.0),
        50: (-1.51e-3, 2.25e-5, 3.32e-10, 0.0, 0.0, 0.0, 0.0),
        51: (2.5826e-2, 1.35e-5, 0.0, -4.4e-7, 1.026e-2, -2.631e-2, 0.0),
        111: (-1.02e-3, 1.35e-5, 4.17e-9, 0.0, 0.0, 0.0, 0.0)}

    # Check supported compounds
    if id not in dat:
        raise NotImplementedError("Compound don't supported")

    # Convert input T in Kelvin to Rankine to use in the correlation
    t = unidades.K2R(T)
    p = unidades.Pressure(P).psi

    # Check input parameter
    if id == 1:
        tmin = 260
        tmax = 2260
        pmin = 15
        pmax = 10000
    elif id in (46, 47):
        tmin = 460
        tmax = 2460
        pmin = 15
        pmax = 10000
    elif id == 48:
        tmin = 460
        tmax = 2460
        pmin = 15
        pmax = 15000
    elif id == 50:
        tmin = 960
        tmax = 2460
        pmin = 15
        pmax = 10000
    elif id in (51, 111):
        tmin = 460
        tmax = 2460
        pmin = 1
        pmax = 100
    if t < tmin or t > tmax or p < pmin or p > pmax:
        raise NotImplementedError("Input out of bound")

    A, B, C, D, E, F, G = dat[id]
    k = A + B*t + C*t**2 + D*p + E*p/t**1.2 + F/(.4*p-.001*t)**.015 + G*log(p)
    return unidades.ThermalConductivity(k, "BtuhftF")


@refDoc(__doi__, [58, 1])
def ThG_StielThodos(T, Tc, Pc, Vc, M, V, ko):
    """Calculate thermal conductivity of compressed gases using the
    Stiel-Thodos correlation.

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tc : float
        Critical temperature, [K]
    Pc : float
        Critical pressure, [Pa]
    Vc : float
        Critical specific volume, [m³/kg]
    M : float
        Molecular weight, [g/mol]
    Vm : float
        Volume at T and P [m³/kg]
    kg : float
        Low-pressure gas thermal conductivity [W/m·K]

    Returns
    -------
    ko : float
        High-pressure gas thermal conductivity [W/m·K]

    Examples
    --------
    Example 10-3 from 1_; nitrous oxide at 105ºC and 138bar

    >>> T = unidades.Temperature(105, "C")
    >>> Vc = 97/44.013/1000
    >>> V = 144/44.013/1000
    >>> "%0.4f" % ThG_StielThodos(T, 309.6, 72.55e5, Vc, 44.013, V, 0.0234)
    '0.0414'
    """
    # Calcualte the Zc internally
    Zc = Pc*1e-3*Vc*M/Tc/R

    # The thermal conductivity in the paper define in cal/s·cm·K
    ko = unidades.ThermalConductivity(ko).calscmK

    rhor = Vc/V
    gamma = (Tc*M**3/(Pc/101325)**4)**(1/6)

    if rhor < 0.5:
        lr = 14*(exp(0.535*rhor) - 1)                                   # Eq 3
    elif rhor < 2:
        lr = 13.1*(exp(0.67*rhor) - 1.069)                              # Eq 4
    else:
        lr = 2.976*(exp(1.155*rhor) + 2.016)                            # Eq 5

    k = ko + lr*1e-8/Zc**5/gamma
    return unidades.ThermalConductivity(k, "calscmK")


@refDoc(__doi__, [49, 1])
def ThG_P_Chung(T, Tc, Vc, M, w, D, k, rho, ko):
    r"""Calculate the thermal conductivity of a compressed gas using the Chung
    correlation

    .. math::
        \lambda = \frac{31.2 \eta^o\Psi}{M}(1/G_2+B_6y)+qB_7y^2T_r^{1/2}G_2

    .. math::
        y = \frac{\rho V_c}{6}

    .. math::
        G_1 = \frac{1-0.5y}{\left(1-y\right)^3

    .. math::
        G_2 = \frac{E_1\left(\left(1-\exp\left(-E_4y\right)\right)/y\right) +
        E_2G_1\expE_5y + E_3G_1}{E_1E_4+E_2+E_3}

    .. math::
        \mu^{**} = E_7y^2G_2\exp(E_8+\frac{E_9}{T^*}+\frac{E_{10}}{T^*^2}

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tc : float
        Critical temperature, [K]
    Vc : float
        Critical volume of the gas [m^3/mol]
    M : float
        Molecular weight, [g/mol]
    w : float
        Acentric factor, [-]
    D : float
        Dipole moment, [Debye]
    k : float, optional
        Corection factor for polar substances, [-]
    rho : float
        Density, [kg/m³]
    ko : float
        Low-pressure gas thermal conductivity[Pa*S]

    Returns
    -------
    k : float
        High-pressure gas thermal conductivity [W/m·K]

    Examples
    --------
    Example 10-4 in 1_, propylene at 473K and 150bar

    >>> Vc = 184.6/42.081/1e3
    >>> rho = 1/172.1*42.081*1e3
    >>> th = ThG_P_Chung(473, 364.9, Vc, 42.081, 0.142, 0.4, 0, rho, 0.0389)
    >>> "%0.3f" % th
    '0.062'
    """
    # units in molar base
    Vc = Vc*M*1000
    rho = rho/M/1000

    # Thermal conductivity in procedure in cal/s·cm·K
    ko = unidades.ThermalConductivity(ko).calscmK

    Tr = T/Tc
    mur = 131.3*D/(Vc*Tc)**0.5                                          # Eq 8

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
        B.append(bo + b1*w + b2*mur**4 + b3*k)
    B1, B2, B3, B4, B5, B6, B7 = B

    Y = rho*Vc/6
    G1 = (1-0.5*Y)/(1-Y)**3
    H2 = (B1*((1-exp(-B4*Y))/Y)+B2*G1*exp(B5*Y)+B3*G1)/(B1*B4+B2+B3)

    # Eq 12
    kk = ko*(1/H2 + B6*Y)
    kp = (3.039e-4*(Tc/M)**0.5/Vc**(2/3))*B7*Y**2*H2*Tr**0.5

    return unidades.ThermalConductivity(kk+kp, "calscmK")


@refDoc(__doi__, [61, 1])
def ThG_TRAPP(T, Tc, Vc, Zc, M, w, rho, ko):
    """Calculate the thermal conductivity of a compressed gas using the TRAPP
    (TRAnsport Property Prediction) method.

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tc : float
        Critical temperature, [K]
    Vc : float
        Critical volume, [m³/kg]
    Zc : float
        Critical pressure, [Pa]
    M : float
        Molecular weight, [g/mol]
    w : float
        Acentric factor, [-]
    rho : float
        Density, [kg/m3]
    ko : float
        Low-pressure gas thermal conductivity, [Pa*S]

    Returns
    -------
    k : float
        High-pressure gas thermal conductivity [W/m·K]

    Examples
    --------
    Example 9-13 in 1_, isobutane at 500K and 100bar

    >>> Vc = 184.6/42.081/1000
    >>> rho = 1/172.1*42.081*1000
    >>> "%0.3f" % ThG_TRAPP(473, 364.9, Vc, 0.2798, 42.081, 0.142, rho, 0.0389)
    '0.061'
    """
    # Reference fluid properties, propane
    TcR = 369.83
    rhocR = 1/200  # mol/cm³
    ZcR = 0.276
    wR = 0.152

    # Convert volume to molar base
    Vc = Vc*M*1000
    rho = rho/M/1000
    rhoc = 1/Vc

    Tr = T/Tc

    f = Tc/TcR*(1+(w-wR)*(0.05203-0.7498*log(Tr)))
    h = rhocR/rhoc*ZcR/Zc*(1-(w-wR)*(0.1436-0.2822*log(Tr)))

    To = T/f
    rho0 = rho*h
    Fl = (44.094/M*f)**0.5/h**(2/3)
    Xl = (1+2.1866*(w-wR)/(1-0.505*(w-wR)))**0.5

    # Coefficients in [53]_, pag 796
    # Density are in mol/dm³
    rho0 *= 1000
    rhocR *= 1000

    rhorR = rho0/rhocR
    TrR = To/TcR
    lR = 15.2583985944*rhorR + 5.29917319127*rhorR**3 + \
        (-3.05330414748+0.450477583739/TrR)*rhorR**4 + \
        (1.03144050679-0.185480417707/TrR)*rhorR**5

    return unidades.ThermalConductivity(Fl*Xl*lR*1e-3 + ko)


# Liquid surface tension
def Tension_Parametric(T, args, Tc):
    r"""Calculates surface tension of fluid using a parametric equation

    .. math::
        $\sigma=A\left(1-T_{r}\right)^{B}$

    .. math::
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


@refDoc(__doi__, [39, 40, 1])
def Tension_BlockBird(T, Tc, Pc, Tb):
    r"""Calculates surface tension of liquid using the Block-Bird correlation
    using the Miller expression for α.

    .. math::
        \frac{\sigma}{P_c^{2/3}T_c^{1/3}} = \left(0.132\alpha_c-0.279\right)
        \left(1-T_r\right)^{11/9}

    .. math::
        \alpha_c = 0.9076\left(1+\frac{T_br\ln{P_c}}{1-T_{br}}\right)

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
    sigma : float
        Surface tension, [N/m]

    Examples
    --------
    Example 12.2 from 1_; ethyl mercaptan at 303K

    >>> "%0.1f" % Tension_BlockBird(303, 499, 54.9e5, 308.15).dyncm
    '22.4'
    """
    Tr = T/Tc
    Trb = Tb/Tc
    Pc_atm = Pc/101325
    Pc_bar = Pc*1e-5

    # Eq 9 in 40_
    alfa = 0.9076*(1+Trb*log(Pc_atm)/(1-Trb))

    # Eq 11 in 39_, Coefficient from 5_
    sr = (0.132*alfa-0.279)*(1-Tr)**(11/9)

    # Eq 5 in 39_
    sigma = sr*Pc_bar**(2/3)*Tc**(1/3)
    return unidades.Tension(sigma, "dyncm")


@refDoc(__doi__, [1])
def Tension_Pitzer(T, Tc, Pc, w):
    r"""Calculates surface tension of liquid using the Pitzer correlation as
    explain in 1_

    .. math::
        \sigma = P_c^{2/3}T_c^{1/3}\frac{1.86 + 1.18\omega}{19.05}
        \left(\frac{3.75 + 0.91\omega}{0.291 - 0.08\omega}\right)^{2/3}
        (1 - T_r)^{11/9}

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
    sigma : float
        Liquid surface tension, [N/m]

    Examples
    --------
    Example 12.2 from 1_; ethyl mercaptan at 303K

    >>> "%0.1f" % Tension_Pitzer(303, 499, 54.9e5, 0.192).dyncm
    '23.5'
    """
    Tr = T/Tc
    Pc_bar = Pc*1e-5

    sigma = Pc_bar**(2/3)*Tc**(1/3)*(1.86+1.18*w)/19.05 * \
        ((3.75+0.91*w)/(0.291-0.08*w))**(2/3)*(1-Tr)**(11/9)
    return unidades.Tension(sigma, "dyncm")


@refDoc(__doi__, [41])
def Tension_ZuoStenby(T, Tc, Pc, w):
    r"""Calculates surface tension of liquid using the Zuo-Stenby correlation

    .. math::
        \sigma_r = \sigma_r^{(1)}+ \frac{\omega - \omega^{(1)}}
        {\omega^{(2)}-\omega^{(1)}} \left(\sigma_r^{(2)}-\sigma_r^{(1)}\right)

    .. math::
        \sigma_r = \ln{\left(\frac{\sigma}{T_c^{1/3}P_c^{2/3}} + 1\right)}

    .. math::
        \sigma^{(1)} = 40.520(1-T_r)^{1.287}
        \sigma^{(2)} = 52.095(1-T_r)^{1.21548}

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
    sigma : float
        Liquid surface tension, [N/m]

    Examples
    --------
    Example 12.2 from 1_; ethyl mercaptan at 303K
    The procedure use the critical properties from meos library, something
    diferent than Poling values, so the last decimal isn't exact

    >>> "%0.0f" % Tension_ZuoStenby(303, 499, 54.9e5, 0.192).dyncm
    '23'
    """
    from lib.mEoS import CH4, nC8

    Tr = T/Tc
    Pc_bar = Pc*1e-5

    s1 = 40.520*(1 - Tr)**1.287                               # Eq 7, methane
    s2 = 52.095*(1 - Tr)**1.21548                             # Eq 8, n-octane

    # Eq 2 for reference fluid
    sr1 = log(s1/(CH4.Tc**(1/3)*CH4.Pc.bar**(2/3)) + 1)
    sr2 = log(s2/(nC8.Tc**(1/3)*nC8.Pc.bar**(2/3)) + 1)

    # Eq 1
    sr = sr1 + (w-CH4.f_acent)/(nC8.f_acent-CH4.f_acent)*(sr2-sr1)

    # Eq 2 for desired fluid
    sigma = Tc**(1/3)*Pc_bar**(2/3)*(exp(sr)-1)
    return unidades.Tension(sigma, "dyncm")


@refDoc(__doi__, [42])
def Tension_SastriRao(T, Tc, Pc, Tb, alcohol=False, acid=False):
    r"""Calculates surface tension of a liquid using the Sastri-Rao correlation

    .. math::
        \sigma = KT_b^xP_c^yT_{br}^z\left(\frac{T_c-T}{T_c-T_b}\right)^m

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tb : float
        Boiling temperature of the fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]

    Returns
    -------
    sigma : float
        Liquid surface tension, [N/m]

    Examples
    --------
    Example 12.2 from 1_; ethyl mercaptan at 303K

    >>> "%0.2f" % Tension_SastriRao(303, 499, 54.9e5, 308.15).dyncm
    '21.92'

    Selected point in Table 3 of [42]_

    >>> from lib.mEoS import Acetone as Ac
    >>> "%0.2f" % Tension_SastriRao(298.16, Ac.Tc, Ac.Pc, Ac.Tb).dyncm
    '22.36'

    >>> from lib.mEoS import Methanol as Met
    >>> "%0.2f" % Tension_SastriRao(333.16, Met.Tc, Met.Pc, Met.Tb, True).dyncm
    '19.34'
    """
    if alcohol:
        K, x, y, z, m = 2.28, 0.175, 0.25, 0, 0.8
    elif acid:
        K, x, y, z, m = 0.125, 0.35, 0.5, -1.85, 11/9
    else:
        K, x, y, z, m = 0.158, 0.35, 0.5, -1.85, 11/9

    Tbr = Tb/Tc
    Pc_bar = Pc*1e-5

    # Eq 3
    sigma = K * Tb**x * Pc_bar**y * Tbr**z * ((Tc-T)/(Tc-Tb))**m
    return unidades.Tension(sigma, "dyncm")


@refDoc(__doi__, [43, 42])
def Tension_Hakim(T, Tc, Pc, w, X):
    r"""Calculates surface tension of a liquid using the Hakim-Steinberg-Stiel
    correlation

    .. math::
        \sigma = P_c^{2/3}T_c^{1/3} \sigma_{r|T_r=0.6}
        \left(\frac{1-T_r}{0.4}\right)^m

    .. math::
        \sigma_{r|T_r=0.6} = 0.1574 + 0.359\omega - 1.769X - 13.69X^2 -
        0.51\omega^2 + 1.298\omega X

    .. math::
        m = 1.21+0.5385\omega-14.61X-32.07X^2-1.65\omega^2+22.03X\omega

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
    X : float
        Stiel Polar Factor, [-]

    Returns
    -------
    sigma : float
        Liquid surface tension, [N/m]

    Examples
    --------
    Selected point in Table 5 of [42]_

    >>> from lib.mEoS import Methanol as Me
    >>> "%0.1f" % Tension_Hakim(313.16, Me.Tc, Me.Pc, Me.f_acent, 0.037).dyncm
    '20.4'
    """
    Tr = T/Tc
    Pc_atm = Pc/101325

    # Eq 6
    sr06 = 0.1574 + 0.359*w - 1.769*X - 13.69*X**2 - 0.510*w**2 + 1.298*X*w
    # Eq 9
    m = 1.210 + 0.5385*w - 14.61*X - 32.07*X**2 - 1.656*w**2 + 22.03*X*w

    # Eq 8
    sigma = Pc_atm**(2/3)*Tc**(1/3)*sr06*((1-Tr)/0.4)**m

    return unidades.Tension(sigma, "dyncm")


@refDoc(__doi__, [44])
def Tension_Miqueu(T, Tc, Vc, M, w):
    r"""Calculates surface tension of a liquid using the Miqueu et al.
    correlation

    .. math::
        \sigma = kT_c\left(\frac{N_A}{V_c}\right)^{2/3}
        (4.35+4.14\omega)t^{1.26}(1+0.19t^{0.5}-0.487t)

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tc : float
        Critical temperature, [K]
    Vc : float
        Critical volume, [m^3/kg]
    M : float
        Molecular weight, [g/mol]
    w : float
        Acentric factor, [-]

    Returns
    -------
    sigma : float
        Liquid surface tension, [N/m]
    """
    t = 1 - T/Tc

    # Eq 13
    sigma = Boltzmann * Tc * (Avogadro/Vc/1000/M)**(2/3) * (4.35+4.14*w) * \
        t**1.26 * (1+0.19*t**0.5-0.25*t)
    return unidades.Tension(sigma, "mNm")


@refDoc(__doi__, [1])
def CpL_Poling(T, Tc, w, Cpgo):
    r"""Calculate liquid isobaric heat capacitiy with the CSP method reported
    in [1]_, Eq 6-6.

    .. math::
        \frac{C_p-C_p^o}{R} = 1.586 + \frac{0.49}{1-T_r} + \omega\left(4.2775 +
        \frac{6.3\left(1-T_r\right)^{1/3}}{T_r} + \frac{0.4355}{1-T_r}\right)

    Parameters
    ----------
    T : float
        Temperature, [K]
    Tc : float
        Critical temperature, [K]
    w : float
        Acentric factor, [-]
    Cpgo : float
        Isobaric ideal gas heat capacity, [J/mol/K]

    Returns
    -------
    Cplm : float
        Liquid constant-pressure heat capacitiy, [J/mol/K]

    Notes
    -----
    This correlation fail with associating compound

    Examples
    --------

    Example 6-3 from [1]_, possible bug in reference

    >>> "%0.1f" % CpL_Poling(350.0, 435.5, 0.203, 91.21)
    '143.8'
    """
    Tr = T/Tc
    Cpr = 1.586 + 0.49/(1-Tr) + \
        w*(4.2775 + 6.3*(1-Tr)**(1/3)/Tr + 0.4355/(1-Tr))
    return Cpgo + R*Cpr


# Acentric factor
@refDoc(__doi__, [4])
def facent_LeeKesler(Tb, Tc, Pc):
    r"""Calculates acentric factor of a fluid using the Lee-Kesler correlation

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
    """
    Tr = Tb/Tc
    Pr = 101325/Pc
    w = (log(Pr) - 5.92714 + 6.09648/Tr + 1.28862*log(Tr) - 0.169347*Tr**6)/(
        15.2518 - 15.6875/Tr - 13.4721*log(Tr) + 0.43577*Tr**6)

    return unidades.Dimensionless(w)


@refDoc(__doi__, [19])
def prop_Edmister(**kwargs):
    """Calculate the missing parameters between Tc, Pc, Tb and acentric factor
    from the Edmister (1958) correlations

    Parameters
    ----------
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


@refDoc(__doi__, [8, 1])
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


# Other properties
@refDoc(__doi__, [17, 18, 5])
def Vc_Riedel(Tc, Pc, w, M):
    """Calculates critical volume of a fluid using the Riedel correlation
    as explain in 5_, procedure 4A3.1, Pag. 302

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
        Critical volume, [m³/kg]

    Examples
    --------
    Example in 5_, n-nonane

    >>> Tc = unidades.Temperature(610.68, "F")
    >>> Pc = unidades.Pressure(331.8, "psi")
    >>> "%0.3f" % Vc_Riedel(Tc, Pc, 0.4368, 128.2551).ft3lb
    '0.068'
    """
    # Eq 1 in [18]
    alfa = 5.811 + 4.919*w

    Vc = R*1000*Tc/Pc/(3.72+0.26*(alfa-7))/M
    return unidades.SpecificVolume(Vc, "lg")


@refDoc(__doi__, [23])
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
    """
    Zra = 0.29056 - 0.08775*w                                           # Eq 1
    return Zra


@refDoc(__doi__, [60, 5])
def Henry(T, args):
    r"""Calculates Henry constant for gases in liquid at low pressure, also
    referenced in API procedure 9A7.1, pag 927

    .. math::
        lnH = A/T + BlnT + CT + D

    Parameters
    ----------
    T : float
        Temperature, [K]
    args : list
        Coefficients for equation

    Returns
    -------
    H : float
        Henry constant, [psi/xmole]

    Notes
    -----
    The parameters for several compound are in database:
        Hydrogen, Helium, Argon, Neon, Krypton, Xenon, Oxygen, Nitrogen,
        Hydrogen sulfide, Carbon monoxide, Carbon dioxide, Sulfur dioxide,
        Nitrous oxide, Chlorine,Bromine, Iodine, Methane, Ethane, Propane,
        Ethylene, Ammonia.
    The Henry constant is returned as unidades.Pressure instance

    Examples
    --------
    Example from 5_; Hydrogen sulfide in water at 77ºF

    >>> T = unidades.Temperature(77, "F")
    >>> "%0.0f" % Henry(T, [-65864.7, -215.127, 0.185874, 1384.15]).psi
    '8257'
    """
    T_R = unidades.K2R(T)
    B1, B2, B3, B4 = args
    H = exp(B1/T_R + B2*log(T_R) + B3*T_R + B4)
    return unidades.Pressure(H, "psi")


class Componente(object):
    """Class to define a chemical compound from the database

    Parameters
    ----------
    id : int
        index of compound in database

    Notes
    -----
    Additionally can define custom calculation method with the parameters:

        * *rhoL*: Liquid density correlation index
        * *RhoLP*: Compressed liquid density correlation index
        * *MuL*: Liquid viscosity correlation index
        * *MuG*: Gas viscosity correlation index
        * *MuGP*: Compressed gas viscosity correlation index
        * *ThCondL*: Liquid thermal conductivity correlation index
        * *ThCondG*: Gas thermal conductivity correlation index
        * *ThCondPG*: Compressed gas thermal conductivity correlation index
        * *Tension*: Surface tension correlation index
        * *facent*: Acentric factor correlation index for missing cases
        * *Pv*: Vapor pressure correlation index

    This option overwrite the project configuration and the user configuration,
    for now only in API usage. Not custom stream property definition in main
    program

    Examples
    --------
    These are several examples of usage of this class with several configuration
    definition, obviously not all correlation return valid values.

    Surface tension methods: Example 12.2 from 1_; ethyl mercaptan at 303K

    >>> cmp1 = Componente(137, Tension=0)
    >>> cmp2 = Componente(137, Tension=2)
    >>> cmp3 = Componente(137, Tension=3)
    >>> cmp4 = Componente(137, Tension=4)
    >>> cmp5 = Componente(137, Tension=5)
    >>> cmp6 = Componente(137, Tension=6)
    >>> "%0.2f %0.2f" % (cmp1.Tension(303).dyncm, cmp2.Tension(303).dyncm)
    '22.23 22.41'
    >>> "%0.2f %0.2f" % (cmp3.Tension(303).dyncm, cmp4.Tension(303).dyncm)
    '23.45 22.65'
    >>> "%0.2f %0.2f" % (cmp5.Tension(303).dyncm, cmp6.Tension(303).dyncm)
    '21.94 22.29'

    Gas viscosity methods: Example 9-1 in 1_, SO2 at 300ºC

    >>> c0 = Componente(51, MuG=0)
    >>> c1 = Componente(51, MuG=1)
    >>> c2 = Componente(51, MuG=2)
    >>> c3 = Componente(51, MuG=3)
    >>> c4 = Componente(51, MuG=4)
    >>> c5 = Componente(51, MuG=5)
    >>> args = (unidades.Temperature(300, "C"), 101325, 0)
    >>> "%0.2f %0.2f" % (c0.Mu_Gas(*args).microP, c1.Mu_Gas(*args).microP)
    '241.51 243.88'
    >>> "%0.2f %0.2f" % (c2.Mu_Gas(*args).microP, c3.Mu_Gas(*args).microP)
    '246.25 250.60'
    >>> "%0.2f %0.2f" % (c4.Mu_Gas(*args).microP, c5.Mu_Gas(*args).microP)
    '243.01 241.95'

    Example 9-9 in 1_, n-pentane at 500K and 101bar

    >>> c0 = Componente(8, MuGP=0)
    >>> c1 = Componente(8, MuGP=6)
    >>> args = (500, 101e5, 0)
    >>> "%0.2f %0.2f" % (c0.Mu_Gas(*args).microP, c1.Mu_Gas(*args).microP)
    '534.84 520.54'

    Example 9-12 in 1_, ammonia at 520K and 600bar

    >>> c0 = Componente(63, MuGP=0)
    >>> c1 = Componente(63, MuGP=1)
    >>> c2 = Componente(63, MuGP=2)
    >>> c3 = Componente(63, MuGP=3)
    >>> c5 = Componente(63, MuGP=5)
    >>> c7 = Componente(63, MuGP=7)
    >>> args = (520, 600e5, 1/48.2*17.031*1000)
    >>> "%0.2f %0.2f" % (c0.Mu_Gas(*args).microP, c1.Mu_Gas(*args).microP)
    '496.31 457.79'
    >>> "%0.2f %0.2f" % (c2.Mu_Gas(*args).microP, c3.Mu_Gas(*args).microP)
    '454.02 506.49'
    >>> "%0.2f %0.2f" % (c5.Mu_Gas(*args).microP, c7.Mu_Gas(*args).microP)
    '454.29 508.31'

    Liquid viscosity correlations: Example 9.196 from 1_; toluene at 383K

    >>> c0 = Componente(40, MuL=0)
    >>> c1 = Componente(40, MuL=1)
    >>> c2 = Componente(40, MuL=2)
    >>> c3 = Componente(40, MuL=3)
    >>> args = (383, 101325)
    >>> "%0.3f %0.3f" % (c0.Mu_Liquido(*args).cP, c1.Mu_Liquido(*args).cP)
    '0.239 0.233'
    >>> "%0.3f %0.3f" % (c2.Mu_Liquido(*args).cP, c3.Mu_Liquido(*args).cP)
    '0.220 0.291'

    Example 9.15 from 1_, methylcyclohexane at 300K 500 bar

    >>> c0 = Componente(39, MuLP=0)
    >>> c2 = Componente(39, MuLP=2)
    >>> args = (300, 500e5)
    >>> "%0.3f %0.3f" % (c0.Mu_Liquido(*args).cP, c2.Mu_Liquido(*args).cP)
    '1.045 1.052'

    Vapor pressure correlations: Example 7-2 from 1_; ethylbenzene at 347.25K

    >>> c0 = Componente(45, Pv=0)
    >>> c2 = Componente(45, Pv=2)
    >>> c3 = Componente(45, Pv=3)
    >>> c4 = Componente(45, Pv=4)
    >>> c5 = Componente(45, Pv=5)
    >>> c6 = Componente(45, Pv=6)
    >>> c7 = Componente(45, Pv=7)
    >>> "%0.3f %0.3f" % (c0.Pv(347.25).kPa, c2.Pv(347.25).kPa)
    '13.342 13.340'
    >>> "%0.3f %0.3f" % (c3.Pv(347.25).kPa, c4.Pv(347.25).kPa)
    '13.311 12.839'
    >>> "%0.3f %0.3f" % (c5.Pv(347.25).kPa, c6.Pv(347.25).kPa)
    '13.104 13.393'
    >>> "%0.3f" % c7.Pv(347.25).kPa
    '13.550'

    Liquid density correlations: Example from 5_; propane at 30ºF

    >>> c0 = Componente(4, RhoL=0)
    >>> c1 = Componente(4, RhoL=1)
    >>> c2 = Componente(4, RhoL=2)
    >>> c3 = Componente(4, RhoL=3)
    >>> c4 = Componente(4, RhoL=4)
    >>> c5 = Componente(4, RhoL=5)
    >>> c6 = Componente(4, RhoL=6)
    >>> c7 = Componente(4, RhoL=7)
    >>> c8 = Componente(4, RhoL=8)
    >>> c9 = Componente(4, RhoL=9)
    >>> args = (unidades.Temperature(30, "F"), 1e5)
    >>> "%0.3f %0.3f" % (c0.RhoL(*args).kgl, c1.RhoL(*args).kgl)
    '0.532 0.532'
    >>> "%0.3f %0.3f" % (c2.RhoL(*args).kgl, c3.RhoL(*args).kgl)
    '0.539 0.528'
    >>> "%0.3f %0.3f" % (c4.RhoL(*args).kgl, c5.RhoL(*args).kgl)
    '0.525 0.530'
    >>> "%0.3f %0.3f" % (c6.RhoL(*args).kgl, c7.RhoL(*args).kgl)
    '0.531 0.529'
    >>> "%0.3f %0.3f" % (c8.RhoL(*args).kgl, c9.RhoL(*args).kgl)
    '0.529 0.525'

    Example from 5_; n-octane at 212ºF and 4410 psi

    >>> T = unidades.Temperature(212, "F")
    >>> P = unidades.Pressure(4410, "psi")
    >>> c0 = Componente(12, RhoLP=0)
    >>> c1 = Componente(12, RhoLP=1)
    >>> c2 = Componente(12, RhoLP=2)
    >>> c3 = Componente(12, RhoLP=3)
    >>> c4 = Componente(12, RhoLP=4)
    >>> c5 = Componente(12, RhoLP=5)
    >>> "%0.3f %0.3f" % (c0.RhoL(T, P).kgl, c1.RhoL(T, P).kgl)
    '0.676 0.705'
    >>> "%0.3f %0.3f" % (c2.RhoL(T, P).kgl, c3.RhoL(T, P).kgl)
    '0.672 0.675'
    >>> "%0.3f %0.3f" % (c4.RhoL(T, P).kgl, c5.RhoL(T, P).kgl)
    '0.678 0.919'

    Liquid thermal conductivity: Example in 5_, n-butylbenzene at 140ºF

    >>> c0 = Componente(78, ThCondL=0)
    >>> c1 = Componente(78, ThCondL=1)
    >>> c2 = Componente(78, ThCondL=2)
    >>> c3 = Componente(78, ThCondL=3)
    >>> c4 = Componente(78, ThCondL=4)
    >>> c5 = Componente(78, ThCondL=5)
    >>> c6 = Componente(78, ThCondL=6)
    >>> c7 = Componente(78, ThCondL=7)
    >>> T = unidades.Temperature(140, "F")
    >>> rho = c0.RhoL(T, 101325)
    >>> args = (unidades.Temperature(140, "F"), 101325, rho)
    >>> "%0.3f %0.3f" % (c0.ThCond_Liquido(*args), c1.ThCond_Liquido(*args))
    '0.120 0.112'
    >>> "%0.3f %0.3f" % (c2.ThCond_Liquido(*args), c3.ThCond_Liquido(*args))
    '0.122 0.145'
    >>> "%0.3f %0.3f" % (c4.ThCond_Liquido(*args), c5.ThCond_Liquido(*args))
    '0.125 0.126'
    >>> "%0.3f %0.3f" % (c6.ThCond_Liquido(*args), c7.ThCond_Liquido(*args))
    '0.112 0.115'

    Example 10-13 from 1_; toluene at 6330bar and 304K

    >>> c0 = Componente(41, ThCondLP=0)
    >>> c1 = Componente(41, ThCondLP=1)
    >>> c2 = Componente(41, ThCondLP=2)
    >>> rho = c0.RhoL(304, 6330e5)
    >>> args = (304, 6330e5, rho)
    >>> "%0.3f %0.3f" % (c0.ThCond_Liquido(*args), c1.ThCond_Liquido(*args))
    '0.356 0.235'
    >>> "%0.3f" % c2.ThCond_Liquido(*args)
    '0.223'

    Gas thermal conductivity: Example in 5_, 2-methylbutane at 212ºF and 1atm

    >>> c0 = Componente(7, ThCondG=0)
    >>> c1 = Componente(7, ThCondG=1)
    >>> c2 = Componente(7, ThCondG=2)
    >>> c3 = Componente(7, ThCondG=3)
    >>> c4 = Componente(7, ThCondG=4)
    >>> c5 = Componente(7, ThCondG=5)
    >>> T = unidades.Temperature(212, "F")
    >>> rho = c0.RhoL(T, 101325)
    >>> args = (T, 101325, rho)
    >>> "%0.3f %0.3f" % (c0.ThCond_Gas(*args), c1.ThCond_Gas(*args))
    '0.022 0.023'
    >>> "%0.3f %0.3f" % (c2.ThCond_Gas(*args), c3.ThCond_Gas(*args))
    '0.024 0.019'
    >>> "%0.3f %0.3f" % (c4.ThCond_Gas(*args), c5.ThCond_Gas(*args))
    '0.024 0.029'

    Example 10-3 from 1_; nitrous oxide at 105ºC and 138bar

    >>> c0 = Componente(110, ThCondGP=0)
    >>> c1 = Componente(110, ThCondGP=1)
    >>> c2 = Componente(110, ThCondGP=2)
    >>> T = unidades.Temperature(105, "C")
    >>> rho = 1/144*c0.M*1000
    >>> args = (T, 138e5, rho)
    >>> "%0.4f %0.4f" % (c0.ThCond_Gas(*args), c1.ThCond_Gas(*args))
    '0.0415 0.0396'
    >>> "%0.4f" % c2.ThCond_Gas(*args)
    '0.0406'
    """

    _bool = False
    kwargs = {"RhoL": None,
              "RhoLP": None,
              "MuL": None,
              "MuLP": None,
              "MuG": None,
              "MuGP": None,
              "ThCondL": None,
              "ThCondLP": None,
              "ThCondG": None,
              "ThCondGP": None,
              "Pv": None,
              "facent": None,
              "Tension": None}

    METHODS_RhoL = ["DIPPR", "Rackett", "Cavett", "COSTALD",
                    "Yen-Woods (1966)", "Yamada-Gun (1973)", "Bhirud (1978)",
                    "Mchaweh (2004)", "Riedel",
                    "Chueh-Prausnitz (1967)"]
    METHODS_RhoLP = ["Tait-COSTALD (1982)", "Chang-Zhao (1990)",
                     "Aalto-Keskinen (1996)", "Aalto-Keskinen (1999)",
                     "Nasrifar (2000)", "API"]

    METHODS_MuL = ["DIPPR", "Parametric", "Letsou-Stiel (1973)",
                   "Przedziecki-Sridhar (1985)"]
    METHODS_MuLP = ["Lucas (1981)", "API", "Kouzel"]
    METHODS_MuG = ["DIPPR", "Chapman-Enskog", "Chung (1988)", "Lucas (1981)",
                   "Stiel-Thodos (1961)", "Gharagheizi (2012)",
                   "Yoon-Thodos (1970)"]
    METHODS_MuGP = ["Lucas", "Chung (1988)", "Brulé", "Jossi", "TRAPP",
                    "Stiel-Thodos", "Reichenberg", "Dean-Stiel", "API"]

    METHODS_ThG = ["DIPPR", "Misic-Thodos", "Chung (1988)", "Eucken",
                   "Modified Eucken", "Riazi-Faghri"]
    METHODS_ThGP = ["Stiel-Thodos", "Chung (1988)", "TRAPP"]
    METHODS_ThL = ["DIPPR", "Pachaiyappan", "Sato-Riedel", "Kanitkar-Thodos",
                   "Riazi-Faghri", "Gharagheizi", "Lakshmi-Prasad", "Nicola"]
    METHODS_ThLP = ["Kanitkar-Thodos", "Lenoir", "Missenard"]

    METHODS_Pv = ["DIPPR", "Wagner", "Antoine", "Ambrose-Walton", "Lee-Kesler",
                  "Riedel", "Sanjari", "Maxwell-Bonnel"]
    METHODS_facent = ["Lee-Kesler", "Edmister", "Ambrose-Walton"]
    METHODS_Tension = ["DIPPR", "Parametric", "Block-Bird", "Pitzer",
                       "Zuo-Stenby", "Sastri-Rao", "Hakim", "Miqueu"]

    def __init__(self, id=None, **kwargs):
        if not id:
            return

        # FIXME DATABASE: Meanwhile check type here

        self._bool = True
        self.id = id
        self.kwargs = Componente.kwargs.copy()
        self.kwargs.update(kwargs)
        self.Config = config.getMainWindowConfig()
        cmp = sql.getElement(id)
        self.formula = cmp[1]
        self.name = cmp[2]
        self.M = cmp[3]
        self.SG = cmp[124]
        self.Tc = unidades.Temperature(cmp[4])
        self.Pc = unidades.Pressure(cmp[5], "atm")
        self.Tb = unidades.Temperature(cmp[129])
        self.Tf = unidades.Temperature(cmp[130])
        if cmp[125] != 0:
            self.f_acent = float(cmp[125])
        elif self.Pc and self.Tc and self.Tb:
            self.f_acent = self._f_acent()
        else:
            self.f_acent = 0
        if cmp[6] != 0:
            # In database the value is in molar base m³/kmol, convert to m³/kg
            self.Vc = unidades.SpecificVolume(cmp[6]/self.M)
        elif self.f_acent != 0 and self.Tc != 0 and self.Pc != 0:
            self.Vc = Vc_Riedel(self.Tc, self.Pc, self.f_acent, self.M)
        else:
            self.Vc = 0

        if self.Vc:
            self.rhoc = 1/self.Vc
        else:
            self.rhoc = 0

        if self.Tc:
            self.Zc = self.Pc*self.Vc*self.M/1000/R/self.Tc
        else:
            self.Zc = 0

        if cmp[7] != 0:
            self.API = cmp[7]
        elif cmp[124] != 0:
            self.API = 141.5/cmp[124]-131.5
        else:
            self.API = 0
        self.cp = cmp[8:14]

        # Parametric parameters
        self.antoine = cmp[14:17]
        self.antoine += cmp[150:154]
        self.wagner = cmp[154:158]
        self._parametricMu = cmp[21:23]
        self._parametricSigma = cmp[23:25]
        self.henry = cmp[17:21]

        # DIPPR parameters
        self._dipprRhoS = cmp[25:33]
        self._dipprRhoL = cmp[33:41]
        self._dipprPv = cmp[41:49]
        self._dipprHv = cmp[49:57]
        self._dipprCpS = cmp[57:65]
        self._dipprCpL = cmp[65:73]
        self._dipprCpG = cmp[73:81]
        self._dipprMuL = cmp[81:89]
        self._dipprMuG = cmp[89:97]
        self._dipprKL = cmp[97:105]
        self._dipprKG = cmp[105:113]
        self._dipprSigma = cmp[113:121]

        self.dipole = unidades.DipoleMoment(cmp[121])
        if cmp[123] != 0.0:
            self.rackett = cmp[123]
        else:
            self.rackett = Rackett(self.f_acent)

        self.Vliq = cmp[122]

        self.SolubilityParameter = unidades.SolubilityParameter(cmp[126])
        self.Kw = cmp[127]
        self.stiel = cmp[128]
        self.CASNumber = str(cmp[131])
        self.alternateFormula = str(cmp[132])
        if cmp[133]:
            self.UNIFAC = eval(cmp[133])
        else:
            self.UNIFAC = []

        self.Dm = cmp[134]
        self.ek = cmp[135]

        self.UNIQUAC_area = cmp[136]
        self.UNIQUAC_volumen = cmp[137]
        self.f_acent_mod = cmp[138]
        self.Hf = unidades.Enthalpy(cmp[139]/self.M)
        self.Gf = unidades.Enthalpy(cmp[140]/self.M)
        self.wilson = unidades.MolarVolume(cmp[141], "ccmol")
        self.NetHeating = unidades.Enthalpy(cmp[142]/self.M)
        self.GrossHeating = unidades.Enthalpy(cmp[143]/self.M)
        self.Synonyms = str(cmp[144])
        self.V_char = cmp[145]
        self.calor_formacion_solido = cmp[146]
        self.energia_formacion_solido = cmp[147]
        self.PolarParameter = cmp[148]
        self.smile = str(cmp[149])

        # Desglosar formula en elementos y átomos de cada elemento
        decmp = atomic_decomposition(self.formula)
        atoms = sum([val for val in decmp.values()])
        self.C = decmp.get("C", 0)
        self.H = decmp.get("H", 0)
        self.O = decmp.get("O", 0)  # noqa
        self.N = decmp.get("N", 0)
        self.S = decmp.get("S", 0)

        if self.C and self.H:
            self.HC = self.H/self.C

        # Chemical class definition
        # Some procedures need the chemical type of compound. This could be
        # hardcoded in database but it try to autodetect below

        # Hydrocarbon definition: Compound only formed by carbon or hydrogen
        if self.C+self.H == atoms:
            self.isHydrocarbon = True
        else:
            self.isHydrocarbon = False

        # Organic compound from name termination
        if self.O >= 1 and self.name[:-2] == "ol":
            self.isAlcohol = True
        else:
            self.isAlcohol = False

        if self.O >= 2 and self.name[:-4].lower() == "acid":
            self.isAcid = True
        else:
            self.isAcid = False

        # TODO: Add other chemical types, this let use other advenced group
        # contribution correlation
        # self.isAldehyde = False
        # self.isAmine = False
        # self.isEster = False
        # self.isEther= False
        # self.isKetone = False
        # self.isHalogen = False
        # self.isCiclico = False

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

    def __bool__(self):
        return self._bool

    # Calculation of undefined properties of compound
    def _f_acent(self):
        """Acentric factor calculation in compounds with undefined property"""
        method = self.kwargs["facent"]
        if method is None or method >= len(Componente.METHODS_facent):
            method = self.Config.getint("Transport", "f_acent")

        if method == 0:
            return facent_LeeKesler(self.Tb, self.Tc, self.Pc)
        elif method == 1:
            return prop_Edmister(Tb=self.Tb, Tc=self.Tc, Pc=self.Pc)
        elif method == 2:
            Pvr = self.Pv(0.7*self.Tc)/self.Pc
            return facent_AmbroseWalton(Pvr)

    @refDoc(__doi__, [47], tab=8)
    def _MuCritical(self):
        """Critical viscosity calculation procedure
        """
        # Define the critical viscosity, use the Table 11A5.4 in 5_
        if 1 < self.id <= 21 and self.id != 7:
            muc = [0, 0, 0.014, 0.02, 0.0237, 0.027, 0.0245, 0, 0.0255, 0.0350,
                   0.0264, 0.0273, 0.0282, 0.0291, 0.0305, 0.0309, 0.0315,
                   0.0328, 0.0337, 0.0348, 0.0355, 0.0362][self.id]
        elif id == 90:
            muc = 0.0370
        elif id == 91:
            muc = 0.0375
        elif id == 92:
            muc = 0.0388
        else:
            xi = self.Tc**(1/6)/self.M**0.5/self.Pc.atm**(2/3)
            muc = 7.7e-4/xi

        return unidades.Viscosity(muc, "cP")

    # Ideal properties
    @refDoc(__doi__, [5], tab=8)
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
        """
        A, B, C, D, E, F = self.cp
        cp = A + B*T + C*T**2 + D*T**3 + E*T**4 + F*T**5
        return unidades.SpecificHeat(cp/self.M, "calgK")

    @refDoc(__doi__, [5], tab=8)
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
        """
        To = 298.15
        A, B, C, D, E, F = self.cp
        H = B*T + C/2*T**2 + D/3*T**3 + E/4*T**4 + F/5*T**5
        Ho = B*To + C/2*To**2 + D/3*To**3 + E/4*To**4 + F/5*To**5
        return unidades.Enthalpy((H-Ho)/self.M, "calg")

    @refDoc(__doi__, [5], tab=8)
    def _so(self, T):
        r"""Ideal gas entropy calculation from polinomial coefficient of
        specific heat saved in database
        Coefficient in database are in the form [A,B,C,D,E,F]
        Explained in procedure 7A1.1, pag 543

        .. math::
            So = A \ln T + BT + C/2T^2 + D/3T^3 + E/4T^4 + F/5T^5

        Parameters
        ----------
        T : float
            Temperature, [K]

        Notes
        -----
        The units in the calculate ideal enthalpy are in cal/mol·K, the
        reference state is set to T=298.15K
        """
        A, B, C, D, E, F = self.cp
        so = A*log(T) + B*T + C/2*T**2 + D/3*T**3 + E/4*T**4 + F/5*T**5
        return unidades.SpecificHeat(so/self.M, "calgK")

    # Physical properties
    def RhoS(self, T):
        """Calculate the density of solid phase using the DIPPR equations"""
        return DIPPR("rhoS", T, self._dipprRhoS[:-2], M=self.M, Tc=self.Tc)

    def RhoL(self, T, P):
        """Calculate the density of liquid phase using any of available
        correlation"""
        method = self.kwargs["RhoL"]
        if method is None or method >= len(Componente.METHODS_RhoL):
            method = self.Config.getint("Transport", "RhoL")
        Pcorr = self.kwargs["RhoLP"]
        if Pcorr is None or method >= len(Componente.METHODS_RhoLP):
            Pcorr = self.Config.getint("Transport", "Corr_RhoL")

        if T > self.Tc:
            T = 0.9*self.Tc

        # Calculate of low pressure viscosity
        if method == 0 and self._dipprRhoL:
            rhos = DIPPR("rhoL", T, self._dipprRhoL[:-2], M=self.M, Tc=self.Tc)
        elif method == 1 and self.rackett != 0 and T < self.Tc:
            rhos = RhoL_Rackett(T, self.Tc, self.Pc, self.rackett, self.M)
        elif method == 2 and self.Vliq != 0:
            rhos = RhoL_Cavett(T, self.Tc, self.M, self.Vliq)
        elif method == 3:
            if self.f_acent_mod != 0:
                w = self.f_acent_mod
            else:
                w = self.f_acent
            rhos = RhoL_Costald(T, self.Tc, w, self.Vc)
        elif method == 4:
            rhos = RhoL_YenWoods(T, self.Tc, self.Vc, self.Zc)
        elif method == 5:
            rhos = RhoL_YamadaGunn(T, self.Tc, self.Pc, self.f_acent, self.M)
        elif method == 6:
            rhos = RhoL_Bhirud(T, self.Tc, self.Pc, self.f_acent, self.M)
        elif method == 7:
            d = Mchaweh_d.get(self.id, 0)/100
            rhos = RhoL_Mchaweh(T, self.Tc, self.Vc, self.f_acent, d)
        elif method == 8:
            rhos = RhoL_Riedel(T, self.Tc, self.Vc, self.f_acent)
        elif method == 9:
            rhos = RhoL_ChuehPrausnitz(T, self.Tc, self.Vc, self.f_acent)
        else:
            if self._dipprRhoL and \
                    self._dipprRhoL[6] <= T <= self._dipprRhoL[7]:
                rhos = DIPPR(
                    "rhoL", T, self._dipprRhoL[:-2], M=self.M, Tc=self.Tc)
            elif self.rackett != 0 and T < self.Tc:
                rhos = RhoL_Rackett(
                    T, self.Tc, self.Pc, self.rackett, self.M)
            elif self.Vliq != 0:
                rhos = RhoL_Cavett(T, self.Tc, self.M, self.Vliq)
            else:
                if self.f_acent_mod != 0:
                    w = self.f_acent_mod
                else:
                    w = self.f_acent
                rhos = RhoL_Costald(T, self.Tc, w, self.Vc)

        # Add correction factor for high pressure
        if P < 1e6:
            rho = rhos
        elif Pcorr == 0:
            if self.f_acent_mod:
                w = self.f_acent_mod
            else:
                w = self.f_acent
            Ps = self.Pv(T)
            rho = RhoL_TaitCostald(
                T, P, self.Tc, self.Pc, w, Ps, rhos)
        elif Pcorr == 1:
            if self.f_acent_mod:
                w = self.f_acent_mod
            else:
                w = self.f_acent
            Ps = self.Pv(T)
            rho = RhoL_ChangZhao(T, P, self.Tc, self.Pc, w, Ps, rhos)
        elif Pcorr == 2:
            if self.f_acent_mod:
                w = self.f_acent_mod
            else:
                w = self.f_acent
            Ps = self.Pv(T)
            rho = RhoL_AaltoKeskinen(T, P, self.Tc, self.Pc, w, Ps, rhos)
        elif Pcorr == 3:
            if self.f_acent_mod:
                w = self.f_acent_mod
            else:
                w = self.f_acent
            Ps = self.Pv(T)
            rho = RhoL_AaltoKeskinen2(T, P, self.Tc, self.Pc, w, Ps, rhos)
        elif Pcorr == 4:
            if self.f_acent_mod:
                w = self.f_acent_mod
            else:
                w = self.f_acent
            Ps = self.Pv(T)
            rho = RhoL_Nasrifar(T, P, self.Tc, self.Pc, w, self.M, Ps, rhos)
        else:
            rho = RhoL_API(T, P, self.Tc, self.Pc, self.SG, rhos)

        return rho

    def Pv(self, T):
        """Vapor pressure calculation procedure using the method defined in
        preferences"""
        method = self.kwargs["Pv"]
        if method is None or method >= len(Componente.METHODS_Pv):
            method = self.Config.getint("Transport", "Pv")

        if method == 0 and self._dipprPv and \
                self._dipprPv[6] <= T <= self._dipprPv[7]:
            return DIPPR("Pv", T, self._dipprPv[:-2], M=self.M, Tc=self.Tc)
        elif method == 1 and self.wagner[0]:
            return Pv_Wagner(T, self.Tc, self.Pc, self.wagner)
        elif method == 2 and self.antoine[0]:
            return Pv_Antoine(T, self.antoine, Tc=self.Tc)
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
                return DIPPR("Pv", T, self._dipprPv[:-2], M=self.M, Tc=self.Tc)
            elif self.wagner[0]:
                return Pv_Wagner(T, self.Tc, self.Pc, self.wagner)
            elif self.antoine[0]:
                return Pv_Antoine(T, self.antoine, Tc=self.Tc)
            elif self.Pc and self.Tc and self.f_acent:
                return Pv_AmbroseWalton(T, self.Tc, self.Pc, self.f_acent)
            elif self.Pc and self.Tc and self.Tb:
                return Pv_Riedel(T, self.Tc, self.Pc, self.Tb)
            elif self.Kw and self.Tb:
                return Pv_MaxwellBonnel(T, self.Tb, self.Kw)

    def ThCond_Liquido(self, T, P, rho):
        """Liquid thermal conductivity procedure using the method defined in
        preferences, use the decision diagram in 5_ Figure 12-0.2 pag 1135"""
        method = self.kwargs["ThCondL"]
        if method is None or method >= len(Componente.METHODS_ThL):
            method = self.Config.getint("Transport", "ThCondL")
        Pcorr = self.kwargs["ThCondLP"]
        if Pcorr is None or method >= len(Componente.METHODS_ThLP):
            Pcorr = self.Config.getint("Transport", "Corr_ThCondL")

        if T > self.Tc:
            T = 0.9*self.Tc

        # Calculate of low pressure viscosity
        if method == 0 and self._dipprKL:
            ko = DIPPR("kL", T, self._dipprKL[:-2], M=self.M, Tc=self.Tc)
        elif method == 1 and T < self.Tc:
            ko = ThL_Pachaiyappan(T, self.Tc, self.M, rho, self.branched)
        elif method == 2 and self.Tb:
            ko = ThL_SatoRiedel(T, self.Tc, self.M, self.Tb)
        elif method == 3:
            ko = ThL_KanitkarThodos(
                T, P, self.Tc, self.Pc, self.Vc, self.M, rho)
        elif method == 4 and self.Tb and self.SG:
            ko = ThL_RiaziFaghri(T, self.Tb, self.SG)
        elif method == 5 and self.Tb:
            ko = ThL_Gharagheizi(T, self.Pc, self.Tb, self.M, self.f_acent)
        elif method == 6:
            ko = ThL_LakshmiPrasad(T, self.M)
        elif method == 7:
            ko = ThL_Nicola(
                T, self.M, self.Tc, self.Pc, self.f_acent, self.dipole.Debye)
        else:
            if self._dipprKL and self._dipprKL[6] <= T <= self._dipprKL[7]:
                ko = DIPPR("kL", T, self._dipprKL[:-2], M=self.M, Tc=self.Tc)
            elif T < self.Tc:
                ko = ThL_Pachaiyappan(T, self.Tc, self.M, rho, self.branched)
            elif self.Tb:
                ko = ThL_SatoRiedel(T, self.Tc, self.M, self.Tb)
            else:
                ko = ThL_KanitkarThodos(
                    T, P, self.Tc, self.Pc, self.Vc, self.M, rho)

        # Add correction factor for high pressure
        if P < 1e6:
            k = ko
        elif Pcorr == 0:
            k = ThL_KanitkarThodos(
                T, P, self.Tc, self.Pc, self.Vc, self.M, rho)
        elif Pcorr == 1:
            k = ThL_Lenoir(T, P, self.Tc, self.Pc, ko)
        elif Pcorr == 2:
            k = ThL_Missenard(T, P, self.Tc, self.Pc, ko)

        return k

    def ThCond_Gas(self, T, P, rho):
        """Vapor thermal conductivity calculation procedure using the method
        defined in preferences, decision diagram in API Databook, pag. 1136"""
        method = self.kwargs["ThCondG"]
        if method is None or method >= len(Componente.METHODS_ThG):
            method = self.Config.getint("Transport", "ThCondG")
        Pcorr = self.kwargs["ThCondGP"]
        if Pcorr is None or method >= len(Componente.METHODS_ThGP):
            Pcorr = self.Config.getint("Transport", "Corr_ThCondG")

        # Calculate of low pressure viscosity
        if method == 0 and self._dipprKG:
            ko = DIPPR("kG", T, self._dipprKG[:-2], M=self.M, Tc=self.Tc)
        elif method == 1:
            cp = self.Cp_Gas_DIPPR(T)
            ko = ThG_MisicThodos(T, self.Tc, self.Pc, self.M, cp)
        elif method == 2:
            cv = self.Cv(T)
            muo = self.Mu_Gas(T, 101325, rho)
            ko = ThG_Chung(T, self.Tc, self.M, self.f_acent, cv, muo)
        elif method == 3:
            cv = self.Cv(T)
            muo = self.Mu_Gas(T, 101325, rho)
            ko = ThG_Eucken(self.M, cv, muo)
        elif method == 4:
            cv = self.Cv(T)
            muo = self.Mu_Gas(T, 101325, rho)
            ko = ThG_EuckenMod(self.M, cv, muo)
        elif method == 5 and self.SG and self.Tb:
            ko = ThG_RiaziFaghri(T, self.Tb, self.SG)
        else:
            if self._dipprKG and self._dipprKG[6] <= T <= self._dipprKG[7]:
                ko = DIPPR("kG", T, self._dipprKG[:-2], M=self.M, Tc=self.Tc)
            else:
                cp = self.Cp_Gas_DIPPR(T)
                ko = ThG_MisicThodos(T, self.Tc, self.Pc, self.M, cp)

        # Add correction factor for high pressure
        if P < 1e7:
            k = ko
        elif self.id in [1, 46, 47, 48, 50, 51, 111]:
            k = ThG_NonHydrocarbon(T, P, self.id)
        elif Pcorr == 0:
            k = ThG_StielThodos(
                T, self.Tc, self.Pc, self.Vc, self.M, 1/rho, ko)
        elif Pcorr == 1:
            K = self._K_Chung()
            k = ThG_P_Chung(T, self.Tc, self.Vc, self.M, self.f_acent,
                            self.dipole.Debye, K, rho, ko)
        elif Pcorr == 2:
            k = ThG_TRAPP(
                T, self.Tc, self.Vc, self.Zc, self.M, self.f_acent, rho, ko)
        return k

    def Mu_Gas(self, T, P, rho):
        """Vapor viscosity calculation procedure using the method defined in
        preferences, decision diagram in API Databook, pag. 1026"""
        method = self.kwargs["MuG"]
        if method is None or method >= len(Componente.METHODS_MuG):
            method = self.Config.getint("Transport", "MuG")
        Pcorr = self.kwargs["MuGP"]
        if Pcorr is None or method >= len(Componente.METHODS_MuGP):
            Pcorr = self.Config.getint("Transport", "Corr_MuG")

        # Calculate of low pressure viscosity
        if method == 0 and self._dipprMuG and \
                self._dipprMuG[6] <= T <= self._dipprMuG[7]:
            muo = DIPPR("muG", T, self._dipprMuG[:-2], M=self.M, Tc=self.Tc)
        elif method == 1 and self.Dm and self.ek:
            omega = self._Collision(T)
            muo = MuG_ChapmanEnskog(T, self.M, self.Dm, omega)
        elif method == 2:
            k = self._K_Chung()
            muo = MuG_Chung(T, self.Tc, self.Vc, self.M, self.f_acent,
                            self.dipole.Debye, k)
        elif method == 3:
            muo = MuG_Lucas(
                T, P, self.Tc, self.Pc, self.Zc, self.M, self.dipole.Debye)
        elif method == 4:
            muo = MuG_StielThodos(T, self.Tc, self.Pc, self.M)
        elif method == 5:
            muo = MuG_Gharagheizi(T, self.Tc, self.Pc, self.M)
        elif method == 6:
            muo = MuG_YoonThodos(T, self.Tc, self.Pc, self.M)
        else:
            if self._dipprMuG and self._dipprMuG[6] <= T <= self._dipprMuG[7]:
                muo = DIPPR(
                    "muG", T, self._dipprMuG[:-2], M=self.M, Tc=self.Tc)
            elif self.Dm and self.ek:
                omega = self._Collision(T)
                muo = MuG_ChapmanEnskog(T, self.M, self.Dm, omega)
            else:
                muo = MuG_StielThodos(T, self.Tc, self.Pc, self.M)

        # Add correction factor for high pressure
        if P < 0.6*self.Pc:
            mu = muo
        elif Pcorr == 0:
            mu = MuG_Lucas(
                T, P, self.Tc, self.Pc, self.Zc, self.M, self.dipole.Debye)
        elif Pcorr == 1:
            k = self._K_Chung()
            mu = MuG_P_Chung(T, self.Tc, self.Vc, self.M, self.f_acent,
                             self.dipole.Debye, k, rho, muo)
        elif Pcorr == 2:
            mu = MuG_Brule(T, self.Tc, self.Vc, self.M, self.f_acent, rho, muo)
        elif Pcorr == 3:
            mu = MuG_Jossi(self.Tc, self.Pc, self.rhoc, self.M, rho, muo)
        elif Pcorr == 4:
            mu = MuG_TRAPP(
                T, self.Tc, self.Vc, self.Zc, self.M, self.f_acent, rho, muo)
        elif Pcorr == 5:
            mu = MuG_P_StielThodos(
                self.Tc, self.Pc, self.rhoc, self.M, rho, muo)
        elif Pcorr == 6:
            mu = MuG_Reichenberg(T, P, self.Tc, self.Pc, self.Vc, self.M,
                                 self.dipole.Debye, muo)
        elif Pcorr == 7:
            mu = MuG_DeanStiel(self.Tc, self.Pc, self.rhoc, self.M, rho, muo)
        elif Pcorr == 8:
            mu = MuG_API(T, P, self.Tc, self.Pc, muo)
        else:
            if self.dipole:
                mu = MuG_Lucas(
                    T, P, self.Tc, self.Pc, self.Zc, self.M, self.dipole.Debye)
            elif self.isHydrocarbon:
                mu = MuG_DeanStiel(
                    self.Tc, self.Pc, self.rhoc, self.M, rho, muo)
            else:
                mu = MuG_TRAPP(T, P, self.Tc, self.Pc, muo)

        return mu

    def _K_Chung(self):
        """Internal procedure to calculate the polar correction factor for
        Chung viscosity correlation

        Chung, T.H., Lee, L.L., Starling, K.E.
        Applications of Kinetic Gas Theories and Multiparameter Correlation for
        Prediction of Dilute Gas Viscosity and Thermal Conductivity
        Ind. Eng. Chem. Fundam. 23(1) (1984) 8-13
        """
        # Table I
        ki = {117: 0.215175,
              134: 0.174823,
              146: 0.143453,
              145: 0.143453,
              160: 0.131671,
              159: 0.131671,
              313: 0.121555,
              335: 0.114230,
              357: 0.108674,
              130: 0.091549,
              62: 0.075908}
        k = ki.get(self.id, 0)

        # Rough alcohol definition using the oxygen atoms count:
        if not k and self.isAlcohol:
            k = 0.0682+0.276659*17*self.O/self.M                       # Eq 10
        return k

    def _Collision(self, T):
        """Internal procudere to calculate the transport collision integral
        necessary for Chapman-Enskog viscosity correlation"""
        T_ = T/self.ek
        omega = Collision_Neufeld(T_)

        # Use the Brokaw Collision integral for polar compounds
        # Brokaw, R.S. Predicting Transport Properties of Dilute Gases. I&EC
        # Process Design and Development 8(22) (1969) 240-253
        if self.PolarParameter:
            omega += 0.2*self.PolarParameter**2/T_
        return omega

    def Mu_Liquido(self, T, P):
        """Liquid viscosity calculation procedure using the method defined in
        preferences, decision diagram in API Databook, pag. 1026"""
        method = self.kwargs["MuL"]
        if method is None or method >= len(Componente.METHODS_MuL):
            method = self.Config.getint("Transport", "MuL")
        Pcorr = self.kwargs["MuLP"]
        if Pcorr is None or method >= len(Componente.METHODS_MuLP):
            Pcorr = self.Config.getint("Transport", "Corr_MuL")

        # Calculate of low pressure viscosity
        if method == 0 and self._dipprMuL and \
                self._dipprMuL[6] <= T <= self._dipprMuL[7]:
            muo = DIPPR("muL", T, self._dipprMuL[:-2], M=self.M, Tc=self.Tc)
        elif method == 1 and self._parametricMu[0]:
            muo = MuL_Parametric(T, self._parametricMu)
        elif method == 2:
            muo = MuL_LetsouStiel(T, self.M, self.Tc, self.Pc, self.f_acent)
        elif method == 3 and self.Tf:
            # Use the critical point as reference volume
            muo = MuL_PrzedzieckiSridhar(
                T, self.Tc, self.Pc, self.Vc, self.f_acent, self.M, self.Tf)
        else:
            if self._dipprMuL and \
                    self._dipprMuL[6] <= T <= self._dipprMuL[7]:
                return DIPPR(
                    "muL", T, self._dipprMuL[:-2], M=self.M, Tc=self.Tc)
            elif self._parametricMu[0]:
                return MuL_Parametric(T, self._parametricMu)
            elif self.Tc and self.Pc and self.f_acent and self.M:
                return MuL_LetsouStiel(
                    T, self.M, self.Tc, self.Pc, self.f_acent)

        # Add correction factor for high pressure
        if P < 0.6*self.Pc:
            mu = muo
        elif Pcorr == 0:
            Ps = self.Pv(T)
            mu = MuL_Lucas(T, P, self.Tc, self.Pc, self.f_acent, Ps, muo)
        elif Pcorr == 1 and self.Pc and self.f_acent:
            muc = self._MuCritical()
            mu = MuL_API(T, P, self.Tc, self.Pc, self.f_acent, muc)
        elif Pcorr == 2:
            mu = MuL_Kouzel(T, P, muo)
        else:
            if self.Pc and self.f_acent:
                Ps = self.Pv(T)
                mu = MuL_Lucas(
                    T, P, self.Tc, self.Pc, self.f_acent, Ps, muo)
            elif self.Tb < 650:
                muc = self._MuCritical()
                mu = MuL_API(
                    T, P, self.Tc, self.Pc, self.f_acent, muc)
            else:
                mu = MuL_Kouzel(T, P, muo)

        return mu

    def Tension(self, T):
        """Liquid surface tension procedure using the method defined in
        preferences"""

        # Calculate this property only if the Temperature is below critical
        # temperature
        if T > self.Tc:
            return 0

        method = self.kwargs["Tension"]
        if method is None or method >= len(Componente.METHODS_Tension):
            method = self.Config.getint("Transport", "Tension")

        if method == 0 and self._dipprSigma and \
                self._dipprSigma[6] <= T <= self._dipprSigma[7]:
            return DIPPR("sigma", T, self._dipprSigma[:-2], M=self.M,
                         Tc=self.Tc)
        elif method == 1 and self._parametricSigma[0]:
            return Tension_Parametric(T, self._parametricSigma, self.Tc)
        elif method == 2 and self.Tb:
            return Tension_BlockBird(T, self.Tc, self.Pc, self.Tb)
        elif method == 3 and self.f_acent:
            return Tension_Pitzer(T, self.Tc, self.Pc, self.f_acent)
        elif method == 4:
            return Tension_ZuoStenby(T, self.Tc, self.Pc, self.f_acent)
        elif method == 5:
            return Tension_SastriRao(T, self.Tc, self.Pc, self.Tb,
                                     alcohol=self.isAlcohol, acid=self.isAcid)
        elif method == 6 and self.stiel:
            return Tension_Hakim(T, self.Tc, self.Pc, self.f_acent, self.stiel)
        elif method == 7 and self.Vc:
            return Tension_Miqueu(T, self.Tc, self.Vc, self.M, self.f_acent)

        else:
            if self._dipprSigma and \
                    self._dipprSigma[6] <= T <= self._dipprSigma[7]:
                return DIPPR("sigma", T, self._dipprSigma[:-2], M=self.M,
                             Tc=self.Tc)
            elif self._parametricSigma[0]:
                return Tension_Parametric(T, self._parametricSigma, self.Tc)
            elif self.stiel:
                return Tension_Hakim(
                    T, self.Tc, self.Pc, self.f_acent, self.stiel)
            elif self.Vc:
                return Tension_Miqueu(
                    T, self.Tc, self.Vc, self.M, self.f_acent)
            elif self.Tb:
                return Tension_BlockBird(T, self.Tc, self.Pc, self.Tb)
            else:
                return Tension_Pitzer(T, self.Tc, self.Pc, self.f_acent)

    def Hv_DIPPR(self, T):
        """Calculate the heat of vaporization using the DIPPR equations"""
        return DIPPR("Hv", T, self._dipprHv[:-2], M=self.M, Tc=self.Tc)

    def Cp_Solido_DIPPR(self, T):
        """Calculate the specific heat of solid using the DIPPR equations"""
        return DIPPR("cpS", T, self._dipprCpS[:-2], M=self.M, Tc=self.Tc)

        return DIPPR("cpL", T, self._dipprCpL[:-2], M=self.M, Tc=self.Tc)

    def Cp_Liquido(self, T):
        if self._dipprCpL:
            return DIPPR("cpL", T, self._dipprCpL[:-2], M=self.M, Tc=self.Tc)
        else:
            Cpo = self._Cpo(T)
            return CpL_Poling(T, self.Tc, self.f_acent, Cpo)

    def Cp_Gas_DIPPR(self, T):
        """Calculate the specific heat of gas using the DIPPR equations"""
        if self._dipprCpG:
            return DIPPR("cpG", T, self._dipprCpG[:-2], M=self.M, Tc=self.Tc)
        else:
            return self._Cpo(T)

    def Cv(self, T):
        """Isochoric specific heat"""
        cp = self.Cp_Gas_DIPPR(T)
        cv = cp/1000-R/self.M
        return unidades.SpecificHeat(cv, "kJkgK")

    def Fase(self, T, P):
        """Método que calcula el estado en el que se encuentra la sustancia"""
        Pv = self.Pv(T)
        if Pv > P:
            return 1
        else:
            return 0
