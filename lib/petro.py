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


from configparser import ConfigParser
import time
import warnings

from PyQt5.QtWidgets import QApplication
from scipy import pi, exp, sqrt, log10, log, sin
from scipy.optimize import fsolve, leastsq, newton
from numpy.linalg import solve
from numpy import array

from lib import unidades
from lib.physics import R_atml, R_Btu
from lib.compuestos import newComponente
from lib.config import conf_dir


__doi__ = {
    1:
        {"autor": "Katz, D. L., and Firoozabadi, A.",
         "title": "Predicting Phase Behavior of Condensate/Crude-oil Systems"
                  "Using Methane Interaction Coefficients",
         "ref": "Journal of Petroleum Technology, Nov. 1978, pp. 1649–1655",
         "doi": "10.2118/6721-pa"},
    2:
        {"autor": "Ahmed, T., G. Cady, and A. Story",
         "title": "A Generalized Correlation for Characterizing the"
                  "Hydrocarbon Heavy Fractions.",
         "ref": "Paper SPE 14266, presented at the 60th Annual Technical"
                "Conference of the Society of Petroleum Engineers, Las Vegas,"
                "September 22–25, 1985.",
         "doi": ""},
    3:
        {"autor": "Riazi, M. R. and Daubert, T. E.",
         "title": "Characterization Parameters for Petroleum Fractions",
         "ref": "Industrial and Engineering Chemistry Research, Vol. 26, 1987,"
                "pp. 755-759.",
         "doi": "10.1021/ie00064a023"},
    4:
        {"autor": "Riazi, M. R., and T. E. Daubert",
         "title": "Simplify Property Predictions",
         "ref": "Hydrocarbon Processing (March 1980): 115–116",
         "doi": ""},
    5:
        {"autor": "Twu, C.H.",
         "title": "An Internally Consistent Correlation for Predicting the"
                  "Critical Properties and Molecular Weights of Petroleum and"
                  "Coal-tar Liquids",
         "ref": "Fluid Phase Equilbria, 16: 137–150.",
         "doi": "10.1016/0378-3812(84)85027-x"},
    6:
        {"autor": "Sim, W.J. and Daubert, T.E.",
         "title": "Prediction of Vapor-Liquid Equilibria of Undefined"
                  "Mixtures",
         "ref": "Industrial and Engineering Chemistry–Process Design and"
                "Development, 19: 386 –393.",
         "doi": "10.1021/i260075a010"},
    7:
        {"autor": "Cavett, R.H.",
         "title": "Physical data for distillation calculations, vapor-liquid"
                  "equilibrium.",
         "ref": "Proceedings of the 27th Meeting, API, San Francisco, Issue 3,"
                "pp. 351–366.",
         "doi": ""},
    8:
        {"autor": "Kesler, M. G., and B. I. Lee",
         "title": "Improve Prediction of Enthalpy of Fractions",
         "ref": "Hydrocarbon Processing Volume 55, Issue 3, March 1976, Pages"
                "153-158.",
         "doi": ""},
    9:
        {"autor": "Tarek Ahmed",
         "title": "Equations of State and PVT Analysis: Applications for"
                  "Improved Reservoir Modeling, 2nd Edition",
         "ref": "Gulf Professional Publishing, 2016, ISBN 9780128015704,",
         "doi": "10.1016/B978-0-12-801570-4.00002-7"},
    10:
        {"autor": "Watansiri, S., Owens, V.H., Starling, K.E.",
         "title": "Correlations for estimating critical constants, acentric"
                  "factor, and dipole moment for undefined coal-fluid"
                  "fractions",
         "ref": "Ind. Eng. Chem. Process. Des. Dev. 24, 294–296",
         "doi": "10.1021/i200029a013"},
    11:
        {"autor": "Rowe, A.M.",
         "title": "Internally Consistent Correlations for Predicting Phase"
                  "Compositions for Use in Reservoir Compositional Simulators",
         "ref": "Paper SPE 7475, In: Presented at the 53rd Annual Society of"
                "Petroleum Engineers Fall Technical Conference and Exhibition,"
                "1978.",
         "doi": "10.2118/7475-MS"},
    12:
        {"autor": "Standing, M.B.",
         "title": "Volumetric and Phase Behavior of Oil Field Hydrocarbon"
                  "Systems.",
         "ref": "Society of Petroleum Engineers, Dallas, TX. 1977",
         "doi": ""},
    13:
        {"autor": "Willman, B., Teja, A.",
         "title": "Prediction of dew points of semicontinuous natural gas and"
                  "petroleum mixtures. 1. Characterization by use of an"
                  "effective carbon number and ideal solution predictions",
         "ref": "Ind. Eng. Chem. Res., 1987, 26 (5), pp 948–952",
         "doi": "10.1021/ie00065a017"},
    14:
        {"autor": "Magoulas, S., Tassios, D.",
         "title": "Predictions of phase behavior of HT-HP reservoir fluids.",
         "ref": "Paper SPE 37294, Society of Petroleum Engineers, Richardson,"
                "TX, 1990.",
         "doi": ""},
    15:
        {"autor": "Sancet, J.,",
         "title": "Heavy Faction C7+ Characterization for PR-EOS.",
         "ref": "SPE 113025, 2007 SPE Annual Conference, November 11–14,"
                "Anaheim, CA 2007.",
         "doi": "10.2118/113026-stu"},
    16:
        {"autor": "Silva, M.B., Rodriguez, F.",
         "title": "Automatic fitting of equations of state for phase behavior"
                  "matching.",
         "ref": "Paper SPE 23703, Society of Petroleum Engineers, Richardson,"
                "TX, 1992.",
         "doi": "10.2118/23703-MS"},
    17:
        {"autor": "Soreide, I.",
         "title": "Improved Phase Behavior Predictions of Petroleum Reservoir"
                  "Fluids From a Cubic Equation of State.",
         "ref": "Doctor of engineering dissertation. Norwegian Institute of"
                "Technology, Trondheim, 1989.",
         "doi": ""},
    18:
        {"autor": "Hall, K. R., and L. Yarborough",
         "title": "New Simple Correlation for Predicting Critical Volume",
         "ref": "Chemical Engineering (November 1971): 76",
         "doi": ""},
    19:
        {"autor": "Edmister, W.C.",
         "title": "Applied Hydrocarbon Thermodynamics, Part 4, Compressibility"
                  "Factors and Equations of State",
         "ref": "Petroleum Refiner. 37 (April, 1958), 173–179",
         "doi": ""},
    20:
        {"autor": "API",
         "title": "Technical Data book: Petroleum Refining 6th Edition",
         "ref": "",
         "doi": ""},
    21:
        {"autor": "Salerno, S, Cascella, M., May, D., Watson, P., Tassios, D.",
         "title": "Prediction of Vapor Pressures and Saturated Volumes with a"
                  "Simple Cubic Equation of State: Part I. A Reliable Data"
                  "Base",
         "ref": "Fluid Phase Equilibria, Volume 27, 1986, Pages 15-34",
         "doi": "10.1016/0378-3812(86)87038-8"},
    22:
        {"autor": "Hougen, O. A., K. M. Watson, and R. A. Ragatz.",
         "title": "Chemical Process Principles, 2nd ed.",
         "ref": "New York: Wiley, 1959, p. 577.",
         "doi": ""},
    23:
        {"autor": "Reid, R., J. M. Prausnitz, and T. Sherwood.",
         "title": "The Properties of Gases and Liquids, 3rd ed. New York:"
                  "McGraw-Hill, 1977, p. 21.",
         "ref": "",
         "doi": ""},
    24:
        {"autor": "Nath, J.",
         "title": "Acentric Factor and the Critical Volumes for Normal Fluids",
         "ref": "Industrial Engineering and Chemical. Fundamentals 21, no. 3"
                "(1985): 325–326.",
         "doi": "10.1021/i100007a023"},
    25:
        {"autor": "Lee, B. I. and Kesler, M. G.",
         "title": "A Generalized Thermodynamic Correlation Based on"
                  "Three-Parameter Corresponding States",
         "ref": "American Institute of Chemical Engineers Journal, 21, 1975",
         "doi": "10.1002/aic.690210313"},
    26:
        {"autor": "Riazi, M. R., Al-Sahhaf, T. A. and Sl-Shammari M. A.",
         "title": "A Generalized Method for Estimation of Critical Constants",
         "ref": "Fluid Phase Equilibria, Vol. 147, 1998, pp. 1-6",
         "doi": "10.1016/s0378-3812(98)00251-9"},
    27:
        {"autor": "Riazi, M. R. and A1-Sahhaf, T. A.",
         "title": "Physical Properties of Heavy Petroleum Fractions and Crude"
                  "Oils",
         "ref": "Fluid Phase Equilibria, Vol. 117, 1996, pp. 217-224.",
         "doi": "10.1016/s0378-3812(98)00251-9"},
    28:
        {"autor": "Riazi, M. R. and A1-Sahhaf, T.",
         "title": "Physical Properties of n-Alkanes and n-Alkyl Hydrocarbons: "
                  "Application to Petroleum Mixtures",
         "ref": "Industrial and Engineering Chemistry Research, Vol. 34, 1995,"
                " pp. 4145-4148.",
         "doi": "10.1021/ie00038a062"},
    29:
        {"autor": "Riazi, M. R.",
         "title": "Characterization and Properties of Petroleum Fractions.",
         "ref": "ASTM manual series MNL50, 2005",
         "doi": ""},
    30:
        {"autor": "Papay, J.A.,",
         "title": "Termelestechnologiai Parameterek Valtozasa a Gazlelepk"
                  "Muvelese",
         "ref": "Soran. OGIL MUSZ, Tud, Kuzl. [Budapest], 1985. pp. 267–273.",
         "doi": ""},
    31:
        {"autor": "Hall, K. R., and L. Yarborough",
         "title": "A New Equation of State for Z-factor Calculations",
         "ref": "Oil and Gas Journal (June 18, 1973): 82–92.",
         "doi": ""},
    32:
        {"autor": "Dranchuk, P. M., and J. H. Abu-Kassem",
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
        {"autor": "Brill, J .P. and Beggs, H .D.",
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
        {"autor": "Dranckuk, P.M., Quon, D.",
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
         "doi": "10.1016/j.fluid.2003.02.003"},
    51:
        {"autor": "Riazi, M. R. and Daubert, T. E.",
         "title": "Analytical Correlations Interconvert Distillation Curve "
                  "Types",
         "ref": "Oil & Gas Journal, Vol. 84, 1986, August 25, pp. 50-57",
         "doi": ""},
    52:
        {"autor": "Daubert, T. E.",
         "title": "Petroleum Fraction Distillation Interconversion",
         "ref": "Hydrocarbon Processing, Vol. 73, No. 9, 1994, pp. 75-78.",
         "doi": ""},
    53:
        {"autor": "Edmister, W. C. and Okamoto, K. K.",
         "title": "Applied Hydrocarbon Thermodynamics, Part 13: Equilibrium "
                  "Flash Vaporization Correlations for Heavy Oils Under "
                  "Subatmospheric Pressures",
         "ref": "Petroleum Refiner, Vol. 38, No. 9, 1959, pp. 271-288.",
         "doi": ""},
    54:
        {"autor": "Riazi, M. R.",
         "title": "Distribution Model for Properties of Hydrocarbon-Plus "
                  "Fractions",
         "ref": "Industrial and Engineering Chemistry Research, Vol. 28, 1989,"
                " pp. 1731-1735.",
         "doi": "10.1021/ie00095a026"},





    55:
        {"autor": "",
         "title": "",
         "ref": "",
         "doi": ""},
    56:
        {"autor": "",
         "title": "",
         "ref": "",
         "doi": ""},
    57:
        {"autor": "",
         "title": "",
         "ref": "",
         "doi": ""},
    58:
        {"autor": "",
         "title": "",
         "ref": "",
         "doi": ""},
    59:
        {"autor": "",
         "title": "",
         "ref": "",
         "doi": ""},
        }


def _unit(key, val, unit=None):
    """Set appropiate unit for value

    Parameters
    ------------
    key : string
        Property name
    val : float
        Property value to unitize
    unit : string, optional
        In case it necessary use a alternate unit

    Returns
    -------
    x : unidades.unidad
        A instance of appropiate unidades.unidad subclass

    Notes
    -----
    Default use the petroleum english standards units, R, psi, ft³/lb...
    """
    if key in ("Tc", "Tb"):
        if unit is None:
            unit = "R"
        x = unidades.Temperature(val, unit)
    elif key == "Pc":
        if unit is None:
            unit = "psi"
        x = unidades.Pressure(val, unit)
    elif key == "Vc":
        if unit is None:
            unit = "ft3lb"
        x = unidades.SpecificVolume(val, unit)
    else:
        x = unidades.Dimensionless(val)
    return x


# Generalized Correlations
def prop_Ahmed(Nc):
    """Calculate petroleum fractions properties with the carbon atom number as
    the only known property

    Parameters
    ------------
    Nc : float
        Carbon atom number [-]

    Returns
    -------
    prop : A dict with the calculated properties
        M: Molecular weight, [-]
        Tc: Critic temperature, [ºR]
        Pc: Critic pressure, [psi]
        Tb: Normal boiling temperature, [ºR]
        w: Acentric factor, [-]
        SG: Specific gravity, [-]
        Vc: Critic volume, [ft³/lb]

    Examples
    --------
    >>> "%0.0f" % prop_Ahmed(6)["Tb"].R
    '604'
    >>> "%0.5f" % prop_Ahmed(45)["Vc"].ft3lb
    '0.06549'

    References
    ----------
    [1] .. Katz, D. L., and Firoozabadi, A. Predicting Phase Behavior of
        Condensate/Crude-oil Systems Using Methane Interaction Coefficients
        Journal of Petroleum Technology, Nov. 1978, pp. 1649–1655
    [2] .. Ahmed, T., G. Cady, and A. Story. A Generalized Correlation for
        Characterizing the Hydrocarbon Heavy Fractions. Paper SPE 14266,
        presented at the 60th Annual Technical Conference of the Society of
        Petroleum Engineers, Las Vegas, September 22–25, 1985.
    """
    a1 = [-131.11375, 915.53747, 275.56275, 434.38878, -0.50862704,
          0.86714949, 5.223458e-2]
    a2 = [24.96156, 41.421337, -12.522269, 50.125279, 8.700211e-2,
          3.41434080e-3, 7.87091369e-4]
    a3 = [-0.34079022, -0.7586849, 0.29926384, -0.9097293, -1.8484814e-3,
          -2.839627e-5, -1.9324432e-5]
    a4 = [2.4941184e-3, 5.8675351e-3, -2.8452129e-3, 7.0280657e-3,
          1.466389e-5, 2.4943308e-8, 1.7547264e-7]
    a5 = [468.32575, -1.3028779e3, 1.7117226e3, -601.856510, 1.8518106,
          -1.1627984, 4.4017952e-2]

    prop = {}
    for i, key in enumerate(["M", "Tc", "Pc", "Tb", "w", "SG", "Vc"]):
        x = a1[i] + a2[i]*Nc + a3[i]*Nc**2 + a4[i]*Nc**3 + a5[i]/Nc
        val = _unit(key, x)
        prop[key] = val
    return prop


def prop_Riazi_Daubert_1980(Tb, SG):
    """Calculate petroleum fractions properties with the Riazi (1980)
    correlation with the boiling temperature and specific gravity as input
    paramters

    Parameters
    ------------
    Tb : float
        Normal boiling temperature, [K]
    SG: float
        Specific gravity, [-]

    Returns
    -------
    prop : A dict with the calculated properties
        M: Molecular weight, [-]
        Tc: Critic temperature, [ºR]
        Pc: Critic pressure, [psi]
        Vc: Critic volume, [ft3/lb]

    Examples
    --------
    Example 2.2 from [9]_: C7+ fraction with Tb=198ºF and SG=0.7365

    >>> T = unidades.Temperature(198, "F")
    >>> p = prop_Riazi_Daubert_1980(T, 0.7365)
    >>> "%.0f %.0f %.0f %.4f" % (p["M"], p["Tc"].R, p["Pc"].psi, p["Vc"].ft3lb)
    '96 990 467 0.0623'

    References
    ----------
    [4] .. Riazi, M. R., and T. E. Daubert. Simplify Property Predictions.
        Hydrocarbon Processing (March 1980): 115–116.
    [9] .. Tarek Ahmed. Equations of State and PVT Analysis: Applications for
        Improved Reservoir Modeling, 2nd Edition. Gulf Professional Publishing,
        2016, ISBN 9780128015704
    """
    # Convert input Tb in Kelvin to Rankine to use in the correlation
    Tb_R = unidades.K2R(Tb)

    a = [4.5673e-5, 24.2787, 3.12281e9, 7.5214e-3]
    b = [2.1962, 0.58848, -2.3125, 0.2896]
    c = [-1.0164, 0.3596, 2.3201, -0.7666]

    prop = {}
    prop["Tb"] = unidades.Temperature(Tb)
    prop["SG"] = unidades.Dimensionless(SG)

    for i, key in enumerate(["M", "Tc", "Pc", "Vc"]):
        x = a[i]*Tb_R**b[i]*SG**c[i]
        val = _unit(key, x)
        prop[key] = val

    return prop


def prop_Riazi_Daubert(tita1, val1, tita2, val2):
    """Calculate petroleum fractions properties known two properties

    Parameters
    ------------
    tita1 : string
        Name of first known property
    tita2 : string
        Name of second known property
    val1: float
        First known property value
    val2: float
        Second known property value

    Returns
    -------
    prop : A dict with the calculated properties
        Tc: Temperatura crítica, [ºR]
        Pc: Presión crítica, [psi]
        Vc: Volumen crítico, [ft³/lb]
        M: Peso molecular, [-]
        Tb: Temperatura fusión, [ºR]
        SG: Gravedad específica, [-]
        I: Huang characterization factor, [-]
        CH: Carbon/hydrogen weight ratio, [-]

    Notes
    -----
    The available input properties for tita are:
        Tb: Temperatura fusión, [ºR]
        SG: Gravedad específica, [-]
        M: Peso molecular, [-]
        CH: Carbon/hydrogen weight ratio, [-]
        I: Huang characterization factor, [-]
        v1: Kinematic viscosity at 100ºF, [m²/s]

    Raises
    ------
    NotImplementedError : If the input pair are unsupported

    NotImplementedError : If input isn't in limit
        * 80ºF ≤ T ≤ 650ºF
        * 70 ≤ M ≤ 300

    Examples
    --------
    Example 2.1 from [9]_: C7+ fraction with M=150 and SG=0.78

    >>> p = prop_Riazi_Daubert("M", 150, "SG", 0.78)
    >>> "%i %i %.4f %.0f" % (p["Tc"].R, p["Pc"].psi, p["Vc"].ft3lb, p["Tb"].R)
    '1160 320 0.0636 825'

    Example 2.2 from [9]_: Petroleum fraction with Tb=198ºF and SG=0.7365

    >>> T = unidades.Temperature(198, "F")
    >>> p = prop_Riazi_Daubert("Tb", T, "SG", 0.7365)
    >>> "%.0f %.0f %.0f %.4f" % (p["M"], p["Tc"].R, p["Pc"].psi, p["Vc"].ft3lb)
    '97 986 466 0.0626'
    >>> "%.3f" % prop_Edmister(Tc=p["Tc"], Pc=p["Pc"], Tb=p["Tb"])["w"]
    '0.287'

    References
    ----------
    [3] .. Riazi, M. R., and T. E. Daubert. “Characterization Parameters for
        Petroleum Fractions.” Industrial Engineering and Chemical Research 26,
        no. 24 (1987): 755–759.
    [9] .. Tarek Ahmed. Equations of State and PVT Analysis: Applications for
        Improved Reservoir Modeling, 2nd Edition. Gulf Professional Publishing,
        2016, ISBN 9780128015704
    """
    # Check correct input parameter
    p1 = ("Tb", "M", "v1")
    p2 = ("SG", "I", "CH")

    if tita1 in p2 and tita2 in p1:
        # Invert input parameters
        tita1, tita2 = tita2, tita1
        val1, val2 = val2, val1

    if tita1 not in p1 or tita2 not in p2:
        raise NotImplementedError(QApplication.translate(
            "pychemqt", "Undefined input pair"))
    elif tita1 == "M" and (val1 < 70 or val1 > 300):
        raise NotImplementedError(QApplication.translate(
            "pychemqt", "Molecular weight input out of bounds"))
    elif tita1 == "Tb" and (val1 < 300 or val1 > 620):
        raise NotImplementedError(QApplication.translate(
            "pychemqt", "Boiling temperature input out of bounds"))

    # Convert input Tb in Kelvin to Rankine to use in the correlation
    if tita1 == "Tb":
        val1 = unidades.K2R(val1)

    par = {
      # Table IV
      "Tc": {
        "Tb-SG": [10.6443, -5.1747e-4, -.54444, 3.5995e-4, .81067, .53691],
        "Tb-I": [5.62596e5, -7.317e-4, -16.9097, 2.5131e-3, .6154, 4.3469],
        "Tb-CH": [2.2452, 1.9152e-4, -0.06487, -6.0192e-4, 0.7699, 0.900],
        "M-SG": [554.4, -1.3478e-4, -0.61641, 0.0, 0.2998, 1.0555],
        "M-I": [2.4254e6, 2.001e-4, -13.049, 0.0, 0.2383, 4.0642],
        "M-CH": [37.332, 1.3848e-3, -0.1379, -2.7e-4, 0.3526, 1.4191],
        "v1-SG": [251.026, -3.177e-2, 1.6587, 0.0, 0.1958, -0.9431],
        "v1-I": [4.414e3, -0.0291, -1.2664, 0.0, 0.1884, 0.7492],
        "v1-CH": [4.939e2, -2.8e-2, -8.91e-2, 0.0, 0.1928, 0.7744]},

      # Table V
      "Pc": {
        "Tb-SG": [6.162e6, -4.725e-3, -4.8014, 3.1939e-3, -0.4844, 4.0846],
        "Tb-I": [2.2337e25, -6.7041e-3, -74.5612, .019, -1.0303, 18.43302],
        "Tb-CH": [158.96, -2.1357e-3, -0.3454, 0.0, -0.1801, 3.2223],
        "M-SG": [4.5203e4, -1.8078e-3, -0.3084, 0.0, -0.8063, 1.6015],
        "M-I": [2.9384e17, -0.01415, -48.5809, 0.0451, -0.8097, 12.9148],
        "M-CH": [815.99, -2.139e-3, -0.265, 0.0, -0.6616, 2.4004],
        "v1-SG": [1.271e5, -0.2523, -5.6028, 0.355, -0.5913, 6.0793],
        "v1-I": [6.1475e22, -0.4586, -71.905, 1.8854, -0.6395, 20.7032],
        "v1-CH": [40.9115, 0.01906, 0.1323, 0.0, 0.471, 1.6306]},

      # Table VI
      "Vc": {
        "Tb-SG": [6.233e-4, -1.4679e-3, -.26404, 1.095e-3, .7506, -1.2028],
        "Tb-I": [1.3077e-3, -1.799e-3, -3.5349, 4.425e-3, .7687, -.72011],
        "Tb-CH": [0.2048, -9.2189e-4, 0.05345, 1.4805e-4, 0.1657, -1.4439],
        "M-SG": [1.206e-2, -2.657e-3, 0.5287, 2.6012e-3, 0.20378, -1.3036],
        "M-I": [1.016e-6, -2.0208e-3, 14.1853, 4.5318e-3, .2556, -4.60413],
        "M-CH": [0.2558, -2.3533e-3, 0.1082, 3.826e-4, 0.0706, -1.3362],
        "v1-SG": [1.64424e-3, -2.04563e-1, 3.513392, 2.12365e-1, 1.19093e-1,
                  -3.801261],
        "v1-I": [3.219523e-12, -1.63181e-1, 36.09011, .4608, .1417, -10.65067],
        "v1-CH": [.245582, -.11261, .086387, .016031, .046004, -1.028488]},

      # Table VII
      "M": {
        "Tb-SG": [581.96, 5.43076e-4, -9.53384, 1.11056e-3, 0.97476,
                  6.51274],
        "Tb-I": [2.606e-6, 8.6574e-6, 4.2376, 0.0, 2.0935, -1.9985],
        "Tb-CH": [3.06584e-3, 5.3305e-4, 7.9113e-2, -2.87657e1, 1.6736,
                  -0.68681],
        "v1-SG": [1.51723e6, -0.195411, -9.63897, .16247, .56370, 6.89383],
        "v1-I": [4.0e-9, -8.9854e-2, 38.106, 0.0, 0.6675, -10.6],
        "v1-CH": [84.1505, -5.976e-2, -0.10741, 0.0, 0.5596, 0.65815],
        "v1-v2": [288.916, 0.1380, -0.7311, -5.704e-3, 0.051, 0.8411]},

      # Table VIII
      "Tb": {
        "M-SG": [6.77857, 3.77409e-3, 2.984036, -4.25288e-3,
                 4.01673e-1, -1.58262],
        "M-I": [136.395, 0.0, 0.0, 0.0, 0.4748, 0.4283],
        "M-CH": [36.45625, -1.57415e-4, -4.5707e-2, 9.22926e-6, 5.12976e-1,
                 4.72372e-1],
        "v1-SG": [4.28375e3, -1.3051e-2, -1.68759, -2.1247e-2, 2.62914e-1,
                  1.34890],
        "v1-I": [9.1133e-2, -6.5236e-2, 14.9371, 6.029e-2, .3228, -3.8798],
        "v1-CH": [444.377, -3.8093e-2, -7.7305e-2, 0.0, 0.2899, 0.6056]},

      # Table IX
      "SG": {
        "Tb-I": [6.9195, -2.33e-4, -23.5535, 2.2152e-3, 2.9806e7, -0.3418],
        "Tb-CH": [7.3238e-1, -1.01845e-3, -8.1635e-2, 3.60649e-5,
                  1.69916e-3, 8.90041e-1],
        "M-I": [6.3028, -1.588e-3, -20.594, 7.344e-3, 1.1284e6, -7.71e-2],
        "M-CH": [9.19255e-1, -1.48844e-3, -7.925e-2, 4.921118e-5,
                 6.84403e-2, 2.89844e-1],
        "v1-I": [8.04224, -6.1406e-2, -26.3934, .2533, 3.8083e7, -.02353],
        "v1-CH": [1.17777, 0.02614, -0.10966, -5.654e-3, .18242, .05245]},

      # Table X
      "I": {
        "Tb-SG": [0.022657, 3.9052e-4, 2.468316, -5.70425e-4, 5.7209e-2,
                  -0.719895],
        "Tb-CH": [4.307e-3, -9.8747e-5, -6.0737e-2, -4.414e-5, 0.4470, 0.9896],
        "M-SG": [0.422375, 3.18857e-4, -0.200996, -4.24514e-4,
                 -8.43271e-3, 1.117818],
        "M-CH": [4.239e-2, -5.6946e-4, -6.836e-2, 0.0, 0.1656, 0.8291],
        "v1-SG": [.26376, 1.7458e-2, .231043, -1.8441e-2, -1.1275e-2, .770779],
        "v1-CH": [0.08716, 6.1396e-3, -7.019e-2, -2.5935e-3, .05166, .84599]},

      # Table XI
      "CH": {
        "Tb-SG": [17.22022, 8.24983e-3, 16.9402, -6.93931e-3, -2.72522,
                  -6.79769],
        "Tb-I": [1.8866e-12, 4.2873e-3, 71.6531, -0.0116, -1.3773, -13.6139],
        "M-SG": [2.35475, 9.3485e-3, 4.74695, -8.01719e-3, -0.68418, -0.7682],
        "M-I": [2.9004e-13, 7.8276e-3, 60.3484, -0.02445, -0.37884, -12.34051],
        "v1-SG": [2.523e-12, .482811, 29.98797, -0.55768, -.146565, -20.31303],
        "v1-I": [2.143e-12, 0.2832, 53.7316, -0.91085, -0.17158, -10.88065]}}

    prop = {}
    prop[tita1] = _unit(tita1, val1)
    prop[tita2] = _unit(tita2, val2)

    for key in par.keys():
        if key not in (tita1, tita2):
            a, b, c, d, e, f = par[key]["%s-%s" % (tita1, tita2)]
            x = a*val1**e*val2**f*exp(b*val1+c*val2+d*val1*val2)         # Eq 3

            val = _unit(key, x)
            prop[key] = val

    return prop


def prop_Cavett(Tb, API):
    """Calculate petroleum fractions properties with the Cavett (1980)
    correlation. Use API as a alternate input parameter to specific gravity

    Parameters
    ------------
    Tb : float
        Normal boiling temperature, [K]
    API : float
        API gravity, [-]

    Returns
    -------
    prop : A dict with the calculated properties
        Tc: Critic temperature, [ªR]
        Pc: Critic pressure, [Pa]

    Examples
    --------
    Example 2.2 from [9]_: Petroleum fraction with Tb=198ºF and SG=0.7365

    >>> T = unidades.Temperature(198, "F")
    >>> API = 141.5/0.7365-131.5
    >>> p = prop_Cavett(T, API)
    >>> "%.1f %.0f" % (p["Tc"].R, p["Pc"].psi)
    '978.1 466'

    References
    ----------
    [7] .. Cavett, R.H., 1962. Physical data for distillation calculations,
        vapor-liquid equilibrium. Proceedings of the 27th Meeting, API, San
        Francisco, 1962, Issue 3, pp. 351–366.
    [9] .. Tarek Ahmed. Equations of State and PVT Analysis: Applications for
        Improved Reservoir Modeling, 2nd Edition. Gulf Professional Publishing,
        2016, ISBN 9780128015704
    """
    # Convert input Tb in Kelvin to Fahrenheit to use in the correlation
    Tb_F = unidades.K2F(Tb)

    a = [768.07121, 1.7133693, -0.0010834003, -0.0089212579, 0.38890584e-6,
         0.5309492000e-5, 0.3271160000e-7]
    Tc = a[0] + a[1]*Tb_F + a[2]*Tb_F**2 + a[3]*API*Tb_F + a[4]*Tb_F**3 + \
        a[5]*API*Tb_F**2 + a[6]*API**2*Tb_F**2

    b = [2.82904060, 0.94120109e-3, -0.30474749e-5, -0.20876110e-4,
         0.15184103e-8, 0.11047899e-7, -0.48271599e-7, 0.13949619e-9]
    Pc = 10**(b[0] + b[1]*Tb_F + b[2]*Tb_F**2 + b[3]*API*Tb_F + b[4]*Tb_F**3 +
              b[5]*API*Tb_F**2+b[6]*API**2*Tb_F+b[7]*API**2*Tb_F**2)

    prop = {}
    prop["Tb"] = unidades.Temperature(Tb)
    prop["API"] = unidades.Dimensionless(API)
    prop["Tc"] = unidades.Temperature(Tc, "R")
    prop["Pc"] = unidades.Pressure(Pc, "psi")
    return prop


def prop_Lee_Kesler(Tb, SG):
    """Calculate petroleum fractions properties with the Lee-Kesker (1976)
    correlation. Use API as a alternate input parameter to specific gravity

    Parameters
    ------------
    Tb : float
        Normal boiling temperature, [K]
    SG: float
        Specific gravity, [-]

    Returns
    -------
    prop : A dict with the calculated properties
        Tc: Critic temperature, [ªR]
        Pc: Critic pressure, [Pa]
        M: Molecular weight, [-]
        w: Acentric factor, [-]

    Examples
    --------
    Example 2.2 from [9]_: Petroleum fraction with Tb=198ºF and SG=0.7365

    >>> T = unidades.Temperature(198, "F")
    >>> p = prop_Lee_Kesler(T, 0.7365)
    >>> "%.0f %.0f %.1f %.3f" % (p["Tc"].R, p["Pc"].psi, p["M"], p["w"])
    '981 470 98.6 0.306'

    References
    ----------
    [8] .. Kesler, M. G., and B. I. Lee. Improve Prediction of Enthalpy of
        Fractions.” Hydrocarbon Processing Volume 55, Issue 3, March 1976,
        Pages 153-158
    [9] .. Tarek Ahmed. Equations of State and PVT Analysis: Applications for
        Improved Reservoir Modeling, 2nd Edition. Gulf Professional Publishing,
        2016, ISBN 9780128015704
    """
    # Convert input Tb in Kelvin to Rankine to use in the correlation
    Tb_R = unidades.K2R(Tb)

    Tc = 341.7+811.1*SG+(0.4244+0.1174*SG)*Tb_R+(0.4669-3.26238*SG)*1e5/Tb_R
    Pc = exp(8.3634 - 0.0566/SG - (0.24244+2.2898/SG+0.11857/SG**2)*1e-3*Tb_R +
             (1.4685+3.648/SG+0.47227/SG**2)*1e-7*Tb_R**2 -
             (0.42019+1.6977/SG**2)*1e-10*Tb_R**3)
    M = -12272.6 + 9486.4*SG + (4.6523-3.3287*SG)*Tb_R + \
        (1-0.77084*SG-0.02058*SG**2)*(1.3437-720.79/Tb_R)*1e7/Tb_R + \
        (1-0.80882*SG+0.02226*SG**2)*(1.8828-181.98/Tb_R)*1e12/Tb_R**3

    tr = Tb_R/Tc
    Kw = Tb_R**(1./3)/SG
    if tr > 0.8:
        w = -7.904+0.1352*Kw-0.007465*Kw**2+8359*tr+(1.408-0.01063*Kw)/tr
    else:
        w = (-log(Pc/14.7)-5.92714+6.09648/tr+1.28862*log(tr)-0.169347*tr**6)/(
            15.2518-15.6875/tr-13.4721*log(tr)+0.43577*tr**6)

    prop = {}
    prop["Tb"] = unidades.Temperature(Tb)
    prop["SG"] = unidades.Dimensionless(SG)
    prop["Tc"] = unidades.Temperature(Tc, "R")
    prop["Pc"] = unidades.Pressure(Pc, "psi")
    prop["M"] = unidades.Dimensionless(M)
    prop["w"] = unidades.Dimensionless(w)
    return prop


def prop_Sim_Daubert(Tb, SG):
    """Calculate petroleum fractions properties with the Sim-Daubert (1980)
    computerized version of Winn (1957) graphical correlations

    Parameters
    ------------
    Tb : float
        Normal boiling temperature, [K]
    SG: float
        Specific gravity, [-]

    Returns
    -------
    prop : A dict with the calculated properties
        M: Molecular weight, [-]
        Tc: Critic temperature, [ªR]
        Pc: Critic pressure, [Pa]

    Examples
    --------
    Example 2.2 from [9]_: Petroleum fraction with Tb=198ºF and SG=0.7365

    >>> T = unidades.Temperature(198, "F")
    >>> p = prop_Sim_Daubert(T, 0.7365)
    >>> "%.0f %.0f %.0f" % (p["Tc"].R, p["Pc"].psi, p["M"])
    '979 479 96'

    References
    ----------
    [6] .. Sim, W.J. and Daubert, T.E., 1980, Prediction of vapor-liquid
        equilibria of undefined mixtures, Industrial and Engineering
        Chemistry–Process Design and Development, 19: 386 –393.
    [9] .. Tarek Ahmed. Equations of State and PVT Analysis: Applications for
        Improved Reservoir Modeling, 2nd Edition. Gulf Professional Publishing,
        2016, ISBN 9780128015704
    """

    M = 5.805e-5*Tb**2.3776/SG**0.9371                                   # Eq 3
    Pc = 6.1483e12*Tb**-2.3177*SG**2.4853                                # Eq 4
    Tc = exp(4.2009*Tb**0.08615*SG**0.04614)                             # Eq 5

    prop = {}
    prop["Tb"] = unidades.Temperature(Tb)
    prop["SG"] = unidades.Dimensionless(SG)
    prop["M"] = unidades.Dimensionless(M)
    prop["Tc"] = unidades.Temperature(Tc, "R")
    prop["Pc"] = unidades.Pressure(Pc)
    return prop


def prop_Watansiri_Owens_Starling(Tb, SG, M):
    """Calculate petroleum fractions properties with the Watansiri-Owens-
    Starling (1985) correlation using boiling temperature and molecular weight
    as input parameters

    Parameters
    ------------
    Tb : float
        Normal boiling temperature, [K]
    SG: float
        Specific gravity, [-]
    M: float
        Molecular weight, [-]

    Returns
    -------
    prop : A dict with the calculated properties
        Tc: Critic temperature, [K]
        Pc: Critic pressure, [atm]
        Vc: Critic volume, [cm3/gr]
        w: Acentric factor, [-]
        Dm: Dipole moment, [Debye]

    Examples
    --------
    Example 2.2 from [9]_: Petroleum fraction with Tb=198ºF and SG=0.7365

    >>> T = unidades.Temperature(198, "F")
    >>> p = prop_Watansiri_Owens_Starling(T, 0.7365, 96)
    >>> "%.0f %.5f" % (p["Tc"].R, p["Vc"].ft3lb)
    '980 0.06548'

    References
    ----------
    [10] .. Watansiri, S., Owens, V.H., Starling, K.E., 1985. Correlations for
        estimating critical constants, acentric factor, and dipole moment for
        undefined coal-fluid fractions. Ind. Eng. Chem. Process. Des. Dev. 24,
        294–296.
    [9] .. Tarek Ahmed. Equations of State and PVT Analysis: Applications for
        Improved Reservoir Modeling, 2nd Edition. Gulf Professional Publishing,
        2016, ISBN 9780128015704
    """

    Tc = exp(-0.00093906*Tb + 0.03095*log(M) + 1.11067*log(Tb) +
             M*(0.078154*SG**0.5 - 0.061061*SG**(1./3.) - 0.016943*SG))  # Eq 1
    Vc = exp(80.4479 - 129.8038*SG + 63.175*SG**2 - 13.175*SG**3 +
             1.10108*log(M) + 42.1958*log(SG))                           # Eq 2
    Pc = exp(3.95431 + 0.70682*(Tc/Vc)**0.8 - 4.84*M/Tc - 0.15919*Tb/M)  # Eq 3
    w = (0.92217e-3*Tb + 0.507288*Tb/M + 382.904/M + 0.242e-5*(Tb/SG)**2 -
         0.2165e-4*Tb*M + 0.1261e-2*SG*M + 0.1265e-4*M**2 + 0.2016e-4*SG*M**2 -
         80.6495*Tb**(1/3.)/M - 0.378e-2*Tb**(2/3.)/SG**2)*Tb/M          # Eq 4

    # Eq 6
    HVNP = -10397.5 + 46.2681*Tb - 1373.91*Tb**0.5 + 4595.81*log(Tb)
    u1 = (197.933/M + 0.039177*M)*Vc/Tc
    u2 = 0.3185e-2*Vc + 0.956247e-2*Tb - 0.5479e-3*HVNP
    u3 = -1.34634*w + 0.906609*log(w)
    u4 = -4.85638 - 0.013548*M + 0.271949e-3*M**2 + 1.04024*log(M)
    Dm = u1/u4 + u2*u3

    prop = {}
    prop["Tb"] = unidades.Temperature(Tb)
    prop["M"] = unidades.Dimensionless(M)
    prop["Tc"] = unidades.Temperature(Tc)
    prop["Pc"] = unidades.Pressure(Pc, "atm")
    prop["Vc"] = unidades.SpecificVolume(Vc/M, "ccg")
    prop["w"] = unidades.Dimensionless(w)
    prop["Dm"] = unidades.DipoleMoment(Dm, "Debye")
    return prop


def prop_Rowe(M):
    """Calculate petroleum fractions properties with the Rowe (1986)
    correlations (1978) using only the molecular weight as input parameter

    Parameters
    ------------
    M: float
        Molecular weight, [-]

    Returns
    -------
    prop : A dict with the calculated properties
        Tc: Critic temperature, [R]
        Tb: Boiling temperature, [R]
        Pc: Critic pressure, [psi]

    Examples
    --------
    Example 2.3 from [9]_: Petroleum fraction with M=216 and SG=0.8605

    >>> p = prop_Rowe(216)
    >>> "%.1f" % p["Tc"].R
    '1279.8'

    References
    ----------
    [11] .. Rowe, A.M. Internally consistent correlations for predicting phase
        compositions for use in reservoir compositional simulators. Paper SPE
        7475, In: Presented at the 53rd Annual Society of Petroleum Engineers
        Fall Technical Conference and Exhibition, 1978.
    [9] .. Tarek Ahmed. Equations of State and PVT Analysis: Applications for
        Improved Reservoir Modeling, 2nd Edition. Gulf Professional Publishing,
        2016, ISBN 9780128015704

    Notes
    -----
    Critical pressure examples in [9]_ don't get result, the equations isn't
    equal in both references
    """
    n = (M-2)/14                                                        # Eq D1
    Tc = (961-10**(2.95597-(0.090597*n**(2/3.))))*1.8                   # Eq D2
    Tb = 4.347e-4*Tc**2+265                                             # Eq D3
    Y = -0.0134426826*n+0.6801481651                                    # Eq D4
    C = 10**Y*1e5                                                       # Eq D5
    Pc = C/Tc                                                           # Eq D6

    prop = {}
    prop["M"] = unidades.Dimensionless(M)
    prop["Tc"] = unidades.Temperature(Tc, "R")
    prop["Tb"] = unidades.Temperature(Tb, "R")
    prop["Pc"] = unidades.Pressure(Pc, "psi")
    return prop


def prop_Standing(M, SG):
    """Calculate petroleum fractions properties with the Standing (1977)
    correlations based in the graphical plot of Mathews et al. (1942) using
    molecular weight and specific gravity as input parameter

    Parameters
    ------------
    M: float
        Molecular weight, [-]
    SG: float
        Specific gravity, [-]

    Returns
    -------
    prop : A dict with the calculated properties
        Tc: Critic temperature, [K]
        Pc: Critic pressure, [atm]

    Examples
    --------
    Example 2.3 from [9]_: Petroleum fraction with M=216 and SG=0.8605

    >>> p = prop_Standing(216, 0.8605)
    >>> "%.1f %.0f" % (p["Tc"].R, p["Pc"].psi)
    '1269.3 270'

    References
    ----------
    [12] .. Standing, M.B., Volumetric and Phase Behavior of Oil Field
        Hydrocarbon Systems. Society of Petroleum Engineers, Dallas, TX. 1977
    [9] .. Tarek Ahmed. Equations of State and PVT Analysis: Applications for
        Improved Reservoir Modeling, 2nd Edition. Gulf Professional Publishing,
        2016, ISBN 9780128015704
    """
    Tc = 608 + 364*log10(M-71.2) + (2450*log10(M)-3800)*log10(SG)       # Eq 25
    Pc = 1188 - 431*log10(M-61.1) + (2319-852*log10(M-53.7))*(SG-0.8)   # Eq 26

    prop = {}
    prop["M"] = unidades.Dimensionless(M)
    prop["SG"] = unidades.Dimensionless(SG)
    prop["Tc"] = unidades.Temperature(Tc, "R")
    prop["Pc"] = unidades.Pressure(Pc, "psi")
    return prop


def prop_Willman_Teja(Tb):
    """Calculate petroleum fractions properties with the Willman-Teja (1987)
    correlations using only the boiling temperature as input parameter

    Parameters
    ------------
    Tb: float
        Boiling temperature, [K]

    Returns
    -------
    prop : A dict with the calculated properties
        Tc: Critic temperature, [K]
        Pc: Critic pressure, [MPa]
        n: Carbon number, [-]

    Examples
    --------
    Example 2.2 from [9]_: Petroleum fraction with Tb=198ºF and SG=0.7365

    >>> Tb = unidades.Temperature(198, "F")
    >>> p = prop_Willman_Teja(Tb)
    >>> "%.0f %.0f" % (p["Tc"].R, p["Pc"].psi)
    '977 442'

    References
    ----------
    [9] .. Tarek Ahmed. Equations of State and PVT Analysis: Applications for
        Improved Reservoir Modeling, 2nd Edition. Gulf Professional Publishing,
        2016, ISBN 9780128015704
    [13] .. Bert Willman, Amyn S. Teja, Prediction of dew points of
        semicontinuous natural gas and petroleum mixtures. 1. Characterization
        by use of an effective carbon number and ideal solution predictions.
        Ind. Eng. Chem. Res., 1987, 26 (5), pp 948–952
    """
    # Eq 11
    A = [95.50441892007, 3.742203001499, 2295.53031513, -1042.57256080,
         -22.66229823925, -1660.893846582, 439.13226915]

    def f(n):
        return A[0] + A[1]*n + A[2]*n**0.667 + A[3]*n**0.5 + A[4]*log10(n) + \
            A[5]*n**0.8 + A[6]*n**0.9 - Tb

    n = fsolve(f, 10)

    Tc = Tb*(1+1/(1.25127+0.137242*n))                                  # Eq 12
    Pc = (2.33761+8.16448*n)/(0.873159+0.54285*n)**2                    # Eq 13

    prop = {}
    prop["Tb"] = unidades.Dimensionless(Tb)
    prop["n"] = unidades.Dimensionless(n)
    prop["Tc"] = unidades.Temperature(Tc)
    prop["Pc"] = unidades.Pressure(Pc, "MPa")
    return prop


def prop_Magoulas_Tassios(M, SG):
    """Calculate petroleum fractions properties with the Magoulas-Tassios(1990)
    correlations using molecular weight and specific gravity as input parameter

    Parameters
    ------------
    M: float
        Molecular weight, [-]
    SG: float
        Specific gravity, [-]

    Returns
    -------
    prop : A dict with the calculated properties
        Tc: Critic temperature, [K]
        Pc: Critic pressure, [MPa]
        w: Acentric factor, [-]

    Examples
    --------
    Example 2.3 from [9]_: Petroleum fraction with M=216 and SG=0.8605

    >>> p = prop_Magoulas_Tassios(216, 0.8605)
    >>> "%.0f %.0f" % (p["Tc"].R, p["Pc"].psi)
    '1317 273'

    References
    ----------
    [9] .. Tarek Ahmed. Equations of State and PVT Analysis: Applications for
        Improved Reservoir Modeling, 2nd Edition. Gulf Professional Publishing,
        2016, ISBN 9780128015704
    [14] .. Magoulas, S., Tassios, D., Predictions of phase behavior of HT-HP
        reservoir fluids. Paper SPE 37294, Society of Petroleum Engineers,
        Richardson, TX, 1990.
    """
    Tc = -1247.4 + 0.792*M + 1971*SG - 27000./M + 707.4/SG
    Pc = exp(0.01901 - 0.0048442*M + 0.13239*SG + 227./M - 1.1663/SG +
             1.2702*log(M))
    w = -0.64235 + 0.00014667*M + 0.021876*SG - 4.539/M + 0.21699*log(M)

    prop = {}
    prop["M"] = unidades.Dimensionless(M)
    prop["SG"] = unidades.Dimensionless(SG)
    prop["Tc"] = unidades.Temperature(Tc, "R")
    prop["Pc"] = unidades.Pressure(Pc, "psi")
    prop["w"] = unidades.Dimensionless(w)
    return prop


def prop_Twu(Tb, SG):
    """Calculate petroleum fractions properties with the Two (1983)
    correlation with the boiling temperature and specific gravity as input
    paramters

    Parameters
    ------------
    Tb : float
        Normal boiling temperature, [K]
    SG: float
        Specific gravity, [-]

    Returns
    -------
    prop : A dict with the calculated properties
        M: Molecular weight, [-]
        Tc: Critic temperature, [ºR]
        Pc: Critic pressure, [psi]
        Vc: Critic volume, [ft3/lb]

    Examples
    --------
    >>> crit = prop_Twu(510, 1.097)
    >>> "%.1f %.1f %.1f" % (crit["Tc"].R, crit["Pc"].psi, crit["M"])
    '1380.3 556.8 130.4'

    References
    ----------
    [5] .. Twu, C.H., 1984, An internally consistent correlation for predicting
        the critical properties and molecular weights of petroleum and coal-tar
        liquids, Fluid Phase Equilbria, 16: 137–150.
    """

    # Convert input Tb in Kelvin to Rankine to use in the correlation
    Tb_R = unidades.K2R(Tb)

    Tco = Tb_R/(0.533272 + 0.191017e-3*Tb_R + 0.779681e-7*Tb_R**2 -
                0.284376e-10*Tb_R**3 + 0.959468e28/Tb_R**13)             # Eq 1
    alfa = 1-Tb_R/Tco                                                    # Eq 5
    Vco = (1-(0.419869 - 0.505839*alfa - 1.56436*alfa**3 -
              9481.7*alfa**14))**-8                                      # Eq 2
    SGo = 0.843593-0.128624*alfa-3.36159*alfa**3-13749.5*alfa**12        # Eq 3

    Pco = (3.83354 + 1.19629*alfa**0.5 + 34.8888*alfa + 36.1952*alfa**2 +
           104.193*alfa**4)**2                                           # Eq 8

    # Eq 4
    def f(M):
        Tita = log(M)
        return exp(5.71419 + 2.71579*Tita - 0.286590*Tita**2 - 39.8544/Tita -
                   0.122488/Tita**2) - 24.7522*Tita + 35.3155*Tita**2 - Tb_R
    Moo = Tb_R/(10.44-0.0052*Tb_R)
    Mo = newton(f, Moo)

    DSGt = exp(5*(SGo-SG))-1                                            # Eq 13
    ft = DSGt*(-0.362456/Tb_R**0.5+(.0398285-0.948125/Tb_R**0.5)*DSGt)  # Eq 12
    Tc = Tco*((1+2*ft)/(1-2*ft))**2                                     # Eq 11

    DSGv = exp(4*(SGo**2-SG**2))-1                                      # Eq 16
    fv = DSGv*(0.46659/Tb_R**0.5+(-0.182421+3.01721/Tb_R**0.5)*DSGv)    # Eq 15
    Vc = Vco*((1+2*fv)/(1-2*fv))**2                                     # Eq 14

    DSGp = exp(0.5*(SGo-SG))-1                                          # Eq 19
    fp = DSGp*(2.53262-46.1955/Tb_R**0.5-0.00127885*Tb_R + (
        -11.4277+252.14/Tb_R**0.5+0.00230535*Tb_R)*DSGp)                # Eq 18
    Pc = Pco*Tc/Tco*Vco/Vc*((1+2*fp)/(1-2*fp))**2                       # Eq 17

    DSGm = exp(5*(SGo-SG))-1                                            # Eq 23
    x = abs(0.012342-0.328086/Tb_R**0.5)                                # Eq 22
    fm = DSGm*(x+(-0.0175691+0.193168/Tb_R**0.5)*DSGm)                  # Eq 21
    M = exp(log(Mo)*((1+2*fm)/(1-2*fm))**2)                             # Eq 20

    prop = {}
    prop["Tb"] = unidades.Temperature(Tb)
    prop["SG"] = unidades.Dimensionless(SG)
    prop["M"] = unidades.Dimensionless(M)
    prop["Tc"] = unidades.Temperature(Tc, "R")
    prop["Pc"] = unidades.Pressure(Pc, "psi")
    prop["Vc"] = unidades.SpecificVolume(Vc/M, "ft3lb")
    return prop


def prop_Sancet(M):
    """Calculate petroleum fractions properties with the Sancet (2007)
    correlations using only the molecular weight as input parameter

    Parameters
    ------------
    M: float
        Molecular weight, [-]

    Returns
    -------
    prop : A dict with the calculated properties
        Tc: Critic temperature, [ºR]
        Pc: Critic pressure, [psi]
        Tb: Boiling temperature, [ºR]

    Examples
    --------
    Example 2.1 from [9]_: C7+ fraction with M=150 and SG=0.78

    >>> p = prop_Sancet(150)
    >>> "%.0f %.0f %.0f" % (p["Tc"].R, p["Pc"].psi, p["Tb"].R)
    '1133 297 828'

    References
    ----------
    [5] .. Sancet, J., Heavy Faction .. math::`C_{7+}` Characterization for
        PR-EOS. In: SPE 113026, 2007 SPE Annual Conference, November 11–14,
        Anaheim, CA 2007.
    [9] .. Tarek Ahmed. Equations of State and PVT Analysis: Applications for
        Improved Reservoir Modeling, 2nd Edition. Gulf Professional Publishing,
        2016, ISBN 9780128015704
    """
    Pc = 82.82 + 653*exp(-0.007427*M)                                   # Eq 16
    Tc = -778.5 + 383.5*log(M-4.075)                                    # Eq 17
    Tb = 194 + 0.001241*Tc**1.869                                       # Eq 18

    prop = {}
    prop["M"] = unidades.Dimensionless(M)
    prop["Pc"] = unidades.Pressure(Pc, "psi")
    prop["Tc"] = unidades.Temperature(Tc, "R")
    prop["Tb"] = unidades.Temperature(Tb, "R")
    return prop


def prop_Silva_Rodriguez(M):
    """Calculate petroleum fractions properties with the Silva-Rodriguez
    (1992) correlations

    Parameters
    ------------
    M: float
        Molecular weight, [-]

    Returns
    -------
    prop : A dict with the calculated properties
        Tb: Boiling temperature, [K]
        SG: Specific gravity, [-]

    Examples
    --------
    Example 2.1 from [9]_: C7+ fraction with M=150 and SG=0.78

    >>> p = prop_Silva_Rodriguez(150)
    >>> "%.0f %.4f" % (p["Tb"].R, p["SG"])
    '839 0.7982'

    References
    ----------
    [16] .. Silva, M.B., Rodriguez, F., Automatic fitting of equations of state
        for phase behavior matching. Paper SPE 23703, Society of Petroleum
        Engineers, Richardson, TX, 1992.
    [9] .. Tarek Ahmed. Equations of State and PVT Analysis: Applications for
        Improved Reservoir Modeling, 2nd Edition. Gulf Professional Publishing,
        2016, ISBN 9780128015704
    """
    Tb = 447.08725*log(M/64.2576)                                       # Eq A4
    SG = 0.132467*log(Tb)+0.0116483                                     # Eq A5

    prop = {}
    prop["M"] = unidades.Dimensionless(M)
    prop["SG"] = unidades.Dimensionless(SG)
    prop["Tb"] = unidades.Temperature(Tb, "F")
    return prop


def Tb_Soreide(M, SG):
    """Calculate petroleum fractions boiling temperature with the Soreide
    (1989) correlations

    Parameters
    ------------
    M: float
        Molecular weight, [-]
    SG: float
        Specific gravity, [-]

    Returns
    -------
    Tb : float
        Boiling temperature, [K]

    References
    ----------
    [17] .. Soreide, I. Improved Phase Behavior Predictions of Petroleum
        Reservoir Fluids From a Cubic Equation of State. Doctor of engineering
        dissertation. Norwegian Institute of Technology, Trondheim, 1989.
    """
    Tb = 1928.3 - 1.695e5*SG**3.266/M**0.03522*exp(
        -4.922e-3*M-4.7685*SG+3.462e-3*M*SG)                          # Eq 3.59
    return unidades.Temperature(Tb, "R")


def vc_Hall_Yarborough(M, SG):
    """Calculate petroleum fractions critic volume with the Hall-Yarborough
    (1971) correlation

    Parameters
    ------------
    M: float
        Molecular weight, [-]
    SG: float
        Specific gravity, [-]

    Returns
    -------
    vc : float
        Boiling temperature, [K]

    References
    ----------
    [18] .. Hall, K. R., and L. Yarborough. New Simple Correlation for
        Predicting Critical Volume. Chemical Engineering (November 1971): 76.
    """
    Vc = 0.025*M**1.15/SG**0.7985
    return unidades.SpecificVolume(Vc/M, "ft3lb")


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


def prop_Riazi(SG, tita, val):
    """Calculate petroleum fractions properties using the Riazi correlation.
    This correlation is recommendated for heavy weight fractions C20-C50.

    Parameters
    ------------
    SG: float
        Specific gravity, [-]
    tita : string
        Name of second known property
    val float
        known property value

    Returns
    -------
    prop : A dict with the calculated properties
        Tc: Temperatura crítica, [ºR]
        Pc: Presión crítica, [psi]
        Vc: Volumen crítico, [ft³/lb]
        Tb: Temperatura fusión, [ºR]
        d20: Liquid density at 20ºC, [g/cm³]
        I: Huang characterization factor, [-]

    Notes
    -----
    The available input properties for tita are:
        Tb: Temperatura fusión, [ºR]
        M: Peso molecular, [-]

    The critic volume is calculate in mole base, so be careful to correct with
    molecular weight

    References
    ----------
    .. [29] Riazi, M. R. Characterization and Properties of Petroleum
        Fractions. ASTM manual series MNL50, 2005
    """

    if tita not in ["M", "Tb"]:
        raise NotImplementedError(QApplication.translate(
            "pychemqt", "Undefined input pair"))

    # Convert input Tb in Kelvin to Rankine to use in the correlation
    if tita == "Tb":
        val = unidades.K2R(val)

    # Table 2.9
    par = {
      "Tc": {
        "Tb": [35.9416, -6.9e-4, -1.4442, 4.91e-4, 0.7293, 1.2771],
        "M": [218.9592, -3.4e-4, -0.40852, -2.5e-5, 0.331, 0.8136]},

      "Pc": {
        "Tb": [6.9575, -0.0135, -0.3129, 9.174e-3, 0.6791, -0.6807],
        "M": [8.2365e4, -9.04e-3, -3.3304, 0.01006, -0.9366, 3.1353]},

      "Vc": {
        "Tb": [6.1677e10, -7.583e-3, -28.5524, 0.01172, 1.20493, 17.2074],
        "M": [9.703e6, -9.512e-3, -15.8092, 0.01111, 1.08283, 10.5118]},

      "Tb": {
        "M": [9.3369, 1.65e-4, 1.4103, -7.5152e-4, 0.5369, -0.7276]},

      "I": {
        "Tb": [3.2709e-3, 8.4377e-4, 4.59487, -1.0617e-3, 0.03201, -2.34887],
        "M": [1.2419e-2, 7.27e-4, 3.3323, -8.87e-4, 6.438e-3, -1.61166]},

      "d20": {
        "Tb": [0.997, 2.9e-4, 5.0425, -3.1e-4, -0.00929, 1.01772],
        "M": [1.04908, 2.9e-4, -7.339e-2, -3.4e-4, 3.484e-3, 1.05015]}}

    prop = {}
    prop[tita] = _unit(tita, val)

    for key in par.keys():
        if key != tita and (tita != "Tb" or key != "M"):
            a, b, c, d, e, f = par[key][tita]
            x = a*val**e*SG**f*exp(b*val+c*SG+d*val*SG)               # Eq 2.46

            unit = ""
            if key == "Pc":
                unit = "bar"

            val = _unit(key, x, unit)
            prop[key] = val

    return prop


def prop_Riazi_Alsahhaf(Tb, M, rho):
    """Calculate petroleum fractions properties with the Riazi-AlSahhaf (1998)
    correlation with the boiling temperature and liquid density at 20ºC as
    input paramters.

    Parameters
    ------------
    Tb : float
        Normal boiling temperature, [K]
    rho : float
        Liquid density at 20ºC and 1 atm, [g/cm³]

    Returns
    -------
    prop : A dict with the calculated properties
        Tc: Critic temperature, [K]
        Pc: Critic pressure, [MPa]
        Vc: Critic volume, [cm³/g]

    Notes
    -----
    This correlation generalized the Riazi-Daubert method to non-polar
    compounds and gases.

    References
    ----------
    [26] .. Riazi, M. R., Al-Sahhaf, T. A. and Sl-Shammari M. A. A Generalized
        Method for Estimation of Critical Constants. Fluid Phase Equilibria,
        Vol. 147, 1998, pp. 1-6
    """
    a = [1.60193, 10.74145, -8.84800]
    b = [0.00558, 0.07434, -0.03632]
    c = [-0.00112, -0.00047, -0.00547]
    d = [-0.52398, -2.10482, 0.16629]
    e = [0.00104, 0.00508, -0.00028]
    f = [-0.06403, -1.18869, 0.04660]
    g = [0.93857, -0.66773, 2.00241]
    h = [-0.00085, -0.01154, 0.00587]
    i = [0.28290, 1.53161, -0.96608]

    prop = {}
    prop["Tb"] = unidades.Temperature(Tb)
    prop["d"] = unidades.Density(d, "gcc")

    for j, key in enumerate(["Tc", "Pc", "Vc"]):
        x = exp(a[j] + b[j]*M + c[j]*Tb + d[j]*rho + e[j]*Tb*rho) * \
            M**f[j] * Tb**(g[j]+h[j]*M) * rho**i[j]
        val = _unit(key, x, "")
        prop[key] = val

    return prop


def prop_Riazi_Alsahhaf_PNA(M, cmp):
    """Calculate petroleum fractions properties with the Riazi-AlSahhaf PNA
    correlation with the molecular weight as input parameter.
    The procedure calculate the phisical properties for paraphins, naphthenes
    and aromatics constituents of fraction.

    Parameters
    ------------
    M : float
        Molecular weight, [-]
    cmp : string
        Hydrocarbon type, P (paraffins), N (naphthenes) or A (aromatics)

    Returns
    -------
    prop : A dict with the calculated properties
        Tf: Freezing temperature, [K]
        Tb: Boiling temperature, [K]
        SG: Specific gravity, [-]
        d20: Density at 20ºC, [g/cm³]
        Tc: Critic temperature, [K]
        Pc: Critic pressure, [bar]
        Vc: Critic volume, [cm³/g]
        w: Acentric factor, [-]
        I: Refractive index parameter, [-]
        sigma: Surface tension, [dyn/cm]

    References
    ----------
    [27] .. Riazi, M. R. and A1-Sahhaf, T. A. Physical Properties of Heavy
        Petroleum Fractions and Crude Oils. Fluid Phase Equilibria, Vol. 117,
        1996, pp. 217-224.
    [28] .. Riazi, M. R. and A1-Sahhaf, T. Physical Properties of n-Alkanes and
        n-Alkyl Hydrocarbons: Application to Petroleum Mixtures. Industrial and
        Engineering Chemistry Research, Vol. 34, 1995, pp. 4145-4148.
    """
    if cmp == "P":
        tita = [397, 1070, 0.85, 0.859, 0.2833, 1.15, 0, 0.26, 0.3, 33.2]
        a = [6.5096, 6.98291, 92.22793, 88.01379, 87.6593, -0.41966, 4.65757,
             -3.50532, -3.06826, 5.29577]
        b = [0.14187, 0.02013, 89.82301, 85.7446, 86.62167, 0.02436, 0.13423,
             1.5e-6, -1.04987, 0.61653]
        c = [0.470, 2/3, 0.01, 0.01, 0.01, 0.58, 0.5, 2.38, 0.2, 0.32]
    elif cmp == "N":
        tita = [370, 1028, 0.853, 0.857, 0.283, 1.2, 0, -0.255, 0.3, 30.6]
        a = [6.52504, 6.95649, 97.72532, 85.1824, 87.55238, 0.06765, 7.25857,
             -3.18846, -8.25682, 14.17595]
        b = [0.04945, 0.02239, 95.73589, 83.65758, 86.97556, 0.13763, 1.13139,
             0.1658, -5.33934, 7.02549]
        c = [2/3, 2/3, 0.01, 0.01, 0.01, 0.35, 0.26, 0.5, 0.08, 0.12]
    elif cmp == "A":
        tita = [375, 1015, -0.8562, -0.854, -0.2829, 1.03, 0, -0.22, 0, 30.4]
        a = [6.53599, 6.91062, 224.7257, 238.791, 137.0918, -0.29875, 9.77968,
             -1.43083, -14.97, 1.98292]
        b = [0.04912, 0.02247, 218.518, 232.315, 135.433, 0.06814, 3.07555,
             0.12744, -9.48345, -0.0142]
        c = [2/3, 2/3, 0.01, 0.01, 0.01, 0.5, 0.15, 0.5, 0.08, 1.0]

    prop = {}
    prop["M"] = unidades.Dimensionless(M)

    properties = ["Tf", "Tb", "SG", "d20", "I", "Tc", "Pc", "Vc", "w", "sigma"]
    for i, key in enumerate(properties):
        x = tita[i]-exp(a[i]-b[i]*M**c[i])

        # Make conversion where necessary
        if key == "Tc":
            x = prop["Tb"]/x
        elif key in ["Pc", "w"]:
            x = -x
        elif key == "Vc":
            x = 1/x

        # Define units
        unit = ""
        if key == "Pc":
            unit = "bar"
        elif key == "Vc":
            unit = "ccg"
        elif key == "sigma":
            unit = "dyncm"

        val = _unit(key, x, unit)
        prop[key] = val

    return prop


# Zc Methods
def Zc_Hougen(w):
    """Calculate the critical compressibility factdor Zc using the Hougen
    correlation (1959)

    Parameters
    ------------
    w : float
        Acentric factor, [-]

    Returns
    -------
    Zc : float
        Critical compressibility factor, [-]

    References
    ----------
    [22] .. Hougen, O. A., K. M. Watson, and R. A. Ragatz. Chemical Process
        Principles, 2nd ed. New York: Wiley, 1959, p. 577.
    """
    Zc = 1/(1.28*w+3.41)

    return unidades.Dimensionless(Zc)


def Zc_Reid(w):
    """Calculate the critical compressibility factdor Zc using the Reid
    Prausnith and Sherwood (1977) correlation

    Parameters
    ------------
    w : float
        Acentric factor, [-]

    Returns
    -------
    Zc : float
        Critical compressibility factor, [-]

    References
    ----------
    [23] .. Reid, R., J. M. Prausnitz, and T. Sherwood. The Properties of Gases
        and Liquids, 3rd ed. New York: McGraw-Hill, 1977, p. 21.
    """
    Zc = 0.2918 - 0.0928*w

    return unidades.Dimensionless(Zc)


def Zc_Salerno(w):
    """Calculate the critical compressibility factdor Zc using the Salerno
    correlation (1985)

    Parameters
    ------------
    w : float
        Acentric factor, [-]

    Returns
    -------
    Zc : float
        Critical compressibility factor, [-]

    References
    ----------
    [21] .. Salerno, S., Cascella, M., May, D., Watson, P., Tassios, D.
        Prediction of Vapor Pressures and Saturated Volumes with a Simple Cubic
        Equation of State: Part I. A Reliable Data Base, Fluid Phase
        Equilibria, Volume 27, 1986, Pages 15-34
    """
    Zc = 0.291 - 0.08*w - 0.016*w**2                                     # Eq 3

    return unidades.Dimensionless(Zc)


def Zc_Nath(w):
    """Calculate the critical compressibility factdor Zc using the Nath
    correlation (1985)

    Parameters
    ------------
    w : float
        Acentric factor, [-]

    Returns
    -------
    Zc : float
        Critical compressibility factor, [-]

    References
    ----------
    [24] .. Nath, J. Acentric Factor and the Critical Volumes for Normal
        Fluids. Industrial Engineering and Chemical. Fundamentals 21, no. 3
        (1985): 325–326.
    """
    Zc = 0.2908 - 0.0825*w                                              # Eq 1

    return unidades.Dimensionless(Zc)


def Zc_Lee_Kesler(w):
    """Calculate the critical compressibility factdor Zc using the Lee-Kesler
    correlation (1975)

    Parameters
    ------------
    w : float
        Acentric factor, [-]

    Returns
    -------
    Zc : float
        Critical compressibility factor, [-]

    References
    ----------
    [25] .. Lee, B. I. and Kesler, M. G., A Generalized Thermodynamic
        Correlation Based on Three-Parameter Corresponding States. American
        Institute of Chemical Engineers Journal, Vot. 21, 1975
    """
    Zc = 0.2905 - 0.085*w                                               # Eq 21

    return unidades.Dimensionless(Zc)


# Gas compresibility factor
def Z_Papay(Tr, Pr):
    """Calculate gas compressibility factor using the correlation of Papay
    (1985)

    Parameters
    ------------
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

    References
    ----------
    .. [30] Papay, J.A. Termelestechnologiai Parameterek Valtozasa a Gazlelepk
        Muvelese Soran. OGIL MUSZ, Tud, Kuzl. [Budapest], 1985. pp. 267–273.
    .. [9] Tarek Ahmed. Equations of State and PVT Analysis: Applications for
        Improved Reservoir Modeling, 2nd Edition. Gulf Professional Publishing,
        2016, ISBN 9780128015704
    """
    Z = 1 - 3.53*Pr/10**(0.9813*Tr) + 0.274*Pr**2/10**(0.8157*Tr)
    return unidades.Dimensionless(Z)


def Z_Hall_Yarborough(Tr, Pr):
    """Calculate gas compressibility factor using the correlation of Hall &
    Yarborough (1973)

    Parameters
    ------------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]

    References
    ----------
    [31] .. Hall, K. R., and L. Yarborough. A New Equation of State for
        Z-factor Calculations. Oil and Gas Journal (June 18, 1973): 82–92.
    [9] .. Tarek Ahmed. Equations of State and PVT Analysis: Applications for
        Improved Reservoir Modeling, 2nd Edition. Gulf Professional Publishing,
        2016, ISBN 9780128015704
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


def Z_Dranchuk_Abu_Kassem(Tr, Pr):
    """Calculate gas compressibility factor using the correlation of Dranchuk-
    Abu-Kassem (1975)

    Parameters
    ------------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]

    Raises
    ------
    NotImplementedError : If input isn't in limit
        * 1 ≤ Tr ≤ 3
        * 0.2 ≤ Pr ≤ 30

    References
    ----------
    [32] .. Dranchuk, P. M., and J. H. Abu-Kassem. Calculate of Z-factors for
        Natural Gases Using Equations-of-State. Journal of Canadian Petroleum
        Technology (July–September 1975): 34–36
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


def Z_Dranchuk_Purvis_Robinson(Tr, Pr):
    """Calculate gas compressibility factor using the correlation of Dranchuk-
    Purvis-Robinson (1974)

    Parameters
    ------------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]

    Raises
    ------
    NotImplementedError : If input isn't in limit
        * 1.05 ≤ Tr ≤ 3
        * 0.2 ≤ Pr ≤ 30

    References
    ----------
    .. [33] Dranchuk, P.M., Purvis, R.A., Robinson, D.B. Computer Calculations
        of Natural Gas Compressibility Factors Using the Standing and Katz
        Correlation. Technical Series, no. IP 74–008. Institute of Petroleum,
        Alberta, Canada, 1974.
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


def Z_Brill_Beggs(Tr, Pr):
    """Calculate gas compressibility factor using the correlation of Brill-
    Beggs (1974)

    Parameters
    ------------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]

    Raises
    ------
    NotImplementedError : If input isn't in limit
        * 1.15 ≤ Tr ≤ 2.4
        * 0.2 ≤ Pr ≤ 15

    References
    ----------
    .. [34] Brill, J .P. and Beggs, H .D. Two-Phase Flow in Pipes. University
        of Tulsa, INTERCOMP Course, The Hague, 1974
    .. [35] Kumar, N. Compressibility factors for natural and sour reservoir
        gases by correlations and cubic equations of state. Thesis of master
        of science in Petroleum Engineering, 2004, Texas Tech University.
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


def Z_Gopal(Tr, Pr):
    """Calculate gas compressibility factor using the correlation of Gopal
    (1974)

    Parameters
    ------------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]

    Raises
    ------
    NotImplementedError : If input isn't in limit
        * 1.05 ≤ Tr ≤ 3
        * 0.2 ≤ Pr ≤ 15

    References
    ----------
    .. [36] Gopal, V.N. Gas Z-Factor Equations Developed for Computer. Oil and
        Gas J. (Aug. 8, 1977) 58-60
    .. [35] Kumar, N. Compressibility factors for natural and sour reservoir
        gases by correlations and cubic equations of state. Thesis of master
        of science in Petroleum Engineering, 2004, Texas Tech University.
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


def Z_ShellOil(Tr, Pr):
    """Calculate gas compressibility factor using the Shell Oil Company
    correlation (2004)

    Parameters
    ------------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]

    References
    ----------
    .. [35] Kumar, N. Compressibility factors for natural and sour reservoir
        gases by correlations and cubic equations of state. Thesis of master
        of science in Petroleum Engineering, 2004, Texas Tech University.
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


def Z_Sarem(Tr, Pr):
    """Calculate gas compressibility factor using the Sarem (1969) correlation

    Parameters
    ------------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]

    Raises
    ------
    NotImplementedError : If input isn't in limit
        * 1.05 ≤ Tr ≤ 2.95
        * 0.1 ≤ Pr ≤ 14.9

    References
    ----------
    .. [37] Sarem, A.M. Z-Factor Equation Developed for Use in Digital
        Computers. Oil and Gas J. (Sept. 18, 1961) 118
    """
    # Check input in range of validity
    if Tr < 1.05 or Tr > 2.95 or Pr < 0.1 or Pr > 14.9:
        raise NotImplementedError("Incoming out of bound")

    x = (2.*Pr-15)/14.8
    y = (2.*Tr-4)/1.9
    Aij = [
        [2.1433504, .0831762, -.0214670, -.0008714, .0042846, -.0016595],
        [.3312352, -.1340361, .0668810,  -.0271743,  .0088512,  -.002152],
        [.1057287,  -.0503937,  .0050925,  .0105513,  -.0073182,  .0026960],
        [.0521840,  .0443121,  -.0193294,  .0058973,  .0015367,  -.0028327],
        [.0197040,  -.0263834, .019262,  -.0115354,  .0042910,  -.0081303],
        [.0053096,  .0089178,  -.0108948,  .0095594,  -.0060114, .0031175]]

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


def Z_Leung(Tr, Pr):
    """Calculate gas compressibility factor using the Leung (1964) correlation

    Parameters
    ------------
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
    The correlation is in cited referencence, the parameters are least square
    fitting by Leung.

    Raises
    ------
    NotImplementedError : If input isn't in limit
        * 1.1 ≤ Tr ≤ 2.6
        * 0.5 ≤ Pr ≤ 11

    References
    ----------
    .. [38] Dranckuk, P.M., Quon, D. A General Solution of the Equations
        Describing Steady State Turbulent Compressible Flow in Circular
        Conduits. Journal of Canadian Petroleum Technology 3(2):60-65, 1964
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


def Z_Burnett(Tr, Pr):
    """Calculate gas compressibility factor using the Burnett (1979)
    correlation

    Parameters
    ------------
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
    The correlation is in cited referencence, the parameters are least square
    fitting by Leung.

    Raises
    ------
    NotImplementedError : If input isn't in limit
        * 1.3 ≤ Tr ≤ 3
        * 0.2 ≤ Pr ≤ 4

    References
    ----------
    .. [39] Burnett, R.R. Calculator gives compressibility factors. Oil & Gas
        Journal, June 11, 1979, pp. 70-74.
    .. [40] Takacs., G. Comparing Methods for Calculating Z-factor. Oil & Gas
        Journal, May 15, 1989, pp. 43-46.
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


def Z_Sanjari_Lay(Tr, Pr):
    """Calculate gas compressibility factor using the Sanjari-Lay (2012)
    correlation

    Parameters
    ------------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]

    Raises
    ------
    NotImplementedError : If input isn't in limit
        * 1.01 ≤ Tr ≤ 3
        * 0.01 ≤ Pr ≤ 15

    References
    ----------
    .. [41] Sanjari E, Lay E.N. An accurate empirical correlation for
        predicting natural gas compressibility factors. Journal of Natural Gas
        Chemistry 21(2012):184-188.
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


def Z_Heidaryan_Salarabadi(Tr, Pr):
    """Calculate gas compressibility factor using the Heidaryan-Salarabadi-
    Moghadasi (2010) correlation

    Parameters
    ------------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]

    Raises
    ------
    NotImplementedError : If input isn't in limit
        * 1.2 ≤ Tr ≤ 3
        * 0.2 ≤ Pr ≤ 15

    References
    ----------
    .. [43] Heidaryan, E., Salarabadi, A., Moghadasi, J. A novel correlation
        approach for prediction of natural gas compressibility factor. J. Nat.
        Gas Chem. 19 (2) 2010, 189–192.
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


def Z_Heidaryan_Moghadasi(Tr, Pr):
    """Calculate gas compressibility factor using the Heidaryan-Moghadasi-
    Rahimi (2010) correlation

    Parameters
    ------------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]

    Raises
    ------
    NotImplementedError : If input isn't in limit
        * 1.2 ≤ Tr ≤ 3
        * 0.2 ≤ Pr ≤ 15

    References
    ----------
    .. [42] Heidaryan E, Moghadasi J, Rahimi M. New correlations to predict
        natural gas viscosity and compressibility factor. Journal of Petroleum
        Science and Engineering 73 (2010):67-72
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


def Z_Azizi(Tr, Pr):
    """Calculate gas compressibility factor using the Azizi-Behbahani-Isazadeh
    (2010) correlation

    Parameters
    ------------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]

    Raises
    ------
    NotImplementedError : If input isn't in limit
        * 1.1 ≤ Tr ≤ 2
        * 0.2 ≤ Pr ≤ 11

    References
    ----------
    .. [44] Azizi N, Behbahani R, Isazadeh M A. An efficient correlation for
        calculating compressibility factor of natural gases. Journal of Natural
        Gas Chemistry 19 (2010) 642-645
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


def Z_Shokir(Tr, Pr):
    """Calculate gas compressibility factor using the Shokir-Awad-Quraishi
    (2012) correlation

    Parameters
    ------------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]

    References
    ----------
    .. [46] Shokir, Eissa M.El-M., El-Awad, Musaed N., Al-Quraishi,
        Adulhrahman A., Al-Mahdy, Osama A. Compressibility factor model of
        sweet, sour, and condensate gases using genetic programming. Chem. Eng.
        Res. Des. 90 (2012), 785-792.
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


def Z_Bahadori(Tr, Pr):
    """Calculate gas compressibility factor using the Bahadori-Mokhatab-Towler
    (2007) correlation

    Parameters
    ------------
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

    Raises
    ------
    NotImplementedError : If input isn't in limit
        * 1.05 ≤ Tr ≤ 2.4
        * 0.2 ≤ Pr ≤ 16

    References
    ----------
    .. [47] Bahadori, A., Mokhatab, S., Towler, B.F. Rapidly estimating natural
        gas compressibility factor. J. Nat. Gas Chem. 16 (4) 2007, 349-353
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


def Z_Londono_DAK(Tr, Pr, pure=False):
    """Calculate gas compressibility factor using the Londono-Archer-Blasingame
    (2007) correlation
    This method implement the Dranckuk-Abu-Kassem optimized version from paper

    Parameters
    ------------
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

    References
    ----------
    .. [48] Londono, F.E., Archer, R.A., Blasingame, T.A. Correlations for
        hydrocarbon-gas viscosity and gas density-validation and correlation
        of behavior using a large-scale database. SPE Reserv. Evalu. Eng. 8 (6)
        2005, 561–572.
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


def Z_Londono_NS(Tr, Pr, pure=False):
    """Calculate gas compressibility factor using the Londono-Archer-Blasingame
    (2007) correlation
    This method implement the Nishiumi-Saito optimized version from paper

    Parameters
    ------------
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

    References
    ----------
    .. [48] Londono, F.E., Archer, R.A., Blasingame, T.A. Correlations for
        hydrocarbon-gas viscosity and gas density-validation and correlation
        of behavior using a large-scale database. SPE Reserv. Evalu. Eng. 8 (6)
        2005, 561–572.
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


def Z_Hall_Iglesias(Tr, Pr):
    """Calculate gas compressibility factor using the correlation of Hall-
    Iglesias (2007). This is a extension of Hall-Yarborough correlation to
    reduced values of temperature.

    Parameters
    ------------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]

    References
    ----------
    .. [45] Hall, K.R., Iglesias-Silva, G.A.. Improved equations for the
        Standing-Katz tables. Hydrocarb. Process 86 (4) 2007, 107-110.
    .. [49] Chankalani, A., Mae'soumi, A., Sameni, A. An Intelligent Approach
        for Optimal Prediction of Gas Deviation Factor Using Particle Swarm
        Optimization and Genetic Algorithm. Journal of Natural Gas Science and
        Engineering 14 (2013) 132-143
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


def Z_Elsharkawy(Tr, Pr):
    """Calculate gas compressibility factor using the Elsharkawy (2003)
    correlation

    Parameters
    ------------
    Tr : float
        Reduced temperature [-]
    Pr : float
        Reduced pressure [-]

    Returns
    -------
    Z : float
        Gas compressibility factor [-]

    References
    ----------
    .. [50] Elsharkawy, A.M. Efficient methods for calculations of
        compressibility, density and viscosity of natural gases. Fluid Phase
        Equilibria 218:1 (2004) 1-13
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


# Distillation curves interconversion methods Ref. [29] Pag 117
def D86_TBP_Riazi(Ti, Xi=None, reverse=False):
    """Interconversion between D86 and TBP at atmospheric pressure using the
    Riazi correlation

    Parameters
    ------------
    Ti : list
        Distillation data with the D86 distillation temperature
    Xi : list
        Volumetric cut point range

    reverse : boolean
        To do the reverse conversion

    Returns
    -------
    TBP : list
        Calculated distillation data

    Example
    -------
    Example 3.2 from [29]_
    >>> X = [0, 0.1, 0.3, 0.5, 0.7, 0.9]
    >>> TBP = [10, 71.1, 143.3, 204.4, 250.6, 291.7]
    >>> TBP_K = [t+273.15 for t in TBP]
    >>> D86 = D86_TBP_Riazi(TBP_K, X, reverse=True)
    >>> D86_C = [t-273.15 for t in D86]
    >>> "%.0f" % D86_C[0]
    '32'

    Example 3.3 from [29]_
    >>> X = [0, 0.1, 0.3, 0.5, 0.7, 0.9]
    >>> D86 = [165.6, 173.7, 193.3, 206.7, 222.8, 242.8]
    >>> D86_K = [t+273.15 for t in D86]
    >>> TDB = D86_TBP_Riazi(D86_K, X)
    >>> TDB_C = [t-273.15 for t in TDB]
    >>> "%.1f %.1f %.1f %.1f %.1f %.1f" % tuple(TDB_C)
    '134.2 157.4 190.3 209.0 230.2 254.7'

    References
    ----------
    .. [51] Riazi, M. R. and Daubert, T. E. Analytical Correlations
        Interconvert Distillation Curve Types. Oil & Gas Journal, Vol. 84,
        1986, August 25, pp. 50-57
    .. [29] Riazi, M. R. Characterization and Properties of Petroleum
        Fractions. ASTM manual series MNL50, 2005
    """
    # Define the default value for volumetric volume
    if Xi is None:
        Xi = [0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95]

    a = [0.9177, 0.5564, 0.76517, 0.9013, 0.8821, 0.9552, 0.8177]
    b = [1.0019, 1.09, 1.0425, 1.0176, 1.0226, 1.011, 1.0355]
    Xmin = [0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95]
    Xmax = [0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 1.0001]

    # Calculate parameters
    A = []
    B = []
    for x in Xi:
        for xmin, xmax, ai, bi in zip(Xmin, Xmax, a, b):
            if xmin <= x < xmax:
                A.append(ai)
                B.append(bi)
                break

    # Calculate desired distillation data
    TBP = []
    if reverse:
        for ai, bi, T in zip(A, B, Ti):
            TBP.append(unidades.Temperature((1/ai)**(1/bi) * T**(1/bi)))
    else:
        for ai, bi, T in zip(A, B, Ti):
            TBP.append(unidades.Temperature(ai*T**bi))
    return TBP


def D86_TBP_Daubert(Ti, Xi=None, T50=None, reverse=False):
    """Interconversion between D86 and TBP at atmospheric pressure using the
    Daubert correlation, API procedure 3A1.1

    Parameters
    ------------
    Ti : list
        Distillation data with the T0, T10, T30, T50, T70, T90, T95
    Xi : list
        Volumetric cut point range
    T50 : float
        D86 distillation at 50% point temperature
    reverse : boolean
        To do the reverse conversion

    Returns
    -------
    TBP : list
        Calculated distillation data

    Example
    -------
    Example from [20]_
    >>> X = [0.1, 0.3, 0.5, 0.7, 0.9]
    >>> D86 = [350, 380, 404, 433, 469]
    >>> D86_K = [unidades.F2K(t) for t in D86]
    >>> TDB = D86_TBP_Daubert(D86_K, X, D86_K[2])
    >>> TDB_F = [t.F for t in TDB]
    >>> "%.1f %.1f %.1f %.1f %.1f" % tuple(TDB_F)
    '316.5 372.6 411.2 451.2 496.7'

    Inverse case
    >>> X = [0.1, 0.3, 0.5, 0.7, 0.9]
    >>> TDB = [321, 371, 409, 447, 469]
    >>> TDB_K = [unidades.F2K(t) for t in TDB]
    >>> D86 = D86_TBP_Daubert(TDB_K, X, TDB_K[2], reverse=True)
    >>> D86_F = [t.F for t in D86]
    >>> "%.1f %.1f" % (D86_F[1], D86_F[2])
    '378.4 401.9'

    Example 3.3 from [29]_
    >>> X = [0, 0.1, 0.3, 0.5, 0.7, 0.9]
    >>> D86 = [165.6, 173.7, 193.3, 206.7, 222.8, 242.8]
    >>> D86_K = [t+273.15 for t in D86]
    >>> TDB = D86_TBP_Daubert(D86_K, X)
    >>> TDB_C = [t-273.15 for t in TDB]
    >>> "%.1f %.1f %.1f %.1f %.1f %.1f" % tuple(TDB_C)
    '133.5 154.2 189.2 210.7 232.9 258.2'

    References
    ----------
    .. [52] Daubert, T. E. Petroleum Fraction Distillation Interconversion.
        Hydrocarbon Processing, Vol. 73, No. 9, 1994, pp. 75-78
    .. [20] API. Technical Data book: Petroleum Refining 6th Edition
    .. [29] Riazi, M. R. Characterization and Properties of Petroleum
        Fractions. ASTM manual series MNL50, 2005
    """
    # Convert to Fahrenheit the input distillation temperature
    Ti = [unidades.K2F(t) for t in Ti]

    # Define the default value for volumetric volume
    if Xi is None:
        Xi = [0, 0.1, 0.3, 0.5, 0.7, 0.9, 1]

    if T50 is None:
        T50 = Ti[3]
    else:
        T50 = unidades.K2F(T50)

    a = [7.4012, 4.9004, 3.0305, 2.5282, 3.0419, 0.11798]
    b = [0.60244, 0.71644, 0.80076, 0.82002, 0.75497, 1.6606]
    Xmin = [0, 0.1, 0.3, 0.5, 0.7, 0.9]
    Xmax = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0001]

    # Calculate parameters
    A = []
    B = []
    for x in Xi:
        for xmin, xmax, ai, bi in zip(Xmin, Xmax, a, b):
            if xmin <= x < xmax:
                A.append(ai)
                B.append(bi)
                break

    # Calculate the desviation coefficient
    if reverse:
        TBP50 = exp(log(T50/0.87180)/1.0258)
        Y = [exp(log((Ti[i+1]-T)/A[i])/B[i]) for i, T in enumerate(Ti[:-1])]
    else:
        TBP50 = 0.87180*T50**1.0258
        Y = [A[i]*(Ti[i+1]-T)**B[i] for i, T in enumerate(Ti[:-1])]

    # Calculate desired distillation data
    TBP = [TBP50]*len(Ti)
    for i, T in enumerate(TBP):
        for y, x in zip(Y[i:], Xi[i:]):
            if Xi[i] < 0.5 and x < 0.5:
                TBP[i] -= y
        for y, x in zip(Y[:i], Xi[:i]):
            if Xi[i] >= 0.5 and x >= 0.5:
                TBP[i] += y

    return [unidades.Temperature(t, "F") for t in TBP]


def SD_D86_Riazi(Ti, Xi=None, F=None):
    """Interconversion between SD and D86

    Parameters
    ------------
    Ti : list
        Distillation data with the D86 distillation temperature
    Xi : list
        Volumetric cut point range
    F : float
        Parameter, :math: $F = 0.01411*SD_{10%}^0.05434*SD_{50%}^0.6147

    Returns
    -------
    D86 : list
        Calculated distillation data

    Example
    -------
    Example 3.5 from [29]_
    >>> X = [0.1, 0.3, 0.5, 0.7, 0.9]
    >>> SD = [33.9, 64.4, 101.7, 140.6, 182.2]
    >>> SD_K = [t+273.15 for t in SD]
    >>> F = 0.01411*SD_K[0]**0.05434*SD_K[2]**0.6147
    >>> D86 = SD_D86_Riazi(SD_K, X, F)
    >>> D86_C = [t-273.15 for t in D86]
    >>> "%.1f %.1f %.1f %.1f %.1f" % tuple(D86_C)
    '53.2 70.9 96.0 131.3 168.3'

    References
    ----------
    .. [51] Riazi, M. R. and Daubert, T. E. Analytical Correlations
        Interconvert Distillation Curve Types. Oil & Gas Journal, Vol. 84,
        1986, August 25, pp. 50-57
    .. [29] Riazi, M. R. Characterization and Properties of Petroleum
        Fractions. ASTM manual series MNL50, 2005
    """
    # Define the default value for volumetric volume
    if Xi is None:
        Xi = [0, 0.1, 0.3, 0.5, 0.7, 0.9, 1]

    if F is None:
        F = 0.01411*Ti[1]**0.05434*Ti[3]**0.6147

    # Possible typo in referencen [20]_, the a_50% parameter with a 10 factor
    a = [5.1764, 3.7452, 4.2749, 18.445, 1.0751, 1.0849, 1.7991]
    b = [0.7445, 0.7944, 0.7719, 0.5425, 0.9867, 0.9834, 0.9007]
    c = [0.2879, 0.2671, 0.345, 0.7132, 0.0486, 0.0354, 0.0625]
    Xmin = [0, 0.1, 0.3, 0.5, 0.7, 0.9, 1]
    Xmax = [0.1, 0.3, 0.5, 0.7, 0.9, 1, 1.0001]

    # Calculate parameters
    A = []
    B = []
    C = []
    for x in Xi:
        for xmin, xmax, ai, bi, ci in zip(Xmin, Xmax, a, b, c):
            if xmin <= x < xmax:
                A.append(ai)
                B.append(bi)
                C.append(ci)
                break

    D86 = []
    for ai, bi, ci, T in zip(A, B, C, Ti):
        D86.append(unidades.Temperature(ai*T**bi*F**ci))
    return D86


def SD_D86_Daubert(Ti, Xi=None, SD50=None):
    """Interconversion between gas chromatography (ASTM D2887) and ASTM D86 at
    atmospheric pressure using the Daubert correlation, API procedure 3A3.2

    Parameters
    ------------
    Ti : list
        Distillation data with the T0, T10, T30, T50, T70, T90, T95
    Xi : list
        Volumetric cut point range
    SD50 : float
        SD86 distillation at 50% point temperature

    Returns
    -------
    TBP : list
        Calculated distillation data

    Example
    -------
    Example from [20]_
    >>> X = [0, 0.1, 0.3, 0.5, 0.7, 0.9, 1]
    >>> SD = [77, 93, 148, 215, 285, 360, 408]
    >>> SD_K = [unidades.F2K(t) for t in SD]
    >>> D86 = SD_D86_Daubert(SD_K, X)
    >>> D86_F = [t.F for t in D86]
    >>> "%.1f %.1f %.1f %.1f %.1f %.1f %.1f" % tuple(D86_F)
    '121.3 128.2 154.8 206.3 270.6 334.0 367.5'

    Example 3.5 from [29]_
    >>> X = [0.1, 0.3, 0.5, 0.7, 0.9]
    >>> SD = [33.9, 64.4, 101.7, 140.6, 182.2]
    >>> SD_K = [t+273.15 for t in SD]
    >>> D86 = SD_D86_Daubert(SD_K, X, SD_K[2])
    >>> D86_C = [t-273.15 for t in D86]
    >>> "%.1f %.1f %.1f %.1f %.1f" % tuple(D86_C)
    '53.5 68.2 96.9 132.6 167.8'

    References
    ----------
    .. [52] Daubert, T. E. Petroleum Fraction Distillation Interconversion.
        Hydrocarbon Processing, Vol. 73, No. 9, 1994, pp. 75-78
    .. [20] API. Technical Data book: Petroleum Refining 6th Edition
    .. [29] Riazi, M. R. Characterization and Properties of Petroleum
        Fractions. ASTM manual series MNL50, 2005
    """
    # Convert to Fahrenheit the input distillation temperature
    Ti = [unidades.K2F(t) for t in Ti]

    # Define the default value for volumetric volume
    if Xi is None:
        Xi = [0, 0.1, 0.3, 0.5, 0.7, 0.9]

    if SD50 is None:
        SD50 = Ti[3]
    else:
        SD50 = unidades.K2F(SD50)

    e = [0.3047, 0.06069, 0.07978, 0.14862, 0.30785, 2.6029]
    f = [1.1259, 1.5176, 1.5386, 1.4287, 1.2341, 0.65962]
    Xmin = [0, 0.1, 0.3, 0.5, 0.7, 0.9]
    Xmax = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0001]

    # Calculate parameters
    E = []
    F = []
    for x in Xi:
        for xmin, xmax, ei, fi in zip(Xmin, Xmax, e, f):
            if xmin <= x < xmax:
                E.append(ei)
                F.append(fi)
                break

    # Calculate the desviation coefficient
    U = [E[i]*(Ti[i+1]-T)**F[i] for i, T in enumerate(Ti[:-1])]

    # Calculate desired distillation data
    D86_50 = 0.77601*SD50**1.0395
    D86 = [D86_50]*len(Ti)
    for i, T in enumerate(D86):
        for y, x in zip(U[i:], Xi[i:]):
            if Xi[i] < 0.5 and x < 0.5:
                D86[i] -= y
        for y, x in zip(U[:i], Xi[:i]):
            if Xi[i] >= 0.5 and x >= 0.5:
                D86[i] += y

    return [unidades.Temperature(t, "F") for t in D86]


def D86_EFV(Ti, Xi=None, SG=None, reverse=False):
    """Interconversion between D86 and EFV

    Parameters
    ------------
    Ti : list
        Distillation data with the D86 distillation temperature
    Xi : list
        Volumetric cut point range
    SG : float
        Specific gravity, [-]
    reverse : boolean
        To do the reverse conversion

    Returns
    -------
    EFV : list
        Calculated distillation data

    Example
    -------
    Example 3.2 from [29]_
    >>> X = [0, 0.1, 0.3, 0.5, 0.7, 0.9]
    >>> TBP = [10, 71.1, 143.3, 204.4, 250.6, 291.7]
    >>> TBP_K = [t+273.15 for t in TBP]
    >>> D86 = D86_TBP_Riazi(TBP_K, X, reverse=True)
    >>> EFV = D86_EFV(D86, X, SG=0.7862)
    >>> EFV_C = [t-273.15 for t in EFV]
    >>> "%.0f" % EFV_C[0]
    '68'

    References
    ----------
    .. [51] Riazi, M. R. and Daubert, T. E. Analytical Correlations
        Interconvert Distillation Curve Types. Oil & Gas Journal, Vol. 84,
        1986, August 25, pp. 50-57
    .. [29] Riazi, M. R. Characterization and Properties of Petroleum
        Fractions. ASTM manual series MNL50, 2005
    """
    # Define the default value for volumetric volume
    if Xi is None:
        Xi = [0, 0.1, 0.3, 0.5, 0.7, 0.9, 1]

    if SG is None:
        SG = 0.08342*Ti[1]**0.10731*Ti[3]**0.26288

    a = [2.9747, 1.4459, 0.8506, 3.268, 8.2873, 10.6266, 7.9952]
    b = [0.8466, 0.9511, 1.0315, 0.8274, 0.6874, 0.6529, 0.6949]
    c = [0.4209, 0.1287, 0.0817, 0.6214, 0.934, 1.1025, 1.0737]
    Xmin = [0, 0.1, 0.3, 0.5, 0.7, 0.9, 1]
    Xmax = [0.1, 0.3, 0.5, 0.7, 0.9, 1, 1.0001]

    # Calculate parameters
    A = []
    B = []
    C = []
    for x in Xi:
        for xmin, xmax, ai, bi, ci in zip(Xmin, Xmax, a, b, c):
            if xmin <= x < xmax:
                A.append(ai)
                B.append(bi)
                C.append(ci)
                break

    # Calculate desired distillation data
    EFV = []
    if reverse:
        for ai, bi, ci, T in zip(A, B, C, Ti):
            EFV.append(unidades.Temperature((T/ai/SG**ci)**(1./bi)))
    else:
        for ai, bi, ci, T in zip(A, B, C, Ti):
            EFV.append(unidades.Temperature(ai*T**bi*SG**ci))
    return EFV


def SD_TBP(Ti, Xi=None, SD50=None):
    """Interconversion between gas chromatography (ASTM D2887) and TBP at
    atmospheric pressure using the Daubert correlation, API procedure 3A3.1

    Parameters
    ------------
    Ti : list
        Distillation data with the T0, T10, T30, T50, T70, T90, T95
    Xi : list
        Volumetric cut point range
    SD50 : float
        SD86 distillation at 50% point temperature

    Returns
    -------
    TBP : list
        Calculated distillation data

    Example
    -------
    Example from [20]_
    >>> X = [0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95]
    >>> SD = [293, 305, 324, 336, 344, 359, 369]
    >>> SD_K = [unidades.F2K(t) for t in SD]
    >>> TDB = SD_TBP(SD_K, X)
    >>> TDB_F = [t.F for t in TDB]
    >>> "%.1f %.1f %.1f %.1f %.1f %.1f %.1f" % tuple(TDB_F)
    '322.2 327.7 332.4 336.0 339.6 350.1 357.4'

    Example 3.4 from [29]_
    >>> X = [0.1, 0.3, 0.5, 0.7, 0.9]
    >>> SD = [151.7, 162.2, 168.9, 173.3, 181.7]
    >>> SD_K = [t+273.15 for t in SD]
    >>> TDB = SD_TBP(SD_K, X, SD_K[2])
    >>> TDB_C = [t-273.15 for t in TDB]
    >>> "%.1f %.1f %.1f %.1f %.1f" % tuple(TDB_C)
    '164.3 166.9 168.9 170.9 176.8'

    References
    ----------
    .. [52] Daubert, T. E. Petroleum Fraction Distillation Interconversion.
        Hydrocarbon Processing, Vol. 73, No. 9, 1994, pp. 75-78
    .. [20] API. Technical Data book: Petroleum Refining 6th Edition
    .. [29] Riazi, M. R. Characterization and Properties of Petroleum
        Fractions. ASTM manual series MNL50, 2005
    """
    # Convert to Fahrenheit the input distillation temperature
    Ti = [unidades.K2F(t) for t in Ti]

    # Define the default value for volumetric volume
    if Xi is None:
        Xi = [0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95]

    if SD50 is None:
        SD50 = Ti[3]
    else:
        SD50 = unidades.K2F(SD50)

    c = [0.15779, 0.011903, 0.05342, 0.19861, 0.31531, 0.97476, 0.02172]
    d = [1.4296, 2.0253, 1.6988, 1.3975, 1.2938, 0.8723, 1.9733]
    Xmin = [0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95]
    Xmax = [0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 1.0001]

    # Calculate parameters
    C = []
    D = []
    for x in Xi:
        for xmin, xmax, ci, di in zip(Xmin, Xmax, c, d):
            if xmin <= x < xmax:
                C.append(ci)
                D.append(di)
                break

    # Calculate the desviation coefficient
    W = [C[i]*(Ti[i+1]-T)**D[i] for i, T in enumerate(Ti[:-1])]

    # Calculate desired distillation data
    TBP = [SD50]*len(Ti)
    for i, T in enumerate(TBP):
        for y, x in zip(W[i:], Xi[i:]):
            if Xi[i] < 0.5 and x < 0.5:
                TBP[i] -= y
        for y, x in zip(W[:i], Xi[:i]):
            if Xi[i] >= 0.5 and x >= 0.5:
                TBP[i] += y

    return [unidades.Temperature(t, "F") for t in TBP]


def D1160_TBP_10mmHg(Ti, Xi=None, reverse=False):
    """Interconversion between ASTM D1160 and TBP distillation curve at 10 mmHg

    Parameters
    ------------
    Ti : list
        Distillation data with the T0, T10, T30, T50, T70, T90, T95
    Xi : list
        Volumetric cut point range

    Returns
    -------
    TBP : list
        Calculated distillation data

    Example
    -------
    Example 3.6 from [29]_
    >>> X = [0.1, 0.3, 0.5, 0.7, 0.9]
    >>> D1160 = [150, 205, 250, 290, 350]
    >>> D1160_K = [t+273.15 for t in D1160]
    >>> TDB = D1160_TBP_10mmHg(D1160_K, X)
    >>> TDB_C = [t-273.15 for t in TDB]
    >>> "%.1f %.1f %.0f %.0f %.0f" % tuple(TDB_C)
    '146.6 200.9 250 290 350'

    References
    ----------
    .. [53] Edmister, W. C. and Okamoto, K. K. Applied Hydrocarbon
        Thermodynamics, Part 13: Equilibrium Flash Vaporization Correlations
        for Heavy Oils Under Subatmospheric Pressures," Petroleum Refiner,
        Vol. 38, No. 9, 1959, pp. 271-288.
    .. [29] Riazi, M. R. Characterization and Properties of Petroleum
        Fractions. ASTM manual series MNL50, 2005
    """
    # Define the default value for volumetric volume
    if Xi is None:
        Xi = [0.1, 0.3, 0.5, 0.7, 0.9, 1]

    TBP = []
    for i, T in enumerate(Ti):
        if Xi[i] >= 0.5:
            TBP.append(T)
        else:
            DT = Ti[i+1]-T
            if Xi[i] >= 0.1:
                F = 0.3 + 1.2775*DT - 5.539e-3*DT**2 + 2.7486e-5*DT**3
            else:
                F = 2.2566*DT - 266.2e-4*DT**2 + 1.4093e-4*DT**3
            TBP.append(Ti[i+1]-F)

    return [unidades.Temperature(t) for t in TBP]


def Tb_Pressure(T, P, Kw=None, reverse=False):
    """Conversion of boiling point at Sub or super atmospheric pressure to the
    normal boiling point, API Procedure 5A1.19

    Parameters
    ------------
    T : float
        Temperature, [K]
    P : float
        Pressure, [Pa]
    Kw: float
        Watson characterization factor, [-]
    reverse : boolean
        Do the inverse calculation from normal boiling point to other pressure

    Returns
    -------
    TBP : list
        Calculated distillation data

    Example
    -------
    Example from [20]_
    >>> T = unidades.Temperature(365, "F")
    >>> P = unidades.Pressure(10, "mmHg")
    >>> TDB = Tb_Pressure(T, P, 12.5)
    >>> "%.0f" % TDB
    '602'

    References
    ----------
    .. [20] API. Technical Data book: Petroleum Refining 6th Edition
    .. [29] Riazi, M. R. Characterization and Properties of Petroleum
        Fractions. ASTM manual series MNL50, 2005
    """
    p = unidades.Pressure(P, "Pa")
    if p.mmHg < 2:
        Q = (6.76156-0.987672*log10(p.mmHg))/(3000.538-43*log10(p.mmHg))
    elif p.mmHg < 760:
        Q = (5.994296-0.972546*log10(p.mmHg))/(2663.129-95.76*log10(p.mmHg))
    else:
        Q = (6.412631-0.989679*log10(p.mmHg))/(2770.085-36.*log10(p.mmHg))

    if reverse:
        if T < 367 or not Kw:
            F = 0
        elif T < 478:
            F = -3.2985+0.009*T*(Kw-12)
        else:
            F = 1
        Tb_ = T-1.3889*F*log10(p.atm)
        Tb = Tb_/(748.1*Q-Tb_*(0.3861*Q-0.00051606))
    else:
        Tb_ = 748.1*Q*T/(1+T*(0.3861*Q-0.00051606))
        if Tb_ < 367 or not Kw:
            F = 0
        elif Tb_ < 478:
            F = -3.2985+0.009*Tb_*(Kw-12)
        else:
            F = 1
        Tb = Tb_+1.3889*F*log10(p.atm)

    return unidades.Temperature(Tb)


def curve_Predicted(x, T):
    """Fill the missing point of a distillation curve

    Parameters
    ------------
    x : list
        Array with mole/volume/weight fraction, [-]
    T : list
        Array with boiling point temperatures, [K]

    Returns
    -------
    TBP : list
        Calculated distillation data

    Example
    -------
    Example 4.7 from [54]_
    >>> xw = [0.261, 0.254, 0.183, 0.14, 0.01, 0.046, 0.042, 0.024, 0.015, \
        0.009, 0.007]
    >>> x = [sum(xw[:i+1]) for i, xi in enumerate(xw)]
    >>> T = [365, 390, 416, 440, 461, 482, 500, 520, 539, 556, 573]
    >>> To, A, B = curve_Predicted(x, T)
    >>> "%.0f %.4f %.4f" % (To, A, B)
    '350 0.1679 1.2586'

    References
    ----------
    .. [54] Riazi, M. R. Distribution Model for Properties of Hydrocarbon-Plus
        Fractions. Industrial and Engineering Chemistry Research, Vol. 28,
        1989, pp. 1731-1735.
    """
    x = array(x)
    T = array(T)
    p0 = [T[0], 0.1, 1]

    def errf(p, xw, Tb):
        return _Tb_Predicted(p, xw) - Tb

    p, cov, info, mesg, ier = leastsq(errf, p0, args=(x, T), full_output=True)

    ss_err = (info['fvec']**2).sum()
    ss_tot = ((T-T.mean())**2).sum()
    r = 1-(ss_err/ss_tot)
    return p, r


def _Tb_Predicted(par, x):
    """Calculate a specific point in the distillation curve"""
    return par[0]+par[0]*(par[1]/par[2]*log(1/(1-x)))**(1./par[2])


# Others properties
def PourPoint(SG, Tb, v100=None):
    """Calculate the pour point of petroleum fractions using the procedure
    2B8 from API technical databook, pag. 235
    The pour point is the lowest temperature at which it will flow or can be
    poured.

    Parameters
    ------------
    SG : float
        Specific gravity, [-]
    Tb : float
        Mean average boiling point, [K]
    v100 : float
        Kinematic viscosity at 100ºF, [cSt]

    Returns
    -------
    PP : float
        Pour point, [ªR]

    Raises
    ------
    NotImplementedError : If input isn't in limit
        * 0.8 ≤ SG ≤ 1
        * 800ºR ≤ Tb ≤ 1500ºR
        * 2cSt ≤ v100 ≤ 960cSt

    Examples
    --------
    >>> T = unidades.Temperature(972, "R")
    >>> "%0.0f" % PourPoint(0.839, T, 3).R
    '458'

    References
    ----------
    [20] .. API. Technical Data book: Petroleum Refining 6th Edition
    """
    # Convert input Tb in Kelvin to Rankine to use in the correlation
    Tb_R = unidades.K2R(Tb)

    # Check input in range of validity
    if SG < 0.8 or SG > 1 or Tb_R < 800 or Tb_R > 1500:
        raise NotImplementedError("Incoming out of bound")
    elif v100 is not None and (v100 < 2 or v100 > 960):
        raise NotImplementedError("Incoming out of bound")

    if v100 is not None:
        # Eq 2B8.1-1
        PP = 753 + 136*(1-exp(-0.15*v100)) - 572*SG + 0.0512*v100 + 0.139*Tb_R
    else:
        # Eq 2B8.1-2
        PP = 3.85e-8*Tb_R**5.49*10**-(0.712*Tb.R**0.315+0.133*SG) + 1.4

    return unidades.Temperature(PP, "R")


def AnilinePoint(SG, Tb):
    """Calculate the aniline point of petroleum fractions using the procedure
    2B9 from API technical databook, pag. 237
    The aniline point is the lowest temperature at which a petroleum fraction
    is completely miscible with an equal volume of distilled aniline.

    Parameters
    ------------
    SG : float
        Specific gravity, [-]
    Tb : float
        Mean average boiling point, [ºR]

    Returns
    -------
    AP : float
        Aniline point, [ªR]

    Raises
    ------
    NotImplementedError : If input isn't in limit
        * 0.7 ≤ SG ≤ 1
        * 200ºF ≤ Tb ≤ 1100ºF

    Examples
    --------
    >>> T = unidades.Temperature(570.2, "F")
    >>> "%0.1f" % AnilinePoint(0.8304, T).R
    '635.5'

    References
    ----------
    [20] .. API. Technical Data book: Petroleum Refining 6th Edition
    """
    # Convert input Tb in Kelvin to Rankine to use in the correlation
    Tb_F = unidades.K2F(Tb)
    Tb_R = unidades.K2R(Tb)
    Kw = Tb_R**(1/3)/SG

    # Check input in range of validity
    if SG < 0.7 or SG > 1 or Tb_F < 200 or Tb_F > 1100:
        raise NotImplementedError("Incoming out of bound")

    # Eq 2B9.1-1
    AP = -1253.7 - 0.139*Tb_R + 107.8*Kw + 868.7*SG

    return unidades.Temperature(AP, "R")


def SmokePoint(SG, Tb):
    """Calculate the smoke point of petroleum fractions using the procedure 2B10
    from API technical databook, pag. 239
    The smoke point is the height in millimeters of the flame that is
    produced in a lamp at standard conditions without causing smoking.

    Parameters
    ------------
    SG : float
        Specific gravity, [-]
    Tb : float
        Mean average boiling point, [ºF]

    Returns
    -------
    SP : float
        Cloud point, [mm]

    Raises
    ------
    NotImplementedError : If input isn't in limit
        * 0.7 ≤ SG ≤ 0.86
        * 200ºF ≤ Tb ≤ 550ºF

    Examples
    --------
    >>> T = unidades.Temperature(414.5, "F")
    >>> "%0.1f" % SmokePoint(0.853, T).mm
    '16.7'

    References
    ----------
    [20] .. API. Technical Data book: Petroleum Refining 6th Edition
    """
    # Convert input Tb in Kelvin to Rankine to use in the correlation
    Tb_F = unidades.K2F(Tb)
    Tb_R = unidades.K2R(Tb)
    Kw = Tb_R**(1/3)/SG

    # Check input in range of validity
    if SG < 0.7 or SG > 0.86 or Tb_F < 200 or Tb_F > 550:
        raise NotImplementedError("Incoming out of bound")

    # Eq 2B10.1-1
    SP = exp(-1.028+0.474*Kw-0.00168*Tb_R)

    return unidades.Length(SP, "mm")


def FreezingPoint(SG, Tb):
    """Calculate freezing point of petroleum fractions using the procedure 2B11
    from API technical databook, pag. 241
    The freezing point is the temperature at which solid crystals formed on
    cooling disappear as the temperature is raised.

    Parameters
    ------------
    SG : float
        Specific gravity, [-]
    Tb : float
        Mean average boiling point, [ºR]

    Returns
    -------
    FP : float
        Cloud point, [ºR]

    Raises
    ------
    NotImplementedError : If input isn't in limit
        * 0.74 ≤ SG ≤ 0.9
        * 725ºR ≤ Tb ≤ 1130ºR

    Examples
    --------
    >>> T = unidades.Temperature(874.5, "R")
    >>> "%0.0f" % CloudPoint(0.799, T).R
    '417'

    References
    ----------
    [20] .. API. Technical Data book: Petroleum Refining 6th Edition
    """
    # Convert input Tb in Kelvin to Rankine to use in the correlation
    Tb_R = unidades.K2R(Tb)
    Kw = Tb_R**(1/3)/SG

    # Check input in range of validity
    if SG < 0.74 or SG > 0.90 or Tb_R < 725 or Tb_R > 1130:
        raise NotImplementedError("Incoming out of bound")

    # Eq 2B11.1-1
    FP = -2390.42 + 1826*SG + 122.49*Kw - 0.135*Tb_R

    return unidades.Temperature(FP, "R")


def CloudPoint(SG, Tb):
    """Calculate cloud point of petroleum fractions using the procedure 2B12
    from API technical databook, pag. 243
    The cloud point of is the temperature at which its solid paraffin content
    begins to solidify and separate in tiny crystals, causing the oil to appear
    cloudy.

    Parameters
    ------------
    SG : float
        Specific gravity, [-]
    Tb : float
        Mean average boiling point, [ºR]

    Returns
    -------
    CP : float
        Cloud point, [ºR]

    Raises
    ------
    NotImplementedError : If input isn't in limit
        * 0.77 ≤ SG ≤ 0.93
        * 800ºR ≤ Tb ≤ 1225ºR

    Examples
    --------
    >>> T = unidades.Temperature(811.5, "R")
    >>> "%0.1f" % CloudPoint(0.787, T).R
    '383.4'

    References
    ----------
    [20] .. API. Technical Data book: Petroleum Refining 6th Edition
    """
    # Convert input Tb in Kelvin to Rankine to use in the correlation
    Tb_R = unidades.K2R(Tb)

    # Check input in range of validity
    if SG < 0.77 or SG > 0.93 or Tb_R < 800 or Tb_R > 1225:
        raise NotImplementedError("Incoming out of bound")

    # Eq 2B12.1-1
    CP = 10**(-7.41+5.49*log10(Tb_R)-0.712*Tb_R**0.315-0.133*SG)

    return unidades.Temperature(CP, "R")


def CetaneIndex(API, Tb):
    """Calculate cetane index of petroleum fractions using the procedure 2B13
    from API technical databook, pag. 245
    Cetane index is the number equal to the porcentage of cetane in a blend of
    cetane and α-methyl naphthalene having the same ignition quality as a
    sample of the petroleum fraction.

    Parameters
    ------------
    API : float
        API gravity, [-]
    Tb : float
        Mean average boiling point, [ºF]

    Returns
    -------
    CI : float
        Cetane Index, [-]

    Raises
    ------
    NotImplementedError : If input isn't in limit
        * 27 ≤ API ≤ 47
        * 360ºF ≤ Tb ≤ 700ºF

    Examples
    --------
    >>> T = unidades.Temperature(617, "F")
    >>> "%0.1f" % CetaneIndex(32.3, T)
    '57.1'

    References
    ----------
    [20] .. API. Technical Data book: Petroleum Refining 6th Edition
    """
    # Convert input Tb in Kelvin to Fahrenheit to use in the correlation
    Tb_F = unidades.K2F(Tb)

    # Check input in range of validity
    if API < 27 or API > 47 or Tb_F < 360 or Tb_F > 700:
        raise NotImplementedError("Incoming out of bound")

    # Eq 2B13.1-1
    CI = 415.26 - 7.673*API + 0.186*Tb_F + 3.503*API*log10(Tb_F) - \
        193.816*log10(Tb.F)

    return unidades.Dimensionless(CI)


class Petroleo(newComponente):
    """Class to define a heavy oil fraction with a unknown composition

    Parámetros:
        nombre: nombre del componente, generalmente una fracción de petróleo
        M: peso_molecular
        Tb: temperaura normal de ebullición, K
        SG: gravedad específica a 60ºF
        API: gravedad API
        CH: relación C/H
        I: parámetro de Huang
        n: indice de refracción
        Nc: número de carbonos de la fracción
        Kw: factor característico de watson
        v100: viscosidad cinemática a 100ºF
        v210: viscosidad cinemática a 210ºF
        H: Porcentaje de contenido en hidrogeno de la fracción
        S: Porcentaje de contenido en azufre de la fracción
        N: Porcentaje de contenido en nitrógeno de la fracción

        Curvas de destilación:
            -D86: Distribución de temperaturas de ebullición de la fracción según la norma ASTM D86
            -TBP: True boiling point
            -SD: Simulated distillation según la norma ASTM D2887
            -EFV: Equilibrium Flash Vaporization
            -D1160: Curva específicas a baja presión
            -P_dist: Presión de la curva de destilación, en mmHg
            -T_dist: array con los valores de las fracciones en porcentaje (peso o volumen) de los datos de la curva de destilación

        Otros parámeros más específicos
        Aplicadas a gasolinas:
            oleffin: porcentaje de contenido de olefinas
            TML: concentración de tetrametilplomo añadido en ml/galón UK
            TEL: concentración de tetraetilplomo añadido en ml/galón UK

    Opciones de definición por prioridad:
        1   -   Tb y SG
        2   -   M y SG
        3   -   Tb y I
        4   -   M y I
        5   -   Tb y CH
        6   -   M y CH
        7   -   v100 y I
        8   -   Nc (opción muy poco precisa)
        9   -   curva de destilación

    API es una forma equivalente alternativa a SG
    Kw es una forma alternativa de indicar SG o Tb
    n es una forma equivalente alternativa a I
    """

    kwargs={"name": "",
                    "M": 0.0,
                    "Tb": 0.0,
                    "SG": 0.0,
                    "API": 0.0,
                    "CH": 0.0,
                    "n": 0.0,
                    "I": 0.0,
                    "Nc": 0,
                    "Kw": 0.0,
                    "v100": 0.0,
                    "v210": 0.0,
                    "H": 0.0,
                    "S": 0.0,
                    "N": 0.0,

                    "P_dist": 0.0,
                    "T_dist": 0.0,
                    "D86": [],
                    "TBP": [],
                    "EFV": [],
                    "SD": [],
                    "D1160": [],

                    "oleffin": 0.0,
                    "TEL": 0.0,
                    "TML": 0.0,
                    }

    status=0
    _bool=False
    msg=""
    definicion=0

    def __init__(self, **kwargs):
        self.Preferences=ConfigParser()
        self.Preferences.read(conf_dir+"pychemqtrc")
        self.__call__(**kwargs)

    def __call__(self, **kwargs):
        self.kwargs.update(kwargs)
        if kwargs:
            self._bool=True
        if self.isCalculable():
            self.calculo()

    def isCalculable(self):
        """
        1   -   Tb y SG
        2   -   M y SG
        3   -   Tb y I
        4   -   M y I
        5   -   Tb y CH
        6   -   M y CH
        7   -   v100 y I
        8   -   Nc (opción muy poco precisa)
        9   -   curva de destilación
        """
        self.hasSG=self.kwargs["SG"] or self.kwargs["API"] or (self.kwargs["Kw"] and self.kwargs["Tb"])
        self.hasRefraction=self.kwargs["n"] or self.kwargs["I"]
        self.hasCurve=self.kwargs["D86"] or self.kwargs["TBP"] or self.kwargs["EFV"] or self.kwargs["D1160"]

        #Tipo Definición
        if self.kwargs["Tb"] and self.hasSG:
            self.definicion=1
        elif self.kwargs["M"] and self.hasSG:
            self.definicion=2
        elif self.kwargs["Tb"] and self.hasRefraction:
            self.definicion=3
        elif self.kwargs["M"] and self.hasRefraction:
            self.definicion=4
        elif self.kwargs["Tb"] and self.kwargs["CH"]:
            self.definicion=5
        elif self.kwargs["M"] and self.kwargs["CH"]:
            self.definicion=6
        elif self.kwargs["v100"] and self.hasRefraction:
            self.definicion=7
        elif self.kwargs["Nc"]:
            self.definicion=8
        elif self.hasCurve:
            self.definicion=9
        else:
            self.status=0
            self.msg=QApplication.translate("pychemqt", "Insufficient input")

        if self.definicion:
            self.status=1
            self.msg=""
            return True


    def calculo(self):
        self.formula=""
        if self.kwargs["name"]:
            self.name=str(self.kwargs["name"])
        else:
            self.name=self.__class__.__name__+"_"+time.strftime("%d/%m/%Y-%H:%M:%S")

        self.cp=[]
        self.Vliq=0
        self.rackett=0
        self.Tf=0
        self.Hf=0
        self.Gf=0

        if self.hasCurve==1:
            curva=D86 or TBP or SD or EFV or D1160
            parameters_Curva=curve_Predicted(T_dist, curva)
            curva_normalizada=[]
            for i in [0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 99]:
                if i in T_dist:
                    curva_normalizada.append(unidades.Temperature(curva[T_dist.index(i)]))
                else:
                    curva_normalizada.append(unidades.Temperature(T_Predicted(parameters_Curva, i)))

            # Si SG no está disponible se puede calcular a partir de la curva de destilación Pag 117 Riazi
            if not self.hasSG:
                if self.kwargs["D86"]:
                    self.SG=0.08342*curva_normalizada[2]**0.10731*curva_normalizada[6]**0.26288
                elif self.kwargs["TBP"]:
                    self.SG=0.10431*curva_normalizada[2]**0.12550*curva_normalizada[6]**0.20862
                else:
                    self.SG=0.09138*curva_normalizada[2]**-0.0153*curva_normalizada[6]**0.36844

            if self.kwargs["D86"]:
                self.D86=curva_normalizada
                self.TBP=[D86_TBP_Riazi, D86_TBP_Daubert][self.Preferences.getint("petro", "curva")](self.D86)
                self.EFV=D86_EFV(self.D86, self.SG)
                self.SD=None
                TBP_10mmHg=[Tb_Presion(t, 10./760, reverse=True) for t in self.TBP]
                D1160_10mmHg=D1160_TBP_10mmHg(TBP_10mmHg, reverse=True)
                self.D1160=[Tb_Presion(t, 10./760) for t in D1160_10mmHg]

            elif self.kwargs["TBP"]:
                if P_dist.mmHg == 760:
                    self.TBP=curva_normalizada
                    TBP_10mmHg=[Tb_Presion(t, 10./760, reverse=True) for t in self.TBP]
                elif P_dist.mmHg == 10:
                    self.TBP=[Tb_Presion(t, P_dist) for t in curva_normalizada]
                    TBP_10mmHg=curva_normalizada
                else:
                    self.TBP=[Tb_Presion(t, P_dist) for t in curva_normalizada]
                    TBP_10mmHg=[Tb_Presion(t, 10./760, reverse=True) for t in self.TBP]
                self.D86=[D86_TBP_Riazi, D86_TBP_Daubert][self.Preferences.getint("petro", "curva")](self.TBP, reverse=True)
                self.EFV=D86_EFV(self.D86, SG)
                self.SD=None
                D1160_10mmHg=D1160_TBP_10mmHg(TBP_10mmHg, reverse=True)
                self.D1160=[Tb_Presion(t, 10./760) for t in D1160_10mmHg]

            elif self.kwargs["SD"]:
                self.SD=curva_normalizada
                self.D86=[SD_D86_Riazi, SD_D86_Daubert][self.Preferences.getint("petro", "curva")](self.SD)
                self.TBP=SD_TBP(self.SD)
                self.EFV=D86_EFV(self.D86, self.SG)
                TBP_10mmHg=[Tb_Presion(t, 10./760, reverse=True) for t in self.TBP]
                D1160_10mmHg=D1160_TBP_10mmHg(TBP_10mmHg, reverse=True)
                self.D1160=[Tb_Presion(t, 10./760) for t in D1160_10mmHg]

            elif self.kwargs["EFV"]:
                if P_dist.mmHg==760:
                    self.EFV=curva_normalizada
                else:
                    self.EFV=[Tb_Presion(t, P_dist/760.) for t in curva_normalizada]
                self.D86=D86_EFV(self.EFV, self.SG, reverse=True)
                self.SD=None
                self.TBP=[D86_TBP_Riazi, D86_TBP_Daubert][self.Preferences.getint("petro", "curva")](self.D86)
                TBP_10mmHg=[Tb_Presion(t, 10./760, reverse=True) for t in self.TBP]
                D1160_10mmHg=D1160_TBP_10mmHg(TBP_10mmHg, reverse=True)
                self.D1160=[Tb_Presion(t, 10./760) for t in D1160_10mmHg]

            else:
                if P_dist.mmHg==10:
                    D1160_10mmHg=curva_normalizada
                    self.D1160=[Tb_Presion(t, 10./760) for t in D1160_10mmHg]
                elif P_dist.mmHg==760:
                    self.D1160=curva_normalizada
                    D1160_10mmHg=[Tb_Presion(t, 10./760, reverse=True) for t in self.D1160]
                else:
                    self.D1160=[Tb_Presion(t, P_dist) for t in curva_normalizada]
                    D1160_10mmHg=[Tb_Presion(t, 10./760, reverse=True) for t in self.D1160]
                TBP_10mmHg=D1160_TBP_10mmHg(D1160_10mmHg)
                self.TBP=[Tb_Presion(t, 10./760) for t in TBP_10mmHg]
                self.D86=[D86_TBP_Riazi, D86_TBP_Daubert][self.Preferences.getint("petro", "curva")](self.TBP, reverse=True)
                self.EFV=D86_EFV(self.D86, self.SG)
                self.SD=None


            self.VABP=unidades.Temperature((self.D86[2]+self.D86[4]+self.D86[6]+self.D86[8]+self.D86[10])/5.)
            SL=(self.D86[-2].F-self.D86[2].F)/80.
            self.WABP=unidades.Temperature(self.VABP.F+exp(-3.062123-0.01829*(self.VABP.F-32)**(2./3)+4.45818*SL**0.25), "F")
            self.MABP=unidades.Temperature(self.VABP.F-exp(-0.56379-0.007981*(self.VABP.F-32)**(2./3)+3.04729*SL**(1./3)), "F")
            self.CABP=unidades.Temperature(self.VABP.F-exp(-0.23589-0.06906*(self.VABP.F-32)**0.45+1.8858*SL**0.45), "F")
            self.MeABP=unidades.Temperature(self.VABP.F-exp(-0.94402-0.00865*(self.VABP.F-32)**(2./3)+2.99791*SL**(1./3)), "F")
            if not Tb:
                Tb=self.MeABP


            self.T5=self.D86[1]
            self.T10=self.D86[2]
            self.T30=self.D86[4]
            self.T50=self.D86[6]
            self.T90=self.D86[10]
            self.ReidVP=unidades.Pressure(3.3922-0.02537*self.T5.C-0.070739*self.T10.C+0.00917*self.T30.C-0.0393*self.T50.C+6.8257e-4*self.T10.C**2, "bar")

            A, B, To=parameters_Curva
            E70=100-100*exp(-B/A*(343.15-To)**B/To**B)
            self.VL_12=88.5-0.19*E70-42.5*self.ReidVP.bar
            self.VL_20=90.6-0.25*E70-39.2*self.ReidVP.bar
            self.VL_36=94.7-0.36*E70-32.3*self.ReidVP.bar
            self.CVLI=4.27+0.24*E70+0.069*self.ReidVP.bar
            self.FVI=1000.*self.ReidVP.bar+7.*E70

            #Cálculo de flash point, API procedure 2B7.1 pag 233
            self.FlashPc=unidades.Temperature(0.69*self.T10-118.2, "F")
            self.FlashPo=unidades.Temperature(0.68*self.T10-109.6, "F")

#            Tav=To*(1+(A/B)**(1./B)*gamma(1+1./B))
#            print Tav, self.VABP, self.WABP

        #Cálculo de la composición PNA
            if self.Preferences.getint("petro", "PNA")==0:
                self.xp, self.xn, self.xa=self.PNA_Peng_Robinson()
            elif self.Preferences.getint("petro", "PNA")==1:
                self.xp, self.xn, self.xa=self.PNA_Bergman()
            elif self.Preferences.getint("petro", "PNA")==2:
                self.xp, self.xn, self.xa=self.PNA_Riazi()
            else:
                self.xp, self.xn, self.xa=self.PNA_van_Nes()

            T=self.Tb.C/100
            RONp=92.809-70.97*T-53.*T**2+20.*T**3+10.*T**4
            RONi=(95.927+92.069+109.38+97.652)/4+(-157.53+57.63-38.83-20.8)/4*T+(561-65-26-58)/4*T**2-800/4.*T**3+300/4.*T**4
            RONn=-77.536+471.59*T-418.*T**2+100.*T**3
            RONa=145.668-54.336*T+16.276*T**2
            self.RON=self.xp/2.*(RONp+RONi)+self.xn*RONn+self.xa*RONa
            self.MON=22.5+0.83*self.RON-20.*self.SG-0.12*self.kwargs["oleffin"]+0.5*self.kwargs["TML"]-0.2*self.kwargs["TML"]


        if self.hasSG:
            if self.kwargs["SG"]:
                self.SG=self.kwargs["SG"]
                self.API=141.5/self.SG-131.5
            elif self.kwargs["API"]:
                self.API=self.kwargs["API"]
                self.SG=141.5/(self.API+131.5)
            elif Kw and Tb:
                self.SG=unidades.Temperature(self.kwargs["Tb"]).R**(1./3)/self.kwargs["Kw"]
                self.API=141.5/self.SG-131.5
        else:
            SG=[self.SG_Riazi_Daubert, self.SG_Riazi_Alsahhaf][self.Preferences.getint("petro", "SG")]
            self.SG=SG()

        if self.kwargs["Tb"]:
            self.Tb=unidades.Temperature(self.kwargs["Tb"])
        if self.kwargs["M"]:
            self.M=self.kwargs["M"]
        if not self.kwargs["Tb"]:
            Tb=[self.tb_Riazi_Daubert_ext, self.tb_Riazi_Adwani, self.tb_Edmister, self.tb_Soreide][self.Preferences.getint("petro", "t_ebull")]
            self.Tb=Tb()

        if not self.kwargs["M"]:
            Peso_Molecular=[self.peso_molecular_Riazi_Daubert, self.peso_molecular_Riazi_Daubert_ext, self.peso_molecular_Lee_Kesler, self.peso_molecular_Sim_Daubert, self.peso_molecular_API, self.peso_molecular_Goossens, self.peso_molecular_Twu][self.Preferences.getint("petro", "molecular_weight")]
            self.M=Peso_Molecular()

        self.Nc=self.kwargs["Nc"]
        if self.definicion==7:
            prop = prop_Ahmed(self.kwargs["Nc"])
            self.M = prop["M"]
            self.Tc = prop["Tc"]
            self.Pc = prop["Pc"]
            self.Vc = prop["Vc"]
            self.Tb = prop["Tb"]
            self.f_acent = prop["w"]
            self.SG = prop["SG"]
            self.API=141.5/self.SG-131.5

        if self.hasRefraction:
            if self.kwargs["n"]:
                self.n=self.kwargs["n"]
                self.I=(self.n**2-1)/(self.n**2+2)
            else:
                self.I=self.kwargs["I"]
                self.n=((1+2*self.I)/(1-self.I))**0.5
        else:
            self.I=self.I_Riazi_Daubert_ext()
            self.n=((1+2*self.I)/(1-self.I))**0.5

        if self.kwargs["Kw"]:
            self.watson=self.kwargs["Kw"]
        else:
            self.watson=self.Tb.R**(1./3)/self.SG

        if self.kwargs["CH"]:
            self.CH=self.kwargs["CH"]
        else:
            self.CH=self.CH_Riazi_Daubert_ext()

        Tc=[self.tc_Riazi_Daubert_ext, self.tc_Riazi_Daubert_ext, self.tc_Riazi_Adwani, self.tc_Riazi_Daubert_ext, self.tc_Riazi_Daubert_ext, self.tc_Sim_Daubert, self.tc_Watansiri_Owens_Starling, self.tc_Edmister, self.tc_Magoulas, self.tc_Twu, self.tc_Tsonopoulos][self.Preferences.getint("petro", "critical")]
        Pc=[self.pc_Riazi_Daubert, self.pc_Riazi_Daubert_ext, self.pc_Riazi_Adwani, self.pc_Lee_Kesler, self.pc_cavett, self.pc_Sim_Daubert, self.pc_Watansiri_Owens_Starling, self.pc_Edmister, self.pc_Magoulas, self.pc_Twu, self.pc_Tsonopoulos][self.Preferences.getint("petro", "critical")]
        Vc=[self.vc_Riazi_Daubert, self.vc_Riazi_Daubert_ext, self.vc_Riazi_Adwani, self.vc_Watansiri_Owens_Starling, self.vc_Twu, self.vc_Hall_Yarborough][self.Preferences.getint("petro", "critical")]
        Factor_acentrico=[self.factor_acentrico_Edmister, self.factor_acentrico_Lee_Kesler, self.factor_acentrico_Watansiri_Owens_Starling, self.factor_acentrico_Magoulas][self.Preferences.getint("petro", "f_acent")]
        self.Tc=Tc()
        self.Pc=Pc()
        self.Vc=Vc()
        self.f_acent=Factor_acentrico()
        Zc=[self.Zc_Lee_Kesler, self.Zc_Haugen, self.Zc_Reid, self.Zc_Salerno, self.Zc_Nath][self.Preferences.getint("petro", "Zc")]
        self.Zc=Zc()

        self.Parametro_solubilidad=prop_Riazi_Alsahhaf(9, self.M), "calcc"

        #Calculo de las viscosidades cinemáticas si no indican
        if self.kwargs["v100"]:
            self.v100=unidades.Diffusivity(self.kwargs["v100"])
        else:
            self.v100=self.V100_API()
        if self.kwargs["v210"]:
            self.v210=unidades.Diffusivity(self.kwargs["v210"])
        else:
            self.v210=self.V210_API()

       #Calculo de Viscosity gravity constant (VGC), preferiblemente usando la viscosidad a 100ºF
        if self.kwargs["v210"] and not self.kwargs["v100"]:
            v99SUS=SUS(99+273.15, self.v210.cSt)
            self.VGC=(self.SG-0.24-0.022*log10(v99SUS-35.5))
        else:
            v38SUS=SUS(38+273.15, self.v100.cSt)
            self.VGC=(10*self.SG-1.0752*log10(v38SUS-38))/(10-log10(v38SUS-38))

        self.d20=self.SG-4.5e-3*(2.34-1.9*self.SG)
        self.Ri=self.n-self.d20/2.
        self.m=self.M*(self.n-1.475)
        self.VI=self.Viscosity_Index()

        if self.hasCurve:
            self.CetaneI=48640/self.VABP+473.7*self.SG-456.8
        else:
            self.CetaneI=48640/self.Tb+473.7*self.SG-456.8

        self.PourP = PourPoint(self.SG, self.Tb, self.kwargs["v100"])
        self.CloudP = CloudPoint(self.SG, self.Tb)
        self.FreezingP = FreezingPoint(self.SG, self.Tb)
        self.AnilineP = AnilinePoint(self.SG, self.Tb)
        self.SmokeP = SmokePoint(self.SG, self.Tb)
        self.CetaneI = CetaneIndex(self.API, self.Tb)

        #Cálculo del indice Diesel
        self.DI=self.API*self.AnilineP.F/100.


        if self.kwargs["S"]:
            self.S = self.kwargs["S"]
        else:
            self.S = self.S_Riazi()
        if self.kwargs["H"]:
            self.H = H
        else:
            H = [self.H_Riazi, self.H_Goossens, self.H_ASTM, self.H_Jenkins_Walsh][self.Preferences.getint("petro", "H")]
            self.H = H()
        self.N = self.kwargs["N"]

    def tr(self, T):
        return T/self.Tc

    def pr(self, P):
        return P/self.Pc.atm

    def SG_Riazi_Daubert(self):
        if self.definicion == 3:
            return prop_Riazi_Daubert("Tb", self.Tb, "I", self.I)["SG"]
        elif self.definicion == 4:
            return prop_Riazi_Daubert("M", self.M, "I", self.I)["SG"]
        elif self.definicion == 5:
            return prop_Riazi_Daubert("Tb", self.Tb, "CH", self.CH)["SG"]
        elif self.definicion == 6:
            return prop_Riazi_Daubert("M", self.M, "CH", self.CH)["SG"]
        elif self.definicion == 7:
            return prop_Riazi_Daubert("v1", self.v100, "I", self.I)["SG"]

    def SG_Riazi_Alsahhaf(self):
        if not self.kwargs["M"]:
            return self.SG_Riazi_Daubert()
        else:
            return prop_Riazi_Alsahhaf(2, self.kwargs["M"])

    def peso_molecular_Riazi_Daubert(self):
        return prop_Riazi_Daubert_1980(self.Tb, self.SG)["M"]

    def tc_Riazi_Daubert(self):
        return prop_Riazi_Daubert_1980(self.Tb, self.SG)["Tc"]

    def pc_Riazi_Daubert(self):
        return prop_Riazi_Daubert_1980(self.Tb, self.SG)["Pc"]

    def vc_Riazi_Daubert(self):
        return prop_Riazi_Daubert_1980(self.Tb, self.SG)["Vc"]

    def peso_molecular_Riazi_Daubert_ext(self):
        if self.definicion == 1:
            return prop_Riazi_Daubert("Tb", self.Tb, "SG", self.SG)["M"]
        elif self.definicion == 3:
            return prop_Riazi_Daubert("Tb", self.Tb, "I", self.I)["M"]
        elif self.definicion == 5:
            return prop_Riazi_Daubert("Tb", self.Tb, "CH", self.CH)["M"]
        elif self.definicion == 7:
            return prop_Riazi_Daubert("v1", self.v100, "I", self.I)["M"]

    def tc_Riazi_Daubert_ext(self):
        if self.definicion == 1:
            return prop_Riazi_Daubert("Tb", self.Tb, "SG", self.SG)["Tc"]
        elif self.definicion == 2:
            return prop_Riazi_Daubert("M", self.M, "SG", self.SG)["Tc"]
        elif self.definicion == 3:
            return prop_Riazi_Daubert("Tb", self.Tb, "I", self.I)["Tc"]
        elif self.definicion == 4:
            return prop_Riazi_Daubert("M", self.M, "I", self.I)["Tc"]
        elif self.definicion == 5:
            return prop_Riazi_Daubert("Tb", self.Tb, "CH", self.CH)["Tc"]
        elif self.definicion == 6:
            return prop_Riazi_Daubert("M", self.M, "CH", self.CH)["Tc"]
        elif self.definicion == 7:
            return prop_Riazi_Daubert("v1", self.v100, "I", self.I)["Tc"]

    def pc_Riazi_Daubert_ext(self):
        if self.definicion == 1:
            return prop_Riazi_Daubert("Tb", self.Tb, "SG", self.SG)["Pc"]
        elif self.definicion == 2:
            return prop_Riazi_Daubert("M", self.M, "SG", self.SG)["Pc"]
        elif self.definicion == 3:
            return prop_Riazi_Daubert("Tb", self.Tb, "I", self.I)["Pc"]
        elif self.definicion == 4:
            return prop_Riazi_Daubert("M", self.M, "I", self.I)["Pc"]
        elif self.definicion == 5:
            return prop_Riazi_Daubert("Tb", self.Tb, "CH", self.CH)["Pc"]
        elif self.definicion == 6:
            return prop_Riazi_Daubert("M", self.M, "CH", self.CH)["Pc"]
        elif self.definicion == 7:
            return prop_Riazi_Daubert("v1", self.v100, "I", self.I)["Pc"]

    def vc_Riazi_Daubert_ext(self):
        if self.definicion == 1:
            return prop_Riazi_Daubert("Tb", self.Tb, "SG", self.SG)["Vc"]
        elif self.definicion == 2:
            return prop_Riazi_Daubert("M", self.M, "SG", self.SG)["Vc"]
        elif self.definicion == 3:
            return prop_Riazi_Daubert("Tb", self.Tb, "I", self.I)["Vc"]
        elif self.definicion == 4:
            return prop_Riazi_Daubert("M", self.M, "I", self.I)["Vc"]
        elif self.definicion == 5:
            return prop_Riazi_Daubert("Tb", self.Tb, "CH", self.CH)["Vc"]
        elif self.definicion == 6:
            return prop_Riazi_Daubert("M", self.M, "CH", self.CH)["Vc"]
        elif self.definicion == 7:
            return prop_Riazi_Daubert("M", self.v100, "I", self.I)["Vc"]

    def tb_Riazi_Daubert_ext(self):
        if self.definicion == 2:
            return prop_Riazi_Daubert("M", self.M, "SG", self.SG)["Tb"]
        elif self.definicion == 4:
            return prop_Riazi_Daubert("M", self.M, "I", self.I)["Tb"]
        elif self.definicion == 6:
            return prop_Riazi_Daubert("M", self.M, "CH", self.CH)["Tb"]
        elif self.definicion == 7:
            return prop_Riazi_Daubert("v1", self.v100, "I", self.I)["Tb"]

    def I_Riazi_Daubert_ext(self):
        if self.definicion == 1:
            return prop_Riazi_Daubert("Tb", self.Tb, "SG", self.SG)["I"]
        elif self.definicion == 2:
            return prop_Riazi_Daubert("M", self.M, "SG", self.SG)["I"]
        elif self.definicion == 5:
            return prop_Riazi_Daubert("Tb", self.Tb, "CH", self.CH)["I"]
        elif self.definicion == 6:
            return prop_Riazi_Daubert("M", self.M, "CH", self.CH)["I"]

    def CH_Riazi_Daubert_ext(self):
        if self.definicion == 1:
            return prop_Riazi_Daubert("Tb", self.Tb, "SG", self.SG)["CH"]
        elif self.definicion == 2:
            return prop_Riazi_Daubert("M", self.M, "SG", self.SG)["CH"]
        elif self.definicion == 3:
            return prop_Riazi_Daubert("Tb", self.Tb, "I", self.I)["CH"]
        elif self.definicion == 4:
            return prop_Riazi_Daubert("M", self.M, "I", self.I)["CH"]
        elif self.definicion == 7:
            return prop_Riazi_Daubert("v1", self.v100, "I", self.I)["CH"]

    def tb_Riazi_Adwani(self):
        return unidades.Temperature(prop_Riazi_Adwani(0, 1, self.M, self.SG))

    def tc_Riazi_Adwani(self):
        if self.definicion==1:
            tc=prop_Riazi_Adwani(1, 0, self.Tb, self.SG)
        elif self.definicion==2:
            tc=prop_Riazi_Adwani(1, 1, self.M, self.SG)
        else:
            tc=self.tc_Riazi_Daubert_ext()
        return unidades.Temperature(tc)

    def pc_Riazi_Adwani(self):
        if self.definicion==1:
            pc=prop_Riazi_Adwani(2, 0, self.Tb, self.SG)
        elif self.definicion==2:
            pc=prop_Riazi_Adwani(2, 1, self.M, self.SG)
        else:
            pc=self.pc_Riazi_Daubert_ext().bar
        return unidades.Pressure(pc, "bar")

    def vc_Riazi_Adwani(self):
        if self.definicion==1:
            vc=prop_Riazi_Adwani(3, 0, self.Tb, self.SG)
        elif self.definicion==2:
            vc=prop_Riazi_Adwani(3, 1, self.M, self.SG)
        else:
            vc=self.vc_Riazi_Daubert_ext().cm3g
        return unidades.SpecificVolume(vc*self.M, "cm3g")

    def I_Riazi_Adwani(self):
        if self.definicion==1:
            i=prop_Riazi_Adwani(4, 0, self.Tb, self.SG)
        elif self.definicion==2:
            i=prop_Riazi_Adwani(4, 1, self.M, self.SG)
        else:
            i=self.I_Riazi_Daubert_ext()
        return i

    def d20_Riazi_Adwani(self):
        if self.definicion==1:
            d20=prop_Riazi_Adwani(5, 0, self.Tb, self.SG)
        elif self.definicion==2:
            d20=prop_Riazi_Adwani(5, 1, self.M, self.SG)
        return unidades.Density(d20, "gcc")


    def tb_Soreide(self):
        """Soreide, I. “Improved Phase Behavior Predictions of Petroleum Reservoir Fluids from a Cubic Equation of State.” Doctor of Engineering dissertation, Norwegian Institute of Technology, Trondheim, 1989.
        """
        return unidades.Temperature(1928.3-1.695e5*self.SG**3.266/self.M**0.03522*exp(-4.922e-3*self.M-4.7685*self.SG+3.462*e-3*self.self.M*self.SG), "R")


    def tc_cavett(self):
        """Cavett, R. H. “Physical Data for Distillation Calculations—Vapor–Liquid Equilibrium.” Proceedings of the 27th Meeting, API, San Francisco, 1962, pp. 351–366."""
        return unidades.Temperature(10**(768.07121+1.7133693*self.Tb.F-0.0010834003*self.Tb.F**2-0.0089212579*self.API*self.Tb.F+0.38890584e-6*self.Tb.F**3+0.5309492e-5*self.API*self.Tb.F**2+0.327116e-7*self.API**2*self.Tb.F**2), "R")

    def pc_cavett(self):
        """Cavett, R. H. “Physical Data for Distillation Calculations—Vapor–Liquid Equilibrium.” Proceedings of the 27th Meeting, API, San Francisco, 1962, pp. 351–366."""
        return unidades.Pressure(10**(2.82904060+0.94120109e-3*self.Tb.F-0.30474749e-5*self.Tb.F**2-0.2087611e-4*self.API*self.Tb.F+0.15184103e-8*self.Tb.F**3+0.11047899e-7*self.API*self.Tb.F**2-0.48271599e-7*self.API**2*self.Tb.F+0.13949619e-9*self.API**2*self.Tb.F**2), "psi")


    def tc_Lee_Kesler(self):
        """Kesler, M. G., and B. I. Lee. “Improve Prediction of Enthalpy of Fractions.” Hydrocarbon Processing (March 1976): 153–158."""
        return unidades.Temperature(341.7+811.1*self.SG+(0.4244+0.1174*self.SG)*self.Tb.R+(0.4669-3.26238*self.SG)*1e5/self.Tb.R, "R")

    def pc_Lee_Kesler(self):
        """Kesler, M. G., and B. I. Lee. “Improve Prediction of Enthalpy of Fractions.” Hydrocarbon Processing (March 1976): 153–158."""
        return unidades.Pressure(exp(8.3634-0.0566/self.SG-(0.24244+2.2898/self.SG+0.11857/self.SG**2)*1e-3*self.Tb.R+(1.4685+3.648/self.SG+0.47227/self.SG**2)*1e-7*self.Tb.R**2-(0.42019+1.6977/self.SG**2)*1e-10*self.Tb.R**3), "psi")

    def peso_molecular_Lee_Kesler(self):
        """Kesler, M. G., and B. I. Lee. “Improve Prediction of Enthalpy of Fractions.” Hydrocarbon Processing (March 1976): 153–158."""
        return -12272.6+9486.4*self.SG+(4.6523-3.3287*self.SG)*self.Tb.R+(1-0.77084*self.SG-0.02058*self.SG**2)*(1.3437-720.79/self.Tb.R)*1e7/self.Tb.R+(1-0.80882*self.SG-0.02226*self.SG**2)*(1.8828-181.98/self.Tb.R)*1e12/self.Tb.R**3

    def factor_acentrico_Lee_Kesler(self):
        """Kesler, M. G., and B. I. Lee. “Improve Prediction of Enthalpy of Fractions.” Hydrocarbon Processing (March 1976): 153–158."""
        tita=self.Tb/self.Tc
        if tita>0.8:
            return -7.904+0.1352*self.watson-0.007465*self.watson**2+8359*tita+(1.408-0.01063*self.watson)/tita
        else:
            return (-log(self.Pc.atm)-5.92714+6.09648/tita+1.58862*log(tita)-0.169347*tita**6)/(15.2518-15.6875/tita-13.4721*log(tita)+0.43577*tita**6)


    def tc_Sim_Daubert(self):
        """Sim, W. J., and T. E. Daubert. “Prediction of Vapor-Liquid Equilibria of Undefined Mixtures.” Ind. Eng. Chem. Process Dis. Dev. 19, no. 3 (1980): 380–393."""
        return unidades.Temperature(exp(3.9934718*self.Tb.R**0.08615*self.SG**0.04614), "R")

    def pc_Sim_Daubert(self):
        """Sim, W. J., and T. E. Daubert. “Prediction of Vapor-Liquid Equilibria of Undefined Mixtures.” Ind. Eng. Chem. Process Dis. Dev. 19, no. 3 (1980): 380–393."""
        return unidades.Pressure(3.48242e9*self.Tb.R**-2.3177*self.SG**2.4853, "psi")

    def peso_molecular_Sim_Daubert(self):
        """Sim, W. J., and T. E. Daubert. “Prediction of Vapor-Liquid Equilibria of Undefined Mixtures.” Ind. Eng. Chem. Process Dis. Dev. 19, no. 3 (1980): 380–393."""
        return 1.4350476e-5*self.Tb.R**2.3776*self.SG**-0.9371


    def tc_Tsonopoulos(self):
        """Tsonopoulos, C., Heidman, J. L., and Hwang, S.-C.,Thermodynamic and Transport Properties of Coal Liquids, An Exxon Monograph, Wiley, New York, 1986."""
        return unidades.Temperature(10**(1.20016-0.61954*log10(self.Tb)+0.48262*log10(self.SG)+0.67365*log10(self.SG)**2))

    def pc_Tsonopoulos(self):
        """Tsonopoulos, C., Heidman, J. L., and Hwang, S.-C.,Thermodynamic and Transport Properties of Coal Liquids, An Exxon Monograph, Wiley, New York, 1986."""
        return unidades.Pressure(10**(7.37498-2.15833*log10(self.Tb)+3.35417*log10(self.SG)+5.64019*log10(self.SG)**2), "bar")


    def tc_Watansiri_Owens_Starling(self):
        """Watansiri, S., V. H. Owens, and K. E. Starling. “Correlations for Estimating Critical Constants, Acentric Factor, and Dipole Moment for Undefined Coal-Fluid Fractions.” Ind. Eng. Chem. Process Des. Dev. 24 (1985): 294–296."""
        return unidades.Temperature(exp(-0.0650504-0.0005217*self.Tb.R+0.03095*log(self.M)+1.11067*log(self.Tb.R)+self.M*(0.078154*self.SG**0.5-0.061061*self.SG**(1./3.)-0.016943*self.SG)), "R")

    def vc_Watansiri_Owens_Starling(self):
        """Watansiri, S., V. H. Owens, and K. E. Starling. “Correlations for Estimating Critical Constants, Acentric Factor, and Dipole Moment for Undefined Coal-Fluid Fractions.” Ind. Eng. Chem. Process Des. Dev. 24 (1985): 294–296."""
        return unidades.SpecificVolume(exp(76.313887-129.8038*self.SG+63.175*self.SG**2-13.175*self.SG**3+1.10108*log(self.M)+42.1958*log(self.SG))/self.M, "ft3lb")

    def pc_Watansiri_Owens_Starling(self):
        """Watansiri, S., V. H. Owens, and K. E. Starling. “Correlations for Estimating Critical Constants, Acentric Factor, and Dipole Moment for Undefined Coal-Fluid Fractions.” Ind. Eng. Chem. Process Des. Dev. 24 (1985): 294–296."""
        return unidades.Pressure(exp(6.6418853+0.01617283*(self.Tc.R/self.Vc.ft3lb)**0.8-8.712*self.M/self.Tc.R-0.08843889*self.Tb.R/self.M), "psi")

    def factor_acentrico_Watansiri_Owens_Starling(self):
        """Watansiri, S., V. H. Owens, and K. E. Starling. “Correlations for Estimating Critical Constants, Acentric Factor, and Dipole Moment for Undefined Coal-Fluid Fractions.” Ind. Eng. Chem. Process Des. Dev. 24 (1985): 294–296."""
        return 5*self.Tb.R/9/self.M*(5.12316667e-4*self.Tb.R+0.281826667*self.Tb.R/self.M+382.904/self.M+0.074691e-5*self.Tb.R**2/self.SG**2-0.120227778e-4*self.Tb.R*self.M+0.001261*self.SG*self.M+0.1265e-4*self.M**2+0.2016e-4*self.SG*self.M**2-66.29959*self.Tb.R**(1/3.)/self.M-0.00255452*self.Tb.R**(2./3)/self.SG**2)


    def tc_Edmister(self):
        """Edmister, W. C. “Applied Hydrocarbon Thermodynamics, Part 4, Compressibility Factors and Equations of State.” Petroleum Refiner 37 (April 1958): 173–179."""
        return unidades.Temperature(self.Tb.R*(3*log(self.Pc.atm)/7/(self.f_acent+1)+1), "R")

    def pc_Edmister(self):
        """Edmister, W. C. “Applied Hydrocarbon Thermodynamics, Part 4, Compressibility Factors and Equations of State.” Petroleum Refiner 37 (April 1958): 173–179."""
        return unidades.Pressure(10**(7/3.*(self.f_acent+1)*(self.Tc.R/self.Tb.R-1)), "atm")

    def tb_Edmister(self):
        """Edmister, W. C. “Applied Hydrocarbon Thermodynamics, Part 4, Compressibility Factors and Equations of State.” Petroleum Refiner 37 (April 1958): 173–179."""
        return unidades.Temperature(self.Tb.R/(3*log(self.Pc.atm)/7/(self.f_acent+1)+1), "R")

    def factor_acentrico_Edmister(self):
        """Edmister, W. C. “Applied Hydrocarbon Thermodynamics, Part 4, Compressibility Factors and Equations of State.” Petroleum Refiner 37 (April 1958): 173–179."""
        return 3*log10(self.Pc.atm)/7./(self.Tc.R/self.Tb.R-1)-1

    def factor_acentrico_Korsten(self):
        """Korsten, H., "Internally Consistent Prediction of Vapor Pressure and Related Properties," Industrial and Engineering Chemistry Research, 2000, Vol. 39, pp. 813-820."""
        tbr=self.Tb/self.Tc
        return 0.5899*tbr**1.3/(1-tbr**1.3)*log10(self.Pc.atm)-1.


    def tc_Magoulas(self):
        """Magoulas, S., and D. Tassios. Predictions of Phase Behavior of HT-HP Reservoir Fluids. Paper SPE no.37294. Richardson, TX: Society of Petroleum Engineers, 1990."""
        return unidades.Temperature(-1247.4+0.792*self.M+1971*self.SG-27000./self.M+707.4/self.SG, "R")

    def pc_Magoulas(self):
        """Magoulas, S., and D. Tassios. Predictions of Phase Behavior of HT-HP Reservoir Fluids. Paper SPE no.37294. Richardson, TX: Society of Petroleum Engineers, 1990."""
        return unidades.Pressure(exp(0.01901-0.0048442*self.M+0.13239*self.SG+227./self.M-1.1663/self.SG+1.2702*log(self.M)), "atm")

    def factor_acentrico_Magoulas(self):
        """Magoulas, S., and D. Tassios. Predictions of Phase Behavior of HT-HP Reservoir Fluids. Paper SPE no.37294. Richardson, TX: Society of Petroleum Engineers, 1990."""
        return -0.64235+0.00014667*self.M+0.021876*self.SG-4.539/self.M+0.21699*log(self.M)

    def tc_Twu(self):
        t, p, v, g=Paraffin_Twu(self.Tb.R)
        f=(exp(5*(g-self.SG))-1)*(-0.362456/self.Tb.R**0.5+(0.0398285-0.948125/self.Tb.R**0.5)*(exp(5*(g-self.SG))-1))
        tc=t*((1+2*f)/(1-2*f))**2
        return unidades.Temperature(tc, "R")

    def pc_Twu(self):
        t, p, v, g=Paraffin_Twu(self.Tb.R)
        f=(exp(0.5*(g-self.SG))-1)*((2.53262-46.19553/self.Tb.R**0.5-0.00127885*self.Tb.R)+(-11.4277+252.14/self.Tb.R**0.5+0.00230535*self.Tb.R)*(exp(0.5*(g-self.SG))-1))
        pc=p*self.Tc/t*v/self.Vc*((1+2*f)/(1-2*f))**2
        return unidades.Pressure(pc, "psi")

    def vc_Twu(self):
        t, p, v, g=Paraffin_Twu(self.Tb.R)
        f=(exp(4*(g**2-self.SG**2))-1)*(0.46659/self.Tb.R**0.5+(-0.182421+3.01721/self.Tb.R**0.5)*(exp(4*(g**2-self.SG**2))-1))
        vc=v*((1+2*f)/(1-2*f))**2
        return unidades.SpecificVolume(vc/self.M, "ft3lb")

    def peso_molecular_Twu(self):
        t, p, v, g=Paraffin_Twu(self.Tb.R)
        mo=self.Tb/(5.8-0.0052*self.Tb)
        X=abs(0.012342-0.244541/self.Tb**0.5)
        DG=exp(5*(g-self.SG))-1.
        f=DG*(X+(-0.0175691+0.143979/self.Tb**0.5)*DG)
        return exp(log(mo)*((1+2*f)/(1-2*f))**2)


    def vc_Hall_Yarborough(self):
        """Hall, K. R., and L. Yarborough. “New Simple Correlation for Predicting Critical Volume.” Chemical Engineering (November 1971): 76."""
        return unidades.SpecificVolume(0.025*self.M**0.15/self.SG**0.7985, "ft3lb")

    def vc_Riedel(self):
        """Estimación del volumen crítico haciendo uso del método de Riedel. API procedure 4A3.1 pag. 302"""
        riedel=5.811+4.919*self.f_acent
        return unidades.SpecificVolume(R_atml*self.Tc/self.Pc.atm/(3.72+0.26*(riedel-7.))/self.M, "lg")

    def peso_molecular_API(self):
        return 20.486*exp(1.165e-4*self.Tb.R-7.78712*self.SG+1.1582e-3*self.Tb.R*self.SG)*self.Tb.R**1.26007*self.SG**4.98308
        return 20.486*self.SG**(1.565*self.SG)

    def peso_molecular_pesado(self):
        """Método alternativo para calcular el peso molecular en el caso de fracciones pesadas de petroleo, API procedure 2B2.3 pag.219
        Riazi, M. R., Daubert, 1". E., "Molecular Weight of Heavy Fractions from Viscosity," Oil and Gas Journal, Vol. 58, No. 52, 1987, pp. 110-113."""
        return 223.56*self.v100**(-1.2435+1.1228*self.SG)*self.v210**(3.4758-3.038*self.SG)*self.SG**-0.6665

    def peso_molecular_ASTM(self):
        """Método histórico, ASTM D 2502
        ASTM, Annual Book of Standards, Section Five: Petroleum Products, Lubricants, and Fossil Fuels (in Five Volumes), ASTM International, West Conshohocken, PA, 2002"""
        H100=870.*log10(log10(self.v100.cSt+0.6))+154.
        H210=870.*log10(log10(self.v210.cSt+0.6))+154.
        VSF=H100-H210
        K=4.145-1.733*log10(VSF-145)
        return 180+K*(H100+60.)

    def peso_molecular_Goossens(self):
        """Goossens, A. G., "Prediction of Molecular Weight of Petroleum Fractions," Industrial and Engineering Chemistry Research, Vol. 35, 1996, pp. 985 988."""
        b=1.52869+0.06486*log(self.Tb/(1078-self.Tb))
        return 0.01077*self.Tb**b/self.d20.gcc


    def Zc(self):
        return self.Pc.atm*self.Vc*self.M/R_atml/self.Tc

    def Zc_Lee_Kesler(self):
        """Lee, B. I. and Kesler, M. G., "A Generalized Thermodynamic Correlation Based on Three- Parameter Corresponding States," American Institute of Chemical Engineers Journal, Vot. 21, 1975,"""
        return 0.2905-0.085*self.f_acent

    def Zc_Haugen(self):
        """Haugen, O. A., K. M. Watson, and R. A. Ragatz. Chemical Process Principles, 2nd ed. New York: Wiley, 1959, p. 577."""
        return 1./(1.28*self.f_acent+3.41)

    def Zc_Reid(self):
        """Reid, R., J. M. Prausnitz, and T. Sherwood. The Properties of Gases and Liquids, 3rd ed. New York: McGraw-Hill, 1977, p. 21."""
        return 0.291-0.08*self.f_acent

    def Zc_Salerno(self):
        """Salerno, S., et al. “Prediction of Vapor Pressures and Saturated Volumes.” Fluid Phase Equilibria 27 (June 10, 1985): 15–34."""
        return 0.291-0.08*self.f_acent-0.016*self.f_acent**2

    def Zc_Nath(self):
        """Nath, J. “Acentric Factor and the Critical Volumes for Normal Fluids.” Industrial Engineering and Chemical. Fundamentals 21, no. 3 (1985): 325–326."""
        return 0.2918-0.0928*self.f_acent





    def nT(self, T):
        """Indice de refracción a una temperatura diferente de los 20º"""
        return self.n-0.0004*(T-293.15)


    def PNA_Riazi(self):
        """Calculo de la composición en parafinas, naftanos y compuestos aromáticos de la fracción petrolífera, API procedure 2B4.1 pag 225"""
        if self.has_v100:
            if self.M<200:
                a, b, c, d, e, f=(-13.359, 14.4591, -1.41344, 23.9825, -23.333, 0.81517)
            else:
                a, b, c, d, e, f=(2.5737, 1.0133, -3.573, 2.464, -3.6701, 1.96312)
            xp=a+b*self.Ri+c*self.VGC
            xn=d+e*self.Ri+f*self.VGC
        elif self.has_CH:
            if self.M<200:
                xp=2.57-2.877*self.SG+0.02876*self.CH
                xn=0.52641-0.7494*xp-0.021811*self.m
            else:
                xp=1.9842-0.27722*self.Ri-0.15643*self.CH
                xn=0.5977-0.761745*self.Ri+0.068048*self.CH
        else:
            if self.M<200:
                xp=3.7387-4.0829*self.SG+0.014772*self.m
                xn=-1.5027+2.10152*self.SG-0.02388*self.m
            else:
                xp=1.9382+0.074855*self.m-0.19966*self.CH
                xn=-0.4226-0.00777*self.m+0.107625*self.CH
        xa=1-xp-xn
        return xp, xn, xa

    def PNA_Peng_Robinson(self):
        """Robinson, D. B., and D. Y. Peng. “The Characterization of the Heptanes and Heavier Fractions.” Research Report 28. Tulsa, OK: GPA, 1978."""
        Tbp=exp(log(1.8)+5.8345183+0.84909035e-1*(self.Nc-6)-0.52635428e-2*(self.Nc-6)**2+0.21252908e-3*(self.Nc-6)**3-0.44933363e-5*(self.Nc-6)**4+0.37285365e-7*(self.Nc-6)**5)
        Tbn=exp(log(1.8)+5.8579332+0.79805995e-1*(self.Nc-6)-0.43098101e-2*(self.Nc-6)**2+0.14783123e-3*(self.Nc-6)**3-0.27095216e-5*(self.Nc-6)**4+0.19907794e-7*(self.Nc-6)**5)
        Tba=exp(log(1.8)+5.867176+0.80436947e-1*(self.Nc-6)-0.47136506e-2*(self.Nc-6)**2+0.18233365e-3*(self.Nc-6)**3-0.38327239e-5*(self.Nc-6)**4+0.32550576e-7*(self.Nc-6)**5)
        Mp=14.026*self.Nc+2.016
        Mn=14.026*self.Nc-14.026
        Ma=14.026*self.Nc-20.074
        a=[[1, 1, 1], [Tbp*Mp, Tbn*Mn, Tba*Ma], [Mp, Mn, Ma]]
        b=[1, self.M*self.WABP.R, self.M]
        Xp, Xn, Xa=solve(a, b)
        pcp=(206.126096*self.Nc+29.67136)/(0.227*self.Nc+0.34)**2
        pcn=(206.126096*self.Nc+206.126096)/(0.227*self.Nc+0.137)**2
        pca=(206.126096*self.Nc+295.007504)/(0.227*self.Nc+0.325)**2
        pc=Xp*pcp+Xn*pcn+Xa*pca
        wp=0.432*self.Nc-0.0457
        wn=0.0432*self.Nc-0.088
        wa=0.0445*self.Nc-0.0995
        S=0.996704+0.00043155*self.Nc
        S1=0.99627245+0.00043155*self.Nc
        tcp=S*(1+(3*log10(pcp)-3.501952)/7/(1+wp))*Tbp
        tcn=S1*(1+(3*log10(pcn)-3.501952)/7/(1+wn))*Tbn
        tca=S1*(1+(3*log10(pca)-3.501952)/7/(1+wa))*Tba
        tc=Xp*tcp+Xn*tcn+Xa*tca
        return Xp, Xn, Xa

    def PNA_Bergman(self):
        """Bergman, D. F., M. R. Tek, and D. L. Katz. “Retrograde Condensation in Natural Gas Pipelines.” Project PR 2-29 of Pipelines Research Committee, AGA, January 1977."""
        Xwa=8.47-self.watson
        gp=0.582486+0.00069481*(self.Tb.R-460)-0.7572818e-6*(self.Tb.R-460)**2+0.3207736e-9*(self.Tb.R-460)**3
        gn=0.694208+0.0004909267*(self.Tb.R-460)-0.659746e-6*(self.Tb.R-460)**2+0.330966e-9*(self.Tb.R-460)**3
        ga=0.916103-0.000250418*(self.Tb.R-460)+0.357967e-6*(self.Tb.R-460)**2-0.166318e-9*(self.Tb.R-460)**3
        a=[[1, 1], [1/gp, 1/gn]]
        b=[1-Xwa, 1/self.SG-Xwa/ga]
        Xwp, Xwn=solve(a, b)
        Tcp=275.23+1.2061*(self.Tb.R-460)-0.00032984*(self.Tb.R-460)**2
        Pcp=573.011-1.13707*(self.Tb.R-460)+0.00131625*(self.Tb.R-460)**2-0.85103e-6*(self.Tb.R-460)**3
        wp=0.14+0.0009*(self.Tb.R-46)+0.233e-6*(self.Tb.R-460)**2
        Tcn=156.8906+2.6077*(self.Tb.R-460)-0.003801*(self.Tb.R-460)**2+0.2544e-5*(self.Tb.R-460)**3
        Pcn=726.414-1.3275*(self.Tb.R-460)+0.9846e-3*(self.Tb.R-460)**2-0.45169e-6*(self.Tb.R-460)**3
        wn=wp-0.075
        Tca=289.535+1.7017*(self.Tb.R-460)-0.0015843*(self.Tb.R-460)**2+0.82358e-6*(self.Tb.R-460)**3
        Pca=1184.514-3.44681*(self.Tb.R-460)+0.0045312*(self.Tb.R-460)**2-0.23416e-5*(self.Tb.R-460)**3
        wa=wp-0.1
        pc=Xwp*Pcp+Xwn*Pcn+Xwa*Pca
        tc=Xwp*Tcp+Xwn*Tcn+Xwa*Tca
        w=Xwp*wp+Xwn*wn+Xwa*wa
        return Xwp, Xwn, Xwa

    def PNA_van_Nes(self):
        """Van Nes, K. and Van Western, H. A., Aspects of the Constitution of Mineral Oils, Elsevier, New York, 1951."""
        v=2.51*(n-1.475)-(self.d20-0.851)
        w=(self.d20-0.851)-1.11*(self.n-1.475)
        if v>0: a=430
        else: a=670
        if w>0:
            Cr=820*w-3*self.S+10000./self.M
            Rt=1.33+0.146*self.M*(w-0.005*self.S)
        else:
            Cr=1440*w-3*self.S+10600./self.M
            Rt=1.33+0.18*self.M*(w-0.005*self.S)
        Ra=0.44+b*v*self.M
        Rn=Rt-Ra
        Ca=a*v+3660/self.M
        Cn=Cr-Ca
        Cp=100-Cr
        return Cp/100., Cn/100., Ca/100.

    def H_Riazi(self):
        return (100.-self.S)/(1.+self.CH)

    def H_Goossens(self):
        return 30.346+(82.952-65.341*self.n)/self.d20-306./self.M

    def H_ASTM(self):
        """ASTM, Annual Book of Standards, ASTM International, West Conshohocken, PA, 2002."""
        Tb=(self.T10+self.T50+self.T90)/3.
        return (5.2407+0.01448*Tb-7.018*self.xa)/self.SG-0.901*self.xa+0.01298*self.xa*Tb-0.01345*Tb-0.01345*Tb+5.6879

    def H_Jenkins_Walsh(self):
        """Jenkins, G. I. and Walsh, R. E, "Quick Measure of Jet Fuel Properties," Hydrocarbon Processing, Vol. 47, No. 5, 1968, pp. 161-164."""
        return 11.17-12.89*self.SG+0.0389*self.AP

    def S_Riazi(self):
        """Riazi, M. R., Nasimi, N., and Roomi, Y., "Estimating Sulfur Content of Petroleum Products and Crude Oils," Industrial and Engineering Chemistry Research, Vol. 38, No. 11, 1999, pp. 4507-4512."""
        if self.M<200:
            S=177.448-170.946*self.Ri+0.2258*self.m+4.054*self.SG
        else:
            S=-58.02+38.463*self.Ri-0.023*self.m+22.4*self.SG
        return S


    def Reid_Blend(self):
        """Método de cálculo de la presión de vapor Reid, API procedure 5B1.3 pag 407"""
        suma=0
        for i in range(len(self.componente)):
            Reid=Pressure(self.componente[i].Pv_DIPPR(Temperature(100, "F")))
            suma+=self.fraccion[i]*Reid.psi**1.2
        return suma**(1/1.2)


    def Reid_simulate(self):
        """Método de cálculo de la presión de vapor Reid simulando un proceso de destilación ASTM D323-94, API procedure 5B1.4"""
        pass


    def tc(self):
        """Método de cálculo de la temperatura crítica de fracciones de petroleo, API procedure 4D1.1, pag 331"""
        delta=self.SG*(self.VABP.F+100)
        return unidades.Temperature(186.16+1.6667*delta-0.7127e-3*delta**2, "F")

    def pc(self):
        """Método de cálculo de la presión cricia de fracciones de petroleo, API procedure 4D2.1 pag 334
        Dificil de implementar ya que se trata de un método gráfico sin ecuaciones asociadas"""
        pass

    def Pv_Maxwell_Bonnell(self, T, mod=False):
        """Maxwell, J. B. and Bonnell, L. S., Vapor Pressure Charts for Petroleum Engineers, Exxon Research and Engineering Company, Florham Park, NJ, 1955. Reprinted in 1974."Deviation and Precision of a New Vapor Pressure Correlation for Petroleum Hydrocarbons," Industrial and Engineering Chemistry, Vol. 49, 1957, pp. 1187-1196.
        modificación de Tsonopoulos para coal liquids:  Tsonopoulos, C., Heidman, J. L., and Hwang, S.-C.,Thermodynamic and Transport Properties of Coal Liquids, An Exxon Monograph, Wiley, New York, 1986."""
        if mod:
            if self.Tb<=366.5:
                F1=0
            else:
                F1=-1+0.009*(self.Tb-255.37)
            F2=(self.watson-12)-0.01304*(self.watson-12)**2
        else:
            if  self.Tb<367:
                F=0
            elif self.Tb<478:
                F=-3.2985+0.0009*self.Tb
            else:
                F=-3.2985+0.009*self.Tb

        pvap=0.
        pvapcalc=1.
        while abs(pvap-pvapcalc)<1e-6:
            pvap=pvapcalc
            if mod:
                if pvap<=760:
                    F3=1.47422*log10(pvap/760.)
                else:
                    F3=1.47422*log10(pvap/760.)+1.190833*(log10(pvap/760.))**2
                DTb=F1*F2*F3
            else:
                DTb=1.3889*F*(self.watson-12)*log10(pvap/760.)
            Tb_=self.Tb-DTb
            Q=(Tb_/T-0.00051606*Tb_)/(748.1-0.3861*Tb_)
            if Q>0.0022:
                pvapcalc=10**((3000.538*Q-6.76156)/(43*Q-0.987672))
            elif Q>=0.0013:
                pvapcalc=10**((2663.129*Q-5.994296)/(95.76*Q-0.972546))
            else:
                pvapcalc=10**((2770.085*Q-6.412631)/(36*Q-0.989679))

        return unidades.Pressure(pvapcalc, "mmHg")


    def Pv_Tsonopoulos(self, T):
        """Tsonopoulos, C., Heidman, J. L., and Hwang, S.-C.,Thermodynamic and Transport Properties of Coal Liquids, An Exxon Monograph, Wiley, New York, 1986."""
        Tr=T/self.tpc
        A=5.671485+12.439604*self.f_acent
        B=5.809839+12.755971*self.f_acent
        C=0.867513+9.654169*self.f_acent
        D=0.1383536+0.316367*self.f_acent
        pr=exp(A-B/Tr-C*log(Tr)+D*Tr**6)
        return unidades.Pressure(pr, "bar")

    def Pv_simple(self, T):
        """eq 7.25"""
        pv=10**(3.2041*(1.-0.998*(self.Tb-41)/(self.Tb-41)*(1393-T)/(1393-self.Tb)))
        return unidades.Pressure(pr, "bar")

    def Pv_gasoline(self, T):
        """Método de cálculo de la presión de vapor de productos terminados de petroleo, gasolinas, naftas, etc, API procedure 5B1.1 pag 404"""
        #FIXME: No da un resultado correcto, posiblemente algún parénteis mal puesto
        t=unidades.Temperature(T)
        SL=(unidades.Temperature(self.D86[-3]).F-unidades.Temperature(self.D86[1]).F)/15.
        p=exp(1/t.R*(21.36412862-6.7769666*sqrt(SL)-0.93213944*log(self.ReidVP.psi)+1.42680425*sqrt(SL)*log(self.ReidVP.psi)-0.29458386*self.ReidVP.psi+(-0.00568374+0.00577103*sqrt(SL)-0.00106045*sqrt(SL)*log(self.ReidVP.psi)+0.00060246*self.ReidVP.psi)*t.R+(-10177.78660360+2306.00561642*sqrt(SL)+1097.68947465*log(self.ReidVP.psi)-463.19014182*sqrt(SL)*log(self.ReidVP.psi)+65.61239475*self.ReidVP.psi+0.13751932*self.ReidVP.psi**2)))
        return unidades.Pressure(p, "psi")

    def Pv_API(self, T):
        """Método de cálculo de la presión de vapor de fracciones de petroleo, API procedure 5B1.2 pag 406"""
        t=unidades.Temperature(T)
        p=exp(7.78511307-1.08100387*log(self.ReidVP.psi)+0.05319502*self.ReidVP.psi+0.00451316*t.R+(-5756.85623050+1104.41248797*log(self.ReidVP.psi)-0.00068023*self.ReidVP.psi**4)/t.R)
        return unidades.Pressure(p, "psi")


    def RhoL(self, T):
        """Método de cálculo de la densidad del líquido de fraciones de petroleo a presión atmosférica, API procedure 6A3.5,pag 494"""
        t=unidades.Temperature(T)
        rho=62.3636*sqrt(self.SG**2-(1.2655*self.SG-0.5098+8.011e-5*self.Tb.R)*(t.R-519.67)/self.Tb.R)
        return unidades.Density(rho, "lbft3")

    def RhoL_Rackett(self, T):
        """Método de cálculo de la densidad de la fase líquida de fracciones de petróleo a presión atmosférica, API procedure 6A3.6 pag 495"""
        rho=self.SG*1000
        t=unidades.Temperature(60, "F")
        Zra=(self.Pc.atm/(rho/self.M*R_atml*self.Tc))**(1/(1+(1-t/self.Tc)**(2./7)))
        inv=R_atml*self.Tc/self.Pc.atm*Zra**(1+(1-T/self.Tc)**(2./7))
        return unidades.Density(1/inv*self.M, "gl")

    def RhoL_Presion(self, T, P):
        """Método de cálculo de la densidad de la fase líquida de fracciones de petróleo a alta presión, API procedure 6A3.7, 6A3.10 pag 497"""
        rho0=self.RhoL_Rackett(T)
        p=unidades.Pressure(P, "atm")
        t=unidades.Temperature(T)
        B20=10**(-6.1e-4*t.F+4.9547+0.7133*rho0.kgl)
        m=21646+0.0734*p.psig+1.4463e-7*p.psig**2
        X=(B20-100000)/23170
        B1=1.52e+4+4.704*p.psig-2.5807e-5*p.psig**2+1.0611e-10*p.psig**3
        Bt=m*X+B1
        return unidades.Density(rho0/(1.0-p.psig/Bt))


    def Reduccion_volumetrica(self, C, G):
        """Cálculo de la reducción volumétrica que se produce al mezclar hidrocarburos de bajo peso molecular con otra fracción pesada de crudo, API procedure 6A3.11, pag 504
        requiere dos parámetros:
        C: fracción volumétria del componente ligero en la mezclar
        G: Diferencia de gravedad en grados API"""
        return 2.14e-3*C**-0.0704*G**1.76


    def Entalpia(self, T, P):
        """Cálculo de la entalpia de fracciones petrolíferas, API procedure 7B4.7, pag 658"""
        T=unidades.Temperature(T)
        if T<self.Tb and self.tr(T)<=0.8 and self.pr(P)<1.:
            A1=1e-3*(-1171.26+(23.722+24.907*self.SG)*self.watson+(1149.82-46.535*self.watson)/self.SG)
            A2=1e-6*((1.+0.82463*self.watson)*(56.086-13.817/self.SG))
            A3=-1e-9*((1.+0.82463*self.watson)*(9.6757-2.3653/self.SG))
            H=A1*(T.R-259.7)+A2*(T.R**2-259.7**2)+A3*(T.R**3-259.7**3)
        else:
            Hl=self.Entalpia(0.8*self.Tc, 0.1*self.Pc.atm)
            if T<self.Tb:
                H_adimensional_presion=Lee_Kesler_Entalpia_lib(self.tr(T), self.pr(P), self.f_acent, fase=1)
            else:
                H_adimensional_presion=Lee_Kesler_Entalpia_lib(self.tr(T), self.pr(P), self.f_acent, fase=0)
            if 10.<self.watson<12.8 and 0.7<self.SG<0.885:
                B4=((12.8/self.watson-1.)*(1.-10./self.watson)*(self.SG-0.885)*(self.SG-0.7)*1e4)**2
            else:
                B4=0
            B1=1e-3*(-356.44+29.72*self.watson+B4*(295.02-248.49/self.SG))
            B2=1e-6*(-146.24+(77.62-2.772*self.watson)*self.watson-B4*(301.42-253.87/self.SG))
            B3=1e-9*(-56.487-2.95*B4)
            H=Hl.Btulb+B1*(T.R-0.8*self.Tc.R)+B2*(T.R**2-0.64*self.Tc.R**2)+B3*(T.R**3-0.512*self.Tc.R**3)+R_Btu*self.Tc.R/self.M*(4.507+5.266*self.f_acent-H_adimensional_presion)
        return unidades.Enthalpy(H, "Btulb")


    def Cp_liquid_API(self, T):
        """Cálculo de la capacidad calorífica isobárica de la fase liquida de fracciones petrolíferas, API procedure 7D2.2, pag 702"""
        t=unidades.Temperature(T)
        A1=-1.17126+(0.023722+0.024907*self.SG)*self.watson+(1.14982-0.046535*self.watson)/self.SG
        A2=1e-4*(1+0.82463*self.watson)*(1.12172-0.27634/self.SG)
        A3=-1e-8*(1+0.82463*self.watson)*(2.9027-0.70958/self.SG)
        return unidades.SpecificHeat(A1+A2*t.R+A3*t.R**2, "BtulbF")

    def Cp_liquid_Tsonopoulos(self, T):
        """Tsonopoulos, C., Heidman, J. L., and Hwang, S.-C.,Thermodynamic and Transport Properties of Coal Liquids, An Exxon Monograph, Wiley, New York, 1986."""
        cp=(0.28299+0.23605*self.watson)*(0.645-0.05959*self.SG+(2.32056-0.94752*self.SG)*(T/1000.-0.25537))
        return unidades.SpecificHeat(cp, "kJkgK")

    def Cp_liquid_Kesler_Lee(self, T):
        """Kesler, M. G. and Lee, B. I., "Improve Prediction of Enthalpy of Fractions," Hydrocarbon Processing, Vol. 55, No. 3, 1976, pp. 153-158."""
        a=1.4651+0.2302*self.watson
        b=0.306469-0.16734*self.SG
        c=0.001467-0.000551*self.SG
        return unidades.SpecificHeat(a*(b+c*T), "kJkgK")

    def Cp_liquid_Bondi(self, T):
        """Ref Eq 7.40 Riazi - Characterization of petroleum fraction, pag 333"""
        Tr=T/self.Tc
        Cp_ideal
        cp_adimensional=1.586+0.49*(1-Tr)+self.f_acent*(4.2775+6.3*(1-Tr)**(1./3)/Tr+0.4355/(1-Tr))
        return unidades.SpecificHeat(cp_adimensional*R+Cp_ideal, "kJkgK")


    def Lee_Kesler_lib_Cp(self, T, P):
        """Librería para el cálculo de capacidades calorificas, usada a continuación en diferentes funciones
        Procedure API 7E1.6 Pag.726"""
        #FIXME: No concuerdan mucho los valores de cp y cv con los valores por DIPPR
        Tr=self.tr(T)
        vr0, vrh=self.Lee_Kesler_lib(T, P.atm, "gas")
        E=0.042724/2/Tr**3/0.060167*(0.65392+1-(0.65392+1+0.060167/vr0**2)*exp(-0.060167/vr0**2))
        Cv0=-2*(0.154790+3*0.030323/Tr)/Tr**2/vr0+6*E

        E=0.041577/2/Tr**3/0.03754*(1.226+1-(1.226+1+0.03754/vrh**2)*exp(-0.03754/vrh**2))
        Cvh=-2*(0.027655+3*0.203488/Tr)/Tr**2/vrh+3*0.016901/Tr**3/vrh**2+6*E

        return Cv0, Cvh

    def Cp_gas(self, T, P):
        """Cálculo de la capacidad calorífica isobárica de la fase vapor de fracciones petrolíferas, API procedure 7D4.2, pag 717"""
        if 10.<=self.watson<=12.8 and 0.70<self.SG<=0.885:
            A4=((12.8/self.watson-1.)*(1.-10/self.watson)*(self.SG-0.885)*(self.SG-0.7)*1e4)**2
        else: A4=0
        A1=-0.35644+0.02972*self.watson+A4*(0.29502-0.24846/self.SG)
        A2=-1e-4*(2.9247-(1.5524-0.05543*self.watson)*self.watson+A4*(6.0283-5.0694/self.SG))
        A3=-1e-7*(1.6946+0.0844*A4)
        #TODO: calculo factor cp adimensional
        Cp0, Cph=self.Lee_Kesler_lib_Cp(T, P)
#        Cp_adimensional=-1.107
        Cp_adimensional=Cp0+self.f_acent/0.3978*(Cph-Cp0)
        return unidades.SpecificHeat(A1+A2*t.R+A3*t.R**2-R_Btu/self.M*Cp_adimensional, "BtulbF")


    def flash(self):
        """Método de cálculo del equilibrio líquido-vapor, API procedure 8D1.5, pag 828"""
        pass

    def flash2(self):
        """Método de cálculo del equilibrio líquido-vapor, cuando parte de la composición sí es conocida, API procedure 8D1.6, pag 832"""
        pass


    def Tension_API(self,T):
        """Método de cálculo de la tensión superficial de fracciones petrolíferas, API procedure 10A3.2, pag 997"""
        t=unidades.Temperature(T)
        sigma=673.7/self.watson*((self.Tc.R-t.R)/self.Tc.R)**1.232
        return unidades.Tension(sigma, "dyncm")

    def Tension_Baker_Swerdloff(self, T, P=1):
        """Baker, O. and Swerdloff, W.: "Finding Surface Tension of Hydrocarbon Liquids," Oil and Gas J. (Jan. 2, 1956) 125"""
        t=unidades.Temperature(T)
        p=unidades.Pressure(P, "atm")
        s68=39.-0.2571*self.API
        s100=37.5-0.2571*self.API
        if t.F<=68.:
            s=s68
        elif t.F<100:
            s=s68-((t.F-68)*(s68-s100))/32
        else:
            s=s100
        F=1.-0.024*p.psi**0.45
        return unidades.Tension(F*s, "dyncm")

    def Tension_Tsonopoulos(self, T):
        """Método alternativo de cálculo de la tensión superficial
        Tsonopoulos, C., Heidman, J. L., and Hwang, S.-C.,Thermodynamic and Transport Properties of Coal Liquids, An Exxon Monograph, Wiley, New York, 1986."""
        Pa=1.7237*self.Tb**0.05873*self.SG**-0.64927
        rhoL=self.RhoL(T)
        rhoV
        return unidades.Tension((Pa*(rhoL-rhoV))**4, "dyncm")

    def Tension_Miqueu(self, T):
        """Método alternativo de cálculo de la tensión superficial
        Miqueu, C., Satherley, J., Mendiboure, B., Lachiase, J., and Graciaa, A., "The Effect of P/N/A Distribution on the Parachors of Petroleum Fractions," Fluid Phase Equilibria, Vol. 180, 2001,pp. 327-344."""
        Pa=(0.85-0.19*self.f_acent)*self.Tc**(12/11.)/(self.Pc.bar/10)**(9./11)
        rhoL=self.RhoL(T)
        rhoV
        return unidades.Tension((Pa/self.M*(rhoL-rhoV))**(11./3), "dyncm")

    def Tension_PNA(self, T):
        """Método alternativo de cálculo de la tensión superficial si se conoce la composición PNA
        Miqueu, C., Satherley, J., Mendiboure, B., Lachiase, J., and Graciaa, A., "The Effect of P/N/A Distribution on the Parachors of Petroleum Fractions," Fluid Phase Equilibria, Vol. 180, 2001,pp. 327-344."""
        PaP=27.503+2.9963*self.M
        PaN=18.384+2.7367*self.M
        PaA=25.511+2.8332*self.M
        Pa=self.xa*PaA+self.xn*PaN+self.xp*PaP
        rhoL=self.RhoL(T)
        rhoV
        return unidades.Tension((Pa/self.M*(rhoL-rhoV))**(11./3), "dyncm")


    def Mu_Singh(self, T, mu):
        """Calculo de la viscosidad cinemática de fracciones petrolíferas a baja presión, conocida la viscosidad cinemática a 100ºF, API procedure 11A4.1 pag 1054
        mu: viscosidad cinematica experimental a 100ºF
        valor obtenido en centistokes"""
        t=unidades.Temperature(T)
        S=0.28008*log10(mu)+1.8616
        B=log10(mu)+0.86960
        return 10**(B*(559.67/t.R)**S-0.86960)


    def V100_API(self):
        """Cálculo de la viscosidad cinemática a 100ºF de fracciones petrolíferas a baja presión, API procedure 11A4.2, pag 1056"""
        A1=34.9310-8.84387e-2*self.Tb.R+6.73513e-5*self.Tb.R**2-1.01394e-8*self.Tb.R**3
        A2=-2.92649+6.98405e-3*self.Tb.R-5.09947e-6*self.Tb.R**2+7.49378e-10*self.Tb.R**3
        mu_ref=10**(-1.35579+8.16059e-4*self.Tb.R+8.38505e-7*self.Tb.R**2)
        mu_cor=10**(A1+A2*self.watson)
        return unidades.Diffusivity(mu_ref+mu_cor, "cSt")

    def V210_API(self):
        """Cálculo de la viscosidad cinemática a 210ºF de fracciones petrolíferas a baja presión, API procedure 11A4.2, pag 1056"""
        if self.v100:
            v100=self.v100
        else:
            v100=self.V100_API()
        return unidades.Diffusivity(10**(-1.92353+2.41071e-4*self.Tb.R+0.5113*log10(self.Tb.R*v100)), "cSt")

    def Viscosidad_ASTM(self, T, T1=unidades.Temperature(100, "F"), T2=unidades.Temperature(210, "F"), mu1=0, mu2=0):
        """Cálculo de la viscosidad cinemática a cualquier temperatura, conociendo otros dos valores de viscosidad a otras temperaturas, API procedure 11A4.4, pag 1063
        Parámetros:
        T:Temperatura a la que se quiere calcular la viscosidad
        T1,T2:opcional, temperatura a la que se conoce la viscosidad
        mu1,mu2:opcionales, valores de viscosidad conocidos
        Si no se suministran los parámetros opcionales se consideran los valores a 100 y 210ºF
        """
        if mu1==0:
            mu1=self.v100.cSt
        if mu2==0:
            mu2=self.v210.cSt
        t=unidades.Temperature(T)
        Z1=mu1+0.7+exp(-1.47-1.84*mu1-0.51*mu1**2)
        Z2=mu2+0.7+exp(-1.47-1.84*mu2-0.51*mu2**2)
        B=(log10(log10(Z1))-log10(log10(Z2)))/(log10(T1.R)-log10(T2.R))
        Z=10**(10**(log10(log10(Z1))+B*(log10(t.R)-log10(T1.R))))
        return unidades.Diffusivity(Z-0.7-exp(-0.7487-3.295*(Z-0.7)+0.6119*(Z-0.7)**2-0.3191*(Z-0.7)**3), "cSt")

    def Viscosidad_liquido_blend(self, T, fraccion_masica, petro1, petro2):
        """Método de cálculo de la viscosidad de líquidos en mezclas de fracciones petrolíferas, API procedure 11A4.5, pag 1066
        Los parámetros petro tienen la estructura [T1,T2,mu1,mu2]"""
        #TODO: de momoento el procedimiento requiere como parámetros petro1 y petro2, matrices con cuatro elementos, dos temperaturas y sus correspondientes viscosidades, cuando se defina correctamente las fracciones petroliferas estos parámetros serán sustituidos por un simple id de fracción petrolífera
        t=unidades.Temperature(T)
        T1=unidades.Temperature(petro1[0])
        T2=unidades.Temperature(petro1[1])

        ml=(log(log(petro1[3]+0.7))-log(log(petro1[2]+0.7)))/(log(T2.R)-log(T1.R))
        bl=log(log(petro1[2]+0.7))-ml*log(T1.R)
        mh=(log(log(petro2[3]+0.7))-log(log(petro2[2]+0.7)))/(log(T2.R)-log(T1.R))
        bh=log(log(petro2[2]+0.7))-mh*log(T1.R)

        Tl=exp((log(log(petro2[2]+0.7))-bl)/ml)
        Tx=exp(fraccion_masica[0]*log(Tl)+fraccion_masica[1]*log(T1.R))
        Th=exp((log(log(petro1[3]+0.7))-bh)/mh)
        Ty=exp(fraccion_masica[0]*log(T2.R)+fraccion_masica[1]*log(Th))

        m=(log(log(petro1[3]+0.7))-log(log(petro2[2]+0.7)))/(log(Ty)-log(Tx))
        b=log(log(petro2[2]+0.7))-m*log(Tx)

        return exp(exp(m*log(t.R)+b))-0.7


    def Viscosity_Index(self):
        """Método de calculo del indice de viscosidad, API procedure 11A6.1, pag 1083"""
        if self.v210.cSt>70:
            L=0.8353*self.v210.cSt**2+14.67*self.v210.cSt-216
            H=0.1684*self.v210.cSt**2+11.85*self.v210.cSt-97
        else: #Ajuste de los datos de la tabla, producción propia
            """Polynomial Fit of dataset: Table1_L, using function: a0+a1*x+a2*x^2+a3*x^3+a4*x^4+a5*x^5
            --------------------------------------------------------------------------------------
            Chi^2/doF = 1,4854228443647e+00
            R^2 = 0,9999990418201
            Adjusted R^2 = 0,9999990229086
            RMSE (Root Mean Squared Error) = 1,218779243491
            RSS (Residual Sum of Squares) = 453,0539675312
            ---------------------------------------------------------------------------------------

            Polynomial Fit of dataset: Table1_H, using function: a0+a1*x+a2*x^2+a3*x^3+a4*x^4+a5*x^5
            --------------------------------------------------------------------------------------
            Chi^2/doF = 1,7949865589785e-01
            R^2 = 0,999998895674
            Adjusted R^2 = 0,9999988738781
            RMSE (Root Mean Squared Error) = 0,4236728170391
            RSS (Residual Sum of Squares) = 54,74709004884
            ---------------------------------------------------------------------------------------
            """
            L=-1.2756691045812e+01+6.1466654146190*self.v210.cSt+9.9774520581931e-01*self.v210.cSt**2-2.5045430263656e-03*self.v210.cSt**3+3.1748553181177e-05*self.v210.cSt**4-1.8264076604682e-07*self.v210.cSt**5
            H=-7.6607640162508+5.4845434135144*self.v210.cSt+0.38222987985934*self.v210.cSt**2-4.6556329076069e-03*self.v210.cSt**3+5.1653200038471e-05*self.v210.cSt**4-2.3274903246922e-07*self.v210.cSt**5
        if self.v100.cSt>H:
            VI=(L-self.v100.cSt)/(L-H)*100
        else:
            N=(log10(H)-log10(self.v100.cSt))/log10(self.v210.cSt)
            VI=(10**N-1)/0.00715+100
        return VI




    def ThCond_Liquido_simple(self, T):
        """Método de cálculo de la conductividad térmica de fracciones petrolíferas líquidas a baja presión, API procedure 12A3.1, pag 1149"""
        t=unidades.Temperature(T)
        k=0.07577-4.1e-5*t.F
        return unidades.Conductividad_termica(k, "BtuhftF")

    def ThCond_Liquido_Tsonopoulos(self, T):
        """Método de cálculo de la conductividad térmica de fracciones petrolíferas líquidas a baja presión,
        Tsonopoulos, C., Heidman, J. L., and Hwang, S.-C.,Thermodynamic and Transport Properties of Coal Liquids, An Exxon Monograph, Wiley, New York, 1986."""
        k=0.05351+0.10177*(1-T/self.Tc)**(2./3)
        return unidades.Conductividad_termica(k)

    def ThCond_Liquido_API(self, T):
        """Método de cálculo de la conductividad térmica de fracciones petrolíferas en líquidas a baja presión, API procedure 12A3.2, pag 1151"""
        t=unidades.Temperature(T)
        k=self.Tb.R**0.2904*(9.961e-3-5.364e-6*t.F)
        return unidades.Conductividad_termica(k, "BtuhftF")

    def ThCond_Liquido_Riazi_Faghri(self, T):
        """Método de cálculo de la conductividad térmica de fracciones petrolíferas a baja presión,
        Riazi, M. R. and Faghri, A., "Thermal Conductivity of Liquid and Vapor Hydrocarbon Systems: Pentanes and Heavier at Low Pressures," Industrial and Engineering Chemistry, Process Design and Development, Vol. 24, No. 2, 1985, pp. 398-401."""
        t=(1.8*T-460)/100
        A=exp(-4.5093-0.6844*t-0.1305*t**2)
        B=0.3003+0.0918*t+0.01195*t**2
        C=0.1029+0.0894*t+0.0292*t**2
        k=1.7307*A*self.Tb.R**B*self.SG**C
        return unidades.Conductividad_termica(k)

    def ThCond_Liquido_Lenoir(self, T, P, ko=0):
        """Método alternativo para el cálculo de la conductividad de líquidos a alta presión, API procedure 12A4.1, pag 1156
        Opcionalmente puede aceptar otro parametros ko que indica un valor experimental de la conductividad térmica (WmK, así como la temperatura y presión a la que se da, en un array [k,T,P]"""
        Tr=self.tr(T)
        if ko==0:
            k1=self.Conductividad_termica_liquido(T)
            C1=17.77+0.065*self.pr(1)-7.764*Tr-2.065*Tr**2/exp(0.2*self.pr(1))
        else:
            k1=ko[0]
            C1=17.77+0.065*self.pr(ko[2])-7.764*self.tr(ko[1])-2.065*self.tr(ko[1])**2/exp(0.2*self.pr(ko[2]))
        C2=17.77+0.065*self.pr(P)-7.764*Tr-2.065*Tr**2/exp(0.2*self.pr(P))
        k=k1*C2/C1
        return unidades.Conductividad_termica(k)


    def ThCond_Gas(self, T):
        """Método de cálculo de la conductividad térmica de vapores de fracciones petrolíferas a baja presión, API procedure 12B3.1, pag 1168"""
        t=unidades.Temperature(T)
        k=0.0013349+0.24628/self.M+1.1493/self.M**2+t.F*(3.2768e-5+4.1881e-5/self.M+0.0018427/self.M**2)
        return unidades.Conductividad_termica(k, "BtuhftF")

    def ThCond_Gas_Riazi_Faghri(self, T):
        """Método de cálculo de la conductividad térmica de fracciones petrolíferas a baja presión,
        Riazi, M. R. and Faghri, A., "Thermal Conductivity of Liquid and Vapor Hydrocarbon Systems: Pentanes and Heavier at Low Pressures," Industrial and Engineering Chemistry, Process Design and Development, Vol. 24, No. 2, 1985, pp. 398-401."""
        t=(1.8*T-460)/100
        A=exp(21.78-8.07986*t+1.12981*t**2-0.05309*t**3)
        B=-4.13948+1.29924*t-0.17813*t**2+0.00833*t**3
        C=0.19876-0.0312*t-0.00567*t**2
        k=1.7307*A*self.Tb.R**B*self.SG**C
        return unidades.Conductividad_termica(k)


    def Calor_combustion_bruto(self):
        """Método de cálculo del calor de combustión bruto de una fracción petrolífera, API procedure 14A1.3, pag 1236"""
        h=17.672+66.6*self.API-0.316*self.API**2-0.0014*self.API**3
        hhv=h-0.01*h*(self.water+self.S+self.ash)+40.5*self.S
        return unidades.Enthalpy(hhv, "Btulb")


    def Calor_combustion_neto(self):
        """Método de cálculo del calor de combustión neto de una fracción petrolífera, API procedure 14A1.3, pag 1236"""
        h=16.796+54.5*self.API-0.217*self.API**2-0.0019*self.API**3
        lhv=h-0.01*h*(self.water+self.S+self.ash)+40.5*self.S-10.53*self.water
        return unidades.Enthalpy(lhv, "Btulb")




class Crudo(Petroleo):
    """Clase que define una fracción de petroleo a partir de la base de datos

    Parámetros:
        indice: indice de la base de datos
        Cplus: numero de fraccioness en las que dividir la fracción para compuestos más pesados que C6
    """
    kwargs=Petroleo.kwargs.copy()
    kwarg={"indice": 0,
                    "Cplus": 0,

                    "Rgo": 0.0,
                    "gas": None,
                    "water": None}
    kwargs.update(kwarg)
    status=0
    _bool=False
    msg=""
    hasCurve=False
    hasSG=True
    hasRefraction=False

    def __call__(self, **kwargs):
        self.kwargs.update(kwargs)
        if kwargs:
            self._bool=True
        if self.isCalculable():
            self.calculo()

    def isCalculable(self):
        if self.kwargs["indice"]:
            self.status=1
            self.msg=""
            return True
        else:
            self.status=0
            self.msg=QApplication.translate("pychemqt", "Undefined petrol")


    def calculo(self):
        propiedades=crudo[self.kwargs["indice"]]

        API=propiedades[3]
        SG=141.5/(API+131.5)
        PP=unidades.Temperature(propiedades[7], "F")
        v100=propiedades[9]
        Tb=unidades.Temperature((-753.-136*(1.-exp(-0.15*v100))+572*SG-0.0512*v100+PP.R)/0.139, "R")

        self.definicion=1
        self.kwargs["name"]=", ".join(propiedades[0:2])
        self.kwargs["SG"]=SG
        self.kwargs["Tb"]=Tb
        self.kwargs["S"]=propiedades[4]
        self.kwargs["N"]=propiedades[5]
        self.kwargs["v100"]=propiedades[9]
        Petroleo.calculo(self)

        self.vanadium=propiedades[10]
        self.nickel=propiedades[11]
        self.carbonResid=propiedades[12]
        self.asphaltene=propiedades[13]
        self.S2=propiedades[5]
        self.H2S=propiedades[17]
        self.nNeutralization=propiedades[18]
        self.ash=propiedades[20]
        self.salt=propiedades[21]
        self.water=propiedades[19]
        self.NPentane=propiedades[14]
        if propiedades[15]:
            self.reidVP=unidades.Pressure(propiedades[15], "psi")
        else:
            self.reidVP=None
        if propiedades[16]:
            self.FlashP=unidades.Temperature(propiedades[16], "F")
        self.PourP=PP

        self.C1=propiedades[22]/100.
        self.C2=propiedades[23]/100.
        self.C3=propiedades[24]/100.
        self.iC4=propiedades[25]/100.
        self.nC4=propiedades[26]/100.
        self.iC5=propiedades[27]/100.
        self.nC5=propiedades[28]/100.

        SGo=0.7
        SG_=(SG-SGo)/SGo
        B=3.
        A=SG_**3/0.619**3
        Cplus=int(self.kwargs["Cplus"])
        Tbi=[unidades.Temperature(1090-exp(6.9955-0.11193*Nc**(2./3))) for Nc in range(6, Cplus)]
        SGi=[1.07-exp(3.65097-3.8864*Nc**0.1) for Nc in range(6, Cplus)]
        x=[1-1/exp(A/B*(SG-SGo)**B/SGo**B) for SG in SGi]
        Mi=[prop_Riazi_Alsahhaf(1, g, reverse=True) for g in SGi]
        APIi=[141.5/SG-131.5 for SG in SGi]
        Kwi=[Tbi[i].R**(1./3)/SGi[i] for i in range(len(Mi))]
        di=[unidades.Density(prop_Riazi_Alsahhaf(2, M), "gcc") for M in Mi]
        Ii=[prop_Riazi_Alsahhaf(3, M) for M in Mi]
        Tci=[unidades.Temperature(Tbi[i]/prop_Riazi_Alsahhaf(4, M)) for i, M in enumerate(Mi)]
        Pci=[unidades.Pressure(prop_Riazi_Alsahhaf(5, M), "bar") for M in Mi]
        Vci=[unidades.SpecificVolume(1/prop_Riazi_Alsahhaf(6, M), "ccg") for M in Mi]
        Wi=[prop_Riazi_Alsahhaf(7, M) for M in Mi]
        Tensioni=[unidades.Tension(prop_Riazi_Alsahhaf(8, M), "dyncm") for M in Mi]
        ParSoli=[unidades.SolubilityParameter(prop_Riazi_Alsahhaf(9, M), "calcc") for M in Mi]


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
        A=self.CO2+self.H2S
        e=120.*(A**0.9-A**1.6)+15*(self.H2S**0.5-self.H2S**4)
        tpc=self.tpc.R-e
        ppc=self.ppc.psi*tpc/(self.tpc.R+self.H2S*(1-self.H2S)*e)
        return unidades.Temperature(tpc, "R"), unidades.Pressure(ppc, "psi")

    def Critical_Carr_Kobayashi_Burrows(self):
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

if __name__ == '__main__':
#    petroleo=Petroleo()
#    print petroleo.VABP
#    print "WABP: ", petroleo.WABP
#    print "MABP: ", petroleo.MABP
#    print "CABP: ", petroleo.CABP
#    print "MeABP: ", petroleo.MeABP.F
#    print "Peso molecular: ", petroleo.M
#    print "TBP:", petroleo.D86_to_TBP([350, 380, 404, 433, 469])
#    print "TBP:", petroleo.TBP_to_D86(petroleo.D86_to_TBP([350, 380, 404, 433, 469]))
#    print "TBP:", petroleo.D2887_to_TBP([293, 305, 324, 336, 344, 359, 369])
#    print "TBP:", petroleo.D2887_to_D86([293, 305, 324, 336, 344, 359, 369])
#    print petroleo.RhoL_altapresion(unidades.Temperature(68, "F"), 368.4482)
#    s=petroleo.Reduccion_volumetrica(5000/95000., 86.5-30.7)
#    print s*5000
#    t=unidades.Temperature(455, "F")
#    print petroleo.Cp_liquido(t).BtulbF
#    print petroleo.Cp_gas(t).BtulbF
#    print petroleo.Tension_superficial(t).dyncm
#    print petroleo.Viscosidad_Singh(t, 1.38)
#    print petroleo.Viscosidad_Fitzgerald_100()
#    print petroleo.Viscosidad_Fitzgerald_210()
#    print petroleo.Viscosidad_ASTM(t, unidades.Temperature(32, "F"), unidades.Temperature(104, "F"), 1.644, 0.925)
#    print petroleo.Viscosidad_liquido_blend(t, [0.6, 0.4], [unidades.Temperature(50, "F"), unidades.Temperature(122, "F"), 14.22, 4.85], [unidades.Temperature(50, "F"), unidades.Temperature(122, "F"), 163.4, 24.98])

#    print petroleo.Conductividad_termica_liquido_simple(t).BtuhftF
#    print petroleo.Conductividad_termica_liquido(t).BtuhftF
#    print petroleo.Conductividad_termica_gas(t).BtuhftF

#    print PNA_Peng_Robinson(7, 655, 94)

#    gas=Natural_Gas(0.699)
#    print gas.ppc.psi, gas.tpc.R
#    print Z_Papay(1.346, 5.603)
#    print Z_Hall_Yarborough(1.346, 5.603)
#    print Z_Dranchuk_Abu_Kassem(1.346, 5.603)
#    print Z_Dranchuk_Purvis_Robinson(1.346, 5.603)
#    print Z_ShellOil(1.346, 5.603)
#    print Z_Beggs_Brill(1.346, 5.603)
#    print Z_Sarem(1.346, 5.603)
#    print Z_Gopal(1.346, 5.603)

#    petroleo=Petroleo(API=44.4, Tb=unidades.Temperature(319., "F"))
#    print petroleo.SG
#    print petroleo.M, petroleo.tc.F, petroleo.pc.psi, petroleo.vc.ft3lb
#    t=unidades.Temperature(325, "F")
#    print petroleo.Cp_liquido(t).BtulbF

#    petroleo=Petroleo(API=22.5, M=339.7)
#    print petroleo.peso_molecular_pesado(55.1, 5.87)
#    print petroleo.watson

#    petroleo=Petroleo(SG=0.9046, Tb=unidades.Temperature(798, "F"))
#    print petroleo.M
#    print petroleo.refractive_index()
#    print petroleo.composicion_molecular(v100=336)

#    petroleo=Petroleo(SG=0.839, Tb=unidades.Temperature(972, "R"), v100=3)
#    print petroleo.pour_point().R

#    petroleo=Petroleo(SG=0.8304, Tb=unidades.Temperature(570.2, "F"))
#    print petroleo.aniline_point().R, petroleo.watson

#    petroleo=Petroleo(SG=0.853, Tb=unidades.Temperature(414.5, "F"))
#    print petroleo.smoke_point().mm, petroleo.watson

#    petroleo=Petroleo(SG=0.799, Tb=unidades.Temperature(874.5, "R"))
#    print petroleo.freezing_point().R, petroleo.watson

#    petroleo=Petroleo(SG=0.787, Tb=unidades.Temperature(811.5, "R"))
#    print petroleo.cloud_point().R

#    petroleo=Petroleo(API=32.3, Tb=unidades.Temperature(617, "F"))
#    print petroleo.cetane_index()

#    d86=[unidades.Temperature(t, "F") for t in [149, 230, 282, 325, 429]]
#    petroleo=Petroleo(D86=d86)
#    print petroleo.VABP.F
#    print petroleo.MABP.F
#    print petroleo.WABP.F
#    print petroleo.CABP.F
#    print petroleo.MeABP.F
#    print petroleo.pv_fraccion(unidades.Temperature(70, "F"), unidades.Pressure(6, "psi").atm).psi

#    t, bool= Hydrates_Sloan("T", T=300, P=10, y=[0.78, 0.06, 0.03, 0.01, 0.02, 0.06, 0.04, 0.0])
#    print t.F

#    petroleo=Petroleo(API=30.6, Tb=unidades.Temperature(538, "F"))
#    print petroleo.RhoL(unidades.Temperature(160, "F")).gml

#    petroleo=Petroleo(API=31.4, Tb=unidades.Temperature(538, "F"))
#    t=unidades.Temperature(68, "F")
#    p=unidades.Pressure(5400, "psig")
#    print petroleo.RhoL_Presion(t, p.atm).gml

#    petroleo=Petroleo(API=43.5, Tb=unidades.Temperature(407.2, "F"))
#    t=unidades.Temperature(600, "F")
#    p=unidades.Pressure(50, "psi")
#    print petroleo.Entalpia(t, p.atm).Btulb

#    d86=[unidades.Temperature(t, "F") for t in [304, 313, 321, 329, 341]]
#    petroleo=Petroleo(API=44.4, D86=d86)
#    t=unidades.Temperature(325, "F")
#    print petroleo.Cp_liquido(t).BtulbF

#    d86=[unidades.Temperature(t, "F") for t in [304, 313, 321, 329, 341]]
#    petroleo=Petroleo(API=44.4, D86=d86)
#    t=unidades.Temperature(885, "F")
#    p=unidades.Pressure(205, "psi")
#    print petroleo.Cp_gas(t, p.atm).BtulbF

#    R=unidades.V2V(675, "ft3bbl")
#    gas=Natural_Gas(SG=0.95, CO2=0.2, H2S=0.1)
#    petroleo=Petroleo(API=31., Rgo=R, gas=gas)
#    t=unidades.Temperature(180, "F")
#    print petroleo.pb_Standing(t).psi
#    print petroleo.pb_Lasater(t).psi
#    print petroleo.pb_Vazquez_Beggs(t, unidades.Temperature(85, "F"), unidades.Pressure(100, "psi").atm).psi
#    print petroleo.pb_Glaso(t).psi
#    print petroleo.pb_Total(t).psi
#    print petroleo.pb_Al_Marhoun(t).psi
#    print petroleo.pb_Dokla_Osman(t).psi
#    print petroleo.pb_Petrosky_Farshad(t).psi
#    print petroleo.pb_Kartoatmodjo_Schmidt(t, unidades.Temperature(85, "F"), unidades.Pressure(100, "psi").atm).psi

#    R=unidades.V2V(675, "ft3bbl")
#    t=unidades.Temperature(180, "F")
#    p=unidades.Pressure(4000, "psi")
#    petroleo=Petroleo(API=31.)
#    print petroleo.Mu_Beal(t).cP
#    print petroleo.Mu_Beggs_Robinson(t).cP
#    print petroleo.Mu_Glaso(t).cP
#    print petroleo.Mu_Egbogad(t).cP
#    print petroleo.Mu_Kartoatmodjo_Schmidt(t).cP
#    print petroleo.Mu_Chew_Connally(t, R).cP
#    print petroleo.Mu_Beggs_Robinson_vivo(t, R).cP
#    print petroleo.Mu_Kartoatmodjo_Schmidt_vivo(t, R).cP

#    petroleo=Petroleo(API=33., M=250)
#    print petroleo.Viscosidad_gas(unidades.Temperature(100, "F")).cP


#    agua=Water()
#    t=unidades.Temperature(200, "F")
#    p=unidades.Pressure(5000, "psi")
#    print agua.Solubilidad_Culberson_McKetta(t, p.atm, 2).ft3bbl
#    print agua.Solubilidad_McCoy(t, p.atm, 2).ft3bbl

#    agua=Water()
#    t=unidades.Temperature(200, "F")
#    p=unidades.Pressure(5000, "psi")
#    print agua.Factor_Volumetrico_McCain(t, p.atm, 2)
#    print agua.Factor_Volumetrico_McCoy(t, p.atm, 2)

#    R=unidades.V2V(17.8, "ft3bbl")
#    print agua.Compresibilidad_Dodson_Standing(t, p.atm, 2, R)
#    print agua.Compresibilidad_Osif(t, p.atm, 2)

#    print agua.Mu_Van_Wingen(t).cP
#    print agua.Mu_Mattews_Russel(t, p.atm, 2).cP
#    print agua.Mu_McCain(t, p.atm, 2).cP
#    print agua.Mu_McCoy(t, 2).cP
#    print agua.Tension_Jennings_Newman(t, p.atm).dyncm
#    print agua.Rho(t, p.atm, 2).lbft3
#    print agua.Rho_McCain(t, p.atm, 2).lbft3

#    t=unidades.Temperature(200, "F")
#    print SUS(t, 53.)
#    t=unidades.Temperature(0, "F")
#    print SUS(t, 90.)
#    print SUF(122, 300)
#    print SUF(210, 100)


#    petroleo=Petroleo(v100=73.3, v210=8.86)
#    print petroleo.VI
#    petroleo=Petroleo(v100=5000, v210=100)
#    print petroleo.VI
#    petroleo=Petroleo(v100=22.83, v210=5.05)
#    print petroleo.VI
#    petroleo=Petroleo(v100=1500, v210=100)
#    print petroleo.VI

#    v=[10, 30, 50, 70, 90]
#    D1160=[i+273.15 for i in [150, 205, 250, 290, 350]]
#    petroleo=Petroleo(P_dist=10, T_dist=v, D1160=D1160)



#        SG=0.867566
#        SGo=0.7
#        SG_=(SG-SGo)/SGo
#        B=3.
#        A=SG_**3/0.619**3
#        x=arange(0.05, 0.96, 0.05)
#        print x
#        SGi=[SGo+SGo*(A/B*log10(1/(1-xi)))**(1./B) for xi in x]
#        Mi=[prop_Riazi_Alsahhaf(1, g, reverse=True) for g in SGi]
#        print SGi
#        print Mi

    petroleo=Petroleo(M=250, API=43.3)
    t=unidades.Temperature(100, "F")
#    print petroleo.Mu_Gas(t).cP
#    print petroleo.Tension_API(t).dyncm
#    print petroleo.Tension_Baker_Swerdloff(t).dyncm
#    print petroleo.API

    petroleo=Petroleo(API=22.5, M=339.7)
    print(petroleo.API, petroleo.M)
