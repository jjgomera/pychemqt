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
along with this program.  If not, see <http://www.gnu.org/licenses/>.


Petroleum fractions hypotethical pseudocomponent definition

:func:`Petroleo`: The main class with all integrated functionality

Correlation for calculation of properties:

    * :func:`prop_Ahmed`
    * :func:`prop_Riazi_Daubert_1980`
    * :func:`prop_Riazi_Daubert`
    * :func:`prop_Cavett`
    * :func:`prop_Lee_Kesler`
    * :func:`prop_Sim_Daubert`
    * :func:`prop_Watansiri_Owens_Starling`
    * :func:`prop_Rowe`
    * :func:`prop_Standing`
    * :func:`prop_Willman_Teja`
    * :func:`prop_Magoulas_Tassios`
    * :func:`prop_Tsonopoulos`
    * :func:`prop_Twu`
    * :func:`prop_Sancet`
    * :func:`prop_Silva_Rodriguez`
    * :func:`Tb_Soreide`
    * :func:`vc_Hall_Yarborough`
    * :func:`M_Goossens`
    * :func:`w_Korsten`
    * :func:`prop_Riazi`
    * :func:`prop_Riazi_Alsahhaf`
    * :func:`prop_Riazi_Alsahhaf_PNA`

Correlations for Zc calculations:
    * :func:`Zc_Hougen`
    * :func:`Zc_Reid`
    * :func:`Zc_Salerno`
    * :func:`Zc_Nath`
    * :func:`Zc_Lee_Kesler`

Distillation curves interconversiono methods:
    * :func:`D86_TBP_Riazi`
    * :func:`D86_TBP_Daubert`
    * :func:`SD_D86_Riazi`
    * :func:`SD_D86_Daubert`
    * :func:`D86_EFV`
    * :func:`SD_TBP`
    * :func:`D1160_TBP_10mmHg`
    * :func:`Tb_Pressure`
    * :func:`curve_Predicted`
    * :func:`_Tb_Predicted`

Advanced petroleum fraction properties:
    * :func:`PourPoint`
    * :func:`AnilinePoint`
    * :func:`SmokePoint`
    * :func:`FreezingPoint`
    * :func:`CloudPoint`
    * :func:`CetaneIndex`
    * :func:`H_Riazi`
    * :func:`H_Goossens`
    * :func:`H_ASTM`
    * :func:`H_Jenkins_Walsh`
    * :func:`viscoAPI`
    * :func:`SUS`
    * :func:`SFS`
    * :func:`S_Riazi`

PNA decomposition procedures:
    * :func:`PNA_Riazi`
    * :func:`PNA_Peng_Robinson`
    * :func:`PNA_Bergman`
    * :func:`PNA_van_Nes`
'''


from configparser import ConfigParser

from PyQt5.QtWidgets import QApplication
from scipy import exp, sqrt, log10, log
from scipy.interpolate import interp1d
from scipy.optimize import fsolve, leastsq, newton
from numpy.linalg import solve
from numpy import array

from lib import unidades
from lib.physics import R_atml, R_Btu
from lib.newComponent import newComponente
from lib.config import conf_dir
from lib.compuestos import prop_Edmister


__doi__ = {
    1:
        {"autor": "Katz, D.L., Firoozabadi, A.",
         "title": "Predicting Phase Behavior of Condensate/Crude-oil Systems"
                  "Using Methane Interaction Coefficients",
         "ref": "Journal of Petroleum Technology 30(11) (1978) 1649-1655",
         "doi": "10.2118/6721-pa"},
    2:
        {"autor": "Ahmed, T., Cady, G., Story, A.",
         "title": "A Generalized Correlation for Characterizing the"
                  "Hydrocarbon Heavy Fractions.",
         "ref": "Paper SPE 14266, presented at the 60th Annual Technical"
                "Conference of the Society of Petroleum Engineers, Las Vegas,"
                "September 22–25, 1985.",
         "doi": ""},
    3:
        {"autor": "Riazi, M.R., Daubert, T.E.",
         "title": "Characterization Parameters for Petroleum Fractions",
         "ref": "Ind. Eng. Chem. Res. 26(4) (1987) 755-759",
         "doi": "10.1021/ie00064a023"},
    4:
        {"autor": "Riazi, M.R., Daubert, T.E.",
         "title": "Simplify Property Predictions",
         "ref": "Hydrocarbon Processing (March 1980): 115–116",
         "doi": ""},
    5:
        {"autor": "Twu, C.H.",
         "title": "An Internally Consistent Correlation for Predicting the"
                  "Critical Properties and Molecular Weights of Petroleum and"
                  "Coal-tar Liquids",
         "ref": "Fluid Phase Equilbria 16 (1984) 137-150",
         "doi": "10.1016/0378-3812(84)85027-x"},
    6:
        {"autor": "Sim, W.J., Daubert, T.E.",
         "title": "Prediction of Vapor-Liquid Equilibria of Undefined"
                  "Mixtures",
         "ref": "Ind. Eng. Chem. Process Des. Dev. 19(3) (1980) 386-393",
         "doi": "10.1021/i260075a010"},
    7:
        {"autor": "Cavett, R.H.",
         "title": "Physical data for distillation calculations, vapor-liquid"
                  "equilibrium.",
         "ref": "Proceedings of the 27th Meeting, API, San Francisco, Issue 3,"
                "351-366.",
         "doi": ""},
    8:
        {"autor": "Kesler, M.G., Lee, B.I.",
         "title": "Improve Prediction of Enthalpy of Fractions",
         "ref": "Hydrocarbon Processing 55(3) (1976) 153-158",
         "doi": ""},
    9:
        {"autor": "Ahmed, T.",
         "title": "Equations of State and PVT Analysis: Applications for"
                  "Improved Reservoir Modeling, 2nd Edition",
         "ref": "Gulf Professional Publishing, 2016, ISBN 9780128015704,",
         "doi": "10.1016/B978-0-12-801570-4.00002-7"},
    10:
        {"autor": "Watansiri, S., Owens, V.H., Starling, K.E.",
         "title": "Correlations for estimating critical constants, acentric"
                  "factor, and dipole moment for undefined coal-fluid"
                  "fractions",
         "ref": "Ind. Eng. Chem. Process. Des. Dev. 24(2) (1985) 294-296",
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
         "ref": "Ind. Eng. Chem. Res. 26(5) (1987) 948-952",
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
        {"autor": "Hall, K. R., Yarborough, L.",
         "title": "New Simple Correlation for Predicting Critical Volume",
         "ref": "Chemical Engineering (November 1971): 76",
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
         "ref": "Fluid Phase Equilibria 27 (1986) 15-34",
         "doi": "10.1016/0378-3812(86)87038-8"},
    22:
        {"autor": "Hougen, O.A., Watson, K.M., Ragatz, R.A.",
         "title": "Chemical Process Principles, 2nd ed.",
         "ref": "New York: Wiley, 1959, p. 577.",
         "doi": ""},
    23:
        {"autor": "Reid, R., Prausnitz, J.M., Sherwood, T.",
         "title": "The Properties of Gases and Liquids, 3rd ed. New York:"
                  "McGraw-Hill, 1977, p. 21.",
         "ref": "",
         "doi": ""},
    24:
        {"autor": "Nath, J.",
         "title": "Acentric Factor and the Critical Volumes for Normal Fluids",
         "ref": "Ind. Eng. Chem. Fundam. 21(3) (1985) 325-326",
         "doi": "10.1021/i100007a023"},
    25:
        {"autor": "Lee, B.I., Kesler, M.G.",
         "title": "A Generalized Thermodynamic Correlation Based on"
                  "Three-Parameter Corresponding States",
         "ref": "AIChE J. 21(3) (1975) 510-527",
         "doi": "10.1002/aic.690210313"},
    26:
        {"autor": "Riazi, M.R., Al-Sahhaf, T.A., Sl-Shammari M.A.",
         "title": "A Generalized Method for Estimation of Critical Constants",
         "ref": "Fluid Phase Equilibria 147 (1998) 1-6",
         "doi": "10.1016/s0378-3812(98)00251-9"},
    27:
        {"autor": "Riazi, M.R., A1-Sahhaf, T.A.",
         "title": "Physical Properties of Heavy Petroleum Fractions and Crude"
                  "Oils",
         "ref": "Fluid Phase Equilibria 117 (1996) 217-224.",
         "doi": "10.1016/s0378-3812(98)00251-9"},
    28:
        {"autor": "Riazi, M.R., A1-Sahhaf, T.",
         "title": "Physical Properties of n-Alkanes and n-Alkyl Hydrocarbons: "
                  "Application to Petroleum Mixtures",
         "ref": "Ind. Eng. Chem. Res. 34(11) (1995) 4145-4148",
         "doi": "10.1021/ie00038a062"},
    29:
        {"autor": "Riazi, M. R.",
         "title": "Characterization and Properties of Petroleum Fractions.",
         "ref": "ASTM manual series MNL50, 2005",
         "doi": ""},
    30:
        {"autor": "ASTM D2161-05",
         "title": "Standard Practice for Conversion of Kinematic Viscosity to "
                  "Saybolt Universal Viscosity or to Saybolt Furol Viscosity",
         "ref": "ASTM International, West Conshohocken, PA 2005, www.astm.org",
         "doi": "10.1520/D2161-05"},
    31:
        {"autor": "Singh, B., Mutyala, S.R., Puttagunta, V.",
         "title": "Viscosity range from one test",
         "ref": "Hydrocarbon Processing 69 (1990) 39-41 ",
         "doi": ""},

    32:
        {"autor": "",
         "title": "",
         "ref": "",
         "doi": ""},

    51:
        {"autor": "Riazi, M.R., Daubert, T.E.",
         "title": "Analytical Correlations Interconvert Distillation Curve "
                  "Types",
         "ref": "Oil & Gas Journal 84 (1986) 50-57",
         "doi": ""},
    52:
        {"autor": "Daubert, T.E.",
         "title": "Petroleum Fraction Distillation Interconversion",
         "ref": "Hydrocarbon Processing 73(9) (1994) 75-78",
         "doi": ""},
    53:
        {"autor": "Edmister, W.C., Okamoto, K.K.",
         "title": "Applied Hydrocarbon Thermodynamics, Part 13: Equilibrium "
                  "Flash Vaporization Correlations for Heavy Oils Under "
                  "Subatmospheric Pressures",
         "ref": "Petroleum Refiner 38(9) (1959) 271-288",
         "doi": ""},
    54:
        {"autor": "Riazi, M.R.",
         "title": "Distribution Model for Properties of Hydrocarbon-Plus "
                  "Fractions",
         "ref": "Ind. Eng. Chem. Res. 28(11) (1989) 1731-1735.",
         "doi": "10.1021/ie00095a026"},
    55:
        {"autor": "Van Nes, K., Van Western, H.A.",
         "title": "Aspects of the Constitution of Mineral Oils",
         "ref": "Elsevier, New York, 1951",
         "doi": ""},
    56:
        {"autor": "Ahmed, T.",
         "title": "Hydrocarbon Phase Behavior",
         "ref": "Gulf Publishing, Houston, TX, 1989.",
         "doi": ""},
    57:
        {"autor": "Bergman, D.F., Tek, M.R., Katz, D.L.",
         "title": "Retrograde Condensation in Natural Gas Pipelines",
         "ref": "Project PR 2-29 of Pipelines Research Committee, AGA, January"
                " 1977",
         "doi": ""},
    58:
        {"autor": "Robinson, D.B., Peng, D.Y.",
         "title": "The characterization of the heptanes and heavier fractions",
         "ref": "Research Report 28. GPA, 1978. Tulsa, OK.",
         "doi": ""},
    59:
        {"autor": "Riazi, M.R., Nasimi, N., Roomi, Y.",
         "title": "Estimating Sulfur Content of Petroleum Products and Crude"
                  " Oils",
         "ref": "Ind. Eng. Chem. Res. 38(11) (1999) 4507-4512",
         "doi": "10.1021/ie990262d"},
    60:
        {"autor": "Goossens, A.G.",
         "title": "Prediction of the Hydrogen Content of Petroleum Fractions",
         "ref": "Ind. Eng. Chem. Res. 36(6) (1997) 2500-2504",
         "doi": "10.1021/ie960772x"},
    61:
        {"autor": "Jenkins, G.I., Walsh, R.E",
         "title": "Quick Measure of Jet Fuel Properties",
         "ref": "Hydrocarbon Processing 47(5) (1968) 161-164",
         "doi": ""},
    62:
        {"autor": "ASTM",
         "title": "Annual Book of Standards",
         "ref": "ASTM International, West Conshohocken, PA, 2002",
         "doi": ""},
    63:
        {"autor": "Goossens, A.G.",
         "title": "Prediction of Molecular Weight of Petroleum Fractions",
         "ref": "Ind. Eng. Chem. Res. 35(3) (1996) 985-988",
         "doi": "10.1021/ie950484l"},
    64:
        {"autor": "Korsten, H.",
         "title": "Internally Consistent Prediction of Vapor Pressure and "
                  "Related Properties",
         "ref": "Ind. Eng. Chem. Res. 39(3) (2000) 813-820",
         "doi": "10.1021/ie990579d"},
    65:
        {"autor": "Tsonopoulos, C., Heidman, J.L., Hwang, S.C.",
         "title": "Thermodynamic and Transport Properties of Coal Liquids",
         "ref": "An Exxon Monograph, Wiley, New York, 1986",
         "doi": ""},


    66:
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
    prop : dict
        A dict with the calculated properties:

            * M: Molecular weight, [-]
            * Tc: Critic temperature, [ºR]
            * Pc: Critic pressure, [psi]
            * Tb: Normal boiling temperature, [ºR]
            * w: Acentric factor, [-]
            * SG: Specific gravity, [-]
            * Vc: Critic volume, [ft³/lb]

    Examples
    --------
    >>> "%0.0f" % prop_Ahmed(6)["Tb"].R
    '604'
    >>> "%0.5f" % prop_Ahmed(45)["Vc"].ft3lb
    '0.06549'

    References
    ----------
    [1]_ Katz, D.L., Firoozabadi, A. Predicting Phase Behavior of
    Condensate/Crude-oil Systems Using Methane Interaction Coefficients.
    Journal of Petroleum Technology 30(11) (1978) 1649-1655

    [2]_ Ahmed, T., Cady, G., Story, A. A Generalized Correlation for
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
    prop : dict
        A dict with the calculated properties:

            * M: Molecular weight, [-]
            * Tc: Critic temperature, [ºR]
            * Pc: Critic pressure, [psi]
            * Vc: Critic volume, [ft3/lb]

    Examples
    --------
    Example 2.2 from [9]_: C7+ fraction with Tb=198ºF and SG=0.7365

    >>> T = unidades.Temperature(198, "F")
    >>> p = prop_Riazi_Daubert_1980(T, 0.7365)
    >>> "%.0f %.0f %.0f %.4f" % (p["M"], p["Tc"].R, p["Pc"].psi, p["Vc"].ft3lb)
    '96 990 467 0.0623'

    References
    ----------
    [4]_ Riazi, M.R., Daubert, T.E. Simplify Property Predictions.
    Hydrocarbon Processing (March 1980): 115–116.

    [9]_ Ahmed, T. Equations of State and PVT Analysis: Applications for
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
    prop : dict
        A dict with the calculated properties:

            * Tc: Temperatura crítica, [ºR]
            * Pc: Presión crítica, [psi]
            * Vc: Volumen crítico, [ft³/lb]
            * M: Peso molecular, [-]
            * Tb: Temperatura fusión, [ºR]
            * SG: Gravedad específica, [-]
            * I: Huang characterization factor, [-]
            * CH: Carbon/hydrogen weight ratio, [-]

    Notes
    -----
    The available input properties for tita are:

        * Tb: Temperatura fusión, [ºR]
        * SG: Gravedad específica, [-]
        * M: Peso molecular, [-]
        * CH: Carbon/hydrogen weight ratio, [-]
        * I: Huang characterization factor, [-]
        * v1: Kinematic viscosity at 100ºF, [m²/s]

    Raise :class:`NotImplementedError` if input pair are unsupported or input
    isn't in limit:

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
    [3]_ Riazi, M.R., Daubert, T.E. Characterization Parameters for
    Petroleum Fractions. Ind. Eng. Chem. Res. 26(4) (1987) 755-759

    [9]_ Ahmed, T. Equations of State and PVT Analysis: Applications for
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
    prop : dict
        A dict with the calculated properties:

            * Tc: Critic temperature, [ªR]
            * Pc: Critic pressure, [Pa]

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
    [7]_ Cavett, R.H., 1962. Physical data for distillation calculations,
    vapor-liquid equilibrium. Proceedings of the 27th Meeting, API, San
    Francisco, 1962, Issue 3, pp. 351–366.

    [9]_ Ahmed, T. Equations of State and PVT Analysis: Applications for
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
    prop : dict
        A dict with the calculated properties:

            * Tc: Critic temperature, [ªR]
            * Pc: Critic pressure, [Pa]
            * M: Molecular weight, [-]
            * w: Acentric factor, [-]

    Examples
    --------
    Example 2.2 from [9]_: Petroleum fraction with Tb=198ºF and SG=0.7365

    >>> T = unidades.Temperature(198, "F")
    >>> p = prop_Lee_Kesler(T, 0.7365)
    >>> "%.0f %.0f %.1f %.3f" % (p["Tc"].R, p["Pc"].psi, p["M"], p["w"])
    '981 470 98.6 0.306'

    References
    ----------
    [8]_ Kesler, M.G., Lee, B.I. Improve Prediction of Enthalpy of
    Fractions. Hydrocarbon Processing 55(3) (1976) 153-158

    [9]_ Ahmed, T. Equations of State and PVT Analysis: Applications for
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
    prop : dict
        A dict with the calculated properties:

            * M: Molecular weight, [-]
            * Tc: Critic temperature, [ªR]
            * Pc: Critic pressure, [Pa]

    Examples
    --------
    Example 2.2 from [9]_: Petroleum fraction with Tb=198ºF and SG=0.7365

    >>> T = unidades.Temperature(198, "F")
    >>> p = prop_Sim_Daubert(T, 0.7365)
    >>> "%.0f %.0f %.0f" % (p["Tc"].R, p["Pc"].psi, p["M"])
    '979 479 96'

    References
    ----------
    [6]_ Sim, W.J., Daubert, T.E. Prediction of vapor-liquid
    equilibria of undefined mixtures. Ind. Eng. Chem. Process Des. Dev. 19(3)
    (1980) 386-393

    [9]_ Ahmed, T. Equations of State and PVT Analysis: Applications for
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
    prop : dict
        A dict with the calculated properties:

            * Tc: Critic temperature, [K]
            * Pc: Critic pressure, [atm]
            * Vc: Critic volume, [cm3/gr]
            * w: Acentric factor, [-]
            * Dm: Dipole moment, [Debye]

    Examples
    --------
    Example 2.2 from [9]_: Petroleum fraction with Tb=198ºF and SG=0.7365

    >>> T = unidades.Temperature(198, "F")
    >>> p = prop_Watansiri_Owens_Starling(T, 0.7365, 96)
    >>> "%.0f %.5f" % (p["Tc"].R, p["Vc"].ft3lb)
    '980 0.06548'

    References
    ----------
    [10]_ Watansiri, S., Owens, V.H., Starling, K.E., 1985. Correlations for
    estimating critical constants, acentric factor, and dipole moment for
    undefined coal-fluid fractions. Ind. Eng. Chem. Process. Des. Dev. 24(2)
    (1985) 294-296

    [9]_ Ahmed, T. Equations of State and PVT Analysis: Applications for
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
    """Calculate petroleum fractions properties with the Rowe
    correlations (1978) using only the molecular weight as input parameter

    Parameters
    ------------
    M: float
        Molecular weight, [-]

    Returns
    -------
    prop : dict
        A dict with the calculated properties:

            * Tc: Critic temperature, [R]
            * Tb: Boiling temperature, [R]
            * Pc: Critic pressure, [psi]

    Examples
    --------
    Example 2.3 from [9]_: Petroleum fraction with M=216 and SG=0.8605

    >>> p = prop_Rowe(216)
    >>> "%.1f" % p["Tc"].R
    '1279.8'

    References
    ----------
    [11]_ Rowe, A.M. Internally consistent correlations for predicting phase
    compositions for use in reservoir compositional simulators. Paper SPE
    7475, In: Presented at the 53rd Annual Society of Petroleum Engineers
    Fall Technical Conference and Exhibition, 1978.

    [9]_ Ahmed, T. Equations of State and PVT Analysis: Applications for
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
    prop : dict
        A dict with the calculated properties:
            * Tc: Critic temperature, [K]
            * Pc: Critic pressure, [atm]

    Examples
    --------
    Example 2.3 from [9]_: Petroleum fraction with M=216 and SG=0.8605

    >>> p = prop_Standing(216, 0.8605)
    >>> "%.1f %.0f" % (p["Tc"].R, p["Pc"].psi)
    '1269.3 270'

    References
    ----------
    [12]_ Standing, M.B., Volumetric and Phase Behavior of Oil Field
    Hydrocarbon Systems. Society of Petroleum Engineers, Dallas, TX. 1977

    [9]_ Ahmed, T. Equations of State and PVT Analysis: Applications for
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
    prop : dict
        A dict with the calculated properties:

            * Tc: Critic temperature, [K]
            * Pc: Critic pressure, [MPa]
            * n: Carbon number, [-]

    Examples
    --------
    Example 2.2 from [9]_: Petroleum fraction with Tb=198ºF and SG=0.7365

    >>> Tb = unidades.Temperature(198, "F")
    >>> p = prop_Willman_Teja(Tb)
    >>> "%.0f %.0f" % (p["Tc"].R, p["Pc"].psi)
    '977 442'

    References
    ----------
    [9]_ Ahmed, T. Equations of State and PVT Analysis: Applications for
    Improved Reservoir Modeling, 2nd Edition. Gulf Professional Publishing,
    2016, ISBN 9780128015704

    [13]_ Willman, B., Teja, A. Prediction of dew points of semicontinuous
    natural gas and petroleum mixtures. 1. Characterization by use of an
    effective carbon number and ideal solution predictions. Ind. Eng. Chem.
    Res. 26(5) (1987) 948-952
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
    prop : dict
        A dict with the calculated properties:

            * Tc: Critic temperature, [K]
            * Pc: Critic pressure, [MPa]
            * w: Acentric factor, [-]

    Examples
    --------
    Example 2.3 from [9]_: Petroleum fraction with M=216 and SG=0.8605

    >>> p = prop_Magoulas_Tassios(216, 0.8605)
    >>> "%.0f %.0f" % (p["Tc"].R, p["Pc"].psi)
    '1317 273'

    References
    ----------
    [9]_ Ahmed, T. Equations of State and PVT Analysis: Applications for
    Improved Reservoir Modeling, 2nd Edition. Gulf Professional Publishing,
    2016, ISBN 9780128015704

    [14]_ Magoulas, S., Tassios, D., Predictions of phase behavior of HT-HP
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


def prop_Tsonopoulos(SG, Tb):
    """Calculate petroleum fractions properties with the Tsonopoulos (1990)
    correlations using molecular weight and specific gravity as input parameter

    Parameters
    ------------
    SG: float
        Specific gravity, [-]
    Tb: float
        Boiling temperature, [K]

    Returns
    -------
    prop : dict
        A dict with the calculated properties:

            * Tc: Critic temperature, [K]
            * Pc: Critic pressure, [MPa]

    References
    ----------
    [65]_ Tsonopoulos, C., Heidman, J.L., Hwang, S.C. Thermodynamic and
    Transport Properties of Coal Liquids. An Exxon Monograph, Wiley, New
    York, 1986
    """
    # TODO: Search reference
    Tc = 10**(1.20016 - 0.61954*log10(Tb) + 0.48262*log10(SG) +
              0.67365*log10(SG)**2)
    Pc = 10**(7.37498 - 2.15833*log10(Tb) + 3.35417*log10(SG) +
              5.64019*log10(SG)**2)

    prop = {}
    prop["Tb"] = unidades.Temperature(Tb)
    prop["SG"] = unidades.Dimensionless(SG)
    prop["Tc"] = unidades.Temperature(Tc)
    prop["Pc"] = unidades.Pressure(Pc, "bar")
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
    prop : dict
        A dict with the calculated properties:

            * M: Molecular weight, [-]
            * Tc: Critic temperature, [ºR]
            * Pc: Critic pressure, [psi]
            * Vc: Critic volume, [ft3/lb]

    Examples
    --------
    >>> crit = prop_Twu(510, 1.097)
    >>> "%.1f %.1f %.1f" % (crit["Tc"].R, crit["Pc"].psi, crit["M"])
    '1380.3 556.8 130.4'

    References
    ----------
    [5]_ Twu, C.H. An internally consistent correlation for predicting the
    critical properties and molecular weights of petroleum and coal-tar
    liquids. Fluid Phase Equilbria 16 (1984) 137-150
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
    prop : dict
        A dict with the calculated properties:

            * Tc: Critic temperature, [ºR]
            * Pc: Critic pressure, [psi]
            * Tb: Boiling temperature, [ºR]

    Examples
    --------
    Example 2.1 from [9]_: C7+ fraction with M=150 and SG=0.78

    >>> p = prop_Sancet(150)
    >>> "%.0f %.0f %.0f" % (p["Tc"].R, p["Pc"].psi, p["Tb"].R)
    '1133 297 828'

    References
    ----------
    [5]_ Sancet, J., Heavy Faction  :math:`C_{7+}` Characterization for
    PR-EOS. In: SPE 113026, 2007 SPE Annual Conference, November 11–14,
    Anaheim, CA 2007.

    [9]_ Ahmed, T. Equations of State and PVT Analysis: Applications for
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
    prop : dict
        A dict with the calculated properties:

            * Tb: Boiling temperature, [K]
            * SG: Specific gravity, [-]

    Examples
    --------
    Example 2.1 from [9]_: C7+ fraction with M=150 and SG=0.78

    >>> p = prop_Silva_Rodriguez(150)
    >>> "%.0f %.4f" % (p["Tb"].R, p["SG"])
    '839 0.7982'

    References
    ----------
    [16]_ Silva, M.B., Rodriguez, F. Automatic fitting of equations of state
    for phase behavior matching. Paper SPE 23703, Society of Petroleum
    Engineers, Richardson, TX, 1992.

    [9]_ Ahmed, T. Equations of State and PVT Analysis: Applications for
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
    [17]_ Soreide, I. Improved Phase Behavior Predictions of Petroleum
    Reservoir Fluids From a Cubic Equation of State. Doctor of engineering
    dissertation. Norwegian Institute of Technology, Trondheim, 1989.
    """
    Tb = 1928.3 - 1.695e5*SG**3.266/M**0.03522*exp(
        -4.922e-3*M-4.7685*SG+3.462e-3*M*SG)                          # Eq 3.59
    return {"Tb": unidades.Temperature(Tb, "R")}


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
    [18]_ Hall, K.R., Yarborough, L. New Simple Correlation for Predicting
    Critical Volume. Chemical Engineering (November 1971): 76.
    """
    Vc = 0.025*M**1.15/SG**0.7985
    return {"Vc": unidades.SpecificVolume(Vc/M, "ft3lb")}


def M_Goossens(Tb, d20):
    """Calculate petroleum fractions molecular weight with the Goossens
    (1971) correlation

    Parameters
    ------------
    Tb : float
        Normal boiling temperature, [K]
    d20 : float
        Liquid density at 20ºC and 1 atm, [g/cm³]

    Returns
    -------
    M: float
        Molecular weight, [-]

    Examples
    --------
    >>> "%.1f" % M_Goossens(306, 0.6258)["M"]
    '77.0'

    References
    ----------
    [63]_ Goossens, A.G. Prediction of Molecular Weight of Petroleum
    Fractions. Ind. Eng. Chem. Res. 35(3) (1996) 985-988
    """
    b = 1.52869 + 0.06486*log(Tb/(1078-Tb))
    M = 0.01077*Tb**b/d20
    return {"M": unidades.Dimensionless(M)}


def w_Korsten(Tb, Tc, Pc):
    """Calculate petroleum fractions acentric factor with the Korsten (2000)
    correlation

    Parameters
    ------------
    Tb : float
        Normal boiling temperature, [K]
    Tc : float
        Critical temperature, [K]
    Pc : float
        Critical pressure, [Pa]

    Returns
    -------
    w: float
        Acentric factor, [-]

    References
    ----------
    [64]_ Korsten, H. Internally Consistent Prediction of Vapor Pressure and
    Related Properties. Ind. Eng. Chem. Res. 39(3) (2000) 813-820
    """
    tbr = Tb/Tc
    # Eq 29
    w = 0.5899*tbr**1.3/(1-tbr**1.3)*log10(Pc/101325)-1.
    return {"w": unidades.Dimensionless(w)}


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
    prop : dict
        A dict with the calculated properties:

            * Tc: Critic temperature, [ºR]
            * Pc: Critic pressure, [psi]
            * Vc: Critic volume, [ft³/lb]
            * Tb: Boiling temperatura, [ºR]
            * d20: Liquid density at 20ºC, [g/cm³]
            * I: Huang Characterization factor, [-]

    Notes
    -----
    The available input properties for tita are:

        * Tb: Boiling temperature, [ºR]
        * M: Molecular weight, [-]

    The critic volume is calculate in mole base, so be careful to correct with
    molecular weight

    References
    ----------
    [29]_ Riazi, M. R. Characterization and Properties of Petroleum Fractions.
    ASTM manual series MNL50, 2005
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


def prop_Riazi_Alsahhaf(Tb, M, d20):
    """Calculate petroleum fractions properties with the Riazi-AlSahhaf (1998)
    correlation with the boiling temperature and liquid density at 20ºC as
    input paramters.

    Parameters
    ------------
    Tb : float
        Normal boiling temperature, [K]
    M : float
        Molecular Weight, [g/mol]
    d20 : float
        Liquid density at 20ºC and 1 atm, [g/cm³]

    Returns
    -------
    prop : dict
        A dict with the calculated properties:

            * Tc: Critic temperature, [K]
            * Pc: Critic pressure, [MPa]
            * Vc: Critic volume, [cm³/g]

    Notes
    -----
    This correlation generalized the Riazi-Daubert method to non-polar
    compounds and gases.

    References
    ----------
    [26]_ Riazi, M.R., Al-Sahhaf, T.A., Sl-Shammari M.A. A Generalized Method
    for Estimation of Critical Constants. Fluid Phase Equilibria 147 (1998) 1-6
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
        x = exp(a[j] + b[j]*M + c[j]*Tb + d[j]*d20 + e[j]*Tb*d20) * \
            M**f[j] * Tb**(g[j]+h[j]*M) * d20**i[j]
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
    prop : dict
        A dict with the calculated properties:

            * Tf: Freezing temperature, [K]
            * Tb: Boiling temperature, [K]
            * SG: Specific gravity, [-]
            * d20: Density at 20ºC, [g/cm³]
            * Tc: Critic temperature, [K]
            * Pc: Critic pressure, [bar]
            * Vc: Critic volume, [cm³/g]
            * w: Acentric factor, [-]
            * I: Refractive index parameter, [-]
            * sigma: Surface tension, [dyn/cm]

    References
    ----------
    [27]_ Riazi, M.R., A1-Sahhaf, T.A. Physical Properties of Heavy Petroleum
    Fractions and Crude Oils. Fluid Phase Equilibria 117 (1996) 217-224.

    [28]_ Riazi, M.R., A1-Sahhaf, T. Physical Properties of n-Alkanes and
    n-Alkyl Hydrocarbons: Application to Petroleum Mixtures. Ind. Eng. Chem.
    Res. 34(11) (1995) 4145-4148
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
    [22]_ Hougen, O.A., Watson, K.M., Ragatz, R.A. Chemical Process Principles,
    2nd ed. New York: Wiley, 1959, p. 577.
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
    [23]_ Reid, R., Prausnitz, J.M., Sherwood, T. The Properties of Gases and
    Liquids, 3rd ed. New York: McGraw-Hill, 1977, p. 21.
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
    [21]_ Salerno, S., Cascella, M., May, D., Watson, P., Tassios, D.
    Prediction of Vapor Pressures and Saturated Volumes with a Simple Cubic
    Equation of State: Part I. A Reliable Data Base. Fluid Phase Equilibria 27
    (1986) 15-34
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
    [24]_ Nath, J. Acentric Factor and the Critical Volumes for Normal Fluids.
    Ind. Eng. Chem. Fundam. 21(3) (1985) 325-326
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
    [25]_ Lee, B.I., Kesler, M.G. A Generalized Thermodynamic Correlation Based
    on Three-Parameter Corresponding States. AIChE J. 21(3) (1975) 510-527
    """
    Zc = 0.2905 - 0.085*w                                               # Eq 21

    return unidades.Dimensionless(Zc)


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
    [51]_ Riazi, M.R., Daubert, T.E. Analytical Correlations
    Interconvert Distillation Curve Types. Oil & Gas Journal 84 (1986) 50-57

    [29]_ Riazi, M.R. Characterization and Properties of Petroleum Fractions.
    ASTM manual series MNL50, 2005
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
    [52]_ Daubert, T.E. Petroleum Fraction Distillation Interconversion.
    Hydrocarbon Processing 73(9) (1994) 75-78

    [20]_ API. Technical Data book: Petroleum Refining 6th Edition

    [29]_ Riazi, M.R. Characterization and Properties of Petroleum Fractions.
    ASTM manual series MNL50, 2005
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
    [51]_ Riazi, M.R., Daubert, T.E. Analytical Correlations
    Interconvert Distillation Curve Types. Oil & Gas Journal 84 (1986) 50-57

    [29]_ Riazi, M.R. Characterization and Properties of Petroleum Fractions.
    ASTM manual series MNL50, 2005
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
    [52]_ Daubert, T.E. Petroleum Fraction Distillation Interconversion.
    Hydrocarbon Processing 73(9) (1994) 75-78

    [20]_ API. Technical Data book: Petroleum Refining 6th Edition

    [29]_ Riazi, M.R. Characterization and Properties of Petroleum Fractions.
    ASTM manual series MNL50, 2005
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
    [51]_ Riazi, M.R., Daubert, T.E. Analytical Correlations
    Interconvert Distillation Curve Types. Oil & Gas Journal 84 (1986) 50-57

    [29]_ Riazi, M.R. Characterization and Properties of Petroleum Fractions.
    ASTM manual series MNL50, 2005
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
    [52]_ Daubert, T.E. Petroleum Fraction Distillation Interconversion.
    Hydrocarbon Processing 73(9) (1994) 75-78

    [20]_ API. Technical Data book: Petroleum Refining 6th Edition

    [29]_ Riazi, M.R. Characterization and Properties of Petroleum Fractions.
    ASTM manual series MNL50, 2005
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
    [53]_ Edmister, W.C., Okamoto, K.K. Applied Hydrocarbon Thermodynamics,
    Part 13: Equilibrium Flash Vaporization Correlations for Heavy Oils Under
    Subatmospheric Pressures. Petroleum Refiner 38(9) (1959) 271-288

    [29]_ Riazi, M.R. Characterization and Properties of Petroleum Fractions.
    ASTM manual series MNL50, 2005
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
    [20]_ API. Technical Data book: Petroleum Refining 6th Edition

    [29]_ Riazi, M.R. Characterization and Properties of Petroleum Fractions.
    ASTM manual series MNL50, 2005
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
    >>> parameter, r = curve_Predicted(x, T)
    >>> "%.0f %.4f %.4f" % tuple(parameter)
    '327 0.2028 1.3802'

    References
    ----------
    [54]_ Riazi, M.R. Distribution Model for Properties of Hydrocarbon-Plus
    Fractions. Ind. Eng. Chem. Res. 28(11) (1989) 1731-1735.
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
        Pour point, [ºR]

    Notes
    ------
    Raise :class:`NotImplementedError` if input isn't in limit:

        * 0.8 ≤ SG ≤ 1
        * 800ºR ≤ Tb ≤ 1500ºR
        * 2cSt ≤ v100 ≤ 960cSt

    Examples
    --------
    >>> T = unidades.Temperature(972, "R")
    >>> v = unidades.Diffusivity(3, "cSt")
    >>> "%0.0f" % PourPoint(0.839, T, v).R
    '458'

    References
    ----------
    [20]_ API. Technical Data book: Petroleum Refining 6th Edition
    """
    # Convert input Tb in Kelvin to Rankine to use in the correlation
    Tb_R = unidades.K2R(Tb)
    v100 = unidades.Diffusivity(v100).cSt

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

    Notes
    ------
    Raise :class:`NotImplementedError` if input isn't in limit:

        * 0.7 ≤ SG ≤ 1
        * 200ºF ≤ Tb ≤ 1100ºF

    Examples
    --------
    >>> T = unidades.Temperature(570.2, "F")
    >>> "%0.1f" % AnilinePoint(0.8304, T).R
    '635.5'

    References
    ----------
    [20]_ API. Technical Data book: Petroleum Refining 6th Edition
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

    Notes
    ------
    Raise :class:`NotImplementedError` if input isn't in limit:

        * 0.7 ≤ SG ≤ 0.86
        * 200ºF ≤ Tb ≤ 550ºF

    Examples
    --------
    >>> T = unidades.Temperature(414.5, "F")
    >>> "%0.1f" % SmokePoint(0.853, T).mm
    '16.7'

    References
    ----------
    [20]_ API. Technical Data book: Petroleum Refining 6th Edition
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

    Notes
    ------
    Raise :class:`NotImplementedError` if input isn't in limit:

        * 0.74 ≤ SG ≤ 0.9
        * 725ºR ≤ Tb ≤ 1130ºR

    Examples
    --------
    >>> T = unidades.Temperature(874.5, "R")
    >>> "%0.0f" % CloudPoint(0.799, T).R
    '417'

    References
    ----------
    [20]_ API. Technical Data book: Petroleum Refining 6th Edition
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

    Notes
    ------
    Raise :class:`NotImplementedError` if input isn't in limit:

        * 0.77 ≤ SG ≤ 0.93
        * 800ºR ≤ Tb ≤ 1225ºR

    Examples
    --------
    >>> T = unidades.Temperature(811.5, "R")
    >>> "%0.1f" % CloudPoint(0.787, T).R
    '383.4'

    References
    ----------
    [20]_ API. Technical Data book: Petroleum Refining 6th Edition
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

    Notes
    ------
    Raise :class:`NotImplementedError` if input isn't in limit:

        * 27 ≤ API ≤ 47
        * 360ºF ≤ Tb ≤ 700ºF

    Examples
    --------
    >>> T = unidades.Temperature(617, "F")
    >>> "%0.1f" % CetaneIndex(32.3, T)
    '57.1'

    References
    ----------
    [20]_ API. Technical Data book: Petroleum Refining 6th Edition
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


# PNA composition procedures
def PNA_Riazi(M, SG, n, d20=None, VGC=None, CH=None):
    """Calculate fractional compositon of paraffins, naphthenes and aromatics
    contained in petroleum fractions using the procedure 2B4.1 from API
    technical databook, pag. 225

    Parameters
    ------------
    M : float
        Molecular weight, [g/mol]
    SG : float
        Specific gravity, [-]
    n : float
        Refractive index, [-]
    d20 : float
        Density at 20ºC, [g/cm³]
    VGC : float
        Viscosity gravity constant, [-]
    CH : float
        Carbon/hydrogen ratio, [-]

    Returns
    -------
    xp : float
        Paraffins mole fraction, [-]
    xn : float
        Naphthenes mole fraction, [-]
    xa : float
        Aromatics mole fraction, [-]

    Notes
    -----
    Density at 20ºC is optional and can be calculated from specific gravity.
    VGC and CH are optional parameters, the procedure need one of them

    Raise :class:`NotImplementedError` if input viscosity or CH ratio are
    undefined

    Examples
    --------
    >>> "%0.3f %0.3f %0.3f" % PNA_Riazi(378, 0.9046, 1.5002, 0.9, VGC=0.8485)
    '0.606 0.275 0.118'

    References
    ----------
    [20]_ API. Technical Data book: Petroleum Refining 6th Edition
    """
    if d20 is None:
        d20 = SG-4.5e-3*(2.34-1.9*SG)
    Ri = n-d20/2.
    m = M*(n-1.475)
    if VGC is not None:
        if M <= 200:
            a, b, c = -13.359, 14.4591, -1.41344
            d, e, f = 23.9825, -23.333, 0.81517
        else:
            a, b, c = 2.5737, 1.0133, -3.573
            d, e, f = 2.464, -3.6701, 1.96312
        xp = a + b*Ri + c*VGC
        xn = d + e*Ri + f*VGC

    elif CH is not None:
        if M <= 200:
            xp = 2.57 - 2.877*SG + 0.02876*CH
            xn = 0.52641 - 0.7494*xp - 0.021811*m
        else:
            xp = 1.9842 - 0.27722*Ri - 0.15643*CH
            xn = 0.5977 - 0.761745*Ri + 0.068048*CH
    else:
        raise NotImplementedError()

    xa = 1-xp-xn
    return xp, xn, xa


def PNA_Peng_Robinson(Nc, M, WABP):
    """Calculate fractional compositon of paraffins, naphthenes and aromatics
    contained in petroleum fractions using the Peng-Robinson procedure

    Parameters
    ------------
    Nc : index
        Carbon number, [-]
    M : float
        Molecular weight, [g/mol]
    WABP : float
        Weight average boiling temperature, [K]

    Returns
    -------
    xp : float
        Paraffins mole fraction, [-]
    xn : float
        Naphthenes mole fraction, [-]
    xa : float
        Aromatics mole fraction, [-]

    Notes
    -----
    This method can too calculate the critical properties of fraction

    References
    ----------
    [56]_ Ahmed, T. Hydrocarbon Phase Behavior. Gulf Publishing, Houston, TX,
    1989.

    [58]_ Robinson, D.B., Peng, D.Y. The characterization of the heptanes and
    heavier fractions. Research Report 28. GPA, 1978. Tulsa, OK.
    """
    # Convert input Tb in Kelvin to Fahrenheit to use in the correlation
    WABP_R = unidades.K2R(WABP)

    Tbp = exp(log(1.8) + 5.8345183 + 0.84909035e-1*(Nc-6) -
              0.52635428e-2*(Nc-6)**2 + 0.21252908e-3*(Nc-6)**3 -
              0.44933363e-5*(Nc-6)**4 + 0.37285365e-7*(Nc-6)**5)
    Tbn = exp(log(1.8) + 5.8579332 + 0.79805995e-1*(Nc-6) -
              0.43098101e-2*(Nc-6)**2 + 0.14783123e-3*(Nc-6)**3 -
              0.27095216e-5*(Nc-6)**4 + 0.19907794e-7*(Nc-6)**5)
    Tba = exp(log(1.8) + 5.867176 + 0.80436947e-1*(Nc-6) -
              0.47136506e-2*(Nc-6)**2 + 0.18233365e-3*(Nc-6)**3 -
              0.38327239e-5*(Nc-6)**4 + 0.32550576e-7*(Nc-6)**5)

    Mp = 14.026*Nc + 2.016
    Mn = 14.026*Nc - 14.026
    Ma = 14.026*Nc - 20.074

    # Solve the system of equation
    a = [[1, 1, 1], [Tbp*Mp, Tbn*Mn, Tba*Ma], [Mp, Mn, Ma]]
    b = [1, M*WABP_R, M]
    Xp, Xn, Xa = solve(a, b)

    # Calculate the critical properties of fraction
    # pcp = (206.126096*Nc+29.67136)/(0.227*Nc+0.34)**2
    # pcn = (206.126096*Nc+206.126096)/(0.227*Nc+0.137)**2
    # pca = (206.126096*Nc+295.007504)/(0.227*Nc+0.325)**2
    # pc = Xp*pcp+Xn*pcn+Xa*pca

    # wp = 0.432*Nc-0.0457
    # wn = 0.0432*Nc-0.088
    # wa = 0.0445*Nc-0.0995

    # S = 0.996704+0.00043155*Nc
    # S1 = 0.99627245+0.00043155*Nc
    # tcp = S*(1+(3*log10(pcp)-3.501952)/7/(1+wp))*Tbp
    # tcn = S1*(1+(3*log10(pcn)-3.501952)/7/(1+wn))*Tbn
    # tca = S1*(1+(3*log10(pca)-3.501952)/7/(1+wa))*Tba
    # tc = Xp*tcp+Xn*tcn+Xa*tca
    return Xp, Xn, Xa


def PNA_Bergman(Tb, SG, Kw):
    """Calculate fractional compositon of paraffins, naphthenes and aromatics
    contained in petroleum fractions using the Bergman procedure

    Parameters
    ------------
    Tb : float
        Boiling temperature, [K]
    SG : float
        Specific gravity, [-]
    Kw : float
        Watson characterization factor, [-]

    Returns
    -------
    xp : float
        Paraffins mole fraction, [-]
    xn : float
        Naphthenes mole fraction, [-]
    xa : float
        Aromatics mole fraction, [-]

    Notes
    -----
    This method can too calculate the critical properties of fraction

    References
    ----------
    [56]_ Ahmed, T. Hydrocarbon Phase Behavior. Gulf Publishing, Houston, TX,
    1989.

    [57]_ Bergman, D.F., Tek, M.R., Katz, D.L. Retrograde Condensation in
    Natural Gas Pipelines. Project PR 2-29 of Pipelines Research Committee,
    AGA, January 1977
    """
    # Convert input Tb in Kelvin to Fahrenheit to use in the correlation
    Tb_F = unidades.K2F(Tb)

    # Estimate gthe weight fraction of the aromatic content in fraction
    Xwa = 8.47 - Kw

    # Solve the linear equation for weight fractions of paraffins and naphtenes
    gp = 0.582486 + 0.00069481*Tb_F - 0.7572818e-6*Tb_F**2 + \
        0.3207736e-9*Tb_F**3
    gn = 0.694208 + 0.0004909267*Tb_F - 0.659746e-6*Tb_F**2 + \
        0.330966e-9*Tb_F**3
    ga = 0.916103 - 0.000250418*Tb_F + 0.357967e-6*Tb_F**2 - \
        0.166318e-9*Tb_F**3
    a = [[1, 1], [1/gp, 1/gn]]
    b = [1-Xwa, 1/SG - Xwa/ga]
    Xwp, Xwn = solve(a, b)

    # Calculate the Tc, Pc and acentric factor for each cut
    # Tcp = 275.23 + 1.2061*Tb_F - 0.00032984*Tb_F**2
    # Pcp = 573.011 - 1.13707*Tb_F + 0.00131625*Tb_F**2 - 0.85103e-6*Tb_F**3
    # wp = 0.14 + 0.0009*Tb_F + 0.233e-6*Tb_F**2
    # Tcn = 156.8906 + 2.6077*Tb_F - 0.003801*Tb_F**2 + 0.2544e-5*Tb_F**3
    # Pcn = 726.414 - 1.3275*Tb_F + 0.9846e-3*Tb_F**2 - 0.45169e-6*Tb_F**3
    # wn = wp-0.075
    # Tca = 289.535 + 1.7017*Tb_F - 0.0015843*Tb_F**2 + 0.82358e-6*Tb_F**3
    # Pca = 1184.514 - 3.44681*Tb_F + 0.0045312*Tb_F**2 - 0.23416e-5*Tb_F**3
    # wa = wp-0.1
    # pc=Xwp*Pcp+Xwn*Pcn+Xwa*Pca
    # tc=Xwp*Tcp+Xwn*Tcn+Xwa*Tca
    # w=Xwp*wp+Xwn*wn+Xwa*wa

    return Xwp, Xwn, Xwa


def PNA_van_Nes(M, n, d20, S):
    """Calculate fractional compositon of paraffins, naphthenes and aromatics
    contained in petroleum fractions using the Van Nes procedure

    Parameters
    ------------
    M : float
        Molecular weight, [g/mol]
    n : float
        Refractive index, [-]
    d20 : float
        Density at 20ºC, [g/cm³]
    S : float
        Sulfur content, [-]

    Returns
    -------
    xp : float
        Paraffins mole fraction, [-]
    xn : float
        Naphthenes mole fraction, [-]
    xa : float
        Aromatics mole fraction, [-]

    Notes
    -----
    This method can too calculate the total number of rings, the number of
    aromatics rings and the naphthenic rings

    References
    ----------
    [55]_ Van Nes, K., Van Western, H.A. Aspects of the Constitution of Mineral
    Oils. Elsevier, New York, 1951

    [29]_ Riazi, M.R. Characterization and Properties of Petroleum
    Fractions. ASTM manual series MNL50, 2005
    """
    v = 2.51*(n-1.475)-(d20-0.851)
    w = (d20-0.851)-1.11*(n-1.475)
    if v > 0:
        a = 430
    else:
        a = 670
    if w > 0:
        CR = 820*w-3*S + 10000./M
    else:
        CR = 1440*w - 3*S + 10600./M
    Ca = a*v + 3660/M
    Cn = CR - Ca
    Cp = 100 - CR
    return Cp/100., Cn/100., Ca/100.


# Viscosity correlation
def viscoAPI(Tb=None, Kw=None, v100=None, v210=None, T=None, T1=None, T2=None):
    """Calculate viscosity of a petroleum fraction using the API procedure
    11A4.2 and 11A4.4

    Parameters
    ------------
    Tb : float
        Boiling temperature, [K]
    Kw : float
        Watson characterization factor, [-]
    v100 : float, optional
        Kinematic viscosity at 100ºF, [m²/s]
    v210 : float, optional
        Kinematic viscosity at 210ºF, [m²/s]
    T : float, optional
        Optional temperature to calculate the viscosity, [K]
    T1 : float, optional
        Temperature for v100 viscosity value, [K]
    T2 : float, optional
        Temperature for v210 viscosity value, [K]

    Returns
    -------
    v100 : float
        Kinematic viscosity at 100ºF, [m²/s]
    v210 : float
        Kinematic viscosity at 210ºF, [m²/s]
    vT : float, optional
        Kinematic viscosity at T input, [m²/s]

    Notes
    -----
    By default the method calculate the viscosity at 100ºF and 210ºF, can too
    calculate the viscosity at input T temperature
    if T1 and T2 are different to 100ºF and 210ºF v100 and v210 input are at
    that temperatures
    if both viscosity are input Tb and Kw input aren't unnecesary

    Examples
    --------
    Example A in procedure 11A4.2

    >>> Tb = unidades.Temperature(598.7, "F")
    >>> T = unidades.Temperature(140, "F")
    >>> mu = viscoAPI(Tb, 11.87, T=T)
    >>> "%0.2f %0.2f %0.2f" % (mu["v100"].cSt, mu["v210"].cSt, mu["vT"].cSt)
    '5.94 1.88 3.57'

    Example B in procedure 11A4.2

    The values differ, in reference only use 2 decimal number in intermediate
    >>> Tb = unidades.Temperature(647.3, "F")
    >>> T = unidades.Temperature(304, "F")
    >>> mu = viscoAPI(Tb, 9.955, T=T)
    >>> "%0.2f %0.2f %0.2f" % (mu["v100"].cSt, mu["v210"].cSt, mu["vT"].cSt)
    '39.41 5.20 2.12'

    Example C in procedure 11A4.2

    >>> Tb = unidades.Temperature(382.1, "F")
    >>> mu = viscoAPI(Tb, 11.93)
    >>> "%0.3f" % mu["v100"].cSt
    '1.249'

    >>> Tb = unidades.Temperature(1113.2, "F")
    >>> mu = viscoAPI(Tb, 11.94)
    >>> "%0.0f" % mu["v100"].cSt
    '2612'

    Example B in procedure 11A4.4

    >>> T = unidades.Temperature(68, "F")
    >>> T1 = unidades.Temperature(32, "F")
    >>> T2 = unidades.Temperature(104, "F")
    >>> v100 = unidades.Diffusivity(1.644, "cSt")
    >>> v210 = unidades.Diffusivity(0.925, "cSt")
    >>> mu = viscoAPI(T1=T1, T2=T2, T=T, v100=v100, v210=v210)
    >>> "%0.3f" % mu["vT"].cSt
    '1.201'

    References
    ----------
    [20]_ API. Technical Data book: Petroleum Refining 6th Edition
    """
    if v100 is None or v210 is None:
        # Convert input Tb in Kelvin to Rankine to use in the correlation
        Tb_R = unidades.K2R(Tb)

    if v100 is None:
        # Calculate the viscosity at 100ºF
        mu_ref = 10**(-1.35579 + 8.16059e-4*Tb_R + 8.38505e-7*Tb_R**2)
        A1 = 34.931-8.84387e-2*Tb_R+6.73513e-5*Tb_R**2-1.01394e-8*Tb_R**3
        A2 = -2.92649+6.98405e-3*Tb_R-5.09947e-6*Tb_R**2+7.49378e-10*Tb_R**3
        mu_cor = 10**(A1+A2*Kw)
        v100 = unidades.Diffusivity(mu_ref+mu_cor, "cSt")
    else:
        v100 = unidades.Diffusivity(v100)

    if v210 is None:
        # Calculate the viscosity at 210ºF
        mu = 10**(-1.92353+2.41071e-4*Tb_R+0.5113*log10(Tb_R*v100.cSt))
        v210 = unidades.Diffusivity(mu, "cSt")
    else:
        v210 = unidades.Diffusivity(v210)

    vs = {}
    vs["v100"] = v100
    vs["v210"] = v210

    if T is not None:
        # Calculate the viscosity at other temperature
        T_R = unidades.K2R(T)
        if T1 is None:
            T1 = unidades.Temperature(100, "F")
        T1_R = T1.R
        if T2 is None:
            T2 = unidades.Temperature(210, "F")
        T2_R = T2.R

        Z1 = v100.cSt+0.7+exp(-1.47-1.84*v100.cSt-0.51*v100.cSt**2)
        Z2 = v210.cSt+0.7+exp(-1.47-1.84*v210.cSt-0.51*v210.cSt**2)
        B = (log10(log10(Z1))-log10(log10(Z2)))/(log10(T1_R)-log10(T2_R))
        Z = 10**(10**(log10(log10(Z1))+B*(log10(T_R)-log10(T1_R))))
        mu = Z - 0.7 - exp(
            -0.7487-3.295*(Z-0.7)+0.6119*(Z-0.7)**2-0.3191*(Z-0.7)**3)
        vs["vT"] = unidades.Diffusivity(mu, "cSt")

    return vs


def SUS(T, v):
    """Calculate the Saybolt Universal Viscosity in Saybolt Universal Seconds
    (SUS) from kinematic viscosity, also referenced in API Procedure 11A1.1,
    pag 1027

    Parameters
    ------------
    T : float
        Temperature, [K]
    v : float
        Kinematic viscosity, [cSt]

    Returns
    -------
    SUS : float
        Saybolt Universal Seconds, [s]

    Examples
    --------
    Example A from [20]_, SUS at 200ºF for ν=53cSt

    >>> T = unidades.Temperature(617, "F")
    >>> "%0.0f" % SUS(T, 53)
    '254'

    Example B from [20]_, SUS at 0ºF for ν=90cSt

    >>> T = unidades.Temperature(0, "F")
    >>> "%0.0f" % SUS(T, 90)
    '415'

    References
    ----------
    [20]_ API. Technical Data book: Petroleum Refining 6th Edition

    [30]_ ASTM D2161-05. Standard Practice for Conversion of Kinematic
    Viscosity to Saybolt Universal Viscosity or to Saybolt Furol Viscosity.
    ASTM International, West Conshohocken, PA, 2005, www.astm.org
    """
    # Convert input temperature to Fahrenheit
    t_F = unidades.K2F(T)

    # Eq 5
    U100 = 4.6324*v + (1+0.03264*v)/(3930.2+262.7*v+23.97*v**2+1.646*v**3)*1e5

    # Eq 6
    Ut = (1+0.000061*(t_F-100))*U100
    return unidades.Time(Ut)


def SFS(T, v):
    """Calculate the Saybolt Furol Viscosity in Saybolt Furol Seconds
    (SFS) from kinematic viscosity, also referenced in API Procedure 11A1.4,
    pag 1031

    Parameters
    ------------
    T : float
        Temperature, [K]
    v : float
        Kinematic viscosity, [cSt]

    Returns
    -------
    SFS : float
        Saybolt Furol Seconds, [s]

    Examples
    --------
    Example A from [20]_, SFS at 122ºF and 210ºF for ν122=300cSt, v210=100cSt

    >>> T = unidades.Temperature(122, "F")
    >>> "%0.1f" % SFS(T, 300)
    '141.7'
    >>> T = unidades.Temperature(210, "F")
    >>> "%0.1f" % SFS(T, 100)
    '48.4'

    References
    ----------
    [20]_ API. Technical Data book: Petroleum Refining 6th Edition

    [30]_ ASTM D2161-05. Standard Practice for Conversion of Kinematic
    Viscosity to Saybolt Universal Viscosity or to Saybolt Furol Viscosity.
    ASTM International, West Conshohocken, PA, 2005, www.astm.org
    """
    if T == 323.15:
        S = 0.4717*v+13924/(v**2-72.59*v+6816)                         # Eq 7
    else:
        S = 0.4792*v+5610/(v**2+2130)                                   # Eq 8
    return unidades.Time(S)


def MuL_Singh(T, v100):
    """Calculate kinematic viscosity of liquid petroleum fractions at low
    pressure by Singh correlation, also referenced in API Procedure 11A4.1,
    pag 1031

    Parameters
    ------------
    T : float
        Temperature, [K]
    v100 : float
        Kinematic viscosity at 100ºF, [cSt]

    Returns
    -------
    v : float
        Kinematic viscosity at T, [cSt]

    Examples
    --------
    Example from [20]_, Sumatran crude at 210ºF

    >>> T = unidades.Temperature(210, "F")
    >>> "%0.3f" % MuL_Singh(T, 1.38).cSt
    '0.705'

    References
    ----------
    [20]_ API. Technical Data book: Petroleum Refining 6th Edition

    [31]_ Singh, B., Mutyala, S.R., Puttagunta, V. Viscosity range from one
    test. Hydrocarbon Processing 69 (1990) 39-41
    """
    # Convert input temperature to Rankine
    t_R = unidades.K2R(T)
    S = 0.28008*log10(v100) + 1.8616
    B = log10(v100) + 0.86960
    mu = 10**(B*(559.67/t_R)**S - 0.86960)
    return unidades.Diffusivity(mu, "cSt")


# Component predition
def H_Riazi(S, CH):
    """Calculate hydrogen content of petroleum fractions

    Parameters
    ------------
    S : float
        Sulfur content, [%]
    CH : float
        Carbon/hydrogen ratio, [-]

    Returns
    -------
    H : float
        Hydrogen content, [%]

    Notes
    -----
    A Simple mass balance discarting other elements composition than hydrogen,
    carbon and sulfur

    References
    ----------
    [29]_ Riazi, M.R. Characterization and Properties of Petroleum Fractions.
    ASTM manual series MNL50, 2005
    """
    H = (100-S)/(1+CH)
    return H


def H_Goossens(M, n, d20):
    """Calculate hydrogen content of petroleum fractions using the correlation
    of Goossens

    Parameters
    ------------
    M : float
        Molecular weight, [-]
    n : float
        Refractive index, [-]
    d20 : float
        Liquid density at 20ºC, [g/cm³]

    Returns
    -------
    H : float
        Hydrogen content, [%]

    Examples
    --------
    Table 1 data for n-decane

    >>> "%0.2f" % H_Goossens(142, 1.4119, 0.7301)
    '15.45'

    References
    ----------
    [60]_ Goossens, A.G. Prediction of the Hydrogen Content of Petroleum
    Fractions. Ind. Eng. Chem. Res. 36(6) (1997) 2500-2504
    """
    H = 30.346 + (82.952-65.341*n)/d20 - 306/M                           # Eq 3
    return H


def H_ASTM(Tb, SG, xa):
    """Calculate hydrogen content of petroleum fractions

    Parameters
    ------------
    Tb : float
        Mean average boiling point, [ºR]
    SG : float
        Specific gravity, [-]
    xa : float
        Aromatics mole fraction, [-]

    Returns
    -------
    H : float
        Hydrogen content, [%]

    References
    ----------
    [29]_ Riazi, M.R. Characterization and Properties of Petroleum Fractions.
    ASTM manual series MNL50, 2005

    [62]_ ASTM, Annual Book of Standards, ASTM International, West
    Conshohocken, PA, 2002
    """

    H = (5.2407 + 0.01448*Tb - 7.018*xa)/SG - 0.901*xa + 0.01298*xa*Tb - \
        0.01345*Tb + 5.6879
    return H


def H_Jenkins_Walsh(SG, anilineP):
    """Calculate hydrogen content of petroleum fractions using the correlation
    of Jenkins-Walsh

    Parameters
    ------------
    SG : float
        Specific gravity, [-]
    anilineP : float
        Aniline point, [K]

    Returns
    -------
    H : float
        Hydrogen content, [%]

    References
    ----------
    [61]_ Jenkins, G.I., Walsh, R.E. Quick Measure of Jet Fuel Properties.
    Hydrocarbon Processing 47(5) (1968) 161-164
    """
    H = 11.17 - 12.89*SG + 0.0389*anilineP
    return H


def S_Riazi(M, SG, Ri, m):
    """Calculate sulfur content of petroleum fractions

    Parameters
    ------------
    M : float
        Molecular weight, [-]
    SG : float
        Specific gravity, [-]
    Ri : float
        Refractivity intercept, [-]
    m : float
        characterization factor, [-]

    Returns
    -------
    S : float
        Sulfur content, [%]

    Examples
    --------
    S predicted in table 3 for a crude nº1

    >>> "%0.2f" % S_Riazi(202, 0.837, 1.0511, -1.4849)
    '1.19'

    crude nº61

    >>> "%0.2f" % S_Riazi(1484, 1.0209, 1.1611, 289.2786)
    '2.85'

    References
    ----------
    [59]_ Riazi, M.R., Nasimi, N., Roomi, Y. Estimating Sulfur Content of
    Petroleum Products and Crude Oils. Ind. Eng. Chem. Res. 38(11) (1999)
    4507-4512
    """
    if M < 200:
        S = 177.448 - 170.946*Ri + 0.2258*m + 4.054*SG
    else:
        S = -58.02 + 38.463*Ri - 0.023*m + 22.4*SG
    return S


def CombustionHeat(API, water=0, ash=0, S=0):
    """Calculate gross and net heat of combustion at 60ºF for petroleum
    fractions, referenced in API procedure 14A1.3, pag 1236

    Parameters
    ------------
    API : float
        API Specific gravity, [-]
    water : float
        Weight fraction of water, [%]
    ash : float
        Weight fraction of inerts, [%]
    S : float
        Weight fraction of inerts, [%]

    Returns
    -------
    Hgross : float
        Gross heat of combustion, [Btu/lb]
    Hnet : float
        Net heat of combustion, [Btu/lb]

    Examples
    --------
    Fuel oil with 11.3ºAPI, 1.49%S, 1.67%ash, 0.3%water

    >>> hg, hn = CombustionHeat(11.3, 0.3, 1.67, 1.49)
    >>> "%0.0f %0.0f" % (hg.Btulb, hn.Btulb)
    '18218 17222'

    References
    ----------
    [20]_ API. Technical Data book: Petroleum Refining 6th Edition
    """
    # Interpolation data from Table 14-0.1 for normal sulfur and inert content
    if API < 35:
        _api = [0, 5, 10, 15, 20, 15, 30, 35]
        _S = [2.95, 2.35, 1.8, 1.35, 1., 0.7, 0.4, 0.3]
        _ash = [1.15, 1., 0.95, 0.85, 0.75, 0.7, 0.65, 0.6]

        S -= interp1d(_api, _S)(API)
        ash -= interp1d(_api, _ash)(API)

    # Gross heating value
    hg = 17672 + 66.6*API - 0.316*API**2 - 0.0014*API**3
    hgc = unidades.Enthalpy(hg - 0.01*hg*(water+S+ash) + 40.5*S, "Btulb")

    # Net heating value
    hn = 16796 + 54.5*API - 0.217*API**2 - 0.0019*API**3
    hnc = hn - 0.01*hn*(water+S+ash) + 40.5*S - 10.53*water
    hnc = unidades.Enthalpy(hnc, "Btulb")

    return hgc, hnc


class Petroleo(newComponente):
    """Class to define a heavy oil fraction with a unknown composition

    Parameters
    ------------
    M : float
        Molecular weight, [-]
    Tb : float
        Mean average boiling point, [K]
    SG : float
        Specific gravity, [-]
    API : float
        API gravity, [-]
    CH : float
        Carbon/hydrogen ratio, [-]
    I : float
        Huang parameter, [-]
    n : float
        Refractive index, [-]
    Kw : float
        Watson characterization factor, [-]
    v100 : float
        Kinematic viscosity at 100ºF, [m²/s]
    v210 : float
        Kinematic viscosity at 210ºF, [m²/s]
    H : float
        Hydrogen content, [-]
    S : float
        Sulfur content, [-]
    N : float
        Nitrogen content, [-]
    Nc : float
        Fraction carbon number, [-]

    curveType : string
        Name of curve data in curve input: D86, TBP, SD, D1186, EFV
    T_curve : list
        Boiling point temperature, [K]
    X_curve : list
        Volume or weight distillation fraction array
    P_curve : float
        Distillation data pressure, [Pa]
    fit_curve : array
        Array of fit parameter of curve distillation [To, A, B]

    Example
    -------
    Example from [20]_, procedure 2B1.1

    >>> X = [0.1, 0.3, 0.5, 0.7, 0.9]
    >>> D86 = [149, 230, 282, 325, 371]
    >>> D86_K = [unidades.F2K(t) for t in D86]
    >>> oil = Petroleo(curveType="D86", X_curve=X, T_curve=D86_K)
    >>> "%.0f %.0f %.0f %.0f %.0f" % (\
        oil.VABP.F, oil.MABP.F, oil.WABP.F, oil.CABP.F, oil.MeABP.F)
    '271 241 279 264 253'

    Notes
    -----
    The calculated instance has the necessary calculated properties to be used
    in the flowsheet as hypotethical component.

    The definition can be with several input parameters:

        * Distillation data (TBP or other). This is the prefered and accurate
          method to define the pseudocomponent.
        * Any combination of physical properties (Tb, SG, M, I, n, CH, v100)
        * Nc as only paramter to calculata all the other propeties, this is the
          less accurate method and it's not recommended.

    Refractive index is equivalent as input parameter to Huang parameter.
    Watson characterization factor is equivalent to SG or boiling temperature
    API gravity is equivalent to SG as input paramter
    """
    kwargs = {"M": 0.0,
              "Tb": 0.0,
              "SG": 0.0,
              "API": 0.0,
              "CH": 0.0,
              "n": 0.0,
              "I": 0.0,
              "Kw": 0.0,
              "v100": 0.0,
              "v210": 0.0,
              "H": 0.0,
              "S": 0.0,
              "N": 0.0,
              "Nc": 0,
              "name": "",

              "curveType": "",
              "T_curve": [],
              "X_curve": [],
              "P_curve": 101325,
              "fit_curve": []}

    status = 0
    _bool = False
    msg = ""
    definicion = 0
    calculatedMethod = {}

    CURVE_TYPE = ["D86", "TBP", "SD", "D1186", "EFV"]
    METHODS_PNA = ["Riazi (API)",
                   "Peng Robinson (1978)",
                   "Bergman (1977)",
                   "Van Nes (1951)"]

    METHODS_M = ["Riazi-Daubert (1987)",
                 "Riazi-Daubert (1980)",
                 "Lee-Kesler (1976)",
                 "Sim Daubert (1980)",
                 "Goossens (1971)",
                 "Twu (1983)"]
    METHODS_Tb = ["Riazi-Daubert (1987)",
                  "Riazi (Heavy fractions)",
                  "Rowe (1978)",
                  "Sancet (2007)",
                  "Silva-Rodríguez (1992)",
                  "Soreide (1989)"]
    METHODS_SG = ["Riazi-Daubert (1987)", "Silva-Rodríguez (1992)"]
    METHODS_w = ["Edmister (1958)",
                 "Lee-Kesler (1976)",
                 "Korsten (2000)",
                 "Watansiri-Owens-Starling (1985)",
                 "Magoulas-Tassios (1990)"]
    METHODS_crit = ["Riazi-Daubert (1987)",
                    "Riazi-Daubert (1980)",
                    "Cavett (1962)",
                    "Lee-Kesler (1976)",
                    "Sim Daubert (1980)",
                    "Watansiri-Owens-Starling (1985)",
                    "Rowe (1978)",
                    "Standing (1977)",
                    "Willman-Teja (1987)",
                    "Magoulas-Tassios (1990)",
                    "Tsonopoulos (1986)",
                    "Twu (1983)",
                    "Sancet (2007)",
                    "Riazi (Heavy fractions)",
                    "Riazi-Alsahhaf (1998)"]
    METHODS_Vc = ["Riazi-Daubert (1987)",
                  "Riazi-Daubert (1980)",
                  "Sim Daubert (1980)",
                  "Watansiri-Owens-Starling (1985)",
                  "Twu (1983)",
                  "Hall-Yarborough (1971)",
                  "Riazi (Heavy fractions)",
                  "Riazi-Alsahhaf (1998)"]
    METHODS_n = ["Riazi-Daubert (1987)",
                 "Riazi (Heavy fractions)",
                 "Willman-Teja (1987)"]
    METHODS_H = ["Riazi", "Goossens", "ASTM", "Jenkins-Walsh"]
    METHODS_Zc = ["Lee-Kesler (1975)",
                  "Salerno (1986)",
                  "Nath (1985)",
                  "Reid (1977)",
                  "Hougen (1952)"]

    def __init__(self, **kwargs):
        self.Preferences = ConfigParser()
        self.Preferences.read(conf_dir+"pychemqtrc")
        self.__call__(**kwargs)

    def __call__(self, **kwargs):
        self.kwargs.update(kwargs)
        if kwargs:
            self._bool = True
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
        self.status = 0
        self.msg = QApplication.translate("pychemqt", "Insufficient input")

        self.hasSG = self.kwargs["SG"] or self.kwargs["API"] or \
            (self.kwargs["Kw"] and self.kwargs["Tb"])
        self.hasRefraction = self.kwargs["n"] or self.kwargs["I"]
        self.hasCurve = self.kwargs["curveType"] and self.kwargs["T_curve"] \
            and self.kwargs["X_curve"]

        # Tipo Definición
        if self.hasCurve or (self.kwargs["Tb"] and self.hasSG):
            self.definicion = 1
        elif self.kwargs["M"] and self.hasSG:
            self.definicion = 2
        elif self.kwargs["Tb"] and self.hasRefraction:
            self.definicion = 3
        elif self.kwargs["M"] and self.hasRefraction:
            self.definicion = 4
        elif self.kwargs["Tb"] and self.kwargs["CH"]:
            self.definicion = 5
        elif self.kwargs["M"] and self.kwargs["CH"]:
            self.definicion = 6
        elif self.kwargs["v100"] and self.hasRefraction:
            self.definicion = 7
        elif self.kwargs["Nc"]:
            self.definicion = 8

        if self.definicion:
            self.status = 1
            self.msg = ""
            return True

    def calculo(self):
        self.formula = ""
        self.cp = []
        self.Vliq = 0
        self.rackett = 0
        self.Tf = 0
        self.Hf = 0
        self.Gf = 0

        if self.hasCurve:
            # Create the complete curve, using fitting procedure o using the
            # input kwargs parameters
            # The curva variable has the values at normaliced volume/weight
            # fractions
            # 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 99
            # 0, 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12
            if self.kwargs["fit_curve"]:
                parCurva = self.kwargs["fit_curve"]
            else:
                x = self.kwargs["X_curve"]
                T = self.kwargs["T_curve"]
                parCurva, r = curve_Predicted(x, T)

            curva = []
            for x in [0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 99]:
                if x/100 in self.kwargs["X_curve"]:
                    index = self.kwargs["X_curve"].index(x/100)
                    T = self.kwargs["T_curve"][index]
                else:
                    T = _Tb_Predicted(parCurva, x/100)
                curva.append(unidades.Temperature(T))

            # If SG isn't known, use the correlation proposed in Riazi, [29]_
            # Eq. 3.17, Pag 102
            if self.kwargs["curveType"] == "D86":
                a, b, c = 0.08342, 0.10731, 0.26288
            elif self.kwargs["curveType"] == "TBP":
                a, b, c = 0.10431, 0.12550, 0.26288
            elif self.kwargs["curveType"] == "EFV":
                a, b, c = 0.09138, -0.0153, 0.36844
            else:
                pass
            self.SG = a*curva[2]**b*curva[6]**c

            # Calculate the distillation curve missing
            P = unidades.Pressure(10, "mmHg")
            Pcurve = self.kwargs["P_curve"]
            X = self.kwargs["X_curve"]

            if self.kwargs["curveType"] == "D86":
                self.D86 = curva
                self.TBP = self._D86_TBP(curva, X)
                self.EFV = D86_EFV(curva, X, self.SG)
                self.SD = None
                TBP_10mHg = [Tb_Pressure(t, P, reverse=True) for t in self.TBP]
                D1160_10mHg = D1160_TBP_10mmHg(TBP_10mHg, reverse=True)
                self.D1160 = [Tb_Pressure(t, P) for t in D1160_10mHg]

            elif self.kwargs["TBP"]:
                if Pcurve == 101325:
                    self.TBP = curva
                    TBP_10mmHg = [Tb_Pressure(t, P, reverse=True)
                                  for t in self.TBP]
                else:
                    self.TBP = [Tb_Pressure(t, Pcurve) for t in curva]
                    if Pcurve == P:
                        TBP_10mmHg = curva
                    else:
                        TBP_10mmHg = [Tb_Pressure(t, P, reverse=True)
                                      for t in self.TBP]
                self.D86 = self._D86_TBP(curva, X, reverse=True)
                self.EFV = D86_EFV(self.D86, X, self.SG)
                self.SD = None
                D1160_10mmHg = D1160_TBP_10mmHg(TBP_10mmHg, reverse=True)
                self.D1160 = [Tb_Pressure(t, 10./760) for t in D1160_10mmHg]

            elif self.kwargs["SD"]:
                self.SD = curva
                self.D86 = self._SD_D86(curva, X)
                self.TBP = SD_TBP(curva, X)
                self.EFV = D86_EFV(self.D86, X, self.SG)
                TBP_10mmHg = [Tb_Pressure(t, P, True) for t in self.TBP]
                D1160_10mmHg = D1160_TBP_10mmHg(TBP_10mmHg, reverse=True)
                self.D1160 = [Tb_Pressure(t, P) for t in D1160_10mmHg]

            elif self.kwargs["EFV"]:
                if Pcurve == 101325:
                    self.EFV = curva
                else:
                    self.EFV = [Tb_Pressure(t, Pcurve) for t in curva]
                self.D86 = D86_EFV(curva, X, self.SG, reverse=True)
                self.SD = None
                self.TBP = self._D86_TBP(self.D86, X)
                TBP_10mmHg = [Tb_Pressure(t, P, True) for t in self.TBP]
                D1160_10mmHg = D1160_TBP_10mmHg(TBP_10mmHg, reverse=True)
                self.D1160 = [Tb_Pressure(t, P) for t in D1160_10mmHg]

            else:
                if Pcurve == P:
                    D1160_10mmHg = curva
                    self.D1160 = [Tb_Pressure(t, P) for t in D1160_10mmHg]
                elif Pcurve == 101325:
                    self.D1160 = curva
                    D1160_10mmHg = [Tb_Pressure(t, P, True) for t in curva]
                else:
                    self.D1160 = [Tb_Pressure(t, Pcurve) for t in curva]
                    D1160_10mmHg = [Tb_Pressure(t, P, reverse=True)
                                    for t in self.D1160]
                TBP_10mmHg = D1160_TBP_10mmHg(D1160_10mmHg)
                self.TBP = [Tb_Pressure(t, P) for t in TBP_10mmHg]
                self.D86 = self._D86_TBP(self.TBP, X, reverse=True)
                self.EFV = D86_EFV(self.D86, X, self.SG)
                self.SD = None

            # Hardcoded key temperatures in distillation curve
            self.T5 = self.D86[1]
            self.T10 = self.D86[2]
            self.T20 = self.D86[3]
            self.T30 = self.D86[4]
            self.T40 = self.D86[5]
            self.T50 = self.D86[6]
            self.T60 = self.D86[7]
            self.T70 = self.D86[8]
            self.T80 = self.D86[9]
            self.T90 = self.D86[10]

            # Calculate average boiling points
            # Ref. [20] Procdedure 2B1.1
            self.VABP = unidades.Temperature((
                self.T10 + self.T30 + self.T50 + self.T70 + self.T90)/5.)
            SL = (self.T90.F-self.T10.F)/80.
            DT1 = -3.062123-0.01829*(self.VABP.F-32)**0.6667+4.45818*SL**0.25
            DT2 = -0.56379-0.007981*(self.VABP.F-32)**0.6667+3.04729*SL**0.333
            DT3 = -0.23589-0.06906*(self.VABP.F-32)**0.45+1.8858*SL**0.45
            DT4 = -0.94402-0.00865*(self.VABP.F-32)**0.6667+2.99791*SL**0.333
            self.WABP = unidades.Temperature(self.VABP.F+exp(DT1), "F")
            self.MABP = unidades.Temperature(self.VABP.F-exp(DT2), "F")
            self.CABP = unidades.Temperature(self.VABP.F-exp(DT3), "F")
            self.MeABP = unidades.Temperature(self.VABP.F-exp(DT4), "F")
            self.Tb = self.MeABP

            # Reid Vapour Pressure, [29]_ Eq 3.100
            ReidVP = 3.3922 - 0.02537*self.T5.C - 0.070739*self.T10.C + \
                .00917*self.T30.C - .0393*self.T50.C + 6.8257e-4*self.T10.C**2
            self.ReidVP = unidades.Pressure(ReidVP, "bar")

            # V/L Ratio, [29]_ Eq 106
            self.VL_12 = 88.5 - 0.19*self.T70 - 42.5*self.ReidVP.bar
            self.VL_20 = 90.6 - 0.25*self.T70 - 39.2*self.ReidVP.bar
            self.VL_36 = 94.7 - 0.36*self.T70 - 32.3*self.ReidVP.bar

            # Critical vapor locking index, [29]_ Eq 109
            self.CVLI = 4.27 + 0.24*self.T70 + 0.069*self.ReidVP.bar

            # Fuel volatility index, [29]_ Eq 110
            self.FVI = 1000*self.ReidVP.bar + 7*self.T70

            # Flash point, API procedure 2B7.1
            # Pensky-Martens Closed Cup (ASTM D93)
            self.FlashPc = unidades.Temperature(0.69*self.T10.F-118.2, "F")
            # Cleveland Open Cup (ASTM D92)
            self.FlashPo = unidades.Temperature(0.68*self.T10.F-109.6, "F")

            # Calculate PNA decomposition
            # self.xp, self.xn, self.xa = self._PNA()

        if self.definicion == 8:
            # Simple definition known only the number of carbon atoms
            prop = prop_Ahmed(self.kwargs["Nc"])
            self.M = prop["M"]
            self.Tc = prop["Tc"]
            self.Pc = prop["Pc"]
            self.Vc = prop["Vc"]
            self.Tb = prop["Tb"]
            self.f_acent = prop["w"]
            self.SG = prop["SG"]
            self.API = 141.5/self.SG-131.5
            self.d20 = self.SG-4.5e-3*(2.34-1.9*self.SG)

        else:
            # Definition with any pair of variables
            # The curve definition is as SG + Tb with aditional properties
            # calculated above

            # Load input parameters
            if self.kwargs["Tb"]:
                self.Tb = unidades.Temperature(self.kwargs["Tb"])
            if self.kwargs["M"]:
                self.M = self.kwargs["M"]

            if self.hasSG:
                if self.kwargs["SG"]:
                    self.SG = self.kwargs["SG"]
                elif self.kwargs["API"]:
                    self.SG = 141.5/(self.kwargs["API"]+131.5)
                elif self.kwargs["Kw"] and self.kwarg["Tb"]:
                    self.SG = self.Tb.R**(1./3)/self.kwargs["Kw"]
                self.API = 141.5/self.SG-131.5
                self.d20 = self.SG-4.5e-3*(2.34-1.9*self.SG)

            if self.hasRefraction:
                if self.kwargs["n"]:
                    self.n = self.kwargs["n"]
                else:
                    I = self.kwargs["I"]
                    self.n = ((1+2*I)/(1-I))**0.5
                self.I = (self.n**2-1)/(self.n**2+2)

            if self.kwargs["CH"]:
                self.CH = self.kwargs["CH"]

            if self.kwargs["v100"]:
                self.v100 = unidades.Diffusivity(self.kwargs["v100"])
            if self.kwargs["v210"]:
                self.v210 = unidades.Diffusivity(self.kwargs["v210"])

        # Calculate the unknown variables
        if not self.kwargs["M"]:
            self.M = self._M()
        if not self.kwargs["Tb"]:
            self.Tb = self._Tb()

        if self.kwargs["Kw"]:
            self.Kw = self.kwargs["Kw"]
        else:
            self.Kw = self.Tb.R**(1./3)/self.SG

        if not self.hasSG:
            self.SG = self._SG()
            self.API = 141.5/self.SG-131.5

        if not self.hasRefraction:
            self.n, self.I = self._n()

        if not self.kwargs["CH"]:
            self.CH = self._CH()

        self.Tc, self.Pc = self._critic()
        self.Vc = self._Vc()
        self.f_acent = self._f_acent()
        self.Zc = self._Zc()

        # self.Parametro_solubilidad = prop_Riazi_Alsahhaf(9, self.M), "calcc"
        # Tc, Pc, Vc, w, Dm prop_Watansiri_Owens_Starling(Tb, SG, M)
        # Tc, Pc, Vc, Tb, d20, I prop_Riazi(SG, tita, val)

        if not self.kwargs["v100"] and not self.kwargs["v210"]:
            visco = viscoAPI(self.Tb, self.Kw)
            self.v100 = visco["v100"]
            self.v210 = visco["v210"]
        elif not self.kwargs["v210"]:
            visco = viscoAPI(self.Tb, self.Kw, v100=self.kwargs["v100"])
            self.v210 = visco["v210"]
        else:
            visco = viscoAPI(self.Tb, self.Kw)
            self.v100 = visco["v100"]

        # Calculate Viscosity gravity constant (VGC)
        if self.kwargs["v210"] and not self.kwargs["v100"]:
            T = unidades.Temperature(210, "F")
            vSUS = SUS(T, self.v210.cSt)
            self.VGC = (self.SG-0.24-0.022*log10(vSUS-35.5))
        else:
            T = unidades.Temperature(100, "F")
            vSUS = SUS(T, self.v100.cSt)
            self.VGC = (10*self.SG-1.0752*log10(vSUS-38))/(10-log10(vSUS-38))

        if not self.hasSG:
            self.d20 = self.SG-4.5e-3*(2.34-1.9*self.SG)
        self.Ri = self.n-self.d20/2.
        self.m = self.M*(self.n-1.475)
        self.VI = self.Viscosity_Index()

        try:
            self.PourP = PourPoint(self.SG, self.Tb, self.v100)
        except NotImplementedError:
            self.PourP = None

        try:
            self.CloudP = CloudPoint(self.SG, self.Tb)
        except NotImplementedError:
            self.CloudP = None

        try:
            self.FreezingP = FreezingPoint(self.SG, self.Tb)
        except NotImplementedError:
            self.FreezingP = None

        try:
            self.AnilineP = AnilinePoint(self.SG, self.Tb)
            self.DieselI = self.API*self.AnilineP.F/100.
        except NotImplementedError:
            self.AnilineP = None
            self.DieselI = None

        try:
            self.SmokeP = SmokePoint(self.SG, self.Tb)
        except NotImplementedError:
            self.SmokeP = None

        try:
            self.CetaneI = CetaneIndex(self.API, self.Tb)
        except NotImplementedError:
            self.CetaneI = None

        if self.kwargs["S"]:
            self.S = self.kwargs["S"]
        else:
            self.S = S_Riazi(self.M, self.SG, self.Ri, self.m)

        if self.kwargs["H"]:
            self.H = self.kwargs["H"]
        else:
            self.H = self._H()

        self.N = self.kwargs["N"]

        # Conradson carbon residue, ref [29]_ Eq 3.141
        # HC: HC atomic ratio Eq. 2.122
        HC = 11.9147/self.CH
        CCR = 148 - 86.96/HC
        if CCR < 0:
            CCR = 0
        elif CCR > 100:
            CCR = 100
        self.CCR = unidades.Dimensionless(CCR)

        self.NetHeating, self.GrossHeating = CombustionHeat(self.API, S=self.S)

        newComponente.calculo(self)

    def tr(self, T):
        return T/self.Tc

    def pr(self, P):
        return P/self.Pc.atm

    def _M(self):
        """Molecular weight calculation"""
        # Method defined to calculate properties
        if self.definicion != 1:
            # The only defined method is the prop_Rizi_Daubert
            methodIndex = 0
        else:
            methodIndex = self.Preferences.getint("petro", "M")

        methods = [prop_Riazi_Daubert, prop_Riazi_Daubert_1980,
                   prop_Lee_Kesler, prop_Sim_Daubert, M_Goossens, prop_Twu]
        method = methods[methodIndex]
        if method.__name__ not in self.calculatedMethod:
            self.calculateMethod(method)
        M = self.calculatedMethod[method.__name__]["M"]
        return unidades.Dimensionless(M)

    def _critic(self):
        """Critic pressure and temperature calculation"""
        methodIndex = self.Preferences.getint("petro", "critical")

        # Discard method with incomplete input
        if methodIndex in [1, 2, 3, 4, 10, 11] and self.definicion != 1:
            methodIndex = 0
        if methodIndex in [7, 9] and self.definicion != 2:
            methodIndex = 0
        if methodIndex in [6, 12] and not self.kwargs["M"]:
            methodIndex = 0
        if methodIndex == 8 and not self.kwargs["Tb"]:
            methodIndex = 0
        if methodIndex in [13, 14] and self.definicion not in [1, 2]:
            methodIndex = 0

        methods = [prop_Riazi_Daubert, prop_Riazi_Daubert_1980, prop_Cavett,
                   prop_Lee_Kesler, prop_Sim_Daubert,
                   prop_Watansiri_Owens_Starling, prop_Rowe, prop_Standing,
                   prop_Willman_Teja, prop_Magoulas_Tassios, prop_Tsonopoulos,
                   prop_Twu, prop_Sancet, prop_Riazi, prop_Riazi_Alsahhaf]
        method = methods[methodIndex]
        if method.__name__ not in self.calculatedMethod:
            self.calculateMethod(method)
        Tc = self.calculatedMethod[method.__name__]["Tc"]
        Pc = self.calculatedMethod[method.__name__]["Pc"]
        return unidades.Temperature(Tc), unidades.Pressure(Pc)

    def _Vc(self):
        """Critic volume calculation"""
        methodIndex = self.Preferences.getint("petro", "vc")

        # Discard method with incomplete input
        if methodIndex in [1, 2, 4] and self.definicion != 1:
            methodIndex = 0
        if methodIndex == 5 and self.definicion != 2:
            methodIndex = 0
        if methodIndex in [3, 6, 7] and self.definicion not in [1, 2]:
            methodIndex = 0

        methods = [prop_Riazi_Daubert, prop_Riazi_Daubert_1980,
                   prop_Watansiri_Owens_Starling, prop_Twu, vc_Hall_Yarborough,
                   prop_Riazi, prop_Riazi_Alsahhaf]
        method = methods[methodIndex]
        if method.__name__ not in self.calculatedMethod:
            self.calculateMethod(method)
        Vc = self.calculatedMethod[method.__name__]["Vc"]
        return unidades.SpecificVolume(Vc)

    def _f_acent(self):
        """Acentric factor calculation"""
        methodIndex = self.Preferences.getint("petro", "f_acent")
        # This procedure has no limitation, the chossen method use calculated
        # properties is isn't in input

        methods = [prop_Edmister, prop_Lee_Kesler, w_Korsten,
                   prop_Watansiri_Owens_Starling, prop_Magoulas_Tassios]
        method = methods[methodIndex]
        if method.__name__ not in self.calculatedMethod:
            self.calculateMethod(method)
        w = self.calculatedMethod[method.__name__]["w"]
        return unidades.Dimensionless(w)

    def _Tb(self):
        """Boiling temperature calculation"""
        methodIndex = self.Preferences.getint("petro", "Tb")
        if self.definicion not in [1, 2] and methodIndex > 1:
            # Set the method to Riaxi Daubert if the choosen method is
            # unavailable for the input parameters
            methodIndex = 0

        methods = [prop_Riazi_Daubert, prop_Riazi, prop_Rowe, prop_Sancet,
                   prop_Silva_Rodriguez, Tb_Soreide]
        method = methods[methodIndex]
        if method.__name__ not in self.calculatedMethod:
            self.calculateMethod(method)
        Tb = self.calculatedMethod[method.__name__]["Tb"]
        return unidades.Temperature(Tb)

    def _SG(self):
        """Specific gravity procedure calculation"""
        methodIndex = self.Preferences.getint("petro", "SG")
        if methodIndex == 1 and not self.kwargs["M"]:
            # If molecular weight isn't in input parameter the Silva-Rodriguez
            # method can't be used
            methodIndex = 0

        methods = [prop_Riazi_Daubert, prop_Silva_Rodriguez]
        method = methods[methodIndex]
        if method.__name__ not in self.calculatedMethod:
            self.calculateMethod(method)
        SG = self.calculatedMethod[method.__name__]["SG"]
        return unidades.Dimensionless(SG)

    def _n(self):
        """Refractive index procedure calculation"""
        methodIndex = self.Preferences.getint("petro", "n")
        if methodIndex == 2 and not self.kwargs["Tb"]:
            # If boiling temperature isn't in input parameter the Willman-Teja
            # method can't be used
            methodIndex = 0
        if self.definicion not in [1, 2] and methodIndex == 1:
            methodIndex = 0

        methods = [prop_Riazi_Daubert, prop_Willman_Teja, prop_Riazi]
        method = methods[methodIndex]
        if method.__name__ not in self.calculatedMethod:
            self.calculateMethod(method)

        if method == prop_Willman_Teja:
            n = self.calculatedMethod[method.__name__]["n"]
            I = (n**2-1)/(n**2+2)
        else:
            I = self.calculatedMethod[method.__name__]["I"]
            n = ((1+2*I)/(1-I))**0.5
        return unidades.Dimensionless(n), unidades.Dimensionless(I)

    def _CH(self):
        """Carbon-Hydrogen procedure calculation"""
        method = prop_Riazi_Daubert
        if method.__name__ not in self.calculatedMethod:
            self.calculateMethod(method)
        CH = self.calculatedMethod[method.__name__]["CH"]
        return unidades.Dimensionless(CH)

    def calculateMethod(self, method):
        """Calculate method of properties and save to avoid other iteration"""
        Tb_SG = [prop_Riazi_Daubert_1980, prop_Lee_Kesler, prop_Sim_Daubert,
                 prop_Twu]
        M_SG = [prop_Standing, prop_Magoulas_Tassios, Tb_Soreide,
                vc_Hall_Yarborough]
        M = [prop_Rowe, prop_Sancet, prop_Silva_Rodriguez]

        if method == prop_Riazi_Daubert:
            # Special case for Riazi_daubert advanced method
            if self.definicion in [1, 8]:
                prop = prop_Riazi_Daubert("Tb", self.Tb, "SG", self.SG)
            elif self.definicion == 2:
                prop = prop_Riazi_Daubert("M", self.M, "SG", self.SG)
            elif self.definicion == 3:
                prop = prop_Riazi_Daubert("Tb", self.Tb, "I", self.I)
            elif self.definicion == 4:
                prop = prop_Riazi_Daubert("M", self.M, "I", self.I)
            elif self.definicion == 5:
                prop = prop_Riazi_Daubert("Tb", self.Tb, "CH", self.CH)
            elif self.definicion == 6:
                prop = prop_Riazi_Daubert("M", self.M, "CH", self.CH)
            elif self.definicion == 7:
                prop = prop_Riazi_Daubert("M", self.v100, "I", self.I)

        elif method == prop_Riazi:
            if self.definicion in [1, 8]:
                prop = prop_Riazi(self.SG, "Tb", self.Tb)
            elif self.definicion == 2:
                prop = prop_Riazi(self.SG, "M", self.M)

        elif method == prop_Cavett:
            prop = prop_Cavett(self.Tb, self.API)
        elif method == M_Goossens:
            prop = M_Goossens(self.Tb, self.d20)
        elif method == prop_Willman_Teja:
            prop = prop_Willman_Teja(self.Tb)
        elif method == prop_Watansiri_Owens_Starling:
            prop = prop_Watansiri_Owens_Starling(self.Tb, self.SG, self.M)
        elif method == prop_Riazi_Alsahhaf:
            prop = prop_Riazi_Alsahhaf(self.Tb, self.M, self.d20)
        elif method == prop_Edmister:
            prop = prop_Edmister(Tb=self.Tb, Tc=self.Tc, Pc=self.Pc)
        elif method in Tb_SG:
            prop = method(self.Tb, self.SG)
        elif method in M_SG:
            prop = method(self.M, self.SG)
        elif method in M:
            prop = method(self.M)
        self.calculatedMethod[method.__name__] = prop

    def _Zc(self):
        methods = [Zc_Lee_Kesler, Zc_Salerno, Zc_Nath, Zc_Reid, Zc_Hougen]
        Zc = methods[self.Preferences.getint("petro", "Zc")]
        return Zc(self.f_acent)

    def _H(self):
        if self.Preferences.getint("petro", "H") == 0:
            H = H_Riazi(self.S, self.CH)
        elif self.Preferences.getint("petro", "H") == 1:
            H = H_Goossens(self.M, self.n, self.d20)
        elif self.Preferences.getint("petro", "H") == 2:
            if self.hasCurve:
                Tb = (self.T10+self.T50+self.T90)/3.
            else:
                Tb = self.Tb
            H = H_ASTM(Tb, self.SG, self.xa)
        else:
            H = H_Jenkins_Walsh(self.SG, self.AnilineP)
        return H

    def _D86_TBP(self, D86, reverse=False):
        """Use the preferences to calcuate the TBP distillation curve from D86
        with the desired method"""
        index = self.Preferences.getint("petro", "curve")
        method = [D86_TBP_Riazi, D86_TBP_Daubert][index]
        TBP = method(D86, reverse)
        return TBP

    def _SD_D86(self, SD):
        """Use the preferences to calcuate the D86 distillation curve from SD
        with the desired method"""
        index = self.Preferences.getint("petro", "curve")
        method = [SD_D86_Riazi, SD_D86_Daubert][index]
        D86 = method(SD)
        return D86

    def _PNA(self):
        """Calculate the Paraffin-Naphthenes-Aromatics group composition"""
        if self.Preferences.getint("petro", "PNA") == 0:
            xp, xn, xa = PNA_Riazi(
                self.M, self.SG, self.n, d20=None, VGC=self.VGC, CH=None)
        elif self.Preferences.getint("petro", "PNA") == 1:
            xp, xn, xa = PNA_Bergman(self.Tb, self.SG, self.Kw)
        elif self.Preferences.getint("petro", "PNA") == 2:
            xp, xn, xa = PNA_Peng_Robinson(self.Nc, self.M, self.WABP)
        else:
            xp, xn, xa = PNA_van_Nes(self.M, self.n, self.d20, self.S)
        return xp, xn, xa

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

    # petroleo=Petroleo(M=250, API=43.3)
    # t=unidades.Temperature(100, "F")
#    print petroleo.Mu_Gas(t).cP
#    print petroleo.Tension_API(t).dyncm
#    print petroleo.Tension_Baker_Swerdloff(t).dyncm
#    print petroleo.API

    # petroleo=Petroleo(API=22.5, M=339.7)
    # print(petroleo.API, petroleo.M)

    v = [10, 30, 50, 70, 90]
    D1160 = [i+273.15 for i in [150, 205, 250, 290, 350]]
    petroleo = Petroleo(P_dist=10, curveType="D1160", X_curve=v, T_curve=D1160)
