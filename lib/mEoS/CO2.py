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


from unittest import TestCase

from scipy import exp, arccosh

from lib.meos import MEoS
from lib import unidades


class CO2(MEoS):
    """Multiparameter equation of state for carbon dioxide"""
    name = "carbon dioxide"
    CASNumber = "124-38-9"
    formula = "CO2"
    synonym = "R-744"
    _refPropName = "CO2"
    _coolPropName = "CarbonDioxide"
    rhoc = unidades.Density(467.6)
    Tc = unidades.Temperature(304.1282)
    Pc = unidades.Pressure(7.3773, "MPa")
    M = 44.0098  # g/mol
    Tt = unidades.Temperature(216.592)
    Tb = unidades.Temperature(194.686)
    f_acent = 0.22394
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 49

    Fi1 = {"ao_log": [1, 2.5],
           "pow": [0, 1],
           "ao_pow": [8.37304456, -3.70454304],
           "ao_exp": [1.99427042, .62105248, .41195293, 1.04028922, .08327678],
           "titao": [3.15163, 6.11190, 6.77708, 11.32384, 27.08792],
           "ao_hyp": [], "hyp": []}

    Fi2 = {"ao_log": [1, 2.50002],
           "pow": [0, 1],
           "ao_pow": [11.925152758, -16.118762264],
           "ao_exp": [], "titao": [],
           "ao_hyp": [2.04452, -1.06044, 2.03366, 0.01393],
           "hyp": [3.022758166, -2.844425476, 1.589964364, 1.12159609]}

    CP3 = {"ao": 3.5,
           "an": [], "pow": [],
           "ao_exp": [2, 1, 1], "exp": [960.11, 1932, 3380.2],
           "ao_hyp": [], "hyp": []}

    span = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for carbon dioxide of Span "
                    "and Wagner (1996)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "A New Equation of State for Carbon Dioxide "
                             "Covering the Fluid Region from the Triple‐Point "
                             "Temperature to 1100K at Pressures up to 800MPa",
                    "ref": "J. Phys. Chem. Ref. Data, 25(6) (1996) 1509-1596",
                    "doi": "10.1063/1.555991"},

        "R": 8.31451,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 2000., "Pmax": 800000.0, "rhomax": 37.24,
        "Pmin": 517.95, "rhomin": 26.777,

        "nr1": [0.38856823203161, 0.29385475942740e1, -0.55867188534934e1,
                -0.76753199592477, 0.31729005580416, 0.54803315897767,
                0.12279411220335],
        "d1": [1, 1, 1, 1, 2, 2, 3],
        "t1": [0.0, 0.75, 1.0, 2.0, 0.75, 2.0, 0.75],

        "nr2": [0.21658961543220e1, 0.15841735109724e1, -0.23132705405503,
                0.58116916431436e-1, -0.55369137205382, 0.48946615909422,
                -0.24275739843501e-1, 0.62494790501678e-1, -0.12175860225246,
                -0.37055685270086, -0.16775879700426e-1, -0.11960736637987,
                -0.045619362508778, 0.35612789270346e-1, -0.74427727132052e-2,
                -0.17395704902432e-2, -0.021810121289527, 0.24332166559236e-1,
                -0.37440133423463e-1, 0.14338715756878, -0.13491969083286,
                -0.23151225053480e-1, 0.12363125492901e-1, 0.21058321972940e-2,
                -0.33958519026368e-3, 0.0055993651771592, -.30335118055646e-3],
        "d2": [1, 2, 4, 5, 5, 5, 6, 6, 6, 1, 1, 4, 4, 4, 7, 8, 2, 3, 3, 5, 5,
               6, 7, 8, 10, 4, 8],
        "t2": [1.5, 1.5, 2.5, 0, 1.5, 2, 0, 1, 2, 3, 6, 3, 6, 8, 6, 0, 7, 12,
               16, 22, 24, 16, 24, 8, 2, 28, 14],
        "c2": [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4,
               4, 4, 4, 4, 5, 6],
        "gamma2": [1]*27,

        "nr3": [-0.21365488688320e3, 0.26641569149272e5, -0.24027212204557e5,
                -0.28341603423999e3, 0.21247284400179e3],
        "d3": [2, 2, 2, 3, 3],
        "t3": [1., 0., 1., 3., 3.],
        "alfa3": [25, 25, 25, 15, 20],
        "beta3": [325, 300, 300, 275, 275],
        "gamma3": [1.16, 1.19, 1.19, 1.25, 1.22],
        "epsilon3": [1.]*5,

        "nr4": [-0.66642276540751, 0.72608632349897, 0.55068668612842e-1],
        "a4": [3.5, 3.5, 3.],
        "b4": [0.875, 0.925, 0.875],
        "beta4": [0.3]*3,
        "A": [0.7]*3,
        "B": [0.3, 0.3, 1.],
        "C": [10., 10., 12.5],
        "D": [275]*3}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for carbon dioxide of Ely (1987)",
        "__doi__": {"autor": "Ely, J.F., Magee, J.W., and Haynes, W.M.",
                    "title": "Thermophysical properties for special high CO2 "
                             "content mixtures",
                    "ref": "Research Report RR-110, Gas Processors "
                           "Association, Tulsa, OK, 1987.",
                    "doi": ""},

        "R": 8.31434,
        "cp": CP3,
        "ref": "OTO",

        "Tmin": 216.58, "Pmin": 518.2, "rhomin": 26.778,
        "Tmax": 440.1, "Pmax": 40000.0, "rhomax": 27.778,

        "b": [None, -0.981851065838e-2, 0.995062267309, -0.228380160313e2,
              0.281827634529e4, -0.347001262699e6, 0.394706709102e-3,
              -0.325550000110, 0.484320083063e1, -0.352181542995e6,
              -0.324053603343e-4, 0.468596684665e-1, -0.754547012075e1,
              -0.381894354016e-4, -0.442192933859e-1, 0.516925168095e2,
              0.212450985237e-2, -0.261009474785e-4, -0.888533388977e-1,
              0.155226179403e-2, 0.415091004940e6, -0.110173967489e8,
              0.291990583344e4, 0.143254606508e8, 0.108574207533e2,
              -0.247799657039e3, 0.199293590763e-1, 0.102749908059e3,
              0.377618865158e-4, -0.332276512346e-2, 0.179196707121e-7,
              0.945076627807e-5, -0.123400943061e-2]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for carbon dioxide of Kunz "
                    "and Wagner (2004)",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for "
                             "Natural Gases and Other Mixtures: An Expansion "
                             "of GERG-2004",
                    "ref": "J. Chem.Eng. Data 57(11) (2012) 3032-3091",
                    "doi":  "10.1021/je300655b"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 1100., "Pmax": 800000.0, "rhomax": 37.24,
        "Pmin": 0.61166, "rhomin": 55.497,

        "nr1": [0.52646564804653, -0.14995725042592e1, 0.27329786733782,
                0.12949500022786],
        "d1": [1, 1, 2, 3],
        "t1": [0, 1.25, 1.625, 0.375],

        "nr2": [0.15404088341841, -0.58186950946814, -0.18022494838296,
                -0.095389904072812, -0.80486819317679e-2, -0.03554775127309,
                -0.28079014882405, -0.82435890081677e-1, 0.10832427979006e-1,
                -0.67073993161097e-2, -0.46827907600524e-2, -0.028359911832177,
                0.19500174744098e-1, -0.21609137507166, 0.43772794926972,
                -0.22130790113593, 0.15190189957331e-1, -0.15380948953300e-1],
        "d2": [3, 3, 4, 5, 6, 6, 1, 4, 1, 1, 3, 3, 4, 5, 5, 5, 5, 5],
        "t2": [0.375, 1.375, 1.125, 1.375, 0.125, 1.625, 3.75, 3.5, 7.5, 8, 6,
               16, 11, 24, 26, 28, 24, 26],
        "c2": [1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 5, 5, 5, 6, 6],
        "gamma2": [1]*18}

    ely = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for carbon dioxide of Ely "
                    "(1987)",
        "__doi__": {"autor": "Ely, J.F., Magee, J.W., and Haynes, W.M.",
                    "title": "Thermophysical properties for special high CO2 "
                             "content mixtures",
                    "ref": "Research Report RR-110, Gas Processors "
                           "Association, Tulsa, OK, 1987.",
                    "doi": ""},

        "R": 8.31434,
        "cp": CP3,
        "ref": "OTO",

        "Tmin": 216.58, "Tmax": 1000., "Pmax": 100000.0, "rhomax": 26.776,
        "Pmin": 518.03, "rhomin": 26.776,

        "nr1": [0.485497428986, -0.191900462349e1, 0.451739876847,
                0.838475229022e-2, 0.310719428397, -0.183619563850,
                0.448878785519e-1, -0.362211893044e-1, -0.169827491865e-1,
                0.803504394396e-3, 0.320223641512e-3, -0.658956249553e-5,
                -0.461991678692e-4],
        "d1": [1, 1, 1, 2, 2, 2, 3, 3, 3, 6, 7, 7, 8],
        "t1": [0, 1.5, 2.5, -0.5, 1.5, 2, 0, 1, 2.5, 0, 2, 5, 2],

        "nr2": [-0.385989029443, 0.131878614095, 0.109639470331,
                -0.310044422115e-1, -0.989797992915e-1, -0.222934996927e-1,
                -0.225488505376e-1, -0.595661202393e-2, -0.219959964099e-1,
                0.140330955537e-1, -0.315424157971e-2, 0.443394060420e-3,
                -0.487628903103e-2, -0.311643343682e-1, 0.226083669848e-1,
                0.186651858191e-1, -0.399277963883, 0.464945130861,
                -0.817090055061e-1],
        "d2": [1, 1, 2, 2, 3, 3, 5, 6, 7, 8, 10, 2, 3, 3, 4, 4, 5, 5, 5, ],
        "t2": [5, 6, 3.5, 5.5, 3, 7, 6, 8.5, 4, 6.5, 5.5, 22, 11, 18, 11, 23,
               17, 18, 23],
        "c2": [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4],
        "gamma2": [1]*19}

    shortSpan = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for carbon dioxide of "
                    "Span and Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of State for Technical Applications. "
                             "III. Results for Polar Fluids",
                    "ref": "Int. J. Thermophys., 24(1) (2003) 111-162",
                    "doi": "10.1023/A:1022362231796"},

        "R": 8.31451,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 600., "Pmax": 100000.0, "rhomax": 37.24,
        "Pmin": 517.86, "rhomin": 26.795,

        "nr1": [0.89875108, -0.21281985e1, -0.6819032e-1, 0.76355306e-1,
                0.22053253e-3],
        "d1": [1, 1, 1, 3, 7],
        "t1": [.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.41541823, 0.71335657, 0.30354234e-3, -0.36643143,
                -0.14407781e-2, -0.89166707e-1, -0.23699887e-1],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    sun = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for carbon dioxide of Sun "
                    "and Ely (2004)",
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering "
                             "application: Algorithm and  application to "
                             "non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223 (2004) 107-118",
                    "doi": "10.1016/j.fluid.2004.06.028"},

        "R": 8.31451,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40.,
        "Pmin": 0.1, "rhomin": 40.,

        "nr1": [-4.71122371e-1, 9.13375599e-1, -1.96793707, 6.89687161e-2,
                2.15658922e-4, 9.51876380e-2],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [-4.91366518e-3, 7.32487713e-1, 8.70918629e-1, -5.35917679e-3,
                -0.403818537, -2.40820897e-2, -1.04239403e-1, -2.16335828e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    # eq = span, MBWR, GERG, ely, shortSpan, sun
    eq = span, GERG, ely, shortSpan, sun

    _surface = {"sigma": [0.07863], "exp": [1.254]}
    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [],  "expt0": [], "expd0": [],
                   "a1": [7.3455, 0.00335], "expt1": [0, 1], "expd1": [1, 1],
                   "a2": [83.93, 145.1, -578.8, -1012.],
                   "expt2": [0, 1, 0, 1], "expd2": [2, 2, 2.55, 2.55]}
    _melting = {"eq": 1, "Tref": Tt, "Pref": 517.95,
                "Tmin": Tt, "Tmax": 1100.0,
                "a1": [1], "exp1": [0],
                "a2": [1955.539, 2055.4593], "exp2": [1, 2],
                "a3": [], "exp3": []}
    _sublimation = {"eq": 3, "Tref": Tt, "Pref": 517.950,
                    "Tmin": Tt, "Tmax": Tt,
                    "a1": [], "exp1": [],
                    "a2": [-14.740846, 2.4327015, -5.3061778],
                    "exp2": [1, 1.9, 2.9],
                    "a3": [], "exp3": []}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-7.0602087, 1.9391218, -1.6463597, -3.2995634],
        "exp": [1, 1.5, 2., 4.]}
    _liquid_Density = {
        "eq": 4,
        "ao": [1.92451080, -0.62385555, -0.32731127, 0.39245142],
        "exp": [1.02, 1.5, 5.0, 5.5]}
    _vapor_Density = {
        "eq": 4,
        "ao": [-1.7074879, -0.8227467, -4.6008549, -10.111178, -29.742252],
        "exp": [1.02, 1.5, 3.0, 7.0, 14.0]}

    visco0 = {"__name__": "Fenghour (1998)",
              "__doi__": {
                  "autor": "Fenghour, A., Wakeham, W.A., Vesovic, V.",
                  "title": "The Viscosity of Carbon Dioxide",
                  "ref": "J. Phys. Chem. Ref. Data 27(1) (1998) 31-44",
                  "doi": "10.1063/1.556013"},

              "eq": 1, "omega": 1,

              "ek": 251.196, "sigma": 1.,
              "n_chapman": 1.00697/M**0.5,
              "collision": [0.235156, -0.491266, 5.211155e-2, 5.347906e-2,
                            -1.537102e-2],

              "Tref_res": 251.196, "rhoref_res": 1,
              "nr": [0.4071119e-2, 0.7198037e-4, 0.2411697e-16, 0.2971072e-22,
                     -0.1627888e-22],
              "tr": [0, 0, 3, 0, 1],
              "dr": [1, 2, 6, 8, 8],
              "gr": [0, 0, 0, 0, 0],
              "cr": [0, 0, 0, 0, 0]}

    visco1 = {"__name__": u"Quiñones-Cisneros (2006)",
              "__doi__": {
                  "autor": "Quiñones-Cisneros, S.E., Deiters, U.K.",
                  "title": "Generalization of the Friction Theory for "
                           "Viscosity Modeling",
                  "ref": "J. Phys. Chem. B, 110(25) (2006) 12820-12834",
                  "doi": "10.1021/jp0618577"},

              "eq": 4, "omega": 0,

              "Toref": 304.1282,
              "no": [69.18424, -215.8618, 210.94362, -49.0494],
              "to": [0, 0.25, 0.5, 0.75],

              "a": [1.19805e-4,  -1.25861e-4, 5.48871e-5],
              "b": [3.15921e-5, -2.60469e-5, 7.09199e-6],
              "c": [1.80689e-5, -7.41742e-6, 0.0],
              "A": [-2.31066e-9, 0.0, 5.42486e-10],
              "B": [1.04558e-8, -2.20758e-9, 0.0],
              "C": [1.03255e-6, -8.56207e-7, 3.84384e-7]}

    visco2 = {"eq": 0,
              "method": "_visco2",
              "__name__": "Vesovic (1990)",
              "__doi__": {
                  "autor": "Vesovic, V., Wakeham, W.A., Olchowy, G.A., "
                           "Sengers, J.V., Watson, J.T.R., Millat, J.",
                  "title": "The Transport Properties of Carbon Dioxide",
                  "ref": "J. Phys. Chem. Ref. Data 19(3) (1990) 763-808",
                  "doi": "10.1063/1.555875"},

              "omega": 1,
              "ek": 251.196, "sigma": 0.3751,
              "n_chapman": 1.00697/M**0.5*0.3751**2,
              "collision": [0.235156, -0.491266, 5.211155e-2, 5.347906e-2,
                            -1.537102e-2]}

    def _visco2(self, rho, T, fase=None):

        # Zero-Density viscosity
        muo = self._Visco0()

        # Gas-phase viscosity
        # Table 8
        ti = [1, 2, 7]
        ei = [3.6350734e-3, 7.209997e-5, 3.00306e-20]
        mug = 0
        for t, e in zip(ti, ei):
            mug += e*rho**t                                             # Eq 67

        # Liquid-phase viscosity
        B = 18.56+0.014*T                                               # Eq 71
        Vo = 7.41e-4-3.3e-7*T                                           # Eq 72
        mul = 1/B/(1/rho-Vo)                                            # Eq 70

        # Critical enhancement
        # TODO: Not implemented
        muc = 0

        # Eq 74, parameters
        rhos = 467.689
        Ts = 302
        Z = 1

        # Eq 76
        mu = muo + mug/(1+exp(-Z*(T-Ts))) + \
            mug/(1+exp(Z*(T-Ts)))/(1+exp(Z*(rho-rhos))) + \
            (mul-muo)/(1+exp(Z*(T-Ts)))/(1+exp(-Z*(rho-rhos))) + muc

        return unidades.Viscosity(mu, "muPas")

    _viscosity = visco0, visco1, visco2

    thermo0 = {"__name__": "Scalabrin (2006)",
               "__doi__": {
                  "autor": "Scalabrin, G., Marchi, P., Finezzo, F.",
                  "title": "A Reference Multiparameter Thermal Conductivity "
                           "Equation for Carbon Dioxide with an Optimized "
                           "Functional Form",
                  "ref": "J. Phys. Chem. Ref. Data 35(4) (2006) 1549-1575",
                  "doi": "10.1063/1.2213631"},

               "eq": 1,

               "Tref_res": 304.1282, "rhoref_res": 467.6,
               "kref_res": 4.81384e-3,
               "nr": [7.69857587, 0.159885811, 1.56918621, -6.73400790,
                      16.3890156, 3.69415242, 22.3205514, 66.1420950,
                      -0.171779133, 0.00433043347],
               "tr": [0, 0, -1.5, 0, -1, -1.5, -1.5, -1.5, -3.5, -5.5],
               "dr": [1, 5, 1, 1, 2, 0, 5, 9, 0, 0],
               "cr": [0, 0, 0, 2, 2, 2, 2, 2, 2, 2],
               "gr": [0, 0, 0, 5, 5, 5, 5, 5, 5, 5],

               "critical": "_tc"}

    def _tc(self, rho, T, fase):
        """Custom method for critical enhancement"""

        Tr = T/304.1282
        rhor = rho/467.6
        nc = 0.775547504

        # Table 4
        a = [0, 3, 6.70697, 0.94604, 0.3, 0.3, 0.39751, 0.33791, 0.77963,
             0.79857, 0.9, 0.02, 0.2]

        # Eq 6
        alfa = 1-a[10]*arccosh(1+a[11]*((1-Tr)**2)**a[12])

        # Eq 5
        num = rhor*exp(-rhor**a[1]/a[1]-(a[2]*(Tr-1))**2-(a[3]*(rhor-1))**2)
        den1 = pow(pow(1-1/Tr+a[4]*pow(pow(rhor-1, 2), 0.5/a[5]), 2), a[6])
        den2 = pow(pow(a[7]*(rhor-alfa), 2), a[8])
        lc = num / (den1+den2)**a[9]

        return lc*nc*4.81384e-3

    thermo1 = {"__name__": "Vesovic (1990)",
               "__doi__": {
                  "autor": "Vesovic, V., Wakeham, W.A., Olchowy, G.A., "
                           "Sengers, J.V., Watson, J.T.R., Millat, J.",
                  "title": "The Transport Properties of Carbon Dioxide",
                  "ref": "J. Phys. Chem. Ref. Data 19(3) (1990) 763-808",
                  "doi": "10.1063/1.555875"},

               "eq": 1,

               "Toref": 251.196, "koref": 1e-3,
               "no_num": [7.5378307, 4.8109652e-2],
               "to_num": [0.5, -99],
               "no_den": [0.4226159, 0.6280115, -0.5387661, 0.6735941,
                          -0.4362677, 0.2255388],
               "to_den": [0, 1, 2, 3, 6, 7],

               "Tref_res": 1., "rhoref_res": 2.272221e-2, "kref_res": 1e-3,
               "nr": [2.447164e-2, 8.705605e-5, -6.547950e-8, 6.594919e-11],
               "tr": [0, 0, 0, 0],
               "dr": [1, 2, 3, 4],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.2415, "R0": 1.01,
               "Xio": 1.5e-10, "gam0": 0.052, "qd": 0.4e-9, "Tcref": 450.}

    _thermal = thermo0, thermo1


class Test(TestCase):

    def test_span(self):
        # Selected point from Table 34, Pag 1560, saturation state
        st = CO2(T=216.592, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 0.51796)
        self.assertEqual(round(st.Liquido.rho, 2), 1178.46)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -426.74)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), -2.2177)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 5), 0.97466)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 1.9532)
        self.assertEqual(round(st.Liquido.w, 2), 975.85)
        self.assertEqual(round(st.Gas.rho, 3), 13.761)
        self.assertEqual(round(st.Gas.h.kJkg, 3), -76.364)
        self.assertEqual(round(st.Gas.s.kJkgK, 5), -0.59999)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 0.6292)
        self.assertEqual(round(st.Gas.cp.kJkgK, 5), 0.90872)
        self.assertEqual(round(st.Gas.w, 2), 222.78)

        st = CO2(T=240, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 1.2825)
        self.assertEqual(round(st.Liquido.rho, 2), 1088.87)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -379.94)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), -2.0155)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 5), 0.94535)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 2.0510)
        self.assertEqual(round(st.Liquido.w, 2), 806.38)
        self.assertEqual(round(st.Gas.rho, 3), 33.295)
        self.assertEqual(round(st.Gas.h.kJkg, 3), -70.293)
        self.assertEqual(round(st.Gas.s.kJkgK, 5), -0.72532)
        self.assertEqual(round(st.Gas.cv.kJkgK, 5), 0.70534)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 1.1033)
        self.assertEqual(round(st.Gas.w, 2), 222.96)

        st = CO2(T=260, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 2.4188)
        self.assertEqual(round(st.Liquido.rho, 2), 998.89)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -337.34)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), -1.8495)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 5), 0.93227)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 2.2554)
        self.assertEqual(round(st.Liquido.w, 2), 652.58)
        self.assertEqual(round(st.Gas.rho, 3), 64.417)
        self.assertEqual(round(st.Gas.h.kJkg, 3), -70.862)
        self.assertEqual(round(st.Gas.s.kJkgK, 5), -0.82456)
        self.assertEqual(round(st.Gas.cv.kJkgK, 5), 0.79426)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 1.4295)
        self.assertEqual(round(st.Gas.w, 2), 218.19)

        st = CO2(T=280, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 4.1607)
        self.assertEqual(round(st.Liquido.rho, 2), 883.58)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -289.48)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), -1.6792)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 5), 0.96046)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 2.8141)
        self.assertEqual(round(st.Liquido.w, 2), 471.54)
        self.assertEqual(round(st.Gas.rho, 2), 121.74)
        self.assertEqual(round(st.Gas.h.kJkg, 3), -80.840)
        self.assertEqual(round(st.Gas.s.kJkgK, 5), -0.93401)
        self.assertEqual(round(st.Gas.cv.kJkgK, 5), 0.92316)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 2.2769)
        self.assertEqual(round(st.Gas.w, 2), 207.72)

        st = CO2(T=300, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 6.7131)
        self.assertEqual(round(st.Liquido.rho, 2), 679.24)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -223.40)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), -1.4631)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 1.1199)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 8.6979)
        self.assertEqual(round(st.Liquido.w, 2), 245.67)
        self.assertEqual(round(st.Gas.rho, 2), 268.58)
        self.assertEqual(round(st.Gas.h.kJkg, 2), -119.70)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), -1.1175)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 1.2476)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 11.921)
        self.assertEqual(round(st.Gas.w, 2), 185.33)

        st = CO2(T=304, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 7.3555)
        self.assertEqual(round(st.Liquido.rho, 2), 530.30)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -188.42)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), -1.3509)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 2.0531)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 2), 386.88)
        self.assertEqual(round(st.Liquido.w, 2), 134.14)
        self.assertEqual(round(st.Gas.rho, 2), 406.42)
        self.assertEqual(round(st.Gas.h.kJkg, 2), -158.84)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), -1.2536)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 2.0679)
        self.assertEqual(round(st.Gas.cp.kJkgK, 2), 555.58)
        self.assertEqual(round(st.Gas.w, 2), 147.62)

        # Selected point from Table 35, Pag 1562, single phase region
        st = CO2(T=240, P=5e4)
        self.assertEqual(round(st.rho, 4), 1.1084)
        self.assertEqual(round(st.u.kJkg, 3), -93.125)
        self.assertEqual(round(st.h.kJkg, 3), -48.014)
        self.assertEqual(round(st.s.kJkgK, 5), -0.04483)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.59473)
        self.assertEqual(round(st.cp.kJkgK, 5), 0.78815)
        self.assertEqual(round(st.w, 2), 243.88)

        st = CO2(T=1100, P=1e5)
        self.assertEqual(round(st.rho, 5), 0.48109)
        self.assertEqual(round(st.u.kJkg, 2), 675.83)
        self.assertEqual(round(st.h.kJkg, 2), 883.69)
        self.assertEqual(round(st.s.kJkgK, 4), 1.3827)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.0702)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.2593)
        self.assertEqual(round(st.w, 2), 494.61)

        st = CO2(T=280, P=5e5)
        self.assertEqual(round(st.rho, 4), 9.7568)
        self.assertEqual(round(st.u.kJkg, 3), -71.840)
        self.assertEqual(round(st.h.kJkg, 3), -20.594)
        self.assertEqual(round(st.s.kJkgK, 5), -0.36761)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.65255)
        self.assertEqual(round(st.cp.kJkgK, 5), 0.87063)
        self.assertEqual(round(st.w, 2), 257.27)

        st = CO2(T=230, P=1e6)
        self.assertEqual(round(st.rho, 2), 1128.97)
        self.assertEqual(round(st.u.kJkg, 2), -401.07)
        self.assertEqual(round(st.h.kJkg, 2), -400.19)
        self.assertEqual(round(st.s.kJkgK, 4), -2.1006)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.95680)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.9959)
        self.assertEqual(round(st.w, 2), 879.82)

        st = CO2(T=290, P=5e6)
        self.assertEqual(round(st.rho, 2), 148.41)
        self.assertEqual(round(st.u.kJkg, 2), -115.58)
        self.assertEqual(round(st.h.kJkg, 3), -81.892)
        self.assertEqual(round(st.s.kJkgK, 5), -0.95959)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.94334)
        self.assertEqual(round(st.cp.kJkgK, 4), 2.5783)
        self.assertEqual(round(st.w, 2), 207.56)

        st = CO2(T=480, P=1e7)
        self.assertEqual(round(st.rho, 2), 119.65)
        self.assertEqual(round(st.u.kJkg, 3), 50.083)
        self.assertEqual(round(st.h.kJkg, 2), 133.66)
        self.assertEqual(round(st.s.kJkgK, 5), -0.48595)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.84672)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.1732)
        self.assertEqual(round(st.w, 2), 328.74)

        st = CO2(T=226.679, P=5e7)
        self.assertEqual(round(st.rho, 2), 1229.78)
        self.assertEqual(round(st.u.kJkg, 2), -430.79)
        self.assertEqual(round(st.h.kJkg, 2), -390.13)
        self.assertEqual(round(st.s.kJkgK, 4), -2.2377)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.99951)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.7735)
        self.assertEqual(round(st.w, 1), 1146.9)

        st = CO2(T=1100, P=1e8)
        self.assertEqual(round(st.rho, 2), 371.36)
        self.assertEqual(round(st.u.kJkg, 2), 617.86)
        self.assertEqual(round(st.h.kJkg, 2), 887.14)
        self.assertEqual(round(st.s.kJkgK, 5), 0.03152)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.0931)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.3544)
        self.assertEqual(round(st.w, 2), 675.48)

        st = CO2(T=327.673, P=8e8)
        self.assertEqual(round(st.rho, 2), 1495.70)
        self.assertEqual(round(st.u.kJkg, 2), -369.91)
        self.assertEqual(round(st.h.kJkg, 2), 164.96)
        self.assertEqual(round(st.s.kJkgK, 4), -2.1926)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.1961)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.5477)
        self.assertEqual(round(st.w, 1), 2052.8)

    def test_shortSpan(self):
        # Table III, Pag 117
        st = CO2(T=500, rho=500, eq="shortSpan")
        self.assertEqual(round(st.cp0.kJkgK, 4), 1.0141)
        self.assertEqual(round(st.P.MPa, 3), 45.164)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.4994)

        st2 = CO2(T=600, rho=100, eq="shortSpan")
        self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 191.33)
        self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 0.60315)

    def test_fenghour(self):
        # Table 13, Pag 44
        self.assertEqual(round(CO2(T=220, rho=2.440).mu.muPas, 2), 11.06)
        self.assertEqual(round(CO2(T=300, rho=1.773).mu.muPas, 2), 15.02)
        self.assertEqual(round(CO2(T=800, rho=0.662).mu.muPas, 2), 35.09)
        self.assertEqual(round(CO2(T=304, rho=254.320).mu.muPas, 2), 20.90)
        self.assertEqual(round(CO2(T=220, rho=1194.86).mu.muPas, 2), 269.37)
        self.assertEqual(round(CO2(T=300, rho=1029.27).mu.muPas, 2), 132.55)
        self.assertEqual(round(CO2(T=800, rho=407.828).mu.muPas, 2), 48.74)

    def test_Scalabrin(self):
        # Selected values from Table 10, Pag 1568, saturation states
        st = CO2(T=218, x=0.5)
        self.assertEqual(round(st.Liquido.k.mWmK, 2), 181.09)
        self.assertEqual(round(st.Gas.k.mWmK, 3), 10.837)

        st = CO2(T=250, x=0.5)
        self.assertEqual(round(st.Liquido.k.mWmK, 2), 138.61)
        self.assertEqual(round(st.Gas.k.mWmK, 3), 14.227)

        st = CO2(T=274, x=0.5)
        self.assertEqual(round(st.Liquido.k.mWmK, 2), 109.11)
        self.assertEqual(round(st.Gas.k.mWmK, 3), 19.417)

        st = CO2(T=284, x=0.5)
        self.assertEqual(round(st.Liquido.k.mWmK, 3), 96.920)
        self.assertEqual(round(st.Gas.k.mWmK, 3), 24.059)

        st = CO2(T=304, x=0.5)
        self.assertEqual(round(st.Liquido.k.mWmK, 2), 140.30)
        self.assertEqual(round(st.Gas.k.mWmK, 2), 217.95)

        # Selected values from Table 11, pag 1569, single phase states
        self.assertEqual(round(CO2(T=225, P=1e4).k.mWmK, 3), 11.037)
        self.assertEqual(round(CO2(T=300, P=1e6).k.mWmK, 3), 17.248)
        self.assertEqual(round(CO2(T=300, P=2e8).k.mWmK, 2), 219.64)
        self.assertEqual(round(CO2(T=400, P=1e7).k.mWmK, 3), 31.504)
        self.assertEqual(round(CO2(T=500, P=1e5).k.mWmK, 3), 33.143)
        self.assertEqual(round(CO2(T=700, P=8e6).k.mWmK, 3), 52.078)
        self.assertEqual(round(CO2(T=1000, P=2e8).k.mWmK, 2), 116.65)

    def test_vesovic(self):
        # FIXME: Vesovic thermal conductivity correlation dont work
        # Appendix IV, Pag 808
        st = CO2(T=220, rho=2.440, visco=2)
        # self.assertEqual(round(st.k.mWmK, 2), 10.90)
        self.assertEqual(round(st.mu.muPas, 2), 11.06)
        st = CO2(T=300, rho=1.773, visco=2)
        # self.assertEqual(round(st.k.mWmK, 2), 16.77)
        self.assertEqual(round(st.mu.muPas, 2), 15.02)
        st = CO2(T=800, rho=0.662, visco=2)
        # self.assertEqual(round(st.k.mWmK, 2), 56.65)
        self.assertEqual(round(st.mu.muPas, 2), 35.09)
        st = CO2(T=304, rho=254.3205, visco=2)
        # self.assertEqual(round(st.k.mWmK, 2), 42.52)
        self.assertEqual(round(st.mu.muPas, 2), 20.80)
        st = CO2(T=220, rho=1194.86, visco=2)
        # self.assertEqual(round(st.k.mWmK, 2), 187.50)
        self.assertEqual(round(st.mu.muPas, 2), 274.22)
        st = CO2(T=300, rho=1029.27, visco=2)
        # self.assertEqual(round(st.k.mWmK, 2), 137.61)
        self.assertEqual(round(st.mu.muPas, 2), 133.15)
        st = CO2(T=800, rho=407.828, visco=2)
        # self.assertEqual(round(st.k.mWmK, 2), 78.47)
        self.assertEqual(round(st.mu.muPas, 2), 48.62)
