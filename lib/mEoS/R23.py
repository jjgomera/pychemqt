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

from scipy import exp

from lib import unidades
from lib.meos import MEoS


class R23(MEoS):
    """Multiparameter equation of state for R23"""
    name = "trifluoromethane"
    CASNumber = "75-46-7"
    formula = "CHF3"
    synonym = "R23"
    _refPropName = "R23"
    _coolPropName = "R23"
    rhoc = unidades.Density(526.504152)
    Tc = unidades.Temperature(299.293)
    Pc = unidades.Pressure(4832.0, "kPa")
    M = 70.01385  # g/mol
    Tt = unidades.Temperature(118.02)
    Tb = unidades.Temperature(191.132)
    f_acent = 0.263
    momentoDipolar = unidades.DipoleMoment(1.649, "Debye")
    id = 643

    Fi1 = {"ao_log": [1, 2.999],
           "pow": [0, 1],
           "ao_pow": [-8.31386064, 6.55087253],
           "ao_exp": [2.371, 3.237, 2.61, 0.8274],
           "titao": [744/Tc, 1459/Tc, 2135/Tc, 4911/Tc]}

    CP1 = {"ao": 3.999509244,
           "an": [], "pow": [],
           "ao_exp": [1.070326018, 1.566866769, 0.848051597, 1.847243699,
                      1.649657530, 2.043965290],
           "exp": [4368.102594, 1607.104940, 1007.138279, 1973.991027,
                   1657.461854, 729.455868]}

    CP2 = {"ao": 4.0101431,
           "an": [-.55274742e-2, .74008258e-4, -.12590943e-6, .69472178e-10],
           "pow": [1, 2, 3, 4],
           "ao_exp": [], "exp": []}

    penoncello = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R23 of Penoncello (2003)",
        "__doi__": {"autor": "Penoncello, S.G., Lemmon, E.W., Jacobsen, R.T, "
                             "Shan, Z.",
                    "title": "A Fundamental Equation for Trifluoromethane "
                             "(R-23)",
                    "ref": "J. Phys. Chem. Ref. Data 32(4) (2003) 1473-1499",
                    "doi":  "10.1063/1.1559671"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 475.0, "Pmax": 120000.0, "rhomax": 24.31,

        "nr1": [.7041529e1, -.8259512e1, .805304e-2, -.8617615e-1, .633341e-2],
        "d1": [1, 1, 1, 2, 5],
        "t1": [0.744, 0.94, 4.3, 1.46, 0.68],

        "nr2": [-0.1863285, 0.3280510, 0.5191023, 0.6916144e-1, -0.5045875e-2,
                -0.1744221e-1, -0.5003972e-1, 0.4729813e-1, -0.6164031e-1,
                0.1583585e-1, -0.1795790e-2, -0.1099007e-2],
        "d2": [1, 2, 3, 5, 1, 2, 2, 4, 4, 4, 2, 2],
        "t2": [4.8, 1.5, 2.07, 0.09, 9.6, 0.19, 11.2, 0.27, 1.6, 10.3, 14, 15],
        "c2": [1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 4],
        "gamma2": [1]*12}

    penoncello2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R23 of Penoncello (2000)",
        "__doi__": {"autor": "Penoncello, S.G., Shan, Z., Jacobsen, R.T.",
                    "title": "A Fundamental Equation for the Calculation of "
                             "the Thermodynamic Properties of "
                             "Trifluoromethane (R23)",
                    "ref": "ASHRAE Trans. 106(Part 1) (2000) 739-756",
                    "doi":  ""},

        "R": 8.31451,
        "cp": CP1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 473.15, "Pmax": 120000.0, "rhomax": 23.0,

        "nr1": [0.350093635099, -0.131185838025e1, -0.254118065769,
                .104275296122, -.205326997924, .256040993750, .118078220087e-1,
                0.532850915621e-3,  0.956700157221e-3, -0.118990410423e-5],
        "d1": [1, 1, 1, 2, 2, 2, 3, 4, 6, 8],
        "t1": [-0.14, 1.49, 2.41, 0.05, 1.59, 2.04, -0.27, 2.76, -0.06, 3.25],

        "nr2": [-0.180609172794, 0.138077199166, 0.507828500811e-1,
                0.439772083175e-1, -0.723557234469e-1, 0.256500006055e-2,
                0.263213487134e-1, 0.139266509424e-1, -0.105325247813e-1,
                0.136475671500e-2, -0.592653649931e-2, -0.644925101471e-1,
                -0.227635186710e-1, 0.122367812706, 0.318153208563e-1,
                0.146725272055e-1, -0.923639585566e-1],
        "d2": [1, 2, 2, 3, 3, 6, 6, 7, 7, 10, 2, 3, 4, 4, 5, 5, 5],
        "t2": [5.36, 5.28, 4.23, 3.35, 6.93, 8.48, 6.01, 3.34, 7.1, 5.46,
               16.06, 19.37, 10.81, 22.79, 34.95, 9.94, 29.16],
        "c2": [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4],
        "gamma2": [1]*17}

    platzer = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R23 of Platzer (1990)",
        "__doi__": {"autor": "Platzer, B., Polt, A., Maurer, G.",
                    "title": "Thermophysical Properties of Refrigerants",
                    "ref": "Berlin:  Springer-Verlag, 1990.",
                    "doi": ""},

        "R": 8.31451,
        "cp": CP2,
        "ref": "NBP",

        "Tmin": 90.0, "Tmax": 475.0, "Pmax": 60000.0, "rhomax": 16.65,

        "nr1": [-0.133234251368e1, 0.210373595421e1, -0.376198728030,
                0.881622087335, -0.272053790906e1, 0.247468024356e1,
                -0.234010064393e1, 0.303959507238, 0.317372750273e-1,
                .329392142221e-1, .20583853186, .133550139894, -.181698216766,
                -0.245123269882e-1, 0.247477874180e-1, 0.589916583383e-2],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.133234251368e1, -0.210373595421e1, 0.376198728030,
                0.574267667948, -0.762218931280, 0.472710395636e-1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [0.70304082]*6}

    eq = penoncello, penoncello2, platzer

    _surface = {"sigma": [-0.32359, 0.37702], "exp": [1.6055, 1.5232]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.2631,  1.3140, -0.78507, -3.1991],
        "t": [1, 1.5, 2.4, 3.9]}
    _liquid_Density = {
        "eq": 1,
        "n": [2.2636, 0.47007, 0.28660],
        "t": [0.37, 0.94, 3.1]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.5136, -7.7491, -24.871, -65.637],
        "t": [0.43, 1.4, 3.7, 8]}

    visco0 = {"__name__": "Shan (2000)",
              "__doi__": {
                  "autor": "Shan, Z., Penoncello, S.G., Jacobsen, R.T.",
                  "title": "A Generalized Model for Viscosity and Thermal "
                           "Conductivity of Trifluoromethane (R-23)",
                  "ref": "ASHRAE Trans. 106(Part 1) (2000) 757-767",
                  "doi": ""},

              "eq": 0, "omega": 1,
              "method": "_visco0",

              "ek": 243.91, "sigma": 0.4278,
              "n_chapman": 0.2233755/M**0.5,
              "collision": [0.4425728, -0.5138403, 0.1547566, -0.02821844,
                            0.001578286]}

    def _visco0(self, rho, T, fase):
        rhol = 32.174
        C1 = 1.3163
        C2 = 0.1832
        deltaG = 771.23
        nmax = 3.967
        R = 8.31451/self.M

        Drho = rhol-rho/self.M
        delta = rho/self.M-7.5114
        tau = T-299.28

        muo = self._Visco0(T)
        mug = muo*(Drho/rhol)**C1
        mur = (rho/self.M/rhol)**C1*C2*rhol**2/Drho*T**0.5*exp(
                rho/self.M/Drho*deltaG/R/self.M/T)
        muc = 4*nmax/(exp(delta)+exp(-delta))/(exp(tau)+exp(-tau))
        return unidades.Viscosity(mug+mur+muc, "muPas")

    _viscosity = visco0,

    thermo0 = {"eq": 0,
               "method": "_thermo0",
               "__name__": "Shan (2000)",
               "__doi__": {
                  "autor": "Shan, Z., Penoncello, S.G., Jacobsen, R.T.",
                  "title": "A Generalized Model for Viscosity and Thermal "
                           "Conductivity of Trifluoromethane (R-23)",
                  "ref": "ASHRAE Trans. 106(Part 1) (2000) 757-767",
                  "doi": ""}}

    def _thermo0(self, rho, T, fase):
        rhol = 68.345
        B1 = -2.5370
        B2 = 0.05366
        C1 = 0.94215
        C2 = 0.14914
        deltaG = 2508.58
        lmax = 25.

        Drho = rhol-rho/self.M
        delta = rho/self.M-7.5114
        tau = T-299.28

        lg = (B1+B2*T)*(Drho/rhol)**C1
        lr = (rho/self.M/rhol)**C1*C2*rhol**2/Drho*T**0.5*exp(
                rho/self.M/Drho*deltaG/self.R.kJkgK/self.M/T)
        lc = 4*lmax/(exp(delta)+exp(-delta))/(exp(tau)+exp(-tau))
        return unidades.ThermalConductivity(lg+lr+lc, "mWmK")

    _thermal = thermo0,


class Test(TestCase):

    def test_penoncello(self):
        # Appendix pag 1489, saturation state
        st = R23(T=R23.Tt, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 0.00006)
        self.assertEqual(round(st.Liquido.rho, 1), 1701.9)
        self.assertEqual(round(st.Liquido.h.kJkg, 4), -3.6734)
        self.assertEqual(round(st.Liquido.s.kJkgK, 5), -0.07052)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 0.8155)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 1.221)
        self.assertEqual(round(st.Liquido.w, 1), 1211.2)
        self.assertEqual(round(st.Gas.rho, 5), 0.00414)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 289.21)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 2.4111)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 0.3802)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 0.4998)
        self.assertEqual(round(st.Gas.w, 1), 135.7)

        st = R23(T=-125+273.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 0.00355)
        self.assertEqual(round(st.Liquido.rho, 1), 1600.7)
        self.assertEqual(round(st.Liquido.h.kJkg, 3), 32.461)
        self.assertEqual(round(st.Liquido.s.kJkgK, 5), 0.20222)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 0.7238)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 1.197)
        self.assertEqual(round(st.Liquido.w, 1), 1062.7)
        self.assertEqual(round(st.Gas.rho, 5), 0.20309)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 304.01)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 2.0352)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 0.4269)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 0.5527)
        self.assertEqual(round(st.Gas.w, 1), 150.1)

        st = R23(T=-100+273.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 0.03157)
        self.assertEqual(round(st.Liquido.rho, 1), 1512.3)
        self.assertEqual(round(st.Liquido.h.kJkg, 3), 62.561)
        self.assertEqual(round(st.Liquido.s.kJkgK, 5), 0.38982)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 0.7152)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 1.213)
        self.assertEqual(round(st.Liquido.w, 1), 920.7)
        self.assertEqual(round(st.Gas.rho, 4), 1.5672)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 315.98)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.8534)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 0.4842)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 0.6246)
        self.assertEqual(round(st.Gas.w, 1), 159.5)

        st = R23(T=-50+273.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 0.47893)
        self.assertEqual(round(st.Liquido.rho, 1), 1315.3)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 125.53)
        self.assertEqual(round(st.Liquido.s.kJkgK, 5), 0.70712)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 0.7254)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 1.325)
        self.assertEqual(round(st.Liquido.w, 1), 647.9)
        self.assertEqual(round(st.Gas.rho, 3), 20.430)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 335.52)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.6482)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 0.6329)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 0.8701)
        self.assertEqual(round(st.Gas.w, 1), 167.3)

        st = R23(T=273.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 2.49469)
        self.assertEqual(round(st.Liquido.rho, 1), 1035.1)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 200.00)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.0000)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 0.7844)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 1.853)
        self.assertEqual(round(st.Liquido.w, 1), 345.4)
        self.assertEqual(round(st.Gas.rho, 2), 118.67)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 337.64)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.5039)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 0.8324)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 1.715)
        self.assertEqual(round(st.Gas.w, 1), 151.0)

        st = R23(T=25+273.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 4.69863)
        self.assertEqual(round(st.Liquido.rho, 2), 680.09)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 261.94)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.2067)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 0.9429)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 2), 18.87)
        self.assertEqual(round(st.Liquido.w, 1), 140.8)
        self.assertEqual(round(st.Gas.rho, 2), 379.91)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 301.55)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.3396)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 1.029)
        self.assertEqual(round(st.Gas.cp.kJkgK, 2), 26.84)
        self.assertEqual(round(st.Gas.w, 1), 124.3)

        # Pressure input, pag 1492
        st = R23(P=1e3, x=0.5)
        self.assertEqual(round(st.T.C, 3), -136.057)
        self.assertEqual(round(st.Liquido.rho, 1), 1638.6)
        self.assertEqual(round(st.Liquido.h.kJkg, 3), 19.242)
        self.assertEqual(round(st.Liquido.s.kJkgK, 5), 0.10950)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 0.7336)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 1.195)
        self.assertEqual(round(st.Liquido.w, 1), 1128.4)
        self.assertEqual(round(st.Gas.rho, 5), 0.06158)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 298.59)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 2.1471)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 0.4065)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 0.5290)
        self.assertEqual(round(st.Gas.w, 1), 145.2)

        st = R23(P=4.8e6, x=0.5)
        self.assertEqual(round(st.T.C, 3), 25.875)
        self.assertEqual(round(st.Liquido.rho, 2), 611.28)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 270.48)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.2348)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 0.9848)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 2), 93.87)
        self.assertEqual(round(st.Liquido.w, 1), 128.0)
        self.assertEqual(round(st.Gas.rho, 2), 444.22)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 292.12)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.3071)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 1.040)
        self.assertEqual(round(st.Gas.cp.kJkgK, 1), 119.1)
        self.assertEqual(round(st.Gas.w, 1), 122.0)

        # Selected points from pag 1494, single phase region
        st = R23(T=-150+273.15, P=1e5)
        self.assertEqual(round(st.rho, 1), 1685.4)
        self.assertEqual(round(st.u.kJkg, 4), 2.5289)
        self.assertEqual(round(st.h.kJkg, 4), 2.5882)
        self.assertEqual(round(st.s.kJkgK, 5), -0.01908)
        self.assertEqual(round(st.cv.kJkgK, 4), 0.7741)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.205)
        self.assertEqual(round(st.w, 1), 1200.5)

        st = R23(T=-60+273.15, P=2e5)
        self.assertEqual(round(st.rho, 4), 8.3497)
        self.assertEqual(round(st.u.kJkg, 2), 311.94)
        self.assertEqual(round(st.h.kJkg, 2), 335.89)
        self.assertEqual(round(st.s.kJkgK, 4), 1.7452)
        self.assertEqual(round(st.cv.kJkgK, 4), 0.5515)
        self.assertEqual(round(st.cp.kJkgK, 4), 0.7160)
        self.assertEqual(round(st.w, 1), 171.2)

        st = R23(T=-50+273.15, P=5e5)
        self.assertEqual(round(st.rho, 1), 1315.4)
        self.assertEqual(round(st.u.kJkg, 2), 125.15)
        self.assertEqual(round(st.h.kJkg, 2), 125.53)
        self.assertEqual(round(st.s.kJkgK, 5), 0.70706)
        self.assertEqual(round(st.cv.kJkgK, 4), 0.7254)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.325)
        self.assertEqual(round(st.w, 1), 648.1)

        st = R23(T=200+273.15, P=1e6)
        self.assertEqual(round(st.rho, 3), 18.027)
        self.assertEqual(round(st.u.kJkg, 2), 487.77)
        self.assertEqual(round(st.h.kJkg, 2), 543.24)
        self.assertEqual(round(st.s.kJkgK, 4), 2.1784)
        self.assertEqual(round(st.cv.kJkgK, 4), 0.8453)
        self.assertEqual(round(st.cp.kJkgK, 4), 0.9758)
        self.assertEqual(round(st.w, 1), 251.5)

        st = R23(T=100+273.15, P=2e6)
        self.assertEqual(round(st.rho, 3), 48.524)
        self.assertEqual(round(st.u.kJkg, 2), 403.75)
        self.assertEqual(round(st.h.kJkg, 2), 444.97)
        self.assertEqual(round(st.s.kJkgK, 4), 1.8656)
        self.assertEqual(round(st.cv.kJkgK, 4), 0.7370)
        self.assertEqual(round(st.cp.kJkgK, 4), 0.9096)
        self.assertEqual(round(st.w, 1), 217.4)

        st = R23(T=273.15, P=5e6)
        self.assertEqual(round(st.rho, 1), 1075.3)
        self.assertEqual(round(st.u.kJkg, 2), 192.86)
        self.assertEqual(round(st.h.kJkg, 2), 197.51)
        self.assertEqual(round(st.s.kJkgK, 5), 0.98220)
        self.assertEqual(round(st.cv.kJkgK, 4), 0.7692)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.643)
        self.assertEqual(round(st.w, 1), 400.2)

        st = R23(T=20+273.15, P=1e7)
        self.assertEqual(round(st.rho, 1), 1008.2)
        self.assertEqual(round(st.u.kJkg, 2), 216.36)
        self.assertEqual(round(st.h.kJkg, 2), 226.28)
        self.assertEqual(round(st.s.kJkgK, 4), 1.0669)
        self.assertEqual(round(st.cv.kJkgK, 4), 0.7807)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.653)
        self.assertEqual(round(st.w, 1), 372.3)

        st = R23(T=200+273.15, P=2e7)
        self.assertEqual(round(st.rho, 2), 412.39)
        self.assertEqual(round(st.u.kJkg, 2), 436.90)
        self.assertEqual(round(st.h.kJkg, 2), 485.40)
        self.assertEqual(round(st.s.kJkgK, 4), 1.7237)
        self.assertEqual(round(st.cv.kJkgK, 4), 0.9049)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.289)
        self.assertEqual(round(st.w, 1), 268.6)

        st = R23(T=-20+273.15, P=5e7)
        self.assertEqual(round(st.rho, 1), 1365.7)
        self.assertEqual(round(st.u.kJkg, 2), 140.36)
        self.assertEqual(round(st.h.kJkg, 2), 176.97)
        self.assertEqual(round(st.s.kJkgK, 5), 0.76824)
        self.assertEqual(round(st.cv.kJkgK, 4), 0.7453)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.196)
        self.assertEqual(round(st.w, 1), 816.2)

        st = R23(T=200+273.15, P=1e8)
        self.assertEqual(round(st.rho, 1), 1042.4)
        self.assertEqual(round(st.u.kJkg, 2), 367.03)
        self.assertEqual(round(st.h.kJkg, 2), 462.96)
        self.assertEqual(round(st.s.kJkgK, 4), 1.4595)
        self.assertEqual(round(st.cv.kJkgK, 4), 0.9594)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.265)
        self.assertEqual(round(st.w, 1), 633.7)
