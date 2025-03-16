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
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


from unittest import TestCase

from lib import unidades
from lib.meos import MEoS


class D2(MEoS):
    """Multiparameter equation of state for deuterium"""
    # Using the corrected critical pressure and density from
    # Richardson, I.A., Leachman, J.W., Lemmon, E.W.",
    # Erratum: "Fundamental Equation of State for Deuterium" [J. Phys. Chem.
    # Ref. Data 43, 013103 (2014)",
    # 10.1063/1.5016519

    name = "deuterium"
    CASNumber = "7782-39-0"
    formula = "D2"
    synonym = ""
    _refPropName = "D2"
    _coolPropName = "Deuterium"
    rhoc = unidades.Density(69.405886)
    Tc = unidades.Temperature(38.34)
    Pc = unidades.Pressure(1679.6, "kPa")
    M = 4.0282  # g/mol
    Tt = unidades.Temperature(18.724)
    Tb = unidades.Temperature(23.661)
    f_acent = -0.136
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")

    Fi1 = {"ao_log": [1, 1.5],
           "pow": [0, 1],
           "ao_pow": [-2.0677351753, 2.4237151502],
           "ao_exp": [-3.54145, 3.0326, -3.52422, -1.73421, -3.57135, 2.14858,
                      6.23107, -3.30425, 6.23098, -3.57137, 3.32901, 0.97782],
           "titao": [7174.1/Tc, 8635/Tc, 902.7/Tc, 181.1/Tc, 438.5/Tc,
                     5034.2/Tc, 269.9/Tc, 229.9/Tc, 666.4/Tc, 452.8/Tc,
                     192.0/Tc, 1187.6/Tc]}

    CP1 = {"ao": 2.4512991,
           "an": [0.43563077e-02, -0.53169470e-03, 0.17067184e-04,
                  -0.53819932e-08, 0.89310438e-12],
           "pow": [1, 1.5, 2, 3, 4],
           "ao_exp": [0.18403263e2, -0.21257617e2, 0.41091635e1],
           "exp": [319, 361, 518]}

    richardson = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for deuterium of Richardson "
                    "et al. (2014).",
        "__doi__": {"autor": "Richardson, I.A., Leachman, J.W., Lemmon, E.W.",
                    "title": "Fundamental Equation of State for Deuterium",
                    "ref": "J. Phys. Chem. Ref. Data 43(1) (2014) 013103",
                    "doi": "10.1063/1.4864752"},

        "R": 8.3144621,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 2000000.0, "rhomax": 43.351,

        "nr1": [0.006267958, 10.53609, -10.14149, 0.356061, 0.1824472,
                -1.129638, -0.0549812, -0.6791329],
        "d1": [4, 1, 1, 2, 3, 1, 3, 2],
        "t1": [1, 0.462, 0.5584, 0.627, 1.201, 0.309, 1.314, 1.1166],

        "nr2": [1.347918, -0.8657582, 1.719146, -1.917977, 0.1233365,
                -0.07936891],
        "d2": [2, 2, 1, 1, 3, 2],
        "t2": [1.25, 1.25, 1.395, 1.627, 1.0, 2.5],
        "c2": [1, 1, 2, 2, 2, 2],
        "gamma2": [1]*6,

        "nr3": [1.686617, -4.240326, 1.857114, -0.5903705, 1.520171,
                2.361373, -2.297315],
        "d3": [1, 1, 2, 3, 3, 1, 3],
        "t3": [0.635, 0.664, 0.7082, 2.25, 1.524, 0.67, 0.709],
        "alfa3": [0.868, 0.636, 0.668, 0.65, 0.745, 0.782, 0.693],
        "beta3": [0.613, 0.584, 0.57, 1.056, 1.01, 1.025, 1.029],
        "gamma3": [0.6306, 0.711, 0.6446, 0.8226, 0.992, 1.2184, 1.203],
        "epsilon3": [1.46, 1.7864, 1.647, 0.541, 0.969, 1.892, 1.076],
        "nr4": []}

    mccarty = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for deuterium of McCarty (1989)",
        "__doi__": {"autor": "McCarty, R.D.",
                    "title": "Correlations for the Thermophysical Properties "
                             "of Deuterium",
                    "ref": "NIST, Boulder, CO, 1989",
                    "doi": ""},

        "R": 8.31434,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 423.0, "Pmax": 320000.0, "rhomax": 43.38,

        "b": [None, 0.4894244053982e-4, 0.5600164604601e-1, -0.6301493491211,
              0.2538329946038e1, 0.1723475985309e3, 0.2956238369436e-4,
              -0.3926317169317e-2, 0.1195764193293e-1, 0.1136916678824e5,
              -0.1916378195727e-6, 0.3153535946452e-3, 0.2122937335070e-1,
              -0.1057999371607e-5, -0.6722062598854e-4, -0.3030166828627,
              0.1980817195099e-5, -0.1453922641871e-7, 0.1780919116891e-3,
              -0.1823145348424e-5, -0.1135358616578e5, -0.1943542941899e4,
              -0.3632847669580e2, 0.1087745118380e3, -0.4078276062687e-1,
              0.6460021864005e-2, -0.4480242189217e-4, -0.2475011206216e-3,
              -0.8834384656760e-8, -0.1081622159862e-8, -0.1478159334303e-10,
              0.7926922356112e-11, 0.5721547329378e-11]}

    eq = richardson, mccarty

    _surface = {"sigma": [0.009376], "exp": [1.258]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-5.5706, 1.7631, -0.5458, 1.2154, -1.1556],
        "t": [1., 1.5, 2.83, 4.06, 5.4]}
    _liquid_Density = {
        "eq": 1,
        "n": [3.3769, -5.3693, 11.943, -17.361, 15.170, -6.3079],
        "t": [0.512, 1.12, 1.8, 2.55, 3.4, 4.4]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.8111, -7.3624, 2.2294, -21.443, 12.796, -31.334],
        "t": [0.528, 2.03, 3.6, 5.0, 6.5, 9.0]}


class pD2(D2):
    """Multiparameter equation of state for paradeuterium"""
    name = "paradeuterium"
    formula = "pD2"
    _refPropName = ""
    _coolPropName = "ParaDeuterium"

    Fi1 = {"ao_log": [1, 1.5],
           "pow": [0, 1],
           "ao_pow": [-2.0683998716, 2.4241000701],
           "ao_exp": [1.28527, 1.11376, -2.49100, 6.38763, 6.17406, -3.13698,
                      -3.14254, -2.29511, -3.37000, 1.13634, 0.72512],
           "titao": [5068/D2.Tc, 1000.8/D2.Tc, 261.5/D2.Tc, 437.2/D2.Tc,
                     312.3/D2.Tc, 382.8/D2.Tc, 356.8/D2.Tc, 294.7/D2.Tc,
                     682.4/D2.Tc, 246/D2.Tc, 277.1/D2.Tc]}

    richardson = D2.richardson.copy()
    richardson["cp"] = Fi1
    eq = (richardson, )


class oD2(D2):
    """Multiparameter equation of state for orthodeuterium"""
    name = "orthodeuterium"
    formula = "oD2"
    _refPropName = ""
    _coolPropName = "OrthoDeuterium"

    Fi1 = {"ao_log": [1, 1.5],
           "pow": [0, 1],
           "ao_pow": [-2.0672670563, 2.4234599781],
           "ao_exp": [4.04482, -4.65391, -4.65342, 3.46313, -4.58637, -4.65030,
                      -4.65124, 2.67024, 15.20455, 0.87164, -4.76080, 4.32447],
           "titao": [1591/D2.Tc, 481.6/D2.Tc, 472.4/D2.Tc, 362.2/D2.Tc,
                     2038/D2.Tc, 463.2/D2.Tc, 491.3/D2.Tc, 2713.4/D2.Tc,
                     618.6/D2.Tc, 8642/D2.Tc, 961.7/D2.Tc, 253.2/D2.Tc]}

    richardson = D2.richardson.copy()
    richardson["cp"] = Fi1
    eq = (richardson, )


class Test(TestCase):
    """Testing"""
    def test_richardson(self):
        """Selected point from Table 7 pag 12, saturation state"""
        st = D2(T=18.724, x=0.5)
        self.assertEqual(round(st.P.kPa, 3), 17.189)
        self.assertEqual(round(st.Liquido.rho, 2), 174.63)
        self.assertEqual(round(st.Gas.rho, 3), 0.455)
        self.assertEqual(round(st.Liquido.h.kJkg, 3), -30.450)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 286.160)
        self.assertEqual(round(st.Liquido.s.kJkgK, 3), -1.415)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), 15.494)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 3.355)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 3.143)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 5.626)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 5.364)
        self.assertEqual(round(st.Liquido.w, 1), 1085.50)
        self.assertEqual(round(st.Gas.w, 2), 250.95)

        st = D2(T=25, x=0.5)
        self.assertEqual(round(st.P.kPa, 2), 146.40)
        self.assertEqual(round(st.Liquido.rho, 3), 158.870)
        self.assertEqual(round(st.Gas.rho, 3), 3.144)
        self.assertEqual(round(st.Liquido.h.kJkg, 3), 9.282)
        self.assertEqual(round(st.Gas.h.kJkg, 3), 307.750)
        self.assertEqual(round(st.Liquido.s.kJkgK, 3), 0.370)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), 12.309)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 3.442)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 3.318)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 6.996)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 6.197)
        self.assertEqual(round(st.Liquido.w, 2), 937.27)
        self.assertEqual(round(st.Gas.w, 2), 278.61)

        st = D2(T=30, x=0.5)
        self.assertEqual(round(st.P.kPa, 2), 445.75)
        self.assertEqual(round(st.Liquido.rho, 2), 143.160)
        self.assertEqual(round(st.Gas.rho, 3), 9.055)
        self.assertEqual(round(st.Liquido.h.kJkg, 3), 49.449)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 313.220)
        self.assertEqual(round(st.Liquido.s.kJkgK, 3), 1.757)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), 10.550)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 3.576)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 3.568)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 9.036)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 8.013)
        self.assertEqual(round(st.Liquido.w, 2), 777.76)
        self.assertEqual(round(st.Gas.w, 2), 288.11)

        st = D2(T=35, x=0.5)
        self.assertEqual(round(st.P.kPa, 1), 1036.7)
        self.assertEqual(round(st.Liquido.rho, 2), 120.480)
        self.assertEqual(round(st.Gas.rho, 3), 23.053)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 105.180)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 297.410)
        self.assertEqual(round(st.Liquido.s.kJkgK, 3), 3.329)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), 8.821)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 3.846)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 3.954)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 16.103)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 15.430)
        self.assertEqual(round(st.Liquido.w, 2), 561.44)
        self.assertEqual(round(st.Gas.w, 2), 286.83)

        st = D2(T=38, x=0.5)
        self.assertEqual(round(st.P.kPa, 1), 1600.6)
        self.assertEqual(round(st.Liquido.rho, 3), 88.550)
        self.assertEqual(round(st.Gas.rho, 3), 50.659)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 172.270)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 244.900)
        self.assertEqual(round(st.Liquido.s.kJkgK, 3), 5.009)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), 6.920)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 4.187)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 4.300)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 2), 132.060)
        self.assertEqual(round(st.Gas.cp.kJkgK, 2), 137.920)
        self.assertEqual(round(st.Liquido.w, 2), 374.12)
        self.assertEqual(round(st.Gas.w, 2), 297.24)
