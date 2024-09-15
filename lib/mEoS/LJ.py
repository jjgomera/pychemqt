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


from unittest import TestCase

from lib import unidades
from lib.meos import MEoS


class LJ(MEoS):
    """Multiparameter equation of state for Lennard-Jones fluid"""
    name = "Lennard-Jones Fluid"
    CASNumber = ""
    formula = ""
    synonym = ""
    _refPropName = ""
    _coolPropName = ""
    rhoc = unidades.Density(0.31)
    Tc = unidades.Temperature(1.32)
    Pc = unidades.Pressure(0.13006, "kPa")
    M = 1  # g/mol
    Tt = unidades.Temperature(0.6)
    Tb = unidades.Temperature(0.8)
    f_acent = 0
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")

    Fi1 = {"ao_log": [1, 1.5],
           "pow": [0, 1],
           "ao_pow": [6.262265814, -1.51515155]}

    thol = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for the Lennard-Jones fluid "
                    "of Thol (2016).",
        "__doi__": {"autor": "Thol, M., Rutkai, G., Köster, A., Lustig, R., "
                             "Span, R., Vrabec, J.",
                    "title": "Equation of State for the Lennard-Jones Fluid",
                    "ref": "J. Phys. Chem. Ref. Data 45(2) (2016) 023101",
                    "doi": "10.1063/1.4945000"},

        "R": 1,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": 0.661, "Tmax": 9, "Pmax": 65, "rhomax": 1.42,

        "nr1": [.005208073, 2.186252, -2.161016, 1.4527, -2.041792, .18695286],
        "d1": [4, 1, 1, 2, 2, 3],
        "t1": [1, 0.32, 0.505, 0.672, 0.843, 0.898],

        "nr2": [-0.090988445, -0.49745610, 0.10901431, -0.80055922,
                -0.568839, -0.6208625],
        "d2": [5, 2, 2, 3, 1, 1],
        "t2": [1.294, 2.59, 1.786, 2.77, 1.786, 1.205],
        "c2": [1, 2, 1, 2, 2, 1],
        "gamma2": [1]*6,

        "nr3": [-1.4667177, 1.891469, -0.1383701, -0.3869645, 0.1265702,
                0.605781, 1.179189, -0.47732679, -9.9218575, -0.5747932,
                0.003772923],
        "d3": [1, 1, 2, 3, 3, 2, 1, 2, 3, 1, 1],
        "t3": [2.83, 2.548, 4.65, 1.385, 1.46, 1.351, 0.66, 1.496, 1.83,
               1.616, 4.97],
        "alfa3": [2.067, 1.522, 8.82, 1.722, 0.679, 1.883, 3.925, 2.461, 28.2,
                  0.753, 0.82],
        "beta3": [0.625, 0.638, 3.91, 0.156, 0.157, 0.153, 1.16, 1.73, 383,
                  0.112, 0.119],
        "gamma3": [0.71, 0.86, 1.94, 1.48, 1.49, 1.945, 3.02, 1.11, 1.17,
                   1.33, 0.24],
        "epsilon3": [0.2053, 0.409, 0.6, 1.203, 1.829, 1.397, 1.39, 0.539,
                     0.934, 2.369, 2.43]}

    eq = (thol, )

    _vapor_Pressure = {
        "eq": 3,
        "n": [-5.4, 0.44704, -1.853, 0.1989, -1.125],
        "t": [1, 1.5, 4.7, 2.5, 21.4]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.1362e1, 0.2093e1, -0.2110e1, 0.3290, 0.1410e1],
        "t": [0.313, 0.94, 1.63, 17, 2.4]}
    _vapor_Density = {
        "eq": 2,
        "n": [-6.9655, -103.31, -2.0325, -44.481, -18.463, -260.70],
        "t": [1.32, 19.24, 0.36, 8.78, 4.04, 41.6]}


class Test(TestCase):
    """Testing"""

    def test_Thol(self):
        """Table 1, pag 1"""
        # FIXME: Speed of sound values fail
        st = LJ(T=0.8, rho=0.005)
        self.assertEqual(round(st.P, 7), 3.8430053)
        self.assertEqual(round(st.u.kJkg-st.u0.kJkg, 9), -0.054597389)
        self.assertEqual(round(st.cv.kJkgK-st.cv0.kJkgK, 9), 0.055672903)
        # self.assertEqual(round(st.w, 1), 1.1324263)
        self.assertEqual(round(st.a.kJkg, 7), 0.2776817)

        st = LJ(T=0.8, rho=0.8)
        self.assertEqual(round(st.P.kPa, 9), 0.015894013)
        self.assertEqual(round(st.u.kJkg-st.u0.kJkg, 7), -5.7174120)
        self.assertEqual(round(st.cv.kJkgK-st.cv0.kJkgK, 8), 0.95995160)
        # self.assertEqual(round(st.w, 7), 5.0522400)
        self.assertEqual(round(st.a.kJkg, 7), 1.1838093)

        st = LJ(T=1, rho=0.02)
        self.assertEqual(round(st.P.kPa, 9), 0.017886470)
        self.assertEqual(round(st.u.kJkg-st.u0.kJkg, 8), -0.18772644)
        self.assertEqual(round(st.cv.kJkgK-st.cv0.kJkgK, 8), 0.13016045)
        # self.assertEqual(round(st.w, 7), 1.2290934)
        self.assertEqual(round(st.a.kJkg, 7), 1.8318141)

        st = LJ(T=1, rho=0.71)
        self.assertEqual(round(st.P.kPa, 9), 0.075247483)
        self.assertEqual(round(st.u.kJkg-st.u0.kJkg, 7), -4.9564222)
        self.assertEqual(round(st.cv.kJkgK-st.cv0.kJkgK, 8), 0.68903536)
        # self.assertEqual(round(st.w, 7), 4.1644650)
        self.assertEqual(round(st.a.kJkg, 7), 2.9792860)

        st = LJ(T=2, rho=0.5)
        self.assertEqual(round(st.P.kPa, 7), 1.0751638)
        self.assertEqual(round(st.u.kJkg-st.u0.kJkg, 7), -3.1525021)
        self.assertEqual(round(st.cv.kJkgK-st.cv0.kJkgK, 8), 0.31068090)
        # self.assertEqual(round(st.w, 7), 3.5186329)
        self.assertEqual(round(st.a.kJkg, 7), 9.5274193)

        st = LJ(T=5, rho=0.6)
        self.assertEqual(round(st.P.kPa, 7), 6.9432008)
        self.assertEqual(round(st.u.kJkg-st.u0.kJkg, 7), -2.6956781)
        self.assertEqual(round(st.cv.kJkgK-st.cv0.kJkgK, 8), 0.31772707)
        # self.assertEqual(round(st.w, 7), 6.8375197)
        self.assertEqual(round(st.a.kJkg, 6), 26.122755)

        st = LJ(T=7, rho=1)
        self.assertEqual(round(st.P.kPa, 6), 41.531352)
        self.assertEqual(round(st.u.kJkg-st.u0.kJkg, 8), -0.62393078)
        self.assertEqual(round(st.cv.kJkgK-st.cv0.kJkgK, 8), 0.73348579)
        # self.assertEqual(round(st.w, 6), 14.201978)
        self.assertEqual(round(st.a.kJkg, 6), 48.074394)
