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

from lib import unidades
from lib.meos import MEoS


class Cis_2_butene(MEoS):
    """Multiparameter equation of state for cis-butene"""
    name = "cis-butene"
    CASNumber = "590-18-1"
    formula = "CH3-CH=CH-CH3"
    synonym = ""
    _refPropName = "C2BUTENE"
    _coolPropName = "cis-2-Butene"
    rhoc = unidades.Density(238.11522208)
    Tc = unidades.Temperature(435.75)
    Pc = unidades.Pressure(4225.5, "kPa")
    M = 56.10632  # g/mol
    Tt = unidades.Temperature(134.3)
    Tb = unidades.Temperature(276.87)
    f_acent = 0.202
    momentoDipolar = unidades.DipoleMoment(0.3, "Debye")
    id = 25

    Fi1 = {"ao_log": [1, 2.9687],
           "pow": [0, 1],
           "ao_pow": [0.2591542, 2.4189888],
           "ao_exp": [3.2375, 7.0437, 11.414, 7.3722],
           "titao": [248/Tc, 1183/Tc, 2092/Tc, 4397/Tc]}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for 1-butene of Lemmon "
                    "and Ihmels (2005)",
        "__doi__": {"autor": "Lemmon, E.W., Ihmels, E.C.",
                    "title": "Thermodynamic properties of the butenes: Part "
                             "II. Short fundamental equations of state",
                    "ref": "Fluid Phase Equilibria 228-229 (2005) 173-187",
                    "doi":  "10.1016/j.fluid.2004.09.004"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 525., "Pmax": 50000.0, "rhomax": 14.09,
        "Pmin": 0.00026, "rhomin": 14.09,

        "nr1":  [0.77827, -2.8064, 1.003, 0.013762, 0.085514, 0.00021268],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.12, 1.3, 1.74, 2.1, 0.28, 0.69],

        "nr2": [0.22962, -0.072442, -0.23722, -0.074071, -0.026547, 0.012032],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.75, 2., 4.4, 4.7, 15., 14.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = lemmon,

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.70022e1, 0.13695e1, -0.30509e1, 0.10012, -0.15577e1],
        "t": [1.0, 1.5, 3.2, 3.46, 6.4]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.46849e1, -0.54614e1, 0.34718e1, 0.50511e1, -0.50389e1],
        "t": [0.402, 0.54, 0.69, 6.6, 7.0]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.28918e1, -0.58582e1, -0.17443e2, -0.24566e2, -0.29413e2,
               -0.11392e3],
        "t": [0.4098, 1.174, 3.11, 6.1, 7.6, 14.8]}


class Test(TestCase):
    def test_shortLemmon(self):
        # Table 9, Pag 186
        st = Cis_2_butene(T=350, rhom=0)
        self.assertEqual(round(st.P.MPa, 4), 0)
        self.assertEqual(round(st.hM.kJkmol, 0), 29735)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 83.593)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 91.907)
        self.assertEqual(round(st.w, 2), 238.80)

        st = Cis_2_butene(T=350, rhom=0.3)
        self.assertEqual(round(st.P.MPa, 5), 0.74661)
        self.assertEqual(round(st.hM.kJkmol, 0), 28294)
        self.assertEqual(round(st.sM.kJkmolK, 3), 84.888)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 89.258)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 106.34)
        self.assertEqual(round(st.w, 2), 209.87)

        st = Cis_2_butene(T=350, rhom=10)
        self.assertEqual(round(st.P.MPa, 4), 5.8051)
        self.assertEqual(round(st.hM.kJkmol, 1), 9632.7)
        self.assertEqual(round(st.sM.kJkmolK, 3), 29.100)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 93.935)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 138.23)
        self.assertEqual(round(st.w, 2), 770.10)

        st = Cis_2_butene(T=440, rhom=4)
        self.assertEqual(round(st.P.MPa, 4), 4.5067)
        self.assertEqual(round(st.hM.kJkmol, 0), 28321)
        self.assertEqual(round(st.sM.kJkmolK, 3), 75.755)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 126.00)
        self.assertEqual(round(st.cpM.kJkmolK, 1), 1686.3)
        self.assertEqual(round(st.w, 2), 130.71)
