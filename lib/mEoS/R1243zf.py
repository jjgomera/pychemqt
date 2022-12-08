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


class R1243zf(MEoS):
    """Multiparameter equation of state for R1243zf"""
    name = "3,3,3-trifluoropropene"
    CASNumber = "677-21-4"
    formula = "CH2=CHCF3"
    synonym = "R-1243zf"
    rhoc = unidades.Density(413.019859)
    Tc = unidades.Temperature(376.93)
    Pc = unidades.Pressure(3517.9, "kPa")
    M = 96.05113  # g/mol
    Tt = unidades.Temperature(0)
    Tb = unidades.Temperature(247.726)
    f_acent = 0.26
    momentoDipolar = unidades.DipoleMoment(0, "Debye")

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1], "ao_pow": [-12.06546081, 8.207216743],
           "ao_exp": [11.247, 8.1391], "titao": [789/Tc, 2196/Tc]}

    akasaka = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R1243zf of Akasaka "
                    "(2016)",
        "__doi__": {"autor": "Akasaka, R., Lemmon, E.W.",
                    "title": "Fundamental Equations of State for cis-1,3,3,3-"
                             "Tetrafluoropropene [R-1234ze(Z)] and 3,3,3-"
                             "Trifluoropropene (R-1243zf)",
                    "ref": "J. Chem. Eng. Data 64(11) (2019) 4679-4691",
                    "doi": "10.1021/acs.jced.9b00007"},

        "R": 8.3144598,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 420.0, "Pmax": 20000.0, "rhomax": 13.26,

        "nr1": [.0383741, 1.551307798, -1.918093187, -.536888401, .116850286],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.43, 1, 1, 0.3],

        "nr2": [-2.004856821, -1.589616187, 0.351405821, -1.110707332,
                -0.014517895],
        "d2": [1, 3, 2, 2, 7],
        "t2": [1.95, 2.19, 0.77, 2.81, 1.03],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [2.521451755, 1.410134214, -0.787246964, -1.103521137,
                -0.692133362],
        "d3": [1, 1, 3, 2, 2],
        "t3": [1.48, 1.12, 1.62, 1.17, 1.6],
        "alfa3": [1.067, 1.498, 1.006, 1.41, 0.98],
        "beta3": [0.698, 2.95, 0.697, 2.94, 0.57],
        "gamma3": [1.35, 1.2, 1.286, 1.18, 1],
        "epsilon3": [0.998, 1.06, 0.655, 0.861, 0.61]}

    eq = (akasaka,)

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.3009, 1.5186, -2.5303, -1.5139],
        "t": [1.0, 1.5, 2.8, 5]}
    _liquid_Density = {
        "eq": 1,
        "n": [2.3870, -1.2213, 4.8105, -5.7706],
        "t": [0.362, 0.83, 1.31, 1.82, 2.4]}
    _vapor_Density = {
        "eq": 2,
        "n": [-2.9588, -7.1173, -21.945, -54.068, -128.95],
        "t": [0.385, 1.25, 3.4, 7.25, 16]}


class Test(TestCase):
    """Test class"""

    def test_akasaka(self):
        """Table 6, pag 4690"""
        st = R1243zf(T=250, rhom=0)
        self.assertEqual(round(st.P.MPa, 6), 0)
        self.assertEqual(round(st.cvM.kJmolK, 7), 0.0690267)
        self.assertEqual(round(st.cpM.kJmolK, 7), 0.0773412)
        self.assertEqual(round(st.w, 3), 155.716)

        st = R1243zf(T=250, rhom=12)
        self.assertEqual(round(st.P.MPa, 5), 21.70009)
        self.assertEqual(round(st.cvM.kJmolK, 7), 0.0862865)
        self.assertEqual(round(st.cpM.kJmolK, 7), 0.1181422)
        self.assertEqual(round(st.w, 3), 852.712)

        st = R1243zf(T=250, rhom=0.05)
        self.assertEqual(round(st.P.MPa, 8), 0.09953467)
        self.assertEqual(round(st.cvM.kJmolK, 7), 0.0707756)
        self.assertEqual(round(st.cpM.kJmolK, 7), 0.0811839)
        self.assertEqual(round(st.w, 3), 150.725)

        st = R1243zf(T=350, rhom=9)
        self.assertEqual(round(st.P.MPa, 6), 7.656805)
        self.assertEqual(round(st.cvM.kJmolK, 6), 0.105027)
        self.assertEqual(round(st.cpM.kJmolK, 6), 0.157970)
        self.assertEqual(round(st.w, 3), 392.414)

        st = R1243zf(T=350, rhom=1)
        self.assertEqual(round(st.P.MPa, 6), 1.939182)
        self.assertEqual(round(st.cvM.kJmolK, 6), 0.105867)
        self.assertEqual(round(st.cpM.kJmolK, 6), 0.158364)
        self.assertEqual(round(st.w, 3), 129.425)

        st = R1243zf(T=377, rhom=4.3)
        self.assertEqual(round(st.P.MPa, 6), 3.522701)
        self.assertEqual(round(st.cvM.kJmolK, 6), 0.129027)
        self.assertEqual(round(st.cpM.kJmolK, 4), 89.7236)
        self.assertEqual(round(st.w, 4), 88.6451)
