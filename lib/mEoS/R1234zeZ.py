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


class R1234zeZ(MEoS):
    """Multiparameter equation of state for R1234ze(Z)"""
    name = "cis-1,3,3,3-tetrafluoropropene"
    CASNumber = "29118-25-0"
    formula = "CHF=CHCF3"
    synonym = "R-1234zeZ"
    rhoc = unidades.Density(456.1664)
    Tc = unidades.Temperature(423.27)
    Pc = unidades.Pressure(3530.6, "kPa")
    M = 114.0416  # g/mol
    Tt = unidades.Temperature(0)
    Tb = unidades.Temperature(282.878)
    f_acent = 0.327
    momentoDipolar = unidades.DipoleMoment(0, "Debye")

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1], "ao_pow": [-2.422442259, 8.1190539844],
           "ao_exp": [4.2365, 13.063], "titao": [20/Tc, 1335/Tc]}

    CP1 = {"ao": -1.6994,
           "an": [24.527/Tc, -9.9249/Tc**2, -1.5158/Tc**3],
           "pow": [1, 2, 3]}

    akasaka = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R1234ze(Z) of Akasaka "
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

        "Tmin": 238, "Tmax": 440.0, "Pmax": 34000.0,

        "nr1": [0.03194509, 1.394592, -2.300799, -0.2556693, 0.1282934],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.333, 1, 1, 0.38],

        "nr2": [-1.335381, -1.366494, 0.2004912, -0.6489709, -0.02220033],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.85, 3.16, 0.607, 2.2, 1],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [1.66538, .3427048, -.6510217, -.5067066, -.1231787, .08828106],
        "d3": [1, 1, 3, 2, 3, 2],
        "t3": [1.83, 3.3, 1.9, 2.6, 2.9, 3.0],
        "alfa3": [1.108, 1.579, 1.098, 0.672, 3.38, 1.6],
        "beta3": [0.563, 1.724, 0.806, 0.505, 26.4, 8.82],
        "gamma3": [1.246, 1.05, 1, 0.677, 1.302, 1.274],
        "epsilon3": [0.933, 0.786, 0.496, 0.327, 0.523, 0.308]}

    akasaka2014 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R1234ze(Z) of Akasaka "
                    "(2014)",
        "__doi__": {
            "autor": "Akasaka, R., Higashi, Y., Miyara, A., Koyama, S.",
            "title": "A fundamental equation of state for cis-1,3,3,3-"
                     "tetrafluoropropene (R-1234ze(Z))",
            "ref": "Int. J. Refrig. 44 (2014) 168-176",
            "doi": "10.1016/j.ijrefrig.2013.12.018"},

        "R": 8.314472,
        "cp": CP1,
        "ref": "IIR",

        "Tmin": 273, "Tmax": 430, "Pmax": 6000,

        "nr1": [7.7652368, -8.7025756, -0.28352251, 0.14534501, 0.0092092105],
        "d1": [4, 1, 1, 2, 5],
        "t1": [0.685, 0.8494, 1.87, 2, 0.142],

        "nr2": [-0.24997382, 0.096674360, 0.024685924, -0.013255083,
                -0.064231330, 0.36638206, -0.25548847, -0.095592361,
                0.086271444, 0.015997412, -0.013127234, 0.0042293990],
        "d2": [1, 3, 5, 7, 1, 2, 2, 3, 4, 2, 3, 5],
        "t2": [4.2, 0.08, 0, 1.1, 5.5, 6.6, 8.4, 7.2, 7.6, 8.5, 23, 18],
        "c2": [1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3],
        "gamma2": [1]*12}

    eq = (akasaka, akasaka2014)

    _surface = {
        "__doi__": {
            "autor": "Kondou, C., Nagata, R., Nii, N., Koyama, S., Higashi, Y",
            "title": "Surface Tension of low GWP refrigerants R1243zf, R1234ze"
                     "(Z) and R1233zd(E)",
            "ref": "Int. J. Refrigeration 53 (2015) 80-89",
            "doi": "10.1016/j.ijrefrig.2015.01.005"},
        "sigma": [0.05657], "exp": [1.220]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.7093, 2.3374, -2.1124, -3.1074],
        "t": [1.0, 1.5, 2, 4.2]}
    _liquid_Density = {
        "eq": 1,
        "n": [1.3241, 2.3135, -1.2904, 0.67545],
        "t": [0.265, 0.75, 1.3, 1.95]}
    _vapor_Density = {
        "eq": 2,
        "n": [-1.9019, -6.4503, -15.73, -47.277],
        "t": [0.3, 0.96, 2.7, 5.8]}


class Test(TestCase):
    """Test"""

    def test_akasaka(self):
        """Table 6, pag 4689"""
        st = R1234zeZ(T=300, rhom=0)
        self.assertEqual(round(st.P.MPa, 6), 0)
        self.assertEqual(round(st.cvM.kJmolK, 7), 0.0858698)
        self.assertEqual(round(st.cpM.kJmolK, 7), 0.0941842)
        self.assertEqual(round(st.w, 3), 154.887)

        st = R1234zeZ(T=300, rhom=11)
        self.assertEqual(round(st.P.MPa, 5), 10.99498)
        self.assertEqual(round(st.cvM.kJmolK, 6), 0.100637)
        self.assertEqual(round(st.cpM.kJmolK, 6), 0.141580)
        self.assertEqual(round(st.w, 3), 704.695)

        st = R1234zeZ(T=300, rhom=0.05)
        self.assertEqual(round(st.P.MPa, 7), 0.1189126)
        self.assertEqual(round(st.cvM.kJmolK, 7), 0.0883387)
        self.assertEqual(round(st.cpM.kJmolK, 7), 0.0991296)
        self.assertEqual(round(st.w, 3), 149.222)

        st = R1234zeZ(T=400, rhom=8)
        self.assertEqual(round(st.P.MPa, 6), 5.085752)
        self.assertEqual(round(st.cvM.kJmolK, 6), 0.119956)
        self.assertEqual(round(st.cpM.kJmolK, 6), 0.191425)
        self.assertEqual(round(st.w, 3), 302.189)

        st = R1234zeZ(T=400, rhom=1)
        self.assertEqual(round(st.P.MPa, 6), 2.115059)
        self.assertEqual(round(st.cvM.kJmolK, 6), 0.121000)
        self.assertEqual(round(st.cpM.kJmolK, 6), 0.181000)
        self.assertEqual(round(st.w, 3), 121.387)

        st = R1234zeZ(T=424, rhom=4)
        self.assertEqual(round(st.P.MPa, 6), 3.578101)
        self.assertEqual(round(st.cvM.kJmolK, 6), 0.145160)
        self.assertEqual(round(st.cpM.kJmolK, 4), 10.1999)
        self.assertEqual(round(st.w, 4), 82.9160)
