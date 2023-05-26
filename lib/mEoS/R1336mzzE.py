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
from lib.mEoS import R134a


class R1336mzzE(MEoS):
    """Multiparameter equation of state for R1336mzzE"""
    name = "trans-1,1,1,4,4,4-hexafluorobutene"
    CASNumber = "66711-86-2"
    formula = "trans-CF3CH=CHCF3"
    synonym = "R-1336mzz(E)"
    rhoc = unidades.Density(513.3096339)
    Tc = unidades.Temperature(403.53)
    Pc = unidades.Pressure(2779, "kPa")
    M = 164.0491  # g/mol
    Tt = unidades.Temperature(200.15)
    Tb = unidades.Temperature(281.02)
    f_acent = 0.413
    momentoDipolar = unidades.DipoleMoment(0, "Debye")

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [-17.5838598852, 11.5244170129],
           "ao_exp": [15.891, 14.143],
           "titao": [550/Tc, 3474/Tc]}

    akasaka = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R1336mzz(E) of Akasaka"
                    "(2023)",
        "__doi__": {"autor": "McLinden, M.O., Akasaka, R.",
                    "title": "A Helmholtz Energy Equation of State for "
                             "trans-1,1,1,4,4,4-Hexafluoro-2-butene "
                             "[R-1336mzz(E)] and an Auxiliary Extended "
                             "Corresponding States Model for the Transport "
                             "Properties",
                    "ref": "Int. J. Thermophys 44(4) (2023) 50",
                    "doi": "10.1007/s10765-022-03143-5"},

        "R": 8.314462618,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": 230, "Tmax": 460, "Pmax": 36000.0,

        "nr1": [0.08297005, 1.1213658, -2.0279038680756, -0.3486706213113,
                0.13227952],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.117, 1, 1, 0.413],


        "nr2": [-1.3751844, -1.3939029, 0.11190839, -0.90635088, -0.050014594],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.44, 2.51, 0.535, 1.927, 1.186],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [1.7953748, 2.0579712, -1.0433081, -1.2459808, -0.49712249],
        "d3": [1, 1, 3, 2, 2],
        "t3": [0.876, 1.23, 0.875, 0.5181, 0.86],
        "alfa3": [1.234, 1.34, 1.13, 1.15, 1.3],
        "beta3": [1.237, 0.384, 1.24, 0.46, 1.326],
        "gamma3": [1.21, 2, 1.216, 2.27, 1.26],
        "epsilon3": [1.23, 0.751, 0.522, 0.32, 1.22]}

    eq = (akasaka, )

    _surface = {
        "__doi__": {
            "autor": "Iwasaki, S., Kondou, C., Higashi, Y.",
            "title": "Correlation Assessment and Temperature Dependence Check "
                     "of Surface Tension and Parachor for New Low-GWP Pure "
                     "Refrigerants",
            "ref": "Trans. JSRAE 37(1) (2020) 73-80",
            "doi": "10.11322/tjsrae.19-36TG_EM_OA"},
        "sigma": [0.05371], "exp": [1.27]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-8.2611, 3.9076, -3.9675, -6.1565],
        "t": [1.0, 1.5, 1.88, 4.47]}
    _liquid_Density = {
        "eq": 1,
        "n": [2.1958, 6.4706, -13.191, 7.7906],
        "t": [0.347, 1.42, 1.82, 2.232]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.1511, -9.2173, -32.976, -86.791],
        "t": [0.387, 1.34, 3.77, 7.85]}


class Test(TestCase):
    """Test class"""

    def test_akasaka(self):
        """Table 7, pag 27 in manuscript"""
        st = R1336mzzE(T=300, rhom=0)
        self.assertEqual(round(st.P.MPa, 6), 0)
        self.assertEqual(round(st.cvM.JmolK, 3), 125.686)
        self.assertEqual(round(st.cpM.JmolK, 3), 134.001)
        self.assertEqual(round(st.w, 3), 127.321)

        st = R1336mzzE(T=300, rhom=8.5)
        self.assertEqual(round(st.P.MPa, 5), 22.73753)
        self.assertEqual(round(st.cvM.JmolK, 3), 146.214)
        self.assertEqual(round(st.cpM.JmolK, 3), 193.851)
        self.assertEqual(round(st.w, 3), 725.684)

        st = R1336mzzE(T=300, rhom=0.05)
        self.assertEqual(round(st.P.MPa, 6), 0.119033)
        self.assertEqual(round(st.cvM.JmolK, 3), 127.121)
        self.assertEqual(round(st.cpM.JmolK, 3), 137.553)
        self.assertEqual(round(st.w, 3), 122.301)

        st = R1336mzzE(T=380, rhom=6)
        self.assertEqual(round(st.P.MPa, 6), 2.279713)
        self.assertEqual(round(st.cvM.JmolK, 3), 154.022)
        self.assertEqual(round(st.cpM.JmolK, 3), 254.161)
        self.assertEqual(round(st.w, 3), 209.403)

        st = R1336mzzE(T=380, rhom=0.5)
        self.assertEqual(round(st.P.MPa, 6), 1.212781)
        self.assertEqual(round(st.cvM.JmolK, 3), 144.095)
        self.assertEqual(round(st.cpM.JmolK, 3), 172.608)
        self.assertEqual(round(st.w, 3), 113.306)

        st = R1336mzzE(T=405, rhom=3)
        self.assertEqual(round(st.P.MPa, 6), 2.859345)
        self.assertEqual(round(st.cvM.JmolK, 3), 165.611)
        self.assertEqual(round(st.cpM.JmolK, 2), 4690.00)
        self.assertEqual(round(st.w, 4), 71.6230)
