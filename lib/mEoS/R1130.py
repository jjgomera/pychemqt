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


class R1130(MEoS):
    """Multiparameter equation of state for trans-1,2-dichloroethene"""
    name = "Trans-1,2-dichlorothene"
    CASNumber = "156-60-5"
    formula = "CHCl=CHCl"
    synonym = "R1130"
    _refPropName = ""
    _coolPropName = ""
    rhoc = unidades.Density(440.1224912)
    Tc = unidades.Temperature(515.69)
    Pc = unidades.Pressure(5255.46, "kPa")
    M = 96.94328  # g/mol
    Tt = unidades.Temperature(223.31)
    Tb = unidades.Temperature(320.367)
    f_acent = 0.194
    momentoDipolar = unidades.DipoleMoment(0, "Debye")
    id = 659

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [-14.91185845478433, 8.646546650620836],
           "ao_exp": [2.697, 6.264, 3.041],
           "titao": [367/Tc, 1246/Tc, 4030/Tc]}

    huber = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-1130 of Huber (2025)",
        "__doi__": {"autor": "Huber, M.L., Kazakov, A.F., Lemmon, E.W.",
                    "title": "Equation of State for the Thermodynamic "
                             "Properties of Trans-1,2-dichloroethene "
                             "[R-1130(E)]",
                    "ref": "Int. J. Thermophys. 46 (2025) 76",
                    "doi": "10.1007/s10765-025-03535-3"},

        "R": 8.314462618,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 625.0, "Pmax": 30000.0, "rhomax": 20,

        "nr1": [0.0256, 0.9, 0.054, 0.153, -0.49787364754, -0.99000024872],
        "d1": [4, 1, 2, 3, 1, 2],
        "t1": [1, 0.188, 0.2496, 0.458, 0.77, 1.0385],

        "nr2": [-0.318, -0.973, -0.098, -0.034, -0.0971],
        "d2": [1, 1, 1, 2, 3],
        "t2": [0.9, 1.21, 4.7, 6.61575, 2.257],
        "c2": [1, 2, 3, 3, 3],
        "gamma2": [2.43, 1.236548, 0.722, 1.483, 0.2255],

        "nr3": [-0.26, -0.135, -0.21, -0.197],
        "d3": [1, 2, 1, 1],
        "t3": [3.2664, 3.59, 1.2, 1.42],
        "alfa3": [5.10245, 4.572, 4.63, 1.785],
        "beta3": [0, 0, 0, 0],
        "gamma3": [1, 1, 1, 1],
        "epsilon3": [-0.31662, 0.065, -0.28568, 0.99]}

    eq = huber,

    # _surface = {"sigma": [0.0556], "exp": [1.24]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.2526, 3.8005, -7.7286, 10.953, -9.7064],
        "t": [1, 1.5, 2.12, 2.84, 3.58]}
    _liquid_Density = {
        "eq": 1,
        "n": [4.577, -6.1201, 7.6048, -4.7251, 1.6728, 1.8053],
        "t": [0.47, 0.83, 1.23, 1.69, 2.27, 9,8]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.063124, -10.541, 23.782, -56.858, 82.754, -76.847, -70.617],
        "t": [0.05, 0.64, 1.13, 1.74, 2.57, 3.25, 8.66]}


class Test(TestCase):
    """Testing"""
    def test_huber(self):
        """Table 4, Pag 6"""
        st = R1130(T=320, rhom=0)
        self.assertEqual(round(st.P.MPa, 5), 0.00000)
        self.assertEqual(round(st.w, 3), 176.452)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 70.151)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 61.836)
        self.assertEqual(round(st.hM.kJkmol, 1), 52896.5)

        st = R1130(T=320, rhom=0.02)
        self.assertEqual(round(st.P.MPa, 5), 0.05229)
        self.assertEqual(round(st.w, 3), 174.174)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 71.732)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 62.647)
        self.assertEqual(round(st.hM.kJkmol, 1), 52732.4)

        st = R1130(T=320, rhom=12.5)
        self.assertEqual(round(st.P.MPa, 5), 3.39671)
        self.assertEqual(round(st.w, 3), 946.434)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 115.586)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 76.665)
        self.assertEqual(round(st.hM.kJkmol, 1), 24883.3)

        st = R1130(T=400, rhom=0.02)
        self.assertEqual(round(st.P.MPa, 5), 0.06584)
        self.assertEqual(round(st.w, 3), 194.273)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 79.550)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 70.843)
        self.assertEqual(round(st.hM.kJkmol, 1), 58759.7)

        st = R1130(T=400, rhom=11)
        self.assertEqual(round(st.P.MPa, 5), 5.40729)
        self.assertEqual(round(st.w, 3), 672.429)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 124.461)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 81.788)
        self.assertEqual(round(st.hM.kJkmol, 1), 34519.1)

        st = R1130(T=520, rhom=9)
        self.assertEqual(round(st.P.MPa, 5), 19.12429)
        self.assertEqual(round(st.w, 3), 454.513)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 136.581)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 88.434)
        self.assertEqual(round(st.hM.kJkmol, 1), 50110.9)
