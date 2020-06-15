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


class SO2(MEoS):
    """Multiparameter equation of state for sulfur dioxide"""
    name = "sulfur dioxide"
    CASNumber = "7446-09-5"
    formula = "SO2"
    synonym = "R-764"
    _refPropName = "SO2"
    _coolPropName = "SulfurDioxide"
    rhoc = unidades.Density(525.002841)
    Tc = unidades.Temperature(430.64)
    Pc = unidades.Pressure(7886.6, "kPa")
    M = 64.0638  # g/mol
    Tt = unidades.Temperature(197.7)
    Tb = unidades.Temperature(263.13)
    f_acent = 0.256
    momentoDipolar = unidades.DipoleMoment(1.6, "Debye")
    id = 51

    Fi1 = {"ao_log": [1, 3.],
           "pow": [0, 1, -1],
           "ao_pow": [-4.5414235721, 4.4732289572, -0.0159272204],
           "ao_exp": [1.0875, 1.916],
           "titao": [783/Tc, 1864/Tc]}

    Fi2 = {"ao_log": [1, 3.],
           "pow": [0, 1, -1],
           "ao_pow": [-4.5328346436, 4.4777967379, -0.01560057996],
           "ao_exp": [1.062, 1.9401],
           "titao": [775/Tc, 1851/Tc]}

    CP2 = {"ao": 0.4021066/8.3143*64.066,
           "an": [0.87348570e-3/8.3143*64.066, -0.45968820e-6/8.3143*64.066,
                  -0.13328400e-11/8.3143*64.066, 0.23785000e-13/8.3143*64.066],
           "pow": [1, 2, 3, 4]}

    gao = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for sulfur dioxide of Gao "
                    "(2016).",
        "__doi__": {
            "autor": "Gao, K., Wu, J., Zhang, P., Lemmon, E.W.",
            "title": "A Helmholtz Energy Equation of State for Sulfur Dioxide",
            "ref": "J. Chem. Eng. Data, 61(6) (2016) 2859-2872",
            "doi":  "10.1021/acs.jced.6b00195"},

        "R": 8.3144621,
        "cp": Fi1,
        "ref": "NBP",

        "rhoc": 8.078,
        "Tmin": Tt, "Tmax": 525.0, "Pmax": 35000.0, "rhomax": 25.40,

        "nr1": [0.01744413, 1.814878, -2.246338, -0.4602906, 0.1097049],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.45, 0.9994, 1, 0.45],

        "nr2": [-0.9485769, -0.8751541, 0.4228777, -0.4174962, -0.002903451],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.907, 2.992, 0.87, 3.302, 1.002],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [1.64041, -0.4103535, -0.08316597, -0.2728126, -0.1075782,
                -0.4348434],
        "d3": [1, 1, 3, 2, 2, 1],
        "t3": [1.15, 0.997, 1.36, 2.086, 0.855, 0.785],
        "alfa3": [1.061, 0.945, 1.741, 1.139, 1.644, 0.647],
        "beta3": [0.967, 2.538, 2.758, 1.062, 1.039, 0.41],
        "gamma3": [1.276, 0.738, 0.71, 0.997, 1.35, 0.919],
        "epsilon3": [0.832, 0.69, 0.35, 0.961, 0.981, 0.333]}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for sulfur dioxide of "
                    "Lemmon and Span (2006).",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 "
                             "Industrial Fluids",
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi":  "10.1021/je050186n"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 525.0, "Pmax": 35000.0, "rhomax": 25.30,

        "nr1": [0.93061, -1.9528, -0.17467, 0.061524, 0.00017711],
        "d1": [1, 1, 1, 3, 7, ],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.21615, 0.51353, 0.010419, -0.25286, -0.054720, -0.059856,
                -0.016523],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    polt = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for sulfur dioxide of Polt "
                    "(1987).",
        "__doi__": {"autor": "Polt, A., Platzer, B., Maurer, G.",
                    "title": "Parameter der thermischen Zustandsgleichung von "
                             "Bender fuer 14 mehratomige reine Stoffe",
                    "ref": "Chem. Technik 22(1992)6 , 216/224",
                    "doi": ""},

        "R": 8.3143,
        "cp": CP2,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 523.0, "Pmax": 32000.0, "rhomax": 22.91,

        "nr1": [0.789407019882, -0.170449580056e1, 0.115984637964e1,
                -0.576307837294, 0.249237283833e1, -0.518115678632e1,
                0.320766081899e1, -0.123636065893e1, 0.144419600938e-1,
                -.15380705504, 0.386324300525, 0.292550313202, -0.372445361392,
                -0.636924333910e-1, 0.986166451596e-1, -0.216993783055e-2],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [-0.789407019882, 0.170449580056e1, -0.115984637964e1,
                -0.480876182378, 0.164910076886e1, -0.133861069604e1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [1]*6}

    eq = gao, lemmon, polt
    _PR = [-0.068, -18.12]

    _surface = {"sigma": [0.0803, 0.0139, -0.0114],
                "exp": [0.928, 1.57, 0.364]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.303, 1.9794, -2.078, -3.5446, 0.51776],
        "t": [1, 1.5, 2.2, 4.7, 6]}
    _liquid_Density = {
        "eq": 1,
        "n": [7.2296, -16.928, 29.832, -27.901, 11.085],
        "t": [0.54, 0.88, 1.23, 1.6, 2]}
    _vapor_Density = {
        "eq": 2,
        "n": [-7.487, 10.118, -13.608, -25.408, -42.04, -38.668],
        "t": [0.545,  0.85, 1.2, 3.7, 7.5, 10]}


class Test(TestCase):
    def test_gao(self):
        # Table 6, Pag M
        st = SO2(T=250, rhom=23.6)
        self.assertEqual(round(st.P.MPa, 7), 12.2955804)
        self.assertEqual(round(st.cvM.kJkmolK, 4), 53.1514)
        self.assertEqual(round(st.cpM.kJkmolK, 4), 86.0982)
        self.assertEqual(round(st.w, 2), 1130.24)

        st = SO2(T=400, rhom=16)
        self.assertEqual(round(st.P.MPa, 7), 8.0793790)
        self.assertEqual(round(st.cvM.kJkmolK, 4), 51.8705)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 117.691)
        self.assertEqual(round(st.w, 3), 449.618)

        st = SO2(T=431, rhom=8.078)
        self.assertEqual(round(st.P.MPa, 7), 7.9343772)
        self.assertEqual(round(st.cvM.kJkmolK, 4), 64.5073)
        self.assertEqual(round(st.cpM.kJkmolK, 1), 19127.4)
        self.assertEqual(round(st.w, 3), 168.147)

        st = SO2(T=250, rhom=0)
        self.assertEqual(round(st.P.MPa, 7), 0)
        self.assertEqual(round(st.cvM.kJkmolK, 4), 29.8406)
        self.assertEqual(round(st.cpM.kJkmolK, 4), 38.1551)
        self.assertEqual(round(st.w, 3), 203.682)

        st = SO2(T=420, rhom=1)
        self.assertEqual(round(st.P.MPa, 7), 2.9365903)
        self.assertEqual(round(st.cvM.kJkmolK, 4), 40.8928)
        self.assertEqual(round(st.cpM.kJkmolK, 4), 59.5297)
        self.assertEqual(round(st.w, 3), 234.103)

        st = SO2(T=450, rhom=11)
        self.assertEqual(round(st.P.MPa, 7), 12.1084452)
        self.assertEqual(round(st.cvM.kJkmolK, 4), 54.7870)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 222.083)
        self.assertEqual(round(st.w, 3), 250.095)

    def test_shortLemmon(self):
        # Table 10, Pag 842
        st = SO2(T=432, rhom=8, eq="lemmon")
        self.assertEqual(round(st.P.kPa, 3), 8052.256)
        self.assertEqual(round(st.hM.kJkmol, 3), 20821.200)
        self.assertEqual(round(st.sM.kJkmolK, 3), 56.819)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 61.478)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 4877.456)
        self.assertEqual(round(st.w, 3), 171.538)
