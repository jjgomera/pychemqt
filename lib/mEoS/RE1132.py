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


class RE1132(MEoS):
    """Multiparameter equation of state for RE1132"""
    name = "trans-1,2-difluroethene"
    CASNumber = "1630-78-0"
    formula = "CHF=CHF"
    synonym = "RE1132"
    _refPropName = ""
    _coolPropName = ""
    rhoc = unidades.Density(434.989755)
    Tc = unidades.Temperature(348.82)
    Pc = unidades.Pressure(5173.7, "kPa")
    M = 64.035  # g/mol
    Tt = unidades.Temperature(184.9)
    Tb = unidades.Temperature(220.51)
    f_acent = 0.245
    momentoDipolar = unidades.DipoleMoment(0, "Debye")

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [-8.5962232475660088, 46.4405339552517979],
           "ao_exp": [4.896, 7.111],
           "titao": [730/Tc, 2259/Tc]}

    akasaka = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-1132 of Marx (1992)",
        "__doi__": {"autor": "Akasaka, R., Lemmon, E.W.",
                    "title": "A Helmholtz Energy Equations of State for "
                             "Calculations of Thermodynamic Properties of "
                             "trans-1,2-Difluoroethene [R-1132(E)]",
                    "ref": "Int. J. Thermophys. 45(12) (2024) 174",
                    "doi": "10.1007/s10765-024-03447-8"},

        "R": 8.314462618,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": 240, "Tmax": 400.0, "Pmax": 6500.0, "rhomax": 49.10,

        "nr1": [0.033144696, 1.8996205, -2.190517288130, -0.8165909562367,
                0.22857476],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.32, 0.71, 1.2, 0.7],

        "nr2": [-1.650708, -1.4409077, 0.60536471, -0.66134087, -0.016787895],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.09, 2.24, 1.17, 2.2, 0.823],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [0.31476822, -0.31991607, 2.2172762, -0.42786753, -0.41335541,
                0.13023087, -0.30592804],
        "d3": [1, 1, 1, 1, 1, 1, 1],
        "t3": [3.14, 4, 1.25, 1.476, 1.5, 0.73, 1.535],
        "alfa3": [34.4, 34.28, 1.303, 1.855, 1.29, 1.17, 2.165],
        "beta3": [1017.0, 1009.0, 1.011, 1.05, 1.33, 1.61, 0.8],
        "gamma3": [1.0617, 1.062, 1.327, 1.05, 1.157, 1.183, 1.2],
        "epsilon3": [0.9637, 0.9636, 0.8535, 1.15, 1.302, 1.54, 0.786]}

    eq = (akasaka, )

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.2916, 2.4572, -2.2495, -5.0401],
        "t": [1.0, 1.5, 1.98, 4.79]}
    _liquid_Density = {
        "eq": 1,
        "n": [1.6128, 2.7201, -8.5516, 7.3398],
        "t": [0.309, 1.08, 1.86, 2.15]}
    _vapor_Density = {
        "eq": 2,
        "n": [-1.9897, -6.0695, -16.288, -59.359],
        "t": [0.312, 0.995, 2.82, 6.32]}


class Test(TestCase):
    """Testing"""
    def test_akasaka(self):
        """Table 6, Pag 17"""
        st = RE1132(T=270, rhom=0)
        self.assertEqual(round(st.P.MPa, 6), 0)
        self.assertEqual(round(st.cvM.kJkmolK, 4), 48.7929)
        self.assertEqual(round(st.cpM.kJkmolK, 4), 57.1074)
        self.assertEqual(round(st.w, 3), 202.562)

        st = RE1132(T=270, rhom=16)
        self.assertEqual(round(st.P.MPa, 6), 4.271529)
        self.assertEqual(round(st.cvM.kJkmolK, 4), 68.2383)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 108.663)
        self.assertEqual(round(st.w, 3), 687.157)

        st = RE1132(T=270, rhom=0.2)
        self.assertEqual(round(st.P.MPa, 7), 0.4171457)
        self.assertEqual(round(st.cvM.kJkmolK, 4), 51.0279)
        self.assertEqual(round(st.cpM.kJkmolK, 4), 62.7398)
        self.assertEqual(round(st.w, 3), 192.403)

        st = RE1132(T=330, rhom=12)
        self.assertEqual(round(st.P.MPa, 6), 3.845082)
        self.assertEqual(round(st.cvM.kJkmolK, 4), 70.9361)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 165.548)
        self.assertEqual(round(st.w, 3), 314.193)

        st = RE1132(T=330, rhom=2)
        self.assertEqual(round(st.P.MPa, 6), 3.369736)
        self.assertEqual(round(st.cvM.kJkmolK, 4), 66.8528)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 137.744)
        self.assertEqual(round(st.w, 3), 159.782)

        st = RE1132(T=349, rhom=7)
        self.assertEqual(round(st.P.MPa, 6), 5.193981)
        self.assertEqual(round(st.cvM.kJkmolK, 4), 91.0108)
        self.assertEqual(round(st.cpM.kJkmolK, 1), 29733.5)
        self.assertEqual(round(st.w, 3), 124.897)
