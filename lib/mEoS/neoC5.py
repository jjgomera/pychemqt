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


class neoC5(MEoS):
    """Multiparameter equation of state for neopentane"""
    name = "neopentane"
    CASNumber = "463-82-1"
    formula = "C(CH3)4"
    synonym = ""
    _refPropName = "NEOPENTN"
    _coolPropName = "Neopentane"
    rhoc = unidades.Density(235.9265106)
    Tc = unidades.Temperature(433.74)
    Pc = unidades.Pressure(3196.0, "kPa")
    M = 72.14878  # g/mol
    Tt = unidades.Temperature(256.6)
    Tb = unidades.Temperature(282.65)
    f_acent = 0.1961
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 9

    Fi1 = {"ao_log": [1, 3.],
           "pow": [0, 1],
           "ao_pow": [0.8702452614, 1.6071746358],
           "ao_exp": [14.422, 12.868, 17.247, 12.663],
           "titao": [710/Tc, 1725/Tc, 3280/Tc, 7787/Tc]}

    f = 72.151/8.3143
    CP1 = {"ao": -0.435375*f,
           "an": [0.96766e-2*f, -0.11533e-4*f, 0.108006e-7*f, -0.44851e-11*f],
           "pow": [1, 2, 3, 4],
           "ao_exp": [], "exp": []}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for neopentane of "
                    "Lemmon and Span (2006).",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 "
                             "Industrial Fluids",
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi":  "10.1021/je050186n"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 550.0, "Pmax": 200000.0, "rhomax": 8.71,

        "nr1": [1.1136, -3.1792, 1.1411, -0.10467, 0.11754, 0.00034058],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.29553, -0.074765, -0.31474, -0.099401, -0.039569, 0.023177],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    polt = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for neopentane of Polt "
                    "(1992)",
        "__doi__": {"autor": "Polt, A., Platzer, B., Maurer, G.",
                    "title": "Parameter der thermischen Zustandsgleichung von "
                             "Bender fuer 14 mehratomige reine Stoffe",
                    "ref": "Chem. Technik 22(1992)6 , 216/224",
                    "doi": ""},

        "R": 8.3143,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": 273.0, "Tmax": 498.0, "Pmax": 20000.0, "rhomax": 8.511,

        "nr1": [-0.146552261671e1, 0.199230626557e1, -0.500821886276,
                0.119809758161e1, -0.363135896710e1, 0.312770556886e1,
                -2.37405105853, 0.473735725047, 0.101500881659, 0.184937708516,
                -0.290527628579e-1, -0.258919377284e-1, 0.748831217999e-1,
                0.216569936506e-1, -0.100375687935, 0.234924630013e-1],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.146552261671e1, -0.199230626557e1, 0.500821886276,
                -0.834410647812, 0.262918341468e1, -0.188136966583e1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [0.968832]*6}

    eq = lemmon, polt

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.70262e1, 0.20090e1, -0.19932e1, -0.28503e1, -0.53760],
        "t": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.56080e1, -0.13549e2, 0.29912e2, -0.28143e2, 0.89021e1],
        "t": [0.45, 0.7, 1.0, 1.25, 1.6]}
    _vapor_Density = {
        "eq": 2,
        "n": [-.25177e1, -.63565e1, -.11985e3, .43740e3, -.10749e4, .74007e3],
        "t": [0.366, 1.14, 4.0, 5.0, 6.0, 6.5]}


class Test(TestCase):
    def test_shortLemmon(self):
        # Table 10, Pag 842
        st = neoC5(T=435, rhom=3)
        self.assertEqual(round(st.P.kPa, 3), 3256.677)
        self.assertEqual(round(st.hM.kJkmol, 3), 34334.720)
        self.assertEqual(round(st.sM.kJkmolK, 3), 92.525)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 184.435)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 5161.767)
        self.assertEqual(round(st.w, 3), 93.352)
