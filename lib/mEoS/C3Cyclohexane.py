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


class C3Cyclohexane(MEoS):
    """Multiparameter equation of state for propylcyclohexane"""
    name = "propylcyclohexane"
    CASNumber = "1678-92-8"
    formula = "C6H11-CH2CH2CH3"
    synonym = ""
    _refPropName = "C2CC6"
    _coolPropName = ""
    rhoc = unidades.Density(260.0527932)
    Tc = unidades.Temperature(630.8)
    Pc = unidades.Pressure(2860.0, "kPa")
    M = 126.23922  # g/mol
    Tt = unidades.Temperature(178.2)
    Tb = unidades.Temperature(429.9)
    f_acent = 0.326
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 184

    CP1 = {"ao": 0.0,
           "an": [9.29427],
           "pow": [0.385871],
           "ao_exp": [1.37051, 106.426, 313.713],
           "exp": [173295, 561.14, 1919.52]}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for propylcyclohexane "
                    "of Lemmon (2007).",
        "__doi__": {"autor": "Lemmon, E.W.",
                    "title": "",
                    "ref": "unpublished equation, 2007",
                    "doi": ""},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 650., "Pmax": 50000.0, "rhomax": 7.03,
        "Pmin": 0.0000007, "rhomin": 7.03,

        "nr1": [1.01911, -2.59762, 0.675152, -0.230891, 0.120966, 0.000309038],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.2, 1.2, 1.8, 1.5, 0.3, 0.9],

        "nr2": [0.526461, -0.0188462, -0.549272, -0.139233, 0.121242],
        "d2": [2, 5, 1, 4, 1],
        "t2": [1.4, 2.2, 3.7, 4.2, 2.4],
        "c2": [1, 1, 2, 2, 1],
        "gamma2": [1]*5}

    eq = lemmon,

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.76296e1, 0.16538e1, -0.28518e1, -0.28205e1, -0.28144e1],
        "t": [1.0, 1.5, 2.7, 4.7, 15.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.39271e-1, 0.38257e2, -0.65743e2, 0.30332e2, 0.17224],
        "t": [0.1, 0.75, 0.87, 1.0, 5.0]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.64572e1, 0.91228e1, -0.25806e2, -0.59044e2, -0.14709e3],
        "t": [0.6, 1.8, 2.2, 6.0, 14.0]}

    thermo0 = {"__name__": "Perkins (2008)",
               "__doi__": {
                   "autor": "Perkins, R.A. Hammerschmidt, U., Huber, M.L.",
                   "title": "Measurement and Correlation of the Thermal "
                            "Conductivity of Methylcyclohexane and "
                            "Propylcyclohexane from 300 to 600 K at Pressures "
                            "to 60 MPa",
                   "ref": "J. Chem. Eng. Data 53(9) (2008) 2120-2127",
                   "doi": "10.1021/je800255r"},

               "eq": 1,

               "Toref": 630.8, "koref": 1,
               "no": [1.07402e-2, -6.09829e-2, 1.38204e-1, -3.81213e-2],
               "to": [0, 1, 2, 3],

               "Tref_res": 630.8, "rhoref_res": 260.05, "kref_res": 1.,
               "nr": [0.116524, -0.102821, -0.113871, 0.126431, 0.0445827,
                      -0.05946, -0.00545736, 0.0098936],
               "tr": [0, -1, 0, -1, 0, -1, 0, -1],
               "dr": [1, 1, 2, 2, 3, 3, 4, 4],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.2415, "R0": 1.01,
               "Xio": 0.15e-9, "gam0": 0.052, "qd": 6.24e-10, "Tcref": 958.725}

    _thermal = thermo0,


class Test(TestCase):

    def test_Perkins(self):
        # The critical enhancement may differ because the correlation used for
        # viscosity in paper isn't implemented in pychemqt, so the values in
        # test could differ

        # Table 5, pag 2125
        st = C3Cyclohexane(T=300, P=1e5)
        self.assertEqual(round(st.rho, 6), 788.480445)
        self.assertEqual(round(st.k, 5), 0.11078)

        st = C3Cyclohexane(T=450, P=1e5)
        self.assertEqual(round(st.rho, 8), 3.52727123)
        self.assertEqual(round(st.k, 5), 0.02438)

        st = C3Cyclohexane(T=450, P=5e7)
        self.assertEqual(round(st.rho, 6), 729.367278)
        self.assertEqual(round(st.k, 5), 0.10543)

        st = C3Cyclohexane(T=600, P=1e5)
        self.assertEqual(round(st.rho, 8), 2.56885356)
        self.assertEqual(round(st.k, 5), 0.04517)

        st = C3Cyclohexane(T=600, P=4.744e6)
        self.assertEqual(round(st.rho, 6), 501.697589)
        # self.assertEqual(round(st.k, 5), 0.07432)

        st = C3Cyclohexane(T=600, P=5e7)
        self.assertEqual(round(st.rho, 6), 644.863706)
        # self.assertEqual(round(st.k, 5), 0.09768)
