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


class R113(MEoS):
    """Multiparameter equation of state for R113"""
    name = "1,1,2-trichloro-1,2,2-trifluoroethane"
    CASNumber = "76-13-1"
    formula = "CCl2FCClF2"
    synonym = "R113"
    _refPropName = "R113"
    _coolPropName = "R113"
    rhoc = unidades.Density(560.)
    Tc = unidades.Temperature(487.21)
    Pc = unidades.Pressure(3392.2, "kPa")
    M = 187.375  # g/mol
    Tt = unidades.Temperature(236.93)
    Tb = unidades.Temperature(320.735)
    f_acent = 0.25253
    momentoDipolar = unidades.DipoleMoment(0.803, "Debye")
    id = 232

    CP1 = {"ao": 3.9999966,
           "an": [], "pow": [],
           "ao_exp": [12.4464495, 2.72181845, 0.692712415, 3.32248298],
           "exp": [5.1143280e2, 1.60676324e3, 4.20292102e3, 1.60618738e3]}

    marx = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-113 of Marx (1992)",
        "__doi__": {"autor": "Marx, V., Pruss, A., and Wagner, W.",
                    "title": "Neue Zustandsgleichungen fuer R 12, R 22, R 11 "
                             "und R 113. Beschreibung des thermodynamishchen "
                             "Zustandsverhaltens bei Temperaturen bis 525 K "
                             "und Druecken bis 200 MPa",
                    "ref": "Düsseldorf: VDI Verlag, Series 19 "
                           "(Waermetechnik/Kaeltetechnik), No. 57, 1992.",
                    "doi": ""},

        "R": 8.314471,
        "cp": CP1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 525.0, "Pmax": 200000.0, "rhomax": 9.10,

        "nr1": [0.8432092286, -0.2019185967e1, 0.2920612996, 0.5323107661e-1,
                0.3214971931e-2, 0.4667858574e-4, -0.1227522799e-5],
        "d1": [1, 1, 2, 3, 4, 8, 8],
        "t1": [0.5, 1.5, 1.5, -0.5, 2, 0, 3],

        "nr2": [0.8167288718, -0.1340790803e1, 0.4065752705, -0.1534754634,
                -0.2414435149e-1, -0.2113056197e-1, -0.3565436205e-1,
                0.1364654968e-2, -0.1251838755e-1, -0.1385761351e-2,
                0.7206335486e-3],
        "d2": [3, 3, 3, 5, 1, 2, 2, 9, 3, 7, 8],
        "t2": [-0.5, 0, 2, 1.5, 6, 2, 10, 6, 18, 15, 33],
        "c2": [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4],
        "gamma2": [1]*11}

    shortSpan = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-113 of Span and "
                    "Wagner (2003).",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of State for Technical Applications. "
                             "III. Results for Polar Fluids",
                    "ref": "Int. J. Thermophys., 24(1) (2003) 111-162",
                    "doi": "10.1023/A:1022362231796"},

        "R": 8.31451,
        "cp": CP1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 9.09,

        "nr1": [0.10519071e1, -0.28724742e1, 0.41983153, 0.87107788e-1,
                0.24105194e-3],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.70738262, 0.93513411, -0.96713512e-2, -0.52595315,
                0.22691984e-1, -0.14556325, -0.2741995e-1],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    eq = marx, shortSpan

    _surface = {"sigma": [0.0556], "exp": [1.24]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.73838e1, 0.32594e1, -0.27761e1, -0.37758e1, -0.19921],
        "t": [1.0, 1.5, 1.8, 4.3, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.15785e1, 0.12404e1, -0.66933, 0.49775e1, -0.55253e1],
        "t": [0.3, 0.7, 2.0, 4.0, 5.0]}
    _vapor_Density = {
        "eq": 2,
        "n": [-2.6225, -6.0753, -15.768, -42.361, -7.9071, -319.66],
        "t": [0.379, 1.13, 2.9, 6.0, 7.0, 15.0]}


class Test(TestCase):

    def test_shortSpan(self):
        # Table III, Pag 117
        st = R113(T=500, rho=500, eq="shortSpan")
        self.assertEqual(round(st.cp0.kJkgK, 4), 0.8055)
        self.assertEqual(round(st.P.MPa, 3), 3.962)
        self.assertEqual(round(st.cp.kJkgK, 4), 4.0344)

        st2 = R113(T=600, rho=100, eq="shortSpan")
        self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 131.00)
        self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 0.26005)
