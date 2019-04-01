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


class R12(MEoS):
    """Multiparameter equation of state for R12"""
    name = "dichlorodifluoromethane"
    CASNumber = "75-69-4"
    formula = "CCl2F2"
    synonym = "R12"
    _refPropName = "R12"
    _coolPropName = "R12"
    rhoc = unidades.Density(565.)
    Tc = unidades.Temperature(385.12)
    Pc = unidades.Pressure(4136.1, "kPa")
    M = 120.913  # g/mol
    Tt = unidades.Temperature(116.099)
    Tb = unidades.Temperature(243.398)
    f_acent = 0.17948
    momentoDipolar = unidades.DipoleMoment(0.510, "Debye")
    id = 216

    CP1 = {"ao": 4.003638529,
           "an": [], "pow": [],
           "ao_exp": [3.160638395, .3712598774, 3.562277099, 2.121533311],
           "exp": [1.4334342e3, 2.4300498e3, 6.8565952e2, 4.1241579e2]}

    marx = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-12 of Marx (1992)",
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

        "Tmin": Tt, "Tmax": 525.0, "Pmax": 200000.0, "rhomax": 15.13,
        "Pmin": 0.000243, "rhomin": 15.1253,

        "nr1": [2.075343402, -2.962525996, 0.1001589616e-1, 0.1781347612e-1,
                0.2556929157e-1, 0.2352142637e-2, -0.8495553314e-4],
        "d1": [1, 1, 1, 2, 4, 6, 8],
        "t1": [0.5, 1, 2, 2.5, -0.5, 0, 0],

        "nr2": [-0.1535945599e-1, -0.2108816776, -0.1654228806e-1,
                -0.1181316130e-1, -0.4160295830e-4, 0.2784861664e-4,
                0.1618686433e-5, -0.1064614686, 0.9369665207e-3,
                0.2590095447e-1, -0.4347025025e-1, 0.1012308449,
                -0.1100003438, -0.3361012009e-2, 0.3789190008e-3],
        "d2": [1, 1, 5, 7, 12, 12, 14, 1, 9, 1, 1, 3, 3, 5, 9],
        "t2": [-0.5, 1.5, 2.5, -0.5, 0, 0.5, -0.5, 4, 4, 2, 4, 12, 14, 0, 14],
        "c2": [1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 4],
        "gamma2": [1]*15}

    shortSpan = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-12 of Span and "
                    "Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of State for Technical Applications. "
                             "III. Results for Polar Fluids",
                    "ref": "Int. J. Thermophys., 24(1) (2003) 111-162",
                    "doi": "10.1023/A:1022362231796"},

        "R": 8.31451,
        "cp": CP1,
        "ref": "IIR",
        "M": 120.914, "Tc": 385.12, "rhoc": 565/120.914,

        "Tmin": 173.0, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 13.9,
        "Pmin": 1.1633, "rhomin": 13.892,

        "nr1": [0.10557228e1, -0.33312001e1, 0.10197244e1, 0.84155115e-1,
                0.28520742e-3],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.39625057, 0.63995721, -0.21423411e-1, -0.36249173,
                0.1934199e-2, -0.92993833e-1, -0.24876461e-1],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    eq = marx, shortSpan

    _surface = {"sigma": [-0.000124, 0.05662], "exp": [0.4318, 1.263]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.70834e1, 0.43562e1, -0.35249e1, -0.28872e1, -0.89926],
        "t": [1.0, 1.5, 1.67, 4.14, 10.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.32983e2, -0.10997e3, 0.17067e3, -0.13342e3, 0.42525e2],
        "t": [0.57, 0.72, 0.89, 1.07, 1.25]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.153, -6.4734, -17.346, -15.918e2, -32.492, -120.72],
        "t": [0.418, 1.32, 3.3, 6.6, 7.0, 15.0]}


class Test(TestCase):

    def test_shortSpan(self):
        # Table III, Pag 117
        st = R12(T=500, rho=500, eq="shortSpan")
        self.assertEqual(round(st.cp0.kJkgK, 4), 0.7421)
        self.assertEqual(round(st.P.MPa, 3), 11.552)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.1052)

        st2 = R12(T=600, rho=100, eq="shortSpan")
        self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 121.54)
        self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 0.28621)
