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
           "ao_exp": [3.160638395, .3712598774, 3.562277099, 2.121533311],
           "exp": [1.4334342e3, 2.4300498e3, 6.8565952e2, 4.1241579e2]}

    marx = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-12 of Marx (1992)",
        "__doi__": {"autor": "Marx, V., Pruss, A., Wagner, W.",
                    "title": "Neue Zustandsgleichungen fuer R 12, R 22, R 11 "
                             "und R 113. Beschreibung des thermodynamishchen "
                             "Zustandsverhaltens bei Temperaturen bis 525 K "
                             "und Druecken bis 200 MPa",
                    "ref": "Düsseldorf: VDI Verlag, Series 19 "
                           "(Waermetechnik/Kaeltetechnik), No. 57, 1992.",
                    "doi": ""},

        "R": 8.314471, "rhoc": 4.672781,
        "cp": CP1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 525.0, "Pmax": 200000.0, "rhomax": 15.13,

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
    _PR = [-0.1824, -16.2971]

    _surface = {"sigma": [-0.000124, 0.05662], "exp": [0.4318, 1.263]}

    _melting = {
        "eq": 2,
        "__doi__": {
            "autor": "Reeves, L.E., Scott, G.J., Babb, S.E. Jr.",
            "title": "Melting Curves of Pressure-Transmitting fluids",
            "ref": "Fluid Phase Equilib., 222-223 (2004) 107-118",
            "doi": "10.1063/1.1725068"},

        "Tmin": Tt, "Tmax": 2000.0,
        "Tref": Tt, "Pref": 0.,
        "a2": [3288e5], "exp2": [2.231]}

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
        "n": [-3.153, -6.4734, -17.346, -15.918, -32.492, -120.72],
        "t": [0.418, 1.32, 3.3, 6.6, 7.0, 15.0]}

    thermo0 = {"__name__": "Krauss (1989)",
               "__doi__": {
                   "autor": "Krauss, R., Stephan, K.",
                   "title": "Thermal Conductivity of Refrigerants in a Wide "
                            "Range of Temperature and Pressure",
                   "ref": "J. Phys. Chem. Ref. Data 18(1) (1989) 43-76",
                   "doi": "10.1063/1.555842"},

               "eq": 1,

               # Typo in paper
               # Temperature reducing value is critical value
               # Thermal conductivity reducing value using the exact value from
               # dimensional analysis with molecular weight in kg/mol
               # R**(5/6)*Pc**(2/3)/Tc**(1/6)/(M/1000)**0.5/Na**(1/3)
               "Toref": 385.15, "koref": 1.898767e-3,
               "no": [-2.8497878, 10.338094],
               "to": [0, 1],

               "rhoref_res": 4.62*120.914, "kref_res": 1.898767e-3,
               "nr": [.82454302, 14.887641, -13.347132, 5.5118177, -.67277383],
               "tr": [0, 0, 0, 0, 0],
               "dr": [1, 2, 3, 4, 5]}

    _thermal = (thermo0, )


class Test(TestCase):
    """Testing"""
    def test_shortSpan(self):
        """Table III, Pag 117"""
        st = R12(T=500, rho=500, eq="shortSpan")
        self.assertEqual(round(st.cp0.kJkgK, 4), 0.7421)
        self.assertEqual(round(st.P.MPa, 3), 11.552)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.1052)

        st2 = R12(T=600, rho=100, eq="shortSpan")
        self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 121.54)
        self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 0.28621)

    def test_krauss(self):
        """Selected point from Table C1 and C2, pag 67"""
        # The values differ because the paper use and old EoS don't
        # implemented in pychemqt
        self.assertEqual(round(R12(T=200, P=1e5, eq=1).k.mWmK, 2), 104.99)
        self.assertEqual(round(R12(T=600, P=5e6, eq=1).k.mWmK, 2), 26.80)
        self.assertEqual(round(R12(T=400, P=1e7, eq=1).k.mWmK, 2), 44.65)
        self.assertEqual(round(R12(T=240, P=6e7, eq=1).k.mWmK, 2), 108.09)

        # Saturation point, Table C2
        st = R12(T=200, x=0.5, eq=1)
        self.assertEqual(round(st.P.bar, 6), 0.099521)
        self.assertEqual(round(st.Gas.k.mWmK, 2), 4.78)
        self.assertEqual(round(st.Liquido.k.mWmK, 2), 104.96)

        st = R12(T=380, x=0.5, eq=1)
        self.assertEqual(round(st.P.bar, 3), 37.739)
        self.assertEqual(round(st.Gas.k.mWmK, 2), 20.57)
        self.assertEqual(round(st.Liquido.k.mWmK, 2), 36.77)
