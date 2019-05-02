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


class R41(MEoS):
    """Multiparameter equation of state for R41"""
    name = "fluoromethane"
    CASNumber = "593-53-3"
    formula = "CH3F"
    synonym = "R41"
    _refPropName = "R41"
    _coolPropName = "R41"
    rhoc = unidades.Density(316.506156)
    Tc = unidades.Temperature(317.28)
    Pc = unidades.Pressure(5897.0, "kPa")
    M = 34.03292  # g/mol
    Tt = unidades.Temperature(129.82)
    Tb = unidades.Temperature(194.84)
    f_acent = 0.2004
    momentoDipolar = unidades.DipoleMoment(1.851, "Debye")
    id = 225

    Fi1 = {"ao_log": [1, 3.],
           "pow": [0, 1, -1],
           "ao_pow": [-4.8676441160, 4.2527951258, -0.0268688568],
           "ao_exp": [5.6936, 2.9351], "titao": [1841/Tc, 4232/Tc]}

    CP1 = {"ao": 4,
           "an": [0.00016937], "pow": [1],
           "ao_exp": [5.6936, 2.9351], "exp": [1841, 4232]}

    CP2 = {"ao": 4.59643,
           "an": [-9.48588e-3, 3.96059e-5, -2.85616e-8], "pow": [1, 2, 3],
           "ao_exp": [], "exp": []}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-41 of Lemmon and "
                    "Span (2006) b.",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 "
                             "Industrial Fluids",
                    "ref": "J. Chem. Eng. Data, 51(3) (2006) 785-850",
                    "doi":  "10.1021/je050186n"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 425.0, "Pmax": 70000.0, "rhomax": 29.66,

        "nr1": [1.6264, -2.8337, 0.0010932, 0.037136, 0.00018724],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.52, 1.12, 4, 0.03, 0.63],

        "nr2": [-.22189, .55021, .0461, -.056405, -.17005, -.032409, -.012276],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [3.4, 2.2, 1.5, 0.1, 4.8, 3.5, 15],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    lemmon2 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-41 of Lemmon and "
                    "Span (2006) a.",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 "
                             "Industrial Fluids",
                    "ref": "J. Chem. Eng. Data, 51(3) (2006) 785-850",
                    "doi":  "10.1021/je050186n"},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 60000.0, "rhomax": 29.6,

        "nr1": [0.85316, -2.6366, 0.69129, 0.054681, 0.00012796],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [-0.37093, 0.33920, -0.0017413, -0.095417, -0.078852, -0.030729,
                -0.011497],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    haynes = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for R-41 of Haynes (1996).",
        "__doi__": {"autor": "Haynes, W.M.",
                    "title": "Thermophysical Properties of HCFC Alternatives",
                    "ref": "NIST, Boulder, Colorado, Final Report for ARTI "
                           "MCLR Project Number 660-50800, 1996. Pag. A-82",
                    "doi":  ""},

        "R": 8.314471,
        "Tc": 317.28, "Pc": 5897, "rhoc": 316.51/M,

        "cp": CP2,
        "ref": "NBP",

        "Tmin": 175.0, "Tmax": 500.0, "Pmax": 60000.0, "rhomax": 27.1006,

        "b": [None, -0.326441485138e-1, 0.338620074694e1, -0.831696847103e2,
              0.139589938388e5, -0.156113972752e7, -0.165160386413e-2,
              0.118821153813e1, -0.137311604695e3, 0.176999573025e6,
              0.164945271187e-4, 0.595329911829e-1, -0.341969857376e2,
              -0.168552064750e-2, -0.758216269071e-2, -0.134800586220e2,
              0.311348265418e-2, -0.651499088798e-4, 0.184033192190e-1,
              -0.281459127843e-3, -0.186344956951e6, 0.110422095705e8,
              -0.147526754027e4, 0.261603025982e8, -0.744431617418e1,
              0.782355157170e3, -0.562784094508e-2, -0.843317187588e3,
              -0.600934897964e-4, 0.145050417148e-1, 0.222324172533e-7,
              -0.204419971811e-4, 0.245556593457e-3]}

    eq = lemmon, lemmon2, haynes

    _surface = {"sigma": [0.05049], "exp": [1.242]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.70970e1, 0.17409e1, -0.11668e1, -0.31830e1, 0.93827],
        "t": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.18181e2, -0.62193e2, 0.85171e2, -0.66958e2, 0.28790e2],
        "t": [0.58, 0.8, 1.0, 1.3, 1.5]}
    _vapor_Density = {
        "eq": 2,
        "n": [-26.966, 54.303, -.36361e2, -.17816e2, -.48535e2, -.86727e2],
        "t": [0.59, 0.72, 0.86, 3.2, 7.0, 15.0]}


class Test(TestCase):
    def test_shortLemmon(self):
        # Table 10, Pag 842, eq b
        st = R41(T=319, rhom=9)
        self.assertEqual(round(st.P.kPa, 3), 6129.100)
        self.assertEqual(round(st.hM.kJkmol, 3), 13670.133)
        self.assertEqual(round(st.sM.kJkmolK, 3), 55.886)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 55.438)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 2796.224)
        self.assertEqual(round(st.w, 3), 189.549)

        st = R41(T=273.15, x=0)
        self.assertEqual(round(st.h.kJkg, 3), 200)
        self.assertEqual(round(st.s.kJkgK, 3), 1)
        # st = R41(T=273.15, x=0, eq="lemmon2")
        # self.assertEqual(round(st.h.kJkg, 3), 200)
        # self.assertEqual(round(st.s.kJkgK, 3), 1)

        # Table 10, Pag 842, eq a
        st = R41(T=319, rhom=9, eq="lemmon2")
        self.assertEqual(round(st.P.kPa, 3), 6130.986)
        # The integration parameters are for the other version
        # self.assertEqual(round(st.hM.kJkmol, 3), 13699.382)
        # self.assertEqual(round(st.sM.kJkmolK, 3), 55.982)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 53.388)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 2724.221)
        self.assertEqual(round(st.w, 3), 194.624)
