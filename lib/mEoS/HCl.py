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

from lib.meos import MEoS
from lib import unidades


class HCl(MEoS):
    """Multiparameter equation of state for hydrogen chloride"""
    name = "hydrogen chloride"
    CASNumber = "7647-01-0 "
    formula = "HCl"
    synonym = ""
    rhoc = unidades.Density(410.96998439164605)
    Tc = unidades.Temperature(324.55)
    Pc = unidades.Pressure(8263.00, "kPa")
    M = 36.460939  # g/mol
    Tt = unidades.Temperature(131.4)
    Tb = unidades.Temperature(188.199)
    f_acent = 0.12875
    momentoDipolar = unidades.DipoleMoment(1.079, "Debye")
    id = 104

    Fi1 = {"ao_log": [1, 2.5],
           "pow": [0, 1, -1, -2],
           "ao_pow": [7.913048, -3.217575, -4.149937e-3, 8.019202e-4],
           "ao_exp": [1.054392],
           "titao": [1.241138e1]}

    thol = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hydrogen chloride of Thol"
                    " et al. (2013)",
        "__doi__": {"autor": "Thol, M., Piazza, L., and Span, R.",
                    "title": "A New Functional Form for Equations of State "
                             "for Some Weakly Associating Fluids",
                    "ref": "Int. J. Thermophys., 35(5):783-811, 2014.",
                    "doi": "10.1007/s10765-014-1633-1"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": 155.0, "Tmax": 330.0, "Pmax": 20000.0, "rhomax": 33.8145,
        "Pmin": 0.7, "rhomin": 33.8145,

        "nr1": [-.40937325, 0.943994574, -1.78830477, 0.128619044,
                4.39018427e-3, 0.0130480908, 1.69387782e-3],
        "d1": [1, 1, 1, 2, 3, 3, 5],
        "t1": [-0.75, -0.25, 1.25, 0.75, -1.0, -0.375, 1.25],

        "nr2": [0.751559060, -0.800007427, 0.430935939, 4.54319457e-3,
                -1.52172259e-1, -4.36174059e-2, -9.70625964e-3, 1.01144098e-2,
                3.76991644e-3],
        "d2": [1, 1, 2, 5, 1, 3, 4, 5, 2],
        "t2": [2.375, 3.0, 2.625, 1.875, 4.5, 5.75, 5.375, 2.75, 14.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 2, 3],
        "gamma2": [1]*9}

    eq = thol,

    _vapor_Pressure = {
        "eq": 6,
        "ao": [-0.01065138, -6.15979914, 1.55860976, -8.42734117],
        "exp": [1.0, 2.0, 6.0, 11.0]}
    _liquid_Density = {
        "eq": 2,
        "ao": [1.89232034, 0.83621066, -0.22094602, 4.70971253, -5.34396174],
        "exp": [1.0, 2.0, 4.0, 11.0, 13.0]}
    _vapor_Density = {
        "eq": 4,
        "ao": [-2.95523223, -8.10448179, -14.78392979, -87.13352586],
        "exp": [1.29, 4.2, 11.1, 24.0]}


class Test(TestCase):
    def test_thol(self):
        # Discard the last 4 number, I'm fairly sure is a problem with the
        # significative figures in the equation parameters in paper

        # Table 9, Pag 26
        st = HCl(T=170, rho=0.01)
        self.assertEqual(round(st.P.MPa, 9), 0.000387586)
        self.assertEqual(round(st.h.kJkg, 5), -102.41433)
        self.assertEqual(round(st.s.kJkgK, 5), 0.82032)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.57116)
        self.assertEqual(round(st.cp.kJkgK, 5), 0.79948)
        self.assertEqual(round(st.w, 5), 232.89853)

        st = HCl(T=170, rho=1230)
        self.assertEqual(round(st.P.MPa, 4), 1.2292)
        self.assertEqual(round(st.h.kJkg, 5), -561.33940)
        self.assertEqual(round(st.s.kJkgK, 5), -2.88832)
        self.assertEqual(round(st.cv.kJkgK, 5), 1.14972)
        self.assertEqual(round(st.cp.kJkgK, 5), 1.55394)
        self.assertEqual(round(st.w, 3), 999.439)

        st = HCl(T=280, rho=0.1)
        self.assertEqual(round(st.P.MPa, 5), 0.00638)
        self.assertEqual(round(st.h.kJkg, 4), -14.6144)
        self.assertEqual(round(st.s.kJkgK, 5), 0.58003)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.57147)
        self.assertEqual(round(st.cp.kJkgK, 5), 0.80008)
        self.assertEqual(round(st.w, 5), 298.83636)

        st = HCl(T=280, rho=900)
        self.assertEqual(round(st.P.MPa, 4), 3.4219)
        self.assertEqual(round(st.h.kJkg, 5), -371.03600)
        self.assertEqual(round(st.s.kJkgK, 5), -2.04410)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.96150)
        self.assertEqual(round(st.cp.kJkgK, 5), 2.15083)
        self.assertEqual(round(st.w, 4), 577.7828)
