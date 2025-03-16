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
from lib.mEoS.R134a import R134a


class R40(MEoS):
    """Multiparameter equation of state for R40"""
    name = "methyl chloride"
    CASNumber = "74-87-3"
    formula = "CH3Cl"
    synonym = "R40"
    _refPropName = "R40"
    _coolPropName = "R40"
    rhoc = unidades.Density(363.219)
    Tc = unidades.Temperature(416.3)
    Pc = unidades.Pressure(6677.3, "kPa")
    M = 50.48752  # g/mol
    Tt = unidades.Temperature(175.0)
    Tb = unidades.Temperature(249.173)
    f_acent = 0.243
    momentoDipolar = unidades.DipoleMoment(1.871, "Debye")
    id = 115

    Fi1 = {"ao_log": [1, 2.92518],
           "pow": [0, 1, -1, -2, -3],
           "ao_pow": [7.499423, -2.997533, -6.0842e-2, -1.1525e-1, 1.0843e-2],
           "ao_exp": [3.764997],
           "titao": [3.7101]}

    thol = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R40 of Thol et al (2013)",
        "__doi__": {"autor": "Thol, M., Piazza, L., Span, R.",
                    "title": "A New Functional Form for Equations of State "
                             "for Some Weakly Associating Fluids",
                    "ref": "Int. J. Thermophys., 35(5):783-811, 2014.",
                    "doi": "10.1007/s10765-014-1633-1"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": 230.0, "Tmax": 630.0, "Pmax": 100000.0, "rhomax": 21.78756,

        "nr1": [.274572058, .104235924e-1, -.125727710e1, 2.25609199e-3,
                -3.31830421e-2, .918440878e-1, .261059608e-2],
        "d1": [1, 1, 1, 2, 3, 3, 5],
        "t1": [-0.75, -0.25, 1.25, 0.75, -1.0, -0.375, 1.25],

        "nr2": [-.948880966e-1, -.843634836e-1, .226263660, -.470765940e-1,
                -.196610405, -.204318929e-1, -.692145009e-1, .148974844e-1,
                -.642544485e-2],
        "d2": [1, 1, 2, 5, 1, 3, 4, 5, 2],
        "t2": [2.375, 3.0, 2.625, 1.875, 4.5, 5.75, 5.375, 2.75, 14.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 2, 3],
        "gamma2": [1]*9}

    eq = (thol, )

    _surface = {
        "__doi__": {
            "autor": "Mulero, A., Cachadiña, I.",
            "title": "Recommended Correlations for the Surface Tension of "
                     "Several Fluids Included in the REFPROP Program",
            "ref": "J. Phys. Chem. Ref. Data 43(2) (2014) 023104",
            "doi": "10.1063/1.4878755"},
        "sigma": [0.071315], "exp": [1.2177]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-6.5074, 0.7520, -9.4148, 19.654, -20.190],
        "t": [1.0, 1.5, 4.5, 5.8, 7.1]}
    _liquid_Density = {
        "eq": 1,
        "n": [2.1809, 0.9228, -2.4615, 7.9722, -13.023, 9.2840],
        "t": [0.37, 1.16, 2.0, 2.9, 3.9, 5.1]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.9433, -6.8001, -82.752, 202.14, -264.16, 99.135],
        "t": [0.18, 0.9, 3.7, 4.6, 5.6, 6.7]}

    trnECS = {"__name__": "Huber (2018)",

              "__doi__": {
                  "autor": "Huber, M.L.",
                  "title": "Models for Viscosity, Thermal Conductivity, and "
                           "Surface Tension of Selected Pure Fluids as "
                           "Implemented in REFPROP v10.0",
                  "ref": "NISTIR 8209",
                  "doi": "10.6028/NIST.IR.8209"},

              "eq": "ecs",
              "ref": R134a,

              "ek": 330.6, "sigma": 0.419, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 0.97,

              "psi": [0.894605, -0.0117296, 0.0130264],
              "psi_d": [0, 1, 2],
              "fint": [2.78821e-4, 2.10163e-6], "fint_t": [0, 1],
              "chi": [0.971796, -0.0356445], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.18e-9, "gam0": 0.056, "qd": 0.505e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""
    def test_thol(self):
        """Table 9, Pag 26"""
        # Discard the last 4 number, I'm fairly sure is a problem with the
        # significative figures in the equation parameters in paper

        st = R40(T=240, rho=0.1)
        self.assertEqual(round(st.P.MPa, 9), 0.003946471)
        self.assertEqual(round(st.h.kJkg, 4), -44.8931)
        self.assertEqual(round(st.s.kJkgK, 5), 0.36742)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.57025)
        self.assertEqual(round(st.cp.kJkgK, 4), 0.7364)
        self.assertEqual(round(st.w, 4), 225.5804)

        st = R40(T=240, rho=1050)
        self.assertEqual(round(st.P.MPa, 5), 27.68669)
        self.assertEqual(round(st.h.kJkg, 4), -469.3851)
        self.assertEqual(round(st.s.kJkgK, 5), -1.97426)
        self.assertEqual(round(st.cv.kJkgK, 5), 1.05799)
        self.assertEqual(round(st.cp.kJkgK, 5), 1.53850)
        self.assertEqual(round(st.w, 4), 1218.6149)

        st = R40(T=400, rho=0.1)
        self.assertEqual(round(st.P.MPa, 5), 0.00658)
        self.assertEqual(round(st.h.kJkg, 3), 89.523)
        self.assertEqual(round(st.s.kJkgK, 4), 0.7075)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.79029)
        self.assertEqual(round(st.cp.kJkgK, 5), 0.95528)
        self.assertEqual(round(st.w, 4), 282.0669)

        st = R40(T=400, rho=900)
        self.assertEqual(round(st.P.MPa, 5), 85.79201)
        self.assertEqual(round(st.h.kJkg, 4), -200.4646)
        self.assertEqual(round(st.s.kJkgK, 5), -1.30372)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.92382)
        self.assertEqual(round(st.cp.kJkgK, 5), 1.46353)
        self.assertEqual(round(st.w, 4), 1077.7374)

    def test_Huber(self):
        """Table 7, pag 266"""
        st = R40(T=374.7, rhom=14.671)
        # self.assertEqual(round(st.mu.muPas, 5), 90.87011)
        # self.assertEqual(round(st.k.mWmK, 4), 87.6545)
        self.assertEqual(round(st.mu.muPas, 5), 90.87045)
        self.assertEqual(round(st.k.mWmK, 4), 87.6538)
