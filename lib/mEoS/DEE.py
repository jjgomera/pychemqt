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
from lib.mEoS import C3


class DEE(MEoS):
    """Multiparameter equation of state for diethyl ether"""
    name = "diethyl ether"
    CASNumber = "60-29-7"
    formula = "C4H10O"
    synonym = ""
    _refPropName = "DEE"
    _coolPropName = "DiethylEther"
    rhoc = unidades.Density(264)
    Tc = unidades.Temperature(466.7)
    Pc = unidades.Pressure(3720.238, "kPa")
    M = 74.1216  # g/mol
    Tt = unidades.Temperature(156.92)
    Tb = unidades.Temperature(307.604)
    f_acent = 0.282
    momentoDipolar = unidades.DipoleMoment(1.151, "Debye")
    id = 162

    Fi1 = {"ao_log": [1, 3.36281],
           "pow": [0, 1, -1, -2, -3],
           "ao_pow": [17.099494, -6.160844, -8.943822, 0.54621, -0.016604]}

    thol = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for diethyl ether of Thol "
                    "et al. (2013)",
        "__doi__": {"autor": "Thol, M., Piazza, L., Span, R.",
                    "title": "A New Functional Form for Equations of State "
                             "for Some Weakly Associating Fluids",
                    "ref": "Int. J. Thermophys., 35(5):783-811, 2014.",
                    "doi": "10.1007/s10765-014-1633-1"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": 270.0, "Tmax": 500.0, "Pmax": 100000.0, "rhomax": 10.6851,

        "nr1": [0.376700499, -0.116630334, -0.73801498, -0.2725701,
                -0.04979231, 0.172267029, 0.0044161891],
        "d1": [1, 1, 1, 2, 3, 3, 5],
        "t1": [-0.75, -0.25, 1.25, 0.75, -1.0, -0.375, 1.25],

        "nr2": [-1.53951612, 1.15606052, -0.0184504019, -0.101800599,
                -0.403598704, 0.00213055571, -0.154741976, 0.0120950552,
                -0.0143106371],
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
        "sigma": [1.2348, -1.1808], "exp": [1.4846, 1.5021]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.3059, 1.1734, 0.7142, -4.3219],
        "t": [1.0, 1.5, 2.2, 3.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.3275, 3.1842, -2.1407, 1.4376],
        "t": [0.12, 0.55, 1.0, 1.4]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.35858, -16.843, 32.476, -33.444, -48.036],
        "t": [0.06, 0.87, 1.3, 1.7, 5.3]}

    trnECS = {"__name__": "Huber (2018)",

              "__doi__": {
                  "autor": "Huber, M.L.",
                  "title": "Models for Viscosity, Thermal Conductivity, and "
                           "Surface Tension of Selected Pure Fluids as "
                           "Implemented in REFPROP v10.0",
                  "ref": "NISTIR 8209",
                  "doi": "10.6028/NIST.IR.8209"},

              "eq": "ecs",
              "ref": C3,
              "visco": "visco1",

              "ek": 370.6, "sigma": 0.53, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 0.9,

              "psi": [1.15039, -0.1535, 0.0330048], "psi_d": [0, 1, 2],
              "fint": [0.00132], "fint_t": [0],
              "chi": [1.16276, -0.0357361], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.196e-9, "gam0": 0.066, "qd": 0.645e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""

    def test_thol(self):
        """Table 9, Pag 26"""
        # Discard the last 4 number, I'm fairly sure is a problem with the
        # significative figures in the equation parameters in paper

        st = DEE(T=280, rho=0.1)
        self.assertEqual(round(st.P.MPa, 9), 0.003134775)
        self.assertEqual(round(st.h.kJkg, 5), -29.23155)
        self.assertEqual(round(st.s.kJkgK, 7), 0.2889437)
        self.assertEqual(round(st.cv.kJkgK, 7), 1.4556700)
        self.assertEqual(round(st.cp.kJkgK, 7), 1.5692554)
        self.assertEqual(round(st.w, 5), 183.65191)

        st = DEE(T=280, rho=750)
        self.assertEqual(round(st.P.MPa, 5), 20.77428)
        self.assertEqual(round(st.h.kJkg, 4), -396.0712)
        self.assertEqual(round(st.s.kJkgK, 6), -1.386035)
        self.assertEqual(round(st.cv.kJkgK, 7), 1.7323676)
        self.assertEqual(round(st.cp.kJkgK, 7), 2.2263812)
        self.assertEqual(round(st.w, 5), 1190.89218)

        st = DEE(T=400, rho=0.1)
        self.assertEqual(round(st.P.MPa, 9), 0.004483768)
        self.assertEqual(round(st.h.kJkg, 4), 182.4662)
        self.assertEqual(round(st.s.kJkgK, 6), 0.873989)
        self.assertEqual(round(st.cv.kJkgK, 6), 1.841364)
        self.assertEqual(round(st.cp.kJkgK, 6), 1.953942)
        self.assertEqual(round(st.w, 4), 218.0490)

        st = DEE(T=400, rho=650)
        self.assertEqual(round(st.P.MPa, 5), 33.84722)
        self.assertEqual(round(st.h.kJkg, 4), -101.4186)
        self.assertEqual(round(st.s.kJkgK, 6), -0.568534)
        self.assertEqual(round(st.cv.kJkgK, 5), 1.94956)
        self.assertEqual(round(st.cp.kJkgK, 6), 2.520106)
        self.assertEqual(round(st.w, 4), 919.5757)

    def test_Huber(self):
        """Table 7, pag 266"""
        st = DEE(T=420, rhom=7.267)
        # self.assertEqual(round(st.mu.muPas, 5), 82.72072)
        self.assertEqual(round(st.mu.muPas, 5), 82.72174)
        self.assertEqual(round(st.k.mWmK, 4), 80.8456)
