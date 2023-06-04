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
from lib.mEoS import R134a


class C4F10(MEoS):
    """Multiparameter equation of state for perfluorobutane"""
    name = "perfluorobutane"
    CASNumber = "355-25-9"
    formula = "C4F10"
    synonym = ""
    _refPropName = "C4F10"
    _coolPropName = ""
    rhoc = unidades.Density(627.677199)
    Tc = unidades.Temperature(386.326)
    Pc = unidades.Pressure(2322.4, "kPa")
    M = 238.027  # g/mol
    Tt = unidades.Temperature(144.0)
    Tb = unidades.Temperature(271.123)
    f_acent = 0.372
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 693

    CP1 = {"ao": 14,
           "ao_exp": [2.164, 15.64], "exp": [368, 810]}

    gao = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for perfluerohexane of Gao "
                    "(2021).",
        "__doi__": {
            "autor": "Gao, K., Köster, A., Thol, M., Wu, J., Lemmon, E.W.",
            "title": "Equations of State for the Thermodynamic Properties of "
                     "n-Perfluorobutane, n-Perfluoropentane, and "
                     "n-Perfluerohexane",
            "ref": "Ind. Eng. Chem. Res. 60(47) (2021) 17207-17227",
            "doi": "10.1021/acs.iecr.1c02969"},
        # Find published paper

        "R": 8.314462618,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 450, "Pmax": 10000,

        "nr1": [0.025377604, 0.97089776, -0.76128126, -1.2517125, 0.28005904],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.135, 1, 1, 0.42],

        "nr2": [-1.7144149, -0.64918553, 1.1662335, -0.35934725],
        "d2": [1, 3, 2, 2],
        "t2": [1.62, 2.35, 1.01, 2.65],
        "c2": [2, 2, 1, 2],
        "gamma2": [1]*4,

        "nr3": [1.4986537, -0.60326234, -0.11389713, -0.2553212, -0.1017598],
        "d3": [1, 1, 3, 2, 2],
        "t3": [0.75, 1.28, 1.5, 1.0, 1.9],
        "alfa3": [1.431, 1.803, 1.608, 1.837, 1.846],
        "beta3": [1.544, 1.366, 0.876, 1.117, 1.29],
        "gamma3": [1.265, 1.156, 0.916, 0.927, 0.926],
        "epsilon3": [0.781, 0.723, 0.842, 0.652, 1.139]}

    eq = (gao,)

    _surface = {"sigma": [0.04429], "exp": [1.242]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-8.2957, 4.5997, -4.4355, -5.0941, -4.1863],
        "t": [1, 1.5, 1.9, 4.3, 15.1]}
    _liquid_Density = {
        "eq": 1,
        "n": [7.2166, -18.074, 32.084, -30.238, 12.446],
        "t": [0.507, 0.824, 1.15, 1.5, 1.9]}
    _vapor_Density = {
        "eq": 2,
        "n": [-6.2029, 7.0601, -11.424, -24.160, -67.136, -182.16],
        "t": [0.496, 0.82, 1.17, 3.3, 6.8, 15.0]}

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

              "ek": 179, "sigma": 0.694, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [1.045], "psi_d": [0],
              "fint": [0.00125], "fint_t": [0],
              "chi": [1.99, -0.33], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.233e-9, "gam0": 0.061, "qd": 0.715e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""

    def test_gao(self):
        """Pag. 14"""
        st = C4F10(T=225, rhom=0)
        self.assertEqual(round(st.P.MPa, 6), 0)
        self.assertEqual(round(st.cvM.JmolK, 4), 171.2260)
        self.assertEqual(round(st.cpM.JmolK, 4), 179.5405)
        self.assertEqual(round(st.w, 5), 90.78029)
        self.assertEqual(round(st.hM.kJmol, 5), 14.88551)

        st = C4F10(T=225, rhom=7.8)
        self.assertEqual(round(st.P.MPa, 5), 37.55722)
        self.assertEqual(round(st.cvM.JmolK, 4), 192.9387)
        self.assertEqual(round(st.cpM.JmolK, 4), 246.5760)
        self.assertEqual(round(st.w, 4), 758.5512)
        self.assertEqual(round(st.hM.kJmol, 6), -8.747458)
        self.assertEqual(round(st.sM.JmolK, 5), -55.99954)

        st = C4F10(T=360, rhom=5.2)
        self.assertEqual(round(st.P.MPa, 6), 3.128110)
        self.assertEqual(round(st.cvM.JmolK, 4), 223.0894)
        self.assertEqual(round(st.cpM.JmolK, 4), 303.2828)
        self.assertEqual(round(st.w, 4), 226.8389)
        self.assertEqual(round(st.hM.kJmol, 5), 24.85197)
        self.assertEqual(round(st.sM.JmolK, 5), 77.45943)

        st = C4F10(T=387, rhom=2.637)
        self.assertEqual(round(st.P.MPa, 6), 2.355390)
        self.assertEqual(round(st.cvM.JmolK, 4), 236.7771)
        self.assertEqual(round(st.cpM.JmolK, 3), 8976.589)
        self.assertEqual(round(st.w, 5), 49.32548)
        self.assertEqual(round(st.hM.kJmol, 5), 37.43821)
        self.assertEqual(round(st.sM.JmolK, 4), 111.3057)

        st = C4F10(T=380, rhom=0.35)
        self.assertEqual(round(st.P.MPa, 7), 0.9312025)
        self.assertEqual(round(st.cvM.JmolK, 4), 218.1084)
        self.assertEqual(round(st.cpM.JmolK, 4), 236.8580)
        self.assertEqual(round(st.w, 5), 99.66618)
        self.assertEqual(round(st.hM.kJmol, 5), 44.83037)
        self.assertEqual(round(st.sM.JmolK, 4), 135.8261)

        st = C4F10(T=400, rhom=3.6)
        self.assertEqual(round(st.P.MPa, 6), 3.513083)
        self.assertEqual(round(st.cvM.JmolK, 4), 233.3552)
        self.assertEqual(round(st.cpM.JmolK, 4), 437.7846)
        self.assertEqual(round(st.w, 5), 89.44035)
        self.assertEqual(round(st.hM.kJmol, 5), 38.40444)
        self.assertEqual(round(st.sM.JmolK, 4), 112.8369)

    def test_Huber(self):
        """Table 7, pag 266"""
        st = C4F10(T=347.7, rhom=5.405)
        # self.assertEqual(round(st.mu.muPas, 4), 137.8044)
        self.assertEqual(round(st.mu.muPas, 4), 137.8051)
        self.assertEqual(round(st.k.mWmK, 4), 56.5125)
