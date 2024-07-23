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


class C6F14(MEoS):
    """Multiparameter equation of state for perflurohexane"""
    name = "perfluorohexane"
    CASNumber = "355-42-0"
    formula = "C6F14"
    synonym = ""
    _refPropName = "C6F14"
    _coolPropName = ""
    rhoc = unidades.Density(616.92665)
    Tc = unidades.Temperature(448)
    Pc = unidades.Pressure(1741.6, "kPa")
    M = 338.042  # g/mol
    Tt = unidades.Temperature(187.07)
    Tb = unidades.Temperature(330.274)
    f_acent = 0.497
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 321

    CP1 = {"ao": 17,
           "ao_exp": [4.902, 23.43, 10.52], "exp": [433, 910, 1982]}

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

        "R": 8.314462618,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 560, "Pmax": 40000,

        "nr1": [0.035689273, 0.89834616, -0.11619207, -1.6558707, 0.3090401],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.146, 1, 1, 0.39],

        "nr2": [-3.0212885, -1.309987, 1.4611604, -0.63849402],
        "d2": [1, 3, 2, 2],
        "t2": [1.56, 2.25, 0.987, 2.602],
        "c2": [2, 2, 1, 2],
        "gamma2": [1]*4,

        "nr3": [-0.40480926, 2.3673483, 0.40072213, -0.43534683, -0.90267664],
        "d3": [1, 1, 1, 3, 2],
        "t3": [1.87, 0.97, 2.22, 1.7, 1.31],
        "alfa3": [0.8775, 1, 1.327, 1.102, 1.274],
        "beta3": [1.171, 1.14, 0.645, 0.658, 0.727],
        "gamma3": [1.254, 1.312, 1.178, 1.326, 0.902],
        "epsilon3": [0.36, 0.755, 1.693, 1.04, 0.646]}

    eq = (gao,)

    _surface = {
        "__doi__": {"autor": "Huber, M.L.",
                    "title": "Models for Viscosity, Thermal Conductivity, and "
                             "Surface Tension of Selected Pure Fluids as "
                             "Implemented in REFPROP v10.0",
                    "ref": "NISTIR 8209",
                    "doi": "10.6028/NIST.IR.8209"},
        "sigma": [0.0230631, 0.0703415], "exp": [0.98534, 2.6579]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-9.0231, 5.1365, -4.8413, -5.2906, -2.937],
        "t": [1, 1.5, 1.97, 3.63, 11.74]}
    _liquid_Density = {
        "eq": 1,
        "n": [4.5417, -5.3549, 5.7116, -2.6333, 1.1928],
        "t": [0.463, 0.85, 1.25, 1.75, 3.45]}
    _vapor_Density = {
        "eq": 2,
        "n": [-4.9298, -13.949, -53.702, -113.62, -262.16],
        "t": [0.476, 1.88, 4.73, 10.4, 21.05]}

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

              "ek": 160, "sigma": 0.805, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [0.673625, 0.35383, -0.0787347], "psi_d": [0, 1, 2],
              "fint": [0.00125], "fint_t": [0],
              "chi": [1.99965, -0.290494], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.254e-9, "gam0": 0.06, "qd": 0.812e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""
    def test_gao(self):
        """Pag. 14"""
        st = C6F14(T=260, rhom=0)
        self.assertEqual(round(st.P.MPa, 6), 0)
        self.assertEqual(round(st.cvM.JmolK, 4), 244.6528)
        self.assertEqual(round(st.cpM.JmolK, 4), 252.9673)
        self.assertEqual(round(st.w, 5), 81.31590)
        self.assertEqual(round(st.hM.kJmol, 6), 9.955560)

        st = C6F14(T=260, rhom=5.5)
        self.assertEqual(round(st.P.MPa, 5), 28.03371)
        self.assertEqual(round(st.cvM.JmolK, 4), 270.3971)
        self.assertEqual(round(st.cpM.JmolK, 4), 329.6023)
        self.assertEqual(round(st.w, 4), 730.8597)
        self.assertEqual(round(st.hM.kJmol, 5), -21.42844)
        self.assertEqual(round(st.sM.JmolK, 5), -91.21752)

        st = C6F14(T=410, rhom=3.7)
        self.assertEqual(round(st.P.MPa, 7), 0.9573522)
        self.assertEqual(round(st.cvM.JmolK, 4), 336.7461)
        self.assertEqual(round(st.cpM.JmolK, 4), 435.6546)
        self.assertEqual(round(st.w, 4), 181.2565)
        self.assertEqual(round(st.hM.kJmol, 5), 31.64675)
        self.assertEqual(round(st.sM.JmolK, 5), 85.06658)

        st = C6F14(T=448.5, rhom=1.825)
        self.assertEqual(round(st.P.MPa, 6), 1.758863)
        self.assertEqual(round(st.cvM.JmolK, 4), 363.5137)
        self.assertEqual(round(st.cpM.JmolK, 2), 16301.72)
        self.assertEqual(round(st.w, 5), 36.59010)
        self.assertEqual(round(st.hM.kJmol, 5), 53.05034)
        self.assertEqual(round(st.sM.JmolK, 4), 133.9994)

        st = C6F14(T=430, rhom=0.23)
        self.assertEqual(round(st.P.MPa, 7), 0.6728015)
        self.assertEqual(round(st.cvM.JmolK, 4), 330.4599)
        self.assertEqual(round(st.cpM.JmolK, 4), 352.0829)
        self.assertEqual(round(st.w, 5), 85.33926)
        self.assertEqual(round(st.hM.kJmol, 5), 58.08506)
        self.assertEqual(round(st.sM.JmolK, 4), 150.4846)

        st = C6F14(T=460, rhom=2.7)
        self.assertEqual(round(st.P.MPa, 6), 2.671262)
        self.assertEqual(round(st.cvM.JmolK, 4), 358.0500)
        self.assertEqual(round(st.cpM.JmolK, 4), 541.4677)
        self.assertEqual(round(st.w, 5), 81.00252)
        self.assertEqual(round(st.hM.kJmol, 5), 53.96219)
        self.assertEqual(round(st.sM.JmolK, 4), 135.1519)

    def test_Huber(self):
        """Table 7, pag 266"""
        st = C6F14(T=403.2, rhom=3.873)
        # self.assertEqual(round(st.mu.muPas, 4), 170.8344)
        self.assertEqual(round(st.mu.muPas, 4), 170.8353)
        self.assertEqual(round(st.k.mWmK, 4), 54.1785)

    def test_Surface(self):
        """Table 10, pag 271"""
        self.assertEqual(round(C6F14(T=403.2, x=0.5).sigma, 7), 0.0025401)
