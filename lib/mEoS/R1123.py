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
from lib.mEoS import R134a


class R1123(MEoS):
    """Multiparameter equation of state for R1123"""
    name = "Trifluoroethene"
    CASNumber = "359-11-5"
    formula = "CF2=CHF"
    synonym = "R-1123"
    rhoc = unidades.Density(504.20490885)
    Tc = unidades.Temperature(331.73)
    Pc = unidades.Pressure(4548.8, "kPa")
    M = 82.02455  # g/mol
    Tt = unidades.Temperature(195.15)
    Tb = unidades.Temperature(211.911)
    f_acent = 0.2429
    momentoDipolar = unidades.DipoleMoment(0, "Debye")

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1], "ao_pow": [-9.899717236, 7.265308393],
           "ao_exp": [7.203, 4.797], "titao": [784/Tc, 2251/Tc]}

    CP1 = {"ao": 3,
           "ao_exp": [5.39533, 7.79874], "exp": [453, 1712]}

    akasaka = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R1123 of Akasaka (2020)",
        "__doi__": {"autor": "Akasaka, R., Higashi, Y., Akoda, N., Fukuda, "
                             "S., Lemmon, E.W.",
                    "title": "Thermodynamic properties of trifluoroethene "
                             "(R1123): (p, ρ, T) behavior and fundamental "
                             "equation of state",
                    "ref": "Int. J. Refrig. 119 (2020) 457-467",
                    "doi": "10.1016/j.ijrefrig.2020.07.011"},

        "R": 8.314462618,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 480, "Pmax": 20000,

        "nr1": [0.03549472, 1.930575, -2.198695, -0.8067075, 0.2329929],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.244, 0.663, 1.157, 0.65],

        "nr2": [-1.62044, -1.457021, 0.58284441, -0.590613, -0.0209976],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2., 2.12, 1.115, 2.286, 0.84],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [0.4400995, -0.4425344, 2.210698, -0.46059068, -0.4125844,
                0.1220793, -0.3193233],
        "d3": [1, 1, 1, 1, 1, 1, 1],
        "t3": [1.26, 1.187, 1.29, 1.45, 1.935, 1.04, 1.92],
        "alfa3": [14.1, 14.1, 1.324, 1.904, 1.224, 1.275, 2.21],
        "beta3": [502.9, 502.6, 1.18, 1.16, 1.37, 1.85, 0.72],
        "gamma3": [1.0615, 1.0616, 1.293, 1.12, 1.05, 1.284, 1.27],
        "epsilon3": [1.0495, 1.0496, 0.8585, 1.145, 1.3, 1.55, 0.768]}

    akasaka2016 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R1123 of Akasaka (2016)",
        "__doi__": {"autor": "Akasaka, R., Fukushima, M., Lemmon, E.W.",
                    "title": "A Helmholtz Energy Equation of State for "
                             "Trifluoroethylene (R-1123)",
                    "ref": "International Refrigeration and Air Conditioning "
                           "Conference at Purdue, July 11-14, 2016, "
                           "http://docs.lib.purdue.edu/iracc/1698.",
                    "doi": ""},

        "R": 8.3144621,
        "cp": CP1,
        "ref": "IIR",
        "rhoc": 6,

        "Tmin": 250, "Tmax": 480, "Pmax": 10000,

        "nr1": [0.044649519, 2.0208994, -2.6417598, -0.41197275, 0.11153993],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.28, 0.782, 1.03, 0.68],

        "nr2": [-1.3190495, -.46464623, -.040932167, .26296637, -.018089075],
        "d2": [1, 3, 2, 2, 7],
        "t2": [1.64, 1.46, 2.23, 1.2, 1.73],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [1.6070681, -0.73580167, -0.26768005, -0.28256773, -0.14045846],
        "d3": [1, 1, 3, 2, 2],
        "t3": [1.05, 1.13, 1.78, 0.96, 1.85],
        "alfa3": [0.721, 1.261, 1.656, 0.804, 1.744],
        "beta3": [2.023, 1.705, 1.81, 3.1, 0.685],
        "gamma3": [1.09, 1.2, 0.9, 1.123, 0.837],
        "epsilon3": [0.557, 0.353, 0.291, 0.736, 1.131]}

    eq = (akasaka, akasaka2016)

    _surface = {
        "__doi__": {"autor": "Huber, M.L.",
                    "title": "Models for Viscosity, Thermal Conductivity, and "
                             "Surface Tension of Selected Pure Fluids as "
                             "Implemented in REFPROP v10.0",
                    "ref": "NISTIR 8209",
                    "doi": "10.6028/NIST.IR.8209"},
        "sigma": [0.0612], "exp": [1.26]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.416, 3.349, -2.855, -4.101, -4.173],
        "t": [1.0, 1.5, 1.84, 4.54, 14.3]}
    _liquid_Density = {
        "eq": 1,
        "n": [1.9578, 2.4060, -5.3582, 4.0164, 0.49887],
        "t": [0.337, 1.38, 2.15, 2.59, 13.6]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.0270, -8.3356, -29.515, -77.370, -222.98, -38.743],
        "t": [0.39, 1.37, 3.96, 8.72, 19.4, 29.4]}

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
#               "Tc": 424.964, "M": 82.02455, "rhoc": 6*82.02455,

              "ek": 263.4, "sigma": 0.4452, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [1], "psi_d": [0],
              "fint": [0.00125], "fint_t": [0],
              "chi": [1], "chi_d": [0],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.153e-9, "gam0": 0.075, "qd": 0.538e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""

    def test_akasaka(self):
        """Table 8, pag 466"""
        st = R1123(T=260, rhom=0)
        self.assertEqual(round(st.P.MPa, 6), 0)
        self.assertEqual(round(st.cvM.JmolK, 4), 54.9842)
        self.assertEqual(round(st.cpM.JmolK, 4), 63.2986)
        self.assertEqual(round(st.w, 3), 174.185)

        st = R1123(T=260, rhom=15)
        self.assertEqual(round(st.P.MPa, 5), 13.11323)
        self.assertEqual(round(st.cvM.JmolK, 4), 67.9165)
        self.assertEqual(round(st.cpM.JmolK, 3), 108.577)
        self.assertEqual(round(st.w, 3), 677.374)

        st = R1123(T=260, rhom=0.2)
        self.assertEqual(round(st.P.MPa, 7), 0.4003891)
        self.assertEqual(round(st.cvM.JmolK, 4), 57.2899)
        self.assertEqual(round(st.cpM.JmolK, 4), 69.1233)
        self.assertEqual(round(st.w, 3), 164.681)

        st = R1123(T=320, rhom=11)
        self.assertEqual(round(st.P.MPa, 6), 5.456590)
        self.assertEqual(round(st.cvM.JmolK, 4), 74.3579)
        self.assertEqual(round(st.cpM.JmolK, 3), 158.839)
        self.assertEqual(round(st.w, 3), 296.996)

        st = R1123(T=320, rhom=2)
        self.assertEqual(round(st.P.MPa, 6), 3.235708)
        self.assertEqual(round(st.cvM.JmolK, 4), 75.2904)
        self.assertEqual(round(st.cpM.JmolK, 3), 148.488)
        self.assertEqual(round(st.w, 3), 134.836)

        st = R1123(T=332, rhom=6.2)
        self.assertEqual(round(st.P.MPa, 6), 4.576478)
        self.assertEqual(round(st.cvM.JmolK, 4), 87.7735)
        self.assertEqual(round(st.cpM.JmolK, 1), 17935.9)
        self.assertEqual(round(st.w, 3), 112.642)

#     def test_Huber(self):
#         """Table 7, pag 266"""
#         st = R1123(T=298.6, rhom=12.527)
#         self.assertEqual(round(st.mu.muPas, 4), 115.7717)
#         self.assertEqual(round(st.k.mWmK, 4), 73.2468)

    def test_Surface(self):
        """Table 10, pag 271"""
        self.assertEqual(round(R1123(T=298.6, x=0.5).sigma, 7), 0.0033577)
