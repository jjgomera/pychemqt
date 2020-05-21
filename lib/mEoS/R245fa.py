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
from lib.mEoS import C3


class R245fa(MEoS):
    """Multiparameter equation of state for R245fa"""
    name = "1,1,1,3,3-pentafluoropropane"
    CASNumber = "460-73-1"
    formula = "CF3CH2CHF2"
    synonym = "R245fa"
    _refPropName = "R245FA"
    _coolPropName = "R245fa"
    rhoc = unidades.Density(519.4357675)
    Tc = unidades.Temperature(427.01)
    Pc = unidades.Pressure(3651.0, "kPa")
    M = 134.04794  # g/mol
    Tt = unidades.Temperature(170.0)
    Tb = unidades.Temperature(288.198)
    f_acent = 0.3783
    momentoDipolar = unidades.DipoleMoment(1.549, "Debye")

    Fi1 = {"ao_log": [1, 3.],
           "pow": [0, 1],
           "ao_pow": [-13.38560883, 9.845374371],
           "ao_exp": [5.5728, 10.385, 12.554],
           "titao": [222/Tc, 1010/Tc, 2450/Tc]}

    Fi2 = {"ao_log": [1, 3.],
           "pow": [0, 1],
           "ao_pow": [-13.4283638514, 9.87236538],
           "ao_exp": [5.5728, 10.385, 12.554],
           "titao": [222/427.16, 1010/427.16, 2450/427.16]}

    akasaka = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-245fa of Akasaka "
                    "(2015)",
        "__doi__": {"autor": "Akasaka, R.; Zhou, Y.; Lemmon, E.W.",
                    "title": "Fundamental Equation of State for 1,1,1,3,3-"
                             "Pentafluoropropane (R-245fa)",
                    "ref": "J. Phys. Chem. Ref. Data 44(1) (2015) 013104",
                    "doi": "10.1063/1.4913493"},

        "R": 8.3144621,
        "cp": Fi1,
        "ref": {"Tref": 273.15, "Pref": 1., "ho": 54490.07, "so": 268.2193},

        "Tmin": Tt, "Tmax": 440.0, "Pmax": 200000.0, "rhomax": 12.3,

        "nr1": [0.057506623, 1.5615975, -2.3614485, -0.51773521, 0.18509788],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.27, 0.9, 1.09, 0.4],

        "nr2": [-.87405626, -.27530955, .57971151, -.39934306, -.033230277],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.9, 1.7, 0.8, 3.6, 1.05],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [0.83210508, -0.33544300, -0.10117801, -0.0091495867],
        "d3": [1, 1, 3, 3],
        "t3": [1.8, 4.0, 4.5, 2.0],
        "alfa3": [1.011, 1.447, 1.079, 7.86],
        "beta3": [1.879, 2.454, 1.256, 21.1],
        "gamma3": [1.081, 0.651, 0.468, 1.293],
        "epsilon3": [0.709, 0.939, 0.703, 0.777]}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-245fa of Lemmon "
                    "and Span (2006).",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 "
                             "Industrial Fluids",
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi":  "10.1021/je050186n"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "NBP",
        "Tc": 427.16, "rhoc": 3.85,

        "Tmin": 171.05, "Tmax": 440.0, "Pmax": 200000.0, "rhomax": 12.3,

        "nr1": [1.2904, -3.2154, 0.50693, 0.093148, 0.00027638],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [.71458, .87252, -.015077, -.40645, -.11701, -.13062, -.022952],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    eq = akasaka, lemmon
    _PR = [-0.0436, -19.1110]

    _surface = {"sigma": [0.073586, 0.0103, -0.02663],
                "exp": [1.0983, 0.60033, 0.72765]}

    # Section 4 in Akasaka
    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.8353, 1.7746, -3.1305,  -3.4216],
        "t": [1, 1.5, 2.5, 5]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.46367, 2.2375, -0.27579, 0.55136],
        "t": [0.17, 0.5, 1.3, 2.5]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.99583, -2.6109, -4.4141, -18.573, -55.961],
        "t": [0.24, 0.61, 1, 2.7, 5.95]}

    trnECS = {"__name__": "Huber (2003)",

              "__doi__": {
                  "autor": "Huber, M.L., Laesecke, A., Perkins, R.A.",
                  "title": "Model for the Viscosity and Thermal Conductivity "
                           "of Refrigerants, Including a New Correlation for "
                           "the Viscosity of R134a",
                  "ref": "Ind. Eng. Chem. Res., 42(13) (2003) 3163-3178",
                  "doi": "10.1021/ie0300880"},

              "eq": "ecs",

              "ref": C3,
              "visco": "visco1",
              "thermo": "thermo0",

              "ek": 329.72, "sigma": 0.5529, "omega": 5,

              "psi": [1.1529, -4.4154e-2], "psi_d": [0, 1],
              "fint": [1.64999e-3, -3.28868e-7], "fint_t": [0, 1],
              "chi": [1.1627, -4.73491e-2], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
              "Xio": 0.194e-9, "gam0": 0.0496, "qd": 5e-10, "Tcref": 1.5*Tc}

    _viscosity = trnECS,
    _thermal = trnECS,


class Test(TestCase):

    def test_Akasaka(self):
        # Table 7, Pag 12
        st = R245fa(T=250, x=0.5)
        self.assertEqual(round(st.P.MPa, 8), 0.01646009)
        self.assertEqual(round(st.Liquido.rhoM, 5), 10.90057)
        self.assertEqual(round(st.Liquido.hM.kJmol, 3), 22.962)
        self.assertEqual(round(st.Liquido.sM.kJmolK, 6), 0.119346)
        self.assertEqual(round(st.Liquido.cvM.kJmolK, 6), 0.114155)
        self.assertEqual(round(st.Liquido.cpM.kJmolK, 6), 0.163322)
        self.assertEqual(round(st.Liquido.w, 3), 873.234)
        self.assertEqual(round(st.Gas.rhoM, 9), 0.008011195)
        self.assertEqual(round(st.Gas.hM.kJmol, 4), 51.9442)
        self.assertEqual(round(st.Gas.sM.kJmolK, 6), 0.235275)
        self.assertEqual(round(st.Gas.cvM.kJmolK, 6), 0.095084)
        self.assertEqual(round(st.Gas.cpM.kJmolK, 5), 0.10397)
        self.assertEqual(round(st.Gas.w, 3), 128.702)

        st = R245fa(T=400, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 2.210563)
        self.assertEqual(round(st.Liquido.rhoM, 5), 7.1518)
        self.assertEqual(round(st.Liquido.hM.kJmol, 4), 51.5591)
        self.assertEqual(round(st.Liquido.sM.kJmolK, 6), 0.206989)
        self.assertEqual(round(st.Liquido.cvM.kJmolK, 6), 0.147755)
        self.assertEqual(round(st.Liquido.cpM.kJmolK, 6), 0.255135)
        self.assertEqual(round(st.Liquido.w, 3), 235.471)
        self.assertEqual(round(st.Gas.rhoM, 6), 1.068455)
        self.assertEqual(round(st.Gas.hM.kJmol, 4), 65.2749)
        self.assertEqual(round(st.Gas.sM.kJmolK, 6), 0.241278)
        self.assertEqual(round(st.Gas.cvM.kJmolK, 6), 0.148819)
        self.assertEqual(round(st.Gas.cpM.kJmolK, 6), 0.234947)
        self.assertEqual(round(st.Gas.w, 3), 106.449)

        st = R245fa(T=250, rhom=11)
        self.assertEqual(round(st.P.MPa, 6), 7.454017)
        self.assertEqual(round(st.hM.kJmol, 4), 23.3683)
        self.assertEqual(round(st.sM.kJmolK, 6), 0.118254)
        self.assertEqual(round(st.cvM.kJmolK, 6), 0.114536)
        self.assertEqual(round(st.cpM.kJmolK, 6), 0.162059)
        self.assertEqual(round(st.w, 3), 908.590)

        st = R245fa(T=250, rhom=0.005)
        self.assertEqual(round(st.P.MPa, 8), 0.01031829)
        self.assertEqual(round(st.hM.kJmol, 4), 51.9792)
        self.assertEqual(round(st.sM.kJmolK, 6), 0.239262)
        self.assertEqual(round(st.cvM.kJmolK, 6), 0.094902)
        self.assertEqual(round(st.cpM.kJmolK, 6), 0.103569)
        self.assertEqual(round(st.w, 3), 129.147)

        st = R245fa(T=400, rhom=9)
        self.assertEqual(round(st.P.MPa, 5), 33.14725)
        self.assertEqual(round(st.hM.kJmol, 4), 50.9308)
        self.assertEqual(round(st.sM.kJmolK, 6), 0.196139)
        self.assertEqual(round(st.cvM.kJmolK, 6), 0.144117)
        self.assertEqual(round(st.cpM.kJmolK, 6), 0.188485)
        self.assertEqual(round(st.w, 3), 593.178)

        st = R245fa(T=400, rhom=0.5)
        self.assertEqual(round(st.P.MPa, 6), 1.352988)
        self.assertEqual(round(st.hM.kJmol, 4), 67.9067)
        self.assertEqual(round(st.sM.kJmolK, 6), 0.250859)
        self.assertEqual(round(st.cvM.kJmolK, 6), 0.138842)
        self.assertEqual(round(st.cpM.kJmolK, 6), 0.162544)
        self.assertEqual(round(st.w, 3), 135.712)

    def test_shortLemmon(self):
        # Table 10, Pag 842
        st = R245fa(T=429, rhom=3, eq="lemmon")
        self.assertEqual(round(st.P.kPa, 3), 3737.844)
        self.assertEqual(round(st.hM.kJkmol, 3), 63909.822)
        self.assertEqual(round(st.sM.kJkmolK, 3), 235.875)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 172.283)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 1891.958)
        self.assertEqual(round(st.w, 3), 78.673)
