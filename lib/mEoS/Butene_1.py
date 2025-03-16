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


class Butene_1(MEoS):
    """Multiparameter equation of state for 1-butene"""
    name = "butene"
    CASNumber = "106-98-9"
    formula = "CH3-CH2-CH=CH2"
    synonym = ""
    _refPropName = "1BUTENE"
    _coolPropName = "1-Butene"
    rhoc = unidades.Density(237.8907968)
    Tc = unidades.Temperature(419.29)
    Pc = unidades.Pressure(4005.1, "kPa")
    M = 56.10632  # g/mol
    Tt = unidades.Temperature(87.8)
    Tb = unidades.Temperature(266.84)
    f_acent = 0.192
    momentoDipolar = unidades.DipoleMoment(0.339, "Debye")
    id = 24

    Fi1 = {"ao_log": [1, 2.9197],
           "pow": [0, 1],
           "ao_pow": [-0.00101126, 2.3869174],
           "ao_exp": [2.9406, 6.5395, 14.535, 5.8971],
           "titao": [274/Tc, 951/Tc, 2127/Tc, 5752/Tc]}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for 1-butene of Lemmon "
                    "and Ihmels (2005)",
        "__doi__": {"autor": "Lemmon, E.W., Ihmels, E.C.",
                    "title": "Thermodynamic properties of the butenes: Part "
                             "II. Short fundamental equations of state",
                    "ref": "Fluid Phase Equilibria 228-229 (2005) 173-187",
                    "doi": "10.1016/j.fluid.2004.09.004"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 525., "Pmax": 70000.0, "rhomax": 14.59,

        "nr1": [0.78084, -2.8258, 0.99403, 0.017951, 0.088889, 0.00024673],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.12, 1.3, 1.74, 2.1, 0.28, 0.69],

        "nr2": [0.22846, -0.074009, -0.22913, -0.062334, -0.025385, 0.011040],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.75, 2., 4.4, 4.7, 15., 14.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = (lemmon, )
    _PR = [-0.2050, -16.0087]

    _surface = {"sigma": [0.05644], "exp": [1.248]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.71727e1, 0.26360e1, -0.20781e1, -0.28860e1, -0.13041e1],
        "t": [1, 1.5, 2, 4.35, 16.]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.16857e2, -0.46280e2, 0.53727e2, -0.23314e2, 0.18889e1],
        "t": [0.547, 0.73, 0.92, 1.14, 2.1]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.1106, -6.3103, -19.272, -48.739, -99.898, -190.01],
        "t": [0.415, 1.27, 3.34, 7.0, 14.5, 28.0]}

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

              "ek": 319, "sigma": 0.5198, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [1.12449, -0.147034, 0.036655],
              "psi_d": [0, 1, 2],
              "fint": [0.000900239, 1.13436e-6], "fint_t": [0, 1],
              "chi": [0.838527, 0.0648013], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.211e-9, "gam0": 0.057, "qd": 0.607e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""
    def test_shortLemmon(self):
        """Table 9, Pag 186"""
        st = Butene_1(T=350, rho=0)
        self.assertEqual(round(st.P.MPa, 3), 0)
        self.assertEqual(round(st.hM.kJkmol, 0), 29617)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 88.208)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 96.522)
        self.assertEqual(round(st.w, 2), 238.24)

        st = Butene_1(T=350, rho=0.3*Butene_1.M)
        self.assertEqual(round(st.P.MPa, 5), 0.75679)
        self.assertEqual(round(st.hM.kJkmol, 0), 28321)
        self.assertEqual(round(st.sM.kJkmolK, 3), 87.626)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 92.719)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 108.45)
        self.assertEqual(round(st.w, 2), 211.38)

        st = Butene_1(T=350, rho=10*Butene_1.M)
        self.assertEqual(round(st.P.MPa, 3), 17.864)
        self.assertEqual(round(st.hM.kJkmol, 0), 11377)
        self.assertEqual(round(st.sM.kJkmolK, 3), 31.563)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 97.760)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 135.36)
        self.assertEqual(round(st.w, 2), 843.31)

        st = Butene_1(T=440, rho=4*Butene_1.M)
        self.assertEqual(round(st.P.MPa, 4), 5.3245)
        self.assertEqual(round(st.hM.kJkmol, 0), 29454)
        self.assertEqual(round(st.sM.kJkmolK, 3), 80.191)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 124.13)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 416.03)
        self.assertEqual(round(st.w, 2), 151.49)

    def test_Huber(self):
        """Table 7, pag 266"""
        st = Butene_1(T=377.4, rhom=8.494)
        # self.assertEqual(round(st.mu.muPas, 5), 73.68086)
        self.assertEqual(round(st.mu.muPas, 5), 73.68090)
        self.assertEqual(round(st.k.mWmK, 4), 72.9207)
