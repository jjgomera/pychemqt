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


class Trans_2_butene(MEoS):
    """Multiparameter equations of state for trans-butene"""
    name = "trans-butene"
    CASNumber = "624-64-6"
    formula = "CH3-CH=CH-CH3"
    synonym = ""
    _refPropName = "T2BUTENE"
    _coolPropName = "trans-2-Butene"
    rhoc = unidades.Density(236.37592616)
    Tc = unidades.Temperature(428.61)
    Pc = unidades.Pressure(4027.3, "kPa")
    M = 56.10632  # g/mol
    Tt = unidades.Temperature(167.6)
    Tb = unidades.Temperature(274.03)
    f_acent = 0.21
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 26

    Fi1 = {"ao_log": [1, 2.9988],
           "pow": [0, 1],
           "ao_pow": [0.5917816, 2.1427758],
           "ao_exp": [5.3276, 13.29, 9.6745, 0.40087],
           "titao": [362/Tc, 1603/Tc, 3729/Tc, 4527/Tc]}

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

        "Tmin": Tt, "Tmax": 525.0, "Pmax": 50000.0, "rhomax": 13.141,

        "nr1": [0.81107, -2.8846, 1.0265, 0.016591, 0.086511, 0.00023256],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.12, 1.3, 1.74, 2.1, 0.28, 0.69],

        "nr2": [0.22654, -0.072182, -0.24849, -0.071374, -0.024737, 0.011843],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.75, 2., 4.4, 4.7, 15., 14.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = (lemmon, )
    _PR = [-0.1594, -16.6119]

    _surface = {
        "__doi__": {
            "autor": "Mulero, A., Cachadiña, I.",
            "title": "Recommended Correlations for the Surface Tension of "
                     "Several Fluids Included in the REFPROP Program",
            "ref": "J. Phys. Chem. Ref. Data 43(2) (2014) 023104",
            "doi": "10.1063/1.4878755"},
        "sigma": [0.0001859, 0.05539], "exp": [0.07485, 1.224]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.6226, 7.9421, -6.9631, -6.5517, 3.9584],
        "t": [1.0, 1.5, 1.65, 4.8, 5.3]}
    _liquid_Density = {
        "eq": 1,
        "n": [12.452, -34.419, 52.257, -42.889, 15.463],
        "t": [0.52, 0.73, 0.97, 1.24, 1.5]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.1276, -6.0548, -18.243, -60.842, 135.95, -182.70],
        "t": [0.412, 1.24, 3.2, 7.0, 10.0, 11.0]}

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

              "ek": 259, "sigma": 0.5508, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [1.12449, -0.147034, 0.036655],
              "psi_d": [0, 1, 2],
              "fint": [0.00102143, 6.64409e-7], "fint_t": [0, 1],
              "chi": [0.838527, 0.0648013], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.21e-9, "gam0": 0.057, "qd": 0.609e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""

    def test_shortLemmon(self):
        """Table 9, Pag 186"""
        st = Trans_2_butene(T=350, rho=0)
        self.assertEqual(round(st.P.MPa, 4), 0)
        self.assertEqual(round(st.hM.kJkmol, 0), 29959)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 89.965)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 98.279)
        self.assertEqual(round(st.w, 2), 238.03)

        st = Trans_2_butene(T=350, rho=0.3*Trans_2_butene.M)
        self.assertEqual(round(st.P.MPa, 5), 0.74692)
        self.assertEqual(round(st.hM.kJkmol, 0), 28521)
        self.assertEqual(round(st.sM.kJkmolK, 3), 86.364)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 95.429)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 112.42)
        self.assertEqual(round(st.w, 2), 208.86)

        st = Trans_2_butene(T=350, rho=10*Trans_2_butene.M)
        self.assertEqual(round(st.P.MPa, 3), 12.844)
        self.assertEqual(round(st.hM.kJkmol, 0), 10494)
        self.assertEqual(round(st.sM.kJkmolK, 3), 29.866)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 99.512)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 139.07)
        self.assertEqual(round(st.w, 2), 821.74)

        st = Trans_2_butene(T=440, rho=4*Trans_2_butene.M)
        self.assertEqual(round(st.P.MPa, 4), 4.7490)
        self.assertEqual(round(st.hM.kJkmol, 0), 29180)
        self.assertEqual(round(st.sM.kJkmolK, 3), 78.589)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 128.04)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 692.14)
        self.assertEqual(round(st.w, 2), 139.25)

    def test_Huber(self):
        """Table 7, pag 266"""
        st = Trans_2_butene(T=385.7, rhom=8.529)
        # self.assertEqual(round(st.mu.muPas, 5), 75.84518)
        self.assertEqual(round(st.mu.muPas, 5), 75.84560)
        self.assertEqual(round(st.k.mWmK, 4), 74.3078)
