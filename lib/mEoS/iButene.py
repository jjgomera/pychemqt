#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2023, Juan José Gómez Romera <jjgomera@gmail.com>

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


class iButene(MEoS):
    """Multiparameter equation of state for isobutene"""
    name = "isobutene"
    CASNumber = "115-11-7"
    formula = "CH2=C(CH3)2"
    synonym = ""
    _refPropName = "IBUTENE"
    _coolPropName = "IsoButene"
    rhoc = unidades.Density(233.9633544)
    Tc = unidades.Temperature(418.09)
    Pc = unidades.Pressure(4009.8, "kPa")
    M = 56.10632  # g/mol
    Tt = unidades.Temperature(132.4)
    Tb = unidades.Temperature(266.15)
    f_acent = 0.193
    momentoDipolar = unidades.DipoleMoment(0.5, "Debye")
    id = 27

    Fi1 = {"ao_log": [1, 3.],
           "pow": [0, 1],
           "ao_pow": [-0.12737888, 2.3125128],
           "ao_exp": [4.8924, 7.832, 7.2867, 8.7293],
           "titao": [399/Tc, 1270/Tc, 2005/Tc, 4017/Tc]}

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

        "Tmin": Tt, "Tmax": 550.0, "Pmax": 50000.0, "rhomax": 13.67,

        "nr1": [0.77111, -2.7971, 1.0118, 0.02073, 0.085086, 0.00021968],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.12, 1.3, 1.74, 2.1, 0.28, 0.69],

        "nr2": [0.20633, -0.078843, -0.23726, -0.080211, -0.027001, 0.013072],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.75, 2., 4.4, 4.7, 15., 14.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = (lemmon, )
    _PR = [-0.1815, -16.3103]

    _surface = {"sigma": [0.0545], "exp": [1.23]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.68973e1, 0.12475e1, -0.25441e1, -0.29282e1, 0.15778e1],
        "t": [1., 1.5, 3.16, 6.2, 7.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.62591e2, -0.20805e3, 0.33243e3, -0.29555e3, 0.11148e3],
        "t": [0.65, 0.8, 0.98, 1.16, 1.3]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.1841, -6.4014, -9.3817, -0.11160e2, -0.52298e2, -0.12195e3],
        "t": [0.431, 1.29, 3.3, 3.54, 7.3, 15.8]}

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

              "ek": 425, "sigma": 0.4776, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [1.12449, -0.147034, 0.036655],
              "psi_d": [0, 1, 2],
              "fint": [0.00102143, 6.64409e-7], "fint_t": [0, 1],
              "chi": [0.838527, 0.0648013], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.212e-9, "gam0": 0.058, "qd": 0.611e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""

    def test_shortLemmon(self):
        """Table 9, Pag 186"""
        st = iButene(T=350, rho=0)
        self.assertEqual(round(st.P.MPa, 4), 0)
        self.assertEqual(round(st.hM.kJkmol, 0), 29966)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 92.121)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 100.44)
        self.assertEqual(round(st.w, 2), 237.80)

        st = iButene(T=350, rho=0.3*iButene.M)
        self.assertEqual(round(st.P.MPa, 5), 0.75754)
        self.assertEqual(round(st.hM.kJkmol, 0), 28666)
        self.assertEqual(round(st.sM.kJkmolK, 3), 88.966)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 96.794)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 112.57)
        self.assertEqual(round(st.w, 2), 211.02)

        st = iButene(T=350, rho=10*iButene.M)
        self.assertEqual(round(st.P.MPa, 3), 17.776)
        self.assertEqual(round(st.hM.kJkmol, 0), 11782)
        self.assertEqual(round(st.sM.kJkmolK, 3), 32.951)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 101.72)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 139.45)
        self.assertEqual(round(st.w, 2), 838.25)

        st = iButene(T=440, rho=4*iButene.M)
        self.assertEqual(round(st.P.MPa, 4), 5.4086)
        self.assertEqual(round(st.hM.kJkmol, 0), 30169)
        self.assertEqual(round(st.sM.kJkmolK, 3), 82.345)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 127.28)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 407.00)
        self.assertEqual(round(st.w, 2), 151.13)

    def test_Huber(self):
        """Table 7, pag 266"""
        st = iButene(T=376.3, rhom=8.515)
        # self.assertEqual(round(st.mu.muPas, 5), 73.80192)
        # self.assertEqual(round(st.k.mWmK, 4), 73.6222)
        self.assertEqual(round(st.mu.muPas, 5), 73.80193)
        self.assertEqual(round(st.k.mWmK, 4), 73.6223)
