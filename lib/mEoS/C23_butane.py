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


class C23_butane(MEoS):
    """Multiparameter equation of state for 2,3-dimethylbutane"""
    name = "2,3-dimethylbutane"
    CASNumber = "79-29-8"
    formula = "CH3-CH(CH3)-CH(CH3)-CH3"
    synonym = ""
    _refPropName = ""
    _coolPropName = ""
    rhoc = unidades.Density(241.291008)
    Tc = unidades.Temperature(500.6)
    Pc = unidades.Pressure(3161.0, "kPa")
    M = 86.17536  # g/mol
    Tt = unidades.Temperature(145.05)
    Tb = unidades.Temperature(331.177)
    f_acent = 0.247
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 55

    Fi1 = {"ao_log": [1, 7.5],
           "pow": [0, 1],
           "ao_pow": [4.13396107278, -1.39038584556],
           "ao_exp": [7.2248, 26.48, 17.893],
           "titao": [535/Tc, 1693/Tc, 4369/Tc]}

    gao = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for 2,3-DiMethylbutane of Gao"
                    "et al. (2021)",
        "__doi__": {
            "autor": "Gao, K., Wu, J., Lemmon, E.W.",
            "title": "Equations of State for the Thermodynamic Properties of "
                     "Three Hexane Isomers: 3-Methylpantane, "
                     "2,2-Dimethylbutane, and 2,3-Dimethylbutane",
            "ref": "J. Phys. Chem. Ref. Data 50(3) (2021) 033103",
            "doi": "10.1063/1.5093644"},

        "R": 8.314462618,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 550, "Pmax": 1100000, "rhomax": 9.12,

        "nr1": [0.007194931, 0.97492236, -0.6880694, -1.0251264, 0.20316871],
        "d1": [5, 1, 1, 2, 3],
        "t1": [1.0, 0.25, 1.0, 1.0, 0.364],

        "nr2": [-1.7247168, -0.75015882, 1.0712536, -0.30179884, -0.029150384],
        "d2": [1, 3, 2, 2, 8],
        "t2": [1.37, 1.77, 1.0, 2.3, 1.5],
        "c2": [2, 2, 1, 2, 2],
        "gamma2": [1]*5,

        "nr3": [0.69062135, -0.010365966, -0.21028711, -.11507614, -.19936642],
        "d3": [1, 1, 3, 2, 1],
        "t3": [0.43, 1.15, 0.4, 0.8, 1.375],
        "alfa3": [1.432, 1.787, 1.412, 1.542, 1.2],
        "beta3": [2.3, 0.46, 1.9, 1.0, 2.86],
        "gamma3": [1.16, 0.91, 1.0, 1.05, 1.06],
        "epsilon3": [0.689, 2.139, 0.314, 0.992, 0.619]}

    eq = (gao, )

    _surface = {
        "__doi__": {"autor": "Huber, M.L.",
                    "title": "Models for Viscosity, Thermal Conductivity, and "
                             "Surface Tension of Selected Pure Fluids as "
                             "Implemented in REFPROP v10.0",
                    "ref": "NISTIR 8209",
                    "doi": "10.6028/NIST.IR.8209"},
        "sigma": [0.05235], "exp": [1.24897]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.6041, 3.4948, -2.6831, -3.0964, -1.2618],
        "t": [1.0, 1.5, 1.9, 3.85, 17.3]}
    _liquid_Density = {
        "eq": 1,
        "n": [10.132, -29.743, 49.275, -39.267, 12.406],
        "t": [0.54, 0.818, 1.1, 1.4, 1.75]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.8824, -8.0209, -25.626, -56.727, -145.5],
        "t": [0.448, 1.55, 3.85, 7.85, 17.15]}

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

              "ek": 397.5, "sigma": 0.574, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [1.25327, -0.178928, 0.0348391], "psi_d": [0, 1, 2],
              "fint": [0.00116], "fint_t": [0],
              "chi": [1.00438, 0.014484], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.238e-9, "gam0": 0.058, "qd": 0.701e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""

    def test_Gao(self):
        """Table 13, pag 22"""
        st = C23_butane(T=260, rhom=8.3)
        self.assertEqual(round(st.P.MPa, 5), 33.95611)
        self.assertEqual(round(st.cvM.JmolK, 4), 137.4007)
        self.assertEqual(round(st.cpM.JmolK, 4), 171.8992)
        self.assertEqual(round(st.w, 3), 1402.545)
        self.assertEqual(round(st.hM.kJmol, 5), -10.45592)
        self.assertEqual(round(st.sM.JmolK, 5), -50.03608)

        st = C23_butane(T=460, rhom=5.5)
        self.assertEqual(round(st.P.MPa, 6), 3.333684)
        self.assertEqual(round(st.cvM.JmolK, 4), 205.0061)
        self.assertEqual(round(st.cpM.JmolK, 4), 274.9322)
        self.assertEqual(round(st.w, 4), 394.4134)
        self.assertEqual(round(st.hM.kJmol, 5), 30.36433)
        self.assertEqual(round(st.sM.JmolK, 5), 75.59279)

        st = C23_butane(T=501, rhom=2.8)
        self.assertEqual(round(st.P.MPa, 6), 3.180216)
        self.assertEqual(round(st.cvM.JmolK, 4), 227.2418)
        self.assertEqual(round(st.cpM.JmolK, 2), 18442.03)
        self.assertEqual(round(st.w, 5), 86.22365)
        self.assertEqual(round(st.hM.kJmol, 5), 46.71702)
        self.assertEqual(round(st.sM.JmolK, 4), 109.3744)

        st = C23_butane(T=260, rhom=0)
        self.assertEqual(round(st.P.MPa, 5), 0)
        self.assertEqual(round(st.cvM.JmolK, 4), 118.9821)
        self.assertEqual(round(st.cpM.JmolK, 4), 127.2965)
        self.assertEqual(round(st.w, 4), 163.8248)
        self.assertEqual(round(st.hM.kJmol, 5), 17.84952)

        st = C23_butane(T=480, rhom=0.37)
        self.assertEqual(round(st.P.MPa, 6), 1.223390)
        self.assertEqual(round(st.cvM.JmolK, 4), 205.9278)
        self.assertEqual(round(st.cpM.JmolK, 4), 224.3395)
        self.assertEqual(round(st.w, 4), 183.0561)
        self.assertEqual(round(st.hM.kJmol, 5), 52.98636)
        self.assertEqual(round(st.sM.JmolK, 4), 127.3124)

        st = C23_butane(T=520, rhom=3.8)
        self.assertEqual(round(st.P.MPa, 6), 4.763216)
        self.assertEqual(round(st.cvM.JmolK, 4), 229.0952)
        self.assertEqual(round(st.cpM.JmolK, 4), 402.1148)
        self.assertEqual(round(st.w, 4), 161.0000)
        self.assertEqual(round(st.hM.kJmol, 5), 48.68907)
        self.assertEqual(round(st.sM.JmolK, 4), 112.3156)

    def test_Huber(self):
        """Table 7, pag 266"""
        st = C23_butane(T=450.5, rhom=5.652)
        # self.assertEqual(round(st.mu.muPas, 5), 90.14569)
        self.assertEqual(round(st.mu.muPas, 5), 90.14598)
        self.assertEqual(round(st.k.mWmK, 4), 67.3762)

    def test_Surface(self):
        """Table 10, pag 271"""
        self.assertEqual(round(C23_butane(T=450.5, x=0.5).sigma, 7), 0.0029538)
