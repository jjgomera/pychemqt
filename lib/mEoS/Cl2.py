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
from lib.mEoS.C3 import C3


class Cl2(MEoS):
    """Multiparamenter equation of state for Chlorine"""
    name = "chlorine"
    CASNumber = "7782-50-5"
    formula = "Cl2"
    synonym = "Chlorine"
    _refPropName = "chlorine"
    _coolPropName = ""
    rhoc = unidades.Density(571.50236)
    Tc = unidades.Temperature(416.8654)
    Pc = unidades.Pressure(7642.4, "kPa")
    M = 70.906  # g/mol
    Tt = unidades.Temperature(172.17)
    Tb = unidades.Temperature(239.198)
    f_acent = 0.07
    momentoDipolar = unidades.DipoleMoment(0, "Debye")
    id = 105

    Fi1 = {"ao_log": [1, 2.5],
           "pow": [0, 1],
           "ao_pow": [-3.95390162183112, 3.83990498815427],
           "ao_exp": [1.0256, 0.067756, 0.14068],
           "titao": [800/Tc, 3000/Tc, 8200/Tc]}

    thol = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for chlorine of Thol (2021).",
        "__doi__": {
            "autor": "Thol, M., Herrig, S., Span, R., Lemmon, E.W.",
            "title": "A fundamental equation of state for the calculation of "
                     "thermodynamic proeprties of chlorine",
            "ref": "AIChE J. 67(9) (2021) 2633-2648",
            "doi": "10.1002/aic.17326"},

        "R": 8.314462618,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 450, "Pmax": 20000, "rhomax": 100,

        "nr1": [0.0245017, 0.9132904, -1.72309, -0.3359344, 0.1200495],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.196, 1, 1.08, 0.39],

        "nr2": [-1.214889, -0.10167, 0.6196819, -0.6578512, -0.009159452],
        "d2": [1, 3, 2, 2, 7],
        "t2": [1.64, 3.2, 1.32, 2.163, 0.93],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [1.909418, -0.07163412, -0.1893345, -0.5698469, -0.8964496],
        "d3": [1, 1, 3, 2, 2],
        "t3": [0.872, 2.08, 1.6, 1.37, 1.05],
        "alfa3": [0.969, 1.89, 1.32, 1.012, 0.98],
        "beta3": [1.22, 6.8, 3.5, 1.276, 1.6],
        "gamma3": [1.142, 1.22, 1.552, 1.135, 0.754],
        "epsilon3": [0.88, 0.73, 0.28, 0.863, 0.554]}

    eq = (thol,)

    _surface = {
        "__doi__": {"autor": "Huber, M.L.",
                    "title": "Models for Viscosity, Thermal Conductivity, and "
                             "Surface Tension of Selected Pure Fluids as "
                             "Implemented in REFPROP v10.0",
                    "ref": "NISTIR 8209",
                    "doi": "10.6028/NIST.IR.8209"},
        "sigma": [0.0783601], "exp": [1.28083]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-6.1289, 1.5112, -1.4523, -5.6038, 3.9923, -1.2651],
        "t": [1, 1.5, 2, 5.94, 7, 14.8]}

    _liquid_Density = {
        "eq": 1,
        "n": [0.9662, 1.7744, -0.23081, 0.47213],
        "t": [0.234, 0.68, 1.3, 3.35]}

    _vapor_Density = {
        "eq": 2,
        "n": [-1.7673, -5.173, -12.539, -37.552, -64.404, -151.49],
        "t": [0.3, 0.994, 2.7, 6.155, 12.4, 24]}

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

              "ek": 257, "sigma": 0.44, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [1.269, -0.08947], "psi_d": [0, 1],
              "fint": [0.0029], "fint_t": [0],
              "chi": [1.24341, -0.0812555], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.179e-9, "gam0": 0.056, "qd": 0.486e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""
    def test_thol(self):
        """Table 4, Pag. 3"""
        st = Cl2(T=300, rhom=0.001)
        self.assertEqual(round(st.P.MPa, 11), 0.00249358615)
        self.assertEqual(round(st.w, 5), 215.78489)
        self.assertEqual(round(st.cpM.JmolK, 7), 33.9883708)
        self.assertEqual(round(st.hM.Jmol, 3), 22546.353)
        self.assertEqual(round(st.sM.JmolK, 6), 123.917492)
        self.assertEqual(round(st.aM.Jmol, 4), -17122.4806)

        st = Cl2(T=300, rhom=20)
        self.assertEqual(round(st.P.MPa, 7), 11.8456782)
        self.assertEqual(round(st.w, 6), 876.623793)
        self.assertEqual(round(st.cpM.JmolK, 7), 67.8943166)
        self.assertEqual(round(st.hM.Jmol, 5), 4358.71206)
        self.assertEqual(round(st.sM.JmolK, 7), 14.1572969)
        self.assertEqual(round(st.aM.Jmol, 6), -480.760931)

        st = Cl2(T=350, rhom=0.01)
        self.assertEqual(round(st.P.MPa, 10), 0.0290386827)
        self.assertEqual(round(st.w, 6), 232.012075)
        self.assertEqual(round(st.cpM.JmolK, 7), 34.8576726)
        self.assertEqual(round(st.hM.Jmol, 4), 24246.9266)
        self.assertEqual(round(st.sM.JmolK, 6), 108.767988)
        self.assertEqual(round(st.aM.Jmol, 4), -16725.7376)

        st = Cl2(T=350, rhom=19)
        self.assertEqual(round(st.P.MPa, 7), 35.7673992)
        self.assertEqual(round(st.w, 6), 841.950732)
        self.assertEqual(round(st.cpM.JmolK, 7), 65.7421046)
        self.assertEqual(round(st.hM.Jmol, 4), 8154.4177)
        self.assertEqual(round(st.sM.JmolK, 6), 22.080586)
        self.assertEqual(round(st.aM.Jmol, 5), -1456.28211)

        st = Cl2(T=450, rhom=10)
        self.assertEqual(round(st.P.MPa, 7), 12.6189188)
        self.assertEqual(round(st.w, 6), 225.643929)
        self.assertEqual(round(st.cpM.JmolK, 6), 171.165598)
        self.assertEqual(round(st.hM.Jmol, 4), 17158.3805)
        self.assertEqual(round(st.sM.JmolK, 6), 47.843764)
        self.assertEqual(round(st.aM.Jmol, 5), -5633.20516)

    def test_Huber(self):
        """Table 7, pag 266"""
        st = Cl2(T=375.2, rhom=15.957)
        self.assertEqual(round(st.mu.muPas, 3), 160.877)
        self.assertEqual(round(st.k.mWmK, 4), 94.8362)

    def test_Surface(self):
        """Table 10, pag 271"""
        self.assertEqual(round(Cl2(T=375.2, x=0.5).sigma, 7), 0.0041019)
