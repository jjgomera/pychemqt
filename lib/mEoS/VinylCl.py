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


class VinylCl(MEoS):
    """Multiparamenter equation of state for vinyl chloride"""
    name = "vinyl chloride"
    CASNumber = "75-01-4"
    formula = "CH2CHCl"
    synonym = "R-1140"
    _refPropName = ""
    _coolPropName = ""
    rhoc = unidades.Density(353.114943)
    Tc = unidades.Temperature(425)
    Pc = unidades.Pressure(5600.3, "kPa")
    M = 62.49822  # g/mol
    Tt = unidades.Temperature(119.31)
    Tb = unidades.Temperature(259.26)
    f_acent = 0.16
    momentoDipolar = unidades.DipoleMoment(0, "Debye")
    id = 122

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [-3.4387095412593389, 3.6032015766362435],
           "ao_exp": [4.51, 4.45, 3.04],
           "titao": [923/Tc, 1907/Tc, 4400/Tc]}

    thol = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for vinyl choride of Thol "
                    "(2022).",
        "__doi__": {
            "autor": "Thol, M., Fenkl, F., Lemmon, E.W.",
            "title": "A Fundamental Equation of State for Chloroethene for "
                     "Temperatures from the Triple Point to 430 K and "
                     "Pressures to 100 MPa",
            "ref": "Int. J. Themophysics 42 (2022) 41",
            "doi": "10.1007/s10765-021-02961-3"},

        "R": 8.314462618,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 530, "Pmax": 100000,

        "nr1": [0.012261627, 0.8604231, -1.0169596, -0.81353519, 0.23529717],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.136, 0.984, 1.02, 0.624],

        "nr2": [-0.87520652, -1.0092224, 0.51890592, -0.15911867, -0.007640016,
                -0.22858771],
        "d2": [1, 3, 2, 2, 7, 1],
        "t2": [1.215, 1.427, 0.9, 4, 0.6, 3.45],
        "c2": [2, 2, 1, 2, 1, 3],
        "gamma2": [1]*6,

        "nr3": [0.9093646, -0.19999899, -0.37519445, 0.97421121, -1.006613],
        "d3": [1, 2, 1, 1, 1],
        "t3": [1.42, 0.9, 1.475, 0.5, 1.156],
        "alfa3": [1.736, 30, 2.425, 0.714, 0.694],
        "beta3": [0.77, 815, 0.77, 0.41, 0.32],
        "gamma3": [1.3, 1.1, 1.59, 1.5, 0.93],
        "epsilon3": [0.82, 0.92, 0.995, 0.71, 0.66]}

    eq = (thol,)

    _surface = {
        "__doi__": {"autor": "Huber, M.L.",
                    "title": "Models for Viscosity, Thermal Conductivity, and "
                             "Surface Tension of Selected Pure Fluids as "
                             "Implemented in REFPROP v10.0",
                    "ref": "NISTIR 8209",
                    "doi": "10.6028/NIST.IR.8209"},
        "Tc": 424.964,
        "sigma": [0.0655789], "exp": [1.16473]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-6.945, 1.7737, -2.0553, -2.813, -2.756],
        "t": [1, 1.5, 2.78, 5.2, 15]}

    _liquid_Density = {
        "eq": 1,
        "n": [1.8041, 5.629, -18.450, 33.235, -66.136, 47.53],
        "t": [0.33, 1.19, 1.75, 2.39, 3.22, 3.46]}

    _vapor_Density = {
        "eq": 2,
        "n": [-0.5868, -5.3182, -14.6323, -46.665, -98.3277, -228.973],
        "t": [0.193, 0.687, 2.540, 6.12, 13, 25.5]}

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

              "ek": 337.46, "sigma": 0.455, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1.06,

              "psi": [0.991393, -0.0190085], "psi_d": [0, 1],
              "fint": [4.68338e-4, 1.55637e-6], "fint_t": [0, 1],
              "chi": [1], "chi_d": [0],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.195e-9, "gam0": 0.059, "qd": 0.551e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""
    def test_thol(self):
        """Table 5, Pag. 8"""
        st = VinylCl(T=250, rhom=0.03)
        self.assertEqual(round(st.P.MPa, 10), 0.0611336026)
        self.assertEqual(round(st.w, 6), 198.051441)
        self.assertEqual(round(st.cpM.JmolK, 7), 49.1470240)
        self.assertEqual(round(st.hM.Jmol, 4), 21830.1713)
        self.assertEqual(round(st.sM.JmolK, 7), 88.2755532)
        self.assertEqual(round(st.aM.Jmol, 5), -2276.50376)

        st = VinylCl(T=250, rhom=16)
        self.assertEqual(round(st.P.MPa, 7), 10.5895545)
        self.assertEqual(round(st.w, 5), 1155.49169)
        self.assertEqual(round(st.cpM.JmolK, 7), 91.1392445)
        self.assertEqual(round(st.hM.Jmol, 6), -477.949658)
        self.assertEqual(round(st.sM.JmolK, 8), -4.48949199)
        self.assertEqual(round(st.aM.Jmol, 7), -17.4238160)

        st = VinylCl(T=300, rhom=0.18)
        self.assertEqual(round(st.P.MPa, 10), 0.413797863)
        self.assertEqual(round(st.w, 6), 204.981105)
        self.assertEqual(round(st.cpM.JmolK, 7), 59.5414183)
        self.assertEqual(round(st.hM.Jmol, 4), 23878.6867)
        self.assertEqual(round(st.sM.JmolK, 7), 80.5342776)
        self.assertEqual(round(st.aM.Jmol, 5), -2580.47364)

        st = VinylCl(T=300, rhom=15)
        self.assertEqual(round(st.P.MPa, 7), 23.0374719)
        self.assertEqual(round(st.w, 5), 1008.04450)
        self.assertEqual(round(st.cpM.JmolK, 7), 91.4066946)
        self.assertEqual(round(st.hM.Jmol, 5), 4518.73773)
        self.assertEqual(round(st.sM.JmolK, 7), 10.8027936)
        self.assertEqual(round(st.aM.Jmol, 6), -257.931823)

        st = VinylCl(T=430, rhom=14)
        self.assertEqual(round(st.P.MPa, 7), 88.2661442)
        self.assertEqual(round(st.w, 6), 982.983659)
        self.assertEqual(round(st.cpM.JmolK, 7), 92.7700338)
        self.assertEqual(round(st.hM.Jmol, 4), 18735.1820)
        self.assertEqual(round(st.sM.JmolK, 7), 37.6281054)
        self.assertEqual(round(st.aM.Jmol, 5), -3749.62785)

    def test_Huber(self):
        """Table 7, pag 266"""
        st = VinylCl(T=382.5, rhom=11.537)
        # self.assertEqual(round(st.mu.muPas, 5), 90.53455)
        # self.assertEqual(round(st.k.mWmK, 4), 84.5575)
        self.assertEqual(round(st.mu.muPas, 5), 90.94743)
        self.assertEqual(round(st.k.mWmK, 4), 84.8306)

    def test_Surface(self):
        """Table 10, pag 271"""
        self.assertEqual(round(VinylCl(T=382.5, x=0.5).sigma, 7), 0.0044838)
