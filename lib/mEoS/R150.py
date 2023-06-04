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


class R150(MEoS):
    """Multiparamenter equation of state for R150"""
    name = "1,2-dichloroethane"
    CASNumber = "107-06-2"
    formula = "CH2ClCH2Cl"
    synonym = "R-150"
    _refPropName = ""
    _coolPropName = ""
    rhoc = unidades.Density(428.49247)
    Tc = unidades.Temperature(561.6)
    Pc = unidades.Pressure(5254.835, "kPa")
    M = 98.959  # g/mol
    Tt = unidades.Temperature(237.52)
    Tb = unidades.Temperature(356.65)
    f_acent = 0.2708
    momentoDipolar = unidades.DipoleMoment(0, "Debye")
    id = 127

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [15.963798537, 0.972870308],
           "ao_exp": [5.35, 10.05],
           "titao": [22.5/Tc, 2015/Tc]}

    thol = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R150 of Thol (2021).",
        "__doi__": {
            "autor": "Thol, M., Rutkai, G., Köster, A., Miroshnichenko, S., "
                     "Wagner, W., Vrabec, J., Span, R.",
            "title": "Equation of state for 1,2-dichloroethane based on a "
                     "hybrid data set",
            "ref": "Molecular Physics 115(9-12) (2017) 1166-1185",
            "doi": "10.1080/00268976.2016.1262557"},

        "R": 8.3144621,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 1000, "Pmax": 1200000,

        "nr1": [0.051, 1.99, -2.595, -0.6653, 0.23595],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.352, 0.89, 0.824, 0.498],

        "nr2": [-0.17e1, -0.4453, 0.672474, -0.21918, -0.03554],
        "d2": [1, 3, 2, 2, 7],
        "t2": [1.63, 4.07, 0.679, 2.85, 1.07],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [0.9765, -0.495179, -0.23291174, -0.01090245, 0.39209],
        "d3": [1, 1, 3, 3, 1],
        "t3": [1.7, 2.09, 1.93, 3.72, 1.58],
        "alfa3": [0.66, 1.36, 0.711, 1.7, 1.11],
        "beta3": [0.574, 1.8, 0.462, 3.22, 2.22],
        "gamma3": [0.995, 0.329, 0.525, 0.85, 0.585],
        "epsilon3": [0.571, 0.862, 0.597, 1.16, 0.208]}

    eq = (thol,)

    _surface = {
        "__doi__": {"autor": "Huber, M.L.",
                    "title": "Models for Viscosity, Thermal Conductivity, and "
                             "Surface Tension of Selected Pure Fluids as "
                             "Implemented in REFPROP v10.0",
                    "ref": "NISTIR 8209",
                    "doi": "10.6028/NIST.IR.8209"},
        "sigma": [0.0785663], "exp": [1.19315]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-8.98372, 15.46, -37.11, 40.852, -20.042],
        "t": [1, 1.5, 1.9, 2.3, 2.8]}

    _liquid_Density = {
        "eq": 1,
        "n": [1.70532, 0.1786, 1.479, -0.62248],
        "t": [0.3, 0.7, 1.1, 1.5]}

    _vapor_Density = {
        "eq": 2,
        "n": [-2.93901, -6.45628, -49.73, 73.273, -83.717, -96.68],
        "t": [0.37, 1.2, 3.5, 4.3, 5.4, 13]}

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

              "ek": 445.96, "sigma": 0.4963, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 0.96,

              "psi": [0.766881, 0.136995, -0.0161687], "psi_d": [0, 1, 2],
              "fint": [9.18633e-4, 7.08996e-7], "fint_t": [0, 1],
              "chi": [1.35752, -0.116398], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.204e-9, "gam0": 0.056, "qd": 0.603e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""
    def test_thol(self):
        """Pag. 2 in suplementary file"""
        st = R150(T=250, rhom=0.0001)
        self.assertEqual(round(st.P.MPa, 11), 2.0782423e-4)
        self.assertEqual(round(st.hM.Jmol, 3), 23536.919)
        self.assertEqual(round(st.sM.JmolK, 5), 112.63617)
        self.assertEqual(round(st.cpM.JmolK, 6), 79.434919)
        self.assertEqual(round(st.w, 5), 153.14712)
        self.assertEqual(round(st.aM.Jmol, 4), -6700.3659)

        st = R150(T=250, rhom=14)
        self.assertEqual(round(st.P.MPa, 5), 131.48464)
        self.assertEqual(round(st.hM.Jmol, 4), -6329.4127)
        self.assertEqual(round(st.sM.JmolK, 6), -54.436721)
        self.assertEqual(round(st.cpM.JmolK, 5), 123.57319)
        self.assertEqual(round(st.w, 4), 1768.2675)
        self.assertEqual(round(st.aM.Jmol, 4), -2111.9925)

        st = R150(T=400, rhom=0.05)
        self.assertEqual(round(st.P.MPa, 8), 0.16082797)
        self.assertEqual(round(st.hM.Jmol, 3), 35896.016)
        self.assertEqual(round(st.sM.JmolK, 6), 96.321504)
        self.assertEqual(round(st.cpM.JmolK, 6), 93.603168)
        self.assertEqual(round(st.w, 5), 187.06153)
        self.assertEqual(round(st.aM.Jmol, 4), -5849.1451)

        st = R150(T=400, rhom=12)
        self.assertEqual(round(st.P.MPa, 6), 72.350760)
        self.assertEqual(round(st.hM.Jmol, 4), 9214.3686)
        self.assertEqual(round(st.sM.JmolK, 7), 8.2503437)
        self.assertEqual(round(st.cpM.JmolK, 5), 130.62037)
        self.assertEqual(round(st.w, 4), 1181.0220)
        self.assertEqual(round(st.aM.Jmol, 5), -114.99891)

        st = R150(T=550, rhom=14)
        self.assertEqual(round(st.P.MPa, 5), 744.15061)
        self.assertEqual(round(st.hM.Jmol, 3), 68017.209)
        self.assertEqual(round(st.sM.JmolK, 6), 24.569494)
        self.assertEqual(round(st.cpM.JmolK, 5), 139.11233)
        self.assertEqual(round(st.w, 4), 2169.7534)
        self.assertEqual(round(st.aM.Jmol, 4), 1350.3729)

    def test_Huber(self):
        """Table 7, pag 266"""
        st = R150(T=505.4, rhom=8.821)
        # self.assertEqual(round(st.mu.muPas, 4), 122.6713)
        # self.assertEqual(round(st.k.mWmK, 4), 83.6726)
        self.assertEqual(round(st.mu.muPas, 4), 122.6731)
        self.assertEqual(round(st.k.mWmK, 4), 83.6747)

    def test_Surface(self):
        """Table 10, pag 271"""
        self.assertEqual(round(R150(T=505.4, x=0.5).sigma, 7), 0.0050403)
