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
from lib.mEoS import R134a


class R13I1(MEoS):
    """Multiparameter equation of state for trifluoroiodomethane"""
    name = "trifluoroiodomethane"
    CASNumber = "2314-97-8"
    formula = "CF3I"
    synonym = "R13I1"
    _refPropName = "CF3I"
    _coolPropName = "R13I1"
    rhoc = unidades.Density(868.00061824)
    Tc = unidades.Temperature(396.44)
    Pc = unidades.Pressure(3953., "kPa")
    M = 195.9104  # g/mol
    Tt = unidades.Temperature(120.)
    Tb = unidades.Temperature(251.3)
    f_acent = 0.18
    momentoDipolar = unidades.DipoleMoment(0.92, "Debye")

    CP1 = {"ao": 4.,
           "ao_exp": [6.2641], "exp": [694]}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R13I1 of Lemmon and "
                    "Span (2015)",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Thermodynamic Properties of R-227ea, R-365mfc, "
                             "R-115, and R-13I1",
                    "ref": "J. Chem. Eng. Data, 60(12) (2015) 3745-3758",
                    "doi": "10.1021/acs.jced.5b00684"},

        "R": 8.3144621,
        "cp": CP1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 420., "Pmax": 20000.0, "rhomax": 14.1,

        "nr1": [1.12191, -3.08087, 1.11307, -0.184885, 0.110971, 0.000325],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.333357, -0.0288288, -0.371554, -0.0997985, -0.0333205,
                0.0207882],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = (lemmon, )

    _surface = {"sigma": [0.05767], "exp": [1.298]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-6.8642, 1.7877, -1.0619, -2.1677],
        "t": [1.0, 1.5, 1.9, 3.8]}
    _liquid_Density = {
        "eq": 1,
        "n": [2.0711, 1.562, -2.599, 1.7177],
        "t": [0.38, 1.3, 1.9, 2.5]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.0987, -6.8771, -19.701, -46.86, -100.02],
        "t": [0.41, 1.33, 3.5, 7.4, 16.0]}

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

              "ek": 314.8, "sigma": 0.4926, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 0.95,

              "psi": [1.22725, -0.0879263], "psi_d": [0, 1],
              "fint": [1.28541e-3, 5.32854e-7], "fint_t": [0, 1],
              "chi": [1], "chi_d": [0],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.21e-9, "gam0": 0.057, "qd": 0.598e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""
    def test_lemmon(self):
        """Table 7, Pag 3754"""
        st = R13I1(T=300, rhom=0)
        self.assertEqual(round(st.P.MPa, 6), 0.0)
        self.assertEqual(round(st.cvM.JmolK, 5), 58.90476)
        self.assertEqual(round(st.cpM.JmolK, 5), 67.21922)
        self.assertEqual(round(st.w, 4), 120.5370)

        st = R13I1(T=300, rhom=10.9)
        self.assertEqual(round(st.P.MPa, 5), 16.52331)
        self.assertEqual(round(st.cvM.JmolK, 5), 66.39803)
        self.assertEqual(round(st.cpM.JmolK, 4), 101.9090)
        self.assertEqual(round(st.w, 4), 529.4695)

        st = R13I1(T=397, rhom=4.4)
        self.assertEqual(round(st.P.MPa, 6), 3.991241)
        self.assertEqual(round(st.cvM.JmolK, 5), 90.32758)
        self.assertEqual(round(st.cpM.JmolK, 2), 11675.81)
        self.assertEqual(round(st.w, 5), 74.67026)

    def test_Huber(self):
        """Table 7, pag 266"""
        st = R13I1(T=356.8, rhom=8.724)
        # self.assertEqual(round(st.mu.muPas, 4), 168.5168)
        self.assertEqual(round(st.mu.muPas, 4), 168.5175)
        self.assertEqual(round(st.k.mWmK, 4), 40.8519)
