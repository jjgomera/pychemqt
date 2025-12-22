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


class ClBenzene(MEoS):
    """Multiparameter equation of state for Chlorobenzene"""
    name = "Chlorobenzene"
    CASNumber = "108-90-7"
    formula = "C6H5Cl"
    synonym = ""
    _refPropName = "CHLOROBENZENE"
    _coolPropName = ""
    rhoc = unidades.Density(364.68468)
    Tc = unidades.Temperature(632.35)
    Pc = unidades.Pressure(4520.6, "kPa")
    M = 112.557  # g/mol
    Tt = unidades.Temperature(227.9)
    Tb = unidades.Temperature(404.87)
    f_acent = 0.2532
    momentoDipolar = unidades.DipoleMoment(0, "Debye")
    id = 172

    CP1 = {"ao": 4,
           "an": [], "pow": [],
           "ao_exp": [6.2566, 16.273, 7.6017],
           "exp": [4160.0, 1580.0, 600.0]}

    refprop = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for Chlorobenzene (2022)",
        "__doi__": {
            "autor": "Huber, M.L., Lemmon, E.W., Bell, I.H., McLinden, M.O.",
            "title": "The NIST REFPROP Database for Highly Accurate Properties"
                     " of Industrially Importants Fluids",
            "ref": "Ind. Eng. Chem. Res. 61(42) (2022) 15449-15472",
            "doi": "10.1021/acs.iecr.2c01427"},

        "R": 8.3144598,
        "cp": CP1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 700.0, "Pmax": 100000.0, "rhomax": 35.57,

        "nr1": [0.03675169, 1.2629, -2.092176, -0.5062699, 0.1826893],
        "d1": [4.0, 1.0, 1.0, 2.0, 3.0],
        "t1": [1.0, 0.25, 0.967, 1.06, 0.527,],

        "nr2": [-0.9710427, -0.3295967, 0.8757209, -0.3980378, -0.02049013],
        "d2": [1.0, 3.0, 2.0, 2.0, 7.0],
        "t2": [1.93, 2.44, 1.28, 3.06, 1.013],
        "c2": [2.0, 2.0, 1.0, 2.0, 1.0],
        "gamma2": [1]*5,

        "nr3": [1.307316, -0.07704369, -0.2117575, -0.5223262, -0.00306347],
        "d3": [1.0, 1.0, 2.0, 2.0, 3.0],
        "t3": [0.768, 1.4, 1.3, 1.16, 1.2],
        "alfa3": [0.815, 1.25, 0.876, 1.034, 3.76],
        "beta3": [1.45, 1.65, 1.51, 1.24, 62.7],
        "gamma3": [1.31, 1.0, 1.29, 1.114, 0.87],
        "epsilon3": [1.042, 1.638, 1.139, 0.799, 0.823]}

    eq = (refprop, )

    _surface = {
        "__doi__": {"autor": "Huber, M.L.",
                    "title": "Models for Viscosity, Thermal Conductivity, and "
                             "Surface Tension of Selected Pure Fluids as "
                             "Implemented in REFPROP v10.0",
                    "ref": "NISTIR 8209",
                    "doi": "10.6028/NIST.IR.8209"},
        "sigma": [0.0610108, 0.0309068], "exp": [1.13941, 3.64067]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.6061, 3.3469, -2.8389, -3.43, -116.4],
        "t": [1.0, 1.5, 1.95, 4.43, 29.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [6.4638, -17.8, 39.155, -47.82, 30.03, -7.079],
        "t": [0.5, 0.88, 1.28, 1.71, 2.18, 2.73]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.991, -39.27, 78.0, -69.23, -60.59, -158.25],
        "t": [0.444, 2.04, 2.55, 3.08, 7.63, 16.8]}

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
              "visco": "visco1",

              "ek": 502.1, "sigma": 0.547, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 0.934,

              "psi": [0.809284, 0.0881819, -0.0147911], "psi_d": [0, 1, 2],
              "fint": [0.002], "fint_t": [0],
              "chi": [1.14085, -0.11208, 0.0189958], "chi_d": [0, 1, 2],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.160e-9, "gam0": 0.098, "qd": 0.666e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""

    def test_Huber(self):
        """Table 7, pag 266"""
        st = ClBenzene(T=569.1, rhom=6.658)
        self.assertEqual(round(st.mu.muPas, 4), 101.6731)
        self.assertEqual(round(st.k.mWmK, 4), 84.6359)

    def test_Surface(self):
        """Table 10, pag 271"""
        self.assertEqual(round(ClBenzene(T=569.1, x=0.5).sigma, 7), 0.0044341)
