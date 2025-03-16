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


class Butadiene13(MEoS):
    """Multiparameter equation of state for 1,3-butadiene"""
    name = "1,3-butadiene"
    CASNumber = "106-99-0"
    formula = "CH2CHCHCH2"
    synonym = ""
    _refPropName = "13BUTADIENE"
    _coolPropName = ""
    rhoc = unidades.Density(245.0296932)
    Tc = unidades.Temperature(425.135)
    Pc = unidades.Pressure(4305.3, "kPa")
    M = 54.09044  # g/mol
    Tt = unidades.Temperature(164.25)
    Tb = unidades.Temperature(268.74)
    f_acent = 0.192
    momentoDipolar = unidades.DipoleMoment(0, "Debye")
    id = 28

    CP1 = {"ao": 4,
           "an": [], "pow": [],
           "ao_exp": [2.3797, 11.213, 8.1824],
           "exp": [308.0, 1210.0, 2631.0]}

    refprop = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for 1,3-butadiene (2022)",
        "__doi__": {
            "autor": "Huber, M.L., Lemmon, E.W., Bell, I.H., McLinden, M.O.",
            "title": "The NIST REFPROP Database for Highly Accurate Properties"
                     " of Industrially Importants Fluids",
            "ref": "Ind. Eng. Chem. Res. 61(42) (2022) 15449-15472",
            "doi": "10.1021/acs.iecr.2c01427"},

        "R": 8.3144598,
        "cp": CP1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 526.0, "Pmax": 10000.0, "rhomax": 35.57,

        "nr1": [0.84958959, -2.028566, 0.03757012, 0.06081274, 0.000204315],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.066043135, 0.42974178, 0.004085874, -0.25316, 0.02068683,
                -0.041815041, -0.02051974],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2.0, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    eq = (refprop, )

    _surface = {
        "__doi__": {"autor": "Huber, M.L.",
                    "title": "Models for Viscosity, Thermal Conductivity, and "
                             "Surface Tension of Selected Pure Fluids as "
                             "Implemented in REFPROP v10.0",
                    "ref": "NISTIR 8209",
                    "doi": "10.6028/NIST.IR.8209"},
        "sigma": [0.045947], "exp": [0.960983]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.2489, 2.8923, -4.1717, 4.3531, -5.5279],
        "t": [1.0, 1.5, 2.2, 3.0, 4.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [8.752, -23.981, 31.388, -20.966, 7.6897],
        "t": [0.495, 0.75, 1.0, 1.4, 1.8]}
    _vapor_Density = {
        "eq": 2,
        "n": [-2.9731, -6.7472, -17.647, -51.073, -114.93],
        "t": [0.392, 1.288, 3.31, 6.9, 15.0]}

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

              "ek": 337.6, "sigma": 0.4889, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [0.569829, 0.169932, -0.00650648], "psi_d": [0, 1, 2],
              "fint": [0.0012], "fint_t": [0],
              "chi": [0.95], "chi_d": [0],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.207e-9, "gam0": 0.057, "qd": 0.593e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""

    def test_Huber(self):
        """Table 7, pag 266"""
        st = Butadiene13(T=382.6, rhom=9.139)
        # self.assertEqual(round(st.mu.muPas, 5), 61.35637)
        self.assertEqual(round(st.mu.muPas, 5), 61.35652)
        self.assertEqual(round(st.k.mWmK, 4), 75.5409)

    def test_Surface(self):
        """Table 10, pag 271"""
        self.assertEqual(round(Butadiene13(T=382.6, x=0.5).sigma, 7), 0.005029)
