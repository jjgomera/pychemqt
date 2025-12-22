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


class DEA(MEoS):
    """Multiparameter equation of state for diethanolamine """
    name = "Diethanolamine"
    CASNumber = "111-42-2"
    formula = "CH2OH-CH2-NH-CH2-CH2OH"
    synonym = "2,2'-Iminodiethanol"
    _refPropName = "DEA"
    _coolPropName = ""
    rhoc = unidades.Density(346.94748)
    Tc = unidades.Temperature(736.5)
    Pc = unidades.Pressure(4950.75, "kPa")
    M = 105.1356  # g/mol
    Tt = unidades.Temperature(301.1)
    Tb = unidades.Temperature(541.54)
    f_acent = 1.013
    momentoDipolar = unidades.DipoleMoment(0.85, "Debye")
    id = 428

    CP1 = {"ao": 4,
           "an": [], "pow": [],
           "ao_exp": [4.25, 37.7],
           "exp": [96.0, 1165.0]}

    refprop = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for diethanolamine (2022)",
        "__doi__": {
            "autor": "Huber, M.L., Lemmon, E.W., Bell, I.H., McLinden, M.O.",
            "title": "The NIST REFPROP Database for Highly Accurate Properties"
                     " of Industrially Importants Fluids",
            "ref": "Ind. Eng. Chem. Res. 61(42) (2022) 15449-15472",
            "doi": "10.1021/acs.iecr.2c01427"},

        "R": 8.3144598,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 840.0, "Pmax": 5000.0, "rhomax": 35.57,

        "nr1": [0.066088158, 6.1059245, -7.0526968, -0.29739545, 0.11592105],
        "d1": [4.0, 1.0, 1.0, 2.0, 3.0],
        "t1": [1.0, 0.507, 0.907, 1.22, 0.649],

        "nr2": [-1.8616953, -0.97392153, .14690655, -0.63284478, -0.037820123],
        "d2": [1.0, 3.0, 2.0, 2.0, 7.0],
        "t2": [2.14, 2.89, 1.54, 3.34, 0.998],
        "c2": [2.0, 2.0, 1.0, 2.0, 1.0],
        "gamma2": [1]*5,

        "nr3": [2.774726, -1.0230468, -0.19552536, -0.14997584, -0.23421918],
        "d3": [1.0, 1.0, 3.0, 2.0, 2.0],
        "t3": [0.92, 1.16, 1.43, 1.24, 0.801],
        "alfa3": [0.8971, 1.499, 1.681, 1.661, 1.245],
        "beta3": [0.9691, 1.518, 1.328, 1.23, 1.112],
        "gamma3": [1.216, 0.6775, 0.7815, 0.8796, 1.357],
        "epsilon3": [0.6694, 0.6466, 0.6669, 1.189, 1.248]}

    eq = (refprop, )

    _surface = {
        "__doi__": {"autor": "Huber, M.L.",
                    "title": "Models for Viscosity, Thermal Conductivity, and "
                             "Surface Tension of Selected Pure Fluids as "
                             "Implemented in REFPROP v10.0",
                    "ref": "NISTIR 8209",
                    "doi": "10.6028/NIST.IR.8209"},
        "sigma": [0.0859443], "exp": [1.15945]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-14.957, 10.22, -5.9784, -4.968, -3.13],
        "t": [1.0, 1.33, 1.8, 3.0, 10.3]}
    _liquid_Density = {
        "eq": 1,
        "n": [-0.1814, 2.3567, 1.81, -3.03, 1.874],
        "t": [0.02, 0.26, 1.43, 2.0, 2.6]}
    _vapor_Density = {
        "eq": 2,
        "n": [0.108, -4.888, -12.7, -47.505, -93.42, -221.0],
        "t": [0.09, 0.395, 1.467, 3.7, 8.05, 16.2]}

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

              "ek": 584.85, "sigma": 0.5434, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 0.74,

              "psi": [0.996593, 0.0399708], "psi_d": [0, 1],
              "fint": [0.00132], "fint_t": [0],
              "chi": [1.47408, -0.123082], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.185e-9, "gam0": 0.068, "qd": 0.662e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""

    def test_Huber(self):
        """Table 7, pag 266"""
        st = DEA(T=662.8, rhom=7.262)
        # self.assertEqual(round(st.mu.muPas, 4), 228.4693)
        self.assertEqual(round(st.mu.muPas, 4), 228.4695)
        self.assertEqual(round(st.k.mWmK, 4), 142.1575)

    def test_Surface(self):
        """Table 10, pag 271"""
        self.assertEqual(round(DEA(T=662.8, x=0.5).sigma, 7), 0.0059581)
