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


class Pentene_1(MEoS):
    """Multiparameter equation of state for 1-Pentene"""
    name = "1-Pentene"
    CASNumber = "109-67-1"
    formula = "CH2=CHCH2CH2CH3"
    synonym = ""
    _refPropName = "1PENTENE"
    _coolPropName = ""
    rhoc = unidades.Density(241.958505)
    Tc = unidades.Temperature(465.74)
    Pc = unidades.Pressure(3598, "kPa")
    M = 70.1329  # g/mol
    Tt = unidades.Temperature(107.797)
    Tb = unidades.Temperature(303.11)
    f_acent = 0.233
    momentoDipolar = unidades.DipoleMoment(0, "Debye")
    id = 29

    CP1 = {"ao": 8,
           "an": [], "pow": [],
           "ao_exp": [9.4073, 16.83, 8.2294],
           "exp": [934.0, 2082.0, 5594.0]}

    refprop = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for 1-Pentene (2022)",
        "__doi__": {
            "autor": "Huber, M.L., Lemmon, E.W., Bell, I.H., McLinden, M.O.",
            "title": "The NIST REFPROP Database for Highly Accurate Properties"
                     " of Industrially Importants Fluids",
            "ref": "Ind. Eng. Chem. Res. 61(42) (2022) 15449-15472",
            "doi": "10.1021/acs.iecr.2c01427"},

        "R": 8.3144598,
        "cp": CP1,
        "ref": "IIR",

        "Tmin": 150, "Tmax": 566.0, "Pmax": 10000.0, "rhomax": 135.57,

        "nr1": [0.772189126, -1.83042132, -0.192259446, 0.115251009, 8.839e-5],
        "d1": [1.0, 1.0, 2.0, 3.0, 8.0,],
        "t1": [0.125, 1.125, 1.25, 0.25, 0.75,],

        "nr2": [.432066965, .190314262, -.24306775, -.030971971, -.011097285],
        "d2": [2.0, 3.0, 1.0, 4.0, 3.0],
        "t2": [0.625, 2.0, 4.125, 4.125, 17.0],
        "c2": [1.0, 1.0, 2.0, 2.0, 3.0],
        "gamma2": [1]*5}

    eq = (refprop, )

    _surface = {
        "__doi__": {"autor": "Huber, M.L.",
                    "title": "Models for Viscosity, Thermal Conductivity, and "
                             "Surface Tension of Selected Pure Fluids as "
                             "Implemented in REFPROP v10.0",
                    "ref": "NISTIR 8209",
                    "doi": "10.6028/NIST.IR.8209"},
        "sigma": [0.050798], "exp": [1.16356]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.48233, 2.89916, -2.8279, -10.6523, -8.26335, 15.6616],
        "t": [1.0, 1.5, 2.15, 5.9, 8.5, 7.15]}
    _liquid_Density = {
        "eq": 1,
        "n": [7.6105, -18.34, 29.764, -23.321, 7.1874, 0.54244],
        "t": [0.546, 0.876, 1.23, 1.6, 2.1, 13.5]}
    _vapor_Density = {
        "eq": 2,
        "n": [-4.4293, -8.8291, -23.467, -50.892, -113.57, -224.82],
        "t": [0.495, 1.743, 4.0, 7.65, 15.5, 30.7]}

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

              "ek": 369.8, "sigma": 0.5354, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 0.96,

              "psi": [0.97], "psi_d": [0],
              "fint": [0.00077636, 1.5844e-6, -1.19304e-9],
              "fint_t": [0, 1, 2],
              "chi": [1.07723, -0.00850457], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.223e-9, "gam0": 0.058, "qd": 0.652e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""

    def test_Huber(self):
        """Table 7, pag 266"""
        st = Pentene_1(T=419.2, rhom=6.898)
        # self.assertEqual(round(st.mu.muPas, 5), 75.38902)
        self.assertEqual(round(st.mu.muPas, 5), 75.38935)
        self.assertEqual(round(st.k.mWmK, 4), 75.3082)

    def test_Surface(self):
        """Table 10, pag 271"""
        self.assertEqual(round(Pentene_1(T=419.2, x=0.5).sigma, 7), 0.0034827)
