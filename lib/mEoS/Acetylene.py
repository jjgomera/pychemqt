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


class Acetylene(MEoS):
    """Multiparameter equation of state for """
    name = "Acetylene"
    CASNumber = "74-86-2"
    formula = "CH≅CH"
    synonym = ""
    _refPropName = ""
    _coolPropName = ""
    rhoc = unidades.Density(229.9091824)
    Tc = unidades.Temperature(308.3)
    Pc = unidades.Pressure(5988.2, "kPa")
    M = 26.03728  # g/mol
    Tt = unidades.Temperature(191.75)
    Tb = unidades.Temperature(189)
    f_acent = 0.178
    momentoDipolar = unidades.DipoleMoment(0, "Debye")
    id = 65

    CP1 = {"ao": 3.5,
           "an": [], "pow": [],
           "ao_exp": [1.8374, 3.4338, 2.6089],
           "exp": [610.0, 1428.0, 6857.0]}

    refprop = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for Acetylene (2022)",
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

        "nr1": [0.8157856, -1.85265, -0.115457, 0.0938171, 6.405e-05],
        "d1": [1.0, 1.0, 2.0, 3.0, 8.0],
        "t1": [0.125, 1.125, 1.25, 0.25, 0.75],

        "nr2": [0.2031037, 0.1417312, -0.1641216, 0.01495196, -0.014536],
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
        "sigma": [0.0615167], "exp": [1.19797]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.2279, 3.8389, -6.5192, 7.0358, -6.3104],
        "t": [1.0, 1.5, 2.1, 2.8, 3.6]}
    _liquid_Density = {
        "eq": 1,
        "n": [1.4146, 1.4232, -0.22333, 0.14913, 2.8826],
        "t": [0.267, 0.84, 1.74, 2.9, 10.0]}
    _vapor_Density = {
        "eq": 2,
        "n": [-1.93926, -4.21026, -10.0508, -28.8004, -82.9808],
        "t": [0.3, 0.883, 2.2, 4.75, 9.8]}

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

              "ek": 244.8, "sigma": 0.3914, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 0.97,

              "psi": [0.98], "psi_d": [0],
              "fint": [6.46532e-4, 9.98935e-7, 1.22755e-10],
              "fint_t": [0, 1, 2],
              "chi": [0.9], "chi_d": [0],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.166e-9, "gam0": 0.056, "qd": 0.47e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""

    def test_Huber(self):
        """Table 7, pag 266"""
        st = Acetylene(T=277.5, rhom=17.772)
        self.assertEqual(round(st.mu.muPas, 5), 70.87286)
        self.assertEqual(round(st.k.mWmK, 4), 105.9641)

    def test_Surface(self):
        """Table 10, pag 271"""
        self.assertEqual(round(Acetylene(T=277.5, x=0.5).sigma, 7), 0.0038951)
