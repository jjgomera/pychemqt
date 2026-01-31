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


class MEA(MEoS):
    """Multiparameter equation of state for Monoethanolamine"""
    name = "Monoethanolamine"
    CASNumber = "141-43-5"
    formula = "CH2OH-CH2-NH2"
    synonym = ""
    _refPropName = "MEA"
    _coolPropName = ""
    rhoc = unidades.Density(329.237909)
    Tc = unidades.Temperature(671.4)
    Pc = unidades.Pressure(8125, "kPa")
    M = 61.0831  # g/mol
    Tt = unidades.Temperature(283.7)
    Tb = unidades.Temperature(443.564)
    f_acent = 0.573
    momentoDipolar = unidades.DipoleMoment(0.7765, "Debye")
    id = 250

    CP1 = {"ao": 4,
           "an": [], "pow": [],
           "ao_exp": [13.7, 11.1],
           "exp": [970.0, 3380.0]}

    herrig = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for Monoethanolamine (2018)",
        "__doi__": {
            "autor": "Herrig, S.",
            "title": "New Helmholtz-Energy Equations of State for Pure "
                     "Fluids and CCS-Relevant Mixtures",
            "ref": "Ph.D. thesis. Bochum: Ruhr-Universität Bochum, 2018.",
            "doi": ""},

        "R": 8.3144598,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 775.0, "Pmax": 10000.0, "rhomax": 35.57,

        "nr1": [0.034371657, 2.804815, -3.5328022, -0.26052106, 0.073728099],
        "d1": [4.0, 1.0, 1.0, 2.0, 3.0],
        "t1": [1.0, 0.53, 1.146, 0.95, 0.35],

        "nr2": [-0.9232864, -0.15243636, .44837938, -0.17517565, -0.012936362],
        "d2": [1.0, 3.0, 2.0, 2.0, 7.0],
        "t2": [1.47, 2.8, 0.9, 3.0, 0.83],
        "c2": [2.0, 2.0, 1.0, 2.0, 1.0],
        "gamma2": [1]*5,

        "nr3": [1.0823719, -0.56755523, -0.38808402, -6.7388446],
        "d3": [1.0, 1.0, 3.0, 3.0],
        "t3": [1.03, 0.76, 0.7, 1.04],
        "alfa3": [0.71, 1.16, 0.733, 4.08],
        "beta3": [1.82, 1.5, 1.74, 57.0],
        "gamma3": [1.04, 1.04, 1.04, 1.37],
        "epsilon3": [0.84, 0.77, 0.6, 0.59]}

    eq = (herrig, )

    _surface = {
        "__doi__": {"autor": "Huber, M.L.",
                    "title": "Models for Viscosity, Thermal Conductivity, and "
                             "Surface Tension of Selected Pure Fluids as "
                             "Implemented in REFPROP v10.0",
                    "ref": "NISTIR 8209",
                    "doi": "10.6028/NIST.IR.8209"},
        "sigma": [0.0776613], "exp": [0.801525]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-9.5739, 5.593, -5.15, -5.047, -4.69],
        "t": [1.0, 1.5, 1.9, 3.7, 13.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [-0.103, 2.116, 3.495, -4.452, 1.795, 0.741],
        "t": [0.06, 0.3, 1.5, 2.0, 2.9, 19.0]}
    _vapor_Density = {
        "eq": 2,
        "n": [0.1, -2.907, -7.405, -25.86, -70.87, -210.8],
        "t": [0.07, 0.34, 1.1, 3.0, 6.34, 14.6]}

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

              "ek": 533.15, "sigma": 0.4614, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 0.88,

              "psi": [1.18676, -0.260695, 0.0789293], "psi_d": [0, 1, 2],
              "fint": [0.00132], "fint_t": [0],
              "chi": [1.61924, -0.210496], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.173e-9, "gam0": 0.065, "qd": 0.559e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""

    def test_Huber(self):
        """Table 7, pag 266"""
        st = MEA(T=604.3, rhom=11.24)
        self.assertEqual(round(st.mu.muPas, 4), 132.2859)
        self.assertEqual(round(st.k.mWmK, 4), 167.1428)

    def test_Surface(self):
        """Table 10, pag 271"""
        self.assertEqual(round(MEA(T=604.3, x=0.5).sigma, 7), 0.0122595)
