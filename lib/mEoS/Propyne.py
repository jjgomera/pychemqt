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


class Propyne(MEoS):
    """Multiparamenter equation of state for propyne"""
    name = "Propyne"
    CASNumber = "74-99-7"
    formula = "CH3-C≡CH"
    synonym = ""
    _refPropName = "PROPYNE"
    _coolPropName = "Propyne"
    rhoc = unidades.Density(244.898798)
    Tc = unidades.Temperature(402.38)
    Pc = unidades.Pressure(5626.0, "kPa")
    M = 40.06  # g/mol
    Tt = unidades.Temperature(170.5)
    Tb = unidades.Temperature(248.0)
    f_acent = 0.204
    momentoDipolar = unidades.DipoleMoment(0.781, "Debye")
    id = 66

    CP1 = {"ao": 0.342418/8.3143*40.06,
           "an": [0.484403e-2/8.3143*40.06, -0.347414e-5/8.3143*40.06,
                  0.144887e-8/8.3143*40.06, -0.26815e-12/8.3143*40.06],
           "pow": [1, 2, 3, 4]}

    polt = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for propyne of Polt (1992)",
        "__doi__": {"autor": "Polt, A., Platzer, B., Maurer, G.",
                    "title": "Parameter der thermischen Zustandsgleichung von "
                             "Bender fuer 14 mehratomige reine Stoffe",
                    "ref": "Chem. Technik 22(1992)6 , 216/224",
                    "doi": ""},

        "R": 8.3143,
        "cp": CP1,
        "ref": "IIR",

        "Tmin": 273, "Tmax": 474.0, "Pmax": 32000.0, "rhomax": 16.28,

        "nr1": [0.102590136933e1, -0.220786016506e1, 0.107889905204e1,
                -0.986950667682, 0.459528109357e1, -0.886063623532e1,
                0.556346955561e1, -0.157450028544e1, -0.159068753573,
                0.235738270184, 0.440755494599, 0.196126150614, -0.36775965033,
                0.792931851008e-2, 0.247509085735e-2, 0.832903610194e-2],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [-0.102590136933e1, 0.220786016506e1, -0.107889905204e1,
                -0.382188466986e1, 0.830345065619e1, -0.448323072603e1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [1.65533788]*6}

    eq = (polt, )

    _surface = {"sigma": [0.05801], "exp": [1.205]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.69162e1, 0.10904e1, -0.74791, 0.75926e1, -0.25926e2],
        "t": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.22754, 0.33173e1, -0.18041e1, 0.22440e1, -0.35823],
        "t": [0.1, 0.53, 1.0, 2.0, 3.0]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.17504, -4.6021, -89.211, 180.02, -243.99, 160.35],
        "t": [0.1, 0.56, 2.5, 3.0, 4.0, 5.0]}

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

              "ek": 246.85, "sigma": 0.478, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [0.98], "psi_d": [0],
              "fint": [0.0012], "fint_t": [0],
              "chi": [0.93], "chi_d": [0],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.186e-9, "gam0": 0.058, "qd": 0.535e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""
    def test_Huber(self):
        """Table 7, pag 266"""
        st = Propyne(T=362.1, rhom=12.571)
        # self.assertEqual(round(st.mu.muPas, 5), 80.62469)
        # self.assertEqual(round(st.k.mWmK, 4), 91.2068)
        self.assertEqual(round(st.mu.muPas, 5), 80.62182)
        self.assertEqual(round(st.k.mWmK, 4), 91.2093)
