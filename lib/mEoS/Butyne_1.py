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


class Butyne_1(MEoS):
    """Multiparameter equation of state for 1-butyne"""
    name = "1-butyne"
    CASNumber = "107-00-6"
    formula = "CH≅CCH2CH3"
    synonym = ""
    _refPropName = "1BUTYNE"
    _coolPropName = ""
    rhoc = unidades.Density(251.520546)
    Tc = unidades.Temperature(432)
    Pc = unidades.Pressure(4141.6, "kPa")
    M = 54.09044  # g/mol
    Tt = unidades.Temperature(147.44)
    Tb = unidades.Temperature(281.22)
    f_acent = 0.28
    momentoDipolar = unidades.DipoleMoment(0, "Debye")
    id = 67

    CP1 = {"ao": 4,
           "an": [], "pow": [],
           "ao_exp": [2.5772, 9.5417, 11.125],
           "exp": [325.0, 991.0, 3349.0]}

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

        "Tmin": Tt, "Tmax": 532.0, "Pmax": 10000.0, "rhomax": 35.57,

        "nr1": [0.806625, -1.75423, -0.219937, 0.109675, 8.9653e-05],
        "d1": [1.0, 1.0, 2.0, 3.0, 8.0],
        "t1": [0.125, 1.125, 1.25, 0.25, 0.75],

        "nr2": [0.277043, 0.159607, -0.348927, -0.0570926, -0.00765105],
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
        "sigma": [0.0564795], "exp": [1.06959]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.788, 4.6419, -4.3658, -3.1948, -6.3609],
        "t": [1.0, 1.5, 1.84, 4.6, 20.5]}
    _liquid_Density = {
        "eq": 1,
        "n": [1.9537, 1.3293, -1.6707, 1.2883, 0.48718],
        "t": [0.314, 1.17, 2.0, 2.75, 12.0]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.08241, -6.39764, -15.8641, -41.4891, -61.7023, -148.468],
        "t": [0.366, 1.243, 2.93, 5.9, 10.8, 19.0]}

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

              "ek": 343.0, "sigma": 0.4847, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [0.98], "psi_d": [0],
              "fint": [0.0012], "fint_t": [0],
              "chi": [0.95], "chi_d": [0],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.199e-9, "gam0": 0.054, "qd": 0.588e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""

    def test_Huber(self):
        """Table 7, pag 266"""
        st = Butyne_1(T=388.8, rhom=9.61)
        # self.assertEqual(round(st.mu.muPas, 5), 86.57655)
        self.assertEqual(round(st.mu.muPas, 5), 86.57846)
        self.assertEqual(round(st.k.mWmK, 4), 81.2199)

    def test_Surface(self):
        """Table 10, pag 271"""
        self.assertEqual(round(Butyne_1(T=388.8, x=0.5).sigma, 7), 0.0048117)
