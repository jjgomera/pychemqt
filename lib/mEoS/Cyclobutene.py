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


class Cyclobutene(MEoS):
    """Multiparameter equation of state for cyclobutene """
    name = "Cyclobutene"
    CASNumber = "822-35-5"
    formula = ""
    synonym = ""
    _refPropName = "CYCLOBUTENE"
    _coolPropName = ""
    rhoc = unidades.Density(279.6475748)
    Tc = unidades.Temperature(448)
    Pc = unidades.Pressure(5149.5, "kPa")
    M = 54.09044  # g/mol
    Tt = unidades.Temperature(150)
    Tb = unidades.Temperature(275.73)
    f_acent = 0.163
    momentoDipolar = unidades.DipoleMoment(0, "Debye")
    # id =

    CP1 = {"ao": 4,
           "an": [], "pow": [],
           "ao_exp": [0.29258, 8.8201, 13.254],
           "exp": [370.0, 1094.0, 2446.0]}

    refprop = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for cyclobutene (2022)",
        "__doi__": {
            "autor": "Huber, M.L., Lemmon, E.W., Bell, I.H., McLinden, M.O.",
            "title": "The NIST REFPROP Database for Highly Accurate Properties"
                     " of Industrially Importants Fluids",
            "ref": "Ind. Eng. Chem. Res. 61(42) (2022) 15449-15472",
            "doi": "10.1021/acs.iecr.2c01427"},

        "R": 8.3144598,
        "cp": CP1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 548.0, "Pmax": 10000.0, "rhomax": 35.57,

        "nr1": [0.80542901, -1.759523, -0.1680518, 0.09838102, 0.000108512],
        "d1": [1.0, 1.0, 2.0, 3.0, 8.0],
        "t1": [0.125, 1.125, 1.25, 0.25, 0.75],

        "nr2": [0.15845628, 0.1815905, -0.2063626, -0.008469521, -0.01158867],
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
        "sigma": [0.0651302], "exp": [1.23574]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.15637, 3.67861, -4.41, 2.67868, -3.44738],
        "t": [1.0, 1.5, 2.0, 2.8, 4.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [1.9824, -0.1502, 2.2258, -3.2595, 1.9792],
        "t": [0.336, 0.911, 1.5, 2.15, 2.9]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.11265, -6.72153, -17.315, -45.5201, -99.9268],
        "t": [0.391, 1.364, 3.35, 7.0, 15.0]}

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

              "ek": 355.8, "sigma": 0.4679, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [1], "psi_d": [0],
              "fint": [0.0012], "fint_t": [0],
              "chi": [0.95], "chi_d": [0],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.2e-9, "gam0": 0.056, "qd": 0.567e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""

    def test_Huber(self):
        """Table 7, pag 266"""
        st = Cyclobutene(T=403.2, rhom=10.301)
        # self.assertEqual(round(st.mu.muPas, 5), 89.63551)
        self.assertEqual(round(st.mu.muPas, 5), 89.63524)
        self.assertEqual(round(st.k.mWmK, 4), 79.4797)

    def test_Surface(self):
        """Table 10, pag 271"""
        self.assertEqual(round(Cyclobutene(T=403.2, x=0.5).sigma, 7), 0.0037848)
