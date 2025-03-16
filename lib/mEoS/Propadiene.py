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


class Propadiene(MEoS):
    """Multiparameter equation of state for Propadiene"""
    name = "Propadiene"
    CASNumber = "463-49-0"
    formula = "CH2=C=CH2"
    synonym = ""
    _refPropName = "PROPADIENE"
    _coolPropName = ""
    rhoc = unidades.Density(236.376774)
    Tc = unidades.Temperature(398)
    Pc = unidades.Pressure(5215.6, "kPa")
    M = 40.06386  # g/mol
    Tt = unidades.Temperature(136.65)
    Tb = unidades.Temperature(238.65)
    f_acent = 0.115
    momentoDipolar = unidades.DipoleMoment(0, "Debye")
    id = 57

    CP1 = {"ao": 4,
           "ao_exp": [2.126, 7.4934, 5.224],
           "exp": [512.0, 1442.0, 3896.0]}

    refprop = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for Propadiene (2022)",
        "__doi__": {
            "autor": "Huber, M.L., Lemmon, E.W., Bell, I.H., McLinden, M.O.",
            "title": "The NIST REFPROP Database for Highly Accurate Properties"
                     " of Industrially Importants Fluids",
            "ref": "Ind. Eng. Chem. Res. 61(42) (2022) 15449-15472",
            "doi": "10.1021/acs.iecr.2c01427"},

        "R": 8.3144598,
        "cp": CP1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 10000.0, "rhomax": 55.57,

        "nr1": [0.7231448, -1.790058, -0.06836828, 0.07947672, 4.0778e-05],
        "d1": [1.0, 1.0, 2.0, 3.0, 8.0],
        "t1": [0.125, 1.125, 1.25, 0.25, 0.75],

        "nr2": [0.1760558, 0.1443484, -0.1494723, 0.008248376, -0.009386559],
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
        "sigma": [0.056], "exp": [1.205]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-6.771, 3.0853, -3.3137, 1.2095, -2.7511],
        "t": [1.0, 1.5, 2.0, 2.8, 4.3]}
    _liquid_Density = {
        "eq": 1,
        "n": [2.5282, -1.5802, 5.4328, -5.5434, 2.351],
        "t": [0.384, 0.9, 1.5, 2.1, 3.0]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.333, -7.2503, -47.109, -18.886, -108.51],
        "t": [0.415, 1.5, 7.5, 3.68, 16.1]}

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

              "ek": 316, "sigma": 0.4477, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [0.99], "psi_d": [0],
              "fint": [0.0012], "fint_t": [0],
              "chi": [0.95], "chi_d": [0],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.195e-9, "gam0": 0.055, "qd": 0.541e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""

    def test_Huber(self):
        """Table 7, pag 266"""
        st = Propadiene(T=358.2, rhom=11.972)
        # self.assertEqual(round(st.mu.muPas, 5), 78.02296)
        self.assertEqual(round(st.mu.muPas, 5), 78.02271)
        self.assertEqual(round(st.k.mWmK, 4), 89.6455)

    def test_Surface(self):
        """Table 10, pag 271"""
        self.assertEqual(round(Propadiene(T=358.2, x=0.5).sigma, 7), 0.0034929)
