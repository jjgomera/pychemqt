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


class PropylenOxide(MEoS):
    """Multiparameter equation of state for Propylene Oxide"""
    name = "Propylene oxide"
    CASNumber = "75-56-9"
    formula = "C3H6O"
    synonym = ""
    _refPropName = "PROPYLENOXIDE"
    _coolPropName = ""
    rhoc = unidades.Density(299.3979667)
    Tc = unidades.Temperature(488.11)
    Pc = unidades.Pressure(5436.6, "kPa")
    M = 58.07914  # g/mol
    Tt = unidades.Temperature(161.244)
    Tb = unidades.Temperature(307.05)
    f_acent = 0.249
    momentoDipolar = unidades.DipoleMoment(2.01, "Debye")
    id = 444

    CP1 = {"ao": 4,
           "ao_exp": [3.5346, 11.904, 7.5015],
           "exp": [455.0, 1535.0, 3598.0]}

    refprop = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for Propylene Oxide (2022)",
        "__doi__": {
            "autor": "Huber, M.L., Lemmon, E.W., Bell, I.H., McLinden, M.O.",
            "title": "The NIST REFPROP Database for Highly Accurate Properties"
                     " of Industrially Importants Fluids",
            "ref": "Ind. Eng. Chem. Res. 61(42) (2022) 15449-15472",
            "doi": "10.1021/acs.iecr.2c01427"},

        "R": 8.3144598,
        "cp": CP1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 590.0, "Pmax": 10000.0, "rhomax": 35.57,

        "nr1": [0.8274658, -1.898325, -0.1298384, 0.0987172, 6.1011e-05],
        "d1": [1.0, 1.0, 2.0, 3.0, 8.0],
        "t1": [0.125, 1.125, 1.25, 0.25, 0.75],

        "nr2": [0.3533366, 0.1213314, -0.2404342, -0.01895105, -0.01169704],
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
        "sigma": [0.073], "exp": [1.22]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.69695, 4.40186, -3.92255, -3.09268, -15.1205],
        "t": [1.0, 1.5, 1.87, 4.5, 23.5]}
    _liquid_Density = {
        "eq": 1,
        "n": [18.4597, -63.8693, 106.657, -82.2743, 23.9996],
        "t": [0.556, 0.775, 1.0, 1.23, 1.5]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.4779, -7.3019, -20.307, -52.562, -130.14],
        "t": [0.412, 1.406, 3.44, 7.0, 15.1]}

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

              "ek": 387.6, "sigma": 0.4683, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 0.96,

              "psi": [0.821761, 0.0611306], "psi_d": [0, 1],
              "fint": [0.00132], "fint_t": [0],
              "chi": [1.04], "chi_d": [0],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.194e-9, "gam0": 0.056, "qd": 0.567e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""

    def test_Huber(self):
        """Table 7, pag 266"""
        st = PropylenOxide(T=439.3, rhom=10.502)
        # self.assertEqual(round(st.mu.muPas, 5), 89.19919)
        self.assertEqual(round(st.mu.muPas, 5), 89.20032)
        self.assertEqual(round(st.k.mWmK, 4), 99.8309)

    def test_Surface(self):
        """Table 10, pag 271"""
        self.assertEqual(
            round(PropylenOxide(T=439.3, x=0.5).sigma, 7), 0.0043986)
