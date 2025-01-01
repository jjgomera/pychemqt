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


class RE347mcc(MEoS):
    """Multiparameter equation of state for RE347mcc"""
    name = "methyl-heptafluoropropyl-ether"
    CASNumber = "375-03-1"
    formula = "CF3CF2CF2OCH3"
    synonym = "HFE-7000"
    _refPropName = "RE347MCC"
    _coolPropName = ""
    rhoc = unidades.Density(528.144783936)
    Tc = unidades.Temperature(437.7)
    Pc = unidades.Pressure(2476.2, "kPa")
    M = 200.0548424  # g/mol
    Tt = unidades.Temperature(250)
    Tb = unidades.Temperature(307.349)
    f_acent = 0.403
    momentoDipolar = unidades.DipoleMoment(3.13, "Debye")

    CP1 = {"ao": 17.916,
           "ao_exp": [0.6505, 0.3794, 21.292], "exp": [21, 7754, 1562]}

    zhou = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for RE347mcc of Zhou (2012)",
        "__doi__": {"autor": "Zhou, Y., Lemmon, E.W., Mahmoud, A.M.",
                    "title": "Equations of state for RE245cb2, RE347mcc, "
                             "RE245fa2 and R1216",
                    "ref": "Preliminary equation",
                    "doi": ""},

        "R": 8.3144598,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 20000.0, "rhomax": 9,

        "nr1": [0.07342, 1.9394, -2.8353, -0.3876, 0.1428],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1.0, 0.185, 0.842, 1.0, 0.5],

        "nr2": [-1.979, -2.0455, 0.3085, -2.166, -0.04225],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.74, 2.74, 0.87, 2.77, 1.24],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [3.0317, -1.0685, -0.3598, -0.4525, 0.9488, 0.3184],
        "d3": [1, 1, 3, 3, 2, 3],
        "t3": [2.3, 1.74, 2.22, 2.1, 1.85, 0.5],
        "alfa3": [1.205, 1.19, 0.94, 1.64, 0.92, 1.71],
        "beta3": [0.53, 2.57, 1.6, 2.56, 1.21, 8.9],
        "gamma3": [0.9, 0.69, 0.87, 1.22, 0.71, 1.11],
        "epsilon3": [0.655, 0.9, 0.655, 0.34, 0.745, 0.28]}

    eq = (zhou, )

    _surface = {
        "__doi__": {
            "autor": "Mulero, A., Cachadiña, I.",
            "title": "Recommended Correlations for the Surface Tension of "
                     "Several Fluids Included in the REFPROP Program",
            "ref": "J. Phys. Chem. Ref. Data 43(2) (2014) 023104",
            "doi": "10.1063/1.4878755"},
        "sigma": [0.05031], "exp": [1.232]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.911, 1.4904, -3.0464, -4.9639, -7.7423],
        "t": [1.0, 1.5, 2.7, 4.8, 13.5]}
    _liquid_Density = {
        "eq": 1,
        "n": [3.1002, -3.1869, 8.0538, -7.5947, 2.8275],
        "t": [0.395, 0.75, 1.15, 1.5, 2.15]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.2144, -7.0853, -23.82, -69.536, -182.42,],
        "t": [0.4, 1.21, 3.19, 6.65, 14.0, 27.85]}

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

              "ek": 347.6, "sigma": 0.5853, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [2.57345, -1.73973, 0.659016, -0.0824925],
              "psi_d": [0, 1, 2, 3]}

    _viscosity = (trnECS, )

    thermo0 = {"__name__": "Huber (2018)",

               "__doi__": {
                   "autor": "Huber, M.L.",
                   "title": "Models for Viscosity, Thermal Conductivity, and "
                            "Surface Tension of Selected Pure Fluids as "
                            "Implemented in REFPROP v10.0",
                   "ref": "NISTIR 8209",
                   "doi": "10.6028/NIST.IR.8209"},

               "eq": 1,

               "Toref": 437.7, "koref": 1,
               "no": [-0.0239098, 0.0960335, -0.060505, 0.012299],
               "to": [1, 2, 3, 4],

               "Tref_res": 437.7, "rhoref_res": 528.144783936,
               "kref_res": 1,
               "nr": [-0.00842403, 0.0545889, -0.0530301, 0.0201447,
                      -0.0025046, 0.00931228, -0.0367016, 0.0392477,
                      -0.0155674, 0.00220816],
               "tr": [0, 0, 0, 0, 0, -1, -1, -1, -1, -1],
               "dr": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
               "Xio": 0.231e-9, "gam0": .058, "qd": 0.5553e-9, "Tcref": 1.5*Tc}

    _thermal = (thermo0, )


class Test(TestCase):
    """Testing"""
    def test_Huber(self):
        """Table 7, pag 266"""
        self.assertEqual(round(
            # RE347mcc(T=393.9, rhom=5.48).mu.muPas, 4), 150.0694)
            RE347mcc(T=393.9, rhom=5.48).mu.muPas, 4), 149.7817)

        # Table 9, pag 271
        self.assertEqual(round(
            # RE347mcc(T=393.93, rhom=5.4797).k.mWmK, 4), 49.3925)
            RE347mcc(T=393.93, rhom=5.4797).k.mWmK, 4), 49.4395)
