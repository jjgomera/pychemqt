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
from lib.mEoS import nC8


class iC8(MEoS):
    """Multiparameter equation of state for isooctane"""
    name = "isooctane"
    CASNumber = "540-84-1"
    formula = "(CH3)2CHCH2C(CH3)3"
    synonym = ""
    _refPropName = "IOCTANE"
    _coolPropName = ""
    rhoc = unidades.Density(242.1644624)
    Tc = unidades.Temperature(544)
    Pc = unidades.Pressure(2572.0, "kPa")
    M = 114.22852  # g/mol
    Tt = unidades.Temperature(165.77)
    Tb = unidades.Temperature(372.358)
    f_acent = 0.303
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 82

    CP1 = {"ao": 10.76,
           "ao_exp": [15.48, 34.42, 21.42],
           "exp": [775, 1900, 5100]}

    blackham = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for isooctane of Blackham "
                    "and Lemmon (2011).",
        "__doi__": {"autor": "Blackham, T.M., Lemmon, E.W.",
                    "title": "",
                    "ref": "to be published in Int. J. Thermophys., 2011.",
                    "doi": ""},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 1000000.0, "rhomax": 6.97,

        "nr1": [0.568901e-1, 0.196155e1, -0.281164e1, -0.815112, 0.326583],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.3, 0.75, 1.11, 0.55],

        "nr2": [-0.160893e1, -0.454734, 0.108306e1, -0.722876, -0.434052e-1],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.2, 3.7, 1.53, 2.1, 0.9],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*6,

        "nr3": [0.196648e1, -0.465082, -0.409398, 0.232131e-1],
        "d3": [1, 1, 3, 3, ],
        "t3": [0.88, 1.1, 2.75, 1.0],
        "alfa3": [0.75, 1.13, 0.87, 4.73],
        "beta3": [0.59, 1.45, 0.5, 10.52],
        "gamma3": [1.44, 0.68, 0.51, 0.8],
        "epsilon3": [0.66, 0.9, 0.54, 0.18]}

    eq = (blackham, )

    _surface = {
        "__doi__": {
            "autor": "Mulero, A., Cachadiña, I.",
            "title": "Recommended Correlations for the Surface Tension of "
                     "Several Fluids Included in the REFPROP Program",
            "ref": "J. Phys. Chem. Ref. Data 43(2) (2014) 023104",
            "doi": "10.1063/1.4878755"},
        "sigma": [0.04794], "exp": [1.209]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.7985, 8.1280, -7.3106, -3.9392, -1.6732],
        "t": [1, 1.5, 1.6, 4.0, 16.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [1.1535, 1.3709, 0.38804],
        "t": [0.286, 0.54, 3.3]}
    _vapor_Density = {
        "eq": 2,
        "n": [-2.5793, -6.4934, -18.631, -54.123, -123.58],
        "t": [0.366, 1.11, 3.0, 6.4, 14.0]}

    trnECS = {"__name__": "Huber (2018)",

              "__doi__": {
                  "autor": "Huber, M.L.",
                  "title": "Models for Viscosity, Thermal Conductivity, and "
                           "Surface Tension of Selected Pure Fluids as "
                           "Implemented in REFPROP v10.0",
                  "ref": "NISTIR 8209",
                  "doi": "10.6028/NIST.IR.8209"},

              "eq": "ecs",
              "ref": nC8,

              "ek": 635.7, "sigma": 0.588, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [1.09755, -0.0223075],
              "psi_d": [0, 1],
              "fint": [0.00115], "fint_t": [0],
              "chi": [0.827544, 0.0391177], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.256e-9, "gam0": 0.059, "qd": 0.771e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""
    def test_Huber(self):
        """Table 7, pag 266"""
        st = iC8(T=489.6, rhom=4.325)
        # self.assertEqual(round(st.mu.muPas, 5), 92.51572)
        # self.assertEqual(round(st.k.mWmK, 4), 59.0987)
        self.assertEqual(round(st.mu.muPas, 5), 92.51538)
        self.assertEqual(round(st.k.mWmK, 4), 59.0985)
