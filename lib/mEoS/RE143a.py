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
from lib.mEoS import R134a


class RE143a(MEoS):
    """Multiparameter equation of state for RE143a"""
    name = "methyl trifluoromethyl ether"
    CASNumber = "421-14-7"
    formula = "CH3-O-CF3"
    synonym = "HFE-143a"
    _refPropName = "RE143A"
    _coolPropName = ""
    rhoc = unidades.Density(465)
    Tc = unidades.Temperature(377.921)
    Pc = unidades.Pressure(3635., "kPa")
    M = 100.0398  # g/mol
    Tt = unidades.Temperature(240)
    Tb = unidades.Temperature(249.572)
    f_acent = 0.289
    momentoDipolar = unidades.DipoleMoment(2.32, "Debye")
    # id = 1817

    f = 8.314472
    CP1 = {"ao": 20.37/f,
           "an": [0.2918/f, -1.950e-4/f, 4.650e-8/f], "pow": [1, 2, 3]}

    akasaka = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for RE143a of Akasaka (2012)",
        "__doi__": {"autor": "Akasaka, R., Kayukawa, Y.",
                    "title": "A fundamental equation of state for "
                             "trifluoromethyl methyl ether (HFE-143m) and its "
                             "application to refrigeration cycle analysis",
                    "ref": "Int. J. Refrig., 35(4) (2012) 1003-1013",
                    "doi": "10.1016/j.ijrefrig.2012.01.003"},

        "R": 8.314472,
        "cp": CP1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 420.0, "Pmax": 7200.0, "rhomax": 12.62,

        "nr1": [0.77715884e1, -0.87042750e1, -0.28095049, 0.14540153,
                0.92291277e-2],
        "d1": [1, 1, 1, 2, 5],
        "t1": [0.682, 0.851, 1.84, 1.87, 0.353],

        "nr2": [-0.2141651, 0.99475155e-1, 0.23247135e-1, -0.12873573e-1,
                -0.57366549e-1, 0.3650465, -0.25433763, -0.90896436e-1,
                0.83503619e-1, 0.15477603e-1, -0.16641941e-1, 0.52410163e-2],
        "d2": [1, 3, 5, 7, 1, 2, 2, 3, 4, 2, 3, 5],
        "t2": [3.92, 1.14, 0.104, 1.19, 6.58, 6.73, 7.99, 7.31, 7.45, 16.5,
               24.8, 10.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3],
        "gamma2": [1]*12}

    eq = (akasaka, )

    _surface = {
        "__doi__": {"autor": "Huber, M.L.",
                    "title": "Models for Viscosity, Thermal Conductivity, and "
                             "Surface Tension of Selected Pure Fluids as "
                             "Implemented in REFPROP v10.0",
                    "ref": "NISTIR 8209",
                    "doi": "10.6028/NIST.IR.8209"},
        "sigma": [0.0371], "exp": [0.98412]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.44314, 1.69164, -2.27778, -4.094],
        "t": [1, 1.5, 2.5, 5]}
    _liquid_Density = {
        "eq": 1,
        "n": [1.20552, 1.33568, 0.0981486, 0.248917],
        "t": [0.33, 0.5, 1.5, 2.5]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.02576, -6.97239, -20.2601, -53.4441],
        "t": [0.38, 1.24, 3.2, 6.9]}

    trnECS = {"__name__": "Huber (2018)",

              "__doi__": {
                  "autor": "Huber, M.L.",
                  "title": "Models for Viscosity, Thermal Conductivity, and "
                           "Surface Tension of Selected Pure Fluids as "
                           "Implemented in REFPROP v10.0",
                  "ref": "NISTIR 8209",
                  "doi": "10.6028/NIST.IR.8209"},

              "eq": "ecs",
              "ref": R134a,

              "ek": 300.104, "sigma": 0.4847, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [1.008], "psi_d": [0],
              "fint": [0.001129], "fint_t": [0],
              "chi": [0.975], "chi_d": [0],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.198e-9, "gam0": 0.054, "qd": 0.588e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Test class"""

    def test_Huber(self):
        """Table 7, pag 266"""
        st = RE143a(T=340.1, rhom=9.488)
        # self.assertEqual(round(st.mu.muPas, 4), 123.6032)
        # self.assertEqual(round(st.k.mWmK, 4), 62.6479)
        self.assertEqual(round(st.mu.muPas, 4), 123.6033)
        self.assertEqual(round(st.k.mWmK, 4), 62.6478)

    def test_Surface(self):
        """Table 10, pag 271"""
        self.assertEqual(round(RE143a(T=340.1, x=0.5).sigma, 7), 0.0038511)
