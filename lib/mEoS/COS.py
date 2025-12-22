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
from lib.mEoS.C3 import C3


class COS(MEoS):
    """Multiparameter equation of state for carbonyl sulfide"""
    name = "carbonyl sulfide"
    CASNumber = "463-58-1"
    formula = "COS"
    synonym = ""
    _refPropName = "COS"
    _coolPropName = "CarbonylSulfide"
    rhoc = unidades.Density(445.1565)
    Tc = unidades.Temperature(378.77)
    Pc = unidades.Pressure(6370.0, "kPa")
    M = 60.0751  # g/mol
    Tt = unidades.Temperature(134.3)
    Tb = unidades.Temperature(222.99)
    f_acent = 0.0978
    momentoDipolar = unidades.DipoleMoment(0.7152, "Debye")
    id = 219

    Fi1 = {"ao_log": [1, 2.5],
           "pow": [0, 1],
           "ao_pow": [-3.6587449805, 3.7349245016],
           "ao_exp": [2.1651, 0.93456, 1.0623, 0.34269],
           "titao": [768/Tc, 1363/Tc, 3175/Tc, 12829/Tc]}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for carbonyl sulfide "
                    "of Lemmon and Span (2006)",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 "
                             "Industrial Fluids",
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi": "10.1021/je050186n"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 650., "Pmax": 50000.0, "rhomax": 25,

        "nr1": [0.94374, -2.5348, 0.59058, -0.021488, 0.082083, 0.00024689],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.21226, -0.041251, -0.22333, -0.050828, -0.028333, 0.016983],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = (lemmon, )

    _surface = {"sigma": [0.07246], "exp": [1.407]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.67055e1, 0.34248e1, -0.26677e1, -0.24717e1],
        "t": [1., 1.5, 1.78, 4.8]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.76592e1, -0.19226e2, 0.27883e2, -0.23637e2, 0.99803e1],
        "t": [0.515, 0.767, 1.034, 1.4, 1.7]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.32494e1, -0.71460e1, 35.026, -34.039, -64.206, -152.25],
        "t": [0.423, 1.464, 5.3, 4.1, 7.0, 17.0]}

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

              "ek": 335, "sigma": 0.413, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [1], "psi_d": [0],
              "fint": [0.00125], "fint_t": [0],
              "chi": [0.95], "chi_d": [0],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.182e-9, "gam0": 0.056, "qd": 0.5e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""
    def test_shortLemmon(self):
        """Table 10, Pag 842"""
        st = COS(T=380, rho=7*COS.M)
        self.assertEqual(round(st.P.kPa, 3), 6498.429)
        self.assertEqual(round(st.hM.kJkmol, 3), 16511.877)
        self.assertEqual(round(st.sM.kJkmolK, 3), 51.563)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 55.861)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 4139.577)
        self.assertEqual(round(st.w, 3), 161.717)

    def test_Huber(self):
        """Table 7, pag 266"""
        st = COS(T=340.9, rhom=14.426)
        self.assertEqual(round(st.mu.muPas, 4), 104.7952)
        self.assertEqual(round(st.k.mWmK, 4), 72.6218)
