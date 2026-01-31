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


class DMC(MEoS):
    """Multiparameter equation of state for dimethyl carbonate"""
    name = "dimethyl carbonate"
    CASNumber = "616-38-6"
    formula = "C3H6O3"
    synonym = ""
    _refPropName = "DMC"
    _coolPropName = "DimethylCarbonate"
    rhoc = unidades.Density(360.3116)
    Tc = unidades.Temperature(557.)
    Pc = unidades.Pressure(4908.8, "kPa")
    M = 90.0779  # g/mol
    Tt = unidades.Temperature(277.06)
    Tb = unidades.Temperature(363.256)
    f_acent = 0.346
    momentoDipolar = unidades.DipoleMoment(0.899, "Debye")
    # id=1798

    Fi1 = {"ao_log": [1, 8.28421],
           "pow": [0, 1],
           "ao_pow": [4.9916462, -0.1709449],
           "ao_exp": [1.48525, 0.822585, 16.2453, 1.15925],
           "titao": [21/Tc, 1340/Tc, 1672/Tc, 7395/Tc]}

    zhou = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for DMC of Zhou (2011).",
        "__doi__": {"autor": "Zhou, Y., Wu, J., Lemmon, E.W.",
                    "title": "Thermodynamic Properties of Dimethyl Carbonate",
                    "ref": "J. Phys. Chem. Ref. Data, Vol. 40, No. 4 2011",
                    "doi": "10.1063/1.3664084"},
        "R": 8.314472,
        "cp": Fi1,
        "ref": {"Tref": 298.15, "Pref": 1.0, "ho": 26712.371, "so": 109.66202},

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 60000.0, "rhomax": 12.107,

        "nr1": [0.52683187e-3, 1.353396, -2.649283, -0.2785412, 0.1742554,
                0.031606252],
        "d1": [5, 1, 1, 2, 3, 4],
        "t1": [1, 0.227, 1.05, 1.06, 0.5, 0.78],

        "nr2": [0.399866, 1.178144, -0.0235281, -1.015, -0.7880436, -0.12696],
        "d2": [1, 2, 7, 1, 2, 3],
        "t2": [1.3, 1.347, 0.706, 2, 2.5, 4.262],
        "c2": [1, 1, 1, 2, 2, 2],
        "gamma2": [1]*6,

        "nr3": [1.2198, -0.4883, -0.0033293, -0.0035387, -0.51172, -0.16882],
        "d3": [1, 1, 2, 2, 3, 3],
        "t3": [1, 2.124, 0.4, 3.5, 0.5, 2.7],
        "alfa3": [0.9667, 1.5154, 1.0591, 1.6642, 12.4856, 0.9662],
        "beta3": [1.24, 0.821, 15.45, 2.21, 437., 0.743],
        "gamma3": [1.2827, 0.4317, 1.1217, 1.1871, 1.1243, 0.4203],
        "epsilon3": [0.6734, 0.9239, 0.8636, 1.0507, 0.8482, 0.7522]}

    eq = (zhou, )

    _surface = {
        "__doi__": {
            "autor": "Mulero, A., Cachadiña, I.",
            "title": "Recommended Correlations for the Surface Tension of "
                     "Several Fluids Included in the REFPROP Program",
            "ref": "J. Phys. Chem. Ref. Data 43(2) (2014) 023104",
            "doi": "10.1063/1.4878755"},
        "sigma": [0.0825], "exp": [1.39]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-8.3197, 3.4260, -3.5905, -3.3194],
        "t": [1.0, 1.5, 2.3, 4.7]}
    _liquid_Density = {
        "eq": 1,
        "n": [1.1572, 4.969, -14.451, 27.569, -26.223, 10.526],
        "t": [0.27, 0.77, 1.29, 1.85, 2.46, 3.16]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.54715, -5.19277, -94.048, 327.21, -676.871, 716.072,
              -379.799],
        "t": [0.197, 0.6, 2.86, 3.65, 4.5, 5.4, 6.4]}

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

              "ek": 442.3, "sigma": 0.51, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [0.811428, 0.0616704], "psi_d": [0, 1],
              "fint": [0.00132], "fint_t": [0],
              "chi": [1.1238, -0.0154353], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.204e-9, "gam0": 0.059, "qd": 0.62e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""

    def test_Huber(self):
        """Table 7, pag 266"""
        st = DMC(T=501.3, rhom=8.3)
        self.assertEqual(round(st.mu.muPas, 4), 102.0731)
        self.assertEqual(round(st.k.mWmK, 4), 92.2091)
