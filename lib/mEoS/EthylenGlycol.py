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


class EthylenGlycol(MEoS):
    """Multiparameter equation of state for Ethylene glycol"""
    name = "Ethylene glycol"
    CASNumber = "107-21-1"
    formula = "CH2OH-CH2OH"
    synonym = ""
    _refPropName = "EGLYCOL"
    _coolPropName = ""
    rhoc = unidades.Density(364.9588992)
    Tc = unidades.Temperature(719)
    Pc = unidades.Pressure(10508.7, "kPa")
    M = 62.06784  # g/mol
    Tt = unidades.Temperature(260.6)
    Tb = unidades.Temperature(470.45)
    f_acent = 0.619
    momentoDipolar = unidades.DipoleMoment(2.4103, "Debye")
    id = 135

    CP1 = {"ao": 4,
           "an": [], "pow": [],
           "ao_exp": [20.86], "exp": [1260.0]}

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

        "Tmin": 300, "Tmax": 850.0, "Pmax": 100000.0, "rhomax": 55.57,

        "nr1": [0.019393376, 1.2215576, 1.2751617, -3.6681302, -1.4660821,
                0.24628603, -0.063217756],
        "d1": [5.0, 1.0, 1.0, 1.0, 2.0, 3.0, 4.0],
        "t1": [1.0, 0.1, 1.27, 1.244, 1.1, 0.32, 1.0],

        "nr2": [1.4131488, 3.5245547, -2.2658015, 0.94972119, -0.13037392,
                0.19881857, -0.022141839],
        "d2": [2.0, 3.0, 3.0, 4.0, 5.0, 6.0, 7.0],
        "t2": [0.89, 1.2, 1.34, 1.3, 1.49, 1.23, 0.18],
        "c2": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        "gamma2": [1]*7,

        "nr3": [1.1273408, -0.12188623, -0.79487875, -0.024231918,
                -0.00574040155, 0.0083087704, -0.041852456],
        "d3": [1.0, 1.0, 2.0, 3.0, 4.0, 5.0, 3.0],
        "t3": [1.1, 0.75, 0.79, 0.77, 0.6, 1.0, 1.0],
        "alfa3": [0.9, 1.35, 0.8, 1.9, 2.0, 1.3, 20.0],
        "beta3": [0.91, 1.25, 0.97, 0.5, 1.0, 0.42, 1000.0],
        "gamma3": [1.17, 1.49, 1.3, 1.5, 1.2, 1.2, 1.07],
        "epsilon3": [1.37, 0.4, 1.0, 2.45, 1.9, 2.0, 0.9]}

    eq = (refprop, )

    _surface = {
        "__doi__": {"autor": "Mulero, A., Cachadiña, I., Sanjuán, E.L.",
                    "title": "Surface Tension of Alcohols. Data Selection and "
                             "Recommended Correlations",
                    "ref": "J. Phys. Chem. Ref. Data 44(3) (2015) 033104",
                    "doi": "10.1063/1.4927858"},
        "sigma": [0.0719], "exp": [0.755]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-10.193, 6.7548, -6.0404, -3.869, -8.4146],
        "t": [1.0, 1.5, 1.9, 3.64, 18.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [1.3902, -7.6323, 3.5326, 10.85, -8.7298, 3.3638],
        "t": [0.29, 1.25, 0.75, 1.8, 2.5, 3.4]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.47502, -5.5419, -16.175, -50.153, -86.316, -169.85],
        "t": [0.173, 0.6, 2.08, 4.8, 10.0, 18.0]}

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

              "ek": 570.95, "sigma": 0.4482, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [0.864168, -0.00230208, -0.00183225, 0.00904878],
              "psi_d": [0, 1, 2, 3],
              "fint": [0.00132], "fint_t": [0],
              "chi": [1.79177, -0.275354], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.166e-9, "gam0": 0.073, "qd": 0.542e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""

    def test_Huber(self):
        """Table 7, pag 266"""
        st = EthylenGlycol(T=647.1, rhom=12.365)
        # self.assertEqual(round(st.mu.muPas, 4), 124.6233)
        self.assertEqual(round(st.mu.muPas, 4), 124.6256)
        self.assertEqual(round(st.k.mWmK, 4), 204.0773)
