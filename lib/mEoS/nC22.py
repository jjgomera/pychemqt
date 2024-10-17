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
from lib.mEoS import nC12


class nC22(MEoS):
    """Multiparameter equation of state for n-hexadecane"""
    name = "n-docosane"
    CASNumber = "629-97-0"
    formula = "C22H46"
    synonym = ""
    _refPropName = ""
    _coolPropName = ""
    rhoc = unidades.Density(224.564523)
    Tc = unidades.Temperature(792.2)
    Pc = unidades.Pressure(1174., "kPa")
    M = 310.601  # g/mol
    Tt = unidades.Temperature(587.6)
    Tb = unidades.Temperature(641.75)
    f_acent = 0.9722
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    # id = 1765

    Fi1 = {"ao_log": [1, 32.9],
           "pow": [0, 1],
           "ao_pow": [66.73, -44.17],
           "ao_exp": [61.6, 77.7],
           "titao": [1000/Tc, 2400/Tc]}

    romeo = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for n-docosane of Romeo "
                    "(2022).",
        "__doi__": {"autor": "Romeo, R., Lemmon, E.W.",
                    "title": "Equations of State for n-Hexadecane and "
                             "n-Docosane",
                    "ref": "Int. J. Thermophys. 43 (2022) 146",
                    "doi": "10.1007/s10765-022-03059-0"},

        "R": 8.314462618,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": 298, "Tmax": 800.0, "Pmax": 200000.0, "rhomax": 10,

        "nr1": [0.0423455, 2.370432, -4.30263, -0.4039603, 0.4005704],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.224, 0.91, 0.95, 0.555],

        "nr2": [-2.643419, -0.9199641, 0.1394402, -1.448862, -0.0547678],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.36, 3.58, 0.5, 1.72, 1.078],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [4.579069, -0.3534636, -0.8217892, -0.2604273, -0.7618884],
        "d3": [1, 1, 3, 2, 2],
        "t3": [1.14, 2.43, 1.75, 1.1, 1.08],
        "alfa3": [0.641, 1.008, 1.026, 1.21, 0.93],
        "beta3": [0.516, 0.669, 0.25, 1.33, 2.1],
        "gamma3": [1.335, 1.187, 1.39, 1.23, 0.763],
        "epsilon3": [0.75, 1.616, 0.47, 1.306, 0.46]}

    eq = (romeo, )

    _surface = {
        "__doi__": {"autor": "Mulero, A., Cachadiña, I., Bautista, D.",
                    "title": "Recommended Correlations for the Surface "
                             "Tension of n-Alkanes",
                    "ref": "J. Phys. Chem. Ref. Data 50(2) (2021) 023104",
                    "doi": "10.1063/5.0048675"},
        "sigma": [0.0522], "exp": [1.25]}

    # FIXME: Vapor Pressure ancillary equation don't fit the values at low T
    _vapor_Pressure = {
        "eq": 3,
        "n": [-12.3834, 2.8818, -11.6292, -2.7357, -7.3103, 1188.9117],
        "t": [1.0, 1.5, 2.7, 5.5, 14.1, 52.1]}
    _liquid_Density = {
        "eq": 1,
        "n": [6.6254, -11.0123, 13.6452, -8.8244, 3.1241],
        "t": [0.5, 0.8, 1.2, 1.8, 2.5]}
    _vapor_Density = {
        "eq": 2,
        "n": [-5.979, -586.6421, -14.3725, -71.0676, -213.3123, 15.7901],
        "t": [0.5, 24, 1.7, 4.2, 10.3, 10.6]}

    trnECS = {"__name__": "Huber (2018)",

              "__doi__": {
                  "autor": "Huber, M.L.",
                  "title": "Models for Viscosity, Thermal Conductivity, and "
                           "Surface Tension of Selected Pure Fluids as "
                           "Implemented in REFPROP v10.0",
                  "ref": "NISTIR 8209",
                  "doi": "10.6028/NIST.IR.8209"},

              "eq": "ecs",
              "ref": nC12,

              "ek": 515.83, "sigma": 1.062, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [1.20571, -0.0689138],
              "psi_d": [0, 1],
              "fint": [0.00132], "fint_t": [0],
              "chi": [1.31627, -0.083506, 0.0127753], "chi_d": [0, 1, 2],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.31e-9, "gam0": 0.067, "qd": 1.114e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""

    def test_Huber(self):
        """Table 7, pag 266"""
        st = nC22(T=713, rhom=1.568)
        # self.assertEqual(round(st.mu.muPas, 4), 121.8525)
        # self.assertEqual(round(st.k.mWmK, 4), 74.8713)
        self.assertEqual(round(st.mu.muPas, 4), 121.6030)
        self.assertEqual(round(st.k.mWmK, 4), 74.8239)
