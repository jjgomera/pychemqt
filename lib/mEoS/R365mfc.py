#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2017, Juan José Gómez Romera <jjgomera@gmail.com>

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

from lib.meos import MEoS
from lib import unidades


class R365mfc(MEoS):
    """Multiparameter equation of state for R365mfc"""
    name = "1,1,1,3,3-pentafluorobutane"
    CASNumber = "406-58-6"
    formula = "CF3CH2CF2CH3"
    synonym = "R365mfc"
    _refPropName = "R365MFC"
    _coolPropName = "R365mfc"
    rhoc = unidades.Density(473.838464)
    Tc = unidades.Temperature(460.0)
    Pc = unidades.Pressure(3266.0, "kPa")
    M = 148.07452  # g/mol
    Tt = unidades.Temperature(239.0)
    Tb = unidades.Temperature(313.3)
    f_acent = 0.377
    momentoDipolar = unidades.DipoleMoment(3.807, "Debye")

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [-16.3423704513, 10.2889710846],
           "ao_exp": [17.47, 16.29],
           "titao": [569/Tc, 2232/Tc]}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-365mfc of Lemmon and "
                    "Span (2013)",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Thermodynamic Properties of R-227ea, R-365mfc, "
                             "R-115, and R-13I1",
                    "ref": "J. Chem. Eng. Data, 60(12) (2015) 3745-3758",
                    "doi": "10.1021/acs.jced.5b00684"},

        "R": 8.3144621,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 35000.0, "rhomax": 9.3,
        "Pmin": 2.478, "rhomin": 9.3,

        "nr1": [2.20027, -2.86240, 0.384559, -0.621227, 0.0665967],
        "d1": [1, 1, 2, 2, 4],
        "t1": [0.24, 0.67, 0.5, 1.25, 1],

        "nr2": [-1.19383, 0.635935, 0.461728, -0.533472, -1.07101, 0.139290],
        "d2": [1, 3, 6, 6, 2, 3],
        "t2": [3.35, 2.5, 0.96, 1.07, 5.6, 6.9],
        "c2": [1, 1, 1, 1, 2, 2],
        "gamma2": [1]*6,

        "nr3": [-0.385506, 0.885653, 0.226303, -0.166116],
        "d3": [1, 1, 1, 2],
        "t3": [3, 3.6, 5, 1.25],
        "alfa3": [0.97, 0.94, 2.15, 2.66],
        "beta3": [1.07, 1.08, 10.9, 22.6],
        "gamma3": [1.48, 1.49, 1.01, 1.16],
        "epsilon3": [1.02, 0.62, 0.53, 0.48]}

    eq = lemmon,

    _surface = {"sigma": [0.0534], "exp": [1.21]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.80955e1, 0.20414e1, -0.13333e2, 0.25514e2, -0.19967e2],
        "t": [1.0, 1.5, 3.4, 4.3, 5.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.17933e1, -0.18792e1, 0.90006e1, -0.11669e2, 0.56329e1],
        "t": [0.31, 0.6, 0.9, 1.2, 1.5]}
    _vapor_Density = {
        "eq": 2,
        "n": [-.1612e1, -.67679e1, -.24499e2, .33398e1, -.2111e3, .25807e3],
        "t": [0.281, 0.91, 3.0, 5.0, 8.0, 10.0]}


class Test(TestCase):

    def test_lemmon(self):
        # Table 7, Pag 3754
        st = R365mfc(T=300, rhom=0)
        self.assertEqual(round(st.P.MPa, 6), 0.0)
        self.assertEqual(round(st.cvM.JmolK, 4), 137.9013)
        self.assertEqual(round(st.cpM.JmolK, 4), 146.2158)
        self.assertEqual(round(st.w, 4), 133.6443)

        st = R365mfc(T=300, rhom=8.5)
        self.assertEqual(round(st.P.MPa, 6), 2.293261)
        self.assertEqual(round(st.cvM.JmolK, 4), 153.5866)
        self.assertEqual(round(st.cpM.JmolK, 4), 203.6507)
        self.assertEqual(round(st.w, 4), 735.3004)

        st = R365mfc(T=461, rhom=3.2)
        self.assertEqual(round(st.P.MPa, 6), 3.324581)
        self.assertEqual(round(st.cvM.JmolK, 4), 224.7732)
        self.assertEqual(round(st.cpM.JmolK, 2), 10987.87)
        self.assertEqual(round(st.w, 5), 68.71050)
