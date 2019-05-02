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

from lib import unidades
from lib.meos import MEoS


class iC6(MEoS):
    """Multiparameter equation of state for isohexane"""
    name = "isohexane"
    CASNumber = "107-83-5"
    formula = "(CH3)2-CH-(CH2)2-CH3"
    synonym = ""
    _refPropName = "IHEXANE"
    _coolPropName = "Isohexane"
    rhoc = unidades.Density(233.966)
    Tc = unidades.Temperature(497.7)
    Pc = unidades.Pressure(3040.0, "kPa")
    M = 86.17536  # g/mol
    Tt = unidades.Temperature(119.6)
    Tb = unidades.Temperature(333.36)
    f_acent = 0.2797
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 52

    Fi1 = {"ao_log": [1, 3.],
           "pow": [0, 1],
           "ao_pow": [6.9259123919, -0.3128629679],
           "ao_exp": [7.9127, 16.871, 19.257, 14.075],
           "titao": [325/Tc, 1150/Tc, 2397/Tc, 5893/Tc]}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for isohexane of "
                    "Lemmon and Span (2006)",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 "
                             "Industrial Fluids",
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi":  "10.1021/je050186n"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 550.0, "Pmax": 1000000.0, "rhomax": 9.38,

        "nr1": [1.1027, -2.9699, 1.0295, -0.21238, 0.11897, 0.00027738],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.40103, -0.034238, -0.43584, -0.11693, -0.019262, 0.0080783],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = lemmon,

    _surface = {"sigma": [0.05024], "exp": [1.194]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.74130e1, 0.16267e1, -0.22311e1, -0.26040e1, -0.29490e1],
        "t": [1.0, 1.5, 2.62, 4.56, 16.3]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.18489e2, -0.43541e2, 0.43985e2, -0.16581e2, 0.64563],
        "t": [0.59, 0.77, 0.96, 1.15, 3.2]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.41180e1, -0.61956e1, -0.21190e2, -0.58972e2, -0.15824e3],
        "t": [0.4824, 1.418, 3.32, 7.1, 16.1]}


class Test(TestCase):
    def test_shortLemmon(self):
        # Table 10, Pag 842
        st = iC6(T=499, rho=2*iC6.M)
        self.assertEqual(round(st.P.kPa, 3), 3058.917)
        self.assertEqual(round(st.hM.kJkmol, 3), 48733.740)
        self.assertEqual(round(st.sM.kJkmolK, 3), 113.316)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 233.627)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 1129.816)
        self.assertEqual(round(st.w, 3), 90.210)
