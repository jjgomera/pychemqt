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


class Acetone(MEoS):
    """Multiparameter equation of state for Acetone"""
    name = "acetone"
    CASNumber = "67-64-1"
    formula = "CH3COCH3"
    synonym = ""
    _refPropName = "ACETONE"
    _coolPropName = "Acetone"
    rhoc = unidades.Density(272.971958)
    Tc = unidades.Temperature(508.1)
    Pc = unidades.Pressure(4700.0, "kPa")
    M = 58.07914  # g/mol
    Tt = unidades.Temperature(178.5)
    Tb = unidades.Temperature(329.22)
    f_acent = 0.3071
    momentoDipolar = unidades.DipoleMoment(2.88, "Debye")
    id = 140

    Fi1 = {"ao_log": [1, 3.],
           "pow": [0, 1],
           "ao_pow": [-9.4883659997, 7.1422719708],
           "ao_exp": [3.7072, 7.0675, 11.012],
           "titao": [310/Tc, 3480/Tc, 1576/Tc]}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for acetone of Lemmon "
                    "and Span (2006)",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 "
                             "Industrial Fluids",
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi": "10.1021/je050186n"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 550.0, "Pmax": 700000.0, "rhomax": 15.73,
        "Pmin": 0.0023, "rhomin": 15.72,

        "nr1":  [0.90041, -2.1267, -0.083409, 0.065683, 0.00016527],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [-0.039663, 0.72085, 0.0092318, -0.17217, -0.14961, -0.076124,
                -0.018166],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    eq = lemmon,

    _surface = {"sigma": [0.0633], "exp": [1.16]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.76214e1, 0.17441e1, -0.20514e1, -0.26644e1, -0.69437],
        "t": [1, 1.5, 2.57, 4.43, 15.]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.11118e2, -0.29507e2, 0.35255e2, -0.14712e2, 0.95560],
        "t": [0.456, 0.626, 0.8, 1.0, 2.47]}
    _vapor_Density = {
        "eq": 2,
        "n": [-.25200e1, -.66065e1, -.25751e2, .78120e1, -.53778e2, -116.84],
        "t": [0.36, 1.05, 3.2, 4.0, 6.5, 14.0]}


class Test(TestCase):
    def test_shortLemmon(self):
        # Table 10, Pag 842
        st = Acetone(T=510, rhom=4)
        self.assertEqual(round(st.P.kPa, 3), 4807.955)
        self.assertEqual(round(st.hM.kJkmol, 3), 51782.004)
        self.assertEqual(round(st.sM.kJkmolK, 3), 157.331)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 138.449)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 3766.619)
        self.assertEqual(round(st.w, 3), 125.351)
