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


class Cyclopentane(MEoS):
    """Multiparameter equation of state for cyclopentane"""
    name = "cyclopropane"
    CASNumber = "287-92-3"
    formula = "C5H10"
    synonym = ""
    _refPropName = "CYCLOPEN"
    _coolPropName = "Cyclopentane"
    rhoc = unidades.Density(274.920968)
    Tc = unidades.Temperature(511.72)
    Pc = unidades.Pressure(4571.2, "kPa")
    M = 70.1329  # g/mol
    Tt = unidades.Temperature(179.7)
    Tb = unidades.Temperature(322.405)
    f_acent = 0.201
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 36

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [-0.3946233253, 2.4918910143],
           "ao_exp": [1.34, 13.4, 17.4, 6.65],
           "titao": [230/Tc, 1180/Tc, 2200/Tc, 5200/Tc],
           "ao_hyp": [], "hyp": []}

    gedanitz = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for cyclopentane of Gedanitz "
                    "et al. (2015)",
        "__doi__": {"autor": "Gedanitz, H., Davila, M.J., Lemmon, E.W.",
                    "title": "Speed of Sound Measurements and a Fundamental "
                             "Equation of State for Cyclopentane",
                    "ref": "J. Chem. Eng. Data, 60(5) (2015) 1331-1337",
                    "doi": "10.1021/je5010164"},

        "R": 8.3144621,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 550.0, "Pmax": 250000.0, "rhomax": 12.11,
        "Pmin": 0.008854, "rhomin": 12.1,

        "nr1": [0.0630928, 1.50365, -2.37099, -0.484886, 0.191843],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.29, 0.85, 1.185, 0.45],

        "nr2": [-0.835582, -0.435929, 0.545607, -0.209741, -0.0387635],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.28, 1.8, 1.5, 2.9, 0.93],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [0.677674, -0.137043, -0.0852862, -0.128085, -0.00389381],
        "d3": [1, 1, 3, 3, 2],
        "t3": [1.05, 4.0, 2.33, 1.5, 1.0],
        "alfa3": [0.86, 0.85, 0.86, 1.53, 5.13],
        "beta3": [0.63, 2.8, 0.5, 0.95, 0.23],
        "gamma3": [1.22, 0.32, 0.22, 1.94, 1.21],
        "epsilon3": [0.684, 0.7, 0.77, 0.625, 0.42]}

    eq = gedanitz,

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.1905, 1.8637, -1.6442, -2.72],
        "t": [1.0, 1.5, 5.5, 2.9]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.0741, 81.968, 173.88, -68.519, -184.74],
        "t": [0.1, 0.9, 1.25, 1.4, 1.05]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.0559, -6.4211, -46.926, 28.082, -70.838],
        "t": [0.1, 0.65, 3.2, 3.55, 7.5]}

    thermo0 = {"__name__": "Vassiliou (2015)",
               "__doi__": {
                   "autor": "Vassiliou, C.-M., Assael, M.J., Huber, M.L., "
                            "Perkins, R.A.",
                   "title": "Reference Correlation of the Thermal Conductivity"
                            " of Cyclopentane, iso-pentane, and n-Pentane",
                   "ref": "J. Phys. Chem. Ref. Data 44(3) (2015) 033102",
                   "doi": "10.1063/1.4927095"},

               "eq": 1,

               "Toref": 511.72, "koref": 1e-3,
               "no_num": [-8.2523346, 76.33654, -217.6154, 312.29877],
               "to_num": [0, 1, 2, 3],
               "no_den": [1, 0.28341478, 2.7890541, 0.32645005],
               "to_den": [0, 1, 2, 3],

               "Tref_res": 511.72, "rhoref_res": 274.921, "kref_res": 1,
               "nr": [0.0920536, -0.172699, 0.126557, -0.0362296, 0.00388718,
                      -0.0435129, 0.112636, -0.0908663, 0.028095, -0.00280368],
               "tr": [0, 0, 0, 0, 0, -1, -1, -1, -1, -1],
               "dr": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
               "Xio": 0.216e-9, "gam0": 0.058, "qd": 0.624e-9, "Tcref": 767.58}

    _thermal = thermo0,


class Test(TestCase):

    def test_gedanitz(self):
        # Table 5, Pag 1336
        st = Cyclopentane(T=330, rhom=0.01)
        self.assertEqual(round(st.P.MPa, 8), 0.02720379)
        self.assertEqual(round(st.cvM.JmolK, 5), 86.25843)
        self.assertEqual(round(st.cpM.JmolK, 5), 94.88857)
        self.assertEqual(round(st.w, 4), 205.6768)

        st = Cyclopentane(T=330, rhom=11)
        self.assertEqual(round(st.P.MPa, 5), 75.55974)
        self.assertEqual(round(st.cvM.JmolK, 4), 100.6003)
        self.assertEqual(round(st.cpM.JmolK, 4), 130.5568)
        self.assertEqual(round(st.w, 3), 1471.842)

        st = Cyclopentane(T=530, rhom=1)
        self.assertEqual(round(st.P.MPa, 6), 3.240235)
        self.assertEqual(round(st.cvM.JmolK, 4), 156.4924)
        self.assertEqual(round(st.cpM.JmolK, 4), 187.7243)
        self.assertEqual(round(st.w, 4), 195.4293)

        st = Cyclopentane(T=512, rhom=4)
        self.assertEqual(round(st.P.MPa, 6), 4.601539)
        self.assertEqual(round(st.cvM.JmolK, 4), 161.5786)
        self.assertEqual(round(st.cpM.JmolK, 2), 22857.91)
        self.assertEqual(round(st.w, 4), 113.0171)

        st = Cyclopentane(T=520, rhom=6)
        self.assertEqual(round(st.P.MPa, 6), 6.522373)
        self.assertEqual(round(st.cvM.JmolK, 4), 159.2304)
        self.assertEqual(round(st.cpM.JmolK, 4), 276.7530)
        self.assertEqual(round(st.w, 4), 234.2660)

    def test_Vassiliou(self):
        # Section 3.1.2, Pag 7
        # Viscosity value different to used in paper
        self.assertEqual(round(Cyclopentane(T=512, rho=400).k.mWmK, 3), 68.661)
