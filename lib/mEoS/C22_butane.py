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


class C22_butane(MEoS):
    """Multiparameter equation of state for 2,2-dimethylbutane"""
    name = "2,2-dimethylbutane"
    CASNumber = "75-83-2"
    formula = "CH3-C(CH2)2-CH2-CH3"
    synonym = ""
    _refPropName = ""
    _coolPropName = ""
    rhoc = unidades.Density(239.5675008)
    Tc = unidades.Temperature(490.0)
    Pc = unidades.Pressure(3138.0, "kPa")
    M = 86.17536  # g/mol
    Tt = unidades.Temperature(174.20)
    Tb = unidades.Temperature(322.846)
    f_acent = 0.230
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 54

    Fi1 = {"ao_log": [1, 6.5],
           "pow": [0, 1],
           "ao_pow": [4.07450146142, -1.02091722997],
           "ao_exp": [7.8764, 26.017, 22.147],
           "titao": [525/Tc, 1620/Tc, 4370/Tc]}

    gao = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for 2,2-DiMethylbutane of Gao"
                    "et al. (2021)",
        "__doi__": {
            "autor": "Gao, K., Wu, J., Lemmon, E.W.",
            "title": "Equations of State for the Thermodynamic Properties of "
                     "Three Hexane Isomers: 3-Methylpantane, "
                     "2,2-Dimethylbutane, and 2,3-Dimethylbutane",
            "ref": "J. Phys. Chem. Ref. Data 50(3) (2021) 033103",
            "doi": "10.1063/1.5093644"},

        "R": 8.314462618,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 575, "Pmax": 1000000, "rhomax": 8.77,

        "nr1": [0.00702066, 0.70134226, -0.3659372, -1.109303, 0.22742868],
        "d1": [5, 1, 1, 2, 3],
        "t1": [1.0, 0.156, 1.0, 1.0, 0.371],

        "nr2": [-1.8603613, -0.65052551, 1.1465612, -0.31514795, -0.028916258],
        "d2": [1, 3, 2, 2, 8],
        "t2": [1.4, 2.0, 1.0, 2.15, 1.5],
        "c2": [2, 2, 1, 2, 2],
        "gamma2": [1]*5,

        "nr3": [0.9153258, -0.010020802, -0.52298297, -0.15308943,
                -0.21698526, -1.1808573],
        "d3": [1, 1, 3, 2, 1, 3],
        "t3": [0.49, 1.4, 0.687, 1.275, 1.48, 1.0],
        "alfa3": [1.35, 1.278, 1.35, 1.724, 1.042, 27.0],
        "beta3": [1.709, 0.218, 1.19, 0.33, 2.18, 1074.0],
        "gamma3": [1.275, 0.91, 1.108, 1.184, 1.174, 1.094],
        "epsilon3": [0.7384, 2.063, 0.239, 1.057, 0.558, 0.926]}

    eq = gao,

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.4088, 2.5218, -1.5652, -3.4318, -1.221],
        "t": [1.0, 1.5, 2.0, 3.85, 15.85]}
    _liquid_Density = {
        "eq": 1,
        "n": [1.852, -1.405, 6.393, -7.718, 3.708],
        "t": [0.33, 0.65, 0.95, 1.3, 1.7]}
    _vapor_Density = {
        "eq": 2,
        "n": [-1.7198, -5.7667, -17.244, -51.992, -126.59],
        "t": [0.32, 0.904, 2.85, 6.4, 14.75]}


class Test(TestCase):

    def test_Gao(self):
        # Table 13, pag 22
        st = C22_butane(T=250, rhom=8.3)
        self.assertEqual(round(st.P.MPa, 5), 36.16117)
        self.assertEqual(round(st.cvM.JmolK, 4), 131.5927)
        self.assertEqual(round(st.cpM.JmolK, 4), 168.4855)
        self.assertEqual(round(st.w, 3), 1433.881)
        self.assertEqual(round(st.hM.kJmol, 5), -10.34456)
        self.assertEqual(round(st.sM.JmolK, 5), -52.29832)

        st = C22_butane(T=450, rhom=5.5)
        self.assertEqual(round(st.P.MPa, 6), 3.848580)
        self.assertEqual(round(st.cvM.JmolK, 4), 205.7891)
        self.assertEqual(round(st.cpM.JmolK, 4), 268.0357)
        self.assertEqual(round(st.w, 4), 392.0591)
        self.assertEqual(round(st.hM.kJmol, 5), 29.46773)
        self.assertEqual(round(st.sM.JmolK, 5), 74.87199)

        st = C22_butane(T=490.5, rhom=2.78)
        self.assertEqual(round(st.P.MPa, 6), 3.161800)
        self.assertEqual(round(st.cvM.JmolK, 4), 246.6004)
        self.assertEqual(round(st.cpM.JmolK, 2), 19167.87)
        self.assertEqual(round(st.w, 5), 82.45625)
        self.assertEqual(round(st.hM.kJmol, 5), 45.45018)
        self.assertEqual(round(st.sM.JmolK, 4), 108.8356)

        st = C22_butane(T=250, rhom=0)
        self.assertEqual(round(st.P.MPa, 5), 0)
        self.assertEqual(round(st.cvM.JmolK, 4), 113.9448)
        self.assertEqual(round(st.cpM.JmolK, 4), 122.2593)
        self.assertEqual(round(st.w, 4), 160.8752)
        self.assertEqual(round(st.hM.kJmol, 5), 16.76637)

        st = C22_butane(T=470, rhom=0.37)
        self.assertEqual(round(st.P.MPa, 6), 1.192083)
        self.assertEqual(round(st.cvM.JmolK, 4), 204.3634)
        self.assertEqual(round(st.cpM.JmolK, 4), 222.5305)
        self.assertEqual(round(st.w, 4), 180.5417)
        self.assertEqual(round(st.hM.kJmol, 5), 51.31643)
        self.assertEqual(round(st.sM.JmolK, 4), 126.2985)

        st = C22_butane(T=510, rhom=3.8)
        self.assertEqual(round(st.P.MPa, 6), 4.816117)
        self.assertEqual(round(st.cvM.JmolK, 4), 229.2782)
        self.assertEqual(round(st.cpM.JmolK, 4), 393.3932)
        self.assertEqual(round(st.w, 4), 162.7817)
        self.assertEqual(round(st.hM.kJmol, 5), 47.55079)
        self.assertEqual(round(st.sM.JmolK, 4), 112.0508)
