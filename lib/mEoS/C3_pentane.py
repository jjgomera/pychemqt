#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2022, Juan José Gómez Romera <jjgomera@gmail.com>

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


class C3_pentane(MEoS):
    """Multiparameter equation of state for 3-methylpentane"""
    name = "3-Methylpentane"
    CASNumber = "96-14-0"
    formula = "CH3-CH2-CH(CH3)-CH2-CH3"
    synonym = ""
    _refPropName = ""
    _coolPropName = ""
    rhoc = unidades.Density(239.5675008)
    Tc = unidades.Temperature(506)
    Pc = unidades.Pressure(3184.5, "kPa")
    M = 86.17536  # g/mol
    Tt = unidades.Temperature(110.263)
    Tb = unidades.Temperature(336.379)
    f_acent = 0.268
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 53

    Fi1 = {"ao_log": [1, 6],
           "pow": [0, 1],
           "ao_pow": [4.64793693054, -0.906514506457],
           "ao_exp": [7.4047, 25.09, 17.741],
           "titao": [470/Tc, 1555/Tc, 3946/Tc]}

    gao = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for 3-Methylpentane of Gao"
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

        "Tmin": Tt, "Tmax": 550, "Pmax": 1200000, "rhomax": 9.66,

        "nr1": [0.006178288, 0.763315017, -0.5546657, -1.0604327, 0.23117181],

        "d1": [5, 1, 1, 2, 3],
        "t1": [1, 0.16, 1.0, 1.0, 0.386],

        "nr2": [-1.8757299, -0.9327912, 1.0720552, -0.2830806, -0.024600061],
        "d2": [1, 3, 2, 2, 8],
        "t2": [1.54, 2.0, 1.0, 2.5, 1.66],
        "c2": [2, 2, 1, 2, 2],
        "gamma2": [1]*5,

        "nr3": [0.87360786, 0.008687374, -0.27160944, -0.12365512, -0.12052593,
                -0.53359397],
        "d3": [1, 1, 3, 2, 1, 3],
        "t3": [0.44, 1.0, 0.55, 0.705, 1.5, 1.0],
        "alfa3": [1.409, 2.53, 1.781, 2.045, 0.688, 20.1],
        "beta3": [1.876, 1.158, 1.808, 1.646, 1.0, 660.0],
        "gamma3": [1.2603, 1.207, 1.045, 1.069, 0.923, 1.109],
        "epsilon3": [0.7065, 2.19, 0.244, 1.014, 0.689, 0.905]}

    eq = gao,

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.3854, 1.5058, -1.3741, -3.1976, -1.3433],
        "t": [1.0, 1.5, 2.75, 4.0, 15.2]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.4911, 4.2494, -3.965, 2.0896, 0.18507],
        "t": [0.19, 0.645, 1.05, 1.6, 7.5]}
    _vapor_Density = {
        "eq": 2,
        "n": [-1.407, -6.4038, -19.274, -53.302, -108.9, -200.54],
        "t": [0.3, 0.871, 2.95, 6.5, 14.0, 27.0]}


class Test(TestCase):

    def test_Gao(self):
        # Table 13, pag 22
        st = C3_pentane(T=265, rhom=8.3)
        self.assertEqual(round(st.P.MPa, 5), 34.55032)
        self.assertEqual(round(st.cvM.JmolK, 4), 137.1533)
        self.assertEqual(round(st.cpM.JmolK, 4), 176.3509)
        self.assertEqual(round(st.w, 3), 1421.979)
        self.assertEqual(round(st.hM.kJmol, 5), -10.80643)
        self.assertEqual(round(st.sM.JmolK, 5), -50.64068)

        st = C3_pentane(T=465, rhom=5.5)
        self.assertEqual(round(st.P.MPa, 6), 3.513147)
        self.assertEqual(round(st.cvM.JmolK, 4), 208.5403)
        self.assertEqual(round(st.cpM.JmolK, 4), 278.1085)
        self.assertEqual(round(st.w, 4), 394.1040)
        self.assertEqual(round(st.hM.kJmol, 5), 30.76257)
        self.assertEqual(round(st.sM.JmolK, 5), 75.55693)

        st = C3_pentane(T=506.5, rhom=2.78)
        self.assertEqual(round(st.P.MPa, 6), 3.207690)
        self.assertEqual(round(st.cvM.JmolK, 4), 256.1725)
        self.assertEqual(round(st.cpM.JmolK, 2), 23439.18)
        self.assertEqual(round(st.w, 5), 80.18092)
        self.assertEqual(round(st.hM.kJmol, 5), 47.39300)
        self.assertEqual(round(st.sM.JmolK, 4), 109.6177)

        st = C3_pentane(T=265, rhom=0)
        self.assertEqual(round(st.P.MPa, 5), 0)
        self.assertEqual(round(st.cvM.JmolK, 4), 118.0127)
        self.assertEqual(round(st.cpM.JmolK, 4), 126.3272)
        self.assertEqual(round(st.w, 4), 165.4369)
        self.assertEqual(round(st.hM.kJmol, 5), 18.44498)

        st = C3_pentane(T=485, rhom=0.37)
        self.assertEqual(round(st.P.MPa, 6), 1.215479)
        self.assertEqual(round(st.cvM.JmolK, 4), 207.5595)
        self.assertEqual(round(st.cpM.JmolK, 4), 226.9415)
        self.assertEqual(round(st.w, 4), 181.4910)
        self.assertEqual(round(st.hM.kJmol, 5), 53.55585)
        self.assertEqual(round(st.sM.JmolK, 4), 127.1686)

        st = C3_pentane(T=525, rhom=3.8)
        self.assertEqual(round(st.P.MPa, 6), 4.791053)
        self.assertEqual(round(st.cvM.JmolK, 4), 230.9643)
        self.assertEqual(round(st.cpM.JmolK, 4), 395.2492)
        self.assertEqual(round(st.w, 4), 165.7838)
        self.assertEqual(round(st.hM.kJmol, 5), 49.29010)
        self.assertEqual(round(st.sM.JmolK, 4), 112.3864)
