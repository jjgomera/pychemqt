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


class oH2(MEoS):
    """Multiparameter equation of state for hydrogen (orto)"""
    name = "ortohydrogen"
    CASNumber = "1333-74-0o"
    formula = "H2"
    synonym = "R-702o"
    _refPropName = "ORTHOHYD"
    _coolPropName = "OrthoHydrogen"
    rhoc = unidades.Density(31.1361933)
    Tc = unidades.Temperature(33.22)
    Pc = unidades.Pressure(1310.65, "kPa")
    M = 2.01594  # g/mol
    Tt = unidades.Temperature(14.008)
    Tb = unidades.Temperature(20.4)
    f_acent = -0.219
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")

    Fi1 = {"ao_log": [1, 1.5],
           "pow": [0, 1],
           "ao_pow": [-1.4675442336, 1.8845068862],
           "ao_exp": [2.54151, -2.3661, 1.00365, 1.22447],
           "titao": [856/Tc, 1444/Tc, 2194/Tc, 6968/Tc]}

    leachman = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ortohydrogen of Leachman "
                    "et al. (2007)",
        "__doi__": {"autor": "Leachman, J.W., Jacobsen, R.T, Penoncello, S.G.,"
                             " Lemmon, E.W.",
                    "title": "Fundamental equations of state for Parahydrogen,"
                             " Normal Hydrogen, and Orthohydrogen",
                    "ref": "J. Phys. Chem. Ref. Data, 38(3) (2009) 721-748",
                    "doi": "10.1063/1.3160306"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 1000.0, "Pmax": 2000000.0, "rhomax": 38.2,

        "nr1": [-6.83148, .01, 2.11505, 4.38353, 0.211292, -1.00939, 0.142086],
        "d1": [1, 4, 1, 1, 2, 2, 3],
        "t1": [0.7333, 1, 1.1372, 0.5136, 0.5638, 1.6248, 1.829],

        "nr2": [-0.87696, 0.804927],
        "d2": [1, 3],
        "t2": [2.404, 2.105],
        "c2": [1, 1],
        "gamma2": [1]*2,

        "nr3": [-0.710775, 0.0639688, 0.0710858, -0.087654, 0.647088],
        "d3": [2, 1, 3, 1, 1],
        "t3": [4.1, 7.658, 1.259, 7.589, 3.946],
        "alfa3": [1.169, 0.894, 0.04, 2.072, 1.306],
        "beta3": [0.4555, 0.4046, 0.0869, 0.4415, 0.5743],
        "gamma3": [1.5444, 0.6627, 0.763, 0.6587, 1.4327],
        "epsilon3": [0.6366, 0.3876, 0.9437, 0.3976, 0.9626],
        "nr4": []}

    eq = (leachman, )

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.488684e1, 0.105310e1, 0.856947, -0.185355],
        "t": [1.0, 1.5, 2.7, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.43911e1, -0.75872e1, 0.10402e2, -0.72651e1, 0.18302e1],
        "t": [0.53, 0.93, 1.35, 1.8, 2.4]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.1463, -16.183, 31.803, -219.61, 431.23, -255.91],
        "t": [0.491, 2.1, 2.9, 4.4, 5.0, 5.5]}


class Test(TestCase):
    """Testing"""

    def test_leachman(self):
        """Selected point from Table 14, Pag 746, saturation states"""
        st = oH2(T=14.008, x=0.5)
        self.assertEqual(round(st.P.kPa, 4), 7.5601)
        self.assertEqual(round(st.Liquido.rho, 3), 77.010)
        self.assertEqual(round(st.Gas.rho, 5), 0.13273)
        self.assertEqual(round(st.Liquido.h.kJkg, 3), -53.820)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 400.77)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), -3.0625)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), 29.390)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 5.1746)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 6.2707)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 7.1448)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 10.557)
        self.assertEqual(round(st.Liquido.w, 1), 1264.7)
        self.assertEqual(round(st.Gas.w, 2), 307.38)

        st = oH2(T=20, x=0.5)
        self.assertEqual(round(st.P.kPa, 3), 90.417)
        self.assertEqual(round(st.Liquido.rho, 3), 71.291)
        self.assertEqual(round(st.Gas.rho, 4), 1.1977)
        self.assertEqual(round(st.Liquido.h.kJkg, 4), -3.7690)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 448.30)
        self.assertEqual(round(st.Liquido.s.kJkgK, 5), -0.17907)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), 22.425)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 5.6050)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 6.4172)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 9.5449)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 11.813)
        self.assertEqual(round(st.Liquido.w, 1), 1137.1)
        self.assertEqual(round(st.Gas.w, 2), 354.41)

        st = oH2(T=30, x=0.5)
        self.assertEqual(round(st.P.kPa, 2), 804.89)
        self.assertEqual(round(st.Liquido.rho, 3), 54.554)
        self.assertEqual(round(st.Gas.rho, 3), 10.481)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 140.19)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 440.48)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 5.0602)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), 15.070)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 6.2748)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 7.6849)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 25.315)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 31.669)
        self.assertEqual(round(st.Liquido.w, 2), 719.08)
        self.assertEqual(round(st.Gas.w, 2), 378.99)

        st = oH2(T=33, x=0.5)
        self.assertEqual(round(st.P.kPa, 2), 1269.00)
        self.assertEqual(round(st.Liquido.rho, 3), 38.828)
        self.assertEqual(round(st.Gas.rho, 3), 23.810)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 251.54)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 349.69)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 8.2561)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), 11.230)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 7.4490)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 8.6394)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 2), 325.68)
        self.assertEqual(round(st.Gas.cp.kJkgK, 2), 406.90)
        self.assertEqual(round(st.Liquido.w, 2), 438.23)
        self.assertEqual(round(st.Gas.w, 2), 370.01)
