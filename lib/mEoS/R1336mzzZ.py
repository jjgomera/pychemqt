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


class R1336mzzZ(MEoS):
    """Multiparameter equation of state for R1336mzzZ"""
    name = "cis-1,1,1,4,4,4-hexafluorobutene"
    CASNumber = "692-49-9"
    formula = "CF3CH=CHCF3"
    synonym = "R-1336mzz(Z)"
    rhoc = unidades.Density(499.386464)
    Tc = unidades.Temperature(444.5)
    Pc = unidades.Pressure(2903, "kPa")
    M = 164.056  # g/mol
    Tt = unidades.Temperature(0)

    # Kontomaris, K.
    # HFO-1336mzz-Z: High Temperature Chemical Stability and Use as A Working
    # Fluid in Organic Rankine Cycles
    # International Refrigeration and Air Condictioning Conference. Paper 1525
    # http://docs.lib.purdue.edu/iracc/1525
    Tb = unidades.Temperature(306.55)

    f_acent = 0
    momentoDipolar = unidades.DipoleMoment(0, "Debye")

    CP1 = {"ao": 4,
           "ao_exp": [20.2, 5.275],
           "exp": [736, 2299]}

    mclinden = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R1336mzz(Z) of McLinden"
                    "(2020)",
        "__doi__": {"autor": "McLinden, M.O., Akasaka, R.",
                    "title": "Thermodynamic Properties of cis-1,1,1,4,4,4-"
                             "Hexafluorobutene [R-1234ze(Z)]: Vapor Pressure, "
                             "(p, ρ, T) Behavior, and Speed of Sound "
                             "Measurements and Equation of State",
                    "ref": "J. Chem. Eng. Data 65(9) (2020) 4201-4214",
                    "doi": "10.1021/acs.jced.9b01198"},

        "R": 8.314462618,
        "cp": CP1,
        "ref": "IIR",

        "Tmin": 23, "Tmax": 460, "Pmax": 36000.0,

        "nr1": [0.036673095, 1.1956619, -1.8462376, -0.60599297, 0.24973833],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.26, 1, 1, 0.515],

        "nr2": [-1.2548278, -1.4389612, 0.35168887, -0.82104051, -0.031747538],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.6, 3.0, 0.74, 2.68, 0.96],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [1.0281388, 0.21094074, 0.701701, 0.24638528, -1.5295034,
                0.33424978, 1.011324, -0.023457179],
        "d3": [1, 1, 3, 2, 3, 2, 2, 1],
        "t3": [1.06, 3.4, 1.617, 1.865, 1.737, 3.29, 1.242, 2],
        "alfa3": [0.746, 2.406, 0.7804, 1.25, 0.6826, 1.677, 1.762, 21],
        "beta3": [1.118, 3.065, 0.7274, 0.8435, 0.6754, 0.436, 3.808, 1888],
        "gamma3": [0.962, 1.111, 1.135, 1.163, 0.969, 1.286, 1.274, 1.056],
        "epsilon3": [1.225, 0.161, 1.231, 1.395, 0.9072, 0.958, 0.412, 0.944]}

    eq = (mclinden, )

    # TODO: Calculate ancillary saturation equations, the paper don't give any
    # Using meanwhile the correlations for R1243zf
    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.3009, 1.5186, -2.5303, -1.5139],
        "t": [1.0, 1.5, 2.8, 5]}
    _liquid_Density = {
        "eq": 1,
        "n": [2.3870, -1.2213, 4.8105, -5.7706],
        "t": [0.362, 0.83, 1.31, 1.82, 2.4]}
    _vapor_Density = {
        "eq": 2,
        "n": [-2.9588, -7.1173, -21.945, -54.068, -128.95],
        "t": [0.385, 1.25, 3.4, 7.25, 16]}


class Test(TestCase):
    """Test class"""

    def test_mclinden(self):
        """Table S.1, pag 2"""
        st = R1336mzzZ(T=350, rhom=0)
        self.assertEqual(round(st.P.MPa, 6), 0)
        self.assertEqual(round(st.cvM.kJmolK, 6), 0.145277)
        self.assertEqual(round(st.cpM.kJmolK, 6), 0.153591)
        self.assertEqual(round(st.w, 3), 136.943)

        st = R1336mzzZ(T=350, rhom=0.1)
        self.assertEqual(round(st.P.MPa, 7), 0.2662916)
        self.assertEqual(round(st.cvM.kJmolK, 6), 0.149683)
        self.assertEqual(round(st.cpM.kJmolK, 6), 0.162997)
        self.assertEqual(round(st.w, 3), 126.681)

        st = R1336mzzZ(T=350, rhom=8)
        self.assertEqual(round(st.P.MPa, 5), 20.88061)
        self.assertEqual(round(st.cvM.kJmolK, 6), 0.160341)
        self.assertEqual(round(st.cpM.kJmolK, 6), 0.206406)
        self.assertEqual(round(st.w, 3), 617.566)

        st = R1336mzzZ(T=400, rhom=0.5)
        self.assertEqual(round(st.P.MPa, 6), 1.201815)
        self.assertEqual(round(st.cvM.kJmolK, 6), 0.170639)
        self.assertEqual(round(st.cpM.kJmolK, 6), 0.210627)
        self.assertEqual(round(st.w, 3), 108.203)

        st = R1336mzzZ(T=400, rhom=7)
        self.assertEqual(round(st.P.MPa, 5), 11.82407)
        self.assertEqual(round(st.cvM.kJmolK, 6), 0.169977)
        self.assertEqual(round(st.cpM.kJmolK, 6), 0.222507)
        self.assertEqual(round(st.w, 3), 425.612)

        st = R1336mzzZ(T=445, rhom=3.1)
        self.assertEqual(round(st.P.MPa, 6), 2.930174)
        self.assertEqual(round(st.cvM.kJmolK, 6), 0.219506)
        self.assertEqual(round(st.cpM.kJmolK, 4), 19.8479)
        self.assertEqual(round(st.w, 4), 60.4525)
