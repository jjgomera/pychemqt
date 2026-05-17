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


class THF(MEoS):
    """Multiparameter equation of state for tetrahydrofuran"""
    name = "tetrahydrofuran"
    CASNumber = "109-99-9"
    formula = "C4H8O"
    synonym = ""
    rhoc = unidades.Density(317.265168)
    Tc = unidades.Temperature(540.2)
    Pc = unidades.Pressure(5304.5, "kPa")
    M = 72.10572  # g/mol
    Tt = unidades.Temperature(164.76)
    Tb = unidades.Temperature(339.075)
    f_acent = 0.234
    momentoDipolar = unidades.DipoleMoment(1.63, "Debye")
    id = 281

    CP1 = {"ao": 4,
           "ao_exp": [18.2, 11.394, 1.05, 2.37],
           "exp": [1460, 3461, 11000, 517]}

    fiedler = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for tetrahydrofuran of "
                    "Fiedler (2023)",
        "__doi__": {"autor": "Fiedler, F., Karog, J., Lemmon, E.W., Thol, M.",
                    "title": "Fundamental Equation of State for Fluid "
                             "Tetrahydrofuran",
                    "ref": "Int. J. Thermophysics 44 (2023) 153",
                    "doi": "10.1007/s10765-023-03258-3"},

        "R": 8.314462618,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 550, "Pmax": 600000,

        "nr1": [0.04386, 0.766, -1.2355036286776, -0.6899995453364, 0.201742],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.12, 0.94, 1.111, 0.41],

        "nr2": [-0.7603, -0.3754, 0.5317, -0.0354, -0.02196],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.25, 2.77, 0.88, 2.71, 0.85],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [-0.0399, -0.0112, -0.4165, 0.6293, -0.03702],
        "d3": [1, 2, 3, 2, 1],
        "t3": [0.87, 1, 1.035, 0.95, 2.26],
        "alfa3": [1.88, 25, 0.85, 0.81, 0.86],
        "beta3": [2.5, 900, 0.8, 0.79, 1.3],
        "gamma3": [0.85, 1.08, 1.34, 1.33, 1.38],
        "epsilon3": [1, 0.93, 0.59, 0.73, 0.56]}

    eq = (fiedler, )

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.82, 4.1666, -3.43, -0.805, -2.417],
        "t": [1, 1.5, 2, 3.45, 5]}
    _liquid_Density = {
        "eq": 1,
        "n": [6.9, -8.7784, 7.87, -5.75, 2.59],
        "t": [0.5254, 0.782, 1.286, 1.94, 2.5]}
    _vapor_Density = {
        "eq": 2,
        "n": [-4.557, -8.9253, -4.585, -27.86, -60.2, 140],
        "t": [0.4897, 1.82, 3.1, 4.7, 8.9, 18]}


class Test(TestCase):
    """Test class"""

    def test_akasaka(self):
        """Table 6, pag 23"""
        st = THF(T=270, rhom=1e-7)
        self.assertEqual(round(st.P.MPa, 10), 0.0000002245)
        self.assertEqual(round(st.cpM.JmolK, 9), 67.962189060)
        self.assertEqual(round(st.w, 9), 188.343575267)
        self.assertEqual(round(st.hM.Jmol, 6), 24850.336600)
        self.assertEqual(round(st.sM.JmolK, 8), 179.51672775)
        self.assertEqual(round(st.aM.Jmol, 6), -25864.084531)

        st = THF(T=350, rhom=0.04)
        self.assertEqual(round(st.P.MPa, 10), 0.1130755690)
        self.assertEqual(round(st.cpM.JmolK, 9), 93.540905193)
        self.assertEqual(round(st.w, 9), 205.594935289)
        self.assertEqual(round(st.hM.Jmol, 6), 30997.192355)
        self.assertEqual(round(st.sM.JmolK, 9), 90.484710729)
        self.assertEqual(round(st.aM.Jmol, 6), -3499.345627)

        st = THF(T=450, rhom=10)
        self.assertEqual(round(st.P.MPa, 9), 12.357974600)
        self.assertEqual(round(st.cpM.JmolK, 8), 167.23826646)
        self.assertEqual(round(st.w, 9), 739.195761440)
        self.assertEqual(round(st.hM.Jmol, 6), 17240.429690)
        self.assertEqual(round(st.sM.JmolK, 9), 40.905432334)
        self.assertEqual(round(st.aM.Jmol, 6), -2402.812320)

        st = THF(T=550, rhom=5)
        self.assertEqual(round(st.P.MPa, 10), 6.1720378363)
        self.assertEqual(round(st.cpM.JmolK, 8), 763.57251979)
        self.assertEqual(round(st.w, 9), 139.994309340)
        self.assertEqual(round(st.hM.Jmol, 6), 40236.271714)
        self.assertEqual(round(st.sM.JmolK, 9), 87.703653509)
        self.assertEqual(round(st.aM.Jmol, 6), -9235.145283)
