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


class MD2M(MEoS):
    """Multiparameter equation of state for decamethyltetrasiloxane"""
    name = "decamethyltetrasiloxane"
    CASNumber = "141-62-8"
    formula = "C10H30Si4O3"
    synonym = "MD2M"
    _refPropName = "MD2M"
    _coolPropName = "MD2M"
    rhoc = unidades.Density(268.4321856)
    Tc = unidades.Temperature(599.4)
    Pc = unidades.Pressure(1144.0, "kPa")
    M = 310.6854  # g/mol
    Tt = unidades.Temperature(205.2)
    Tb = unidades.Temperature(467.59)
    f_acent = 0.635
    momentoDipolar = unidades.DipoleMoment(1.12, "Debye")
    # id = 1837

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [131.089725, -26.3839138],
           "ao_exp": [28.59, 56.42, 50.12],
           "titao": [20/Tc, 1180/Tc, 4240/Tc]}

    CP2 = {"ao": 4,
           "ao_exp": [28.59, 56.42, 50.12],
           "exp": [20, 1180, 4240]}

    f = 8.314472
    CP1 = {"ao": 331.9/f,
           "ao_sinh": [777.100/f], "sinh": [1813.8],
           "ao_cosh": [521.400/f], "cosh": [795.1]}

    thol = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for "
                    "decamethyltetrasiloxane of Thol (2017).",
        "__doi__": {
            "autor": "Thol, M., Dubberke, F.H., Baumhögger, E., Vrabec, J., "
                     "Span, R.",
            "title": "Speed of Sound Measuements and Fundamental Equations "
                     "State for Octamethyltrisiloxane and "
                     "Decamethyltetrasiloxane",
            "ref": "J. Chem. Eng. Data 62(9) (2017) 2633-2648",
            "doi": "10.1021/acs.jced.7b00092"},

        "R": 8.3144598,
        "cp": CP2,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 600, "Pmax": 130000,

        "nr1": [0.01458333, 3.227554, -3.503565, -2.017391, 0.8606129],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.319, 0.829, 0.78, 0.687],

        "nr2": [-2.196015, -0.9289014, 2.027740, -0.9168439, -0.06383507],
        "d2": [1, 3, 2, 2, 7],
        "t2": [1.29, 3.91, 0.77, 3.055, 1.013],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [2.674255, 0.04662529, -0.3835361, -0.4273462, -1.148009],
        "d3": [1, 1, 3, 2, 2],
        "t3": [1.07, 1.89, 1.133, 0.826, 0.83],
        "alfa3": [0.982, 2.7, 1.347, 0.864, 1.149],
        "beta3": [0.7323, 0.543, 1.26, 0.878, 2.22],
        "gamma3": [1.042, 1.1, 1.146, 1.085, 0.6844],
        "epsilon3": [0.874, 1.43, 0.855, 0.815, 0.491]}

    colonna = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for MD2M of Colonna (2006)",
        "__doi__": {"autor": "Colonna, P., Nannan, N.R., Guardone, A.",
                    "title": "Multiparameter equations of state for siloxanes:"
                             " [(CH3)3-Si-O1/2]2-[O-Si-(CH3)2]i=1,…,3, and "
                             "[O-Si-(CH3)2]6",
                    "ref": "Fluid Phase Equilibria 263:115-130, 2008",
                    "doi": "10.1016/j.fluid.2007.10.001"},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 673.0, "Pmax": 30000.0, "rhomax": 3.033,

        "nr1": [1.33840331, -2.62939393, 0.4398383, -0.53496715, 0.1818844,
                0.40774609e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [1.13444506, 0.5774631e-1, -0.5917498, -0.11020225,
                -0.34942635e-1, 0.7646298e-2],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.0],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = thol, colonna

    _surface = {
        "__doi__": {
            "autor": "Mulero, A., Cachadiña, I.",
            "title": "Recommended Correlations for the Surface Tension of "
                     "Several Fluids Included in the REFPROP Program",
            "ref": "J. Phys. Chem. Ref. Data 43(2) (2014) 023104",
            "doi": "10.1063/1.4878755"},
        "sigma": [0.0456], "exp": [1.41]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-10.174, 9.607, -10.08, -7.242, -30.56],
        "t": [1, 1.5, 1.83, 4.15, 17.8]}
    _liquid_Density = {
        "eq": 1,
        "n": [8.215, -24.65, 47.23, -42.44, 15.18],
        "t": [0.498, 0.855, 1.22, 1.6, 2.04]}
    _vapor_Density = {
        "eq": 2,
        "n": [-4.5483, -101.989, 224.06, -182.79, -110.45, -330.87],
        "t": [0.428, 2.32, 2.8, 3.3, 8.5, 17.5]}


class Test(TestCase):

    # Reference state fail with last decimal number
    def test_thol(self):
        # Table 13, Pag. 2646
        st = MD2M(T=350, rhom=0.0001)
        self.assertEqual(round(st.P.MPa, 10), 0.0002908572)
        # self.assertEqual(round(st.cpM.JmolK, 4), 467.6281)
        self.assertEqual(round(st.w, 5), 97.60521)
        self.assertEqual(round(st.hM.Jmol, 2), -19315.07)
        # self.assertEqual(round(st.sM.JmolK, 5), -12.27574)
        self.assertEqual(round(st.aM.Jmol, 2), -17927.15)

        st = MD2M(T=350, rhom=2.9)
        self.assertEqual(round(st.P.MPa, 5), 97.64338)
        # self.assertEqual(round(st.cpM.JmolK, 4), 568.1542)
        self.assertEqual(round(st.w, 3), 1295.433)
        # self.assertEqual(round(st.hM.Jmol, 2), -47644.93)
        self.assertEqual(round(st.sM.JmolK, 4), -208.0510)
        # self.assertEqual(round(st.aM.Jmol, 3), -8497.211)

        st = MD2M(T=500, rhom=0.05)
        self.assertEqual(round(st.P.MPa, 7), 0.1848863)
        # self.assertEqual(round(st.cpM.JmolK, 4), 587.2108)
        self.assertEqual(round(st.w, 4), 103.4085)
        # self.assertEqual(round(st.hM.Jmol, 2), 58103.75)
        self.assertEqual(round(st.sM.JmolK, 4), 118.4514)
        # self.assertEqual(round(st.aM.Jmol, 3), -4819.665)

        st = MD2M(T=500, rhom=2.7)
        self.assertEqual(round(st.P.MPa, 4), 115.7619)
        # self.assertEqual(round(st.cpM.JmolK, 4), 649.4756)
        self.assertEqual(round(st.w, 3), 1168.896)
        self.assertEqual(round(st.hM.Jmol, 2), 48958.28)
        # self.assertEqual(round(st.sM.JmolK, 6), 5.385695)
        self.assertEqual(round(st.aM.Jmol, 3), 3390.654)

        st = MD2M(T=590, rhom=2.5)
        self.assertEqual(round(st.P.MPa, 5), 89.78835)
        # self.assertEqual(round(st.cpM.JmolK, 4), 685.6283)
        self.assertEqual(round(st.w, 4), 987.7431)
        # self.assertEqual(round(st.hM.Jmol, 1), 101922.4)
        # self.assertEqual(round(st.sM.JmolK, 4), 121.0574)
        # self.assertEqual(round(st.aM.Jmol, 3), -5416.819)
