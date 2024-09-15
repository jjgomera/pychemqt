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
from lib.meos import MEoSBlend


class R404a(MEoSBlend):
    """Multiparameter equation of state for R404A
    (44% R125, 4% R134a, 52% R143a)"""

    name = "R404A"
    CASNumber = ""
    formula = "R125+R134a+R143a"
    synonym = "R404A"
    _refPropName = "R404A"
    _coolPropName = "R404A"
    rhoc = unidades.Density(482.162772)
    Tc = unidades.Temperature(345.27)
    Pc = unidades.Pressure(3734.8, "kPa")
    M = 97.6038  # g/mol
    Tt = unidades.Temperature(200.0)
    Tb = unidades.Temperature(226.93)
    f_acent = 0.293
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")

    Fi1 = {"ao_log": [1, -1],
           "pow": [0, 1, -0.3],
           "ao_pow": [7.00407, 7.98695, -18.8664],
           "ao_exp": [0.63078, 3.5979, 5.0335],
           "titao": [413/Tc, 804/Tc, 1727/Tc]}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-404A of Lemmon (2003)",
        "__doi__": {"autor": "Lemmon, E.W.",
                    "title": "Pseudo-Pure Fluid Equations of State for the "
                             "Refrigerant Blends R-410A, R-404A, R-507A, and "
                             "R-407C",
                    "ref": "Int. J. Thermophys., 24(4) (2003) 991-1006",
                    "doi": "10.1023/A:1025048800563"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": 250, "Tmax": 500.0, "Pmax": 50000.0, "rhomax": 14.5,

        "Tj": 345.270, "Pj": 3.7348,
        "dew": {"i": [0.1*2, 0.972*2, 3.8*2, 9.0*2],
                "n": [-0.00026863, -6.5757, -4.1802, -7.9102]},
        "bubble": {"i": [0.54*2, 0.965*2, 3.7*2, 9.0*2],
                   "n": [0.061067, -6.5646, -3.6162, 3.9771]},

        "nr1": [6.10984, -7.79453, 0.0183377, 0.262270, -0.00351688, 0.0116181,
                0.00105992],
        "d1": [1, 1, 1, 2, 2, 4, 6],
        "t1": [0.67, 0.91, 5.96, 0.7, 6, 0.3, 0.7],

        "nr2": [0.850922, -0.520084, -0.0464225, 0.62119, -0.195505, 0.336159,
                -0.0376062, -0.00636579, -0.0758262, -0.0221041, 0.0310441,
                0.0132798, 0.0689437, -0.0507525, 0.0161382],
        "d2": [1, 1, 1, 2, 2, 3, 4, 7, 2, 3, 4, 4, 2, 3, 5],
        "t2": [1.7, 3.3, 7, 2.05, 4.3, 2.7, 1.8, 1.25, 12, 6, 8.7, 11.6, 13,
               17, 16],
        "c2": [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3],
        "gamma2": [1]*15}

    eq = (lemmon, )

    _surface = {
        "__doi__": {
            "autor": "Okada, M., Shibata, T., Sato. Y., Higashi, Y.",
            "title": "Surface Tension of HFC Refrigerant Mixtures",
            "ref": "Int. J. Thermophysics 20(1) (1999) 119-127",
            "doi": "10.1023/a:1021482231102"},
        "sigma": [0.0536], "exp": [1.249]}

    _liquid_Density = {
        "eq": 1,
        "n": [-0.011777381218307114, 2.4984547682678704, 5367.932468376006,
              -3732982.2139762687, 33474569.015986394, -80262239.93679684],
        "t": [0.066, 0.377, 8.99, 13.98, 15.982, 17.954]}
    _vapor_Density = {
        "eq": 3,
        "n": [-5.211849919504925, 6.798109425092868, -13.038945243933242,
              9.654798311425768, -320.35292894913766, 80023.04183470446],
        "t": [0.485, 0.989, 1.431, 2.495, 7.812, 15.329]}

    thermo0 = {"__name__": "Geller (2001)",
               "__doi__": {
                   "autor": "Geller, V.Z., Nemzer, B.V., Cheremnykh, U.V.",
                   "title": "Thermal  Conductivity of the Refrigerant "
                            "Mixtures R404A, R407C, R410A and R507A",
                   "ref": "Int. J. Termophysics 22(4) (2001) 1035-1043",
                   "doi": "10.1023/a_1010691504352"},

               "eq": 1,

               "Toref": 1, "koref": 1e-3,
               "no": [-8.624, 7.36e-2],
               "to": [0, 1],

               "rhoref_res": 1, "kref_res": 1e-3,
               "nr": [3.222e-2, 2.569e-5, -2.693e-8, 2.007e-11],
               "tr": [0, 0, 0, 0],
               "dr": [1, 2, 3, 4],

               "critical": 0}

    _thermal = (thermo0, )


class Test(TestCase):
    """Testing"""
    def test_lemmon(self):
        """Table V, Pag 998"""
        st = R404a(T=300, rhom=0)
        self.assertEqual(round(st.P.MPa, 3), 0)
        self.assertEqual(round(st.cvM.JmolK, 3), 76.219)
        self.assertEqual(round(st.cpM.JmolK, 3), 84.533)
        self.assertEqual(round(st.w, 2), 168.36)

        st = R404a(T=300, P=R404a._bubbleP(300))
        self.assertEqual(round(st.P.MPa, 4), 1.3169)
        self.assertEqual(round(st.rhoM, 5), 10.60497)
        self.assertEqual(round(st.cvM.JmolK, 3), 90.653)
        self.assertEqual(round(st.cpM.JmolK, 2), 152.11)
        self.assertEqual(round(st.w, 2), 365.11)

        st = R404a(T=300, P=R404a._dewP(300))
        self.assertEqual(round(st.P.MPa, 4), 1.3034)
        self.assertEqual(round(st.rhoM, 5), 0.70599)
        self.assertEqual(round(st.cvM.JmolK, 3), 87.917)
        self.assertEqual(round(st.cpM.JmolK, 2), 121.86)
        self.assertEqual(round(st.w, 2), 132.92)

        st = R404a(T=250, rhom=13)
        self.assertEqual(round(st.P.MPa, 3), 10.435)
        self.assertEqual(round(st.cvM.JmolK, 3), 83.062)
        self.assertEqual(round(st.cpM.JmolK, 2), 123.27)
        self.assertEqual(round(st.w, 2), 688.27)
