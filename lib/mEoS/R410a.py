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


class R410a(MEoSBlend):
    """Multiparameter equation of state for R410A (50% R32, 50% R125)"""
    name = "R410A"
    CASNumber = ""
    formula = "R32+R125"
    synonym = "R410A"
    _refPropName = "R410A"
    _coolPropName = "R410A"
    rhoc = unidades.Density(459.0300696)
    Tc = unidades.Temperature(344.494)
    Pc = unidades.Pressure(4901.2, "kPa")
    M = 72.5854  # g/mol
    Tt = unidades.Temperature(200.0)
    Tb = unidades.Temperature(221.71)
    f_acent = 0.296
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")

    Fi1 = {"ao_log": [1, -1],
           "pow": [0, 1, -0.1],
           "ao_pow": [36.8871, 7.15807, -46.87575],
           "ao_exp": [2.0623, 5.9751, 1.5612],
           "titao": [697/Tc, 1723/Tc, 3875/Tc]}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-410A of Lemmon (2003)",
        "__doi__": {"autor": "Lemmon, E.W.",
                    "title": "Pseudo-Pure Fluid Equations of State for the "
                             "Refrigerant Blends R-410A, R-404A, R-507A, and "
                             "R-407C",
                    "ref": "Int. J. Thermophys., 24(4) (2003) 991-1006",
                    "doi": "10.1023/A:1025048800563"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 50000.0, "rhomax": 19.51,

        "Tj": 344.494, "Pj": 4.9012,
        "dew": {"i": [1*2, 1.6*2, 2.4*2, 5*2],
                "n": [-7.4411, 1.9883, -2.4925, -3.2633]},
        "bubble": {"i": [1*2, 1.8*2, 2.4*2, 4.9*2],
                   "n": [-7.2818, 2.5093, -3.2695, -2.8022]},

        "nr1": [0.987252, -0.103017e1, 0.117666e1, -0.138991, 0.302373e-2],
        "d1": [1, 1, 1, 2, 5],
        "t1": [0.44, 1.2, 2.97, 2.95, 0.2],

        "nr2": [-0.253639e1, -0.196680e1, -0.830480, 0.172477, -0.261116,
                -0.745473e-1, 0.679757, -0.652431, 0.553849e-1, -0.710970e-1,
                -0.875332e-3, 0.200760e-1, -0.139761e-1, -0.185110e-1,
                0.171939e-1, -0.482049e-2],
        "d2": [1, 2, 3, 5, 5, 5, 1, 1, 4, 4, 9, 2, 2, 4, 5, 6],
        "t2": [1.93, 1.78, 3., 0.2, 0.74, 3, 2.1, 4.3, 0.25, 7, 4.7, 13, 16,
               25, 17, 7.4],
        "c2": [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3],
        "gamma2": [1]*16}

    eq = lemmon,

    _surface = {
        "__doi__": {
            "autor": "Okada, M., Shibata, T., Sato. Y., Higashi, Y.",
            "title": "Surface Tension of HFC Refrigerant Mixtures",
            "ref": "Int. J. Thermophysics 20(1) (1999) 119-127",
            "doi": "10.1023/A:1021482231102"},
        "sigma": [0.06238], "exp": [1.246]}

    _liquid_Density = {
        "eq": 1,
        "n": [17.62376356934577, -34.49248423230781, 229.4790425881445,
              -332.5139608181856, 23058076.978355706, -31889422.797179773],
        "t": [0.732, 1.165, 2.775, 3.462, 15.531, 15.944]}
    _vapor_Density = {
        "eq": 3,
        "n": [-25.00709732482156, 34.83593078562833, -33.0203056956001,
              61.19129776625034, -63.32155175461899, 149.8912653231937],
        "t": [0.676, 0.844, 1.394, 2.411, 3.149, 7.808]}

    thermo0 = {"__name__": "Geller (2001)",
               "__doi__": {
                   "autor": "Geller, V.Z., Nemzer, B.V., Cheremnykh, U.V.",
                   "title": "Thermal  Conductivity of the Refrigerant "
                            "Mixtures R404A, R407C, R410A and R507A",
                   "ref": "Int. J. Thermophysics 22(4) (2001) 1035-1043",
                   "doi": "10.1023/a_1010691504352"},

               "eq": 1,

               "Toref": 1, "koref": 1e-3,
               "no": [-8.872, 7.410e-2],
               "to": [0, 1],

               "rhoref_res": 1, "kref_res": 1e-3,
               "nr": [3.576e-2, -9.045e-6, 4.343e-8, -3.705e-12],
               "tr": [0, 0, 0, 0],
               "dr": [1, 2, 3, 4],

               "critical": 0}

    _thermal = thermo0,


class Test(TestCase):

    def test_lemmon(self):
        # Table V, Pag 998
        st = R410a(T=300, rhom=0)
        self.assertEqual(round(st.P.MPa, 3), 0)
        self.assertEqual(round(st.cvM.JmolK, 3), 50.400)
        self.assertEqual(round(st.cpM.JmolK, 3), 58.714)
        self.assertEqual(round(st.w, 2), 200.08)

        st = R410a(T=300, P=R410a._bubbleP(300))
        self.assertEqual(round(st.P.MPa, 4), 1.7404)
        self.assertEqual(round(st.rhoM, 5), 14.45917)
        self.assertEqual(round(st.cvM.JmolK, 3), 67.147)
        self.assertEqual(round(st.cpM.JmolK, 2), 125.50)
        self.assertEqual(round(st.w, 2), 418.60)

        st = R410a(T=300, P=R410a._dewP(300))
        self.assertEqual(round(st.P.MPa, 4), 1.7351)
        self.assertEqual(round(st.rhoM, 5), 0.95997)
        self.assertEqual(round(st.cvM.JmolK, 3), 67.335)
        self.assertEqual(round(st.cpM.JmolK, 2), 107.60)
        self.assertEqual(round(st.w, 2), 160.94)

        st = R410a(T=250, rhom=18)
        self.assertEqual(round(st.P.MPa, 3), 17.651)
        self.assertEqual(round(st.cvM.JmolK, 3), 62.521)
        self.assertEqual(round(st.cpM.JmolK, 3), 98.401)
        self.assertEqual(round(st.w, 2), 800.83)
