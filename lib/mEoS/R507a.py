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
from lib.meos import MEoSBlend


class R507a(MEoSBlend):
    """Multiparameter equation of state for R507A (50% R125, 50% R143a)"""
    name = "R507A"
    CASNumber = ""
    formula = "R125+R143a"
    synonym = "R507A"
    _refPropName = "R507A"
    _coolPropName = "R507A"
    rhoc = unidades.Density(490.7370688)
    Tc = unidades.Temperature(343.765)
    Pc = unidades.Pressure(3704.9, "kPa")
    M = 98.8592  # g/mol
    Tt = unidades.Temperature(200.0)
    Tb = unidades.Temperature(226.41)
    f_acent = 0.286
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")

    Fi1 = {"ao_log": [1, -1],
           "pow": [0, 1, -0.25],
           "ao_pow": [9.93541, 7.9985, -21.6054],
           "ao_exp": [0.95006, 4.1887, 5.5184],
           "titao": [364/Tc, 815/Tc, 1768/Tc]}

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

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 50000.0, "rhomax": 14.13,
        "Pmin": 23.23, "rhomin": 14.13,

        "Tj": 343.765, "Pj": 3.7049,
        "dew": {"i": [1*2, 1.5*2, 2.1*2, 4.7*2],
                "n": [-7.5459, 2.338, -2.237, -4.1535]},
        "bubble": {"i": [1*2, 1.5*2, 2.2*2, 4.6*2],
                   "n": [-7.4853, 2.0115, -2.0141, -3.7763]},

        "nr1": [0.624982e1, -0.807855e1, 0.264843e-1, 0.286215, -0.507076e-2,
                0.109552e-1, 0.116124e-2],
        "d1": [1, 1, 1, 2, 2, 4, 6],
        "t1": [0.692, 0.943, 5.8, 0.77, 5.84, 0.24, 0.69],

        "nr2": [0.138469e1, -0.922473, -0.503562e-1, 0.822098, -0.277727,
                0.358172, -0.126426e-1, -0.607010e-2, -.815653e-1, -.233323e-1,
                .352952e-1, .159566e-1, .755927e-1, -.542007e-1, .170451e-1],
        "d2": [1, 1, 1, 2, 2, 3, 4, 7, 2, 3, 4, 4, 2, 3, 5],
        "t2": [2, 3, 7, 2.2, 4.3, 2.7, 1.2, 1.23, 12, 6, 8.5, 11.5, 13, 17,
               16.2],
        "c2": [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3],
        "gamma2": [1]*15}

    eq = lemmon,

    _surface = {"sigma": [0.06701, -0.04297], "exp": [1.3066, 2.3145]}

    _liquid_Density = {
        "eq": 1,
        "n": [-1.1909272837212852, 2.4301888224201074, 3.512171521080571,
              -3.7013727137258634, 2.358426436855713, -3.9469718132522256],
        "t": [0.093, 0.162, 0.914, 1.425, 2.269, 7.339]}
    _vapor_Density = {
        "eq": 3,
        "n": [3.975317800427674, -74.1368924572824, 75.98726158051952,
              -11.56733447604713, 3.084500064878586, -6.722628041295302],
        "t": [0.11, 0.26, 0.286, 0.563, 2.03, 2.811]}

    thermo0 = {"__name__": "Geller (2001)",
               "__doi__": {
                   "autor": "Geller, V.Z., Nemzer, B.V., Cheremnykh, U.V.",
                   "title": "Thermal  Conductivity of the Refrigerant "
                            "Mixtures R404A, R407C, R410A and R507A",
                   "ref": "Int. J. Termophysics 22(4) (2001) 1035-1043",
                   "doi": "10.1023/a_1010691504352"},

               "eq": 1,

               "Toref": 1, "koref": 1e-3,
               "no": [-8.656, 7.383e-2],
               "to": [0, 1],

               "rhoref_res": 1, "kref_res": 1e-3,
               "nr": [2.799e-2, 3.065e-5, -3.644e-8, 2.609e-11],
               "tr": [0, 0, 0, 0],
               "dr": [1, 2, 3, 4],

               "critical": 0}

    _thermal = thermo0,


class Test(TestCase):

    def test_lemmon(self):
        # Table V, Pag 998
        st = R507a(T=300, rhom=0)
        self.assertEqual(round(st.P.MPa, 3), 0)
        self.assertEqual(round(st.cvM.JmolK, 3), 76.838)
        self.assertEqual(round(st.cpM.JmolK, 3), 85.152)
        self.assertEqual(round(st.w, 2), 167.22)

        st = R507a(T=300, P=R507a._bubbleP(300))
        self.assertEqual(round(st.P.MPa, 4), 1.3462)
        self.assertEqual(round(st.rhoM, 5), 10.50670)
        self.assertEqual(round(st.cvM.JmolK, 3), 91.290)
        self.assertEqual(round(st.cpM.JmolK, 2), 153.79)
        self.assertEqual(round(st.w, 2), 356.72)

        st = R507a(T=300, P=R507a._dewP(300))
        self.assertEqual(round(st.P.MPa, 4), 1.3450)
        self.assertEqual(round(st.rhoM, 5), 0.73552)
        self.assertEqual(round(st.cvM.JmolK, 3), 88.742)
        self.assertEqual(round(st.cpM.JmolK, 2), 124.15)
        self.assertEqual(round(st.w, 2), 130.95)

        st = R507a(T=250, rhom=13)
        self.assertEqual(round(st.P.MPa, 3), 12.916)
        self.assertEqual(round(st.cvM.JmolK, 3), 83.541)
        self.assertEqual(round(st.cpM.JmolK, 2), 123.21)
        self.assertEqual(round(st.w, 2), 697.89)
