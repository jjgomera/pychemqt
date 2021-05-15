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
from lib.meos import MEoS


class nC12(MEoS):
    """Multiparameter equation of state for n-dodecane"""
    name = "dodecane"
    CASNumber = "112-40-3"
    formula = "CH3-(CH2)10-CH3"
    synonym = ""
    _refPropName = "C12"
    _coolPropName = "n-Dodecane"
    rhoc = unidades.Density(226.5453372)
    Tc = unidades.Temperature(658.1)
    Pc = unidades.Pressure(1817.0, "kPa")
    M = 170.33484  # g/mol
    Tt = unidades.Temperature(263.6)
    Tb = unidades.Temperature(489.3)
    f_acent = 0.574
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 16

    CP1 = {"ao": 23.085,
           "ao_exp": [37.776, 29.369, 12.461, 7.7733],
           "exp": [1280, 2399, 5700, 13869]}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for dodecane of Lemmon "
                    "(2004).",
        "__doi__": {"autor": "Lemmon, E.W., Huber, M.L.",
                    "title": "Thermodynamic Properties of n-Dodecane",
                    "ref": "Energy & Fuels, 18(4) (2004) 960-967",
                    "doi": "10.1021_ef0341062"},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 700., "Pmax": 700000.0, "rhomax": 4.53,

        "nr1": [1.38031, -2.85352, .288897, -0.165993, .0923993, .000282772],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.32, 1.23, 1.5, 1.4, 0.07, 0.8],

        "nr2": [.956627, .0353076, -0.445008, -0.118911, -0.0366475, .0184223],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [2.16, 1.1, 4.1, 5.6, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = lemmon,
    _PR = [0.1099, -26.8035]

    _surface = {
        "__doi__": {"autor": "Mulero, A., Cachadiña, I., Bautista, D.",
                    "title": "Recommended Correlations for the Surface "
                             "Tension of n-Alkanes",
                    "ref": "J. Phys. Chem. Ref. Data 50(2) (2021) 023104",
                    "doi": "10.1063/5.0048675"},
        "sigma": [0.04733, 0.01131], "exp": [1.192, 2.977]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.94217e1, -0.41890e1, 0.54999e1, -0.67789e1, -0.17161e1],
        "t": [1.0, 1.5, 1.359, 3.56, 9.2]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.92236, 0.92047, 0.55713e1, -0.92253e1, 0.51763e1],
        "t": [0.21, 0.49, 1.08, 1.49, 1.9]}
    _vapor_Density = {
        "eq": 2,
        "n": [-1.7859, -7.5436, -22.848, -81.355, 92.283, -217.25],
        "t": [0.298, 0.91, 2.8, 6., 9., 11.]}

    visco0 = {"__name__": "Huber (2004)",
              "__doi__": {
                  "autor": "Huber, M.L., Laesecke, A., Perkins, R.A.",
                  "title": "Transport Properties of n-Dodecane",
                  "ref": "Energy & Fuels 18(4) (2004) 968-975.",
                  "doi": "10.1021/ef034109e"},

              "eq": 1, "omega": 1,

              "ek": 522.592, "sigma": 0.735639,
              "n_chapman": 0.021357,
              "collision": [0.382987, -0.561050, 0.0313962],

              "Tref_virial": 522.592,
              "n_virial": [-19.572881, 219.73999, -1015.3226, 2471.0125,
                           -3375.1717, 2491.6597, -787.26086, 14.085455,
                           -0.34664158],
              "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],

              "Tref_res": 658.1, "rhoref_res": 1.33*M, "muref_res": 1000,
              "nr": [-0.0471703, 0.00827816, 0.0298429, -0.0134156],
              "tr": [1, 1, 2, 2],
              "dr": [2, 3, 2, 3],

              "CPf": 503.109,
              "CPg1": 2.32661,
              "CPgi": [2.23089/2.32661],
              "CPti": [-0.5]}

    _viscosity = visco0,

    thermo0 = {"eq": 1,
               "__name__": "Huber (2004)",
               "__doi__": {
                   "autor": "Huber, M.L., Laesecke, A., Perkins, R.A.",
                   "title": "Transport Properties of n-Dodecane",
                   "ref": "Energy & Fuels 18(4) (2004) 968-975.",
                   "doi": "10.1021/ef034109e"},

               "Toref": 658.1, "koref": 1.,
               "no": [0.436343e-2, -0.264054e-1, 0.922394e-1, -0.291756e-1],
               "to": [0, 1, 2, 3],

               "Tref_res": 658.1, "rhoref_res": 1.33*M, "kref_res": 1.,
               "nr": [0.693347e-1, -0.280792e-1, -0.331695e-1, 0.173922e-2,
                      0.676165e-2, 0.309558e-2],
               "tr": [0, -1, 0, -1, 0, -1],
               "dr": [1, 1, 2, 2, 3, 3],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
               "Xio": 0.194e-9, "gam0": 0.0496, "qd": 1.52e-9, "Tcref": 987.15}

    _thermal = thermo0,


class Test(TestCase):

    def test_lemmon(self):
        # Table 5, Pag 967
        st = nC12(T=300, rho=0)
        self.assertEqual(round(st.P.MPa, 3), 0)
        self.assertEqual(round(st.cvM.JmolK, 4), 271.3952)
        self.assertEqual(round(st.cpM.JmolK, 4), 279.7096)
        self.assertEqual(round(st.w, 4), 122.8511)

        st = nC12(T=300, rhom=4.4)
        self.assertEqual(round(st.P.MPa, 6), 7.161946)
        self.assertEqual(round(st.cvM.JmolK, 4), 311.1937)
        self.assertEqual(round(st.cpM.JmolK, 4), 377.1307)
        self.assertEqual(round(st.w, 3), 1316.614)

        st = nC12(T=600, rhom=0.2)
        self.assertEqual(round(st.P.MPa, 7), 0.7228018)
        self.assertEqual(round(st.cvM.JmolK, 4), 488.1173)
        self.assertEqual(round(st.cpM.JmolK, 4), 526.2322)
        self.assertEqual(round(st.w, 4), 123.5851)

        st = nC12(T=600, rhom=4.0)
        self.assertEqual(round(st.P.MPa, 4), 120.2532)
        self.assertEqual(round(st.cvM.JmolK, 4), 499.9857)
        self.assertEqual(round(st.cpM.JmolK, 4), 539.3864)
        self.assertEqual(round(st.w, 3), 1255.247)

        st = nC12(T=658.2, rhom=1.33)
        self.assertEqual(round(st.P.MPa, 6), 1.820055)
        self.assertEqual(round(st.cvM.JmolK, 4), 547.6852)
        self.assertEqual(round(st.cpM.JmolK, 1), 128753.9)
        self.assertEqual(round(st.w, 5), 49.76424)

    def test_Huber(self):
        # Viscosity test point, Pag 972
        self.assertEqual(round(nC12(T=300, rhom=4.4115).mu.muPas, 1), 1484.8)
        self.assertEqual(round(nC12(T=500, P=1e6).mu.muPas, 2), 183.76)

        # Thermal conductivity test point, Pag 974
        self.assertEqual(round(nC12(T=300, P=1e7).k.mWmK, 2), 138.29)
        self.assertEqual(round(nC12(T=500, P=1e6).k.mWmK, 2), 93.84)
        self.assertEqual(round(nC12(T=660, P=1.8714e6).k.mWmK, 3), 90.346)
