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


class Novec649(MEoS):
    """Multiparameter equation of state for Novec649"""
    name = "Novec649"
    CASNumber = "756-13-8"
    formula = "C6F12O"
    synonym = "1,1,1,2,2,4,5,5,5-nonafluoro-4-(trifluromethyl)-3-pentanone"
    _refPropName = "NOVEC649"
    _coolPropName = "Novec649"
    rhoc = unidades.Density(606.805248)
    Tc = unidades.Temperature(441.81)
    Pc = unidades.Pressure(1.869, "MPa")
    M = 316.0444  # g/mol
    Tt = unidades.Temperature(165)
    Tb = unidades.Temperature(322.202)
    f_acent = 0.471
    momentoDipolar = unidades.DipoleMoment(0.43, "Debye")

    Fi1 = {"ao_log": [1, 29.8],
           "pow": [0, 1],
           "ao_pow": [-30.6610503233, 6.8305296372],
           "ao_exp": [29.8],
           "titao": [1940/Tc]}

    mclinden = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for Novec649 of McLinden "
                    "(2015).",
        "__doi__": {
            "autor": "McLinden, M.O., Perkins, R.A., Lemmon, E.W., Fortin, "
                     "T.J.",
            "title": "Thermodynamic Properties of 1,1,1,2,2,4,5,5,5-nonafluoro"
                     "-4-(trifluoromethyl)-3-pentanone: Vapor Pressure, (p, "
                     "ρ, T) Behavior, and Speed of Sound Measurements, and an "
                     "Equation of State",
            "ref": "J. Chem. Eng. Data 60(12) (2015) 3646-3659",
            "doi": "10.1021/acs.jced.5b00623"},

        "R": 8.3144621,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 50000.0, "rhomax": 6.24,

        "nr1": [0.05623648, 2.973616, -6.126970, 3.440240, 1.451737,
                -2.837857, 0.2077767],
        "d1": [4, 1, 1, 1, 2, 2, 3],
        "t1": [1, 0.25, 0.793, 1.16, 0.75, 1.09, 0.75],

        "nr2": [2.168307, -2.124648, -1.296704],
        "d2": [2, 1, 2],
        "t2": [1.3, 2.25, 1.9],
        "c2": [1, 2, 2],
        "gamma2": [1]*3,

        "nr3": [-1.010569, 2.701505, 0.8167202, -1.814579, 0.2075389,
                -1.009347, -0.04848043],
        "d3": [1, 1, 2, 2, 3, 3, 1],
        "t3": [0.88, 1.63, 1.3, 2.0, 1.15, 1.66, 1.5],
        "alfa3": [0.32, 1.32, 1.35, 1.48, 0.51, 1.30, 5.15],
        "beta3": [0.12, 0.83, 0.19, 0.95, 0.1, 0.11, 65.0],
        "gamma3": [1.10, 1.04, 1.15, 0.9, 0.8, 1.2, 1.19],
        "epsilon3": [1.16, 0.793, 1.13, 0.527, 1.19, 0.83, 0.82]}

    eq = mclinden,

    _vapor_Pressure = {
        "eq": 3,
        "n": [-8.4411, 2.711, -3.6354, -5.3872, -8.1641],
        "t": [1.0, 1.5, 2.2, 4.4, 15]}
    _liquid_Density = {
        "eq": 1,
        "n": [1.5545, 1.149, 0.51565],
        "t": [0.297, 0.7, 4.4]}
    _vapor_Density = {
        "eq": 2,
        "n": [-1.6073, -5.8095, -17.824, -61.012, -151.3],
        "t": [0.291, 0.82, 2.45, 5.5, 12]}

    def _Visco0(self, T, coef=None):
        """Special dilute gas viscosity correlation"""
        muo = MEoS._Visco0(self, T)

        # Correction parameter
        mur = 131.3*self.momentoDipolar/(self.Tc/self.rhoc)**0.5         # Eq 7
        Fc = 1 - 0.2756*self.f_acent + 0.059035*mur**4                   # Eq 6

        muo *= Fc
        return muo

    visco0 = {"__name__": "Wen (2017)",
              "__doi__": {
                  "autor": "Wen, C., Meng, X., Huber, M.L., Wu, J.",
                  "title": "Measurement and Correlation of the Viscosity of "
                           "1,1,1,2,2,4,5,5,5-Nanofluoro-4-(trifluromethyl)-"
                           "3-pentanone",
                  # TODO: Search final reference
                  "ref": "J. Chem. Eng. Data ",
                  "doi": "10.1021/acs.jced.7b00572"},

              "eq": 1, "omega": 5,

              "ek": 350.84, "sigma": 0.6509,

              "Tref_virial": 350.84,
              "n_virial": [-19.572881, 219.73999, -1015.3226, 2471.01251,
                           -3375.1717, 2491.6597, -787.26086, 14.085455,
                           -0.34664158],
              "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],

              "Tref_res": Tc, "rhoref_res": rhoc,
              "nr": [22.0057],
              "tr": [-0.5],
              "dr": [2/3],
              "nr_num": [231.063],
              "tr_num": [-0.5],
              "dr_num": [2/3],
              "nr_den": [0.423359, -0.122057, 18.4610, -11.1393, 1.67777],
              "tr_den": [0, 0, -1, -1, -1],
              "dr_den": [0, 1, 0, 1, 2]}

    _viscosity = visco0,

    thermo0 = {"__name__": "Perkins (2018)",
               "__doi__": {
                   "autor": "Perkins, R.A., Huber, M.L.",
                   "title": "Measurement and Correlation of the Thermal "
                            "Conducitivity of 1,1,1,2,2,4,5,5,5-Nanofluoro-4-"
                            "(trifluromethyl)-3-pentanone",
                   "ref": "J. Phys. Chem. Ref. Data 47(4) (2018) 043101",
                   # TODO: Search final reference
                   "ref": "J. Chem. Eng. Data ",
                   "doi": "10.1021/acs.jced.8b00132"},

               "eq": 1,

               "Toref": Tc, "koref": 1,
               "no_num": [1.54022e-3, -15.0745e-3, 49.0451e-3, -60.7192e-3,
                          46.2647e-3, -3.16935e-3],
               "to_num": [0, 1, 2, 3, 4, 5],
               "no_den": [-0.211741, 1.16696, -1.15574, 1.0],
               "to_den": [0, 1, 2, 3],

               "Tref_res": Tc, "rhoref_res": rhoc, "kref_res": 1,
               "nr": [-0.0234542, 0.0158544, 0.0418017, -0.0334181, -0.0274745,
                      0.0282766, 0.00932188, -0.0103009, -0.00111766,
                      0.00147852],
               "tr": [0, -1, 0, -1, 0, -1, 0, -1, 0, -1],
               "dr": [1, 1, 2, 2, 3, 3, 4, 4, 5, 5],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.02, "Xio": 2.51e-10,
               "gam0": 0.061, "qd": 0.334e-9, "Tcref": 1.5*Tc}

    _thermal = thermo0,


class Test(TestCase):

    def test_mclinden(self):
        # Table 9, pag J
        st = Novec649(T=250, rhom=5.6)
        self.assertEqual(round(st.P.MPa, 6), 11.459869)
        self.assertEqual(round(st.cvM.JmolK, 3), 277.136)
        self.assertEqual(round(st.cpM.JmolK, 3), 341.656)
        self.assertEqual(round(st.w, 3), 716.250)

        st = Novec649(T=400, rhom=4)
        self.assertEqual(round(st.P.MPa, 7), 2.9272365)
        self.assertEqual(round(st.cvM.JmolK, 3), 308.183)
        self.assertEqual(round(st.cpM.JmolK, 3), 386.271)
        self.assertEqual(round(st.w, 3), 245.654)

        st = Novec649(T=442, rhom=1.92)
        self.assertEqual(round(st.P.MPa, 7), 1.8757729)
        self.assertEqual(round(st.cvM.JmolK, 3), 351.688)
        self.assertEqual(round(st.cpM.JmolK, 1), 45430.0)
        self.assertEqual(round(st.w, 4), 37.0299)

        st = Novec649(T=250, rhom=0.001)
        self.assertEqual(round(st.P.MPa, 10), 0.0020718017)
        self.assertEqual(round(st.cvM.JmolK, 3), 254.272)
        self.assertEqual(round(st.cpM.JmolK, 3), 262.742)
        self.assertEqual(round(st.w, 4), 82.1676)

        st = Novec649(T=450, rhom=4.6)
        self.assertEqual(round(st.P.MPa, 6), 42.305705)
        self.assertEqual(round(st.cvM.JmolK, 3), 319.294)
        self.assertEqual(round(st.cpM.JmolK, 3), 365.304)
        self.assertEqual(round(st.w, 3), 512.603)

    def test_Wen(self):
        # Table 4, pag 6
        self.assertEqual(round(Novec649(T=250, rho=0).mu.muPas, 2), 8.09)
        self.assertEqual(round(Novec649(T=250, rho=0.41).mu.muPas, 2), 8.33)
        self.assertEqual(round(
            Novec649(T=250, rho=1809.77).mu.muPas, 1), 2377.5)
        self.assertEqual(round(Novec649(T=300, rho=0).mu.muPas, 2), 9.77)
        self.assertEqual(round(Novec649(T=300, rho=3.89).mu.muPas, 2), 10.85)
        self.assertEqual(round(
            Novec649(T=300, rho=1701.48).mu.muPas, 1), 1059.7)
        self.assertEqual(round(Novec649(T=350, rho=0).mu.muPas, 2), 11.43)
        self.assertEqual(round(Novec649(T=350, rho=4.42).mu.muPas, 2), 12.65)
        self.assertEqual(round(
            Novec649(T=350, rho=1595.99).mu.muPas, 2), 587.87)

    def test_perkins(self):
        # Table 3, pag 5
        self.assertEqual(round(Novec649(T=300, rho=0).k.WmK, 6), 0.011876)
        self.assertEqual(round(Novec649(T=300, rho=5.5).k.WmK, 6), 0.011813)
        self.assertEqual(round(Novec649(T=300, rho=1673.3).k.WmK, 6), 0.065209)
        self.assertEqual(round(Novec649(T=445, rho=0).k.WmK, 6), 0.022632)
        self.assertEqual(round(Novec649(T=445, rho=685).k.WmK, 6), 0.036507)
