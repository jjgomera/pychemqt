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

from scipy import exp

from lib import unidades
from lib.meos import MEoS


class EthylBenzene(MEoS):
    """Multiparameter equation of state for ethylbenzene"""
    name = "ethylbenzene"
    CASNumber = "100-41-4"
    formula = "C8H10"
    synonym = ""
    _refPropName = "EBENZENE"
    _coolPropName = "EthylBenzene"
    rhoc = unidades.Density(291.)
    Tc = unidades.Temperature(617.12)
    Pc = unidades.Pressure(3622.4, "kPa")
    M = 106.165  # g/mol
    Tt = unidades.Temperature(178.2)
    Tb = unidades.Temperature(409.314)
    f_acent = 0.305
    momentoDipolar = unidades.DipoleMoment(0.6, "Debye")
    id = 45

    Fi1 = {"ao_log": [1, 4.2557889],
           "pow": [0, 1],
           "ao_pow": [5.70409, -0.52414353],
           "ao_exp": [9.7329909, 11.201832, 25.440749],
           "titao": [585/Tc, 4420/Tc, 1673/Tc]}

    zhou = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ethylbenzene of Zhou et "
                    "al. (2012).",
        "__doi__": {"autor": "Zhou, Y., Lemmon, E.W., and Wu, J.",
                    "title": "Thermodynamic Properties of o-Xylene, m-Xylene, "
                             "p-Xylene, and Ethylbenzene",
                    "ref": "J. Phys. Chem. Ref. Data 41, 023103 (2012).",
                    "doi": "10.1063/1.3703506"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 700.0, "Pmax": 60000.0, "rhomax": 9.124,

        "nr1": [0.0018109418, -0.076824284, 0.041823789, 1.5059649, -2.4122441,
                -0.47788846, 0.18814732],
        "d1": [5, 1, 4, 1, 1, 2, 3],
        "t1": [1, 1, 0.92, 0.27, 0.962, 1.033, 0.513],

        "nr2": [-1.0657412, -0.20797007, 1.1222031, -0.99300799, -0.027300984],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.31, 3.21, 1.26, 2.29, 1.0],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [1.3757894, -0.44477155, -0.07769742, -2.16719],
        "d3": [1, 1, 3, 3],
        "t3": [0.6, 3.6, 2.1, 0.5],
        "alfa3": [1.178, 1.07, 1.775, 15.45],
        "beta3": [2.437, 1.488, 4, 418.6],
        "gamma3": [1.2667, 0.4237, 0.8573, 1.15],
        "epsilon3": [0.5494, 0.7235, 0.493, 0.8566]}

    eq = zhou,

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.8411, 2.5921, -3.502, -2.7613],
        "t": [1.0, 1.5, 2.5, 5.4]}
    _liquid_Density = {
        "eq": 1,
        "n": [3.5146, -3.7537, 5.476, -3.4724, 1.2141],
        "t": [0.43, 0.83, 1.3, 1.9, 3.1]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.2877, -3.6071, -15.878, -53.363, -128.57],
        "t": [0.42, 0.98, 2.48, 5.9, 13.4]}

    visco0 = {"__name__": "Meng (2017)",
              "__doi__": {
                  "autor": "Meng, X.Y., Cao, F.L., Wu, J.T., Vesovic, V.",
                  "title": "Reference Correlation of the Viscosity of "
                           "Ethylbenzene from the Triple Point to 673 K and "
                           "up to 110 MPa",
                  "ref": "J. Phys. Chem. Ref. Data 46(1) (2017) 013101",
                  "doi": "10.1063/1.4973501"},

              "eq": 1, "omega": 3,
              "collision": [-1.4933, 473.2, -57033],

              "sigma": 1,
              "n_chapman": 0.22115/M**0.5,

              "Tref_res": 617.12, "rhoref_res": 2.741016*M,
              "nr": [-0.0376893, 0.168877, 17.9684, 3.57702e-11, 29.996,
                     -8.00082, -25.7468],
              "dr": [209/30, 209/30, 29/30, 731/30, 59/30, 29/30, 455/300],
              "tr": [-0.5, 0.6, -0.5, 2.9, -0.5, -1.5, -0.5],

              "special": "_vir"}

    def _vir(self, rho, T, fase):
        # The initial density dependence has a different expresion, without muo
        # and other normal method calculation so hardcoded here
        muB = 0
        if rho:
            for i, n in enumerate([13.2814, -10862.4, 1664060]):
                muB += n/T**i

        # Special exponential term for residual viscosity, Eq 5
        Ei = [-3.29316e-13, -2.92665e-13, 2.97768e-13, 1.76186e-18]
        ni = [4.6, 11.1, 5.6, 12.4]
        ki = [20.8, 10.6, 19.7, 21.9]
        Tr = T/617.12
        rhor = rho/self.M/2.741016

        # Eq 7
        g = 0
        for E, n, k in zip(Ei, ni, ki):
            g += E*rhor**n/Tr**k

        mur = g*exp(rhor**2)

        return muB*rho/self.M + mur

    _viscosity = visco0,

    thermo0 = {"__name__": "Mylona (2014)",
               "__doi__": {
                   "autor": "Mylona, S.K., Antoniadis, K.D., Assael, M.J., "
                            "Huber, M.L., Perkins, R.A.",
                   "title": "Reference Correlations of the Thermal "
                            "Conductivity of o-Xylene, m-Xylene, p-Xylene, "
                            "and Moderate Pressures",
                   "ref": "J. Phys. Chem. Ref. Data 43(4) (2014) 043104",
                   "doi": "10.1063/1.4901166"},

               "eq": 1,

               "Toref": 617.12, "koref": 1e-3,
               "no_num": [-1.10708, 10.8026, -28.9015, 41.9227, 20.9133,
                          -4.01492],
               "to_num": [0, 1, 2, 3, 4, 5],
               "no_den": [0.259475, -0.343879, 1],
               "to_den": [0, 1, 2],

               "Tref_res": 617.12, "rhoref_res": 291, "kref_res": 1e-3,
               "nr": [-4.97837e1, 1.06739e2, -6.85137e1, 2.26133e1, -2.79455,
                      6.63073e1, -1.46279e2, 1.21439e2, -4.62245e1, 6.58554],
               "tr": [0, 0, 0, 0, 0, -1, -1, -1, -1, -1],
               "dr": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.02, "Xio": 0.235e-9,
               "gam0": 0.056, "qd": 0.706e-9, "Tcref": 925.7}

    _thermal = thermo0,


class Test(TestCase):

    def test_Meng(self):
        # Table 5, saturation state properties, include basic test for Zhou EoS
        st = EthylBenzene(T=273.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 0.0003)
        self.assertEqual(round(st.Gas.rhoM, 4), 0.0001)
        self.assertEqual(round(st.Liquido.rhoM, 4), 8.3287)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 879.4)

        st = EthylBenzene(T=333.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 0.0074)
        self.assertEqual(round(st.Gas.rhoM, 4), 0.0027)
        self.assertEqual(round(st.Liquido.rhoM, 4), 7.8329)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 435.1)

        st = EthylBenzene(T=393.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 0.0643)
        self.assertEqual(round(st.Gas.rhoM, 4), 0.0203)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 8.48)
        self.assertEqual(round(st.Liquido.rhoM, 4), 7.3104)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 263.5)

        st = EthylBenzene(T=453.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 0.2884)
        self.assertEqual(round(st.Gas.rhoM, 4), 0.0842)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 9.79)
        self.assertEqual(round(st.Liquido.rhoM, 4), 6.7317)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 171.0)

        st = EthylBenzene(T=553.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 1.5899)
        self.assertEqual(round(st.Gas.rhoM, 4), 0.4922)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 12.52)
        self.assertEqual(round(st.Liquido.rhoM, 4), 5.4711)
        self.assertEqual(round(st.Liquido.mu.muPas, 2), 85.66)

        # Table 6, Pag 9
        self.assertEqual(round(
            EthylBenzene(T=250, rhom=8.9814).mu.muPas, 3), 2948.109)
        self.assertEqual(round(
            EthylBenzene(T=300, rhom=8.1093).mu.muPas, 3), 616.814)
        self.assertEqual(round(
            EthylBenzene(T=300, rhom=8.4082).mu.muPas, 3), 875.361)
        self.assertEqual(round(
            EthylBenzene(T=300, rhom=8.6762).mu.muPas, 3), 1262.986)
        self.assertEqual(round(EthylBenzene(T=350, rhom=0).mu.muPas, 3), 7.591)
        self.assertEqual(round(EthylBenzene(T=400, rhom=0).mu.muPas, 3), 8.616)
        self.assertEqual(round(
            EthylBenzene(T=400, rhom=7.2481).mu.muPas, 3), 250.283)
        self.assertEqual(round(
            EthylBenzene(T=400, rhom=8.1196).mu.muPas, 3), 509.868)
        self.assertEqual(round(
            EthylBenzene(T=500, rhom=0).mu.muPas, 3), 10.734)
        self.assertEqual(round(
            EthylBenzene(T=500, rhom=0.02).mu.muPas, 3), 10.776)
        self.assertEqual(round(
            EthylBenzene(T=600, rhom=0).mu.muPas, 3), 12.841)
        self.assertEqual(round(
            EthylBenzene(T=600, rhom=6.4831).mu.muPas, 3), 155.940)
        self.assertEqual(round(
            EthylBenzene(T=600, rhom=7.1427).mu.muPas, 3), 229.686)

    def test_Mylona(self):
        # The critical enchancement use a innacurate ecs viscosity correlation
        # This viscosity with that correlation is fairly diferent of Meng
        # correlation, this is the cause of testing error

        # Table 21, the point with critical enhancement differ
        self.assertEqual(round(EthylBenzene(T=200, rho=0).k.mWmK, 2), 3.96)
        self.assertEqual(round(EthylBenzene(T=300, rho=0).k.mWmK, 2), 9.71)
        self.assertEqual(round(EthylBenzene(T=400, rho=0).k.mWmK, 2), 18.39)
        self.assertEqual(round(EthylBenzene(T=500, rho=0).k.mWmK, 2), 29.15)
        self.assertEqual(round(EthylBenzene(T=600, rho=0).k.mWmK, 2), 41.13)
        self.assertEqual(round(EthylBenzene(T=700, rho=0).k.mWmK, 2), 53.83)
        self.assertEqual(round(EthylBenzene(T=200, P=1e5).k.mWmK, 1), 151.0)
        self.assertEqual(round(EthylBenzene(T=300, P=1e5).k.mWmK, 1), 127.3)
        self.assertEqual(round(EthylBenzene(T=400, P=1e5).k.mWmK, 1), 103.8)
        self.assertEqual(round(EthylBenzene(T=500, P=1e5).k.mWmK, 2), 29.19)
        self.assertEqual(round(EthylBenzene(T=600, P=1e5).k.mWmK, 2), 41.24)
        self.assertEqual(round(EthylBenzene(T=700, P=1e5).k.mWmK, 2), 53.98)
        self.assertEqual(round(EthylBenzene(T=200, P=2e7).k.mWmK, 1), 154.4)
        self.assertEqual(round(EthylBenzene(T=300, P=2e7).k.mWmK, 1), 133.0)
        self.assertEqual(round(EthylBenzene(T=400, P=2e7).k.mWmK, 1), 111.8)
        self.assertEqual(round(EthylBenzene(T=500, P=2e7).k.mWmK, 2), 96.94)
        self.assertEqual(round(EthylBenzene(T=300, P=6e7).k.mWmK, 1), 143.0)

        # Critical enhancement point
        self.assertEqual(round(EthylBenzene(T=617, rho=316).k.mWmK, 1), 147.9)
