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


class pXylene(MEoS):
    """Multiparameter equation of state for p-xylene"""
    name = "p-xylene"
    CASNumber = "106-42-3"
    formula = "C8H10"
    synonym = "1,4-dimethylbenzene"
    _refPropName = "PXYLENE"
    _coolPropName = "p-Xylene"
    rhoc = unidades.Density(286)
    Tc = unidades.Temperature(616.168)
    Pc = unidades.Pressure(3531.5, "kPa")
    M = 106.165  # g/mol
    Tt = unidades.Temperature(286.4)
    Tb = unidades.Temperature(411.47)
    f_acent = 0.324
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 44

    Fi1 = {"ao_log": [1, 4.2430504],
           "pow": [0, 1],
           "ao_pow": [5.9815241, -0.52477835],
           "ao_exp": [5.2291378, 19.549862, 16.656178, 5.9390291],
           "titao": [414/Tc, 1256/Tc, 2649/Tc, 6681/Tc]}

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

        "Tmin": Tt, "Tmax": 700.0, "Pmax": 200000.0, "rhomax": 8.166,

        "nr1": [0.0010786811, -0.103161822, 0.0421544125, 1.47865376, -2.4266,
                -0.46575193, 0.190290995],
        "d1": [5, 1, 4, 1, 1, 2, 3],
        "t1": [1.0, 0.83, 0.83, 0.281, 0.932, 1.1, 0.443],

        "nr2": [-1.06376565, -0.209934069, 1.25159879, -0.951328356,
                -0.0269980032],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.62, 2.5, 1.2, 3.0, 0.778],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [1.37103180, -0.494160616, -0.0724317468, -3.69464746],
        "d3": [1, 1, 3, 3],
        "t3": [1.13, 4.5, 2.2, 2.0],
        "alfa3": [1.179, 1.065, 1.764, 13.675],
        "beta3": [2.445, 1.483, 4.971, 413.0],
        "gamma3": [1.267, 0.4242, 0.864, 1.1465],
        "epsilon3": [0.54944, 0.7234, 0.4926, 0.8459]}

    eq = zhou,

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.7221, 1.5789, -13.035, 18.453, -11.345],
        "t": [1.0, 1.5, 3.8, 4.6, 5.5]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.1783, 3.4488, -2.3906, 1.5933],
        "t": [0.15, 0.5, 0.9, 1.3]}
    _vapor_Density = {
        "eq": 2,
        "n": [-6.17784, -0.38825, -19.0575, -541.124, 1251.55, -920.22],
        "t": [0.653, 0.17, 2.6, 7.8, 8.9, 10.]}

    visco0 = {"__name__": "Balogun (2015)",
              "__doi__": {
                  "autor": "Balogun, B., Riesco, N., Vesovic, V.",
                  "title": "Reference Correlation of the Viscosity of "
                           "para-Xylene from the Triple Point to 673K and up "
                           "to 110 MPa",
                  "ref": "J. Phys. Chem. Ref. Data 44(1) (2015) 013103",
                  "doi": "10.1063/1.4908048"},

              "eq": 1, "omega": 3,
              "collision": [-1.4933, 473.2, -57033],

              "sigma": 1,
              "n_chapman": 0.22005/M**0.5,

              "Tref_res": 616.168, "rhoref_res": 2.69392*M,
              "nr": [122.919, -282.329, 279.348, -146.776, 28.361, -0.004585,
                     15.337, -0.0004382, 0.00002307],
              "dr": [13/6, 8/3, 11/3, 14/3, 17/3, 35/3, 13/6, 35/3, 47/3],
              "tr": [0, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5],

              "special": "_vir"}

    def _vir(self, rho, T, fase):
        # The initial density dependence has a different expresion, without muo
        # and other normal method calculation so hardcoded here
        muB = 0
        if rho:
            for i, n in enumerate([13.2814, -10862.4, 1664060]):
                muB += n/T**i
        return muB*rho/self.M

    visco1 = {"__name__": "Meng (2017)",
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
              "n_chapman": 0.22005/M**0.5,

              "Tref_res": 616.168, "rhoref_res": 2.69392*M,
              "nr": [-0.0258436, 0.0323284, 15.4854, 4.30133e-11, 13.8162,
                     -19.9037, -0.000112615],
              "dr": [269/30, 269/30, 19/15, 145/6, 43/15, 28/15, 35/3],
              "tr": [-0.5, 0.2, -0.5, 2.5, -0.5, -1.5, -0.5],

              "special": "_vir"}

    _viscosity = visco0, visco1

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

               "Toref": 616.168, "koref": 1e-3,
               "no_num": [-3.88568, 29.4648, -81.5299, 77.1534, 7.55487,
                          -3.8897, 0.406892],
               "to_num": [0, 1, 2, 3, 4, 5, 6],
               "no_den": [0.00404188, -0.424893, 1],
               "to_den": [0, 1, 2],

               "Tref_res": 616.168, "rhoref_res": 286, "kref_res": 1e-3,
               "nr": [-1.01022e2, 2.24828e2, -1.59100e2, 4.99490e1, -5.62422,
                      1.07531e2, -2.05499e2, 1.50348e2, -5.02584e1, 6.44051],
               "tr": [0, 0, 0, 0, 0, -1, -1, -1, -1, -1],
               "dr": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.02, "Xio": 0.235e-9,
               "gam0": 0.056, "qd": 0.71e-9, "Tcref": 924.3}

    _thermal = thermo0,


class Test(TestCase):

    def test_balogun(self):
        # Table 7, saturation state properties, include basic test for Zhou EoS
        st = pXylene(T=293.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 0.0009)
        self.assertEqual(round(st.Gas.rhoM, 4), 0.0004)
        self.assertEqual(round(st.Liquido.rhoM, 4), 8.1100)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 644.2)

        st = pXylene(T=323.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 0.0043)
        self.assertEqual(round(st.Gas.rhoM, 4), 0.0016)
        self.assertEqual(round(st.Liquido.rhoM, 4), 7.8637)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 458.5)

        st = pXylene(T=353.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 0.0156)
        self.assertEqual(round(st.Gas.rhoM, 4), 0.0054)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 7.59)
        self.assertEqual(round(st.Liquido.rhoM, 4), 7.6118)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 345.4)

        st = pXylene(T=453.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 0.2746)
        self.assertEqual(round(st.Gas.rhoM, 4), 0.0803)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 9.53)
        self.assertEqual(round(st.Liquido.rhoM, 4), 6.6834)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 164.5)

        st = pXylene(T=553.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 1.5487)
        self.assertEqual(round(st.Gas.rhoM, 4), 0.4883)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 12.30)
        self.assertEqual(round(st.Liquido.rhoM, 4), 5.3990)
        self.assertEqual(round(st.Liquido.mu.muPas, 2), 84.23)

        # Table 8, Pag 9
        self.assertEqual(round(pXylene(T=300, rhom=0).mu.muPas, 3), 6.604)

        # This point isn't real, it's in two phases region so need force
        # calculation
        st = pXylene(T=300, rhom=0.0490)
        mu = st._Viscosity(0.0490*st.M, 300, None)
        self.assertEqual(round(mu.muPas, 3), 6.405)

        self.assertEqual(round(
            pXylene(T=300, rhom=8.0548).mu.muPas, 3), 593.272)
        self.assertEqual(round(
            pXylene(T=300, rhom=8.6309).mu.muPas, 3), 1266.337)
        self.assertEqual(round(pXylene(T=400, rhom=0).mu.muPas, 3), 8.573)
        self.assertEqual(round(
            pXylene(T=400, rhom=7.1995).mu.muPas, 3), 239.202)
        self.assertEqual(round(
            pXylene(T=400, rhom=8.0735).mu.muPas, 3), 484.512)
        self.assertEqual(round(pXylene(T=600, rhom=0).mu.muPas, 3), 12.777)
        self.assertEqual(round(pXylene(
            T=600, rhom=7.0985).mu.muPas, 3), 209.151)

    def test_meng(self):
        # Table 8, Pag 9
        self.assertEqual(round(
            pXylene(T=300, rhom=0, visco=1).mu.muPas, 3), 6.604)

        # This point isn't real, it's in two phases region so need force
        # calculation
        st = pXylene(T=300, rhom=0.0490, visco=1)
        mu = st._Viscosity(0.0490*st.M, 300, None)
        self.assertEqual(round(mu.muPas, 3), 6.450)

        self.assertEqual(round(
            pXylene(T=300, rhom=8.0548, visco=1).mu.muPas, 3), 593.513)
        self.assertEqual(round(
            pXylene(T=300, rhom=8.6309, visco=1).mu.muPas, 3), 1257.494)
        self.assertEqual(round(
            pXylene(T=400, rhom=0, visco=1).mu.muPas, 3), 8.573)
        self.assertEqual(round(
            pXylene(T=400, rhom=7.1995, visco=1).mu.muPas, 3), 239.081)
        self.assertEqual(round(
            pXylene(T=400, rhom=8.0735, visco=1).mu.muPas, 3), 488.777)
        self.assertEqual(round(
            pXylene(T=600, rhom=0, visco=1).mu.muPas, 3), 12.777)
        self.assertEqual(round(
            pXylene(T=600, rhom=7.0985, visco=1).mu.muPas, 3), 199.160)

    def test_Mylona(self):
        # The critical enchancement use a innacurate ecs viscosity correlation
        # This viscosity with that correlation is fairly diferent of Balogun
        # or Meng correlations, this is the cause of testing error

        # Table 16, the point with critical enhancement differ
        self.assertEqual(round(pXylene(T=300, rho=0).k.mWmK, 2), 10.57)
        self.assertEqual(round(pXylene(T=400, rho=0).k.mWmK, 2), 19.48)
        self.assertEqual(round(pXylene(T=500, rho=0).k.mWmK, 2), 30.17)
        self.assertEqual(round(pXylene(T=600, rho=0).k.mWmK, 2), 41.73)
        self.assertEqual(round(pXylene(T=700, rho=0).k.mWmK, 2), 53.68)
        self.assertEqual(round(pXylene(T=300, P=1e5).k.mWmK, 1), 126.1)
        self.assertEqual(round(pXylene(T=400, P=1e5).k.mWmK, 1), 102.2)
        self.assertEqual(round(pXylene(T=500, P=1e5).k.mWmK, 2), 30.06)
        self.assertEqual(round(pXylene(T=600, P=1e5).k.mWmK, 2), 41.76)
        self.assertEqual(round(pXylene(T=700, P=1e5).k.mWmK, 2), 53.82)
        self.assertEqual(round(pXylene(T=700, P=2e8).k.mWmK, 1), 143.2)

        # Critical enhancement point
        self.assertEqual(round(pXylene(T=620, rho=287).k.mWmK, 2), 98.04)
