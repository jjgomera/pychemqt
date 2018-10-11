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

from lib.meos import MEoS
from lib import unidades


class oXylene(MEoS):
    """Multiparameter equation of state for o-xylene """
    name = "o-xylene"
    CASNumber = "95-47-6"
    formula = "C8H10"
    synonym = "1,2-dimethylbenzene"
    _refPropName = "OXYLENE"
    _coolPropName = "o-Xylene"
    rhoc = unidades.Density(285)
    Tc = unidades.Temperature(630.259)
    Pc = unidades.Pressure(3737.5, "kPa")
    M = 106.165  # g/mol
    Tt = unidades.Temperature(247.985)
    Tb = unidades.Temperature(417.521)
    f_acent = 0.312
    momentoDipolar = unidades.DipoleMoment(0.63, "Debye")
    id = 42

    Fi1 = {"ao_log": [1, 2.748798],
           "pow": [0, 1],
           "ao_pow": [10.137376, -0.91282993],
           "ao_exp": [4.754892, 6.915052, 15.84813, 10.93886],
           "titao": [225/Tc, 627/Tc, 1726/Tc, 4941/Tc]}

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

        "Tmin": Tt, "Tmax": 700.0, "Pmax": 70000.0, "rhomax": 8.648,
        "Pmin": 0.0228, "rhomin": 8.647,

        "nr1": [0.0036765156, -0.13918171, 0.014104203, 1.5398899, -2.3600925,
                -0.44359159, 0.19596977],
        "d1": [5, 1, 4, 1, 1, 2, 3],
        "t1": [1.0, 0.6, 0.91, 0.3, 0.895, 1.167, 0.435],

        "nr2": [-1.0909408, -0.21890801, 1.1179223, -0.93563815, -0.018102996],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.766, 3.8, 1.31, 3.0, 0.77],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [1.4172368, -0.57134695, -0.081944041, -40.682878],
        "d3": [1, 1, 3, 3],
        "t3": [1.41, 4.8, 1.856, 2.0],
        "alfa3": [1.1723, 1.095, 1.6166, 20.4],
        "beta3": [2.442, 1.342, 3.0, 450.0],
        "gamma3": [1.2655, 0.3959, 0.7789, 1.162],
        "epsilon3": [0.552, 0.728, 0.498, 0.894]}

    eq = zhou,

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.2834, -1.5813, 7.6516, -7.9953, -2.2277],
        "t": [1.0, 1.5, 1.9, 2.4, 6.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.9743, 16.511, -52.934, 87.962, -71.719, 22.569],
        "t": [0.3, 0.96, 1.4, 1.9, 2.4, 3.0]}
    _vapor_Density = {
        "eq": 2,
        "n": [-1.29038, -33.3428, 142.046, -292.211, 293.950, -159.504,
              -88.2170],
        "t": [0.32, 1.14, 1.7, 2.2, 2.8, 3.5, 9.8]}

    visco0 = {"__name__": "Cao (2016)",
              "__doi__": {
                  "autor": "Cao, F.L., Meng, X.Y., Wu, J.T., Vesovic, V.",
                  "title": "Reference Correlation of the Viscosity of "
                           "ortho-Xylene from 273 to 673 K and up to 110 MPa",
                  "ref": "J. Phys. Chem. Ref. Data 45(2) (2016) 023102",
                  "doi": "10.1063/1.4945663"},

              "eq": 1, "omega": 3,
              "collision": [-1.4933, 473.2, -57033],

              "sigma": 1,
              "n_chapman": 0.22225/M**0.5,

              "Tref_res": 630.259, "rhoref_res": 2.6845*M,
              "nr": [-2.05581e-3, 2.65651e-3, 2.38762, 1.77616e-12, 10.4497,
                     -18.2446, 15.9587],
              "dr": [329/30, 329/30, 119/30, 77/3, 71/30, 41/30, 16/15],
              "tr": [-0.5, 0.3, -0.5, 3.9, -0.5, -1.5, -0.5],

              "special": "_vir"}

    def _vir(self, rho, T, fase):
        # The initial density dependence has a different expresion, without muo
        # and other normal method calculation so hardcoded here
        muB = 0
        if rho:
            for i, n in enumerate([13.2814, -10862.4, 1664060]):
                muB += n/T**i
        return muB*rho/self.M

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

               "Toref": 630.259, "koref": 1e-3,
               "no_num": [-0.837488, 12.7856, -37.1925, 63.9548, -4.43443],
               "to_num": [0, 1, 2, 3, 4],
               "no_den": [0.262226, -0.490519, 1],
               "to_den": [0, 1, 2],

               "Tref_res": 630.259, "rhoref_res": 285, "kref_res": 1e-3,
               "nr": [-3.46292e1, 7.57735e1, -6.74378e1, 2.76950e1, -3.74238,
                      4.55879e1, -5.94473e1, 5.50012e1, -2.55522e1, 4.18805],
               "tr": [0, 0, 0, 0, 0, -1, -1, -1, -1, -1],
               "dr": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.02, "Xio": 0.236e-9,
               "gam0": 0.058, "qd": 0.711e-9, "Tcref": 945.4}

    _thermal = thermo0,


class Test(TestCase):

    def test_Cao(self):
        # Table 5, saturation state properties, include basic test for Zhou EoS
        st = oXylene(T=273.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 0.0002)
        self.assertEqual(round(st.Gas.rhoM, 4), 0.0001)
        self.assertEqual(round(st.Liquido.rhoM, 4), 8.4485)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 1107.7)

        st = oXylene(T=333.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 0.0055)
        self.assertEqual(round(st.Gas.rhoM, 4), 0.002)
        self.assertEqual(round(st.Liquido.rhoM, 4), 7.9701)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 499.8)

        st = oXylene(T=393.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 0.0507)
        self.assertEqual(round(st.Gas.rhoM, 4), 0.0160)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 8.50)
        self.assertEqual(round(st.Liquido.rhoM, 4), 7.4658)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 294.5)

        st = oXylene(T=453.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 0.2389)
        self.assertEqual(round(st.Gas.rhoM, 4), 0.0692)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 9.81)
        self.assertEqual(round(st.Liquido.rhoM, 4), 6.9087)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 195.6)

        st = oXylene(T=553.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 1.3884)
        self.assertEqual(round(st.Gas.rhoM, 4), 0.4174)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 12.53)
        self.assertEqual(round(st.Liquido.rhoM, 4), 5.6991)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 103.8)

        # Table 6, Pag 9
        self.assertEqual(round(oXylene(T=300, rhom=0).mu.muPas, 3), 6.670)

        # This point isn't real, it's in two phases region so need force
        # calculation
        st = oXylene(T=300, rhom=0.04)
        mu = st._Viscosity(0.04*st.M, 300, None)
        self.assertEqual(round(mu.muPas, 3), 6.598)

        self.assertEqual(round(
            oXylene(T=300, rhom=8.2369).mu.muPas, 3), 738.286)
        self.assertEqual(round(
            oXylene(T=300, rhom=8.7845).mu.muPas, 3), 1645.436)
        self.assertEqual(round(oXylene(T=400, rhom=0).mu.muPas, 3), 8.658)

        # This point isn't real, it's in two phases region so need force
        # calculation
        st = oXylene(T=400, rhom=0.04)
        mu = st._Viscosity(0.04*st.M, 400, None)
        self.assertEqual(round(mu.muPas, 3), 8.634)

        self.assertEqual(round(
            oXylene(T=400, rhom=7.4060).mu.muPas, 3), 279.954)
        self.assertEqual(round(
            oXylene(T=400, rhom=8.2291).mu.muPas, 3), 595.652)
        self.assertEqual(round(oXylene(T=600, rhom=0).mu.muPas, 3), 12.904)
        self.assertEqual(round(oXylene(T=600, rhom=0.04).mu.muPas, 3), 13.018)
        self.assertEqual(round(
            oXylene(T=600, rhom=7.2408).mu.muPas, 3), 253.530)

    def test_Mylona(self):
        # The critical enchancement use a innacurate ecs viscosity correlation
        # This viscosity with that correlation is fairly diferent of Cao
        # correlation, this is the cause of testing error

        # Table 6, the point with critical enhancement differ
        self.assertEqual(round(oXylene(T=250, rho=0).k.mWmK, 2), 10.06)
        self.assertEqual(round(oXylene(T=300, rho=0).k.mWmK, 2), 13.67)
        self.assertEqual(round(oXylene(T=400, rho=0).k.mWmK, 1), 22.4)
        self.assertEqual(round(oXylene(T=500, rho=0).k.mWmK, 1), 32.0)
        self.assertEqual(round(oXylene(T=600, rho=0).k.mWmK, 1), 41.6)
        self.assertEqual(round(oXylene(T=700, rho=0).k.mWmK, 1), 50.9)
        self.assertEqual(round(oXylene(T=250, P=1e5).k.mWmK, 1), 141.5)
        self.assertEqual(round(oXylene(T=300, P=1e5).k.mWmK, 1), 130.2)
        self.assertEqual(round(oXylene(T=400, P=1e5).k.mWmK, 1), 104.9)
        self.assertEqual(round(oXylene(T=500, P=1e5).k.mWmK, 1), 32.0)
        self.assertEqual(round(oXylene(T=600, P=1e5).k.mWmK, 1), 41.6)
        self.assertEqual(round(oXylene(T=700, P=1e5).k.mWmK, 1), 51.0)

        # Critical enhancement point
        self.assertEqual(round(oXylene(T=635, rho=270).k.mWmK, 2), 101.76)
