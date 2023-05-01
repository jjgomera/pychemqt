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


class R161(MEoS):
    """Multiparameter equation of state for R161"""
    name = "fluoroethane"
    CASNumber = "353-36-6"
    formula = "C2H5F"
    synonym = "R161"
    _refPropName = "R161"
    _coolPropName = "R161"
    rhoc = unidades.Density(302.0010921)
    Tc = unidades.Temperature(375.25)
    Pc = unidades.Pressure(5046.0, "kPa")
    M = 48.0595  # g/mol
    Tt = unidades.Temperature(130.0)
    Tb = unidades.Temperature(235.6)
    f_acent = 0.216
    momentoDipolar = unidades.DipoleMoment(1.9397, "Debye")
    id = 247

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [-6.9170707460, 5.4837900434],
           "ao_exp": [1.08888, 1.80842, 8.72417, 5.67715],
           "titao": [329/Tc, 742/Tc, 1644/Tc, 3922/Tc]}

    Fi2 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [-6.9187, 5.4788],
           "ao_exp": [2.059, 9.253, 6.088],
           "titao": [420/Tc, 1548/Tc, 3882/Tc]}

    CP1 = {"ao": 3.985,
           "ao_exp": [2.077, 9.265, 6.054], "exp": [420, 1548, 3882]}

    qi = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-161 of Qi (2012)",
        "__doi__": {"autor": "Qi, H., Fang, D., Gao, K., Meng, X., Wu, J.",
                    "title": "Compressed Liquid Densities and Helmholtz Energy"
                             " Equation of State for Fluoroethane (R161)",
                    "ref": "Int. J. Thermophys. 37(3) (2016) 55",
                    "doi": "10.1007/s10765-016-2061-1"},

        "R": 8.3144621,
        "Tc": 375.25, "rhoc": 6.2839,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 420.0, "Pmax": 100000.0, "rhomax": 18,

        "nr1": [0.005145283, -0.001882274, 1.884722, -3.1819965, -0.24432415,
                0.27792467],
        "d1": [5, 4, 1, 1, 2, 3],
        "t1": [1, 0.68, 0.32, 0.92, 1.23, 0.846],

        "nr2": [-0.4414064, -0.402065, 0.24171113, -0.16603585, -0.03440867,
                -0.000099185],
        "d2": [1, 3, 2, 2, 7, 5],
        "t2": [4.208, 3.06, 1.85, 4.28, 1.003, 1.12],
        "c2": [2, 2, 1, 2, 1, 1],
        "gamma2": [1]*6,

        "nr3": [1.0146668, -0.03542609, -0.006038245, -0.025437558,
                -0.00515678, 0.006396804],
        "d3": [1, 1, 3, 3, 2, 2],
        "t3": [1.055, 0.8, 4.08, 1.6, 3.85, 0.57],
        "alfa3": [0.96212, 3.2147, 2.6288, 0.8657, 2.3839, 1.7814],
        "beta3": [0.62848, 4.5968, 4.9696, 0.239, 0.788, 7.0874],
        "gamma3": [1.9363, 1.5054, 1.3691, 2.3594, 0.5581, 0.6326],
        "epsilon3": [0.70192, 1.23824, 0.73324, 0.6258, 1.564, 1.4861]}

    wu = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-161 of Wu and "
                    "Zhou (2012)",
        "__doi__": {"autor": "Wu, J., Zhou, Y.",
                    "title": "An Equation of State for Fluoroethane (R161)",
                    "ref": "Int. J. Thermophys. 33(2) (2012) 220-234",
                    "doi": "10.1007/s10765-011-1151-3"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": {"Tref": 273.15, "Pref": 1., "ho": 28559.6, "so": 167.205},

        "Tmin": Tt, "Tmax": 450.0, "Pmax": 5000.0, "rhomax": 20.0,

        "nr1": [1.511, -2.3, -0.457, 0.1683, 0.04133],
        "d1": [1, 1, 2, 3, 4],
        "t1": [0.37, 0.97, 1.14, 0.744, 1.],

        "nr2": [0.62187, -0.0265, -1.03, -0.285, -0.476],
        "d2": [2, 7, 1, 2, 3],
        "t2": [1.26, 1., 1.8, 3., 2.25],
        "c2": [1, 1, 2, 2, 2],
        "gamma2": [1]*5,

        "nr3": [0.82, -0.3532, -0.116, -0.0220583, -1.63148],
        "d3": [11, 1, 3, 3, 3],
        "t3": [1, 1.2, 5.3, 1, 4],
        "alfa3": [0.96, 1.35, 1.26, 1.23, 16.8],
        "beta3": [2.7, 5.2, 3.9, 4.7, 413],
        "gamma3": [0.9, 0.69, 0.67, 0.67, 1.15],
        "epsilon3": [0.683, 0.892, 0.785, 1.33, 0.86]}

    refprop = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-161 of Lemmon "
                    "(2005).",
        "__doi__": {"autor": "Lemmon, E.W.",
                    "title": "preliminary equation, 2005.",
                    "ref": "",
                    "doi": ""},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 400.0, "Pmax": 50000.0, "rhomax": 20.0,

        "nr1": [0.75688, -1.4110, -0.63922, 0.055685, 0.00028395],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [.73357, .67596, .011369, -.56406, -.094362, -.1678, .00034215],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2.0, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 2],
        "gamma2": [1]*7}

    eq = qi, wu, refprop

    _surface = {"sigma": [0.05385], "exp": [1.111]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-8.977955, 25.64713, -80.33162, 89.22478, -34.33593],
        "t": [1.0, 1.5, 1.84, 2.11, 2.47]}
    _liquid_Density = {
        "eq": 1,
        "n": [2.11404, 6.25555, -16.12805, 116.00294, -233.6455, 149.54793],
        "t": [0.46, 1.08, 1.71, 3.56, 4.44, 5.51]}
    _vapor_Density = {
        "eq": 2,
        "n": [-2.26823, 1.39859, -27.93344, 49.91536, -50.51966],
        "t": [0.3, 0.88, 1.46, 2.1, 2.82]}

    visco0 = {"__name__": "Tsolakidou (2017)",
              "__doi__": {
                  "autor": "Tsolakidou, C.M., Assael, M.J., Huber, M.L.,"
                           "Perkins, R.A.",
                  "title": "Correlations for the Viscosity and Thermal "
                           "Conductivity of Ethyl Fluoride (R161)",
                  "ref": "J. Phys. Chem. Ref. Data 46(2) (2017) 023103",
                  "doi": "10.1063/1.4983027"},

              "eq": 1, "omega": 1,

              "n_chapman": 0.021357,
              "ek": 320.39, "sigma": 0.4457,
              "collision": [0.24130, -0.45],

              "Tref_virial": 320.39,
              "n_virial": [-19.572881, 219.73999, -1015.3226, 2471.0125,
                           -3375.1717, 2491.6597, -787.26086, 14.085455,
                           -0.34664158],
              "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],

              "Tref_res": Tc, "rhoref_res": 302.001,
              "nr": [-10.28373, 7.65563, 4.842, 0.42223],
              "tr": [-0.5, -0.5, -2.5, -0.5],
              "dr": [2/3, 5/3, 2/3, 14/3],

              "nr_num": [64.34983, 64.34983],
              "tr_num": [-1.5, -0.5],
              "dr_num": [2/3, 5/3],
              # The Eq 8 in paper have a typo, the denominator of last term
              # has a plus but the correlation work with a minus
              "nr_den": [10.99213, -1],
              "tr_den": [-2, -2],
              "dr_den": [0, 2]}

    _viscosity = (visco0, )

    thermo0 = {"__name__": "Tsolakidou (2017)",
               "__doi__": {
                   "autor": "Tsolakidou, C.M., Assael, M.J., Huber, M.L.,"
                            "Perkins, R.A.",
                   "title": "Correlations for the Viscosity and Thermal "
                            "Conductivity of Ethyl Fluoride (R161)",
                   "ref": "J. Phys. Chem. Ref. Data 46(2) (2017) 023103",
                   "doi": "10.1063/1.4983027"},

               "eq": 1,

               "Toref": Tc, "koref": 1e-3,
               "no_num": [7.96804, -12.5874, -26.3743, 16.9894, 127.545,
                          -32.548],
               "to_num": [0, 1, 2, 3, 4, 5],
               "no_den": [5.406, -18.8331, 24.868, -9.14139, 1],
               "to_den": [0, 1, 2, 3, 4],

               "Tref_res": Tc, "rhoref_res": 302.001, "kref_res": 1e-3,
               "nr": [-8.41553, 7.41456, -39.7744, 44.0586, 106.179, -81.9833,
                      -53.2351, 37.6052, 8.23094, -4.90293],
               "tr": [0, -1, 0, -1, 0, -1, 0, -1, 0, -1],
               "dr": [1, 1, 2, 2, 3, 3, 4, 4, 5, 5],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.02, "Xio": 0.183e-9,
               "gam0": 0.055, "qd": 3.104e-10, "Tcref": 562.88}

    _thermal = (thermo0, )


class Test(TestCase):

    def test_Tsolakidou(self):
        # Table 9, pag 10, saturation state for mEoS testing
        st = R161(T=250, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 0.1880)
        self.assertEqual(round(st.Liquido.rho, 2), 789.54)
        self.assertEqual(round(st.Gas.rho, 2), 4.63)
        self.assertEqual(round(st.Liquido.mu.muPas, 2), 204.34)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 8.15)
        self.assertEqual(round(st.Liquido.k.mWmK, 2), 140.98)
        self.assertEqual(round(st.Gas.k.mWmK, 2), 9.92)

        st = R161(T=275, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 0.4639)
        self.assertEqual(round(st.Liquido.rho, 2), 745.02)
        self.assertEqual(round(st.Gas.rho, 2), 10.96)
        self.assertEqual(round(st.Liquido.mu.muPas, 2), 152.87)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 8.86)
        self.assertEqual(round(st.Liquido.k.mWmK, 2), 125.46)
        self.assertEqual(round(st.Gas.k.mWmK, 2), 13.03)

        st = R161(T=300, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 0.9716)
        self.assertEqual(round(st.Liquido.rho, 2), 693.43)
        self.assertEqual(round(st.Gas.rho, 2), 22.76)
        self.assertEqual(round(st.Liquido.mu.muPas, 2), 115.89)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 9.64)
        self.assertEqual(round(st.Liquido.k.mWmK, 2), 110.40)
        self.assertEqual(round(st.Gas.k.mWmK, 2), 16.56)

        st = R161(T=325, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 1.8072)
        self.assertEqual(round(st.Liquido.rho, 2), 631.73)
        self.assertEqual(round(st.Gas.rho, 2), 43.92)
        self.assertEqual(round(st.Liquido.mu.muPas, 2), 87.59)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 10.68)
        self.assertEqual(round(st.Liquido.k.mWmK, 2), 95.79)
        self.assertEqual(round(st.Gas.k.mWmK, 2), 21.36)

        st = R161(T=350, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 3.0880)
        self.assertEqual(round(st.Liquido.rho, 2), 551.03)
        self.assertEqual(round(st.Gas.rho, 2), 84.51)
        self.assertEqual(round(st.Liquido.mu.muPas, 2), 63.78)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 12.54)
        self.assertEqual(round(st.Liquido.k.mWmK, 2), 81.95)
        self.assertEqual(round(st.Gas.k.mWmK, 2), 30.41)

        # Table 10, pag 10, basic ρ,P,T point for mEoS testing
        st = R161(T=250, P=1e7)
        self.assertEqual(round(st.rho, 1), 803.6)
        self.assertEqual(round(st.mu.muPas, 1), 223.3)
        self.assertEqual(round(st.k.mWmK, 1), 147.9)

        st = R161(T=275, P=1e7)
        self.assertEqual(round(st.rho, 1), 764.0)
        self.assertEqual(round(st.mu.muPas, 1), 169.6)
        self.assertEqual(round(st.k.mWmK, 1), 133.1)

        st = R161(T=300, P=1e7)
        self.assertEqual(round(st.rho, 1), 719.8)
        self.assertEqual(round(st.mu.muPas, 1), 131.8)
        self.assertEqual(round(st.k.mWmK, 1), 118.9)

        st = R161(T=325, P=1e7)
        self.assertEqual(round(st.rho, 1), 670.0)
        self.assertEqual(round(st.mu.muPas, 1), 103.7)
        self.assertEqual(round(st.k.mWmK, 1), 105.4)

        st = R161(T=350, P=1e7)
        self.assertEqual(round(st.rho, 1), 612.6)
        self.assertEqual(round(st.mu.muPas, 1), 81.80)
        self.assertEqual(round(st.k.mWmK, 1), 92.7)

        st = R161(T=250, P=2e7)
        self.assertEqual(round(st.rho, 1), 815.9)
        self.assertEqual(round(st.mu.muPas, 1), 242.0)
        self.assertEqual(round(st.k.mWmK, 1), 154.3)

        st = R161(T=275, P=2e7)
        self.assertEqual(round(st.rho, 1), 780.1)
        self.assertEqual(round(st.mu.muPas, 1), 185.9)
        self.assertEqual(round(st.k.mWmK, 1), 140.2)

        st = R161(T=300, P=2e7)
        self.assertEqual(round(st.rho, 1), 741.5)
        self.assertEqual(round(st.mu.muPas, 1), 147.1)
        self.assertEqual(round(st.k.mWmK, 1), 126.9)

        st = R161(T=325, P=2e7)
        self.assertEqual(round(st.rho, 1), 699.9)
        self.assertEqual(round(st.mu.muPas, 1), 118.9)
        self.assertEqual(round(st.k.mWmK, 1), 114.4)

        # Table 11, pag 11
        st = R161(T=250, rho=0)
        self.assertEqual(round(st.mu.muPas, 3), 8.280)
        self.assertEqual(round(st.k.mWmK, 3), 9.892)

        st = R161(T=250, rho=1)
        self.assertEqual(round(st.mu.muPas, 3), 8.255)
        self.assertEqual(round(st.k.mWmK, 3), 9.884)

        st = R161(T=250, rho=850)
        self.assertEqual(round(st.mu.muPas, 2), 308.22)
        self.assertEqual(round(st.k.mWmK, 2), 175.48)

        st = R161(T=375, rho=0)
        self.assertEqual(round(st.mu.muPas, 3), 12.171)
        self.assertEqual(round(st.k.mWmK, 3), 24.517)

        st = R161(T=375, rho=229)
        self.assertEqual(round(st.mu.muPas, 3), 20.859)
        self.assertEqual(round(st.k.mWmK, 3), 81.297)
