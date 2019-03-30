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


class R1234ze(MEoS):
    """Multiparameter equation of state for R1234ze"""
    name = "trans-1,3,3,3-tetrafluoropropene"
    CASNumber = "29118-24-9"
    formula = "CHF=CHCF3"
    synonym = "R-1234ze"
    _refPropName = "R1234ZE"
    rhoc = unidades.Density(489.238464)
    Tc = unidades.Temperature(382.513)
    Pc = unidades.Pressure(3634.9, "kPa")
    M = 114.0416  # g/mol
    Tt = unidades.Temperature(168.62)
    Tb = unidades.Temperature(254.177)
    f_acent = 0.313
    momentoDipolar = unidades.DipoleMoment(1.27, "Debye")
    id = 671

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [-12.558347537, 8.7912297624],
           "ao_exp": [9.3575, 10.717],
           "titao": [513/Tc, 1972/Tc],
           "ao_hyp": [], "hyp": []}

    Fi2 = {"ao_log": [1, 5.8887],
           "pow": [], "ao_pow": [],
           "ao_exp": [7.0804, 9.3371, 2.5577],
           "titao": [620/Tc, 1570/Tc, 3953/Tc],
           "ao_hyp": [], "hyp": []}

    Fi3 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [-10.8724711, -30.1326538],
           "ao_exp": [6.07536, 9.95795],
           "titao": [289/Tc, 1303/Tc],
           "ao_hyp": [], "hyp": []}

    thol = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R1234ze of Thol (2016)",
        "__doi__": {"autor": "Thol, M. and Lemmon, E.W.",
                    "title": "Equation of State for the Thermodynamic"
                             "Properties of trans-1,3,3,3-Tetrafluoroporpene"
                             "[R-1234ze(E)]",
                    "ref": "Int. J. Thermophys. 37(3) (2016) 28",
                    "doi": "10.1007/s10765-016-2040-6"},

        "R": 8.3144621,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 420.0, "Pmax": 20000.0, "rhomax": 13.26,
        "Pmin": 0.2187, "rhomin": 13.26,

        "nr1": [0.03982797, 1.812227, -2.537512, -0.5333254, 0.1677031],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1.0, 0.223, 0.755, 1.24, 0.44],

        "nr2": [-1.323801, -0.6694654, 0.8072718, -0.7740229, -0.01843846],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2., 2.2, 1.2, 1.5, 0.9],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [1.407916, -0.4237082, -0.2270068, -0.805213, 0.00994318,
                -0.008798793],
        "d3": [1, 1, 3, 3, 2, 1],
        "t3": [1.33, 1.75, 2.11, 1.0, 1.5, 1.0],
        "alfa3": [1.0, 1.61, 1.24, 9.34, 5.78, 3.08],
        "beta3": [1.21, 1.37, 0.98, 171, 47.4, 15.4],
        "gamma3": [0.943, 0.642, 0.59, 1.2, 1.33, 0.64],
        "epsilon3": [0.728, 0.87, 0.855, 0.79, 1.3, 0.71],
        "nr4": []}

    mclinden = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R1234ze of McLinden "
                    "(2010)",
        "__doi__": {"autor": "McLinden, M.O., Thol, M., Lemmon, E.W.",
                    "title": "Thermodynamic Properties of trans-1,3,3,3-"
                             "Tetrafluoropropene [R1234ze(E)]: Measurements "
                             "of Density and Vapor Pressure and a "
                             "Comprehensive Equation of State",
                    "ref": "International Refrigeration and Air Conditioning "
                           "Conference at Purdue, July 12-15, 2010.",
                    "doi": "10.0000_docs.lib.purdue.edu_generic-99DA7EA2C877"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "IIR",
        "Tc": 382.52, "rhoc": 4.29, "M": 114.0415928,

        "Tmin": Tt, "Tmax": 420.0, "Pmax": 20000.0, "rhomax": 13.20,
        "Pmin": 0.23, "rhomin": 13.19,

        "nr1": [0.055563, 1.66927, -2.53408, -0.475075, 0.190055],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1.0, 0.34, 0.91, 1.23, 0.46],

        "nr2": [-1.25154, -0.742195, 0.537902, -0.741246, -0.0355064],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.26, 2.50, 2.0, 2.24, 0.9],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [1.58506, -0.502086, -0.19136, -0.975576],
        "d3": [1, 1, 3, 3],
        "t3": [1.06, 1.79, 3.75, 0.92],
        "alfa3": [1.02, 1.34, 1.08, 6.41],
        "beta3": [1.19, 2.29, 1.15, 131.8],
        "gamma3": [1.14, 0.667, 0.505, 1.22],
        "epsilon3": [0.711, 0.914, 0.694, 0.731]}

    akasaka = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R1234yf of Akasaka "
                    "(2011)",
        "__doi__": {"autor": "Akasaka, R.",
                    "title": "New Fundamental Equations of State with a Common"
                             " Functional Form for 2,3,3,3-Tetrafluoropropene "
                             "(R-1234yf) and trans-1,3,3,3-Tetrafluoropropene "
                             "(R-1234ze(E))",
                    "ref": "Int. J. Thermophys. 32(6) (2011) 1125-1147",
                    "doi": "10.1007/s10765-011-0992-0"},

        "R": 8.314472,
        "cp": Fi3,
        "ref": "IIR",

        "Tmin": 240., "Tmax": 420.0, "Pmax": 15000.0, "rhomax": 13.20,
        "Pmin": 0.23, "rhomin": 13.19,

        "nr1": [8.5579765, -9.4701332, -0.25013623, 0.1378987, 0.012177113],
        "d1": [1, 1, 1, 2, 5],
        "t1": [0.66886, 0.83392, 1.6982, 1.8030, 0.36657],

        "nr2": [-0.14227996, 0.10096648, 0.017504319, -0.017627303,
                -0.014705120, 0.37202269, -0.30138266, -0.092927274,
                0.087051177, 0.01811377, -0.016018424, 0.005380986],
        "d2": [1, 3, 5, 7, 1, 2, 2, 3, 4, 2, 3, 5],
        "t2": [3.8666, 1.0194, 0, 1.1655, 8.3101, 6.1459, 8.3495, 6.0422,
               7.444, 15.433, 21.543, 15.499],
        "c2": [1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3],
        "gamma2": [1]*12}

    eq = thol, mclinden, akasaka

    _surface = {
        "__doi__": {
            "autor": "Tanaka, K., Higashi, Y.",
            "title": "Surface Tensions of trans-1,3,3,3-Tetrafluoropropene "
                     "and trans-1,3,3,3-Tetrafluoropropene + Difluoromethane "
                     "Mixture",
            "ref": "J. Chem. Eng. Japan 46(6) (2013) 371-375",
            "doi": "10.1252/jcej.13we021"},
        "sigma": [0.05681], "exp": [1.23]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.5888, 1.9696, -2.0827, -4.1238],
        "t": [1.0, 1.5, 2.2, 4.6]}
    _liquid_Density = {
        "eq": 1,
        "n": [1.1913, 2.2456, -1.7747, 1.3096],
        "t": [0.27, 0.7, 1.25, 1.9]}
    _vapor_Density = {
        "eq": 2,
        "n": [-1.0308, -5.0422, -11.5, -37.499, -77.945],
        "t": [0.24, 0.72, 2.1, 4.8, 9.5]}

    thermo0 = {"__name__": "Perkins (2011)",
               "__doi__": {
                   "autor": "Perkins, R.A., Huber, M.L.",
                   "title": "Measurement and Correlation of the Thermal "
                            "Conductivity of 2,3,3,3-Tetrafluoroprop-1-ene "
                            "(R1234yf) and trans-1,3,3,3-Tetrafluoropropene "
                            "(R1234ze(E))",
                   "ref": "J. Chem. Eng. Data 56(12) (2011) 4868-4874",
                   "doi": "10.1021/je200811n"},

               "eq": 1,
               "Pc": 3.6363e6,

               "Toref": 382.52, "koref": 1,
               "no": [-0.0103589, 0.0308929, 0.000230348],
               "to": [0, 1, 2],

               "Tref_res": 382.52, "rhoref_res": 489.24, "kref_res": 1.,
               "nr": [-0.0428296, 0.0927099, -0.0702107, 0.0249708,
                      -0.00301838, 0.0434288, -0.0605844, 0.0440187,
                      -0.0155082, 0.0021019],
               "tr": [0, 0, 0, 0, 0, -1, -1, -1, -1, -1],
               "dr": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03, "Xio": 0.194e-9,
               "gam0": 0.0496, "qd": 5.835e-10, "Tcref": 573.78}

    _thermal = thermo0,


class Test(TestCase):

    def test_thol(self):
        # Table 3, pag 15
        st = R1234ze(T=200, rhom=12.6)
        self.assertEqual(round(st.P.MPa, 6), 2.161400)
        self.assertEqual(round(st.cvM.JmolK, 4), 94.9481)
        self.assertEqual(round(st.cpM.JmolK, 3), 139.161)
        self.assertEqual(round(st.w, 3), 982.461)

        st = R1234ze(T=350, rhom=11.4)
        self.assertEqual(round(st.P.MPa, 5), 95.04156)
        self.assertEqual(round(st.cvM.JmolK, 3), 112.495)
        self.assertEqual(round(st.cpM.JmolK, 3), 144.400)
        self.assertEqual(round(st.w, 3), 902.540)

        st = R1234ze(T=383, rhom=4.29)
        self.assertEqual(round(st.P.MPa, 6), 3.670105)
        self.assertEqual(round(st.cvM.JmolK, 3), 151.107)
        self.assertEqual(round(st.cpM.JmolK, 1), 17431.1)
        self.assertEqual(round(st.w, 4), 79.9438)

        st = R1234ze(T=200, rhom=0)
        self.assertEqual(round(st.P.MPa, 6), 0)
        self.assertEqual(round(st.cvM.JmolK, 4), 71.6045)
        self.assertEqual(round(st.cpM.JmolK, 4), 79.9190)
        self.assertEqual(round(st.w, 3), 127.572)

        st = R1234ze(T=360, rhom=1)
        self.assertEqual(round(st.P.MPa, 6), 2.039103)
        self.assertEqual(round(st.cvM.JmolK, 3), 113.342)
        self.assertEqual(round(st.cpM.JmolK, 3), 161.395)
        self.assertEqual(round(st.w, 3), 122.236)

        st = R1234ze(T=420, rhom=8)
        self.assertEqual(round(st.P.MPa, 5), 19.95447)
        self.assertEqual(round(st.cvM.JmolK, 3), 120.029)
        self.assertEqual(round(st.cpM.JmolK, 3), 169.513)
        self.assertEqual(round(st.w, 3), 371.320)

    def test_Perkins(self):
        # Table 2, Pag 4872
        # The used EoS is other version of mclinden referenced in REFPROP
        st = R1234ze(T=250, P=5e4, eq="mclinden")
        # self.assertEqual(round(st.rho, 5), 2.80451)
        self.assertEqual(round(st.k, 7), 0.0098506)

        st = R1234ze(T=300, P=1e5, eq="mclinden")
        # self.assertEqual(round(st.rho, 5), 4.67948)
        self.assertEqual(round(st.k, 6), 0.013934)

        st = R1234ze(T=250, P=2e7, eq="mclinden")
        # self.assertEqual(round(st.rho, 2), 1349.37)
        self.assertEqual(round(st.k, 5), 0.10068)

        st = R1234ze(T=300, P=2e7, eq="mclinden")
        # self.assertEqual(round(st.rho, 2), 1233.82)
        self.assertEqual(round(st.k, 6), 0.085574)
