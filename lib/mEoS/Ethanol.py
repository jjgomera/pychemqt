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


class Ethanol(MEoS):
    """Multiparameter equation of state for ethanol"""
    name = "ethanol"
    CASNumber = "64-17-5"
    formula = "C2H6O"
    synonym = ""
    _refPropName = "ETHANOL"
    _coolPropName = "Ethanol"
    rhoc = unidades.Density(273.1858492)
    Tc = unidades.Temperature(514.71)
    Pc = unidades.Pressure(6268., "kPa")
    M = 46.06844  # g/mol
    Tt = unidades.Temperature(159)
    Tb = unidades.Temperature(351.57)
    f_acent = 0.646
    momentoDipolar = unidades.DipoleMoment(1.6909, "Debye")
    id = 134

    Fi1 = {"ao_log": [1, 3.43069],
           "pow": [0, 1],
           "ao_pow": [-12.7531, 9.39094],
           "ao_exp": [2.14326, 5.09206, 6.60138, 5.70777],
           "titao": [420.4/Tc, 1334/Tc, 1958/Tc, 4420/Tc]}

    CP1 = {"ao": 6.4112,
           "an": [], "pow": [],
           "ao_exp": [1.95988750679, 7.60084166080, 3.89583440622,
                      4.23238091363],
           "exp": [694, 1549, 2911, 4659]}

    schroeder = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ethanol of Schroeder "
                    "(2011).",
        "__doi__": {
            "autor": "Schroeder, J.A.; Penoncello, S.G.; Schroeder, J.S.",
            "title": "A Fundamental Equation of State for Ethanol",
            "ref": "J. Phys. Chem. Ref. Data 43(4) (2014) 043102",
            "doi": "10.1063/1.4895394"},

        # The paper report a diferent value for R, but the program verification
        # table work with this ancient value
        "R": 8.314472,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": 159.0, "Tmax": 650.0, "Pmax": 280000.0, "rhomax": 19.74,
        "Pmin": 0.00000088, "rhomin": 19.731,

        "nr1": [0.58200796e-1, 0.94391227, -0.80941908, 0.55359038,
                -0.14269032e1, 0.13448717],
        "d1": [4, 1, 1, 2, 2, 3],
        "t1": [1, 1.04, 2.72, 1.174, 1.329, 0.195],

        "nr2": [0.42671978, -0.11700261e1, -0.92405872, 0.34891808,
                -0.9132772, 0.22629481e-1, -0.15513423, 0.21055146,
                -0.2199769, -0.65857238e-2],
        "d2": [1, 1, 1, 3, 3, 2, 2, 6, 6, 8],
        "t2": [2.43, 1.274, 4.16, 3.3, 4.177, 2.5, 0.81, 2.02, 1.606, 0.86],
        "c2": [1, 1, 2, 1, 2, 1, 2, 1, 1, 1],
        "gamma2": [1]*10,

        "nr3": [.75564749, .1069411, -.69533844e-1, -.24947395, .27177891e-1,
                -0.9053953e-3, -0.12310953, -0.8977971e-1, -0.39512601],
        "d3": [1, 1, 2, 3, 3, 2, 2, 2, 1],
        "t3": [2.5, 3.72, 1.19, 3.25, 3, 2, 2, 1, 1],
        "alfa3": [1.075, .463, .876, 1.108, .741, 4.032, 2.453, 2.3, 3.143],
        "beta3": [1.207, .0895, .581, .947, 2.356, 27.01, 4.542, 1.287, 3.09],
        "gamma3": [1.194, 1.986, 1.583, 0.756, 0.495, 1.002, 1.077, 1.493,
                   1.542],
        "epsilon3": [.779, .805, 1.869, .694, 1.312, 2.054, .441, .793, .313],
        "nr4": []}

    dillon = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ethanol of Dillon and "
                    "Penoncello (2004)",
        "__doi__": {"autor": "Dillon, H.E., Penoncello, S.G.",
                    "title": "A Fundamental Equation for Calculation of the "
                             "Thermodynamic Properties of Ethanol",
                    "ref": "Int. J. Thermophys., 25(2) (2004) 321-335",
                    "doi": "10.1023/B:IJOT.0000028470.49774.14"},

        "R": 8.314472,
        "Tc": 513.9, "rhoc": 5.991,
        "cp": CP1,

        # Changing reference state cited in paper for first point in sample
        # table to fix h-s values in testing, last decimal added to fix last
        # decimal
        # "ref": {"Tref": 273.15, "Pref": 1., "ho": 45800, "so": 180},
        "ref": {"Tref": 350, "Pref": 100., "ho": 259.9514*M, "so": 1.03631*M},

        "Tmin": 250.0, "Tmax": 650.0, "Pmax": 280000.0, "rhomax": 19.4,
        "Pmin": 0.00000088, "rhomin": 19.4,

        "nr1": [0.114008942201e2, -0.395227128302e2, 0.413063408370e2,
                -0.188892923721e2, 0.472310314140e1, -0.778322827052e-2,
                0.171707850032, -0.153758307602e1, 0.142405508571e1,
                0.132732097050, -0.114231649761, 0.327686088736e-5,
                0.495699527725e-3, -0.701090149558e-4, -0.225019381648e-5],
        "d1": [1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 6, 7, 8, 8],
        "t1": [-0.5, 0, 0.5, 1.5, 2, 5, -0.5, 1, 2, 0, 2.5, 6, 2, 2, 4],

        "nr2": [-0.255406026981, -0.632036870646e-1, -0.314882729522e-1,
                0.256187828185e-1, -0.308694499382e-1, 0.722046283076e-2,
                0.299286406225e-2, 0.972795913095e-3],
        "d2": [1, 3, 3, 6, 7, 8, 2, 7],
        "t2": [5, 3, 7, 5.5, 4, 1, 22, 23],
        "c2": [2, 2, 2, 2, 2, 2, 4, 4],
        "gamma2": [1]*8}

    sun = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ethanol of Sun and Ely "
                    "(2004).",
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering "
                             "application: Algorithm and  application to "
                             "non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223 (2004) 107-118",
                    "doi": "10.1016/j.fluid.2004.06.028"},

        "R": 8.314472,
        "cp": CP1,
        "ref": {"name": "CUSTOM",
                "Tref": 273.15, "Pref": 1., "ho": 45800, "so": 180},

        "Tmin": Tt, "Tmax": 650.0, "Pmax": 280000.0, "rhomax": 19.6,
        "Pmin": 0.00000064, "rhomin": 19.55,

        "nr1":  [-2.95455387, 1.95055493, -1.31612955, -1.47547651e-2,
                 1.39251945e-4, 5.04178939e-1],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [2.52310166e-1, 1.97074652, 8.73146115e-1, 4.27767205e-2,
                9.68966545e-2, -8.39632113e-1, -7.71828521e-2, 1.63430744e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = schroeder, dillon, sun
    _PR = 0.0043733

    _surface = {"sigma": [0.05], "exp": [0.952]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.91043e1, 0.47263e1, -0.97145e1, 0.41536e1, -0.20739e1],
        "t": [1, 1.5, 2.0, 2.55, 4.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.11632e2, -0.21866e3, 0.82694e3, -0.13512e4, 0.10517e4,
              -0.31809e3],
        "t": [0.66, 1.5, 1.9, 2.3, 2.7, 3.1]}
    _vapor_Density = {
        "eq": 2,
        "n": [0.22543e1, -0.24734e2, 0.48993e2, -0.41689e2, -0.45104e2,
              -0.10732e3],
        "t": [0.18, 0.44, 0.68, 0.95, 4.0, 10.0]}

    visco0 = {"__name__": "Kiselev (2005)",
              "__doi__": {
                  "autor": "Kiselev, S. B., Ely, J. F., Abdulagatov, I. M., "
                           "Huber, M. L.",
                  "title": "Generalized SAFT-DFT/DMT Model for the "
                           "Thermodynamic, Interfacial, and Transport "
                           "Properties of Associating Fluids: Application for "
                           "n-Alkanols",
                  "ref": "Ind. Eng. Chem. Res. 44(17) (2005) 6916-6927",
                  "doi": "10.1021/ie050010e"},

              "eq": 1, "omega": 0,

              "ek": 362.6, "sigma": 0.453,
              "no": [-1.03116, 3.48379e-2, -6.50264e-6],
              "to": [0, 1, 2],

              "Tref_virial": 362.6,
              "n_virial": [-19.572881, 219.73999, -1015.3226, 2471.0125,
                           -3375.1717, 2491.6597, -787.26086, 14.085455,
                           -0.34664158],
              "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],

              "Tref_res": 513.9, "rhoref_res": 5.991*M, "muref_res": 1e3,
              "nr": [0.131194057, -0.0805700894, -0.382240694, 0.153811778,
                     0.0, -0.110578307],
              "tr": [0, 0, 1, 1, 2, 2],
              "dr": [2, 3, 2, 3, 2, 3],

              "CPf": 23.7222995e3,
              "CPg1": -3.38264465,
              "CPgi": [12.7568864/-3.38264465],
              "CPti": [-0.5]}

    _viscosity = visco0,

    thermo0 = {"__name__": "Assael (2013)",
               "__doi__": {
                   "autor": "Assael, M.J., Sykioti, E.A., Huber, M.L., "
                            "Perkins, R.A.",
                   "title": "Reference Correlation of the Thermal Conductivity"
                            " of Ethanol from the Triple Point to 600 K and "
                            "up to 245 MPa",
                   "ref": "J. Phys. Chem. Ref. Data 42(2) (2013) 023102",
                   "doi": "10.1063/1.4797368"},

               "eq": 1,

               "Toref": 514.71, "koref": 1e-3,
               "no_num": [-2.09575, 19.9045, -53.964, 82.1223, -1.98864,
                          -0.495513],
               "to_num": [0, 1, 2, 3, 4, 5],
               "no_den": [0.17223, -0.078273, 1.0],
               "to_den": [0, 1, 2],

               "Tref_res": 514.71, "rhoref_res": 273.186, "kref_res": 1.,
               "nr": [.267222E-01, .148279, -.130429, .346232E-01,
                      -.244293E-02, .0, .177166E-01, -.893088E-01,
                      .684664E-01, -.145702E-01, .809189E-03, .0],
               "tr": [0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1],
               "dr": [1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.02, "Xio": 0.164296e-9,
               "gam0": 0.05885, "qd": 0.53e-9, "Tcref": 770.85}

    thermo1 = {"__name__": "Kiselev (2005)",
               "__doi__": {
                   "autor": "Kiselev, S. B., Ely, J. F., Abdulagatov, I. M., "
                            "Huber, M. L.",
                   "title": "Generalized SAFT-DFT/DMT Model for the "
                            "Thermodynamic, Interfacial, and Transport "
                            "Properties of Associating Fluids: Application for"
                            " n-Alkanols",
                   "ref": "Ind. Eng. Chem. Res. 44(17) (2005) 6916-6927",
                   "doi": "10.1021/ie050010e"},

               "eq": 1,

               "Toref": 1, "koref": 1,
               "no_num": [-10.109e-3],
               "to_num": [0.6475],
               "no_den": [1.0, -7332, -2.68e5],
               "to_den": [0, -1, -2],

               "Tref_res": 514.45, "rhoref_res": 5.988*M, "kref_res": 1.,
               "nr": [1.06917458e-1, -5.95897870e-2, -8.65012441e-2,
                      6.14073818e-2, 2.12220237e-2, -1.00317135e-2],
               "tr": [0, -1, 0, -1, 0, -1],
               "dr": [1, 1, 2, 2, 3, 3],

               # TODO: Add critical crossover model from paper
               "critical": 0}

    _thermal = thermo0, thermo1


class Test(TestCase):

    def test_schroeder(self):
        # Table 30, Pag 38
        st = Ethanol(T=300, rhom=18)
        self.assertEqual(round(st.P.MPa, 6), 65.640781)
        self.assertEqual(round(st.cvM.JmolK, 6), 94.211739)
        self.assertEqual(round(st.cpM.JmolK, 5), 110.34295)
        self.assertEqual(round(st.w, 4), 1454.8360)

        st = Ethanol(T=450, rhom=13.2)
        self.assertEqual(round(st.P.MPa, 7), 2.9215460)
        self.assertEqual(round(st.cvM.JmolK, 5), 137.36790)
        self.assertEqual(round(st.cpM.JmolK, 5), 191.28901)
        self.assertEqual(round(st.w, 5), 588.67833)

        st = Ethanol(T=450, rhom=0.5)
        self.assertEqual(round(st.P.MPa, 7), 1.5540771)
        self.assertEqual(round(st.cvM.JmolK, 6), 97.346133)
        self.assertEqual(round(st.cpM.JmolK, 5), 125.53109)
        self.assertEqual(round(st.w, 5), 263.57231)

        st = Ethanol(T=550, rhom=6)
        self.assertEqual(round(st.P.MPa, 6), 10.315532)
        self.assertEqual(round(st.cvM.JmolK, 5), 151.00169)
        self.assertEqual(round(st.cpM.JmolK, 5), 440.13187)
        self.assertEqual(round(st.w, 5), 207.74032)

        st = Ethanol(T=514.8, rhom=6)
        self.assertEqual(round(st.P.MPa, 7), 6.2784176)
        self.assertEqual(round(st.cvM.JmolK, 5), 163.95041)
        self.assertEqual(round(st.cpM.JmolK, 2), 115623.88)
        self.assertEqual(round(st.w, 5), 159.34583)

        # Enthalpy reference state
        st = Ethanol(T=273.15, x=0)
        self.assertEqual(round(st.h.kJkg, 2), 200)
        self.assertEqual(round(st.s.kJkgK, 4), 1)

    def test_dillon(self):
        # Table IV, Pag 329
        st = Ethanol(T=350, P=1e5, eq="dillon")
        self.assertEqual(round(st.rho, 2), 737.85)
        self.assertEqual(round(st.h.kJkg, 2), 259.95)
        self.assertEqual(round(st.s.kJkgK, 4), 1.0363)
        self.assertEqual(round(st.cv.kJkgK, 4), 2.6542)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.1707)
        self.assertEqual(round(st.w, 2), 964.32)

        st = Ethanol(P=1e5, x=0.5, eq="dillon")
        self.assertEqual(round(st.T, 2), 351.05)
        self.assertEqual(round(st.Liquido.rho, 2), 736.78)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 263.30)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.0459)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 2.6614)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 3.1818)
        self.assertEqual(round(st.Liquido.w, 2), 960.72)
        self.assertEqual(round(st.Gas.rho, 4), 1.6269)
        self.assertEqual(round(st.Gas.h.kJkg, 1), 1113.8)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 3.4686)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 1.5894)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 1.8044)
        self.assertEqual(round(st.Gas.w, 2), 260.21)

        st = Ethanol(P=1e6, x=0.5, eq="dillon")
        self.assertEqual(round(st.T, 2), 423.85)
        self.assertEqual(round(st.Liquido.rho, 2), 647.20)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 521.56)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.7090)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 2.9808)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 3.9104)
        self.assertEqual(round(st.Liquido.w, 2), 694.41)
        self.assertEqual(round(st.Gas.rho, 3), 15.246)
        self.assertEqual(round(st.Gas.h.kJkg, 1), 1207.2)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 3.3267)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 1.9648)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 2.4276)
        self.assertEqual(round(st.Gas.w, 2), 260.65)

        st = Ethanol(T=550, P=1e6, eq="dillon")
        self.assertEqual(round(st.rho, 3), 10.439)
        self.assertEqual(round(st.h.kJkg, 1), 1510.2)
        self.assertEqual(round(st.s.kJkgK, 4), 3.9523)
        self.assertEqual(round(st.cv.kJkgK, 4), 2.1985)
        self.assertEqual(round(st.cp.kJkgK, 4), 2.4426)
        self.assertEqual(round(st.w, 2), 320.26)

        st = Ethanol(T=250, P=1e7, eq="dillon")
        self.assertEqual(round(st.rho, 2), 831.84)
        self.assertEqual(round(st.h.kJkg, 4), 9.4714)
        self.assertEqual(round(st.s.kJkgK, 4), 0.1625)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.6724)
        self.assertEqual(round(st.cp.kJkgK, 4), 2.0325)
        self.assertEqual(round(st.w, 1), 1372.8)

        st = Ethanol(T=600, P=1e8, eq="dillon")
        self.assertEqual(round(st.rho, 2), 626.59)
        self.assertEqual(round(st.h.kJkg, 1), 1204.7)
        self.assertEqual(round(st.s.kJkgK, 4), 2.7521)
        self.assertEqual(round(st.cv.kJkgK, 4), 2.7296)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.5932)
        self.assertEqual(round(st.w, 1), 1015.1)

    def test_assael(self):
        # Table 4, Pag 8
        self.assertEqual(round(Ethanol(T=300, rho=850).k.mWmK, 2), 209.68)
        self.assertEqual(round(Ethanol(T=400, rho=2).k.mWmK, 3), 26.108)
        self.assertEqual(round(Ethanol(T=400, rho=690).k.mWmK, 2), 149.21)
        self.assertEqual(round(Ethanol(T=500, rho=10).k.mWmK, 3), 39.594)
