#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2016, Juan José Gómez Romera <jjgomera@gmail.com>

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


from lib.meos import MEoS
from lib import unidades


class Ethanol(MEoS):
    """Multiparameter equation of state for ethanol"""
    name = "ethanol"
    CASNumber = "64-17-5"
    formula = "C2H6O"
    synonym = ""
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
           "ao_exp": [1.95988750679, 7.60084166080, 3.89583440622, 4.23238091363],
           "exp": [694, 1549, 2911, 4659],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ethanol of Schroeder (2011).",
        "__doi__": {"autor": "Schroeder, J. A.; Penoncello, S. G.; Schroeder, J. S.",
                    "title": "A Fundamental Equation of State for Ethanol",
                    "ref": "J. Phys. Chem. Ref. Data 43, 043102 (2014)",
                    "doi": "10.1063/1.4895394"},
        "__test__": """
            >>> st=Ethanol(T=300, rhom=18)
            >>> print "%0.1f %0.1f %0.8g %0.8g %0.8g %0.8g" % ( \
                st.T, st.rhoM, st.P.MPa, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            300.0 18.0 65.640781 94.211739 110.34295 1454.836
            >>> st=Ethanol(T=450, rhom=13.2)
            >>> print "%0.1f %0.1f %0.8g %0.8g %0.8g %0.8g" % ( \
                st.T, st.rhoM, st.P.MPa, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            450.0 13.2 2.921546 137.3679 191.28901 588.67833
            >>> st=Ethanol(T=450, rhom=0.5)
            >>> print "%0.1f %0.1f %0.8g %0.8g %0.8g %0.8g" % ( \
                st.T, st.rhoM, st.P.MPa, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            450.0 0.5 1.5540771 97.346133 125.53109 263.57231
            >>> st=Ethanol(T=550, rhom=6)
            >>> print "%0.1f %0.1f %0.8g %0.8g %0.8g %0.8g" % ( \
                st.T, st.rhoM, st.P.MPa, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            550.0 6.0 10.315533 151.00169 440.13187 207.74032
            >>> st=Ethanol(T=514.8, rhom=6)
            >>> print "%0.1f %0.1f %0.8g %0.8g %0.8g %0.8g" % ( \
                st.T, st.rhoM, st.P.MPa, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            514.8 6.0 6.2784176 163.95041 115623.88 159.34583
            """, # Table 30, Pag 38

        "R": 8.314472,
        "cp": Fi1,
        "ref": {"Tref": 273.15, "Pref": 1., "ho": 52811.79, "so": 209.583},

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
        "gamma3": [1.194, 1.986, 1.583, .756, .495, 1.002, 1.077, 1.493, 1.542],
        "epsilon3": [.779, .805, 1.869, .694, 1.312, 2.054, .441, .793, .313],
        "nr4": []}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ethanol of Dillon and Penoncello (2004)",
        "__doi__": {"autor": "Dillon, H.E. and Penoncello, S.G.",
                    "title": "A Fundamental Equation for Calculation of the Thermodynamic Properties of Ethanol",
                    "ref": "Int. J. Thermophys., 25(2):321-335, 2004.",
                    "doi": "10.1023/B:IJOT.0000028470.49774.14"},
        "__test__": """
            >>> st=Ethanol(T=350, P=1e5, eq=1)
            >>> print "%0.1f %0.0f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (st.P.MPa, st.T, st.rho, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            0.1 350 737.85 259.95 1.0363 2.6542 3.1707 964.32
            >>> st=Ethanol(T=650, P=1e5, eq=1)
            >>> print "%0.1f %0.0f %0.4g %0.5g %0.5g %0.5g %0.5g %0.5g" % (st.P.MPa, st.T, st.rho, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            0.1 650 0.8534 1771.7 4.8011 2.3679 2.5512 355.1
            >>> st=Ethanol(T=450, P=1e6, eq=1)
            >>> print "%i %i %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (st.P.MPa, st.T, st.rho, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            1 450 13.736 1270.3 3.4711 2.0217 2.3954 276.5
            >>> st=Ethanol(T=250, P=1e7, eq=1)
            >>> print "%i %i %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (st.P.MPa, st.T, st.rho, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            10 250 831.84 9.4714 0.1625 1.6724 2.0325 1372.8
            >>> st=Ethanol(T=600, P=1e7, eq=1)
            >>> print "%i %i %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (st.P.MPa, st.T, st.rho, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            10 600 129 1474.5 3.5258 2.7347 4.0688 273.66
            >>> st=Ethanol(T=400, P=1e8, eq=1)
            >>> print "%i %i %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (st.P.MPa, st.T, st.rho, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            100 400 784.52 500.06 1.3279 2.8005 3.2736 1348.2
            """, # Table IV, Pag 329

        "R": 8.314472,
        "cp": CP1,
        "ref": {"Tref": 273.15, "Pref": 1., "ho": 45800, "so": 180},
        "Tc": 513.9, "rhoc": 5.991,

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

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ethanol of Sun and Ely (2004).",
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering application: Algorithm and  application to non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223:107-118, 2004.",
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

    eq = helmholtz1, helmholtz2, helmholtz3
    _PR = 0.0043733

    _surface = {"sigma": [0.05], "exp": [0.952]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.91043e1, 0.47263e1, -0.97145e1, 0.41536e1, -0.20739e1],
        "exp": [1, 1.5, 2.0, 2.55, 4.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.11632e2, -0.21866e3, 0.82694e3, -0.13512e4, 0.10517e4,
               -0.31809e3],
        "exp": [0.66, 1.5, 1.9, 2.3, 2.7, 3.1]}
    _vapor_Density = {
        "eq": 3,
        "ao": [0.22543e1, -0.24734e2, 0.48993e2, -0.41689e2, -0.45104e2,
               -0.10732e3],
        "exp": [0.18, 0.44, 0.68, 0.95, 4.0, 10.0]}

    visco0 = {"eq": 1, "omega": 1,
              "__name__": "Kiselev (2005)",
              "__doi__": {"autor": "Kiselev, S. B., Ely, J. F., Abdulagatov, I. M., Huber, M. L.",
                          "title": "Generalized SAFT-DFT/DMT Model for the Thermodynamic, Interfacial, and Transport Properties of Associating Fluids: Application for n-Alkanols",
                          "ref": "Ind. Eng. Chem. Res., 2005, 44 (17), pp 6916–6927",
                          "doi": "10.1021/ie050010e"},

              "ek": 362.6, "sigma": 0.453,
              "Tref": 1., "rhoref": 1.*M,
              "n_chapman": 0,
              "n_ideal": [-1.03116, 3.48379e-2, -6.50264e-6],
              "t_ideal": [0, 1, 2],

              "n_virial": [-0.19572881e2, 0.21973999e3, -0.10153226e4,
                           0.24710125e4, -0.33751717e4, 0.24916597e4,
                           -0.78726086e3, 0.14085455e2, -0.34664158],
              "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],
              "Tref_virial": 362.6, "etaref_virial": 0.0559816,

              "Tref_res": 513.9, "rhoref_res": 5.991*M, "etaref_res": 1000,
              "n_packed": [-3.38264465, 1.27568864e1],
              "t_packed": [0, 0.5],
              "n_poly": [1.31194057e-1, -8.05700894e-2, -3.82240694e-1,
                         1.53811778e-1, 0.0, -1.10578307e-1, -2.37222995e1],
              "t_poly": [0, 0, -1, -1, -2, -2, 0],
              "d_poly": [2, 3, 2, 3, 2, 3, 1],
              "g_poly": [0, 0, 0, 0, 0, 0, -1],
              "c_poly": [0, 0, 0, 0, 0, 0, 1],
              "n_num": [2.37222995e1],
              "t_num": [0],
              "d_num": [1],
              "g_num": [0],
              "c_num": [0],
              "n_den": [1, -1],
              "t_den": [0, 0],
              "d_den": [0, 1],
              "g_den": [1, 0],
              "c_den": [0, 0]}

    _viscosity = visco0,

    thermo0 = {"eq": 1,
               "__name__": "Assael (2013)",
               "__doi__": {"autor": "M. J. Assael, E. A. Sykioti, M. L. Huber, and R. A. Perkins",
                           "title": "Reference Correlation of the Thermal Conductivity of Ethanol from the Triple Point to 600 K and up to 245 MPa",
                           "ref": "J. Phys. Chem. Ref. Data 42, 023102 (2013)",
                           "doi": "10.1063/1.4797368"},
               "__test__": """
                   >>> st=Ethanol(T=300, rho=850)
                   >>> print "%0.5g" % st.k.mWmK
                   209.68
                   >>> st=Ethanol(T=400, rho=2)
                   >>> print "%0.5g" % st.k.mWmK
                   26.108
                   >>> st=Ethanol(T=400, rho=690)
                   >>> print "%0.5g" % st.k.mWmK
                   149.21
                   >>> st=Ethanol(T=500, rho=10)
                   >>> print "%0.5g" % st.k.mWmK
                   39.594
                   >>> st=Ethanol(T=500, rho=10)
                   >>> print "%0.5g" % st.k.mWmK
                   40.755
                   """, # Table 4, Pag 8

               "Tref": 514.71, "kref": 1e-3,
               "no": [-2.09575, 1.99045e1, -5.39640e1, 8.21223e1, -1.98864, -0.495513],
               "co": [0, 1, 2, 3, 4, 5],
               "noden": [0.17223, -0.078273, 1.0],
               "coden": [0, 1, 2],

               "Trefb": 514.71, "rhorefb": 5.93, "krefb": 1.,
               "nb": [.267222E-01, .148279, -.130429, .346232E-01,
                      -.244293E-02, .0, .177166E-01, -.893088E-01,
                      .684664E-01, -.145702E-01, .809189E-03, .0],
               "tb": [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
               "db": [1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6],
               "cb": [0]*12,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
               "Xio": 0.164296e-9, "gam0": 0.05885, "qd": 0.53e-9, "Tcref": 770.85}

    thermo1 = {"eq": 1,
               "__name__": "Kiselev (2005)",
               "__doi__": {"autor": "Kiselev, S. B., Ely, J. F., Abdulagatov, I. M., Huber, M. L.",
                           "title": "Generalized SAFT-DFT/DMT Model for the Thermodynamic, Interfacial, and Transport Properties of Associating Fluids: Application for n-Alkanols",
                           "ref": "Ind. Eng. Chem. Res., 2005, 44 (17), pp 6916–6927",
                           "doi": "10.1021/ie050010e"},

               "Tref": 1, "kref": 1,
               "no": [-10.109e-3],
               "co": [0.6475],
               "noden": [1.0, -7.332e3, -2.68e5],
               "coden": [0, -1, -2],

               "Trefb": 513.9, "rhorefb": 5.991, "krefb": 1.,
               "nb": [1.06917458e-1, -5.95897870e-2, -8.65012441e-2,
                      6.14073818e-2, 2.12220237e-2, -1.00317135e-2, 0, 0, 0, 0],
               "tb": [0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
               "db": [1, 1, 2, 2, 3, 3, 4, 4, 5, 5],
               "cb": [0]*10,

               # TODO: Add critical crossover model from paper
               "critical": 0}

    _thermal = thermo0, thermo1
