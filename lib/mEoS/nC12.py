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


class nC12(MEoS):
    """Multiparameter equation of state for n-dodecane"""
    name = "dodecane"
    CASNumber = "112-40-3"
    formula = "CH3-(CH2)10-CH3"
    synonym = ""
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
           "an": [], "pow": [],
           "ao_exp": [37.776, 29.369, 12.461, 7.7733],
           "exp": [1280, 2399, 5700, 13869],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for dodecane of Lemmon (2004).",
        "__doi__": {"autor": "Lemmon, E.W. and Huber, M.L.",
                    "title": "Thermodynamic Properties of n-Dodecane",
                    "ref": "Energy & Fuels, 18:960-967, 2004.",
                    "doi": "10.1021_ef0341062"},
        "__test__": """
            >>> st=nC12(T=300, rho=0)
            >>> print "%0.1f %0.1f %0.1f %0.7g %0.7g %0.7g" % (st.T, st.P.MPa, st.rhoM, st.cvM.JmolK, st.cpM.JmolK, st.w)
            300.0 0.0 0.0 271.3952 279.7096 122.8511
            >>> st=nC12(T=300, rhom=4.4)
            >>> print "%0.1f %0.7g %0.1f %0.7g %0.7g %0.7g" % (st.T, st.P.MPa, st.rhoM, st.cvM.JmolK, st.cpM.JmolK, st.w)
            300.0 7.161946 4.4 311.1937 377.1307 1316.614
            >>> st=nC12(T=600, rhom=0.2)
            >>> print "%0.1f %0.7g %0.1f %0.7g %0.7g %0.7g" % (st.T, st.P.MPa, st.rhoM, st.cvM.JmolK, st.cpM.JmolK, st.w)
            600.0 0.7228018 0.2 488.1173 526.2322 123.5851
            >>> st=nC12(T=600, rhom=4.0)
            >>> print "%0.1f %0.7g %0.1f %0.7g %0.7g %0.7g" % (st.T, st.P.MPa, st.rhoM, st.cvM.JmolK, st.cpM.JmolK, st.w)
            600.0 120.2532 4.0 499.9857 539.3864 1255.247
            >>> st=nC12(T=658.2, rhom=1.33)
            >>> print "%0.1f %0.7g %0.2f %0.7g %0.7g %0.7g" % (st.T, st.P.MPa, st.rhoM, st.cvM.JmolK, st.cpM.JmolK, st.w)
            658.2 1.820055 1.33 547.6852 128753.9 49.76424
            """, # Table 5, Pag 967

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 700., "Pmax": 700000.0, "rhomax": 4.53,
        "Pmin": 0.00063, "rhomin": 4.53,

        "nr1": [1.38031, -2.85352, 0.288897, -0.165993, 0.0923993, 0.000282772],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.32, 1.23, 1.5, 1.4, 0.07, 0.8],

        "nr2": [0.956627, 0.0353076, -0.445008, -0.118911, -0.0366475, 0.0184223],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [2.16, 1.1, 4.1, 5.6, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = helmholtz1,

    _surface = {"sigma": [0.0154, 0.0480],
                "exp": [4.18, 1.17]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.94217e1, -0.41890e1, 0.54999e1, -0.67789e1, -0.17161e1],
        "exp": [1.0, 1.5, 1.359, 3.56, 9.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.92236, 0.92047, 0.55713e1, -0.92253e1, 0.51763e1],
        "exp": [0.21, 0.49, 1.08, 1.49, 1.9]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.17859e1, -0.75436e1, -0.22848e2, -0.81355e2, 0.92283e2, -0.21725e3],
        "exp": [0.298, 0.91, 2.8, 6., 9., 11.]}

    visco0 = {"eq": 1, "omega": 1,
              "collision": [0.382987, -0.561050, 0.313962e-1],
              "__name__": "Huber (2004)",
              "__doi__": {"autor": "Huber, M.L., Laesecke, A. and Perkins, R.A.",
                          "title": "Transport Properties of n-Dodecane",
                          "ref": "Energy Fuels, 2004, 18 (4), pp 968–975.",
                          "doi": "10.1021/ef034109e"},
              "__test__": """
                  >>> st=nC12(T=300, rhom=4.4115)
                  >>> print "%0.1f" % st.mu.muPas
                  1484.8
                  >>> st=nC12(T=500, rhom=3.4447)
                  >>> print "%0.1f" % st.mu.muPas
                  183.76
                  """, # Pag 972

              "ek": 522.592, "sigma": 0.735639,
              "Tref": 1., "rhoref": 1.*M,
              "n_chapman": 0.2787353/M**0.5,
              "n_virial": [-0.19572881e2, 0.21973999e3, -0.10153226e4,
                           0.24710125e4, -0.33751717e4, 0.24916597e4,
                           -0.78726086e3, 0.14085455e2, -0.34664158],
              "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],
              "Tref_virial": 522.592, "etaref_virial": 0.2397238,

              "Tref_res": 658.1, "rhoref_res": 1.33*M, "etaref_res": 1000,
              "n_packed": [0.232661e1, 0.223089e1],
              "t_packed": [0, 0.5],
              "n_poly": [-0.471703e-1, 0.827816e-2, 0.298429e-1, -0.134156e-1,
                         -0.503109],
              "t_poly": [-1, -1, -2, -2, 0],
              "d_poly": [2, 3, 2, 3, 1],
              "g_poly": [0, 0, 0, 0, -1],
              "c_poly": [0, 0, 0, 0, 1],
              "n_num": [0.503109],
              "t_num": [0],
              "d_num": [1],
              "g_num": [0],
              "c_num": [0],
              "n_den": [1, -1],
              "t_den": [0, 0],
              "d_den": [0, 1],
              "g_den": [1, 0],
              "c_den": [1, 0]}

    _viscosity = visco0,

    thermo0 = {"eq": 1,
               "__name__": "Huber (2004)",
               "__doi__": {"autor": "Huber, M.L., Laesecke, A. and Perkins, R.A.",
                           "title": "Transport Properties of n-Dodecane",
                           "ref": "Energy Fuels, 2004, 18 (4), pp 968–975.",
                           "doi": "10.1021/ef034109e"},
               "__test__": """
                   >>> st=nC12(T=300, rhom=4.4115)
                   >>> print "%0.2f" % st.k.mWmK
                   138.29
                   >>> st=nC12(T=500, rhom=3.4447)
                   >>> print "%0.2f" % st.k.mWmK
                   93.84
                   >>> st=nC12(T=660, rhom=1.50098)
                   >>> print "%0.3f" % st.k.mWmK
                   90.346
                   """, # Pag 974

               "Tref": 658.1, "kref": 1.,
               "no": [0.436343e-2, -0.264054e-1, 0.922394e-1, -0.291756e-1],
               "co": [0, 1, 2, 3],

               "Trefb": 658.1, "rhorefb": 1.33, "krefb": 1.,
               "nb": [0.693347e-1, -0.280792e-1, -0.331695e-1, 0.173922e-2,
                      0.676165e-2, 0.309558e-2, 0.0, 0.0, 0.0, 0.0],
               "tb": [0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
               "db": [1, 1, 2, 2, 3, 3, 4, 4, 5, 5],
               "cb": [0]*10,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
               "Xio": 0.194e-9, "gam0": 0.0496, "qd": 1.52e-9, "Tcref": 987.15}

    _thermal = thermo0,
