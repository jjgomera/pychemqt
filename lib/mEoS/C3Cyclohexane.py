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


from lib.meos import MEoS
from lib import unidades


class C3Cyclohexane(MEoS):
    """Multiparameter equation of state for propylcyclohexane"""
    name = "propylcyclohexane"
    CASNumber = "1678-92-8"
    formula = "C6H11-CH2CH2CH3"
    synonym = ""
    rhoc = unidades.Density(260.0527932)
    Tc = unidades.Temperature(630.8)
    Pc = unidades.Pressure(2860.0, "kPa")
    M = 126.23922  # g/mol
    Tt = unidades.Temperature(178.2)
    Tb = unidades.Temperature(429.9)
    f_acent = 0.326
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 184

    CP1 = {"ao": 0.0,
           "an": [9.29427],
           "pow": [0.385871],
           "ao_exp": [1.37051, 106.426, 313.713],
           "exp": [173295, 561.14, 1919.52],
           "ao_hyp": [], "hyp": []}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for propylcyclohexane of Lemmon (2007).",
        "__doi__": {"autor": "Lemmon, E.W.",
                    "title": "",
                    "ref": "unpublished equation, 2007",
                    "doi": ""},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 650., "Pmax": 50000.0, "rhomax": 7.03,
        "Pmin": 0.0000007, "rhomin": 7.03,

        "nr1": [1.01911, -2.59762, 0.675152, -0.230891, 0.120966, 0.000309038],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.2, 1.2, 1.8, 1.5, 0.3, 0.9],

        "nr2": [0.526461, -0.0188462, -0.549272, -0.139233, 0.121242],
        "d2": [2, 5, 1, 4, 1],
        "t2": [1.4, 2.2, 3.7, 4.2, 2.4],
        "c2": [1, 1, 2, 2, 1],
        "gamma2": [1]*5}

    eq = lemmon,

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.76296e1, 0.16538e1, -0.28518e1, -0.28205e1, -0.28144e1],
        "exp": [1.0, 1.5, 2.7, 4.7, 15.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.39271e-1, 0.38257e2, -0.65743e2, 0.30332e2, 0.17224],
        "exp": [0.1, 0.75, 0.87, 1.0, 5.0]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.64572e1, 0.91228e1, -0.25806e2, -0.59044e2, -0.14709e3],
        "exp": [0.6, 1.8, 2.2, 6.0, 14.0]}

    visco0 = {"eq": 4, "omega": 1,
              "__name__": "Quiñones-Cisneros (2006)",
              "__doi__": {"autor": "S.E.Quiñones-Cisneros and U.K. Deiters",
                          "title": "Generalization of the Friction Theory for Viscosity Modeling",
                          "ref": "J. Phys. Chem. B, 2006, 110 (25), pp 12820–12834",
                          "doi": "10.1021/jp0618577"},

              "Tref": 630.8, "muref": 1.0,
              "ek": 507.54, "sigma": 0.6321,
              "n_ideal": [0.528175e2, -0.170572e3,  0.171218e3, -0.402745e2],
              "t_ideal": [0, 0.25, 0.5, 0.75],

              "a": [-0.132691e-3, 0.0, 0.469322e-6],
              "b": [-0.121616e-3, 0.157511e-4, 0.487973e-6],
              "c": [0.160622e-2, -0.500143e-3, 0.0],
              "A": [-0.158302e-7, 0.223800e-9, 0.0],
              "B": [0.252822e-7, 0.0, 0.0],
              "C": [0.0, 0.0, 0.0],
              "D": [0.0, 0.0, 0.0]}

    _viscosity = visco0,

    thermo0 = {"eq": 1,
               "__name__": "Perkins (2008)",
               "__doi__": {"autor": "Perkins, R.A. Hammerschmidt, U. and Huber, M.L.",
                           "title": "Measurement and Correlation of the Thermal Conductivity of Methylcyclohexane and Propylcyclohexane from 300 K to 600 K at Pressures to 60 MPa",
                           "ref": "J. Chem. Eng. Data, 2008, 53 (9), pp 2120–2127",
                           "doi": "10.1021/je800255r"},
               "__test__":
                    """
                    >>> st=C3Cyclohexane(T=300, P=1e5)
                    >>> print "%0.2f %0.3f %0.9g %0.9g " % ( \
                        st.T, st.P.MPa, st.rho, st.k.WmK)
                    300.00 0.100 788.480445 0.110772385
                    >>> st=C3Cyclohexane(T=450, P=1e5)
                    >>> print "%0.2f %0.3f %0.9g %0.9g " % ( \
                        st.T, st.P.MPa, st.rho, st.k.WmK)
                    450.00 0.100 3.52727123 0.0243185705
                    >>> st=C3Cyclohexane(T=450, P=5e7)
                    >>> print "%0.2f %0.3f %0.9g %0.9g " % ( \
                        st.T, st.P.MPa, st.rho, st.k.WmK)
                    450.00 50.000 729.367278 0.105375558
                    >>> st=C3Cyclohexane(T=600, P=1e5)
                    >>> print "%0.2f %0.3f %0.9g %0.9g " % ( \
                        st.T, st.P.MPa, st.rho, st.k.WmK)
                    600.00 0.100 2.56885356 0.0451537166
                    >>> st=C3Cyclohexane(T=600, P=4.744e6)
                    >>> print "%0.2f %0.3f %0.9g %0.9g " % ( \
                        st.T, st.P.MPa, st.rho, st.k.WmK)
                    600.00 4.744 501.697589 0.0743155503
                    >>> st=C3Cyclohexane(T=600, P=5e7)
                    >>> print "%0.2f %0.3f %0.9g %0.9g " % ( \
                        st.T, st.P.MPa, st.rho, st.k.WmK)
                    300.00 50.000 644.863706 0.0976809646
                    """, # Table 5, pag 2125

               "Tref": 630.8, "kref": 1,
               "no": [0.107402e-1, -0.609829e-1, 0.138204, -0.381213e-1],
               "co": [0, 1, 2, 3],

               "Trefb": 630.8, "rhorefb": 2.06, "krefb": 1.,
               "nb": [0.116524, -0.102821, -0.113871, 0.126431, 0.445827e-1,
                      -0.5946e-1, -0.54573600e-2, 0.98936000e-2, 0.0, 0.0],
               "tb": [0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
               "db": [1, 1, 2, 2, 3, 3, 4, 4, 5, 5],
               "cb": [0]*10,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.2415, "R0": 1.01,
               "Xio": 0.15e-9, "gam0": 0.052, "qd": 6.24e-10, "Tcref": 958.725}

    _thermal = thermo0,
