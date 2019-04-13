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


class R1233zd(MEoS):
    """Multiparameter equation of state for R1233zd"""
    name = "1-chloro-3,3,3-trifluoroprop-1-ene"
    CASNumber = "102687-65-0"
    formula = "CHCl=CH-CF3"
    synonym = "R1233zd"
    _refPropName = "R1233ZD"
    _coolPropName = "R1233zd(E)"
    rhoc = unidades.Density(480.219392)
    Tc = unidades.Temperature(439.6)
    Pc = unidades.Pressure(3623.7, "kPa")
    M = 130.4944  # g/mol
    Tt = unidades.Temperature(195.15)
    Tb = unidades.Temperature(291.47)
    f_acent = 0.305
    momentoDipolar = unidades.DipoleMoment(1.44, "Debye")

    CP1 = {"ao": 4.0,
           "an": [], "pow": [],
           "ao_exp": [11.765, 8.6848],
           "exp": [630, 2230]}

    mondejar = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R1233zd(E) of Mondejar "
                    "(2013).",
        "__doi__": {"autor": "Mondejar, M.E., McLinden, M.O., Lemmon, E.W.",
                    "title": "Thermodynamic Properties of trans-1-Chloro-3,3,3"
                             "-trifluoropropene (R1233zd(E)): Vapor Pressure, "
                             "(p-ρ-T) Behavior, and Spped of Sound "
                             "Measurements, and Equation of State",
                    "ref": "J. Chem. Eng. Data 60(8) (2015) 2477-2489",
                    "doi": "10.1021/acs.jced.5b00348"},

        "R": 8.3144621,
        "cp": CP1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 550.0, "Pmax": 100000.0, "rhomax": 11.41,

        "nr1": [0.0478487, 1.60644, -2.27161, -0.530687, 0.169641],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1., 0.26, 1.02, 0.7, 0.4],

        "nr2": [-1.85458, -0.321916, 0.636411, -0.121482, -0.0262755],
        "d2": [1, 3, 2, 2, 7],
        "t2": [1.46, 2.3, 0.66, 2.7, 1.19],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [2.37362, -0.901771, -0.455962, -0.602941, -0.0594311],
        "d3": [1, 1, 3, 2, 2],
        "t3": [1.62, 1.13, 1.7, 1.35, 1.5],
        "alfa3": [0.748, 1.473, 1.39, 0.86, 1.8],
        "beta3": [1.29, 1.61, 0.8, 1.34, 0.49],
        "gamma3": [0.89, 1.13, 0.7, 0.91, 1.2],
        "epsilon3": [0.508, 0.366, 0.38, 0.773, 1.17],
        "nr4": []}

    eq = mondejar,

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.6021, 2.3265, -1.9771, -4.8451, -4.8762],
        "t": [1.0, 1.5, 2.0, 4.3, 14.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [2.13083, 0.583568, 0.247871, 0.472173],
        "t": [0.355, 0.9, 3.5, 8.0]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.0152, -6.5621, -19.427, -62.650, -181.64],
        "t": [0.397, 1.2, 3.1, 6.6, 15.0]}

    thermo0 = {"__name__": "Perkins (2017)",
               "__doi__": {
                   "autor": "Perkins, R.A., Huber, M.L.",
                   "title": "Measurement and Correlation of the Thermal "
                            "Conductivity of trans-1-Chloro-3,3,3-"
                            "trifluoropropene (R1233zd(E))",
                   "ref": "J. Chem. Eng. Data 62(9) (2017) 2659-2665",
                   "doi": "10.1021/acs.jced.7b00106"},

               "eq": 1,

               "Toref": 439.6, "koref": 1,
               "no": [-0.140033e-1, 0.378160e-1, -0.245832e-2],
               "to": [0, 1, 2],

               "Tref_res": 439.6, "rhoref_res": 480.219, "kref_res": 1,
               "nr": [0.862816e-2, 0.914709e-3, -0.208988e-1, -0.407914e-2,
                      0.511968e-1, 0.845668e-2, -0.349076e-1, -0.108985e-1,
                      0.975727e-2, 0.538262e-2, -0.926484e-3, -0.806009e-3],
               "tr": [0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1],
               "dr": [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.02, "Xio": 0.213e-9,
               "gam0": 0.059, "qd": 0.598e-9, "Tcref": 1.5*439.6}

    _thermal = thermo0,


class Test(TestCase):

    def test_Perkins(self):
        # Table 2, pag E
        # Include basic testing for Mondejar mEoS
        # The paper use a undocumented ecs viscosity correlation, so critical
        # enhancement differ
        st = R1233zd(T=300, rho=0)
        self.assertEqual(round(st.P.MPa, 2), 0)
        self.assertEqual(round(st.k.WmK, 6), 0.010659)

        st = R1233zd(T=300, rho=5.4411)
        self.assertEqual(round(st.P.MPa, 2), 0.10)
        self.assertEqual(round(st.k.WmK, 6), 0.010766)

        st = R1233zd(T=300, rho=1308.8)
        self.assertEqual(round(st.P.MPa, 2), 20.02)
        self.assertEqual(round(st.k.WmK, 6), 0.091409)

        st = R1233zd(T=445, rho=0)
        self.assertEqual(round(st.P.MPa, 2), 0)
        self.assertEqual(round(st.k.WmK, 6), 0.021758)

        st = R1233zd(T=445, rho=168.52)
        self.assertEqual(round(st.P.MPa, 2), 3.00)
        self.assertEqual(round(st.k.WmK, 6), 0.026044)
