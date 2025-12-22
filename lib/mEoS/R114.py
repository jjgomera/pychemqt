#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2025, Juan José Gómez Romera <jjgomera@gmail.com>

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
from lib.mEoS import C3


class R114(MEoS):
    """Multiparameter equation of state for R114"""
    name = "1,2-dichloro-1,1,2,2-tetrafluoroethane"
    CASNumber = "76-14-2"
    formula = "CClF2CClF2"
    synonym = "R114"
    _refPropName = "R114"
    _coolPropName = "R114"
    rhoc = unidades.Density(579.969)
    Tc = unidades.Temperature(418.83)
    Pc = unidades.Pressure(3257.0, "kPa")
    M = 170.921  # g/mol
    Tt = unidades.Temperature(180.63)
    Tb = unidades.Temperature(276.741)
    f_acent = 0.25253
    momentoDipolar = unidades.DipoleMoment(0.658, "Debye")
    id = 231

    f = 1/8.31451*170.93
    CP1 = {"ao": 0.97651380e-1*f,
           "an": [0.3240861e-2*f, -0.5895364e-5*f, 0.6737929e-8*f,
                  -0.3546364e-11*f],
           "pow": [1, 2, 3, 4]}

    platzer = {
        "__type__": "Helmholtz",
        "__name__": "Bender equation of state for R-114 of Platzer (1990)",
        "__doi__": {"autor": "Platzer, B., Polt, A., Maurer, G.",
                    "title": "Thermophysical Properties of Refrigerants",
                    "ref": "Berlin:  Springer-Verlag, 1990.",
                    "doi": ""},

        "R": 8.31451,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": 273.15, "Tmax": 507.0, "Pmax": 21000.0, "rhomax": 10,

        "nr1": [-0.340776521414, 0.323001398420, -0.424950537596e-1,
                0.107938879710e1, -0.199243619673e1, -0.155135133506,
                -0.121465790553, -0.165038582393e-1, -0.186915808643,
                0.308074612567, 0.115861416115, 0.276358316589e-1,
                0.108043243088, 0.460683793064e-1, -0.174821616881,
                0.317530854287e-1],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.340776521414, -0.323001398420, 0.424950537596e-1,
                -0.166940100976e1, 0.408693082002e1, -0.241738963889e1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [1.21103865]*6}

    eq = (platzer, )
    _PR = [-0.1804, -16.3839]

    _surface = {"sigma": [0.05239], "exp": [1.258]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.72195e1, 0.16357e1, -0.14576e1, -0.69580e1, 0.57181e1],
        "t": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.43023, 0.22722e2, -0.27118e2, 0.13247e2, -0.90529e1],
        "t": [0.095, 0.93, 1.1, 2.0, 3.0]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.46609, -6.8355, -167.15, 1.5805e4, -3.1859e4, 2.1548e4],
        "t": [0.09, 0.76, 4.0, 6.5, 7.0, 8.0]}

    trnECS = {"__name__": "Huber (2018)",

              "__doi__": {
                  "autor": "Huber, M.L.",
                  "title": "Models for Viscosity, Thermal Conductivity, and "
                           "Surface Tension of Selected Pure Fluids as "
                           "Implemented in REFPROP v10.0",
                  "ref": "NISTIR 8209",
                  "doi": "10.6028/NIST.IR.8209"},

              "eq": "ecs",
              "ref": C3,
              "visco": "visco1",

              "ek": 174, "sigma": 0.648, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [1.36002, -0.209356, 0.0373222],
              "psi_d": [0, 1, 2],
              "fint": [0.00132], "fint_t": [0],
              "chi": [1.2005, -0.0533827], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.223e-9, "gam0": 0.059, "qd": 0.656e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )

    thermo0 = {"__name__": "Krauss (1989)",
               "__doi__": {
                   "autor": "Krauss, R., Stephan, K.",
                   "title": "Thermal Conductivity of Refrigerants in a Wide "
                            "Range of Temperature and Pressure",
                   "ref": "J. Phys. Chem. Ref. Data 18(1) (1989) 43-76",
                   "doi": "10.1063/1.555842"},

               "eq": 1,

               # Typo in paper
               # Temperature reducing value is critical value
               # Thermal conductivity reducing value using the exact value from
               # dimensional analysis with molecular weight in kg/mol
               # R**(5/6)*Pc**(2/3)/Tc**(1/6)/(M/1000)**0.5/Na**(1/3)
               "Toref": 418.85, "koref": 1.3812067e-3,
               "no": [-4.4141934, 16.535328],
               "to": [0, 1],

               "rhoref_res": 3.4*170.922, "kref_res": 1.3812067e-3,
               "nr": [7.1429545, -2.7453806, 2.453072],
               "tr": [0, 0, 0],
               "dr": [1, 2, 3]}

    _thermal = (trnECS, thermo0)


class Test(TestCase):
    """Testing"""
    def test_Huber(self):
        """Table 7, pag 266"""
        st = R114(T=376.9, rhom=6.897)
        self.assertEqual(round(st.mu.muPas, 4), 152.2319)
        self.assertEqual(round(st.k.mWmK, 4), 45.8090)

    def test_krauss(self):
        """Selected point from Table C5 and C6, pag 72"""
        # The values differ because the paper use and old EoS don't
        # implemented in pychemqt
        self.assertEqual(round(R114(T=280, P=1e5, thermal=1).k.mWmK, 2), 9.30)
        self.assertEqual(round(R114(T=400, P=5e6, thermal=1).k.mWmK, 2), 43.75)
        self.assertEqual(round(R114(T=300, P=1e7, thermal=1).k.mWmK, 2), 67.64)
        self.assertEqual(round(R114(T=480, P=2e7, thermal=1).k.mWmK, 2), 45.59)

        # Saturation point, Table C6
        st = R114(T=300, x=0.5, thermal=1)
        self.assertEqual(round(st.P.bar, 4), 2.2746)
        self.assertEqual(round(st.Liquido.k.mWmK, 2), 63.80)

        st = R114(T=400, x=0.5, thermal=1)
        self.assertEqual(round(st.P.bar, 3), 23.530)
        self.assertEqual(round(st.Liquido.k.mWmK, 2), 38.68)
