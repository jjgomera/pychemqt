#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2022, Juan José Gómez Romera <jjgomera@gmail.com>

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


class R115(MEoS):
    """Multiparameter equation of state for R115"""
    name = "chloropentafluoroethane"
    CASNumber = "76-15-3"
    formula = "CClF2CF3"
    synonym = "R115"
    _refPropName = "R115"
    _coolPropName = "R115"
    rhoc = unidades.Density(614.77633568)
    Tc = unidades.Temperature(353.1)
    Pc = unidades.Pressure(3129.0, "kPa")
    M = 154.466416  # g/mol
    Tt = unidades.Temperature(173.75)
    Tb = unidades.Temperature(233.9)
    f_acent = 0.248
    momentoDipolar = unidades.DipoleMoment(0.52, "Debye")
    id = 229

    CP1 = {"ao": 4,
           "ao_exp": [7.142, 10.61], "exp": [289, 1301]}

    CP2 = {"ao": 2.4409547,
           "an": [0.053544743, -0.81861429e-4, 0.10410538e-6, -0.71645701e-10],
           "pow": [1, 2, 3, 4]}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-115 of Lemmon "
                    "and Span (2013)",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Thermodynamic Properties of R-227ea, R-365mfc, "
                             "R-115, and R-13I1",
                    "ref": "J. Chem. Eng. Data, 60(12) (2015) 3745-3758",
                    "doi": "10.1021/acs.jced.5b00684"},

        "R": 8.3144621,
        "cp": CP1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 550.0, "Pmax": 60000.0, "rhomax": 11.3,

        "nr1": [1.20873, -3.54460, 0.745302, 0.114128, 0.000436572],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.988385, 1.13878, -0.0215633, -0.630230, 0.0167901,
                -0.149412, -0.0271153],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    platzer = {
        "__type__": "Helmholtz",
        "__name__": "Bender equation of state for R-115 of Platzer (1990)",
        "__doi__": {"autor": "Platzer, B., Polt, A., Maurer, G.",
                    "title": "Thermophysical Properties of Refrigerants",
                    "ref": "Berlin:  Springer-Verlag, 1990.",
                    "doi": ""},

        "R": 8.31451,
        "cp": CP2,
        "ref": "NBP",

        "Tmin": 200.0, "Tmax": 450.0, "Pmax": 7000.0, "rhomax": 10.7,

        "nr1": [-0.377294477051, -0.695891789165e-1, 0.206972205161,
                0.266609543946, -0.117158857583e1, 0.817521154071,
                -0.978729789251, -0.174482448760, 0.143598704796e1,
                -0.265460417723e1, 0.165212655822e1, -0.588257570097,
                0.738774518022, 0.296779702685, -0.534330750773,
                0.659766160237e-1],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.377294477051, 0.695891789165e-1, -0.206972205161,
                -0.350603135603, 0.108682541098e1, -0.619304197853],
        "d2": [2]*6,
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [0, 0, 0, 2, 2, 2],
        "gamma2": [1.50553819]*6}

    eq = lemmon, platzer
    _PR = [-0.1373, -19.0]

    _surface = {"sigma": [0.04771], "exp": [1.246]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.77016e1, 0.43462e1, -0.40020e1, -0.65510e1, 0.39278e1],
        "t": [1.0, 1.5, 1.9, 5.2, 6.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.18245e2, -0.57373e2, 0.78511e2, -0.50979e2, 0.14361e2],
        "t": [0.556, 0.75, 0.95, 1.2, 1.5]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.5696, -8.3593, -340.07, 401.14, 84.442, -221.37],
        "t": [0.421, 1.5, 4.7, 5.0, 5.4, 6.0]}

    trnECS = {"__name__": "Huber (2003)",

              "__doi__": {
                  "autor": "Huber, M.L., Laesecke, A., Perkins, R.A.",
                  "title": "Model for the Viscosity and Thermal Conductivity "
                           "of Refrigerants, Including a New Correlation for "
                           "the Viscosity of R134a",
                  "ref": "Ind. Eng. Chem. Res., 42(13) (2003) 3163-3178",
                  "doi": "10.1021/ie0300880"},

              "eq": "ecs",

              "ref": C3,
              "visco": "visco1",
              "thermo": "thermo0",

              "ek": 201.9, "sigma": 0.5876, "omega": 5,

              "psi": [1.1838, -5.91896e-2], "psi_d": [0, 1],
              "fint": [1.25079e-3, 2.96636e-7], "fint_t": [0, 1],
              "chi": [1.0343, -2.16614e-3], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.03, "Xio": 0.194e-9,
              "gam0": 0.0496, "qd": 3.72933e-10, "Tcref": 1.5*Tc}

    _viscosity = trnECS,
    _thermal = trnECS,


class Test(TestCase):

    def test_lemmon(self):
        # Table 7, Pag 3754
        st = R115(T=300, rhom=0)
        self.assertEqual(round(st.P.MPa, 6), 0.0)
        self.assertEqual(round(st.cvM.JmolK, 4), 102.2181)
        self.assertEqual(round(st.cpM.JmolK, 4), 110.5326)
        self.assertEqual(round(st.w, 4), 132.1423)

        st = R115(T=300, rhom=10)
        self.assertEqual(round(st.P.MPa, 5), 57.40521)
        self.assertEqual(round(st.cvM.JmolK, 4), 111.0519)
        self.assertEqual(round(st.cpM.JmolK, 4), 146.0220)
        self.assertEqual(round(st.w, 4), 728.9583)

        st = R115(T=354, rhom=4)
        self.assertEqual(round(st.P.MPa, 6), 3.190363)
        self.assertEqual(round(st.cvM.JmolK, 4), 137.4032)
        self.assertEqual(round(st.cpM.JmolK, 3), 5990.358)
        self.assertEqual(round(st.w, 5), 70.43849)
