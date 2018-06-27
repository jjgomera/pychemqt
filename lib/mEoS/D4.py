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


class D4(MEoS):
    """Multiparameter equation of state for octamethylcyclotetrasiloxane"""
    name = "octamethylcyclotetrasiloxane"
    CASNumber = "556-67-2"
    formula = "C8H24O4Si4"
    synonym = "D4"
    _refPropName = "D4"
    rhoc = unidades.Density(307.0335906736056)
    Tc = unidades.Temperature(586.49127187)
    Pc = unidades.Pressure(1332.0, "kPa")
    M = 296.61576  # g/mol
    Tt = unidades.Temperature(290.25)
    Tb = unidades.Temperature(448.504)
    f_acent = 0.592
    momentoDipolar = unidades.DipoleMoment(1.090, "Debye")
    # id=1430

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],

           # The paper parameter are incorrect, using alternate parameter
           # to meet reference state OTO
           # "ao_pow": [71.1636049792958, -21.6743650975623],
           "ao_pow": [44.72889170669655, -4.9687471148991968],

           "ao_exp": [0.292757, 38.2456, 58.975],
           "titao": [40/Tc, 200/Tc, 1800/Tc]}

    CP1 = {"ao": -18.256,
           "an": [1427.2e-3, -990.20e-6, 300.0e-9], "pow": [1, 2, 3],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    thol = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for D4 of Thol (2006).",
        "__doi__": {"autor": "Thol, M.",
                    "title": "Empirical Multiparameter Equations of State "
                             "Based on Molecular Simulation and Hybrid Data "
                             "Sets",
                    "ref": "PhD thesis, Ruhr-Universität Bochum, 2015.",
                    "doi":  ""},

        "R": 8.3144621,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 1200.0, "Pmax": 520000.0, "rhomax": 5.266,
        # "Pmin": 0.00269, "rhomin": 5.2,
        "Tc": 586.5, "rhoc": 1.043, "Pc": 1347,

        "nr1": [5.273743e-2, 4.176401, -4.737070, -1.289588, 5.272749e-1],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.27, 0.51, 0.998, 0.56],

        "nr2": [-2.558391, -0.9726737, 0.7208209, -4.789456e-1, -5.563239e-2],
        "d2": [1, 3, 2, 2, 7],
        "t2": [1.75, 3.09, 0.79, 2.71, 0.998],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [3.766589, 8.786997e-2, -1.267646e-1, -1.004246, -1.641887],
        "d3": [1, 1, 3, 2, 2],
        "t3": [0.93, 3.17, 1.08, 1.41, 0.89],
        "alfa3": [0.861, 1.114, 1.01, 1.11, 1.032],
        "beta3": [0.75, 0.55, 1.0, 0.47, 1.36],
        "gamma3": [1.124, 1.388, 1.148, 1.197, 0.817],
        "epsilon3": [0.926, 1.3, 1.114, 0.996, 0.483],
        "nr4": []}

    colonna = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for D4 of Colonna (2006).",
        "__doi__": {"autor": "Colonna, P., Nannan, N.R., Guardone, A., "
                             "Lemmon, E.W.",
                    "title": "Multiparameter Equations of State for Selected "
                             "Siloxanes",
                    "ref": "Fluid Phase Equilibria, 244:193-211, 2006.",
                    "doi":  "10.1016/j.fluid.2006.04.015"},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": 300.0, "Tmax": 673.0, "Pmax": 30000.0, "rhomax": 3.21,
        "Pmin": 0.0696, "rhomin": 3.2,

        "nr1": [1.05392408, -2.22981918, 0.77573923, -0.6937405, 0.18721557,
                0.42193330e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.70301835, 0.47851888e-1, -0.8025348, -0.18968872,
                -0.22211781e-1, 0.60103354e-2],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.0],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = thol, colonna

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.87935e1, 0.27204e1, -0.48174e1, -0.69086e1],
        "exp": [1.0, 1.5, 2.2, 4.4]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.14563e1, -0.94215, 0.45065e1, -0.27688e1, 0.8745],
        "exp": [0.24, 0.5, 0.75, 1.0, 2.0]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.16204e1, -0.57888e1, -0.24291e2, 0.53567e2, -0.12135e3,
               -0.10976e4],
        "exp": [0.31, 0.78, 2.5, 4.4, 5.0, 15.0]}


class Test(TestCase):

    def test_thol(self):
        # Appendix A, Pag 259

        # The values in table are not good, the entropy and enthalpy reference
        # state don't is OTO, using the custom integration parameters check
        # reference state values.
        st = D4(T=298.15, P=101325, eq="thol")
        self.assertEqual(round(st.h, 10), 0)
        self.assertEqual(round(st.s, 10), 0)

        st = D4(T=350, rhom=0.001, eq="thol")
        self.assertEqual(round(st.P.MPa, 10), 2.8952836e-3)
        self.assertEqual(round(st.cpM.JmolK, 2), 422.18)
        self.assertEqual(round(st.w, 4), 99.5573)

        st = D4(T=350, rhom=3.2, eq="thol")
        self.assertEqual(round(st.P.MPa, 6), 38.313384)
        self.assertEqual(round(st.cpM.JmolK, 2), 515.36)
        self.assertEqual(round(st.w, 4), 1003.2309)

        st = D4(T=500, rhom=0.08, eq="thol")
        self.assertEqual(round(st.P.MPa, 8), 0.28137478)
        self.assertEqual(round(st.cpM.JmolK, 2), 549.71)
        self.assertEqual(round(st.w, 4), 100.549)

        st = D4(T=500, rhom=2.5, eq="thol")
        self.assertEqual(round(st.P.MPa, 7), 8.0209605)
        self.assertEqual(round(st.cpM.JmolK, 1), 612.9)
        self.assertEqual(round(st.w, 3), 475.065)

        st = D4(T=600, rhom=3, eq="thol")
        self.assertEqual(round(st.P.MPa, 5), 126.85572)
        self.assertEqual(round(st.cpM.JmolK, 2), 650.74)
        self.assertEqual(round(st.w, 3), 1071.662)
