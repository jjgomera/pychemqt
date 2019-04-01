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


class MM(MEoS):
    """Multiparameter equation of state for hexamethyldisiloxane"""
    name = "hexamethyldisiloxane"
    CASNumber = "107-46-0"
    formula = "C6H18OSi2"
    synonym = "MM"
    _refPropName = "MM"
    _coolPropName = "MM"
    rhoc = unidades.Density(304.4043888253152)
    Tc = unidades.Temperature(518.69997204)
    Pc = unidades.Pressure(1939.39, "kPa")
    M = 162.3768  # g/mol
    Tt = unidades.Temperature(204.93)
    Tb = unidades.Temperature(373.401)
    f_acent = 0.418
    momentoDipolar = unidades.DipoleMoment(0.801, "Debye")
    id = 1376

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [72.110754, -10.431499],
           "ao_exp": [18.59, 29.58, 19.74, 4.87],
           "titao": [20/Tc, 1400/Tc, 3600/Tc, 6300/Tc]}

    f = 8.314472
    CP1 = {"ao": 51.894/f,
           "an": [741.34e-3/f, -416e-6/f, 70e-9/f],
           "pow": [1, 2, 3],
           "ao_exp": [], "exp": []}

    thol = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hexamethyldisiloxane of "
                    "Thol (2015).",
        "__doi__": {"autor": "Thol, M., Dubberke, F.H., Rutkai, G., Windmann, "
                             "T., Köster, A., Span, R., Vrabec, J.",
                    "title": "Fundamental equation of state correlation for "
                             "hexamethyldisiloxane based on experimental and "
                             "molecular simulation data",
                    "ref": "Fluid Phase Equilibria 418 (2016) 133-151",
                    "doi":  "10.1016/j.fluid.2015.09.047"},

        "R": 8.3144621,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": 220.0, "Tmax": 1200.0, "Pmax": 600000.0, "rhomax": 5.266,
        "M": 162.3768, "Tc": 518.7, "rhoc": 1.653, "Pc": 1931.1,

        "nr1": [0.5063651e-1, 8.604724, -9.179684, -1.146325, 0.4878559],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.346, 0.46, 1.01, 0.59],

        "nr2": [-2.434088, -1.621326, 0.6239872, -2.306057, -0.5555096e-1],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.600, 3.330, 0.750, 2.950, 0.930],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [9.385015, -2.493508, -3.308032, -0.1885803, -0.9883865e-1,
                0.1111090, 0.1061928, -0.1452454e-1],
        "d3": [1, 1, 3, 3, 1, 2, 3, 1],
        "t3": [1.33, 1.68, 1.7, 3.08, 5.41, 1.4, 1.1, 5.3],
        "alfa3": [1.0334, 1.544, 1.113, 1.113, 1.11, 7.2, 1.45, 4.73],
        "beta3": [0.4707, 0.32, 0.404, 0.517, 0.432, 7.2, 1.2, 35.8],
        "gamma3": [1.7754, 0.692, 1.242, 0.421, 0.406, 0.163, 0.795, 0.88],
        "epsilon3": [0.8927, 0.5957, 0.559, 1.056, 1.3, 0.106, 0.181, 0.525],
        "nr4": []}

    colonna = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hexamethyldisiloxane of "
                    "Colonna (2006).",
        "__doi__": {"autor": "Colonna, P., Nannan, N.R., Guardone, A., "
                             "Lemmon, E.W.",
                    "title": "Multiparameter Equations of State for Selected "
                             "Siloxanes",
                    "ref": "Fluid Phase Equilibria, 244:193-211, 2006.",
                    "doi":  "10.1016/j.fluid.2006.04.015"},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": 273.0, "Tmax": 673.0, "Pmax": 30000.0, "rhomax": 5.21,
        "Pmin": 0.00269, "rhomin": 5.2,

        "nr1": [1.01686012, -2.19713029, 0.75443188, -0.68003426, 0.19082162,
                0.10530133e-2],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.6284595, 0.30903042e-1, -0.83948727, -0.20262381,
                -0.35131597e-1, 0.25902341e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.0],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = thol, colonna

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.850230e1, 0.380300e1, -0.341500e1, -0.467900e1, -0.310600e1],
        "t": [1.0, 1.5, 1.98, 3.86, 14.6]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.4003e1, -0.6406e1, 0.115e2, -0.1004e2, 0.4e1],
        "t": [0.436, 0.827, 1.24, 1.7, 2.23]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.37421e1, -0.37087e2, 0.7546e2, -0.7167e2, -68.69, -178.4],
        "t": [0.428, 1.79, 2.28, 2.8, 7, 15.4]}


class Test(TestCase):

    def test_thol(self):
        # Test in thesis
        # Thol, M.
        # Empirical Multiparameter Equations of State Based on Molecular
        # Simulation and Hybrid Data Sets
        # PhD thesis, Ruhr-Universität Bochum, 2015.

        # Appendix A, Pag 259

        # The values in table are not good, trying to use the values report by
        # CoolProp

        # from CoolProp. CoolProp import PropsSI
        # P = PropsSI("P", "T", 250, "Dmolar", 0.0001, "MM")
        P = 0.20785
        # Cp = PropsSI("Cpmolar", "T", 250, "Dmolar", 0.0001, "MM")
        Cp = 216.5

        st = MM(T=250, rhom=0.0001, eq="thol")
        self.assertEqual(round(st.P.kPa, 5), P)
        self.assertEqual(round(st.cpM.JmolK, 1), Cp)

        # st = MM(T=250, rhom=0.0001, eq="thol")
        # self.assertEqual(round(st.P.MPa, 8), 2.3550378)
        # self.assertEqual(round(st.cpM.JmolK, 5), 290.08362)
        # self.assertEqual(round(st.w, 4), 1068.3855)
        # self.assertEqual(round(st.hM.Jmol, 3), -38660.059)
        # self.assertEqual(round(st.sM.JmolK, 6), -126.50073)
        # self.assertEqual(round(st.aM.Jmol, 6), -7505.8829)

        # st = MM(T=250, rhom=5, eq="thol")
        # self.assertEqual(round(st.P.MPa, 8), 2.0772979e-4)
        # self.assertEqual(round(st.cpM.JmolK, 5), 2.1658262e-2)
        # self.assertEqual(round(st.w, 4), 115.31572)
        # self.assertEqual(round(st.hM.Jmol, 3), 1715.1940)
        # self.assertEqual(round(st.sM.JmolK, 6), 38.943471)
        # self.assertEqual(round(st.aM.Jmol, 6), -10097.972)

        # st = MM(T=400, rhom=0.05, eq="thol")
        # self.assertEqual(round(st.P.MPa, 8), 0.15367468)
        # self.assertEqual(round(st.cpM.JmolK, 5), 293.72934)
        # self.assertEqual(round(st.w, 4), 134.70433)
        # self.assertEqual(round(st.hM.Jmol, 3), 38493.817)
        # self.assertEqual(round(st.sM.JmolK, 6), 99.143201)
        # self.assertEqual(round(st.aM.Jmol, 6), -4236.9572)

        # st = MM(T=400, rhom=4.5, eq="thol")
        # self.assertEqual(round(st.P.MPa, 8), 40.937214)
        # self.assertEqual(round(st.cpM.JmolK, 5), 339.40134)
        # self.assertEqual(round(st.w, 4), 930.21218)
        # self.assertEqual(round(st.hM.Jmol, 3), 1367.2106)
        # self.assertEqual(round(st.sM.JmolK, 6), 11.063887)
        # self.assertEqual(round(st.aM.Jmol, 6), 149.39229)

        # st = MM(T=560, rhom=4.5, eq="thol")
        # self.assertEqual(round(st.P.MPa, 8), 123.02530)
        # self.assertEqual(round(st.cpM.JmolK, 5), 387.27688)
        # self.assertEqual(round(st.w, 4), 1132.8991)
        # self.assertEqual(round(st.hM.Jmol, 3), 83661.459)
        # self.assertEqual(round(st.sM.JmolK, 6), 119.31485)
        # self.assertEqual(round(st.aM.Jmol, 6), -10493.815)
