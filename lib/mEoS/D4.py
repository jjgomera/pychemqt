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
from lib.mEoS import N2


class D4(MEoS):
    """Multiparameter equation of state for octamethylcyclotetrasiloxane"""
    name = "octamethylcyclotetrasiloxane"
    CASNumber = "556-67-2"
    formula = "C8H24O4Si4"
    synonym = "D4"
    _refPropName = "D4"
    _coolPropName = "D4"
    rhoc = unidades.Density(307.0335906736056)
    Tc = unidades.Temperature(586.49127187)
    Pc = unidades.Pressure(1332.0, "kPa")
    M = 296.61576  # g/mol
    Tt = unidades.Temperature(290.25)
    Tb = unidades.Temperature(448.504)
    f_acent = 0.592
    momentoDipolar = unidades.DipoleMoment(1.090, "Debye")
    # id=1430

    # Integration constant given in paper don't return good values, using Cp
    # Fi1 = {"ao_log": [1, 3],
    #        "pow": [0, 1], "ao_pow": [71.163605, -21.674365],
    #        "ao_exp": [0.292757, 38.2456, 58.975],
    #        "titao": [40/Tc, 200/Tc, 1800/Tc]}

    CP2 = {"ao": 4,
           "ao_exp": [0.292757, 38.2456, 58.975],
           "exp": [40, 200, 1800]}

    f = 8.314472
    CP1 = {"ao": -18.256/f,
           "an": [1427.2e-3/f, -990.20e-6/f, 300.0e-9/f], "pow": [1, 2, 3]}

    thol = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for D4 of Thol (2016).",
        "__doi__": {
            "autor": "Thol, M., Rutkai, G., Köster, A., Dubberke, F.H., "
                     "Windmann, T., Span, R., Vrabec, J.",
            "title": "Thermodynamic Properties of "
                     "Octamiethylciyclotetrasilosane",
            "ref": "J. Chem. Eng. Data 61(7) (2016) 2580-2595",
            "doi": "10.1021/acs.jced.6b00261"},

        "R": 8.3144621,
        "Tc": 586.5, "rhoc": 1.043, "Pc": 1347.2,
        "cp": CP2,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 1200.0, "Pmax": 520000.0,

        "nr1": [0.05273743, 4.176401, -4.73707, -1.289588, 0.5272749],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.27, 0.51, 0.998, 0.56],

        "nr2": [-2.558391, -0.9726737, 0.7208209, -0.4789456, -0.05563239],
        "d2": [1, 3, 2, 2, 7],
        "t2": [1.75, 3.09, 0.79, 2.71, 0.998],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [3.766589, 0.08786997, -0.1267646, -1.004246, -1.641887],
        "d3": [1, 1, 3, 2, 2],
        "t3": [0.93, 3.17, 1.08, 1.41, 0.89],
        "alfa3": [0.861, 1.114, 1.01, 1.11, 1.032],
        "beta3": [0.75, 0.55, 1.0, 0.47, 1.36],
        "gamma3": [1.124, 1.388, 1.148, 1.197, 0.817],
        "epsilon3": [0.926, 1.3, 1.114, 0.996, 0.483]}

    colonna = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for D4 of Colonna (2006).",
        "__doi__": {
            "autor": "Colonna, P., Nannan, N.R., Guardone, A., Lemmon, E.W.",
            "title": "Multiparameter Equations of State for Selected "
                     "Siloxanes",
            "ref": "Fluid Phase Equilibria, 244 (2006) 193-211",
            "doi": "10.1016/j.fluid.2006.04.015"},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": 300.0, "Tmax": 673.0, "Pmax": 30000.0, "rhomax": 3.21,

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

    _surface = {
        "__doi__": {
            "autor": "Mulero, A., Cachadiña, I.",
            "title": "Recommended Correlations for the Surface Tension of "
                     "Several Fluids Included in the REFPROP Program",
            "ref": "J. Phys. Chem. Ref. Data 43(2) (2014) 023104",
            "doi": "10.1063/1.4878755"},
        "sigma": [0.04246], "exp": [1.207]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-9.2842, 3.8173, -4.4415, -7.7628, -6.9289],
        "t": [1.0, 1.5, 2.1, 15, 3.9]}
    _liquid_Density = {
        "eq": 1,
        "n": [2.7216, -1.5754, 3.9887, -3.7683, 1.9445],
        "t": [0.38, 0.89, 1.44, 2.06, 2.78]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.745, -9.2075, -71.786, 108.85, -141.61, -227.19],
        "t": [0.416, 1.35, 3.8, 4.8, 5.8, 14]}

    trnECS = {"__name__": "Huber (2018)",

              "__doi__": {
                  "autor": "Huber, M.L.",
                  "title": "Models for Viscosity, Thermal Conductivity, and "
                           "Surface Tension of Selected Pure Fluids as "
                           "Implemented in REFPROP v10.0",
                  "ref": "NISTIR 8209",
                  "doi": "10.6028/NIST.IR.8209"},

              "eq": "ecs",
              "ref": N2,
              "visco": "visco1",

              "ek": 465.7, "sigma": 0.798, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [2.42579, -0.871777, 0.137283], "psi_d": [0, 1, 2],
              "fint": [0.00132], "fint_t": [0],
              "chi": [1.43353, 0.0407501], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.298e-9, "gam0": 0.064, "qd": 0.983e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""
    def test_thol(self):
        """Table S7, Pag 3"""

        st = D4(T=350, rhom=0.001, eq="thol")
        self.assertEqual(round(st.P.MPa, 10), 0.0028952836)
        self.assertEqual(round(st.cpM.JmolK, 9), 422.187662524)
        self.assertEqual(round(st.w, 6), 99.557286)
        self.assertEqual(round(st.hM.Jmol, 4), -5606.6834)
        self.assertEqual(round(st.sM.JmolK, 7), 3.7361111)
        self.assertEqual(round(st.aM.Jmol, 4), -9809.6059)

        st = D4(T=350, rhom=3.2, eq="thol")
        self.assertEqual(round(st.P.MPa, 6), 38.313384)
        self.assertEqual(round(st.cpM.JmolK, 9), 515.367386403)
        self.assertEqual(round(st.w, 4), 1003.2303)
        self.assertEqual(round(st.hM.Jmol, 3), -47367.713)
        self.assertEqual(round(st.sM.JmolK, 5), -151.31880)
        self.assertEqual(round(st.aM.Jmol, 4), -6379.0666)

        st = D4(T=500, rhom=0.08, eq="thol")
        self.assertEqual(round(st.P.MPa, 8), 0.28137478)
        self.assertEqual(round(st.cpM.JmolK, 9), 549.715457833)
        self.assertEqual(round(st.w, 5), 100.54897)
        self.assertEqual(round(st.hM.Jmol, 3), 63803.842)
        self.assertEqual(round(st.sM.JmolK, 5), 131.62219)
        self.assertEqual(round(st.aM.Jmol, 4), -5524.4399)

        st = D4(T=500, rhom=2.5, eq="thol")
        self.assertEqual(round(st.P.MPa, 7), 8.0209605)
        self.assertEqual(round(st.cpM.JmolK, 8), 612.92148189)
        self.assertEqual(round(st.w, 5), 475.06521)
        self.assertEqual(round(st.hM.Jmol, 3), 31400.054)
        self.assertEqual(round(st.sM.JmolK, 6), 59.639266)
        self.assertEqual(round(st.aM.Jmol, 4), -1627.9631)

        st = D4(T=600, rhom=3, eq="thol")
        self.assertEqual(round(st.P.MPa, 5), 126.85572)
        self.assertEqual(round(st.cpM.JmolK, 9), 650.740004998)
        self.assertEqual(round(st.w, 4), 1071.6620)
        self.assertEqual(round(st.hM.Jmol, 2), 118708.24)
        self.assertEqual(round(st.sM.JmolK, 5), 141.34253)
        self.assertEqual(round(st.aM.Jmol, 4), -8382.5142)

    def test_Huber(self):
        """Table 7, pag 266"""
        st = D4(T=527.8, rhom=2.176)
        self.assertEqual(round(st.mu.muPas, 3), 193.593)
        # self.assertEqual(round(st.k.mWmK, 4), 63.1859)
        self.assertEqual(round(st.k.mWmK, 4), 63.1860)
