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


class MD3M(MEoS):
    """Multiparameter equation of state for dodecamethylpentasiloxane"""
    name = "dodecamethylpentasiloxane"
    CASNumber = "141-63-9"
    formula = "C12H36Si5O4"
    synonym = "MD3M"
    _refPropName = "MD3M"
    _coolPropName = "MD3M"
    rhoc = unidades.Density(269.3873)
    Tc = unidades.Temperature(628)
    Pc = unidades.Pressure(953.95, "kPa")
    M = 384.839  # g/mol
    Tt = unidades.Temperature(192.0)
    Tb = unidades.Temperature(503.02)
    f_acent = 0.73
    momentoDipolar = unidades.DipoleMoment(1.223, "Debye")

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [68.1167204166, -29.8091965426],
           "ao_exp": [81.2386, 61.191, 51.1798],
           "titao": [610/Tc, 2500/Tc, 7500/Tc]}

    f = 8.314472
    CP1 = {"ao": 463.2/f,
           "ao_sinh": [957.2/f, ], "sinh": [2117.1],
           "ao_cosh": [738.3/f], "cosh": [908.5]}

    thol = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for "
                    "decamethylcyclopentasiloxane of Thol (2019).",
        "__doi__": {
            "autor": "Thol, M., Javed, M.A., Baumhögger, E., Span, R., "
                     "Vrabec, J.",
            "title": "Thermodynamic Properties of Dodecamethylpentasiloxane, "
                     "Tetradecamethylhexasiloxane, and "
                     "Decamethylcyclopentasiloxane",
            "ref": "Ind. Eng. Chem. Res. 58(22) (2019) 9617-9635",
            "doi": "10.1021/acs.iecr.9b00608"},

        "R": 8.3144598,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": 200, "Tmax": 730, "Pmax": 125000,

        "nr1": [0.040674325, 4.4936509, -6.0327468, -1.0842396, 0.65985153],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.37, 0.718, 0.79, 0.59],

        "nr2": [-2.3011802, -1.5022099, 0.5051725, -2.2363839, -0.071582853],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.38, 3.14, 0.62, 2.08, 1.042],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [4.7053488, -0.774783117, -0.68302991, 0.41657104, -1.1441135],
        "d3": [1, 1, 3, 2, 2],
        "t3": [0.9, 0.86, 2.06, 0.55, 0.69],
        "alfa3": [1.043, 20, 1.08, 0.47, 1.085],
        "beta3": [0.86, 1099, 0.95, 0.1, 1.85],
        "gamma3": [1.357, 1.097, 1.030, 1.020, 0.800],
        "epsilon3": [0.725, 0.940, 0.546, 0.680, 0.495]}

    colonna = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for MD3M of Colonna (2006).",
        "__doi__": {"autor": "Colonna, P., Nannan, N.R., Guardone, A.",
                    "title": "Multiparameter equations of state for siloxanes:"
                             " [(CH3)3-Si-O1/2]2-[O-Si-(CH3)2]i=1,…,3, and "
                             "[O-Si-(CH3)2]6",
                    "ref": "Fluid Phase Equilibria 263:115-130, 2008",
                    "doi": "10.1016/j.fluid.2007.10.001"},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",
        "Tc": 628, "Pc": 945, "rhoc": 0.685798,

        "Tmin": Tt, "Tmax": 673.0, "Pmax": 30000.0, "rhomax": 2.54,

        "nr1": [1.20540386, -2.42914797, 0.69016432, -0.69268041, 0.18506046,
                0.31161436e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.99862519, 0.74229034e-1, -0.80259136, -0.20865337,
                -0.36461791e-1, 0.19174051e-1],
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
        "sigma": [0.03972], "exp": [1.254]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-9.6774, 3.973, -3.701, -9.1232, -5.467],
        "t": [1.0, 1.5, 1.83, 3.54, 11.9]}
    _liquid_Density = {
        "eq": 1,
        "n": [3.326, 3.889, -2.363, -2.709, 1.325],
        "t": [0.42, 1.46, 0.9, 2.15, 3.15]}
    _vapor_Density = {
        "eq": 2,
        "n": [-4.0084, -7.913, -79.392, -28.572, -211.86, -1800],
        "t": [0.441, 1.244, 5.88, 3.03, 12.7, 32]}

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
              "Tc": 628.96, "rhoc": 0.7*384.839,

              "ek": 499.5, "sigma": 0.911, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [1.45796, -0.15796], "psi_d": [0, 1],
              "fint": [0.00132], "fint_t": [0],
              "chi": [1.72213], "chi_d": [0],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.33e-9, "gam0": 0.066, "qd": 1.127e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""

    def test_thol(self):
        """Table 6, Pag. 9621"""
        st = MD3M(T=300, rhom=2.4)
        self.assertEqual(round(st.P.MPa, 7), 56.5643398)
        self.assertEqual(round(st.hM.Jmol, 3), -133761.828)
        self.assertEqual(round(st.sM.JmolK, 6), -403.543152)
        self.assertEqual(round(st.w, 5), 1241.26649)
        self.assertEqual(round(st.aM.Jmol, 4), -36267.3571)

        st = MD3M(T=390, rhom=0.0005)
        self.assertEqual(round(st.P.MPa, 10), 0.0016139843)
        self.assertEqual(round(st.hM.Jmol, 4), -31595.5295)
        self.assertEqual(round(st.sM.JmolK, 7), -48.7220551)
        self.assertEqual(round(st.w, 7), 92.0127172)
        self.assertEqual(round(st.aM.Jmol, 4), -15821.8965)

        st = MD3M(T=450, rhom=0.003)
        self.assertEqual(round(st.P.MPa, 10), 0.0110320958)
        self.assertEqual(round(st.hM.Jmol, 5), 7104.19531)
        self.assertEqual(round(st.sM.JmolK, 7), 27.6618505)
        self.assertEqual(round(st.w, 7), 97.5667091)
        self.assertEqual(round(st.aM.Jmol, 5), -9021.00267)

        st = MD3M(T=450, rhom=2)
        self.assertEqual(round(st.P.MPa, 7), 18.3032601)
        self.assertEqual(round(st.hM.Jmol, 4), -38524.8786)
        self.assertEqual(round(st.sM.JmolK, 6), -101.224710)
        self.assertEqual(round(st.w, 6), 728.435652)
        self.assertEqual(round(st.aM.Jmol, 5), -2125.38904)

        st = MD3M(T=600, rhom=2)
        self.assertEqual(round(st.P.MPa, 7), 70.6395352)
        self.assertEqual(round(st.hM.Jmol, 3), 100062.177)
        self.assertEqual(round(st.sM.JmolK, 6), 113.577309)
        self.assertEqual(round(st.w, 6), 902.133167)
        self.assertEqual(round(st.aM.Jmol, 5), -3403.97539)

    # def test_Huber(self):
    #     """Table 7, pag 266"""
    #     st = MD3M(T=566.1, rhom=1.47)
    #     self.assertEqual(round(st.mu.muPas, 4), 142.9225)
    #     self.assertEqual(round(st.k.mWmK, 4), 63.0474)
