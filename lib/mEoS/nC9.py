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


class nC9(MEoS):
    """Multiparameter equation of state for n-nonane"""
    name = "nonane"
    CASNumber = "111-84-2"
    formula = "CH3-(CH2)7-CH3"
    synonym = ""
    _refPropName = "NONANE"
    _coolPropName = "n-Nonane"
    rhoc = unidades.Density(232.1417)
    Tc = unidades.Temperature(594.55)
    Pc = unidades.Pressure(2281.0, "kPa")
    M = 128.2551  # g/mol
    Tt = unidades.Temperature(219.7)
    Tb = unidades.Temperature(423.91)
    f_acent = 0.4433
    momentoDipolar = unidades.DipoleMoment(0.07, "Debye")
    id = 13

    Fi1 = {"ao_log": [1, 16.349],
           "pow": [0, 1],
           "ao_pow": [10.7927224829, -8.2418318753],
           "ao_exp": [24.926, 24.842, 11.188, 17.483],
           "titao": [1221/Tc, 2244/Tc, 5008/Tc, 11724/Tc]}

    Fi2 = {"ao_log": [1, 3.0],
           "pow": [0, 1],
           "ao_pow": [16.313913248, -102.160247463],
           "ao_sinh": [18.0241, 53.3415], "sinh": [0.263819696, 2.848860483],
           "ao_cosh": [38.1235], "cosh": [1.370586158]}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for nonane of Lemmon "
                    "and Span (2006)",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 "
                             "Industrial Fluids",
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi":  "10.1021/je050186n"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 800000.0, "rhomax": 6.06,

        "nr1": [1.1151, -2.7020, 0.83416, -0.38828, 0.1376, 0.00028185],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.62037, 0.015847, -0.61726, -0.15043, -0.012982, 0.0044325],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for nonane of Kunz and "
                    "Wagner (2008).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for "
                             "Natural Gases and Other Mixtures: An Expansion "
                             "of GERG-2004",
                    "ref": "J. Chem.Eng. Data 57(11) (2012) 3032-3091",
                    "doi":  "10.1021/je300655b"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 800000.0, "rhomax": 6.06,

        "nr1": [0.11151e1, -0.27020e1, 0.83416, -0.38828, 0.13760, 0.28185e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.62037, 0.015847, -0.61726, -0.15043, -0.012982, 0.0044325],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = lemmon, GERG
    _PR = [0.0755, -22.2752]

    _surface = {
        "__doi__": {"autor": "Mulero, A., Cachadiña, I., Bautista, D.",
                    "title": "Recommended Correlations for the Surface "
                             "Tension of n-Alkanes",
                    "ref": "J. Phys. Chem. Ref. Data 50(2) (2021) 023104",
                    "doi": "10.1063/5.0048675"},
        "sigma": [0.05351], "exp": [1.252]}
    _dielectric = {
        "eq": 1,
        "a": [44.53, 0.045], "b": [286.27, 529.31], "c": [-83471, -90493],
        "Au": 29.84, "D": 2}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.84804e1, 0.28640e1, -0.37414e1, -0.57479e1, 0.18799e1],
        "t": [1.0, 1.5, 2.3, 4.6, 5.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [-0.43785, 0.37240e1, -0.23029e1, 0.18270e1, 0.38664],
        "t": [0.116, 0.32, 0.54, 0.8, 3.5]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.33199e1, -0.23900e1, -0.15307e2, -0.51788e2, -0.11133e3],
        "t": [0.461, 0.666, 2.12, 5.1, 11.0]}

    visco0 = {"__name__": "Huber (2004)",
              "__doi__": {
                  "autor": "Huber, M.L., Laesecke, A. Xiang, H.W.",
                  "title": "Viscosity correlations for minor constituent "
                           "fluids in natural gas: n-octane, n-nonane and "
                           "n-decane",
                  "ref": "Fluid Phase Equilibria 224 (2004) 263-270",
                  "doi": "10.1016/j.fluid.2004.07.012"},

              "eq": 1, "omega": 1,

              "ek": 472.127, "sigma": 0.66383,
              "n_chapman": 0.021357,
              "collision": [0.340344, -0.466455],

              "Tref_virial": 472.127,
              "n_virial": [-19.572881, 219.73999, -1015.3226, 2471.0125,
                           -3375.1717, 2491.6597, -787.26086, 14.085455,
                           -0.34664158],
              "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],

              "Tref_res": 594.55, "rhoref_res": 232.1417, "muref_res": 1000,
              "nr": [-0.314367e-1, 0.639384e-2, 0.326258e-1, -0.108922e-1],
              "tr": [1, 1, 2, 2],
              "dr": [2, 3, 2, 3],

              "CPf": 192.935,
              "CPg1": 2.66987,
              "CPgi": [1.32137/2.66987],
              "CPti": [-0.5]}

    _viscosity = visco0,

    thermo0 = {"__name__": "Huber (2005)",
               "__doi__": {
                   "autor": "Huber, M.L., Perkins, R.A.",
                   "title": "Thermal conductivity correlations for minor "
                            "constituent fluids in natural gas: n-octane, "
                            "n-nonane and n-decane",
                   "ref": "Fluid Phase Equilibria 227 (2005) 47-55",
                   "doi": "10.1016/j.fluid.2004.10.031"},

               "eq": 1,

               "Toref": 594.55, "koref": 1,
               "no": [0.878765e-2, -0.413510e-1, 0.104791, -0.320032e-1],
               "to": [0, 1, 2, 3],

               "Tref_res": 594.55, "rhoref_res": 232.1417, "kref_res": 1,
               "nr": [0.490088e-2, 0.996486e-2, -0.807305e-2, 0.557431e-2],
               "tr": [0, -1, 0, 0],
               "dr": [1, 1, 2, 3],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03, "Xio": 0.194e-9,
               "gam0": 0.0496, "qd": 1.043054e-9, "Tcref": 891.825}

    _thermal = thermo0,


class Test(TestCase):
    def test_shortLemmon(self):
        # Table 10, Pag 842
        st = nC9(T=596, rhom=1)
        self.assertEqual(round(st.P.kPa, 3), 2200.687)
        self.assertEqual(round(st.hM.kJkmol, 3), 81692.218)
        self.assertEqual(round(st.sM.kJkmolK, 3), 156.217)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 379.897)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 715.553)
        self.assertEqual(round(st.w, 3), 85.318)

    def test_viscoHuber(self):
        # Section 3.2 pag 267
        self.assertEqual(round(nC9(T=300, P=1e7).mu.muPas, 2), 709.84)

    def test_thermoHuber(self):
        # Section 3.2 pag 53
        self.assertEqual(round(nC9(T=300, P=1e7).k.mWmK, 2), 130.31)
