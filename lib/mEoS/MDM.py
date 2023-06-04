#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2023, Juan José Gómez Romera <jjgomera@gmail.com>

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


class MDM(MEoS):
    """Multiparamenter equation of state for octamethyltrisiloxane"""
    name = "octamethyltrisiloxane"
    CASNumber = "107-51-7"
    formula = "C8H24O2Si3"
    synonym = "MDM"
    _refPropName = "MDM"
    _coolPropName = "MDM"
    rhoc = unidades.Density(268.22667564)
    Tc = unidades.Temperature(565.3609)
    Pc = unidades.Pressure(1437.5, "kPa")
    M = 236.53146  # g/mol
    Tt = unidades.Temperature(187.2)
    Tb = unidades.Temperature(425.63)
    f_acent = 0.524
    momentoDipolar = unidades.DipoleMoment(1.079, "Debye")
    # id = 1893

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [117.9946064218, -19.6600754238],
           "ao_exp": [28.817, 46.951, 31.054],
           "titao": [20/Tc, 1570/Tc, 4700/Tc]}

    CP2 = {"ao": 4,
           "ao_exp": [28.817, 46.951, 31.054],
           "exp": [20, 1570, 4700]}

    f = 8.314472
    CP1 = {"ao": 275.1/f,
           "ao_sinh": [612.9/f], "sinh": [1829.6],
           "ao_cosh": [413.0/f], "cosh": [802.6]}

    thol = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for "
                    "octamethyltrisiloxane of Thol (2017).",
        "__doi__": {
            "autor": "Thol, M., Dubberke, F.H., Baumhögger, E., Vrabec, J., "
                     "Span, R.",
            "title": "Speed of Sound Measuements and Fundamental Equations "
                     "State for Octamethyltrisiloxane and "
                     "Decamethyltetrasiloxane",
            "ref": "J. Chem. Eng. Data 62(9) (2017) 2633-2648",
            "doi": "10.1021/acs.jced.7b00092"},

        "R": 8.3144598,
        "cp": CP2,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 600, "Pmax": 130000,

        "nr1": [0.05039724, 1.189992, -2.468723, -0.743856, 0.4434056],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.188, 1.03, 0.7, 0.464],

        "nr2": [-1.371359, -1.529621, 0.4445898, -1.009921, -0.05903694],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.105, 1.376, 0.8, 1.8, 1.005],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [3.515188, 0.08367608, 1.646856, -0.2851917, -2.457571],
        "d3": [1, 1, 3, 2, 2],
        "t3": [0.7, 0.66, 1.138, 1.56, 1.31],
        "alfa3": [0.986, 1.715, 0.837, 1.312, 1.191],
        "beta3": [0.966, 0.237, 0.954, 0.861, 0.909],
        "gamma3": [1.25, 1.438, 0.894, 0.9, 0.899],
        "epsilon3": [0.928, 2.081, 0.282, 1.496, 0.805]}

    colonna = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for MDM of Colonna (2006).",
        "__doi__": {"autor": "Colonna, P., Nannan, N.R., Guardone, A.",
                    "title": "Multiparameter equations of state for siloxanes:"
                             " [(CH3)3-Si-O1/2]2-[O-Si-(CH3)2]i=1,…,3, and "
                             "[O-Si-(CH3)2]6",
                    "ref": "Fluid Phase Equilibria 263:115-130, 2008",
                    "doi": "10.1016/j.fluid.2007.10.001"},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 673.0, "Pmax": 30000.0, "rhomax": 3.94,

        "nr1": [1.19735372, -2.40380622, 0.3256564, -0.19971259, 0.11206277,
                0.15893999e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.51234323, -0.20660361e-1, -0.38978114, -0.1186931,
                -0.37203537e-1, 0.18359984e-1],
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
        "sigma": [0.04992], "exp": [1.465]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-8.8192, 4.0952, -4.062, -6.208, -3.212],
        "t": [1, 1.5, 1.9, 3.71, 14.6]}
    _liquid_Density = {
        "eq": 1,
        "n": [7.016, -13.924, 20.84, -16.64, 5.906],
        "t": [0.54, 0.9, 1.3, 1.73, 2.2]}
    _vapor_Density = {
        "eq": 2,
        "n": [-5.3686, -11.85, -16.64, -52.26, -125.6, -235.7],
        "t": [0.515, 4.58, 2.06, 5.25, 11.3, 21.6]}

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

              "ek": 448.9, "sigma": 0.776, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [1.46043, -0.161196], "psi_d": [0, 1],
              "fint": [0.00132], "fint_t": [0],
              "chi": [3.47746, -1.50335, 0.27515], "chi_d": [0, 1, 2],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.295e-9, "gam0": 0.064, "qd": 0.956e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""

    def test_thol(self):
        """Table 13, Pag. 2646"""
        # Reference state fail with last decimal number
        st = MDM(T=320, rhom=0.0001)
        self.assertEqual(round(st.P.MPa, 10), 0.0002659237)
        # self.assertEqual(round(st.cpM.JmolK, 4), 343.4178)
        self.assertEqual(round(st.w, 4), 107.3143)
        # self.assertEqual(round(st.hM.Jmol, 3), -2908.074)
        # self.assertEqual(round(st.sM.JmolK, 5), 28.42609)
        self.assertEqual(round(st.aM.Jmol, 2), -14663.66)

        st = MDM(T=320, rhom=3.5)
        self.assertEqual(round(st.P.MPa, 5), 25.48165)
        # self.assertEqual(round(st.cpM.JmolK, 4), 427.6187)
        self.assertEqual(round(st.w, 3), 1025.009)
        # self.assertEqual(round(st.hM.Jmol, 2), -44465.80)
        self.assertEqual(round(st.sM.JmolK, 4), -140.7280)
        # self.assertEqual(round(st.aM.Jmol, 3), -6713.301)

        st = MDM(T=500, rhom=0.05)
        self.assertEqual(round(st.P.MPa, 7), 0.1919677)
        # self.assertEqual(round(st.cpM.JmolK, 4), 464.1575)
        self.assertEqual(round(st.w, 4), 123.7605)
        # self.assertEqual(round(st.hM.Jmol, 2), 68249.50)
        self.assertEqual(round(st.sM.JmolK, 4), 149.4223)
        self.assertEqual(round(st.aM.Jmol, 2), -10301.03)

        st = MDM(T=500, rhom=3)
        self.assertEqual(round(st.P.MPa, 5), 34.98485)
        # self.assertEqual(round(st.cpM.JmolK, 4), 516.4433)
        self.assertEqual(round(st.w, 4), 752.6067)
        # self.assertEqual(round(st.hM.Jmol, 2), 42566.70)
        # self.assertEqual(round(st.sM.JmolK, 5), 66.33294)
        self.assertEqual(round(st.aM.Jmol, 3), -2261.391)

        st = MDM(T=575, rhom=2)
        self.assertEqual(round(st.P.MPa, 6), 3.511622)
        # self.assertEqual(round(st.cpM.JmolK, 4), 597.2631)
        self.assertEqual(round(st.w, 4), 204.4220)
        # self.assertEqual(round(st.hM.Jmol, 2), 81168.03)
        self.assertEqual(round(st.sM.JmolK, 4), 160.0291)
        self.assertEqual(round(st.aM.Jmol, 2), -12604.53)

    def test_Huber(self):
        """Table 7, pag 266"""
        st = MDM(T=508.8, rhom=2.366)
        self.assertEqual(round(st.mu.muPas, 4), 131.0288)
        self.assertEqual(round(st.k.mWmK, 4), 68.7734)
