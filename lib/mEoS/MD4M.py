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


class MD4M(MEoS):
    """Multiparameter equation of state for tetradecamethylhexasiloxane"""
    name = "tetradecamethylhexasiloxane"
    CASNumber = "107-52-8"
    formula = "C14H42O5Si6"
    synonym = "MD4M"
    _refPropName = "MD4M"
    _coolPropName = "MD4M"
    rhoc = unidades.Density(261.6261696)
    Tc = unidades.Temperature(653.2)
    Pc = unidades.Pressure(828.56, "kPa")
    M = 458.99328  # g/mol
    Tt = unidades.Temperature(214.15)
    Tb = unidades.Temperature(532.845)
    f_acent = 0.8
    momentoDipolar = unidades.DipoleMoment(1.308, "Debye")

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [88.1018724545, -39.5537611892],
           "ao_exp": [97.16, 69.73, 38.43],
           "titao": [610/Tc, 2480/Tc, 6400/Tc]}

    f = 8.314472
    CP1 = {"ao": -20.071/f,
           "an": [2228.5e-3/f, -1311.4e-6/f, 286.2e-9/f],
           "pow": [1, 2, 3]}

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

        "Tmin": 300, "Tmax": 750, "Pmax": 125000,

        "nr1": [0.053362183, 2.8527871, -3.8108356, -0.95254215, 0.44739021],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.3, 0.68, 0.913, 0.434],

        "nr2": [-2.5194015, -1.2945338, 0.43538523, -0.92015738, -0.054299195],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.33, 2.7, 0.61, 2.12, 1.121],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [4.6112643, -0.58630821, -0.7391977, -0.14001997, -1.8085327],
        "d3": [1, 1, 3, 2, 2],
        "t3": [1.13, 0.7, 2.55, 2.59, 1.07],
        "alfa3": [0.81, 17.3, 0.892, 0.82, 0.847],
        "beta3": [0.526, 700, 0.72, 0.056, 1.3],
        "gamma3": [1.34, 1.108, 1.19, 1.68, 0.86],
        "epsilon3": [0.977, 0.92, 0.65, 1.06, 0.659]}

    colonna = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for MD4M of Colonna  (2006).",
        "__doi__": {
            "autor": "Colonna, P., Nannan, N.R., Guardone, A., Lemmon, E.W.",
            "title": "Multiparameter Equations of State for Selected "
                     "Siloxanes",
            "ref": "Fluid Phase Equilibria, 244 (2006) 193-211",
            "doi": "10.1016/j.fluid.2006.04.015"},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",
        "Tc": 653.2, "Pc": 877, "rhoc": 1/1.65,

        "Tmin": Tt, "Tmax": 673.0, "Pmax": 30000.0, "rhomax": 2.09,

        "nr1": [1.18492421, -1.87465636, -0.65713510e-1, -0.61812689,
                0.19535804, 0.50678740e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [1.23544082, 0.49462708e-1, -0.73685283, -0.19991438,
                -0.55118673e-1, 0.28325885e-1],
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
        "sigma": [0.040798], "exp": [1.3323]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-9.6745, 1.3742, 0.71467, -12.967, -7.201],
        "t": [1.0, 1.5, 2.2, 3.2, 11.4]}
    _liquid_Density = {
        "eq": 1,
        "n": [-0.034315, 4.0324, -2.7980, 2.2450, 0.41085],
        "t": [0.14, 0.45, 0.84, 1.28, 7.92]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.316, -8.324, -184.27, 350.11, -287.08, -267.25],
        "t": [0.432, 1.085, 3.78, 4.42, 5.1, 13]}

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

              "ek": 518.7, "sigma": 0.976, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [1.45542, -0.154807], "psi_d": [0, 1],
              "fint": [0.00132], "fint_t": [0],
              "chi": [1.91993], "chi_d": [0],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.347e-9, "gam0": 0.07, "qd": 1.208e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""

    def test_thol(self):
        """Table 6, Pag. 9621"""
        st = MD4M(T=280, rhom=2.1)
        self.assertEqual(round(st.P.MPa, 7), 70.8719158)
        self.assertEqual(round(st.hM.Jmol, 3), -199382.667)
        self.assertEqual(round(st.sM.JmolK, 6), -595.997052)
        self.assertEqual(round(st.w, 5), 1346.73495)
        self.assertEqual(round(st.aM.Jmol, 4), -66252.0236)

        st = MD4M(T=420, rhom=0.0005)
        self.assertEqual(round(st.P.MPa, 10), 0.0017375886)
        self.assertEqual(round(st.hM.Jmol, 4), -46438.0435)
        self.assertEqual(round(st.sM.JmolK, 7), -75.4013910)
        self.assertEqual(round(st.w, 7), 87.2885442)
        self.assertEqual(round(st.aM.Jmol, 4), -18244.6365)

        st = MD4M(T=500, rhom=0.01)
        self.assertEqual(round(st.P.MPa, 10), 0.0391881825)
        self.assertEqual(round(st.hM.Jmol, 4), 17451.7844)
        self.assertEqual(round(st.sM.JmolK, 7), 38.3367562)
        self.assertEqual(round(st.w, 7), 90.1871428)
        self.assertEqual(round(st.aM.Jmol, 5), -5635.41188)

        st = MD4M(T=500, rhom=1.8)
        self.assertEqual(round(st.P.MPa, 6), 67.169626)
        self.assertEqual(round(st.hM.Jmol, 4), -10913.2743)
        self.assertEqual(round(st.sM.JmolK, 7), -99.7222037)
        self.assertEqual(round(st.w, 6), 982.279594)
        self.assertEqual(round(st.aM.Jmol, 5), 1631.36868)

        st = MD4M(T=650, rhom=1.5)
        self.assertEqual(round(st.P.MPa, 7), 31.6991170)
        self.assertEqual(round(st.hM.Jmol, 3), 127131.508)
        self.assertEqual(round(st.sM.JmolK, 6), 178.458722)
        self.assertEqual(round(st.w, 6), 637.354433)
        self.assertEqual(round(st.aM.Jmol, 5), -9999.40628)

    # def test_Huber(self):
    #     """Table 7, pag 266"""
    #     st = MD4M(T=587.9, rhom=1.229)
    #     self.assertEqual(round(st.mu.muPas, 4), 142.4179)
    #     self.assertEqual(round(st.k.mWmK, 4), 66.3927)
