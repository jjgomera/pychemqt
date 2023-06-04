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


class nC11(MEoS):
    """Multiparameter equation of state for n-undecane"""
    name = "undecane"
    CASNumber = "1120-21-4"
    formula = "CH3-9(CH2)-CH3"
    synonym = ""
    _refPropName = "C11"
    _coolPropName = "n-Undecane"
    rhoc = unidades.Density(236.791383074)
    Tc = unidades.Temperature(638.8)
    Pc = unidades.Pressure(1990.4, "kPa")
    M = 156.30826  # g/mol
    Tt = unidades.Temperature(247.541)
    Tb = unidades.Temperature(468.934)
    f_acent = 0.539
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 15

    Fi1 = {"ao_log": [1, -120.4274],
           "pow": [-3, -2, -1, 0, 1, 2],
           "ao_pow": [-3.515339, 28.27708, -136.8378, -46.40384, 107.1876,
                      1.419929],
           "tau*logtau": -31.81246}

    aleksandrov = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for undecane of "
                    "Aleksandrov et al. (2011)",
        "__doi__": {
            "autor": "Aleksandrov, I.S., Gerasimov, A.A., Grigor’ev, B.A.",
            "title": "Using Fundamental Equations of State for Calculating "
                     "the Thermodynamic Properties of Normal Undecane",
            "ref": "Thermal Engineering, 58(8) (2011) 691-698",
            "doi": "10.1134/S0040601511080027"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 700., "Pmax": 500000.0, "rhomax": 4.97,

        "nr1": [-0.66172706, 1.3375396, -2.5608399, 0.1067891, 0.28873614e-3,
                0.49587209e-1],
        "d1": [1, 1, 1, 3, 7, 2],
        # Typo in paper, the three first τ factor are disordered
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [0.55407101e-7, 0.99754712, 1.5774025, 0.13108354e-2,
                -0.59326961, -0.93001876e-1, -0.17960228, -0.22560853e-1],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = (aleksandrov, )
    _PR = [0.1099, -26.8035]

    _surface = {
        "__doi__": {"autor": "Mulero, A., Cachadiña, I., Bautista, D.",
                    "title": "Recommended Correlations for the Surface "
                             "Tension of n-Alkanes",
                    "ref": "J. Phys. Chem. Ref. Data 50(2) (2021) 023104",
                    "doi": "10.1063/5.0048675"},
        "sigma": [0.0556], "exp": [1.32]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-9.3961, 4.4531, -5.2658, -4.7352],
        "t": [1, 1.5, 2.2, 4.5]}
    _liquid_Density = {
        "eq": 1,
        "n": [4.5273, -7.5714, 13.920, -13.464, 5.8411],
        "t": [0.46, 0.84, 1.25, 1.7, 2.2]}
    _vapor_Density = {
        "eq": 2,
        "n": [-4.3093, -3.4358, -17.473, -58.573, -133.83],
        "t": [0.466, 1.02, 2.4, 5.3, 11.4]}

    visco0 = {"__name__": "Assael (2017)",
              "__doi__": {
                  "autor": "Assael, M.J., Papalas, T.B., Huber, M.L.",
                  "title": "Reference Correlations for the Viscosity and "
                           "Thermal Conductivity of n-Undecane",
                  "ref": "J. Phys. Chem. Ref. Data 46(3) (2017) 033103",
                  "doi": "10.1063/1.4996885"},

              "eq": 1, "omega": 0,

              "Toref": Tc,
              "no_num": [0.773488, -1.53641, 19.9976, -7.58148, 2.15143,
                         -0.261065],
              "to_num": [0, 1, 2, 3, 4, 5],
              "no_den": [0.313626, 1],
              "to_den": [0, 1],

              "Tref_res": Tc, "rhoref_res": 236.7914,
              "nr_num": [256.66394],
              "tr_num": [-0.5],
              "dr_num": [2/3],
              "nr_den": [10.351826, 6.4977736, 1, 1, -1.968383, -6.4530492],
              "tr_den": [0, -1, 0, -2, -1, 0],
              "dr_den": [0, 0, 2, 0, 1, 1]}

    _viscosity = (visco0, )

    thermo0 = {"__name__": "Assael (2017)",
               "__doi__": {
                   "autor": "Assael, M.J., Papalas, T.B., Huber, M.L.",
                   "title": "Reference Correlations for the Viscosity and "
                            "Thermal Conductivity of n-Undecane",
                   "ref": "J. Phys. Chem. Ref. Data 46(3) (2017) 033103",
                   "doi": "10.1063/1.4996885"},

               "eq": 1,

               "Toref": Tc, "koref": 1e-3,
               "no_num": [-37.3793, 767.377, -3043.34, 9056.43, -5922.11,
                          1527.46],
               "to_num": [0, 1, 2, 3, 4, 5],
               "no_den": [27.743, 27.1621, 1],
               "to_den": [0, 1, 2],

               "Tref_res": Tc, "rhoref_res": 236.7914, "kref_res": 1e-3,
               "nr": [-57.3413, 64.6731, 81.5949, -44.3965, -35.4049, 1.53679,
                      8.31716, 3.20177, -0.723814, -0.308355],
               "tr": [0, -1, 0, -1, 0, -1, 0, -1, 0, -1],
               "dr": [1, 1, 2, 2, 3, 3, 4, 4, 5, 5],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
               "Xio": 0.267e-9, "gam0": 0.059, "qd": 8.66e-10, "Tcref": 958.2}

    _thermal = (thermo0, )


class Test(TestCase):
    """Testing"""

    def test_Assael(self):
        """Table 9, pag 10, saturation state for mEoS testing"""
        st = nC11(T=300, x=0.5)
        self.assertEqual(round(st.P.MPa, 8), 6.639e-5)
        self.assertEqual(round(st.Liquido.rho, 2), 734.99)
        self.assertEqual(round(st.Gas.rho, 5), 0.00416)
        self.assertEqual(round(st.Liquido.mu.muPas, 0), 1047)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 4.83)
        self.assertEqual(round(st.Liquido.k.mWmK, 1), 133.9)
        self.assertEqual(round(st.Gas.k.mWmK, 2), 8.27)

        st = nC11(T=350, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 0.001484)
        self.assertEqual(round(st.Liquido.rho, 2), 696.87)
        self.assertEqual(round(st.Gas.rho, 4), 0.0800)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 556.6)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 5.72)
        self.assertEqual(round(st.Liquido.k.mWmK, 1), 120.6)
        self.assertEqual(round(st.Gas.k.mWmK, 2), 11.66)

        st = nC11(T=400, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 0.01289)
        self.assertEqual(round(st.Liquido.rho, 2), 657.80)
        self.assertEqual(round(st.Gas.rho, 3), 0.616)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 346.8)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 6.75)
        self.assertEqual(round(st.Liquido.k.mWmK, 1), 108.8)
        self.assertEqual(round(st.Gas.k.mWmK, 2), 15.68)

        st = nC11(T=450, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 0.06218)
        self.assertEqual(round(st.Liquido.rho, 2), 616.19)
        self.assertEqual(round(st.Gas.rho, 2), 2.73)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 234.2)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 8.03)
        self.assertEqual(round(st.Liquido.k.mWmK, 2), 98.31)
        self.assertEqual(round(st.Gas.k.mWmK, 2), 20.26)

        st = nC11(T=500, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 0.2051)
        self.assertEqual(round(st.Liquido.rho, 2), 569.74)
        self.assertEqual(round(st.Gas.rho, 2), 8.68)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 164.6)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 9.72)
        self.assertEqual(round(st.Liquido.k.mWmK, 2), 89.37)
        self.assertEqual(round(st.Gas.k.mWmK, 2), 25.43)

        st = nC11(T=550, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 0.5268)
        self.assertEqual(round(st.Liquido.rho, 2), 514.09)
        self.assertEqual(round(st.Gas.rho, 1), 23.0)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 116.5)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 12.10)
        self.assertEqual(round(st.Liquido.k.mWmK, 2), 82.19)
        self.assertEqual(round(st.Gas.k.mWmK, 2), 31.71)

        st = nC11(T=600, x=0.5)
        self.assertEqual(round(st.P.MPa, 3), 1.150)
        self.assertEqual(round(st.Liquido.rho, 2), 436.84)
        self.assertEqual(round(st.Gas.rho, 1), 59.2)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 78.7)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 16.16)
        self.assertEqual(round(st.Liquido.k.mWmK, 2), 77.25)
        self.assertEqual(round(st.Gas.k.mWmK, 2), 41.98)

        # Table 10, pag 10, basic ρ,P,T point for mEoS testing
        st = nC11(T=300, P=1e7)
        self.assertEqual(round(st.rho, 2), 742.37)
        self.assertEqual(round(st.mu.muPas, 1), 1168.8)
        self.assertEqual(round(st.k.mWmK, 1), 137.7)

        st = nC11(T=350, P=1e7)
        self.assertEqual(round(st.rho, 2), 706.72)
        self.assertEqual(round(st.mu.muPas, 1), 619.6)
        self.assertEqual(round(st.k.mWmK, 1), 125.3)

        st = nC11(T=400, P=1e7)
        self.assertEqual(round(st.rho, 2), 671.21)
        self.assertEqual(round(st.mu.muPas, 1), 389.5)
        self.assertEqual(round(st.k.mWmK, 1), 114.4)

        st = nC11(T=450, P=1e7)
        self.assertEqual(round(st.rho, 2), 634.96)
        self.assertEqual(round(st.mu.muPas, 1), 268.2)
        self.assertEqual(round(st.k.mWmK, 1), 105.1)

        st = nC11(T=500, P=1e7)
        self.assertEqual(round(st.rho, 2), 597.18)
        self.assertEqual(round(st.mu.muPas, 1), 195.1)
        self.assertEqual(round(st.k.mWmK, 1), 97.5)

        st = nC11(T=550, P=1e7)
        self.assertEqual(round(st.rho, 2), 557.06)
        self.assertEqual(round(st.mu.muPas, 1), 146.9)
        self.assertEqual(round(st.k.mWmK, 1), 91.8)

        st = nC11(T=600, P=1e7)
        self.assertEqual(round(st.rho, 2), 513.80)
        self.assertEqual(round(st.mu.muPas, 1), 113.1)
        self.assertEqual(round(st.k.mWmK, 1), 88.2)

        st = nC11(T=300, P=2.5e7)
        self.assertEqual(round(st.rho, 2), 752.25)
        self.assertEqual(round(st.mu.muPas, 1), 1366.2)
        self.assertEqual(round(st.k.mWmK, 1), 142.9)

        st = nC11(T=350, P=2.5e7)
        self.assertEqual(round(st.rho, 2), 719.42)
        self.assertEqual(round(st.mu.muPas, 1), 716.6)
        self.assertEqual(round(st.k.mWmK, 1), 131.5)

        st = nC11(T=400, P=2.5e7)
        self.assertEqual(round(st.rho, 2), 687.53)
        self.assertEqual(round(st.mu.muPas, 1), 452.0)
        self.assertEqual(round(st.k.mWmK, 1), 121.6)

        st = nC11(T=450, P=2.5e7)
        self.assertEqual(round(st.rho, 2), 656.08)
        self.assertEqual(round(st.mu.muPas, 1), 315.1)
        self.assertEqual(round(st.k.mWmK, 1), 113.4)

        st = nC11(T=500, P=2.5e7)
        self.assertEqual(round(st.rho, 2), 624.74)
        self.assertEqual(round(st.mu.muPas, 1), 233.8)
        self.assertEqual(round(st.k.mWmK, 1), 106.9)

        st = nC11(T=550, P=2.5e7)
        self.assertEqual(round(st.rho, 2), 593.38)
        self.assertEqual(round(st.mu.muPas, 1), 181.3)
        self.assertEqual(round(st.k.mWmK, 1), 102.1)

        st = nC11(T=600, P=2.5e7)
        self.assertEqual(round(st.rho, 2), 562.02)
        self.assertEqual(round(st.mu.muPas, 1), 145.2)
        self.assertEqual(round(st.k.mWmK, 1), 99.2)

        st = nC11(T=300, P=5e7)
        self.assertEqual(round(st.rho, 2), 766.49)
        self.assertEqual(round(st.mu.muPas, 1), 1744.2)
        self.assertEqual(round(st.k.mWmK, 1), 150.6)

        st = nC11(T=350, P=5e7)
        self.assertEqual(round(st.rho, 2), 736.95)
        self.assertEqual(round(st.mu.muPas, 1), 889.3)
        self.assertEqual(round(st.k.mWmK, 1), 140.5)

        st = nC11(T=400, P=5e7)
        self.assertEqual(round(st.rho, 2), 708.95)
        self.assertEqual(round(st.mu.muPas, 1), 557.2)
        self.assertEqual(round(st.k.mWmK, 1), 132.0)

        st = nC11(T=450, P=5e7)
        self.assertEqual(round(st.rho, 2), 682.06)
        self.assertEqual(round(st.mu.muPas, 1), 389.7)
        self.assertEqual(round(st.k.mWmK, 1), 124.9)

        st = nC11(T=500, P=5e7)
        self.assertEqual(round(st.rho, 2), 656.08)
        self.assertEqual(round(st.mu.muPas, 1), 291.8)
        self.assertEqual(round(st.k.mWmK, 1), 119.3)

        st = nC11(T=550, P=5e7)
        self.assertEqual(round(st.rho, 2), 630.88)
        self.assertEqual(round(st.mu.muPas, 1), 229.2)
        self.assertEqual(round(st.k.mWmK, 1), 115.3)

        st = nC11(T=600, P=5e7)
        self.assertEqual(round(st.rho, 2), 606.47)
        self.assertEqual(round(st.mu.muPas, 1), 186.5)
        self.assertEqual(round(st.k.mWmK, 1), 112.7)

        # Table 11, pag 10
        st = nC11(T=550, rho=0)
        self.assertEqual(round(st.mu.muPas, 3), 8.935)
        self.assertEqual(round(st.k.mWmK, 3), 31.153)

        st = nC11(T=550, rho=10)
        self.assertEqual(round(st.mu.muPas, 3), 10.702)
        self.assertEqual(round(st.k.mWmK, 3), 31.211)

        st = nC11(T=550, rho=600)
        self.assertEqual(round(st.mu.muPas, 2), 188.68)
        self.assertEqual(round(st.k.mWmK, 2), 104.24)

        st = nC11(T=635, rho=0)
        self.assertEqual(round(st.mu.muPas, 3), 10.252)
        self.assertEqual(round(st.k.mWmK, 3), 41.522)

        st = nC11(T=635, rho=325)
        self.assertEqual(round(st.mu.muPas, 3), 49.077)
        self.assertEqual(round(st.k.mWmK, 3), 78.669)
