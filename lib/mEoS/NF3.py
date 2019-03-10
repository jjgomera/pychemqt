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


class NF3(MEoS):
    """Multiparameter equation of state for nitrogen trifluoride"""
    name = "nitrogen trifluoride"
    CASNumber = "7783-54-2"
    formula = "NF3"
    synonym = ""
    _refPropName = "NF3"
    _coolPropName = ""
    rhoc = unidades.Density(562.47)
    Tc = unidades.Temperature(234.0)
    Pc = unidades.Pressure(4460.7, "kPa")
    M = 71.019  # g/mol
    Tt = unidades.Temperature(66.36)
    Tb = unidades.Temperature(144.138)
    f_acent = 0.126
    momentoDipolar = unidades.DipoleMoment(0.235, "Debye")
    # id = 951

    CP1 = {"ao": -7.140693612211,
           "an": [0.7427518245951e6, -0.4389825372134e5, 0.1012629224351e4,
                  0.5481339146452e-1, -0.7677196006769e-4, 0.4203630864340e-7],
           "pow": [-3, -2, -1, 1, 2, 3],
           "ao_exp": [-0.6328752997967],
           "exp": [3000],
           "ao_hyp": [], "hyp": []}

    younglove = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for nitrogen trifluoride of "
                    "Younglove (1982)",
        "__doi__": {"autor": "Younglove, B.A.",
                    "title": "Thermophysical Properties of Fluids. I. Argon, "
                             "Ethylene, Parahydrogen, Nitrogen, Nitrogen "
                             "Trifluoride, and Oxygen",
                    "ref": "J. Phys. Chem. Ref. Data, 11(Suppl. 1) (1982)",
                    "doi": ""},

        "R": 8.31441,
        "rhoc": 7.92, "Tc": 234, "Pc": 4460.7, "M": 71.019,

        "cp": CP1,
        "ref": {"Tref": 300, "Pref": 101.325, "ho": 11900, "so": 260.9},

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 50000.0, "rhomax": 26.4,
        "Pmin": 0.000186, "rhomin": 26.32,

        "gamma": -0.0056,
        "b": [None, 0.1774353868e-1, -0.5409379418, 0.3976634466e1,
              -0.5209476694e3, -0.3286322888e5, -0.5990517411e-3, 0.9217525601,
              -0.4848977075e3, -0.4235892691e7, -0.9824248063e-5, .05432235989,
              -0.1462388500e2, -0.3366180440e-2, 0.2801374599, 0.8435288597e1,
              -0.1324421452e-1, 0.1875604377e-3, 0.2959643991, -0.700997687e-2,
              0.4365820912e7, -0.1111397536e8, 0.2411866612e5, 0.3179136276e7,
              0.6166849090e2, 0.4260854720e2, 0.1090598789, -0.3340951059e2,
              0.8597429644e-4, 0.1240544214e-2, 0.1286224248e-6,
              -0.8941104276e-6, 0.3353054595e-4]}

    eq = younglove,

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.66672e1, 0.33864e1, -0.28222e1, -0.50602e1, 0.32481e1],
        "t": [1.0, 1.5, 1.7, 5.5, 7.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.22080e1, 0.35709e2, -0.92868e2, 0.66666e2, -0.93589e1],
        "t": [0.35, 2.4, 2.7, 3.0, 4.0]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.061, -8.0541, -19.619, -13.432, -32.76, -67.907],
        "t": [0.421, 1.48, 3.9, 7.0, 8.0, 15.0]}


class Test(TestCase):

    def test_younglove(self):
        # The saturation state use ancillary equation for saturation pressure
        # and densities calculated values so differ of equation values
        # Selected point from Appendix J, Pag 267, single phase region
        st = NF3(T=120, P=2e4)
        self.assertEqual(round(st.rho, 0), 1649)
        self.assertEqual(round(st.rhoM, 2), 23.22)
        self.assertEqual(round(st.uM.Jmol, 0), -8547)
        self.assertEqual(round(st.hM.Jmol, 0), -8546)
        self.assertEqual(round(st.sM.JmolK, 1), 134.6)
        self.assertEqual(round(st.cvM.JmolK, 2), 41.82)
        self.assertEqual(round(st.cpM.JmolK, 2), 70.90)
        self.assertEqual(round(st.w, 1), 894.1)

        st = NF3(T=200, P=4e4)
        self.assertEqual(round(st.rho, 3), 1.717)
        self.assertEqual(round(st.rhoM, 5), 0.02418)
        self.assertEqual(round(st.uM.Jmol, 0), 5420)
        self.assertEqual(round(st.hM.Jmol, 0), 7074)
        self.assertEqual(round(st.sM.JmolK, 1), 249.2)
        self.assertEqual(round(st.cvM.JmolK, 2), 34.57)
        self.assertEqual(round(st.cpM.JmolK, 2), 43.07)
        self.assertEqual(round(st.w, 1), 169.9)

        st = NF3(T=500, P=6e4)
        self.assertEqual(round(st.rho, 3), 1.025)
        self.assertEqual(round(st.rhoM, 5), 0.01443)
        self.assertEqual(round(st.uM.Jmol, 0), 20051)
        self.assertEqual(round(st.hM.Jmol, 0), 24208)
        self.assertEqual(round(st.sM.JmolK, 1), 296.4)
        self.assertEqual(round(st.cvM.JmolK, 2), 59.25)
        self.assertEqual(round(st.cpM.JmolK, 2), 67.59)
        self.assertEqual(round(st.w, 1), 258.4)

        st = NF3(T=68, P=1e5)
        self.assertEqual(round(st.rho, 0), 1863)
        self.assertEqual(round(st.rhoM, 2), 26.23)
        self.assertEqual(round(st.uM.Jmol, 0), -12230)
        self.assertEqual(round(st.hM.Jmol, 0), -12226)
        self.assertEqual(round(st.sM.JmolK, 2), 94.41)
        self.assertEqual(round(st.cvM.JmolK, 2), 46.54)
        self.assertEqual(round(st.cpM.JmolK, 2), 73.05)
        self.assertEqual(round(st.w, 0), 1429)

        # Reference state
        st = NF3(T=300, P=101325)
        self.assertEqual(round(st.rho, 3), 2.895)
        self.assertEqual(round(st.rhoM, 5), 0.04077)
        self.assertEqual(round(st.uM.Jmol, 0), 9415)
        self.assertEqual(round(st.hM.Jmol, 0), 11900)
        self.assertEqual(round(st.sM.JmolK, 1), 260.9)
        self.assertEqual(round(st.cvM.JmolK, 2), 45.28)
        self.assertEqual(round(st.cpM.JmolK, 2), 53.74)
        self.assertEqual(round(st.w, 1), 203.4)

        st = NF3(T=155, P=2e5)
        self.assertEqual(round(st.rho, 2), 11.73)
        self.assertEqual(round(st.rhoM, 4), 0.1652)
        self.assertEqual(round(st.uM.Jmol, 0), 3826)
        self.assertEqual(round(st.hM.Jmol, 0), 5037)
        self.assertEqual(round(st.sM.JmolK, 1), 224.5)
        self.assertEqual(round(st.cvM.JmolK, 2), 31.20)
        self.assertEqual(round(st.cpM.JmolK, 2), 42.08)
        self.assertEqual(round(st.w, 1), 146.7)

        st = NF3(T=140, P=3e5)
        self.assertEqual(round(st.rho, 0), 1558)
        self.assertEqual(round(st.rhoM, 2), 21.94)
        self.assertEqual(round(st.uM.Jmol, 0), -7129)
        self.assertEqual(round(st.hM.Jmol, 0), -7115)
        self.assertEqual(round(st.sM.JmolK, 1), 145.6)
        self.assertEqual(round(st.cvM.JmolK, 2), 40.40)
        self.assertEqual(round(st.cpM.JmolK, 2), 71.73)
        self.assertEqual(round(st.w, 1), 781.8)

        st = NF3(T=420, P=4e5)
        self.assertEqual(round(st.rho, 3), 8.164)
        self.assertEqual(round(st.rhoM, 4), 0.1150)
        self.assertEqual(round(st.uM.Jmol, 0), 15431)
        self.assertEqual(round(st.hM.Jmol, 0), 18911)
        self.assertEqual(round(st.sM.JmolK, 1), 269.1)
        self.assertEqual(round(st.cvM.JmolK, 2), 54.89)
        self.assertEqual(round(st.cpM.JmolK, 2), 63.46)
        self.assertEqual(round(st.w, 1), 237.6)

        st = NF3(T=170, P=5e5)
        self.assertEqual(round(st.rho, 0), 1407)
        self.assertEqual(round(st.rhoM, 2), 19.81)
        self.assertEqual(round(st.uM.Jmol, 0), -4914)
        self.assertEqual(round(st.hM.Jmol, 0), -4889)
        self.assertEqual(round(st.sM.JmolK, 1), 159.9)
        self.assertEqual(round(st.cvM.JmolK, 2), 40.40)
        self.assertEqual(round(st.cpM.JmolK, 2), 77.55)
        self.assertEqual(round(st.w, 1), 615.2)

        st = NF3(T=220, P=6e5)
        self.assertEqual(round(st.rho, 2), 24.77)
        self.assertEqual(round(st.rhoM, 4), 0.3487)
        self.assertEqual(round(st.uM.Jmol, 0), 5915)
        self.assertEqual(round(st.hM.Jmol, 0), 7635)
        self.assertEqual(round(st.sM.JmolK, 1), 229.9)
        self.assertEqual(round(st.cvM.JmolK, 2), 37.81)
        self.assertEqual(round(st.cpM.JmolK, 2), 48.56)
        self.assertEqual(round(st.w, 1), 170.9)

        st = NF3(T=185, P=8e5)
        self.assertEqual(round(st.rho, 2), 43.39)
        self.assertEqual(round(st.rhoM, 4), 0.6109)
        self.assertEqual(round(st.uM.Jmol, 0), 4432)
        self.assertEqual(round(st.hM.Jmol, 0), 5741)
        self.assertEqual(round(st.sM.JmolK, 1), 218.3)
        self.assertEqual(round(st.cvM.JmolK, 2), 36.33)
        self.assertEqual(round(st.cpM.JmolK, 2), 52.63)
        self.assertEqual(round(st.w, 1), 149.0)

        st = NF3(T=68, P=1e6)
        self.assertEqual(round(st.rho, 0), 1864)
        self.assertEqual(round(st.rhoM, 2), 26.24)
        self.assertEqual(round(st.uM.Jmol, 0), -12235)
        self.assertEqual(round(st.hM.Jmol, 0), -12197)
        self.assertEqual(round(st.sM.JmolK, 1), 94.3)
        self.assertEqual(round(st.cvM.JmolK, 2), 47.25)
        self.assertEqual(round(st.cpM.JmolK, 2), 73.01)
        self.assertEqual(round(st.w, 0), 1412)

        st = NF3(T=205, P=2e6)
        self.assertEqual(round(st.rho, 0), 1180)
        self.assertEqual(round(st.rhoM, 2), 16.62)
        self.assertEqual(round(st.uM.Jmol, 0), -2005)
        self.assertEqual(round(st.hM.Jmol, 0), -1884)
        self.assertEqual(round(st.sM.JmolK, 1), 175.5)
        self.assertEqual(round(st.cvM.JmolK, 2), 43.09)
        self.assertEqual(round(st.cpM.JmolK, 2), 99.23)
        self.assertEqual(round(st.w, 1), 387.4)

        st = NF3(T=500, P=3e6)
        self.assertEqual(round(st.rho, 2), 51.57)
        self.assertEqual(round(st.rhoM, 4), 0.7261)
        self.assertEqual(round(st.uM.Jmol, 0), 19732)
        self.assertEqual(round(st.hM.Jmol, 0), 23864)
        self.assertEqual(round(st.sM.JmolK, 1), 263.2)
        self.assertEqual(round(st.cvM.JmolK, 2), 59.23)
        self.assertEqual(round(st.cpM.JmolK, 2), 68.87)
        self.assertEqual(round(st.w, 1), 259.5)

        st = NF3(T=230, P=4e6)
        self.assertEqual(round(st.rho, 1), 864.4)
        self.assertEqual(round(st.rhoM, 2), 12.17)
        self.assertEqual(round(st.uM.Jmol, 1), 848.9)
        self.assertEqual(round(st.hM.Jmol, 0), 1178)
        self.assertEqual(round(st.sM.JmolK, 1), 188.9)
        self.assertEqual(round(st.cvM.JmolK, 2), 50.45)
        self.assertEqual(round(st.cpM.JmolK, 1), 300.6)
        self.assertEqual(round(st.w, 1), 171.9)

        st = NF3(T=68, P=5e6)
        self.assertEqual(round(st.rho, 0), 1867)
        self.assertEqual(round(st.rhoM, 2), 26.29)
        self.assertEqual(round(st.uM.Jmol, 0), -12255)
        self.assertEqual(round(st.hM.Jmol, 0), -12065)
        self.assertEqual(round(st.sM.JmolK, 2), 94.04)
        self.assertEqual(round(st.cvM.JmolK, 2), 50.28)
        self.assertEqual(round(st.cpM.JmolK, 2), 72.82)
        self.assertEqual(round(st.w, 0), 1338)

        st = NF3(T=272, P=6e6)
        self.assertEqual(round(st.rho, 1), 280.4)
        self.assertEqual(round(st.rhoM, 3), 3.948)
        self.assertEqual(round(st.uM.Jmol, 0), 6213)
        self.assertEqual(round(st.hM.Jmol, 0), 7732)
        self.assertEqual(round(st.sM.JmolK, 1), 214.4)
        self.assertEqual(round(st.cvM.JmolK, 2), 47.86)
        self.assertEqual(round(st.cpM.JmolK, 2), 89.87)
        self.assertEqual(round(st.w, 1), 160.7)

        st = NF3(T=70, P=7e6)
        self.assertEqual(round(st.rho, 0), 1861)
        self.assertEqual(round(st.rhoM, 2), 26.21)
        self.assertEqual(round(st.uM.Jmol, 0), -12122)
        self.assertEqual(round(st.hM.Jmol, 0), -11855)
        self.assertEqual(round(st.sM.JmolK, 2), 95.98)
        self.assertEqual(round(st.cvM.JmolK, 2), 50.28)
        self.assertEqual(round(st.cpM.JmolK, 2), 70.96)
        self.assertEqual(round(st.w, 0), 1249)

        st = NF3(T=284, P=8e6)
        self.assertEqual(round(st.rho, 1), 373.9)
        self.assertEqual(round(st.rhoM, 3), 5.265)
        self.assertEqual(round(st.uM.Jmol, 0), 6210)
        self.assertEqual(round(st.hM.Jmol, 0), 7729)
        self.assertEqual(round(st.sM.JmolK, 1), 212.8)
        self.assertEqual(round(st.cvM.JmolK, 2), 49.06)
        self.assertEqual(round(st.cpM.JmolK, 2), 98.58)
        self.assertEqual(round(st.w, 1), 169.4)

        st = NF3(T=500, P=1e7)
        self.assertEqual(round(st.rho, 1), 172.3)
        self.assertEqual(round(st.rhoM, 3), 2.426)
        self.assertEqual(round(st.uM.Jmol, 0), 19004)
        self.assertEqual(round(st.hM.Jmol, 0), 23127)
        self.assertEqual(round(st.sM.JmolK, 1), 251.8)
        self.assertEqual(round(st.cvM.JmolK, 2), 59.17)
        self.assertEqual(round(st.cpM.JmolK, 2), 71.81)
        self.assertEqual(round(st.w, 1), 266.9)

        st = NF3(T=72, P=2e7)
        self.assertEqual(round(st.rho, 0), 1868)
        self.assertEqual(round(st.rhoM, 2), 26.30)
        self.assertEqual(round(st.uM.Jmol, 0), -12047)
        self.assertEqual(round(st.hM.Jmol, 0), -11286)
        self.assertEqual(round(st.sM.JmolK, 2), 97.02)
        self.assertEqual(round(st.cvM.JmolK, 2), 53.87)
        self.assertEqual(round(st.cpM.JmolK, 2), 68.70)
        self.assertEqual(round(st.w, 0), 1047)

        st = NF3(T=300, P=3e7)
        self.assertEqual(round(st.rho, 1), 969.3)
        self.assertEqual(round(st.rhoM, 2), 13.65)
        self.assertEqual(round(st.uM.Jmol, 0), 3606)
        self.assertEqual(round(st.hM.Jmol, 0), 5804)
        self.assertEqual(round(st.sM.JmolK, 1), 198.8)
        self.assertEqual(round(st.cvM.JmolK, 2), 49.48)
        self.assertEqual(round(st.cpM.JmolK, 2), 82.92)
        self.assertEqual(round(st.w, 1), 379.0)

        st = NF3(T=80, P=4e7)
        self.assertEqual(round(st.rho, 0), 1863)
        self.assertEqual(round(st.rhoM, 2), 26.23)
        self.assertEqual(round(st.uM.Jmol, 0), -11626)
        self.assertEqual(round(st.hM.Jmol, 0), -10101)
        self.assertEqual(round(st.sM.JmolK, 1), 102.6)
        self.assertEqual(round(st.cvM.JmolK, 2), 45.38)
        self.assertEqual(round(st.cpM.JmolK, 2), 69.45)
        self.assertEqual(round(st.w, 0), 1068)

        st = NF3(T=500, P=5e7)
        self.assertEqual(round(st.rho, 1), 687.9)
        self.assertEqual(round(st.rhoM, 3), 9.686)
        self.assertEqual(round(st.uM.Jmol, 0), 16264)
        self.assertEqual(round(st.hM.Jmol, 0), 21426)
        self.assertEqual(round(st.sM.JmolK, 1), 234.2)
        self.assertEqual(round(st.cvM.JmolK, 2), 60.63)
        self.assertEqual(round(st.cpM.JmolK, 2), 78.20)
        self.assertEqual(round(st.w, 1), 380.2)
