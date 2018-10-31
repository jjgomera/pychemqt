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


class EthyOxide(MEoS):
    """Multiparameter equation of state for ethylene oxide"""
    name = "ethylene oxide"
    CASNumber = "75-21-8"
    formula = "C2H4O"
    synonym = ""
    _refPropName = ""
    _coolPropName = "EthyleneOxide"
    rhoc = unidades.Density(315.8568552)
    Tc = unidades.Temperature(468.92)
    Pc = unidades.Pressure(7.3047, "MPa")
    M = 44.05256  # g/mol
    Tt = unidades.Temperature(160.65)
    Tb = unidades.Temperature(283.6)
    f_acent = 0.1974
    momentoDipolar = unidades.DipoleMoment(1.8887, "Debye")
    id = 129

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [7.2881975, -2.782872],
           "ao_exp": [6.79, 4.53, 3.68],
           "titao": [1330/Tc, 2170/Tc, 4470/Tc]}

    thol = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ethylene oxide of Thol "
                    "et al. (Corrigendum) (2015)",
        "__doi__": {"autor": "Thol, M., Rutkai, G., Köster, A., Kortmann, M., "
                             "Span, R., Vrabec, J.",
                    "title": "Corrigendum to 'Fundamental equation of state "
                             "for ethylene oxide based on a hybrid dataset",
                    "ref": "Chem. Eng. Sci. 134 (2015) 887-890",
                    "doi": "10.1016/j.ces.2015.06.020"},

        "R": 8.3144621,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 1000.0, "Pmax": 700000.0, "rhomax": 100,

        "nr1": [0.0300676, 2.1629, -2.72041, -0.53931, 0.181051],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1.0, 0.41, 0.79, 1.06, 0.58],

        "nr2": [-2.61292, -2.08004, 0.3169968, -1.6532, -0.01981719],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.0, 2.2, 0.73, 2.4, 0.97],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [3.34387, -0.950671, -0.445528, -0.005409938, -0.0638824,
                -0.093912],
        "d3": [1, 1, 3, 3, 2, 1],
        "t3": [1.87, 2.08, 2.8, 0.97, 3.15, 0.7],
        "alfa3": [1.02, 1.55, 1.44, 14, 1.63, 1.9],
        "beta3": [0.62, 1.11, 0.62, 368, 0.66, 1.87],
        "gamma3": [0.847, 0.34, 0.265, 1.13, 0.36, 1.05],
        "epsilon3": [0.705, 0.821, 0.791, 1.08, 1.64, 1.51]}

    thol2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ethylene oxide of Thol "
                    "et al. (2015)",
        "__doi__": {"autor": "Thol, M., Rutkai, G., Köster, A., Kortmann, M., "
                             "Span, R., Vrabec, J.",
                    "title": "Fundamental equation of state for ethylene oxide"
                             " based on a hybrid dataset",
                    "ref": "Chem. Eng. Sci. 121 (2015) 87-99",
                    "doi": "10.1016/j.ces.2014.07.051"},

        "R": 8.3144621,
        "cp": Fi1,
        "ref": "OTO",
        "rhoc": 7.32,

        "Tmin": Tt, "Tmax": 1000.0, "Pmax": 700000.0, "rhomax": 100,

        "nr1": [0.03805675, 1.359482, -1.83337, -0.575445, 0.153649],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1.0, 0.312, 0.86, 1.114, 0.5],

        "nr2": [-1.59813, -0.682609, 0.643696, -0.535307, -0.0187222],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.1, 1.7, 0.754, 2.5, 0.9],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [1.23884, -0.431546, -0.229587, -19.3128, -0.0528359],
        "d3": [1, 1, 3, 3, 2],
        "t3": [2.18, 3.5, 2.34, 4.33, 3.9],
        "alfa3": [1.01, 1.65, 0.896, 22, 1.73],
        "beta3": [1.12, 2.16, 0.91, 196, 0.13],
        "gamma3": [0.874, 0.617, 0.476, 1.24, 0.562],
        "epsilon3": [0.7202, 0.9110, 0.688, 0.91, 1.21]}

    eq = thol, thol2

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.74136e1, 0.19870e1, -0.66330e1, 0.71500e1, -0.47200e1],
        "t": [1, 1.5, 3.5, 4.3, 5.2]}

    # The reported liquid density ancillary equation in paper don't work, using
    # the equation used in CoolProp
    _liquid_Density = {
        "eq": 1,
        "n": [2.3014, -0.08549, 2.055, -2.883, 1.686],
        "t": [0.382, 0.93, 1.48, 2.1, 2.95]}
    # _liquid_Density = {
    #     "eq": 2,
    #     "n": [0.6610, 0.4045e1, -0.4488e1, 0.3445e1, -0.9230],
    #     "t": [0.25, 0.7, 1.2, 1.75, 2.4]}

    _vapor_Density = {
        "eq": 2,
        "n": [-0.10592e1, -0.10712e2, 0.16812e2, -0.27664e2, -0.54968e2,
              -0.20428e3],
        "t": [0.3, 0.91, 1.6, 2.26, 6.6, 17]}


class Test(TestCase):
    def test_thol(self):
        # Table 5, Pag 97
        st = EthyOxide(T=200, x=0.5, eq=1)
        self.assertEqual(round(st.P.MPa, 10), 0.0007171788)
        self.assertEqual(round(st.Liquido.rhoM, 10), 22.4762797391)
        self.assertEqual(round(st.Liquido.hM.Jmol, 5), -33442.98983)
        self.assertEqual(round(st.Liquido.sM.JmolK, 7), -122.0751209)
        self.assertEqual(round(st.Liquido.cvM.JmolK, 10), 54.1084845523)
        self.assertEqual(round(st.Liquido.cpM.JmolK, 10), 81.5266043376)
        self.assertEqual(round(st.Liquido.w, 8), 1794.54046849)
        self.assertEqual(round(st.Liquido.aM.Jmol, 8), -9027.99755817)
        self.assertEqual(round(st.Gas.rhoM, 10), 0.0004315688)
        self.assertEqual(round(st.Gas.hM.Jmol, 8), -4103.02312657)
        self.assertEqual(round(st.Gas.sM.JmolK, 10), 24.6247126168)
        self.assertEqual(round(st.Gas.cvM.JmolK, 10), 28.2762101332)
        self.assertEqual(round(st.Gas.cpM.JmolK, 10), 36.6153026835)
        self.assertEqual(round(st.Gas.w, 9), 220.943064557)
        self.assertEqual(round(st.Gas.aM.Jmol, 7), -10689.7605167)

        st = EthyOxide(T=300, x=0.5, eq=1)
        self.assertEqual(round(st.P.MPa, 10), 0.1852431635)
        self.assertEqual(round(st.Liquido.rhoM, 10), 19.5606827885)
        self.assertEqual(round(st.Liquido.hM.Jmol, 7), -25005.6597985)
        self.assertEqual(round(st.Liquido.sM.JmolK, 10), -88.0098778295)
        self.assertEqual(round(st.Liquido.cvM.JmolK, 10), 58.0568818566)
        self.assertEqual(round(st.Liquido.cpM.JmolK, 9), 89.697506934)
        self.assertEqual(round(st.Liquido.w, 8), 1152.98334771)
        self.assertEqual(round(st.Liquido.aM.Jmol, 8), 1387.83337152)
        self.assertEqual(round(st.Gas.rhoM, 10), 0.0776886235)
        self.assertEqual(round(st.Gas.hM.Jmol, 7), -298.7845167)
        self.assertEqual(round(st.Gas.sM.JmolK, 9), -5.653626890)
        self.assertEqual(round(st.Gas.cvM.JmolK, 9), 41.442653701)
        self.assertEqual(round(st.Gas.cpM.JmolK, 9), 51.838824193)
        self.assertEqual(round(st.Gas.w, 9), 254.127483231)
        self.assertEqual(round(st.Gas.aM.Jmol, 8), -987.12746629)

        st = EthyOxide(T=400, x=0.5, eq=1)
        self.assertEqual(round(st.P.MPa, 10), 2.3448898851)
        self.assertEqual(round(st.Liquido.rhoM, 10), 15.5640200379)
        self.assertEqual(round(st.Liquido.hM.Jmol, 7), -14928.2462421)
        self.assertEqual(round(st.Liquido.sM.JmolK, 10), -59.5392920534)
        self.assertEqual(round(st.Liquido.cvM.JmolK, 9), 69.046404868)
        self.assertEqual(round(st.Liquido.cpM.JmolK, 9), 117.352380777)
        self.assertEqual(round(st.Liquido.w, 9), 590.414507617)
        self.assertEqual(round(st.Liquido.aM.Jmol, 8), 8736.80963955)
        self.assertEqual(round(st.Gas.rhoM, 10), 0.9448808588)
        self.assertEqual(round(st.Gas.hM.Jmol, 7), 2699.8256174)
        self.assertEqual(round(st.Gas.sM.JmolK, 9), -15.469112405)
        self.assertEqual(round(st.Gas.cvM.JmolK, 9), 62.639070990)
        self.assertEqual(round(st.Gas.cpM.JmolK, 9), 93.32073484)
        self.assertEqual(round(st.Gas.w, 9), 238.903280941)
        self.assertEqual(round(st.Gas.aM.Jmol, 8), 6405.79274306)

        st = EthyOxide(T=500, P=1e6, eq=1)
        self.assertEqual(round(st.rhoM, 10), 0.2509683066)
        self.assertEqual(round(st.hM.Jmol, 7), 11943.4908181)
        self.assertEqual(round(st.sM.JmolK, 10), 11.6066851141)
        self.assertEqual(round(st.cvM.JmolK, 10), 67.9588531666)
        self.assertEqual(round(st.cpM.JmolK, 10), 78.0665039036)
        self.assertEqual(round(st.w, 9), 315.413932985)
        self.assertEqual(round(st.aM.Jmol, 8), 2155.58138992)

        st = EthyOxide(T=500, P=1e7, eq=1)
        self.assertEqual(round(st.rhoM, 10), 5.5466493279)
        self.assertEqual(round(st.hM.Jmol, 7), 2602.9531350)
        self.assertEqual(round(st.sM.JmolK, 10), -22.6269845211)
        self.assertEqual(round(st.cvM.JmolK, 10), 81.9472541880)
        self.assertEqual(round(st.cpM.JmolK, 10), 256.331691752)
        self.assertEqual(round(st.w, 9), 214.249497552)
        self.assertEqual(round(st.aM.Jmol, 7), 12113.5551443)
