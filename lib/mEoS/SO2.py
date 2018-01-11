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


class SO2(MEoS):
    """Multiparameter equation of state for sulfur dioxide"""
    name = "sulfur dioxide"
    CASNumber = "7446-09-5"
    formula = "SO2"
    synonym = "R-764"
    rhoc = unidades.Density(525.002841)
    Tc = unidades.Temperature(430.64)
    Pc = unidades.Pressure(7884.0, "kPa")
    M = 64.0638  # g/mol
    Tt = unidades.Temperature(197.7)
    Tb = unidades.Temperature(263.13)
    f_acent = 0.2557
    momentoDipolar = unidades.DipoleMoment(1.6, "Debye")
    id = 51

    Fi1 = {"ao_log": [1, 3.],
           "pow": [0, 1, -1],
           "ao_pow": [-4.5328346436, 4.4777967379, -0.01560057996],
           "ao_exp": [1.062, 1.9401],
           "titao": [775/Tc, 1851/Tc]}

    CP2 = {"ao": 0.4021066/8.3143*64.066,
           "an": [0.87348570e-3/8.3143*64.066, -0.45968820e-6/8.3143*64.066,
                  -0.13328400e-11/8.3143*64.066, 0.23785000e-13/8.3143*64.066],
           "pow": [1, 2, 3, 4],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for sulfur dioxide of "
                    "Lemmon and Span (2006).",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 "
                             "Industrial Fluids",
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi":  "10.1021/je050186n"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 525.0, "Pmax": 35000.0, "rhomax": 25.30,
        "Pmin": 1.66, "rhomin": 25.29,

        "nr1": [0.93061, -1.9528, -0.17467, 0.061524, 0.00017711],
        "d1": [1, 1, 1, 3, 7, ],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.21615, 0.51353, 0.010419, -0.25286, -0.054720, -0.059856,
                -0.016523],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for sulfur dioxide of Polt (1987).",
        "__doi__": {"autor": "Polt, A., Platzer, B., and Maurer, G.",
                    "title": "Parameter der thermischen Zustandsgleichung von Bender fuer 14 mehratomige reine Stoffe",
                    "ref": "Chem. Technik 22(1992)6 , 216/224",
                    "doi": ""},
        "R": 8.3143,
        "cp": CP2,
        "ref": "NBP",

        "Tmin": 273.0, "Tmax": 523.0, "Pmax": 32000.0, "rhomax": 22.91,
        "Pmin": 11.82, "rhomin": 23.0,

        "nr1": [0.789407019882, -0.170449580056e1, 0.115984637964e1,
                -0.576307837294, 0.249237283833e1, -0.518115678632e1,
                0.320766081899e1, -0.123636065893e1, 0.144419600938e-1,
                -0.15380705504, 0.386324300525, 0.292550313202, -0.372445361392,
                -0.636924333910e-1, 0.986166451596e-1, -0.216993783055e-2],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [-0.789407019882, 0.170449580056e1, -0.115984637964e1,
                -0.480876182378, 0.164910076886e1, -0.133861069604e1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [1]*6}

    eq = helmholtz1, helmholtz2

    _surface = {"sigma": [0.0803, 0.0139, -0.0114],
                "exp": [0.928, 1.57, 0.364]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.73845e1, 0.22867e1, -0.24669e1, -0.32217e1, 0.23109],
        "exp": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.17156e2, -0.60441e2, 0.81407e2, -0.51871e2, 0.16754e2],
        "exp": [0.57, 0.8, 1.0, 1.3, 1.6]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-3.3832, -7.6873, -23.614, -137.20, 1866.4, -2446.9],
        "exp": [0.424, 1.4, 3.6, 8.5, 13.0, 14.0]}


class Test(TestCase):
    def test_shortLemmon(self):
        # Table 10, Pag 842
        st = SO2(T=432, rhom=8)
        self.assertEqual(round(st.P.kPa, 3), 8052.256)
        self.assertEqual(round(st.hM.kJkmol, 3), 20821.200)
        self.assertEqual(round(st.sM.kJkmolK, 3), 56.819)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 61.478)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 4877.456)
        self.assertEqual(round(st.w, 3), 171.538)
