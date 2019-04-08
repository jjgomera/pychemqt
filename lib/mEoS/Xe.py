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


class Xe(MEoS):
    """Multiparameter equation of state for xenon"""
    name = "xenon"
    CASNumber = "7440-63-3"
    formula = "Xe"
    synonym = ""
    _refPropName = "XENON"
    _coolPropName = "Xenon"
    rhoc = unidades.Density(1102.8612)
    Tc = unidades.Temperature(289.733)
    Pc = unidades.Pressure(5842.0, "kPa")
    M = 131.293  # g/mol
    Tt = unidades.Temperature(161.405)
    Tb = unidades.Temperature(165.05)
    f_acent = 0.00363
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    # id = 994

    Fi1 = {"ao_log": [1, 1.5],
           "pow": [0, 1],
           "ao_pow": [-3.8227178129, 3.8416395351],
           "ao_exp": [], "titao": []}

    CP1 = {"ao": 2.5,
           "an": [], "pow": [],
           "ao_exp": [], "exp": []}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for xenon of Lemmon "
                    "and Span (2006).",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 "
                             "Industrial Fluids",
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi":  "10.1021/je050186n"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 750.0, "Pmax": 700000.0, "rhomax": 28.78,
        "Pmin": 81.77, "rhomin": 22.59,

        "nr1": [0.83115, -2.3553, 0.53904, 0.014382, 0.066309, 0.00019649],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.14996, -0.035319, -0.15929, -0.027521, -0.023305, 0.0086941],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6,

        "nr3": [],
        "nr4": []}

    eq = lemmon,

    _surface = {"sigma": [-0.11538, 0.16598], "exp": [1.0512, 1.098]}
    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [],  "expt0": [], "expd0": [],
                   "a1": [10.122], "expt1": [0], "expd1": [1],
                   "a2": [31.97, 46.97, -82.51, -948.4],
                   "expt2": [0, 1, 0], "expd2": [2, 2, 2.7]}

    _melting = {
        "eq": 2,
        "__doi__": {"autor": "Michels, A., Prins, C.",
                    "title": "The Melting Lines of Argon, Krypton and Xenon "
                             "up to 1500 Atm; Representation of the Results "
                             "by a Law of Corresponding States",
                    "ref": "Physica 28 (1962) 101-116",
                    "doi": "10.1016/0031-8914(62)90096-4"},

        "Tmin": Tt, "Tmax": 1300.0,
        "Tref": 1, "Pref": -2576*101325,
        "a1": [0.79932770279654*101325], "exp1": [1.589165]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.60231e1, 0.14989e1, -0.74906, -0.12194e1, -0.44905],
        "t": [1., 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.13570e2, -0.47545e2, 0.63876e2, -0.39983e2, 0.12701e2],
        "t": [0.56, 0.8, 1.0, 1.3, 1.6]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.0026, -6.056, -0.60339e2, 0.48838e3, -0.81974e3, 0.47287e3],
        "t": [0.435, 1.4, 4.4, 6.2, 7.0, 8.6]}


class Test(TestCase):
    def test_shortLemmon(self):
        # Table 10, Pag 842
        st = Xe(T=291, rhom=8)
        self.assertEqual(round(st.P.kPa, 3), 5986.014)
        self.assertEqual(round(st.hM.kJkmol, 3), 9193.668)
        self.assertEqual(round(st.sM.kJkmolK, 3), 36.895)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 28.692)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 3063.309)
        self.assertEqual(round(st.w, 3), 125.648)

    def test_Michels(self):
        # Table I, pag 105
        self.assertEqual(round(Xe._Melting_Pressure(161.554).atm, 2), 3.66)
        self.assertEqual(round(Xe._Melting_Pressure(167.154).atm, 2), 147.21)
        self.assertEqual(round(Xe._Melting_Pressure(191.144).atm, 2), 794.08)
        self.assertEqual(round(Xe._Melting_Pressure(215.264).atm, 2), 1494.59)
