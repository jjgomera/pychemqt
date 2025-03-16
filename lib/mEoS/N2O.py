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
from lib.mEoS.N2 import N2


class N2O(MEoS):
    """Multiparameter equation of state for nitrous oxide"""
    name = "nitrous oxide"
    CASNumber = "10024-97-2"
    formula = "N2O"
    synonym = "R-744A"
    _refPropName = "N2O"
    _coolPropName = "NitrousOxide"
    rhoc = unidades.Density(452.011456)
    Tc = unidades.Temperature(309.52)
    Pc = unidades.Pressure(7245.0, "kPa")
    M = 44.0128  # g/mol
    Tt = unidades.Temperature(182.33)
    Tb = unidades.Temperature(184.68)
    f_acent = 0.1613
    momentoDipolar = unidades.DipoleMoment(0.1608, "Debye")
    id = 110

    Fi1 = {"ao_log": [1, 2.5],
           "pow": [0, 1],
           "ao_pow": [-4.4262736272, 4.3120475243],
           "ao_exp": [2.1769, 1.6145, 0.48393],
           "titao": [879/Tc, 2372/Tc, 5447/Tc]}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for nitrous oxide of "
                    "Lemmon and Span (2006)",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 "
                             "Industrial Fluids",
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi": "10.1021/je050186n"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 525.0, "Pmax": 50000.0, "rhomax": 30,

        "nr1": [0.88045, -2.4235, 0.38237, 0.068917, 0.00020367],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.13122, 0.46032, -0.0036985, -0.23263, -0.00042859, -0.042810,
                -0.023038],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    eq = (lemmon, )

    _surface = {"sigma": [0.07087], "exp": [1.204]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.69078e1, 0.26620e1, -0.22386e1, -0.38002e1, 0.76922],
        "t": [1.0, 1.5, 1.9, 4.8, 5.8]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.67919e1, -0.16069e2, 0.25632e2, -0.20755e2, 0.71963e1],
        "t": [0.47, 0.72, 1.0, 1.3, 1.6]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.1287, -77.651, .21442e3, -.47809e3, .75185e3, -.46279e3],
        "t": [0.409, 1.91, 2.33, 3.0, 3.6, 4.0]}

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

              "ek": 232.4, "sigma": 0.3828, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [0.88769, 0.0214265], "psi_d": [0, 1],
              "fint": [5.15648e-4, 2.85508e-6, -2.46391e-9],
              "fint_t": [0, 1, 2],
              "chi": [0.923824, 0.03315898], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.159e-9, "gam0": 0.057, "qd": 0.446e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""
    def test_shortLemmon(self):
        """Table 10, Pag 842"""
        st = N2O(T=311, rhom=10)
        self.assertEqual(round(st.P.kPa, 3), 7474.778)
        self.assertEqual(round(st.hM.kJkmol, 3), 13676.531)
        self.assertEqual(round(st.sM.kJkmolK, 3), 52.070)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 50.336)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 2997.404)
        self.assertEqual(round(st.w, 3), 185.945)

    def test_Huber(self):
        """Table 7, pag 266"""
        st = N2O(T=278.6, rhom=20.524)
        self.assertEqual(round(st.mu.muPas, 5), 90.96945)
        self.assertEqual(round(st.k.mWmK, 4), 100.6284)
