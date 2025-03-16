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
from lib.mEoS import C3


class R142b(MEoS):
    """Multiparameter equation of state for R142b"""
    name = "1-chloro-1,1-difluoroethane"
    CASNumber = "75-68-3"
    formula = "CClF2CH3"
    synonym = "R142b"
    _refPropName = "R142B"
    _coolPropName = "R142b"
    rhoc = unidades.Density(445.997)
    Tc = unidades.Temperature(410.26)
    Pc = unidades.Pressure(4055.0, "kPa")
    M = 100.49503  # g/mol
    Tt = unidades.Temperature(142.72)
    Tb = unidades.Temperature(264.03)
    f_acent = 0.2321
    momentoDipolar = unidades.DipoleMoment(2.14, "Debye")
    id = 241

    Fi1 = {"ao_log": [1, 3.],
           "pow": [0, 1],
           "ao_pow": [-12.6016527149, 8.3160183265],
           "ao_exp": [5.0385, 6.8356, 4.0591, 2.8136],
           "titao": [473/Tc, 1256/Tc, 2497/Tc, 6840/Tc]}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-142b of Lemmon "
                    "and Span (2006)",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 "
                             "Industrial Fluids",
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi": "10.1021/je050186n"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 470.0, "Pmax": 60000.0, "rhomax": 14.44,

        "nr1": [1.0038, -2.7662, 0.42921, 0.081363, 0.00024174],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.48246, 0.75542, -0.007430, -0.41460, -0.016558, -0.10644,
                -0.021704],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    eq = (lemmon, )
    _PR = [-0.0378, -18.7463]

    _surface = {"sigma": [0.05685], "exp": [1.237]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.73074e1, 0.23186e1, -0.23278e1, -0.32761e1, 0.42103],
        "t": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.17162e2, -0.47495e2, 0.57171e2, -0.25404e2, 0.15855e1],
        "t": [0.53, 0.71, 0.9, 1.1, 2.3]}
    _vapor_Density = {
        "eq": 2,
        "n": [-.3146e1, -.65221e1, -.18006e2, -.46694e2, -.26087e1, -.1102e3],
        "t": [0.408, 1.28, 3.2, 6.6, 7.0, 14.0]}

    trnECS = {"__name__": "Huber (2003)",

              "__doi__": {
                  "autor": "Huber, M.L., Laesecke, A., Perkins, R.A.",
                  "title": "Model for the Viscosity and Thermal Conductivity "
                           "of Refrigerants, Including a New Correlation for "
                           "the Viscosity of R134a",
                  "ref": "Ind. Eng. Chem. Res., 42(13) (2003) 3163-3178",
                  "doi": "10.1021/ie0300880"},

              "eq": "ecs",

              "ref": C3,
              "visco": "visco1",

              "ek": 278.2, "sigma": 0.5362, "omega": 5,

              "psi": [0.9716, 1.9181e-2], "psi_d": [0, 1],
              "fint": [9.40725e-4, 9.88196e-7], "fint_t": [0, 1],
              "chi": [1.0749, -1.77916e-2], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.03, "Xio": 0.194e-9,
              "gam0": 0.0496, "qd": 6.15654e-10, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""
    def test_shortLemmon(self):
        """Table 10, Pag 842"""
        st = R142b(T=412, rhom=4)
        self.assertEqual(round(st.P.kPa, 3), 4165.653)
        self.assertEqual(round(st.hM.kJkmol, 3), 44982.401)
        self.assertEqual(round(st.sM.kJkmolK, 3), 170.029)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 117.705)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 3187.213)
        self.assertEqual(round(st.w, 3), 96.468)
