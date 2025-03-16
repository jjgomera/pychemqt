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


class R116(MEoS):
    """Multiparameter equation of state for R116"""
    name = "hexafluoroethane"
    CASNumber = "76-16-4"
    formula = "CF3CF3"
    synonym = "R116"
    _refPropName = "R116"
    _coolPropName = "R116"
    rhoc = unidades.Density(613.3245)
    Tc = unidades.Temperature(293.03)
    Pc = unidades.Pressure(3048.0, "kPa")
    M = 138.01182  # g/mol
    Tt = unidades.Temperature(173.1)
    Tb = unidades.Temperature(195.06)
    f_acent = 0.2566
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 236

    CP1 = {"ao": 4,
           "ao_exp": [2.4818, 7.0622, 7.9951],
           "exp": [190, 622, 1470]}

    CP2 = {"ao": 27.4009901,
           "an": [-2.6057376855e-6, 9.7501305219e-10, -6559.250418,
                  787904.9649, -34166787.86],
           "pow": [2, 3, -1.001, -2, -3]}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-116 of Lemmon "
                    "and Span (2006)",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 "
                             "Industrial Fluids",
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi": "10.1021/je050186n"},

        "R": 8.314472,
        "rhoc": 4.444,
        "cp": CP1,
        "ref": "IIR",

        "Tmin": 200, "Tmax": 425.0, "Pmax": 50000.0, "rhomax": 12.31,

        "nr1": [1.1632, -2.8123, 0.77202, -0.14331, 0.10227, 0.00024629],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.30893, -0.028499, -0.30343, -0.068793, -0.027218, 0.010665],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    kozlov = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-116 of Kozlov (1996).",
        "__doi__": {"autor": "Kozlov A.D.",
                    "title": "Private communication with Dr. Alexander D. "
                             "Kozlov, Director, VNITs SMV Russian Research "
                             "Center for Standartization Information and "
                             "Certification of Materials",
                    "ref": "",
                    "doi": ""},

        "R": 8.31451,
        "cp": CP2,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 425.0, "Pmax": 50000.0, "rhomax": 12.23,

        "nr1": [2.1775273, -5.5052198, -1.3675742, -8.1284229e-1,
                -4.0207525e-1, 2.5890073, 1.4500537, -1.0445036, 9.8965288e-1,
                -8.6794888e-1, 2.8240917e-1, 4.5154220e-2, -3.0294024e-2,
                -1.7668398e-2, 2.0592774e-3],
        "d1": [1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 5, 6, 6, 7, 8],
        "t1": [0.25, 1, 3, 4, 0.25, 1, 3.5, 1.5, 2.5, 3, 3, 1, 3, 1, 1],

        "nr2": [4.2059839, 2.1500380e-1, -1.6449561e-1, -1.2396086e-1,
                1.5814552e-1, -1.4362345e-1, 1.8637877e-2, 1.6342835e-2],
        "d2": [1, 1, 4, 4, 5, 5, 8, 4],
        "t2": [2, 5, 2, 4, 8, 10, 10, 18],
        "c2": [1, 2, 2, 2, 3, 3, 3, 4],
        "gamma2": [1]*8}

    eq = lemmon, kozlov
    _PR = [-0.2951, -15.3806]

    _surface = {"sigma": [0.047593, -0.0073402], "exp": [1.2666, 1.9892]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.73997e1, 0.22554e1, -0.23385e1, -0.35244e1, 0.40350],
        "t": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.68490e2, -0.24772e3, 0.35824e3, -0.25290e3, 0.76880e2],
        "t": [0.64, 0.79, 0.95, 1.14, 1.33]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.4135, -.14529e3, .23651e3, -.22276e3, .23103e3, -.17433e3],
        "t": [0.428, 2.0, 2.24, 3.0, 4.0, 5.0]}

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
              "thermo": "thermo0",

              "ek": 226.16, "sigma": 0.5249, "omega": 5,

              "psi": [1.21996, -6.47835e-2], "psi_d": [0, 1],
              "fint": [1.32e-3], "fint_t": [0],
              "chi": [1.1804, -5.39975e-2], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
              "Xio": 0.194e-9, "gam0": 0.0496, "qd": 5e-10, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""
    def test_shortLemmon(self):
        """Table 10, Pag 842"""
        st = R116(T=295, rhom=4)
        self.assertEqual(round(st.P.kPa, 3), 3180.336)
        self.assertEqual(round(st.hM.kJkmol, 3), 34509.528)
        self.assertEqual(round(st.sM.kJkmolK, 3), 161.389)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 120.218)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 2189.730)
        self.assertEqual(round(st.w, 3), 73.317)
