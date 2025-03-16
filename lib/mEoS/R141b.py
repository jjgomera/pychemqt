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


class R141b(MEoS):
    """Multiparameter equation of state for R141b"""
    name = "1,1-dichloro-1-fluoroethane"
    CASNumber = "1717-00-6"
    formula = "CCl2FCH3"
    synonym = "R141b"
    _refPropName = "R141B"
    _coolPropName = "R141b"
    rhoc = unidades.Density(458.55946)
    Tc = unidades.Temperature(477.5)
    Pc = unidades.Pressure(4212.0, "kPa")
    M = 116.94962  # g/mol
    Tt = unidades.Temperature(169.68)
    Tb = unidades.Temperature(305.20)
    f_acent = 0.2195
    momentoDipolar = unidades.DipoleMoment(2.014, "Debye")
    # id = 1633

    Fi1 = {"ao_log": [1, 3.],
           "pow": [0, 1],
           "ao_pow": [-15.5074814985, 9.1871858933],
           "ao_exp": [6.8978, 7.8157, 3.2039],
           "titao": [502/Tc, 1571/Tc, 4603/Tc]}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-141b of Lemmon "
                    "and Span (2006)",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 "
                             "Industrial Fluids",
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi": "10.1021/je050186n"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 400000.0, "rhomax": 12.56,

        "nr1": [1.1469, -3.6799, 1.3469, 0.083329, 0.00025137],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.32720, 0.46946, -0.029829, -0.31621, -0.026219, -0.078043,
                -0.020498],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    eq = (lemmon, )
    _PR = [-0.1122, -17.5406]

    _surface = {"sigma": [7.3958e-5, 0.059941], "exp": [0.066331, 1.2214]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.73784e1, 0.52955e1, -0.46639e1, -0.31122e1, -0.18972e1],
        "t": [1.0, 1.5, 1.7, 4.2, 9.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.10443e2, -0.24726e2, 0.27718e2, -0.11220e2, 0.75848],
        "t": [0.49, 0.68, 0.88, 1.1, 2.9]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.31177e1, -0.68872e1, -0.18566e2, -0.40311e2, -0.95472e1,
              -0.12482e3],
        "t": [0.398, 1.33, 3.3, 6.7, 7.0, 14.0]}

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

              "ek": 370.44, "sigma": 0.5493, "omega": 5,

              "psi": [0.921345, 4.1091e-2], "psi_d": [0, 1],
              "fint": [5.21722e-4, 2.92456e-6], "fint_t": [0, 1],
              "chi": [1.0867, -2.16469e-2], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
              "Xio": 0.194e-9, "gam0": 0.0496, "qd": 5e-10, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""
    def test_shortLemmon(self):
        """Table 10, Pag 842"""
        st = R141b(T=479, rhom=3)
        self.assertEqual(round(st.P.kPa, 3), 4267.596)
        self.assertEqual(round(st.hM.kJkmol, 3), 61757.274)
        self.assertEqual(round(st.sM.kJkmolK, 3), 214.894)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 126.963)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 1485.482)
        self.assertEqual(round(st.w, 3), 93.674)
