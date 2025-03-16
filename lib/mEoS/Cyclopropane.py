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


class Cyclopropane(MEoS):
    """Multiparameter equation of state for cyclopropane"""
    name = "cyclopropane"
    CASNumber = "75-19-4"
    formula = "cyclo(CH2)3"
    synonym = ""
    _refPropName = "CYCLOPRO"
    _coolPropName = "CycloPropane"
    rhoc = unidades.Density(258.5)
    Tc = unidades.Temperature(398.3)
    Pc = unidades.Pressure(5579.7, "kPa")
    M = 42.081  # g/mol
    Tt = unidades.Temperature(145.7)
    Tb = unidades.Temperature(241.67)
    f_acent = 0.1305
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 258

    CP1 = {"ao": 1.26016/8.3143*42.081,
           "an": [-0.90530700e-2/8.3143*42.081, 0.50550400e-4/8.3143*42.081,
                  -0.77223700e-7/8.3143*42.081, 0.40538000e-10/8.3143*42.081],
           "pow": [1, 2, 3, 4]}

    # FIXME: Bad saturation point
    polt = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for cyclopropane of Polt et "
                    "al. (1992)",
        "__doi__": {"autor": "Polt, A., Platzer, B., Maurer, G.",
                    "title": "Parameter der thermischen Zustandsgleichung von "
                             "Bender fuer 14 mehratomige reine Stoffe",
                    "ref": "Chem. Technik 22(1992)6 , 216/224",
                    "doi": ""},

        "R": 8.3143,
        "cp": CP1,
        "ref": "IIR",

        "Tmin": 273, "Tmax": 473.0, "Pmax": 28000.0, "rhomax": 15.595,

        "nr1": [-0.137016097588e1, 0.212444673002e1, -0.578908942724,
                -0.115633726379e1, 0.252574014413e1, -0.282265442929e1,
                0.283576113255, -0.842718450726e-1, 0.931086305879,
                -0.105296584292e1, 0.432020532920, -0.251108254803,
                0.127725582443, 0.483621161849e-1, -0.116473795607e-1,
                0.334005754773e-3],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.137016097588e1, -0.212444673002e1, 0.578908942724,
                0.304945770499, -0.184276165165, -0.292111460397],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [1]*6}

    eq = (polt, )

    _surface = {
        "__doi__": {
            "autor": "Mulero, A., Cachadiña, I.",
            "title": "Recommended Correlations for the Surface Tension of "
                     "Several Fluids Included in the REFPROP Program",
            "ref": "J. Phys. Chem. Ref. Data 43(2) (2014) 023104",
            "doi": "10.1063/1.4878755"},
        "sigma": [0.06812], "exp": [1.314]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.73438e1, 0.17584e2, -0.34265e2, 0.20155e2, -0.77259e1],
        "t": [1.0, 1.5, 1.71, 1.95, 4.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.16998, 0.35101e1, -0.27092e1, 0.17644e1],
        "t": [0.11, 0.5, 0.8, 1.1]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.33232, -29.566, 57.762, -142.21, 325.73, -244.39],
        "t": [0.1, 0.87, 1.14, 1.78, 2.32, 2.6]}

    trnECS = {"__name__": "Huber (2018)",

              "__doi__": {
                  "autor": "Huber, M.L.",
                  "title": "Models for Viscosity, Thermal Conductivity, and "
                           "Surface Tension of Selected Pure Fluids as "
                           "Implemented in REFPROP v10.0",
                  "ref": "NISTIR 8209",
                  "doi": "10.6028/NIST.IR.8209"},

              "eq": "ecs",
              "ref": C3,
              "visco": "visco1",

              "ek": 316.29, "sigma": 0.442, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 0.95,

              "psi": [1], "psi_d": [0],
              "fint": [0.00122], "fint_t": [0],
              "chi": [1], "chi_d": [0],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.191e-9, "gam0": 0.057, "qd": 0.534e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""
    def test_Huber(self):
        """Table 7, pag 266"""
        st = Cyclopropane(T=358.5, rhom=12.343)
        # self.assertEqual(round(st.mu.muPas, 5), 82.90981)
        # self.assertEqual(round(st.k.mWmK, 4), 95.1952)
        self.assertEqual(round(st.mu.muPas, 5), 82.90939)
        self.assertEqual(round(st.k.mWmK, 4), 95.345)
