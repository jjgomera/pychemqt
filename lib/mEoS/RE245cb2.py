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
from lib.mEoS import R134a


class RE245cb2(MEoS):
    """Multiparameter equation of state for RE245cb2"""
    name = "methyl-pentafluoroethyl-ether"
    CASNumber = "22410-44-2"
    formula = "CF3CF2OCH3"
    synonym = "HFE-245cb2"
    _refPropName = "RE245CB2"
    _coolPropName = ""
    rhoc = unidades.Density(499.507581544)
    Tc = unidades.Temperature(406.813)
    Pc = unidades.Pressure(2886.4, "kPa")
    M = 150.047336  # g/mol
    Tt = unidades.Temperature(250)
    Tb = unidades.Temperature(278.76)
    f_acent = 0.354
    momentoDipolar = unidades.DipoleMoment(2.785, "Debye")
    # id = 1817

    CP1 = {"ao": 10.196438,
           "ao_exp": [10.214789, 10.503071, 0.98682562],
           "exp": [814, 2031, 3040]}

    zhou = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for RE245cb2 of Zhou (2010)",
        "__doi__": {"autor": "Zhou, Y., Lemmon, E.W., Mahmoud, A.M.",
                    "title": "Equations of state for RE245cb2, RE347mcc, "
                             "RE245fa2 and R1216",
                    "ref": "Preliminary equation",
                    "doi": ""},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 400000.0, "rhomax": 10.02,

        "nr1": [0.041453162, 1.5010352, -2.3142144, -0.471412, 0.17182],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.25, 0.786, 1.32, 0.338],

        "nr2": [-0.98793, -0.392049, 0.6848583, -0.32413816, -0.02414796],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.82, 2., 1., 3., 0.766],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [0.82792487, -0.31833343, -0.11929747, -0.65010212],
        "d3": [1, 1, 3, 3],
        "t3": [1.75, 3.5, 3.86, 2.75],
        "alfa3": [1.023, 1.384, 0.998, 6.9],
        "beta3": [1.727, 1.543, 1.075, 88],
        "gamma3": [1.1, 0.64, 0.5, 1.26],
        "epsilon3": [0.713, 0.917, 0.69, 0.743]}

    eq = (zhou, )

    _surface = {
        "__doi__": {
            "autor": "Mulero, A., Cachadiña, I.",
            "title": "Recommended Correlations for the Surface Tension of "
                     "Several Fluids Included in the REFPROP Program",
            "ref": "J. Phys. Chem. Ref. Data 43(2) (2014) 023104",
            "doi": "10.1063/1.4878755"},
        "sigma": [0.04534], "exp": [1.237]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.8026, 1.8804, -2.8375, -4.3077],
        "t": [1, 1.5, 2.5, 5]}
    _liquid_Density = {
        "eq": 1,
        "n": [1.8378, 2.5311, -7.084, 18.678, -30.228, 22.985],
        "t": [0.32, 1.08, 1.9, 2.8, 3.8, 4.9]}
    _vapor_Density = {
        "eq": 2,
        "n": [-1.5224, -5.7245, -15.972, -50.473, -6.8916],
        "t": [0.286, 0.82, 2.5, 5.6, 7.3]}

    trnECS = {"__name__": "Huber (2018)",

              "__doi__": {
                  "autor": "Huber, M.L.",
                  "title": "Models for Viscosity, Thermal Conductivity, and "
                           "Surface Tension of Selected Pure Fluids as "
                           "Implemented in REFPROP v10.0",
                  "ref": "NISTIR 8209",
                  "doi": "10.6028/NIST.IR.8209"},

              "eq": "ecs",
              "ref": R134a,

              "ek": 323.05, "sigma": 0.5418, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [1.0692, -0.023595], "psi_d": [0, 1],
              "fint": [0.001129], "fint_t": [0],
              "chi": [0.96324, 0.027265], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.217e-9, "gam0": 0.057, "qd": 0.66e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""
    def test_Huber(self):
        """Table 7, pag 266"""
        st = RE245cb2(T=366.1, rhom=6.908)
        # self.assertEqual(round(st.mu.muPas, 4), 131.9492)
        self.assertEqual(round(st.mu.muPas, 4), 131.9498)
        self.assertEqual(round(st.k.mWmK, 4), 52.1205)
