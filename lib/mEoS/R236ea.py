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


class R236ea(MEoS):
    """Multiparameter equation of state for R236ea"""
    name = "1,1,1,2,3,3-hexafluoropropane"
    CASNumber = "431-63-0"
    formula = "CF3CHFCHF2"
    synonym = "R236ea"
    _refPropName = "R236EA"
    _coolPropName = "R236EA"
    rhoc = unidades.Density(565.)
    Tc = unidades.Temperature(412.44)
    Pc = unidades.Pressure(3420.0, "kPa")
    M = 152.0384  # g/mol
    Tt = unidades.Temperature(170.0)
    Tb = unidades.Temperature(279.322)
    f_acent = 0.369
    momentoDipolar = unidades.DipoleMoment(1.129, "Debye")
    # id = 1873

    # Using the integration constants from the corrigendum file
    # Fluid Phase Equilibria 348 (2013) 83",
    # doi: 10.1016/j.fluid.2012.12.026
    Fi1 = {"ao_log": [1, 2.762],
           "pow": [0, 1],
           "ao_pow": [-14.121424135, 10.2355589225],
           "ao_exp": [0.7762, 10.41, 12.18, 3.332],
           "titao": [144/Tc, 385/Tc, 1536/Tc, 7121/Tc]}

    CP1 = {"ao": 5.30694,
           "an": [0.03973, -1.859e-5], "pow": [1, 2]}

    rui = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R236ea of Rui (2013)",
        "__doi__": {"autor": "Rui, X., Pan, J., Wang, Y.",
                    "title": "An Equation of State for Thermodynamic "
                             "Properties of 1,1,1,2,3,3-Hexafluoropropane "
                             "(R236ea)",
                    "ref": "Fluid Phase Equilibria 341 (2013) 78-85",
                    "doi": "10.1016/j.fluid.2012.12.026"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": {"Tref": 273.15, "Pref": 1.,
                "ho": 56317.4970978844, "so": 282.8465334259},

        "Tmin": Tt, "Tmax": 420.0, "Pmax": 6000.0, "rhomax": 11,

        "nr1": [0.051074, 2.5584, -2.9180, -0.71485, 0.15534],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1., 0.264, 0.5638, 1.306, 0.2062],

        "nr2": [-1.5894, -0.784, 0.85767, -0.67235, -0.017953],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.207, 2.283, 1.373, 2.33, 0.6376],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [1.3165, -0.42023, -0.28053, -1.4134, -0.0000062617],
        "d3": [1, 1, 3, 3, 2],
        "t3": [1.08, 1.67, 3.502, 4.357, 0.6945],
        "alfa3": [1.019, 1.341, 1.034, 5.264, 24.44],
        "beta3": [1.3, 2.479, 1.068, 79.85, 49.06],
        "gamma3": [1.13, 0.6691, 0.465, 1.28, 0.8781],
        "epsilon3": [0.7119, 0.9102, 0.678, 0.7091, 1.727]}

    eq = (rui, )
    _PR = [-0.1998, -19.0730]

    _surface = {"sigma": [0.306974, -0.247277], "exp": [1.12614, 1.09899]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.9095, 2.3374, -2.6453, -5.7058],
        "t": [1, 1.5, 2.15, 4.75]}
    _liquid_Density = {
        "eq": 1,
        "n": [1.6074, 1.5021, -1.106, 0.91146],
        "t": [0.31, 0.75, 1.3, 1.9]}
    _vapor_Density = {
        "eq": 2,
        "n": [-2.7426, -6.2268, -15.109, -49.524, -114.11],
        "t": [0.376, 1.1, 2.7, 5.5, 11]}

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
              "visco": "visco1",

              "ek": 318.33, "sigma": 0.5604, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [1.19985, -0.0906827, 0.0128243], "psi_d": [0, 1, 2],
              "fint": [.0054277, -2.33425e-5, 3.46098e-8], "fint_t": [0, 1, 2],
              "chi": [0.961712, 0.0337897], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.208e-9, "gam0": 0.06, "qd": 0.636e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""
    def test_Huber(self):
        """Table 7, pag 266"""
        st = R236ea(T=371.2, rhom=7.645)
        self.assertEqual(round(st.mu.muPas, 4), 157.3024)
        self.assertEqual(round(st.k.mWmK, 4), 57.4729)
