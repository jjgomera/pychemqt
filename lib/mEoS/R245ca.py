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


class R245ca(MEoS):
    """Multiparameter equation of state for R245ca"""
    name = "1,1,2,2,3-pentafluoropropane"
    CASNumber = "679-86-7"
    formula = "CHF2CF2CH2F"
    synonym = "R245ca"
    _refPropName = "R245CA"
    _coolPropName = "R245ca"
    rhoc = unidades.Density(525.4679248)
    Tc = unidades.Temperature(447.57)
    Pc = unidades.Pressure(3940.7, "kPa")
    M = 134.04794  # g/mol
    Tt = unidades.Temperature(191.5)
    Tb = unidades.Temperature(298.412)
    f_acent = 0.355
    momentoDipolar = unidades.DipoleMoment(1.740, "Debye")

    Fi1 = {"ao_log": [1, 7.888],
           "pow": [0, 1],
           "ao_pow": [-18.09410031, 8.996084665],
           "ao_exp": [0.8843, 5.331, 14.46],
           "titao": [865/Tc, 2830/Tc, 1122/Tc]}

    CP1 = {"ao": 8.888,
           "ao_exp": [0.8843, 5.331, 14.46], "exp": [865, 2830, 1122]}

    CP2 = {"ao": -3.8444,
           "an": [5.24008e-1, -3.74976e-4], "pow": [1, 2]}

    zhou = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-245ca of Zhou (2013)",
        "__doi__": {"autor": "Zhou, Y., Lemmon, E.W.",
                    "title": "Equation of State for the Thermodynamic "
                             "Properties of 1,1,2,2,3-Pentafluoropropane "
                             "(R-245ca)",
                    "ref": "Int. J. Thermophys., 37(3) (2016) 27",
                    "doi": "10.1007/s10765-016-2039-z"},

        "R": 8.3144621,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 550.0, "Pmax": 10000.0, "rhomax": 13,

        "nr1": [0.04489247, 1.526476, -2.408320, -0.5288088, 0.18222346],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.26, 1., 1.2, 0.67],

        "nr2": [-1.063228, -0.223149, 1.18738, -0.9772383, -0.02296938],
        "d2": [1, 3, 2, 2, 7],
        "t2": [1.92, 2., 1.5, 1.93, 1.06],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [1.364444, -0.5080666, -0.06649496, -1.128359],
        "d3": [1, 1, 3, 3],
        "t3": [0.17, 3.9, 1., 1.],
        "alfa3": [1.16, 1.1, 1.64, 13.8],
        "beta3": [2.4, 1.5, 4.2, 379],
        "gamma3": [1.265, 0.42, 0.864, 1.15],
        "epsilon3": [0.55, 0.724, 0.524, 0.857]}

    eq = (zhou, )
    _PR = [-0.0856, -17.8089]

    _surface = {"sigma": [0.069297, -0.022419], "exp": [1.2795, 3.1368]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.8807, 2.1026, -3.0768, -4.9894],
        "t": [1.0, 1.5, 2.5, 4.95]}
    _liquid_Density = {
        "eq": 1,
        "n": [4.0176, -4.7916, 7.8662, -7.1049, 3.1949],
        "t": [0.48, 1.0, 1.62, 2.3, 3.1]}
    _vapor_Density = {
        "eq": 2,
        "n": [-4.65885, -1.03328, -13.5047, -48.4225, -104.097],
        "t": [0.5, 1.09, 2.1, 5.1, 10.4]}

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

              "ek": 355.41, "sigma": 0.5131, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [0.850332, 0.185603, -0.0415622], "psi_d": [0, 1, 2],
              "fint": [7.33395e-4, 1.62265e-6], "fint_t": [0, 1],
              "chi": [1.04155, -0.0118606], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.205e-9, "gam0": 0.06, "qd": 0.624e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""

    def test_zhou(self):
        """Table 3, Pag 11"""
        st = R245ca(T=250, rhom=11.3)
        self.assertEqual(round(st.P.MPa, 6), 7.998072)
        self.assertEqual(round(st.cvM.JmolK, 3), 127.995)
        self.assertEqual(round(st.cpM.JmolK, 3), 176.671)
        self.assertEqual(round(st.w, 3), 971.213)

        st = R245ca(T=400, rhom=8)
        self.assertEqual(round(st.P.MPa, 6), 2.017058)
        self.assertEqual(round(st.cvM.JmolK, 3), 151.601)
        self.assertEqual(round(st.cpM.JmolK, 3), 226.229)
        self.assertEqual(round(st.w, 3), 326.167)

        st = R245ca(T=448, rhom=3.92)
        self.assertEqual(round(st.P.MPa, 6), 3.970007)
        self.assertEqual(round(st.cvM.JmolK, 3), 189.859)
        self.assertEqual(round(st.cpM.JmolK, 1), 23758.6)
        self.assertEqual(round(st.w, 4), 73.1890)

        st = R245ca(T=250, rhom=0)
        self.assertEqual(round(st.P.MPa, 6), 0)
        self.assertEqual(round(st.cvM.JmolK, 4), 96.4517)
        self.assertEqual(round(st.cpM.JmolK, 3), 104.766)
        self.assertEqual(round(st.w, 3), 129.781)

        st = R245ca(T=420, rhom=1)
        self.assertEqual(round(st.P.MPa, 6), 2.282553)
        self.assertEqual(round(st.cvM.JmolK, 3), 155.316)
        self.assertEqual(round(st.cpM.JmolK, 3), 218.413)
        self.assertEqual(round(st.w, 3), 112.579)

        st = R245ca(T=450, rhom=7)
        self.assertEqual(round(st.P.MPa, 6), 8.304540)
        self.assertEqual(round(st.cvM.JmolK, 3), 160.706)
        self.assertEqual(round(st.cpM.JmolK, 3), 237.245)
        self.assertEqual(round(st.w, 3), 256.123)

    def test_Huber(self):
        """Table 7, pag 266"""
        st = R245ca(T=402.8, rhom=8.056)
        self.assertEqual(round(st.mu.muPas, 4), 154.9094)
        self.assertEqual(round(st.k.mWmK, 4), 65.7106)
