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
from lib.mEoS import R134a


class C5F12(MEoS):
    """Multiparameter equation of state for perfluropentane"""
    name = "perfluoropentane"
    CASNumber = "678-26-2"
    formula = "C5F12"
    synonym = ""
    _refPropName = "C5F12"
    _coolPropName = ""
    rhoc = unidades.Density(625.03378)
    Tc = unidades.Temperature(421)
    Pc = unidades.Pressure(2063.0, "kPa")
    M = 288.034  # g/mol
    Tt = unidades.Temperature(148.21)
    Tb = unidades.Temperature(302.453)
    f_acent = 0.436
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")

    CP1 = {"ao": 15,
           "ao_exp": [5.761, 19.37, 7.096], "exp": [485, 1026, 2009]}

    gao = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for perfluerohexane of Gao "
                    "(2021).",
        "__doi__": {
            "autor": "Gao, K., Köster, A., Thol, M., Wu, J., Lemmon, E.W.",
            "title": "Equations of State for the Thermodynamic Properties of "
                     "n-Perfluorobutane, n-Perfluoropentane, and "
                     "n-Perfluerohexane",
            "ref": "Ind. Eng. Chem. Res. 60(47) (2021) 17207-17227",
            "doi": "10.1021/acs.iecr.1c02969"},

        "R": 8.314462618,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 500, "Pmax": 10000,

        "nr1": [0.036359237, 1.1222855, -0.63721964, -1.3599906, 0.24995429],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.19, 1, 1, 0.33],

        "nr2": [-2.291332, -0.95997986, 1.2463966, -0.41785522],
        "d2": [1, 3, 2, 2],
        "t2": [1.51, 2.22, 0.97, 2.54],
        "c2": [2, 2, 1, 2],
        "gamma2": [1]*4,

        "nr3": [-0.49360387, 1.8336733, -0.27152337, -0.11731069, -0.488446],
        "d3": [1, 1, 1, 3, 2],
        "t3": [0.75, 0.375, 0.5, 1, 0.815],
        "alfa3": [1.168, 0.944, 2.28, 1.486, 1.78],
        "beta3": [2.13, 2.13, 2.13, 1.865, 1.59],
        "gamma3": [1.14, 1.185, 0.95, 1.02, 1.025],
        "epsilon3": [0.317, 0.702, 0.686, 1.257, 0.947]}

    eq = (gao,)

    _surface = {"sigma": [0.04394], "exp": [1.254]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-8.4733, 3.5899, -3.3162, -4.1966, -1.6897],
        "t": [1, 1.5, 1.97, 3.63, 11.74]}
    _liquid_Density = {
        "eq": 1,
        "n": [3.9956, -4.6464, 8.1411, -7.7694, 3.4916],
        "t": [0.464, 0.9, 1.4, 1.9, 2.6]}
    _vapor_Density = {
        "eq": 2,
        "n": [-4.456, -7.6886, -24.64, -62.48, -163.57],
        "t": [0.479, 1.51, 3.35, 6.71, 14.76]}

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
              "visco": "visco0",

              "ek": 195, "sigma": 0.736, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [0.663775, 0.38206, -0.0882706], "psi_d": [0, 1, 2],
              "fint": [0.00125], "fint_t": [0],
              "chi": [1.99279, -0.308118], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.244e-9, "gam0": 0.062, "qd": 0.765e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""
    def test_gao(self):
        """Pag. 14"""
        st = C5F12(T=250, rhom=0)
        self.assertEqual(round(st.P.MPa, 6), 0)
        self.assertEqual(round(st.cvM.JmolK, 4), 199.2576)
        self.assertEqual(round(st.cpM.JmolK, 4), 207.5721)
        self.assertEqual(round(st.w, 5), 86.70462)
        self.assertEqual(round(st.hM.kJmol, 5), 15.00556)

        st = C5F12(T=250, rhom=6.5)
        self.assertEqual(round(st.P.MPa, 5), 45.74829)
        self.assertEqual(round(st.cvM.JmolK, 4), 219.3241)
        self.assertEqual(round(st.cpM.JmolK, 4), 268.4391)
        self.assertEqual(round(st.w, 4), 769.7973)
        self.assertEqual(round(st.hM.kJmol, 5), -10.43608)
        self.assertEqual(round(st.sM.JmolK, 5), -64.99865)

        st = C5F12(T=390.0, rhom=4.2)
        self.assertEqual(round(st.P.MPa, 6), 1.496384)
        self.assertEqual(round(st.cvM.JmolK, 4), 273.9917)
        self.assertEqual(round(st.cpM.JmolK, 4), 375.3160)
        self.assertEqual(round(st.w, 4), 182.6921)
        self.assertEqual(round(st.hM.kJmol, 5), 29.10127)
        self.assertEqual(round(st.sM.JmolK, 5), 83.36270)

        st = C5F12(T=421.5, rhom=2.17)
        self.assertEqual(round(st.P.MPa, 6), 2.083314)
        self.assertEqual(round(st.cvM.JmolK, 4), 302.6768)
        self.assertEqual(round(st.cpM.JmolK, 2), 15207.46)
        self.assertEqual(round(st.w, 5), 41.76442)
        self.assertEqual(round(st.hM.kJmol, 5), 44.91926)
        self.assertEqual(round(st.sM.JmolK, 4), 121.6715)

        st = C5F12(T=410, rhom=0.3)
        self.assertEqual(round(st.P.MPa, 6), 0.841555)
        self.assertEqual(round(st.cvM.JmolK, 4), 273.0757)
        self.assertEqual(round(st.cpM.JmolK, 4), 294.6374)
        self.assertEqual(round(st.w, 5), 91.41883)
        self.assertEqual(round(st.hM.kJmol, 5), 51.90776)
        self.assertEqual(round(st.sM.JmolK, 4), 143.3826)

        st = C5F12(T=450, rhom=3)
        self.assertEqual(round(st.P.MPa, 6), 4.159190)
        self.assertEqual(round(st.cvM.JmolK, 4), 297.3112)
        self.assertEqual(round(st.cpM.JmolK, 4), 431.0872)
        self.assertEqual(round(st.w, 5), 99.09973)
        self.assertEqual(round(st.hM.kJmol, 5), 51.36937)
        self.assertEqual(round(st.sM.JmolK, 4), 134.6778)

    def test_Huber(self):
        """Table 7, pag 266"""
        st = C5F12(T=378.9, rhom=4.483)
        # self.assertEqual(round(st.mu.muPas, 4), 159.7795)
        # self.assertEqual(round(st.k.mWmK, 4), 54.6119)
        self.assertEqual(round(st.mu.muPas, 4), 159.7803)
        self.assertEqual(round(st.k.mWmK, 4), 54.6120)
