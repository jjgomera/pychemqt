#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2023, Juan José Gómez Romera <jjgomera@gmail.com>

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


class R236fa(MEoS):
    """Multiparameter equation of state for R236fa"""
    name = "1,1,1,3,3,3-hexafluoropropane"
    CASNumber = "690-39-1"
    formula = "CF3CH2CF3"
    synonym = "R236fa"
    _refPropName = "R236FA"
    _coolPropName = "R236FA"
    rhoc = unidades.Density(551.2912384)
    Tc = unidades.Temperature(398.07)
    Pc = unidades.Pressure(3200.0, "kPa")
    M = 152.0384  # g/mol
    Tt = unidades.Temperature(179.6)
    Tb = unidades.Temperature(271.66)
    f_acent = 0.377
    momentoDipolar = unidades.DipoleMoment(1.982, "Debye")

    Fi1 = {"ao_log": [1, 9.175],
           "pow": [0, 1],
           "ao_pow": [-17.5983849, 8.87150449],
           "ao_exp": [9.8782, 18.236, 49.934],
           "titao": [962/Tc, 2394/Tc, 5188/Tc]}

    CP1 = {"ao": 53.4662555/8.314471,
           "an": [0.228092134/8.314471, 0.352999168e-4/8.314471],
           "pow": [1, 2]}

    pan = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R236fa of Pan (2012)",
        "__doi__": {"autor": "Pan, J., Rui, X., Zhao, X., Qiu, L.",
                    "title": "An equation of state for the thermodynamic "
                             "properties of 1,1,1,3,3,3-hexafluoropropane "
                             "(HFC-236fa)",
                    "ref": "Fluid Phase Equilib., 321 (2012) 10-16",
                    "doi": "10.1016/j.fluid.2012.02.012"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 400.0, "Pmax": 70000.0, "rhomax": 11.235,

        "nr1": [0.04453255, 1.777017, -2.230519, -0.6708606, 0.1587907],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1.07, 0.222, 0.66, 1.33, 0.227],

        "nr2": [-1.425119, -0.6461628, 0.8469985, -0.5635356, -0.01535611],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.33, 1.94, 1.53, 2.65, 0.722],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [1.156362, -0.407031, -0.2172753, -1.007176, -0.00006902909],
        "d3": [1, 1, 3, 3, 2],
        "t3": [1.11, 2.31, 3.68, 4.23, 0.614],
        "alfa3": [1.02, 1.336, 1.055, 5.84, 16.2],
        "beta3": [1.42, 2.31, 0.89, 80, 108],
        "gamma3": [1.13, 0.67, 0.46, 1.28, 1.2],
        "epsilon3": [0.712, 0.91, 0.677, 0.718, 1.64]}

    outcalt = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for R-236fa of Outcalt and "
                    "McLinden (1995)",
        "__doi__": {"autor": "Outcalt, S.L. McLinden, M.O.",
                    "title": "An Equation of State for the Thermodynamic "
                             "Properties of R236fa",
                    "ref": "NIST report to sponsor under contract "
                           "N61533-94-F-0152, 1995.",
                    "doi": ""},

        "R": 8.314471,
        "cp": CP1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 40000.0, "rhomax": 11.30,

        "b": [None, -0.661121874831e-1, 0.861763902745e1, -0.233732255968e3,
              0.437486232843e5, -0.539677761508e7, -0.757588552002e-2,
              0.107379563512e2, -0.106626588551e5, -0.103047455432e6,
              -0.194868091617e-2, 0.438365228107e1, -0.111207843880e4,
              -0.263710051508, 0.477521163113e2, 0.197804035098e4,
              -0.485710898935e1, 0.144821196401, -0.221059322936e2,
              0.926270169913, 0.577920666161e7, -0.985511065626e9,
              0.197199808018e6, 0.319420123094e10, 0.792946107314e4,
              -0.693606295610e6, 0.849836259084e2, 0.209702051124e7,
              0.110600369167e1, 0.953714711849e2, -0.881815206562e-2,
              0.973194908842e1, -0.935516922205e3]}

    eq = pan, outcalt
    _PR = [-0.0413, -18.8297]

    _surface = {"sigma": [0.05389], "exp": [1.249]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.78785e1, 0.15884e1, -0.48864e1, -0.50273e1, 0.89900e1],
        "t": [1.0, 1.5, 3.1, 8.0, 10.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.12320e2, -0.27579e2, 0.40114e2, -0.35461e2, 0.13769e2],
        "t": [0.579, 0.77, 0.97, 1.17, 1.4]}
    _vapor_Density = {
        "eq": 2,
        "n": [-4.4507, -3.7583, -20.279, -268.01, 501.71, -349.17],
        "t": [0.506, 1.16, 2.8, 7.0, 8.0, 9.0]}

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

              "ek": 307.24, "sigma": 0.5644, "omega": 5,

              "psi": [1.10195, -2.94253e-2], "psi_d": [0, 1],
              "fint": [1.00946e-3, 1.21255e-6], "fint_t": [0, 1],
              "chi": [1.1627, -4.37246e-2], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
              "Xio": 0.194e-9, "gam0": 0.0496, "qd": 5e-10, "Tcref": 1.5*Tc}

    _viscosity = trnECS,
    _thermal = trnECS,


class Test(TestCase):

    def test_Pan(self):
        # State point in paragraph 5. Conclusion
        st = R236fa(T=398.1, rhom=3.63)
        self.assertEqual(round(st.P.MPa, 5), 3.19277)
        self.assertEqual(round(st.hM.Jmol, 1), 59895.4)
        self.assertEqual(round(st.sM.JmolK, 3), 237.269)
        self.assertEqual(round(st.cvM.JmolK, 3), 188.426)
        self.assertEqual(round(st.w, 4), 65.0881)
