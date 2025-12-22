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


class R1233zd(MEoS):
    """Multiparameter equation of state for R1233zd"""
    name = "1-chloro-3,3,3-trifluoroprop-1-ene"
    CASNumber = "102687-65-0"
    formula = "CHCl=CH-CF3"
    synonym = "R1233zd"
    _refPropName = "R1233ZD"
    _coolPropName = "R1233zd(E)"
    rhoc = unidades.Density(483.3579248)
    Tc = unidades.Temperature(438.86)
    Pc = unidades.Pressure(3582.8, "kPa")
    M = 130.4962  # g/mol
    Tt = unidades.Temperature(165.75)
    Tb = unidades.Temperature(291.28)
    f_acent = 0.304
    momentoDipolar = unidades.DipoleMoment(1.44, "Debye")

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [-17.291229931888, 10.404947446884],
           "ao_exp": [13.7, 7.0974],
           "titao": [761/Tc, 2870/Tc]}

    CP1 = {"ao": 4.0,
           "ao_exp": [11.765, 8.6848],
           "exp": [630, 2230]}

    akasaka = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R1233zd(E) of Akasaka"
                    "(2022).",
        "__doi__": {"autor": "Akasaka, R., Lemmon, E.W.",
                    "title": "An International Standard Formulation for "
                             "trans-1-Chloro-3,3,3-trifluroprop-1-ene "
                             "[R1233zd(E)] Covering Temperatures from the"
                             "Triple-Pont Temperature to 450K and Pressures"
                             "up to 100 MPa",
                    "ref": "J. Phys. Chem. Ref. Data 51(2) (2022) 023101",
                    "doi": "10.1063/5.0083026"},

        "R": 8.314462618,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": 200, "Tmax": 550.0, "Pmax": 100000.0, "rhomax": 11.85,

        "nr1": [0.04394257, 1.062919, -1.2879140374971, -0.8088618920845,
                0.2372427],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1.0, 0.182, 0.865, 1.0924, 0.49],

        "nr2": [-1.9403, -2.831967, 0.373421, -1.515798, -0.02755627],
        "d2": [1, 3, 2, 2, 7],
        "t2": [1.958, 2.05, 0.658, 2.051, 0.862],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [4.24023, -0.03152671, -1.366494, 2.647143, -2.325463,
                -0.2541521, 0.1330834, -0.1569217],
        "d3": [1, 2, 3, 1, 1, 1, 1, 1],
        "t3": [1.852, 1.92, 1.936, 1.515, 2.668, 1.755, 0.526, 2.98],
        "alfa3": [1.532, 0.635, 1.4056, 1.451, 1.395, 2.259, 24.3, 23.6],
        "beta3": [0.2912, 0.6245, 0.669, 0.5798, 0.4643, 2.449, 1061.4, 917.8],
        "gamma3": [1.7171, 0.63, 0.7852, 2.251, 1.821, 2.074, 1.0797, 1.084],
        "epsilon3": [0.8834, 1.386, 0.5196, 1.133, .9788, 1.166, .9244, .9372]}

    mondejar = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R1233zd(E) of Mondejar "
                    "(2015).",
        "__doi__": {"autor": "Mondejar, M.E., McLinden, M.O., Lemmon, E.W.",
                    "title": "Thermodynamic Properties of trans-1-Chloro-3,3,3"
                             "-trifluoropropene (R1233zd(E)): Vapor Pressure, "
                             "(p-ρ-T) Behavior, and Spped of Sound "
                             "Measurements, and Equation of State",
                    "ref": "J. Chem. Eng. Data 60(8) (2015) 2477-2489",
                    "doi": "10.1021/acs.jced.5b00348"},

        "R": 8.3144621,
        "cp": CP1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 550.0, "Pmax": 100000.0, "rhomax": 11.41,
        "M": 130.4944, "Tc": 439.6, "rhoc": 3.68, "Pc": 3623.7,

        "nr1": [0.0478487, 1.60644, -2.27161, -0.530687, 0.169641],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1., 0.26, 1.02, 0.7, 0.4],

        "nr2": [-1.85458, -0.321916, 0.636411, -0.121482, -0.0262755],
        "d2": [1, 3, 2, 2, 7],
        "t2": [1.46, 2.3, 0.66, 2.7, 1.19],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [2.37362, -0.901771, -0.455962, -0.602941, -0.0594311],
        "d3": [1, 1, 3, 2, 2],
        "t3": [1.62, 1.13, 1.7, 1.35, 1.5],
        "alfa3": [0.748, 1.473, 1.39, 0.86, 1.8],
        "beta3": [1.29, 1.61, 0.8, 1.34, 0.49],
        "gamma3": [0.89, 1.13, 0.7, 0.91, 1.2],
        "epsilon3": [0.508, 0.366, 0.38, 0.773, 1.17]}

    eq = akasaka, mondejar

    _surface = {
        "__doi__": {
            "autor": "Kondou, C., Nagata, R., Nii, N., Koyama, S., Higashi, Y",
            "title": "Surface Tension of low GWP refrigerants R1243zf, R1234ze"
                     "(Z) and R1233zd(E)",
            "ref": "Int. J. Refrigeration 53 (2015) 80-89",
            "doi": "10.1016/j.ijrefrig.2015.01.005"},
        "sigma": [0.06195], "exp": [1.277]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.4798, 1.5791, -1.7959, -3.6716],
        "t": [1, 1.5, 2.39, 4.53]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.61448, 2.6345, -1.316, 0.86885, 0.34071],
        "t": [0.22, 0.55, 0.94, 1.4, 6.8]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.80785, -3.2355, -2.7567, -10.863, -33.456, -69.31],
        "t": [0.24, 0.59, 1, 2.2, 4.7, 9]}

    trnECS = {"__name__": "Huber (2003)",

              "__doi__": {
                  "autor": "Huber, M.L., Laesecke, A., Perkins, R.A.",
                  "title": "Model for the Viscosity and Thermal Conductivity "
                           "of Refrigerants, Including a New Correlation for "
                           "the Viscosity of R134a",
                  "ref": "Ind. Eng. Chem. Res., 42(13) (2003) 3163-3178",
                  "doi": "10.1021/ie0300880"},

              "eq": "ecs",
              "ref": R134a,
              "visco": "visco1",

              "ek": 349.1, "sigma": 0.524, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 0.93,

              "psi": [-0.0848988, 1.22693, -0.463275, 0.0568798],
              "psi_d": [0, 1, 2, 3]}

    _viscosity = (trnECS, )

    thermo0 = {"__name__": "Perkins (2017)",
               "__doi__": {
                   "autor": "Perkins, R.A., Huber, M.L.",
                   "title": "Measurement and Correlation of the Thermal "
                            "Conductivity of trans-1-Chloro-3,3,3-"
                            "trifluoropropene (R1233zd(E))",
                   "ref": "J. Chem. Eng. Data 62(9) (2017) 2659-2665",
                   "doi": "10.1021/acs.jced.7b00106"},

               "eq": 1,

               "Toref": 439.6, "koref": 1,
               "no": [-0.140033e-1, 0.378160e-1, -0.245832e-2],
               "to": [0, 1, 2],

               "Tref_res": 439.6, "rhoref_res": 480.219, "kref_res": 1,
               "nr": [0.862816e-2, 0.914709e-3, -0.208988e-1, -0.407914e-2,
                      0.511968e-1, 0.845668e-2, -0.349076e-1, -0.108985e-1,
                      0.975727e-2, 0.538262e-2, -0.926484e-3, -0.806009e-3],
               "tr": [0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1],
               "dr": [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.02, "Xio": 0.213e-9,
               "gam0": 0.059, "qd": 0.598e-9, "Tcref": 1.5*439.6}

    _thermal = (thermo0, )


class Test(TestCase):
    """Testing"""

    def test_Akasaka(self):
        """Table 9, pag 16"""
        st = R1233zd(T=300, rhom=0)
        self.assertEqual(round(st.P.MPa, 5), 0)
        self.assertEqual(round(st.cvM.JmolK, 4), 93.7166)
        self.assertEqual(round(st.cpM.JmolK, 3), 102.031)
        self.assertEqual(round(st.w, 3), 144.257)

        st = R1233zd(T=300, rhom=10)
        self.assertEqual(round(st.P.MPa, 5), 18.20558)
        self.assertEqual(round(st.cvM.JmolK, 3), 107.988)
        self.assertEqual(round(st.cpM.JmolK, 3), 150.070)
        self.assertEqual(round(st.w, 3), 797.711)

        st = R1233zd(T=300, rhom=0.05)
        self.assertEqual(round(st.P.MPa, 7), 0.1192035)
        self.assertEqual(round(st.cvM.JmolK, 4), 94.9935)
        self.assertEqual(round(st.cpM.JmolK, 3), 105.252)
        self.assertEqual(round(st.w, 3), 138.949)

        st = R1233zd(T=400, rhom=8)
        self.assertEqual(round(st.P.MPa, 5), 10.79073)
        self.assertEqual(round(st.cvM.JmolK, 3), 122.693)
        self.assertEqual(round(st.cpM.JmolK, 3), 176.124)
        self.assertEqual(round(st.w, 3), 441.123)

        st = R1233zd(T=400, rhom=0.8)
        self.assertEqual(round(st.P.MPa, 6), 1.791900)
        self.assertEqual(round(st.cvM.JmolK, 3), 121.820)
        self.assertEqual(round(st.cpM.JmolK, 3), 171.027)
        self.assertEqual(round(st.w, 3), 116.943)

        st = R1233zd(T=439, rhom=3.8)
        self.assertEqual(round(st.P.MPa, 6), 3.591512)
        self.assertEqual(round(st.cvM.JmolK, 3), 151.871)
        self.assertEqual(round(st.cpM.JmolK, 1), 56925.2)
        self.assertEqual(round(st.w, 4), 77.5936)

    def test_Perkins(self):
        """Table 2, pag E"""
        # Include basic testing for Mondejar mEoS
        # The paper use a undocumented ecs viscosity correlation, so critical
        # enhancement differ
        st = R1233zd(T=300, rho=0, eq="mondejar")
        self.assertEqual(round(st.P.MPa, 2), 0)
        self.assertEqual(round(st.k.WmK, 6), 0.010659)

        st = R1233zd(T=300, rho=5.4411, eq="mondejar")
        self.assertEqual(round(st.P.MPa, 2), 0.10)
        self.assertEqual(round(st.k.WmK, 6), 0.010766)

        st = R1233zd(T=300, rho=1308.8, eq="mondejar")
        self.assertEqual(round(st.P.MPa, 2), 20.02)
        self.assertEqual(round(st.k.WmK, 6), 0.091403)

        st = R1233zd(T=445, rho=0, eq="mondejar")
        self.assertEqual(round(st.P.MPa, 2), 0)
        self.assertEqual(round(st.k.WmK, 6), 0.021758)

        st = R1233zd(T=445, rho=168.52, eq="mondejar")
        self.assertEqual(round(st.P.MPa, 2), 3.00)
        self.assertEqual(round(st.k.WmK, 6), 0.026587)

    def test_Huber(self):
        """Table 7, pag 266"""
        self.assertEqual(round(
            # R1233zd(T=395.6, rhom=7.561).mu.muPas, 5), 113.0913)
            R1233zd(T=395.6, rhom=7.561).mu.muPas, 5), 112.9085)
