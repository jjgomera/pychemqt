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


class R1336mzzZ(MEoS):
    """Multiparameter equation of state for R1336mzzZ"""
    name = "cis-1,1,1,4,4,4-hexafluorobutene"
    CASNumber = "692-49-9"
    formula = "CF3CH=CHCF3"
    synonym = "R-1336mzz(Z)"
    rhoc = unidades.Density(499.386464)
    Tc = unidades.Temperature(444.5)
    Pc = unidades.Pressure(2903, "kPa")
    M = 164.056  # g/mol
    Tt = unidades.Temperature(0)

    # Kontomaris, K.
    # HFO-1336mzz-Z: High Temperature Chemical Stability and Use as A Working
    # Fluid in Organic Rankine Cycles
    # International Refrigeration and Air Condictioning Conference. Paper 1525
    # http://docs.lib.purdue.edu/iracc/1525
    Tb = unidades.Temperature(306.55)

    f_acent = 0
    momentoDipolar = unidades.DipoleMoment(0, "Debye")

    CP1 = {"ao": 4,
           "ao_exp": [20.2, 5.275],
           "exp": [736, 2299]}

    mclinden = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R1336mzz(Z) of McLinden"
                    "(2020)",
        "__doi__": {"autor": "McLinden, M.O., Akasaka, R.",
                    "title": "Thermodynamic Properties of cis-1,1,1,4,4,4-"
                             "Hexafluorobutene [R-1336mzz(Z)]: Vapor Pressure,"
                             " (p, ρ, T) Behavior, and Speed of Sound "
                             "Measurements and Equation of State",
                    "ref": "J. Chem. Eng. Data 65(9) (2020) 4201-4214",
                    "doi": "10.1021/acs.jced.9b01198"},

        "R": 8.314462618,
        "cp": CP1,
        "ref": "IIR",

        "Tmin": 230, "Tmax": 560, "Pmax": 36000.0,

        "nr1": [0.036673095, 1.1956619, -1.8462376, -0.60599297, 0.24973833],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.26, 1, 1, 0.515],

        "nr2": [-1.2548278, -1.4389612, 0.35168887, -0.82104051, -0.031747538],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.6, 3.0, 0.74, 2.68, 0.96],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [1.0281388, 0.21094074, 0.701701, 0.24638528, -1.5295034,
                0.33424978, 1.011324, -0.023457179],
        "d3": [1, 1, 3, 2, 3, 2, 2, 1],
        "t3": [1.06, 3.4, 1.617, 1.865, 1.737, 3.29, 1.242, 2],
        "alfa3": [0.746, 2.406, 0.7804, 1.25, 0.6826, 1.677, 1.762, 21],
        "beta3": [1.118, 3.065, 0.7274, 0.8435, 0.6754, 0.436, 3.808, 1888],
        "gamma3": [0.962, 1.111, 1.135, 1.163, 0.969, 1.286, 1.274, 1.056],
        "epsilon3": [1.225, 0.161, 1.231, 1.395, 0.9072, 0.958, 0.412, 0.944]}

    eq = (mclinden, )

    _surface = {
        "__doi__": {"autor": "Huber, M.L.",
                    "title": "Models for Viscosity, Thermal Conductivity, and "
                             "Surface Tension of Selected Pure Fluids as "
                             "Implemented in REFPROP v10.0",
                    "ref": "NISTIR 8209",
                    "doi": "10.6028/NIST.IR.8209"},
        "sigma": [0.06], "exp": [1.22]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.9009, 1.5186, -2.5303, -1.5139],
        "t": [1.0, 1.5, 2.8, 5]}
    _liquid_Density = {
        "eq": 1,
        "n": [2.215826, 0.80344, -1.5691755, 1.677364],
        "t": [0.3517, 1, 2, 2.24163]}
    _vapor_Density = {
        "eq": 2,
        "n": [-2.5907, -6.55386, -82.16925, -22.4756],
        "t": [0.3608, 1.0607, 6.9789, 3.00187]}

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

              "ek": 352.97, "sigma": 0.5582, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 0.94,

              "psi": [0.615513, 0.281144, -0.0527921], "psi_d": [0, 1, 2],
              "fint": [0.00109396, 0.675562e-6], "fint_t": [0, 1],
              "chi": [1.09323, -0.0316036], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.221e-9, "gam0": 0.058, "qd": 0.681e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )

    thermo0 = {"__name__": "Perkins (2020)",
               "__doi__": {
                   "autor": "Perkins, R.A., Huber, M.L.",
                   "title": "Measurement and Correlation of the Thermal "
                            "Conductivity of cis-1,1,1,4,4,4-hexafluoro-2-"
                            "butene",
                   "ref": "Int. J. Thermophysics 41 (2020) 103",
                   "doi": "10.1007/s10765-020-02681-0"},

               "eq": 1,

               "special": "_thermo0",

               "Tref_res": Tc, "rhoref_res": rhoc, "kref_res": 1,
               "nr": [0.0182670, 0.105672, -0.119719, 0.0431845, -0.00497023,
                      -0.0138014, -0.0756253, .0919892, -0.0327931, .00384035],
               "tr": [0, 0, 0, 0, 0, -1, -1, -1, -1, -1],
               "dr": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
               "Xio": 0.221e-9, "gam0": 0.058, "qd": 0.955e-9, "Tcref": 666.75}

    _thermal = (thermo0, trnECS)

    def _thermo0(self, rho, T, fase):
        """Custom dilute-gas limit thermal conductivity for Perkins method"""
        Tr = T/self.Tc
        ai = [0.033207, -0.089513, 0.101275]

        # Eq 4
        ko = Tr**0.5/(sum(a/Tr**(0.5*j) for j, a in enumerate(ai)))
        return ko*1e-3


class Test(TestCase):
    """Test class"""

    def test_mclinden(self):
        """Table S.1, pag 2"""
        st = R1336mzzZ(T=350, rhom=0)
        self.assertEqual(round(st.P.MPa, 6), 0)
        self.assertEqual(round(st.cvM.kJmolK, 6), 0.145277)
        self.assertEqual(round(st.cpM.kJmolK, 6), 0.153591)
        self.assertEqual(round(st.w, 3), 136.943)

        st = R1336mzzZ(T=350, rhom=0.1)
        self.assertEqual(round(st.P.MPa, 7), 0.2662916)
        self.assertEqual(round(st.cvM.kJmolK, 6), 0.149683)
        self.assertEqual(round(st.cpM.kJmolK, 6), 0.162997)
        self.assertEqual(round(st.w, 3), 126.681)

        st = R1336mzzZ(T=350, rhom=8)
        self.assertEqual(round(st.P.MPa, 5), 20.88061)
        self.assertEqual(round(st.cvM.kJmolK, 6), 0.160341)
        self.assertEqual(round(st.cpM.kJmolK, 6), 0.206406)
        self.assertEqual(round(st.w, 3), 617.566)

        st = R1336mzzZ(T=400, rhom=0.5)
        self.assertEqual(round(st.P.MPa, 6), 1.201815)
        self.assertEqual(round(st.cvM.kJmolK, 6), 0.170639)
        self.assertEqual(round(st.cpM.kJmolK, 6), 0.210627)
        self.assertEqual(round(st.w, 3), 108.203)

        st = R1336mzzZ(T=400, rhom=7)
        self.assertEqual(round(st.P.MPa, 5), 11.82407)
        self.assertEqual(round(st.cvM.kJmolK, 6), 0.169977)
        self.assertEqual(round(st.cpM.kJmolK, 6), 0.222507)
        self.assertEqual(round(st.w, 3), 425.612)

        st = R1336mzzZ(T=445, rhom=3.1)
        self.assertEqual(round(st.P.MPa, 6), 2.930174)
        self.assertEqual(round(st.cvM.kJmolK, 6), 0.219506)
        self.assertEqual(round(st.cpM.kJmolK, 4), 19.8479)
        self.assertEqual(round(st.w, 4), 60.4525)

    def test_Huber(self):
        """Table 7, pag 266"""
        st = R1336mzzZ(T=400, rhom=6.422, thermal=1)
        self.assertEqual(round(st.mu.muPas, 4), 121.0337)
        self.assertEqual(round(st.k.mWmK, 4), 53.5265)

    def test_Perkins(self):
        """Table 3, pag 103"""
        self.assertEqual(round(R1336mzzZ(T=300, rho=0).k.WmK, 6), 0.011056)
        self.assertEqual(round(R1336mzzZ(T=300, rho=1363).k.WmK, 6), 0.078701)
        self.assertEqual(round(R1336mzzZ(T=450, rho=0).k.WmK, 6), 0.022723)
        self.assertEqual(round(R1336mzzZ(T=450, rho=500).k.WmK, 6), 0.053990)

    def test_Surface(self):
        """Table 10, pag 271"""
        self.assertEqual(round(R1336mzzZ(T=400, x=0.5).sigma, 7), 0.0036203)
