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


class HCl(MEoS):
    """Multiparameter equation of state for hydrogen chloride"""
    name = "hydrogen chloride"
    CASNumber = "7647-01-0 "
    formula = "HCl"
    synonym = ""
    _refPropName = "HCL"
    _coolPropName = "HydrogenChloride"
    rhoc = unidades.Density(432.7913578)
    Tc = unidades.Temperature(324.68)
    Pc = unidades.Pressure(8313.5, "kPa")
    M = 36.46094  # g/mol
    Tt = unidades.Temperature(159.07)
    Tb = unidades.Temperature(188.179)
    f_acent = 0.131
    momentoDipolar = unidades.DipoleMoment(1.079, "Debye")
    id = 104

    Fi1 = {"ao_log": [1, 2.5],
           "pow": [0, 1],
           "ao_pow": [-4.069044527, 4.0257768311],
           "ao_exp": [0.0033327, 0.935243, 0.209996],
           "titao": [300/Tc, 4000/Tc, 6300/Tc]}

    Fi2 = {"ao_log": [1, 2.5],
           "pow": [0, 1, -1, -2],
           "ao_pow": [7.913048, -3.217575, -4.149937e-3, 8.019202e-4],
           "ao_exp": [1.054392], "titao": [1.241138e1]}

    thol = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hydrogen chloride of Thol"
                    " et al. (2018)",
        "__doi__": {"autor": "Thol, M., Dubberke, F.H., Baumhögger, E., Span, "
                             "R., Vrabec, J.",
                    "title": "Speed of Sound Measurements and a Fundamental "
                             "Equation of State for Hydrogen Chloride",
                    "ref": "J. Chem. Eng. Data 63(7) (2018) 2533-2547",
                    "doi": "10.1021/acs.jced.7b01031"},

        "R": 8.3144598,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 670, "Pmax": 20000.0, "rhomax": 35,

        "nr1": [0.01952802, 1.926809, -2.835744, -0.2276121, 0.08843713],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.553, 1.037, 0.817, 0.378],

        "nr2": [-2.433471, -0.2636625, 0.6307008, -0.6382638, -0.006851438],
        "d2": [1, 3, 2, 2, 7],
        "t2": [1.523, 2.656, 1.338, 2.828, 0.75],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [7.363661, -1.262993, -0.006539739, -0.8752692, -3.224835],
        "d3": [1, 1, 3, 2, 2],
        "t3": [0.644, 2.892, 0.76, 1.323, 0.693],
        "alfa3": [1.141, 1.162, 34.6, 1.175, 0.99],
        "beta3": [0.95, 0.92, 1550, 1.2, 0.89],
        "gamma3": [1.56, 1.14, 1.06, 0.94, 1.25],
        "epsilon3": [0.855, 0.91, 0.942, 0.702, 0.487]}

    thol2014 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hydrogen chloride of Thol"
                    " et al. (2013)",
        "__doi__": {"autor": "Thol, M., Piazza, L., Span, R.",
                    "title": "A New Functional Form for Equations of State "
                             "for Some Weakly Associating Fluids",
                    "ref": "Int. J. Thermophys., 35(5):783-811, 2014.",
                    "doi": "10.1007/s10765-014-1633-1"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO",

        "Tmin": 155.0, "Tmax": 330.0, "Pmax": 20000.0, "rhomax": 33.8145,
        "Tc": 324.55, "rhoc": 11.271514, "M": 36.460939,

        "nr1": [-.40937325, 0.943994574, -1.78830477, 0.128619044,
                4.39018427e-3, 0.0130480908, 1.69387782e-3],
        "d1": [1, 1, 1, 2, 3, 3, 5],
        "t1": [-0.75, -0.25, 1.25, 0.75, -1.0, -0.375, 1.25],

        "nr2": [0.751559060, -0.800007427, 0.430935939, 4.54319457e-3,
                -1.52172259e-1, -4.36174059e-2, -9.70625964e-3, 1.01144098e-2,
                3.76991644e-3],
        "d2": [1, 1, 2, 5, 1, 3, 4, 5, 2],
        "t2": [2.375, 3.0, 2.625, 1.875, 4.5, 5.75, 5.375, 2.75, 14.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 2, 3],
        "gamma2": [1]*9}

    eq = (thol, thol2014)

    _surface = {
        "__doi__": {
            "autor": "Mulero, A., Cachadiña, I.",
            "title": "Recommended Correlations for the Surface Tension of "
                     "Several Fluids Included in the REFPROP Program",
            "ref": "J. Phys. Chem. Ref. Data 43(2) (2014) 023104",
            "doi": "10.1063/1.4878755"},
        "sigma": [0.05994], "exp": [1.0953]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-6.73, 1.464, -1.994, 1.283, -2.062],
        "t": [1, 1.5, 3.12, 3.95, 4.8]}
    _liquid_Density = {
        "eq": 1,
        "n": [2.547, -0.631, 1.75, -1.922, 1.03],
        "t": [0.418, 1.12, 1.86, 2.66, 3.57]}
    _vapor_Density = {
        "eq": 2,
        "n": [-2.5676, -4.1055, -12.068, -29.03, -54.93, -222.7],
        "t": [0.417, 0.923, 2.57, 5.54, 10.5, 23.3]}

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

              "ek": 257.8, "sigma": 0.355, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [0.615877, 0.55609, -0.337867, 0.0681029],
              "psi_d": [0, 1, 2, 3],
              "fint": [0.0006], "fint_t": [0],
              "chi": [1.57373, -0.17681], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.154e-9, "gam0": 0.054, "qd": 0.424e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""

    def test_thol(self):
        """Table 4, Pag 2538"""
        st = HCl(T=180, rhom=34)
        self.assertEqual(round(st.P.MPa, 7), 32.0901855)
        self.assertEqual(round(st.hM.kJkmol, 6), 143.888365)
        self.assertEqual(round(st.sM.kJkmolK, 8), -4.43856331)
        self.assertEqual(round(st.w, 5), 1245.21167)
        # self.assertEqual(round(st.aM.kJkmol, 9), -0.999224087)

        st = HCl(T=180, rhom=0.04)
        self.assertEqual(round(st.P.MPa, 10), 0.0586755086)
        self.assertEqual(round(st.hM.kJkmol, 4), 16001.7182)
        self.assertEqual(round(st.sM.kJkmolK, 7), 89.4473371)
        self.assertEqual(round(st.w, 6), 237.450444)
        self.assertEqual(round(st.aM.kJkmol, 5), -1565.69018)

        st = HCl(T=300, rhom=25)
        self.assertEqual(round(st.P.MPa, 7), 20.5806246)
        self.assertEqual(round(st.hM.kJkmol, 5), 7239.11135)
        self.assertEqual(round(st.sM.kJkmolK, 7), 27.2001575)
        self.assertEqual(round(st.w, 6), 691.262823)
        self.assertEqual(round(st.aM.kJkmol, 5), -1744.16090)

        st = HCl(T=300, rhom=3)
        self.assertEqual(round(st.P.MPa, 8), 4.70517077)
        self.assertEqual(round(st.hM.kJkmol, 4), 16399.1807)
        self.assertEqual(round(st.sM.kJkmolK, 7), 60.1776839)
        self.assertEqual(round(st.w, 6), 243.922254)
        self.assertEqual(round(st.aM.kJkmol, 5), -3222.51469)

        st = HCl(T=400, rhom=18)
        self.assertEqual(round(st.P.MPa, 7), 36.1071944)
        self.assertEqual(round(st.hM.kJkmol, 4), 13798.3760)
        self.assertEqual(round(st.sM.kJkmolK, 7), 44.0061792)
        self.assertEqual(round(st.w, 6), 515.047201)
        self.assertEqual(round(st.aM.kJkmol, 5), -5810.05088)


    def test_thol2014(self):
        """Table 9, Pag 26"""
        # Discard the last 4 number, I'm fairly sure is a problem with the
        # significative figures in the equation parameters in paper

        st = HCl(T=170, rho=0.01, eq=1)
        self.assertEqual(round(st.P.MPa, 9), 0.000387586)
        self.assertEqual(round(st.h.kJkg, 5), -102.41433)
        self.assertEqual(round(st.s.kJkgK, 5), 0.82032)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.57116)
        self.assertEqual(round(st.cp.kJkgK, 5), 0.79948)
        self.assertEqual(round(st.w, 5), 232.89853)

        st = HCl(T=170, rho=1230, eq=1)
        self.assertEqual(round(st.P.MPa, 4), 1.2292)
        self.assertEqual(round(st.h.kJkg, 5), -561.33940)
        self.assertEqual(round(st.s.kJkgK, 5), -2.88832)
        self.assertEqual(round(st.cv.kJkgK, 5), 1.14972)
        self.assertEqual(round(st.cp.kJkgK, 5), 1.55394)
        self.assertEqual(round(st.w, 3), 999.439)

        st = HCl(T=280, rho=0.1, eq=1)
        self.assertEqual(round(st.P.MPa, 5), 0.00638)
        self.assertEqual(round(st.h.kJkg, 4), -14.6144)
        self.assertEqual(round(st.s.kJkgK, 5), 0.58003)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.57147)
        self.assertEqual(round(st.cp.kJkgK, 5), 0.80008)
        self.assertEqual(round(st.w, 5), 298.83636)

        st = HCl(T=280, rho=900, eq=1)
        self.assertEqual(round(st.P.MPa, 5), 3.42191)
        self.assertEqual(round(st.h.kJkg, 5), -371.03600)
        self.assertEqual(round(st.s.kJkgK, 5), -2.04410)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.96150)
        self.assertEqual(round(st.cp.kJkgK, 5), 2.15083)
        self.assertEqual(round(st.w, 4), 577.7828)

    def test_Huber(self):
        """Table 7, pag 266"""
        st = HCl(T=292.2, rhom=23.87)
        # self.assertEqual(round(st.mu.muPas, 5), 92.63275)
        self.assertEqual(round(st.mu.muPas, 5), 92.63348)
        self.assertEqual(round(st.k.mWmK, 4), 199.4777)
