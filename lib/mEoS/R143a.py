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
from lib.mEoS import R134a


class R143a(MEoS):
    """Multiparameter equation of state for R143a"""
    name = "1,1,1-trifluoroethane"
    CASNumber = "420-46-2"
    formula = "CF3CH3"
    synonym = "R143a"
    _refPropName = "R143A"
    _coolPropName = "R143a"
    rhoc = unidades.Density(431.00006645)
    Tc = unidades.Temperature(345.857)
    Pc = unidades.Pressure(3761.0, "kPa")
    M = 84.041  # g/mol
    Tt = unidades.Temperature(161.34)
    Tb = unidades.Temperature(225.909)
    f_acent = 0.2615
    momentoDipolar = unidades.DipoleMoment(2.34, "Debye")
    id = 243

    Fi1 = {"ao_log": [1, -1],
           "pow": [0, 1, -0.33],
           "ao_pow": [5.903087, 7.307253, -16.59105],
           "ao_exp": [4.4402, 3.7515],
           "titao": [1791/Tc, 823/Tc]}

    Fi2 = {"R": 8.31451,
           "ao_log": [1, -0.8999794],
           "pow": [0, 1, -1.5, 1.25, -1],
           "ao_pow": [-5.556942, 8.93748, 1.652398, -0.6827433, -8.113464]}

    CP1 = {"ao": 1.838736,
           "an": [3.01994e-2, -1.78455e-5, 4.42442e-9],
           "pow": [1, 2, 3]}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-143a of Lemmon and "
                    "Jacobsen (2000)",
        "__doi__": {"autor": "Lemmon, E.W., Jacobsen, R.T.",
                    "title": "An International Standard Formulation for the "
                             "Thermodynamic Properties of 1,1,1-"
                             "Trifluoroethane (HFC-143a) for Temperatures From"
                             " 161 to 450 K and Pressures to 50 MPa",
                    "ref": "J. Phys. Chem. Ref. Data 29(4) (2000) 521-552",
                    "doi": "10.1063/1.1318909"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": {"Tref": 273.15, "Pref": 1., "ho": 33936.4, "so": 198.961},

        "Tmin": Tt, "Tmax": 650.0, "Pmax": 100000.0, "rhomax": 15.85,

        "nr1": [.77736443e1, -.870185e1, -.27779799, .1460922, .89581616e-2],
        "d1": [1, 1, 1, 2, 5],
        "t1": [0.67, 0.833, 1.7, 1.82, 0.35],

        "nr2": [-0.20552116, 0.10653258, 0.23270816e-1, -0.13247542e-1,
                -0.42793870e-1, 0.36221685, -0.25671899, -0.92326113e-1,
                0.83774837e-1, 0.17128445e-1, -0.17256110e-1, 0.49080492e-2],
        "d2": [1, 3, 5, 7, 1, 2, 2, 3, 4, 2, 3, 5],
        "t2": [3.9, 0.95, 0, 1.19, 7.2, 5.9, 7.65, 7.5, 7.45, 15.5, 22, 19],
        "c2": [1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3],
        "gamma2": [1]*12}

    outcalt = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for R-143a of Outcalt and "
                    "McLinden (1996)",
        "__doi__": {"autor": "Outcalt, S.L., McLinden, M.O.",
                    "title": "An Equation of State for the Thermodynamic "
                             "Properties of R143a (1,1,1-Trifluoroethane)",
                    "ref": "Int. J. Thermophys., 18(6) (1997) 1445-1463",
                    "doi": "10.1007/BF02575344"},

        "R": 8.314471,
        "M": 84.041, "Tc": 346.04, "Pc": 3775.6, "rhoc": 5.151118,

        "cp": CP1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 60000.0, "rhomax": 15.84,

        "b": [None, -0.240561786316e-1, 0.262345913719e1, -0.650858041394e2,
              0.995952053681e4, -0.147536464961e7, 0.135498153308e-2,
              -0.281726617426e1, 0.134371062574e4, 0.850286316514e6,
              -0.180516636446e-3, 0.618889066246, -0.223083798271e3,
              -0.119095922349e-1, -0.173933336877e1, -0.420847601180e3,
              0.213502079796, -0.565708555185e-2, 0.185442296800e1,
              -0.520377059921e-1, -0.846735696108e6, -0.207964483848e8,
              -0.349977290513e5, 0.576427827667e9, -0.389131863941e3,
              0.103074054089e5, -0.447627052215e1, -0.106673161101e6,
              -0.219511369081e-1, 0.642186519493e1, -0.938317030843e-4,
              -0.478594713528e-1, -0.206555883874e1]}

    li = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-143a of Li (1999)",
        "__doi__": {
            "autor": "Li, J., Tillner-Roth, R., Sato, H., Watanabe, K.",
            "title": "An Equation of State for 1,1,1-Trifluoroethane (R-143a)",
            "ref": "Int. J. Thermophys., 20(6) (1999) 1639-1651",
            "doi": "10.1023/A:1022645626800"},

        "R": 8.31451,
        "cp": Fi2,
        "ref": "IIR",
        "Tc": 345.86, "rhoc": 434/M,

        "Tmin": Tt, "Tmax": 650.0, "Pmax": 50000.0, "rhomax": 15.84,

        "nr1": [0.01606645, 4.163515, -5.031058, -0.01920208, 0.001470093],
        "d1": [5, 1, 1, 2, 4],
        "t1": [0, 0.5, 0.75, 2.5, 2.5],

        "nr2": [0.1775429, -0.7316069e-2, -0.9555916e-1, -0.5822518,
                -0.4211022e-3, -0.2059847e-1, 0.3711325e-1, 0.1799723e-3,
                -0.4145922e-1, 0.7682566e-4, -0.2089695e-2, 0.1958633e-2,
                -0.3198325e-5, -0.5376561e-2],
        "d2": [3, 8, 3, 1, 10, 1, 4, 8, 2, 12, 8, 2, 5, 3],
        "t2": [0.25, 0.25, 2, 3, 3, 8, 8, 8, 10, 8, 17, 20, 35, 27],
        "c2": [1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4],
        "gamma2": [1]*14}

    shortSpan = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-143a of Span and "
                    "Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of State for Technical Applications. "
                             "III. Results for Polar Fluids",
                    "ref": "Int. J. Thermophys., 24(1) (2003) 111-162",
                    "doi": "10.1023/A:1022362231796"},

        "R": 8.31451,
        "cp": Fi2,
        "ref": "IIR",
        "M": 80.04, "Tc": 345.86, "Pc": 3764, "rhoc": 434.1/80.04,

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 15.82,

        "nr1": [1.0306886, -2.9497307, 0.6943523, 0.071552102, 1.9155982e-4],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.79764936e-1, 0.56859424, -0.90946566e-2, -0.24199452,
                -0.70610813e-1, -0.75041709e-1, -0.16411241e-1],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    eq = lemmon, outcalt, li, shortSpan
    _PR = [0.0889, -22.5000]

    _surface = {
        "__doi__": {
            "autor": "Mulero, A., Cachadiña, I.",
            "title": "Recommended Correlations for the Surface Tension of "
                     "Several Fluids Included in the REFPROP Program",
            "ref": "J. Phys. Chem. Ref. Data 43(2) (2014) 023104",
            "doi": "10.1063/1.4878755"},
        "sigma": [0.0537], "exp": [1.25]}

    _surface = {"sigma": [0.05416], "exp": [1.255]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.73938e1, 0.19948e1, -0.18487e1, -0.41927e1, 0.14862e1],
        "t": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.21135e1, 0.10200e2, -0.30836e2, 0.39909e2, -0.18557e2],
        "t": [0.348, 1.6, 2.0, 2.4, 2.7]}
    _vapor_Density = {
        "eq": 2,
        "n": [-2.8673, -6.3818, -16.314, -45.947, -1.3956, -246.71],
        "t": [0.384, 1.17, 3.0, 6.2, 7.0, 15.0]}

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

              "ek": 301.76, "sigma": 0.4827, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 0.992,

              "psi": [0.942896, 0.0142114], "psi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.193e-9, "gam0": 0.055, "qd": 0.23e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )

    thermo0 = {"__name__": "Huber (2018)",

               "__doi__": {
                   "autor": "Huber, M.L.",
                   "title": "Models for Viscosity, Thermal Conductivity, and "
                            "Surface Tension of Selected Pure Fluids as "
                            "Implemented in REFPROP v10.0",
                   "ref": "NISTIR 8209",
                   "doi": "10.6028/NIST.IR.8209"},

               "eq": 1,

               "Toref": 1, "koref": 1,
               "no": [-7.00852e-3, 6.56307e-5, 2.62499e-8],
               "to": [0, 1, 2],

               "Tref_res": 345.857, "rhoref_res": 431.00006645, "kref_res": 1,
               "nr": [-0.0812212, -0.0166652, 0.0874477, -0.0351468, 0.0039957,
                      0.0762355, -0.0227662, -0.0175726, 0.00379467,
                      0.000776919],
               "tr": [0, 0, 0, 0, 0, -1, -1, -1, -1, -1],
               "dr": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
               "Xio": 0.198e-9, "gam0": 0.054, "qd": 0.588e-9, "Tcref": Tc*1.5}

    _thermal = (thermo0, )


class Test(TestCase):
    """Testing"""

    def test_lemmon(self):
        """Selected points from Table, Pag 541, saturation state"""
        st = R143a(T=R143a.Tt, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 0.00107)
        self.assertEqual(round(st.Liquido.rho, 1), 1330.5)
        self.assertEqual(round(st.Liquido.h.kJkg, 3), 52.520)
        self.assertEqual(round(st.Liquido.s.kJkgK, 5), 0.31417)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 0.8138)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 1.211)
        self.assertEqual(round(st.Liquido.w, 1), 1058.1)
        self.assertEqual(round(st.Gas.rho, 5), 0.06754)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 319.59)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.9695)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 0.5283)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 0.6299)
        self.assertEqual(round(st.Gas.w, 1), 137.6)

        st = R143a(T=-50+273.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 0.08874)
        self.assertEqual(round(st.Liquido.rho, 1), 1173.9)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 130.05)
        self.assertEqual(round(st.Liquido.s.kJkgK, 5), 0.71968)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 0.8608)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 1.318)
        self.assertEqual(round(st.Liquido.w, 1), 763.9)
        self.assertEqual(round(st.Gas.rho, 4), 4.2098)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 358.58)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.7438)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 0.7046)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 0.8331)
        self.assertEqual(round(st.Gas.w, 1), 154.2)

        st = R143a(T=273.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 0.61967)
        self.assertEqual(round(st.Liquido.rho, 1), 1024.3)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 200.00)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.0000)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 0.9408)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 1.495)
        self.assertEqual(round(st.Liquido.w, 1), 524.3)
        self.assertEqual(round(st.Gas.rho, 3), 27.306)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 387.81)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.6876)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 0.8756)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 1.109)
        self.assertEqual(round(st.Gas.w, 1), 153.1)

        st = R143a(T=50+273.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 2.30735)
        self.assertEqual(round(st.Liquido.rho, 2), 802.97)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 283.90)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.2748)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 1.051)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 2.118)
        self.assertEqual(round(st.Liquido.w, 1), 262.7)
        self.assertEqual(round(st.Gas.rho, 2), 120.31)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 402.43)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.6416)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 1.093)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 2.070)
        self.assertEqual(round(st.Gas.w, 1), 127.9)

        st = R143a(T=70+273.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 3.55268)
        self.assertEqual(round(st.Liquido.rho, 2), 600.85)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 333.20)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.4172)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 1.198)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 7.720)
        self.assertEqual(round(st.Liquido.w, 1), 122.4)
        self.assertEqual(round(st.Gas.rho, 2), 270.10)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 385.42)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.5694)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 1.272)
        self.assertEqual(round(st.Gas.cp.kJkgK, 2), 11.50)
        self.assertEqual(round(st.Gas.w, 1), 104.2)

        st = R143a(P=1e5, x=0.5)
        self.assertEqual(round(st.T.C, 3), -47.518)
        self.assertEqual(round(st.Liquido.rho, 1), 1167.1)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 133.33)
        self.assertEqual(round(st.Liquido.s.kJkgK, 5), 0.73426)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 0.8644)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 1.324)
        self.assertEqual(round(st.Liquido.w, 1), 752.1)
        self.assertEqual(round(st.Gas.rho, 4), 4.7099)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 360.15)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.7395)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 0.7123)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 0.8433)
        self.assertEqual(round(st.Gas.w, 1), 154.5)

        st = R143a(P=3.6e6, x=0.5)
        self.assertEqual(round(st.T.C, 3), 70.629)
        self.assertEqual(round(st.Liquido.rho, 2), 585.69)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 335.85)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.4247)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 1.213)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 9.788)
        self.assertEqual(round(st.Liquido.w, 1), 116.5)
        self.assertEqual(round(st.Gas.rho, 2), 283.61)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 383.21)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.5624)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 1.284)
        self.assertEqual(round(st.Gas.cp.kJkgK, 2), 14.81)
        self.assertEqual(round(st.Gas.w, 1), 103.0)

        # Selected point from Table, pag 544, Single phase region
        st = R143a(T=-100+273.15, P=1e5)
        self.assertEqual(round(st.rho, 1), 1302.0)
        self.assertEqual(round(st.u.kJkg, 3), 66.841)
        self.assertEqual(round(st.h.kJkg, 3), 66.918)
        self.assertEqual(round(st.s.kJkgK, 5), 0.39984)
        self.assertEqual(round(st.cv.kJkgK, 4), 0.8116)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.220)
        self.assertEqual(round(st.w, 1), 1002.1)

        st = R143a(T=-30+273.15, P=2e5)
        self.assertEqual(round(st.rho, 4), 8.9809)
        self.assertEqual(round(st.u.kJkg, 2), 349.27)
        self.assertEqual(round(st.h.kJkg, 2), 371.54)
        self.assertEqual(round(st.s.kJkgK, 4), 1.7237)
        self.assertEqual(round(st.cv.kJkgK, 4), 0.7650)
        self.assertEqual(round(st.cp.kJkgK, 4), 0.9139)
        self.assertEqual(round(st.w, 1), 156.6)

        st = R143a(T=300+273.15, P=5e5)
        self.assertEqual(round(st.rho, 4), 8.8753)
        self.assertEqual(round(st.u.kJkg, 2), 690.02)
        self.assertEqual(round(st.h.kJkg, 2), 746.36)
        self.assertEqual(round(st.s.kJkgK, 4), 2.5782)
        self.assertEqual(round(st.cv.kJkgK, 3), 1.274)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.378)
        self.assertEqual(round(st.w, 1), 246.1)

        st = R143a(T=10+273.15, P=1e6)
        self.assertEqual(round(st.rho, 2), 990.23)
        self.assertEqual(round(st.u.kJkg, 2), 214.20)
        self.assertEqual(round(st.h.kJkg, 2), 215.21)
        self.assertEqual(round(st.s.kJkgK, 4), 1.0533)
        self.assertEqual(round(st.cv.kJkgK, 4), 0.9584)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.548)
        self.assertEqual(round(st.w, 1), 477.5)

        st = R143a(T=50+273.15, P=2e6)
        self.assertEqual(round(st.rho, 3), 91.768)
        self.assertEqual(round(st.u.kJkg, 2), 391.16)
        self.assertEqual(round(st.h.kJkg, 2), 412.95)
        self.assertEqual(round(st.s.kJkgK, 4), 1.6833)
        self.assertEqual(round(st.cv.kJkgK, 3), 1.030)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.565)
        self.assertEqual(round(st.w, 1), 140.6)

        st = R143a(T=-100+273.15, P=3e6)
        self.assertEqual(round(st.rho, 1), 1306.3)
        self.assertEqual(round(st.u.kJkg, 3), 66.129)
        self.assertEqual(round(st.h.kJkg, 3), 68.426)
        self.assertEqual(round(st.s.kJkgK, 5), 0.39571)
        self.assertEqual(round(st.cv.kJkgK, 4), 0.8137)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.216)
        self.assertEqual(round(st.w, 1), 1013.2)

        st = R143a(T=300+273.15, P=5e6)
        self.assertEqual(round(st.rho, 3), 93.689)
        self.assertEqual(round(st.u.kJkg, 2), 675.52)
        self.assertEqual(round(st.h.kJkg, 2), 728.89)
        self.assertEqual(round(st.s.kJkgK, 4), 2.3254)
        self.assertEqual(round(st.cv.kJkgK, 3), 1.289)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.452)
        self.assertEqual(round(st.w, 1), 238.9)

        st = R143a(T=273.15, P=1e7)
        self.assertEqual(round(st.rho, 1), 1068.3)
        self.assertEqual(round(st.u.kJkg, 2), 192.37)
        self.assertEqual(round(st.h.kJkg, 2), 201.73)
        self.assertEqual(round(st.s.kJkgK, 5), 0.97356)
        self.assertEqual(round(st.cv.kJkgK, 4), 0.9401)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.411)
        self.assertEqual(round(st.w, 1), 620.1)

        st = R143a(T=-100+273.15, P=2e7)
        self.assertEqual(round(st.rho, 1), 1329.3)
        self.assertEqual(round(st.u.kJkg, 3), 62.387)
        self.assertEqual(round(st.h.kJkg, 3), 77.432)
        self.assertEqual(round(st.s.kJkgK, 5), 0.37323)
        self.assertEqual(round(st.cv.kJkgK, 4), 0.8229)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.199)
        self.assertEqual(round(st.w, 1), 1074.5)

        st = R143a(T=40+273.15, P=5e7)
        self.assertEqual(round(st.rho, 1), 1098.8)
        self.assertEqual(round(st.u.kJkg, 2), 226.53)
        self.assertEqual(round(st.h.kJkg, 2), 272.04)
        self.assertEqual(round(st.s.kJkgK, 4), 1.0876)
        self.assertEqual(round(st.cv.kJkgK, 3), 1.013)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.373)
        self.assertEqual(round(st.w, 1), 761.6)

        st = R143a(T=300+273.15, P=1e8)
        self.assertEqual(round(st.rho, 2), 890.50)
        self.assertEqual(round(st.u.kJkg, 2), 569.53)
        self.assertEqual(round(st.h.kJkg, 2), 681.83)
        self.assertEqual(round(st.s.kJkgK, 4), 1.9149)
        self.assertEqual(round(st.cv.kJkgK, 3), 1.368)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.595)
        self.assertEqual(round(st.w, 1), 671.2)

    def test_shortSpan(self):
        """Table III, Pag 117"""
        # The values don't work, I think paper problem

        # st = R143a(T=500, rho=500, eq="shortSpan")
        # self.assertEqual(round(st.cp0.kJkgK, 4), 1.2785)
        # self.assertEqual(round(st.P.MPa, 3), 20.152)
        # self.assertEqual(round(st.cp.kJkgK, 4), 1.6702)

        # st2 = R143a(T=600, rho=100, eq="shortSpan")
        # self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 201.13)
        # self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 0.47846)

    def test_Huber(self):
        """Table 7, pag 266"""
        self.assertEqual(round(
            # R143a(T=311.3, rhom=10.627).mu.muPas, 4), 103.2252)
            R143a(T=311.3, rhom=10.627).mu.muPas, 4), 103.2256)

        # Table 9, pag 271
        # Critical enhancement give a really high value in paper, maybe a typo
        # there or but in pychemqt
        # self.assertEqual(round(
        #     R143a(T=311.27, rhom=10.6289).k.mWmK, 4), 65.5116)
