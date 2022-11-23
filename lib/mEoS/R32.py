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


class R32(MEoS):
    """Multiparameter equation of state for R32"""
    name = "difluoromethane"
    CASNumber = "75-10-5"
    formula = "CH2F2"
    synonym = "R32"
    _refPropName = "R32"
    _coolPropName = "R32"
    rhoc = unidades.Density(424.)
    Tc = unidades.Temperature(351.255)
    Pc = unidades.Pressure(5782., "kPa")
    M = 52.024  # g/mol
    Tt = unidades.Temperature(136.34)
    Tb = unidades.Temperature(221.499)
    f_acent = 0.2769
    momentoDipolar = unidades.DipoleMoment(1.978, "Debye")
    id = 645

    Fi1 = {"R": 8.314471,
           "ao_log": [1, 3.004486],
           "pow": [0, 1],
           "ao_pow": [-8.258096, 6.353098],
           "ao_exp": [1.160761, 2.645151, 5.794987, 1.129475],
           "titao": [2.2718538, 11.9144210, 5.1415638, 32.7682170]}

    Fi2 = {"ao_log": [1, 2.999660],
           "pow": [0, 1],
           "ao_pow": [-8.253834, 6.351918],
           "ao_exp": [3.12115, 0.9994221, 2.412721, 3.055435],
           "titao": [4.559777, 2.164788, 1.234687e1, 5.877902]}

    CP2 = {"ao": 36.79959/8.314471,
           "an": [-0.06304821/8.314471, 3.757936e-4/8.314471,
                  -3.219812e-7/8.314471],
           "pow": [1, 2, 3]}

    # Expression in tau term, dividing by Tc in all terms
    CP3 = {"ao": 4.3914,
           "an": [-2.5143/351.35, 5.3885/351.35**2, -1.6057/351.35*3],
           "pow": [1, 2, 3]}

    tillner = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-32 of Tillner-Roth "
                    "and Yokozeki (1997)",
        "__doi__": {"autor": "Tillner-Roth, R., Yokozeki, A.",
                    "title": "An International Standard Equation of State for "
                             "Difluoromethane (R-32) for Temperatures from "
                             "the Triple Point at 136.34 K to 435 K at "
                             "Pressures up to 70 MPa",
                    "ref": "J. Phys. Chem. Ref. Data 26(6) (1997) 1273-1328",
                    "doi": "10.1063/1.556002"},

        "R": 8.314471,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 435.0, "Pmax": 70000.0, "rhomax": 27.4734,

        "nr1": [0.1046634e1, -0.5451165, -0.2448595e-2, -0.4877002e-1,
                0.3520158e-1, 0.1622750e-2, 0.2377225e-4, 0.2914900e-1],
        "d1": [1, 2, 5, 1, 1, 3, 8, 4],
        "t1": [0.25, 1., -0.25, -1., 2., 2., 0.75, 0.25],

        "nr2": [0.3386203e-2, -0.4202444e-2, 0.4782025e-3, -0.5504323e-2,
                -0.2418396e-1, 0.4209034, -0.4616537, -0.1200513e1,
                -0.2591550e1, -0.1400145e1,  0.8263017],
        "d2": [4, 4, 8, 3, 5, 1, 1, 3, 1, 2, 3],
        "t2": [18., 26., -1., 25., 1.75, 4., 5., 1., 1.5, 1., 0.5],
        "c2": [4, 3, 1, 4, 1, 2, 2, 1, 1, 1, 1],
        "gamma2": [1]*11}

    shortSpan = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-32 of Span and "
                    "Wagner (2003).",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of State for Technical Applications. "
                             "III. Results for Polar Fluids",
                    "ref": "Int. J. Thermophys., 24(1) (2003) 111-162",
                    "doi": "10.1023/A:1022362231796"},

        "R": 8.31451,
        "cp": Fi1,
        "ref": "IIR",
        "Tc": 351.35, "rhoc": 427/M, "Pc": 5795,

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 27.41,

        "nr1": [0.92876414, -2.4673952, 0.40129043, 0.055101049, 1.1559754e-4],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [-0.25209758, 0.42091879, 0.0037071833, -0.10308607,
                -0.11592089, -0.044350855, -0.012788805],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    astina = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-32 of Astina (2003)",
        "__doi__": {"autor": "Astina, I.M., Sato, H.",
                    "title": "A Rational Helmholtz Fundamental Equation of "
                             "State for Difluoromethane with an Intermolecular"
                             " Potential Background",
                    "ref": "Int. J. Thermophys. 24(4) (2003) 963-990",
                    "doi": "10.1023/A:1025096716493"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 450.0, "Pmax": 72000.0, "rhomax": 27.48,

        "nr1": [2.118688, -4.531096, 1.442456, 2.053906e-1, -1.311675e-1,
                1.022272e-2],
        "d1": [1, 1, 1, 3, 3, 4],
        "t1": [0.5, 1.125, 1.625, 0.875, 1.5, 1.75],

        "nr2": [4.873982e-1, -1.062213, -4.542051e-3, -6.933347e-4,
                -3.510307e-2, -5.606161e-2, 8.849625e-2, -1.850758e-2,
                7.878071e-3, -3.384115e-2, 1.641979e-4, -1.459172e-3],
        "d2": [1, 1, 5, 5, 6, 1, 2, 5, 6, 2, 2, 8],
        "t2": [1.75, 2.75, 0.25, 3.75, 1, 6.5, 2.5, 7.5, 7.5, 11, 16, 13],
        "c2": [1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3],
        "gamma2": [1]*12}

    outcalt = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for R-32 of Outcalt and McLinden "
                    "(1995)",
        "__doi__": {"autor": "Outcalt, S.L., McLinden, M.O.",
                    "title": "Equations of State for the Thermodynamic "
                             "Properties of R32 (Difluoromethane) and R125 "
                             "(Pentafluoroethane)",
                    "ref": "Int. J. Thermophysics 16(1) (1995) 79-89.",
                    "doi": "10.1007/BF01438959"},

        "R": 8.314471,
        "Tc": 351.35, "Pc": 5795, "rhoc": 8.2078,

        "cp": CP2,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 60000.0, "rhomax": 27.48,

        "b": [None, -0.131275405202e-3, 0.899927934911, -0.281400805178e2,
              0.436091182784e4, -0.837235280004e6, -0.782176408963e-6,
              -0.111226606825e1, 0.539331431878e3, 0.288600276863e6,
              -0.352264609289e-4, 0.189661830119, -0.686549003993e2,
              -0.349007064245e-2, -0.749983559476e-1, -0.321524283063e2,
              0.913057921906e-2, -0.171082181849e-3, 0.503986984347e-1,
              -0.830354867752e-3, -0.245522676708e6, -0.107859056038e8,
              -0.429514279646e4, 0.808724729567e8, -0.125945229993e2,
              -0.105735009761e4, -0.904064745354e-1, -0.183578733048e4,
              -0.169690612464e-3, 0.639250820631e-1, -0.204925767440e-6,
              -0.165629700870e-3, -0.932607493424e-2]}

    sun = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-32 of Sun and Ely "
                    "(2004)",
        "__doi__": {"autor": "Sun, L., Ely, J.F.",
                    "title": "Universal equation of state for engineering "
                             "application: Algorithm and  application to "
                             "non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223 (2004) 107-118",
                    "doi": "10.1016/j.fluid.2004.06.028"},

        "R": 8.314471,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40.,

        "nr1": [2.75866232e-1, 9.26526641e-1, -2.44296579, 5.34289357e-2,
                1.06739638e-4, 3.46487335e-2],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [9.07435007e-2, -1.93104843e-1, 5.11370826e-1, 3.09453923e-3,
                -1.53328967e-1, -1.03816916e-1, -3.8066998e-2, -1.16075825e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    vasserman = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-32 of Vasserman and"
                    "Fominsky (2001)",
        "__doi__": {"autor": "Vasserman A.A., Fominsky D.V.",
                    "title": "Equations of State for the Ozone-Safe "
                             "Refrigerants R32 and R125",
                    "ref": "Int. J. Thermophysics 22(4) (2001) 1089-1098",
                    "doi": "10.1023/a_1010699806169"},

        "R": 0.159821*M,
        "Tc": 351.35, "rhoc": 427/M,
        "cp": CP3,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40.,

        "nr1": [1.183486, -2.430934, -1.472179e-2, -4.506743e-1, 1.721527,
                -1.349166, -6.052212e-1, 9.265910e-1, 8.081905e-2,
                -1.999587e-1, 3.655934e-3, 8.217181e-3, -3.230880e-3,
                5.778584e-3, -2.536027e-6],
        "d1": [1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 6, 6, 10],
        "t1": [0, 1, 4, 0, 1, 2, 1, 2, 0, 2, 3, 2, 0, 1, 0],

        "nr2": [-6.546357e-2, -2.784785e-1, 1.113400, -2.954417, 4.898234,
                -2.354906, -7.709682e-1, 6.502963e-1, 2.168338e-1,
                -5.499117e-1, 1.978099e-2, 9.535163e-2, -1.425744e-2,
                3.921874e-3],
        "d2": [1, 1, 2, 2, 2, 2, 3, 4, 5, 5, 6, 6, 8, 9],
        "t2": [4, 5, 1, 2, 4, 5, 5, 5, 3, 4, 3, 5, 4, 2],
        "c2": [2]*14,
        "gamma2": [1]*14}

    eq = tillner, outcalt, shortSpan, astina, vasserman, sun
    _PR = [0.18086, -25.5000]

    _surface = {
        "__doi__": {
            "autor": "Tanaka, K., Higashi, Y.",
            "title": "Surface Tensions of trans-1,3,3,3-Tetrafluoropropene "
                     "and trans-1,3,3,3-Tetrafluoropropene + Difluoromethane "
                     "Mixture",
            "ref": "J. Chem. Eng. Japan 46(6) (2013) 371-375",
            "doi": "10.1252/jcej.13we021"},
        "sigma": [0.07216], "exp": [1.252]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.74883e1, 0.19697e1, -0.17496e1, -0.40224e1, 0.15209e1],
        "t": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.12584e1, 0.46410e1, -0.54870e1, 0.33115e1, -0.61370],
        "t": [0.27, 0.8, 1.1, 1.5, 1.8]}
    _vapor_Density = {
        "eq": 2,
        "n": [-.22002e1, -.5972e1, -.14571e2, -.42598e2, .42686e1, -.73373e2],
        "t": [0.336, 0.98, 2.7, 5.7, 6.5, 11.0]}

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

              "ek": 289.65, "sigma": 0.4098, "omega": 5,

              "psi": [0.7954, 5.42658e-2], "psi_d": [0, 1],
              "fint": [4.36654e-4, 1.78134e-6], "fint_t": [0, 1],
              "chi": [1.2942, -9.24549e-2], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
              "Xio": 0.194e-9, "gam0": 0.0496, "qd": 5e-10, "Tcref": 1.5*Tc}

    _viscosity = trnECS,
    _thermal = trnECS,


class Test(TestCase):

    def test_tillner(self):
        # Selected point from Table 12, Pag 1293, saturation state

        st = R32(T=R32.Tt, x=0.5)
        self.assertEqual(round(st.P.kPa, 2), 0.05)
        self.assertEqual(round(st.Liquido.rho, 1), 1429.3)
        self.assertEqual(round(st.Gas.rho, 4), 0.0022)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -19.07)
        self.assertEqual(round(st.Hvap.kJkg, 2), 463.38)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 444.31)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), -0.105)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 3.2937)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 1.592)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 0.660)

        st = R32(T=-100+273.15, x=0.5)
        self.assertEqual(round(st.P.kPa, 2), 3.81)
        self.assertEqual(round(st.Liquido.rho, 1), 1339.0)
        self.assertEqual(round(st.Gas.rho, 4), 0.1385)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 38.83)
        self.assertEqual(round(st.Hvap.kJkg, 2), 429.48)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 468.31)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 0.2711)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 2.7515)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 1.560)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 0.703)

        st = R32(T=-50+273.15, x=0.5)
        self.assertEqual(round(st.P.kPa, 2), 110.14)
        self.assertEqual(round(st.Liquido.rho, 1), 1208.4)
        self.assertEqual(round(st.Gas.rho, 4), 3.2316)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 117.22)
        self.assertEqual(round(st.Hvap.kJkg, 2), 380.06)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 497.27)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 0.6683)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 2.3714)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 1.589)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 0.883)

        st = R32(T=273.15, x=0.5)
        self.assertEqual(round(st.P.kPa, 2), 813.10)
        self.assertEqual(round(st.Liquido.rho, 1), 1055.3)
        self.assertEqual(round(st.Gas.rho, 4), 22.091)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 200.00)
        self.assertEqual(round(st.Hvap.kJkg, 2), 315.30)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 515.30)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.0000)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 2.1543)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 1.745)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 1.251)

        st = R32(T=50+273.15, x=0.5)
        self.assertEqual(round(st.P.kPa, 1), 3141.2)
        self.assertEqual(round(st.Liquido.rho, 2), 839.26)
        self.assertEqual(round(st.Gas.rho, 3), 98.550)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 297.49)
        self.assertEqual(round(st.Hvap.kJkg, 2), 209.62)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 507.10)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.3183)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.9670)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 2.439)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 2.477)

        st = R32(T=74+273.15, x=0.5)
        self.assertEqual(round(st.P.kPa, 1), 5304.6)
        self.assertEqual(round(st.Liquido.rho, 2), 624.57)
        self.assertEqual(round(st.Gas.rho, 2), 240.12)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 367.53)
        self.assertEqual(round(st.Hvap.kJkg, 2), 98.88)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 466.41)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.5179)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.8027)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 8.052)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 12.094)

        st = R32(T=78+273.15, x=0.5)
        self.assertEqual(round(st.P.kPa, 1), 5769.7)
        self.assertEqual(round(st.Liquido.rho, 2), 484.61)
        self.assertEqual(round(st.Gas.rho, 2), 367.24)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 400.38)
        self.assertEqual(round(st.Hvap.kJkg, 2), 28.52)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 428.90)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.6095)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.6907)

        # Table 13, pag. 1295
        st = R32(P=100e3, x=0.5)
        self.assertEqual(round(st.T.C, 2), -51.91)
        self.assertEqual(round(st.Liquido.rho, 1), 1213.6)
        self.assertEqual(round(st.Gas.rho, 4), 2.9512)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 114.18)
        self.assertEqual(round(st.Hvap.kJkg, 2), 382.14)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 496.32)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 0.6546)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 2.3819)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 1.587)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 0.873)

        st = R32(P=5000e3, x=0.5)
        self.assertEqual(round(st.T.C, 2), 71.18)
        self.assertEqual(round(st.Liquido.rho, 2), 666.31)
        self.assertEqual(round(st.Gas.rho, 2), 207.51)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 356.00)
        self.assertEqual(round(st.Hvap.kJkg, 2), 120.27)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 476.27)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.4859)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.8352)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 5.428)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 7.607)

        # Selected point from Table 14, Pag 1297, single phase region
        st = R32(T=-85+273.15, P=10e3)
        self.assertEqual(round(st.rho, 4), 0.3354)
        self.assertEqual(round(st.h.kJkg, 2), 478.09)
        self.assertEqual(round(st.s.kJkgK, 4), 2.6526)

        st = R32(T=-75+273.15, P=20e3)
        self.assertEqual(round(st.rho, 4), 0.6401)
        self.assertEqual(round(st.h.kJkg, 2), 484.34)
        self.assertEqual(round(st.s.kJkgK, 4), 2.5754)

        st = R32(T=160+273.15, P=50e3)
        self.assertEqual(round(st.rho, 4), 0.7234)
        self.assertEqual(round(st.h.kJkg, 2), 687.70)
        self.assertEqual(round(st.s.kJkgK, 4), 3.0905)

        st = R32(T=-55+273.15, P=100e3)
        self.assertEqual(round(st.rho, 1), 1222.2)
        self.assertEqual(round(st.h.kJkg, 2), 109.29)
        self.assertEqual(round(st.s.kJkgK, 4), 0.6324)

        st = R32(T=273.15, P=200e3)
        self.assertEqual(round(st.rho, 4), 4.7455)
        self.assertEqual(round(st.h.kJkg, 2), 536.37)
        self.assertEqual(round(st.s.kJkgK, 4), 2.4378)

        st = R32(T=-25+273.15, P=300e3)
        self.assertEqual(round(st.rho, 4), 8.2165)
        self.assertEqual(round(st.h.kJkg, 2), 509.91)
        self.assertEqual(round(st.s.kJkgK, 4), 2.2747)

        st = R32(T=105+273.15, P=500e3)
        self.assertEqual(round(st.rho, 4), 8.4945)
        self.assertEqual(round(st.h.kJkg, 2), 627.61)
        self.assertEqual(round(st.s.kJkgK, 4), 2.5771)

        st = R32(T=5+273.15, P=1000e3)
        self.assertEqual(round(st.rho, 1), 1038.0)
        self.assertEqual(round(st.h.kJkg, 2), 208.80)
        self.assertEqual(round(st.s.kJkgK, 4), 1.0313)

        st = R32(T=-85+273.15, P=2000e3)
        self.assertEqual(round(st.rho, 1), 1303.6)
        self.assertEqual(round(st.h.kJkg, 2), 63.18)
        self.assertEqual(round(st.s.kJkgK, 4), 0.3976)

        st = R32(T=160+273.15, P=3000e3)
        self.assertEqual(round(st.rho, 3), 47.926)
        self.assertEqual(round(st.h.kJkg, 2), 662.45)
        self.assertEqual(round(st.s.kJkgK, 4), 2.3927)

        st = R32(T=70+273.15, P=5000e3)
        self.assertEqual(round(st.rho, 2), 689.76)
        self.assertEqual(round(st.h.kJkg, 2), 350.24)
        self.assertEqual(round(st.s.kJkgK, 4), 1.4692)

        st = R32(T=45+273.15, P=10000e3)
        self.assertEqual(round(st.rho, 2), 930.32)
        self.assertEqual(round(st.h.kJkg, 2), 280.50)
        self.assertEqual(round(st.s.kJkgK, 4), 1.2414)

        st = R32(T=-85+273.15, P=20000e3)
        self.assertEqual(round(st.rho, 1), 1323.9)
        self.assertEqual(round(st.h.kJkg, 2), 72.12)
        self.assertEqual(round(st.s.kJkgK, 4), 0.3723)

        st = R32(T=160+273.15, P=40000e3)
        self.assertEqual(round(st.rho, 2), 760.12)
        self.assertEqual(round(st.h.kJkg, 2), 475.86)
        self.assertEqual(round(st.s.kJkgK, 4), 1.6700)

        # Table 15, Pag 1319, isobaric heat capacity
        st = R32(T=-85+273.15, P=20e3)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.559)

        st = R32(T=160+273.15, P=6000e3)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.412)

        st = R32(T=273.15, P=50000e3)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.518)

        # Table 16, Pag 1322, isochoric heat capacity
        st = R32(T=-85+273.15, P=20e3)
        self.assertEqual(round(st.cv.kJkgK, 3), 0.970)

        st = R32(T=160+273.15, P=6000e3)
        self.assertEqual(round(st.cv.kJkgK, 3), 0.976)

        st = R32(T=273.15, P=50000e3)
        self.assertEqual(round(st.cv.kJkgK, 3), 0.949)

        # Table 17, Pag 1325, speed of sound
        st = R32(T=-85+273.15, P=20e3)
        self.assertEqual(round(st.w, 1), 1143.4)

        st = R32(T=160+273.15, P=6000e3)
        self.assertEqual(round(st.w, 2), 254.16)

        st = R32(T=273.15, P=50000e3)
        self.assertEqual(round(st.w, 2), 989.78)

    def test_shortSpan(self):
        # Table III, Pag 117
        st = R32(T=500, rho=500, eq="shortSpan")
        self.assertEqual(round(st.cp0.kJkgK, 4), 1.1419)
        self.assertEqual(round(st.P.MPa, 3), 30.358)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.8390)

        st2 = R32(T=600, rho=100, eq="shortSpan")
        self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 235.82)
        self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 0.59788)
